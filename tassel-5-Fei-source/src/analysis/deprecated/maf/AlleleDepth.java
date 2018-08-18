/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.maf;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author Fei Lu
 */
public class AlleleDepth implements Comparable<AlleleDepth> {
    int baseCnt = 4;
    DepthRecord[] depthInfo = null;
    
    public AlleleDepth () {
        
    }
   
    public AlleleDepth (String infileS, IOFileFormat format) {
        this.readAlleleDepth(infileS, format);
    }
    
    public AlleleDepth (DepthRecord[] info) {
        this.depthInfo = info;
    }
    
    public int getSiteNumber () {
        if (depthInfo == null) return 0;
        return depthInfo.length;
    }
    
    public byte getChromosome (int index) {
        return depthInfo[index].chr;
    }
    
    public int getPosition (int index) {
        return depthInfo[index].pos;
    }
    
    public short[] getDepthArrray (int index) {
        return depthInfo[index].depth;
    }
    
    public int getSiteDepth (int index) {
        return depthInfo[index].getSiteDepth();
    }
    
    public short getMajorDepth (int index) {
        return depthInfo[index].depth[0];
    }
    
    public short getMinorDepth (int index) {
        return depthInfo[index].depth[1];
    }
    
    public short getThirdDepth (int index) {
        return depthInfo[index].depth[2];
    }
    
    public short getForthDepth (int index) {
        return depthInfo[index].depth[3];
    }
    
    public byte[] getQualityArrray (int index) {
        return depthInfo[index].qual;
    }
    
    public byte getMajorQuality (int index) {
        return depthInfo[index].qual[0];
    }
    
    public byte getMinorQuality (int index) {
        return depthInfo[index].qual[1];
    }
    
    public byte getThirdQuality (int index) {
        return depthInfo[index].qual[2];
    }
    
    public byte getForthQuality (int index) {
        return depthInfo[index].qual[3];
    }
    
    public double getMinorAlleleFrequency (int index) {
        return (double)this.getMinorDepth(index)/(this.getMinorDepth(index)+this.getMajorDepth(index));
    }
    
    public double getCorrectedMinorAlleleFrequency (int index) {
        double v = (double)this.getMinorDepth(index)-(double)(this.getThirdDepth(index)+this.getForthDepth(index))/2;
        return v/(this.getMajorDepth(index)+v);
    }
    
    /**
     * Return the index of closest site n bp downstream of a position
     * The returned index exclusive
     * @param chr
     * @param pos
     * @param nBp
     * @return 
     */
    public int getDownstreamSiteIndexWithinNbp (int chr, int pos, int nBp) {
        DepthRecord query = new DepthRecord(chr, pos+nBp);
        int index = Arrays.binarySearch(depthInfo, query);
        if (index < 0) {
            index = -index-1;
        }
        return index;
    }
    
    /**
     * Return the index of closest site n bp upstream of a position
     * The returned index is inclusive
     * @param chr
     * @param pos
     * @param nBp
     * @return 
     */
    public int getUpstreamSiteIndexWithinNbp (int chr, int pos, int nBp) {
        DepthRecord query = new DepthRecord(chr, pos-nBp);
        int index = Arrays.binarySearch(depthInfo, query);
        if (index < 0) {
            index = -index-1;
        }
        return index;
    }
    
    /**
     * Return the site number of nbp interval around a position
     * @param chr
     * @param pos
     * @param intervelLength
     * @return 
     */
    public int getSiteNumberOfNbpInterval (int chr, int pos, int intervalLength) {
        int half = intervalLength/2;
        return this.getDownstreamSiteIndexWithinNbp(chr, pos, half)-this.getUpstreamSiteIndexWithinNbp(chr, pos, half);
    }
 
    /**
     * Return the site number of nbp interval around a position in MAF interval
     * @param chr
     * @param pos
     * @param intervalLength
     * @param mafLow
     * @param mafUp
     * @return 
     */
    public int getSiteNumberOfNbpIntervalByMAF (int chr, int pos, int intervalLength, double mafLow, double mafUp) {
        int half = intervalLength/2;
        int cnt  = 0;
        int start = this.getUpstreamSiteIndexWithinNbp(chr, pos, half);
        int end = this.getDownstreamSiteIndexWithinNbp(chr, pos, half);
        for (int i = start; i < end; i++) {
            double maf = this.getCorrectedMinorAlleleFrequency(i);
            if (maf > mafLow && maf < mafUp) cnt++;
        }
        return cnt;
    }
    
    /**
     * Return the index of a position
     * @param chr
     * @param pos
     * @return 
     */
    public int getSiteIndex (int chr, int pos) {
        DepthRecord query = new DepthRecord(chr, pos);
        int index = Arrays.binarySearch(depthInfo, query);
        return index;
    }
    
    public AlleleDepth getSubset (int[] indices) {
        DepthRecord[] dInfo = new DepthRecord[indices.length];
        for (int i = 0; i < dInfo.length; i++) {
            dInfo[i] = this.depthInfo[indices[i]];
        }
        return new AlleleDepth(dInfo);
    }
    
    public AlleleDepth merge (AlleleDepth another) {
        DepthRecord[] dInfo = new DepthRecord[this.getSiteNumber()+another.getSiteNumber()];
        for (int i = 0; i < this.getSiteNumber(); i++) dInfo[i] = this.depthInfo[i];
        for (int i = 0; i < another.getSiteNumber(); i++) {
            dInfo[i+this.getSiteNumber()] = another.depthInfo[i];
        }
        System.out.println("Another AlleleDepth merged");
        AlleleDepth newAD = new AlleleDepth(dInfo);
        newAD.sortByPosition();
        return new AlleleDepth(dInfo);       
    }
    
    public void readFromHapMap (String infileS) {
        System.out.println("Reading file from " + infileS);
        try {
            BufferedReader br = null;
            if (infileS.endsWith("gz")) {
                br = IoUtils.getTextGzipReader(infileS);
            }
            else {
                br = IoUtils.getTextReader(infileS);
            }
            ArrayList<DepthRecord> dl = new ArrayList();
            String temp = br.readLine();
            String[] tem;
            DepthQual[] info = new DepthQual[baseCnt];
            for (int i = 0; i < info.length; i++) info[i] = new DepthQual();
            int cnt = 0;
            
            while ((temp = br.readLine()) != null) {
                tem = temp.split("\t");
                for (int i = 0; i < baseCnt; i++) {
                    info[i].setValue(Integer.valueOf(tem[3+i*2]), Integer.valueOf(tem[3+i*2+1]));
                }
                Arrays.sort(info, Collections.reverseOrder());
                short[] depth = new short[baseCnt];
                byte[] qual = new byte[baseCnt];
                for (int i = 0; i < baseCnt; i++) {
                    depth[i] = info[i].depth;
                    qual[i] = info[i].qual;
                }
                dl.add(new DepthRecord(Integer.valueOf(tem[0]), Integer.valueOf(tem[1]), depth, qual));
                cnt++;
                if (cnt%10000000 == 0) System.out.println("Read in " + String.valueOf(cnt) + " sites");
            }
            depthInfo = dl.toArray(new DepthRecord[dl.size()]);
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private class DepthQual implements Comparable<DepthQual> {
        short depth = 0;
        byte qual = 0;
        
        public DepthQual () {}
        
        public void setValue (int depth, int qual) {
            if (depth > Short.MAX_VALUE) depth = Short.MAX_VALUE;
            this.depth = (short)depth;
            this.qual = (byte)qual;
        }
        
        @Override
        public int compareTo(DepthQual o) {
            return this.depth-o.depth;
        }
    }
    
    public void readAlleleDepth (String infileS, IOFileFormat format) {
        System.out.println("Reading file from " + infileS);
        if (format == IOFileFormat.Binary) {
            this.readBinaryAlleleDepth(infileS);
        }
        else if (format == IOFileFormat.Text) {
            this.readTextAlleleDepth(infileS);
        }
        else {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
        System.out.println("File read from " + infileS);
    }
    
    private void readBinaryAlleleDepth (String infileS) {
        try {
            DataInputStream dis = IoUtils.getBinaryReader(infileS);
            depthInfo = new DepthRecord[dis.readInt()];
            for (int i = 0; i < this.getSiteNumber(); i++) {
                byte chr = dis.readByte();
                int pos = dis.readInt();
                short[] depth = new short[baseCnt];
                byte[] qual = new byte[baseCnt];
                for (int j = 0; j < baseCnt; j++) depth[j] = dis.readShort();
                for (int j = 0; j < baseCnt; j++) qual[j] = dis.readByte();
                depthInfo[i] = new DepthRecord(chr, pos, depth, qual);
            }
            dis.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void readTextAlleleDepth (String infileS) {
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            ArrayList<DepthRecord> dl = new ArrayList();
            String temp = br.readLine();
            String[] tem;
            while ((temp = br.readLine()) != null) {
                tem = temp.split("\t");
                byte chr = Byte.valueOf(tem[0]);
                int pos = Integer.valueOf(tem[1]);
                short[] depth = new short[baseCnt];
                byte[] qual = new byte[baseCnt];
                for (int j = 0; j < baseCnt; j++) depth[j] = Short.valueOf(tem[2+j]);
                for (int j = 0; j < baseCnt; j++) qual[j] = Byte.valueOf(tem[2+baseCnt+j]);
                dl.add(new DepthRecord(chr, pos, depth, qual));
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void writeAlleleDepth (String outfileS, IOFileFormat format) {
        System.out.println("Writing file to " + outfileS);
        if (format == IOFileFormat.Binary) {
            this.writeBinaryAlleleDepth(outfileS);
        }
        else if (format == IOFileFormat.Text) {
            this.writeTextAlleleDepth(outfileS);
        }
        else {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
        System.out.println("File written to " + outfileS);
    }
    
    private void writeBinaryAlleleDepth (String outfileS) {
        try {
            DataOutputStream dos = IoUtils.getBinaryWriter(outfileS);
            dos.writeInt(this.getSiteNumber());
            for (int i = 0; i < this.getSiteNumber(); i++) {
                dos.writeByte(this.getChromosome(i));
                dos.writeInt(this.getPosition(i));
                for (int j = 0; j < this.getDepthArrray(i).length; j++) dos.writeShort(depthInfo[i].depth[j]);
                for (int j = 0; j < this.getQualityArrray(i).length; j++) dos.writeByte(depthInfo[i].qual[j]);
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void writeTextAlleleDepth (String outfileS) {
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Chr\tPos\tMajorDepth\tMinorDepth\tThirdDepth\tForthDepth\tMajorQual\tMinorQual\tThirdQual\tForthQual");
            bw.newLine();
            for (int i = 0; i < this.getSiteNumber(); i++) {
                bw.write(this.depthInfo[i].getRecordString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void sortByPosition () {
        Arrays.sort(this.depthInfo);
    }
    
    @Override
    public int compareTo(AlleleDepth o) {
        if (this.getChromosome(0) == o.getChromosome(0)) return 0;
        else if (this.getChromosome(0) < o.getChromosome(0)) return -1;
        else return 1;
    }
    
    class DepthRecord implements Comparable<DepthRecord> {
        byte chr = Byte.MIN_VALUE;
        int pos = Integer.MIN_VALUE;
        short[] depth = new short[baseCnt];
        byte[] qual = new byte[baseCnt];
        
        public DepthRecord (int chr, int pos) {
            this.chr = (byte)chr;
            this.pos = pos;
        }
        
        public DepthRecord (int chr, int pos, short[] depth, byte[] qual) {
            this(chr, pos);
            this.depth = depth;
            this.qual = qual;
        }

        public String getRecordString () {
            StringBuilder sb = new StringBuilder();
            sb.append(chr).append("\t").append(pos);
            for (int i = 0; i < depth.length; i++) {
                sb.append("\t").append(depth[i]);
            }
            for (int i = 0; i < qual.length; i++) {
                sb.append("\t").append(qual[i]);
            }
            return sb.toString();
        }
        
        public int getSiteDepth () {
            int dep = 0;
            for (int i = 0; i < depth.length; i++) dep+=depth[i];
            return dep;
        }
        
        @Override
        public int compareTo(DepthRecord t) {
            if (chr == t.chr) {
                if (pos == t.pos) return 0;
                else if (pos < t.pos) return -1;
                else return 1;
            }
            return chr-t.chr;
        }
    }
    
    public static AlleleDepth[] getAlleleDepthByChrome () {
        String infileDirS = "M:\\production\\maf\\wgs\\alleleDepth\\";
        File[] fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        ArrayList<AlleleDepth> adList = new ArrayList();
        fsList.parallelStream().forEach(e -> {
            adList.add(new AlleleDepth(e.getAbsolutePath(), IOFileFormat.Binary));
        });
        AlleleDepth[] ads = adList.toArray(new AlleleDepth[adList.size()]);
        Arrays.sort(ads);
        return ads;
    }
    
    public static File[] getAlleleDepthFileByChrome () {
        String infileDirS = "M:\\production\\maf\\wgs\\alleleDepth\\";
        File[] fs = new File(infileDirS).listFiles();
        return fs;
    }
}
