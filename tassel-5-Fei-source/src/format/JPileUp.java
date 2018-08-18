/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package format;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author Fei Lu
 */
public class JPileUp {
    PileUpRecord[] pileInfo;
    
    public JPileUp () {
        
    }
    
    public JPileUp (String pileupFileS, IOFileFormat format) {
        this.ReadFile(pileupFileS, format);
    }
    
    public void readFromPileup (String infileS) {
        BufferedReader br = null;
        ArrayList<PileUpRecord> pileList = new ArrayList();
        int cnt = 0;
        try {
            if (infileS.endsWith(".gz")) br = IoUtils.getTextGzipReader(infileS);
            else br = IoUtils.getTextReader(infileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                PileUpRecord pr = new PileUpRecord(temp);
                if (pr.allele == null) continue;
                pileList.add(pr);
                cnt++;
                if (cnt%1000000 == 0) System.out.println("Processed "+ String.valueOf(cnt)+" sites from " + infileS);
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        pileInfo = pileList.toArray(new PileUpRecord[pileList.size()]);
    }
    
    public int getSiteNumber () {
        return this.pileInfo.length;
    }
    
    public PileUpRecord getSiteRecord (int index) {
        return this.pileInfo[index];
    }
    
    public void writeFile (String outfileS, IOFileFormat format) {
        if (format == IOFileFormat.Binary) {
            this.writeBinaryFile(outfileS);
        }
        else if (format == IOFileFormat.Text) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
        else {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
    }
    
    private void writeBinaryFile (String outfileS) {
        try {
            DataOutputStream dos = IoUtils.getBinaryWriter(outfileS);
            dos.writeInt(this.getSiteNumber());
            for (int i = 0; i < this.getSiteNumber(); i++) {
                dos.writeInt(this.getSiteRecord(i).chromosome);
                dos.writeInt(this.getSiteRecord(i).position);
                dos.writeByte(this.getSiteRecord(i).refBase);
                dos.writeInt(this.getSiteRecord(i).totalDepth);
                dos.writeInt(this.getSiteRecord(i).getRefReadCount());
                for (int j = 0; j < this.getSiteRecord(i).getRefReadCount(); j++) {
                    dos.writeByte(this.getSiteRecord(i).getRefBaseQuality(j));
                }
                dos.writeInt(this.getSiteRecord(i).getAlleleNumber());
                for (int j = 0; j < this.getSiteRecord(i).getAlleleNumber(); j++) {
                    dos.writeUTF(this.getSiteRecord(i).allele[j]);
                    dos.writeShort(this.getSiteRecord(i).alleleCount[j]);
                    for (int k = 0; k < this.getSiteRecord(i).alleleCount[j]; k++) {
                        dos.writeByte(this.getSiteRecord(i).alleleBaseQ[j][k]);
                    }
                }
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void ReadFile (String infileS, IOFileFormat format) {
        if (format == IOFileFormat.Binary) {
            this.readBinaryFile(infileS);
        }
        else if (format == IOFileFormat.Text) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
        else {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
    }
    
    private void readBinaryFile (String infileS) {
        try {
            DataInputStream dis = IoUtils.getBinaryReader(infileS);
            int siteNumber = dis.readInt();
            this.pileInfo = new PileUpRecord[siteNumber];
            for (int i = 0; i < siteNumber; i++) {
                int chromosome = dis.readInt();
                int position = dis.readInt();
                byte refBase = dis.readByte();
                int totalDepth = dis.readInt();
                int refDepth = dis.readInt();
                byte[] refBaseQ = new byte[refDepth];
                for (int j = 0; j < refBaseQ.length; j++) {
                    refBaseQ[j] = dis.readByte();
                }
                int alleleNumber = dis.readInt();
                String[] allele = new String[alleleNumber];
                short[] alleleCount = new short[alleleNumber];
                byte[][] alleleBaseQ = new byte[alleleNumber][];
                for (int j = 0; j < alleleNumber; j++) {
                    allele[j] = dis.readUTF();
                    alleleCount[j] = dis.readShort();
                    alleleBaseQ[j] = new byte[alleleCount[j]];
                    for (int k = 0; k < alleleBaseQ[j].length; k++) {
                        alleleBaseQ[j][k] = dis.readByte();
                    }
                }
                pileInfo[i] = new PileUpRecord(chromosome, position, refBase, totalDepth, refBaseQ, allele, alleleCount, alleleBaseQ);
                if ((i+1)%5000000 == 0) System.out.println("Read in " + String.valueOf(i+1) + " sites from " + infileS);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void sortByPosition () {
        Arrays.sort(pileInfo);
    }
    
    public int getChromosome (int siteIndex) {
        return this.pileInfo[siteIndex].getChromosome();
    }
    
    public int getPosition (int siteIndex) {
        return this.pileInfo[siteIndex].getPosition();
    }
    
    public int getSiteDepth (int siteIndex) {
        return this.pileInfo[siteIndex].totalDepth;
    }
    
    public int getAlleleNumber (int siteIndex) {
        return this.pileInfo[siteIndex].getAlleleNumber();
    }
    
    public short getAlleleDepth (int siteIndex, int alleleIndex)  {
        return this.pileInfo[siteIndex].alleleCount[alleleIndex];
    }
    
    public byte getAlleleBaseQuality (int siteIndex, int alleleIndex, int readIndex) {
        return this.pileInfo[siteIndex].getAlleleBaseQuality(alleleIndex, readIndex);
    }
    
    public byte[] getAlleleBaseQuality (int siteIndex, int alleleIndex) {
        return this.pileInfo[siteIndex].getAlleleBaseQuality(alleleIndex);
    }
    
    public int getRefDepth (int siteIndex) {
        return this.pileInfo[siteIndex].refBaseQ.length;
    }
    
    public byte getRefBaseQuality (int siteIndex, int readIndex) {
        return this.pileInfo[siteIndex].refBaseQ[readIndex];
    }
    
    public byte[] getRefBaseQuality (int siteIndex) {
        return this.pileInfo[siteIndex].refBaseQ;
    }
    
    /**
     * Return site index of a position
     * @param chromosome
     * @param position
     * @return minus value if the position does not exist
     */
    public int getSiteIndex (int chromosome, int position) {
        int index = Arrays.binarySearch(pileInfo, new PileUpRecord(chromosome, position));
        return index;
    }
    
    public int getNearbySiteIndex (int chromosome, int position) {
        int index = this.getSiteIndex(chromosome, position);
        if (index > -1) {
            return index;
        }
        else if (index == -1) {
            return 0;
        }
        else {
            if (index < -this.getSiteNumber()) {
                return -index-2;
            }
            else {
                int d1 = Math.abs(position-this.getPosition(-index-1));
                int d2 = Math.abs(position-this.getPosition(-index-2));
                if (d1 <= d2) return -index-1;
                else return -index-2;
            }
        }
    }
    
    class PileUpRecord implements Comparable<PileUpRecord>{
        int chromosome = Integer.MIN_VALUE;
        int position = Integer.MIN_VALUE;
        byte refBase = Byte.MIN_VALUE;
        int totalDepth = Integer.MIN_VALUE;
        String[] allele = null;
        short[] alleleCount = null;
        byte[][] alleleBaseQ = null;
        byte[] refBaseQ = null;
        
        public PileUpRecord (int chromosome, int position, byte refBase, int totalDepth, byte[] refBaseQ, String[] allele, short[] alleleCount, byte[][] alleleBaseQ) {
            this.chromosome = chromosome;
            this.position = position;
            this.refBase = refBase;
            this.totalDepth = totalDepth;
            this.refBaseQ = refBaseQ;
            this.allele = allele;
            this.alleleCount = alleleCount;
            this.alleleBaseQ = alleleBaseQ;
        }
        
        public PileUpRecord (int chromosome, int position) {
            this.chromosome = chromosome;
            this.position = position;
        }
        
        public PileUpRecord (String inputline) {
            this.process(inputline);
        }
        
        public int getChromosome () {
            return this.chromosome;
        }
        
        public int getPosition () {
            return this.position;
        }
        
        public int getAlleleNumber () {
            return allele.length;
        }
        
        public int getRefReadCount () {
            return this.refBaseQ.length;
        }
        
        public byte getRefBaseQuality (int readIndex) {
            return this.refBaseQ[readIndex];
        }
        
        public byte[] getRefBaseQualituy () {
            return this.refBaseQ;
        }
        
        public short getAlleleReadCount (int alleleIndex) {
            return alleleCount[alleleIndex];
        }
        
        public byte getAlleleBaseQuality (int alleleIndex, int readIndex) {
            return alleleBaseQ[alleleIndex][readIndex];
        }
        
        public byte[] getAlleleBaseQuality (int alleleIndex) {
            return alleleBaseQ[alleleIndex];
        }
        
        private void process (String inputline) {
            String[] tem = inputline.split("\t");
            if (tem[2].startsWith("N") || tem[2].startsWith("n")) return;
            int nBam = (tem.length-3)/3;
            this.chromosome = Integer.valueOf(tem[0]);
            this.position = Integer.valueOf(tem[1]);
            this.refBase = tem[2].getBytes()[0];
            int cnt = 0;
            for (int i = 0; i < nBam; i++) cnt += Integer.valueOf(tem[i*3+3]);
            this.totalDepth = cnt;
            StringBuilder sbb = new StringBuilder();
            StringBuilder sbq = new StringBuilder();
            for (int i = 0; i < nBam; i++) {
                if (Integer.valueOf(tem[i*3+3]) == 0) {
                    tem[i*3+4] = tem[i*3+4].replaceAll("\\*", "");
                    tem[i*3+5] = tem[i*3+5].replaceAll("\\*", "");
                }
                sbb.append(tem[i*3+4]);
                sbq.append(tem[i*3+5]);
            }
            String bases = sbb.toString();
            String quals = sbq.toString();
            bases = bases.replaceAll("\\^.|\\$|>|<", "");
            bases = bases.replaceAll("a", "A");
            bases = bases.replaceAll("g", "G");
            bases = bases.replaceAll("t", "T");
            bases = bases.replaceAll("c", "C");
            bases = bases.replaceAll("n", "N");
            bases = bases.replaceAll(",", ".");
            bases = bases.replaceAll("\\*", ".");
            HashSet<String> alleleSet = new HashSet();
            ArrayList<String> alleleList = new ArrayList();
            try {
                char[] baseChars = bases.toCharArray();
                for (int i = 0; i < baseChars.length; i++) {
                    //A,T,G,C,N
                    if (baseChars[i]> 64 & baseChars[i] < 91) {
                        if (i+1 < baseChars.length && (baseChars[i+1] == 43 || baseChars[i+1] == 45)) {
                            String s = "";
                            for (int j = i+2; j < baseChars.length; j++) {
                                if (baseChars[j] > 47 && baseChars[j] < 58) {
                                    s=s+Character.toString(baseChars[j]);
                                }
                                else {
                                    break;
                                }
                            }
                            int ns = Integer.valueOf(s);
                            int endIndex = i+2+s.length()+ns;
                            alleleList.add(bases.substring(i, endIndex));
                            i = endIndex-1;
                        }
                        else {
                            alleleList.add(Character.toString(baseChars[i]));
                        }
                    }
                    //+,-
                    else if (baseChars[i] == 43 || baseChars[i] == 45) {

                    }
                    //.
                    else if (baseChars[i] == 46) {
                        if (i+1 < baseChars.length && (baseChars[i+1] == 43 || baseChars[i+1] == 45)) {
                            String s = "";
                            for (int j = i+2; j < baseChars.length; j++) {
                                if (baseChars[j] > 47 && baseChars[j] < 58) {
                                    s=s+Character.toString(baseChars[j]);
                                }
                                else {
                                    break;
                                }
                            }
                            int ns = Integer.valueOf(s);
                            int endIndex = i+2+s.length()+ns;
                            alleleList.add(bases.substring(i+1, endIndex));
                            i = endIndex-1;
                        }
                        else {
                            alleleList.add(Character.toString(baseChars[i]));
                        }
                    }
                }
                String[] alleleArray = alleleList.toArray(new String[alleleList.size()]);
                int refCnt = 0;
                for (int i = 0; i < alleleArray.length; i++) {
                    if (alleleArray[i].equals(".")) {
                        refCnt++;
                        continue;
                    }
                    alleleSet.add(alleleArray[i]);
                }
                this.refBaseQ = new byte[refCnt];
                char[] qualC = quals.toCharArray();
                refCnt = 0;
                for (int i = 0; i < alleleArray.length; i++) {
                    if (alleleArray[i].equals(".")) {
                        refBaseQ[refCnt] = (byte)(qualC[i]-33);
                        refCnt++;
                    }
                }
                allele = alleleSet.toArray(new String[alleleSet.size()]);
                if (allele.length == 0) return;
                Arrays.sort(allele);
                this.alleleCount = new short[allele.length];
                int[] tempCnt = new int[allele.length];
                for (int i = 0; i < alleleArray.length; i++) {
                    if (alleleArray[i].equals(".")) continue;
                    int index = Arrays.binarySearch(allele, alleleArray[i]);
                    tempCnt[index]++;
                }
                for (int i = 0; i < allele.length; i++) {
                    if (tempCnt[i] > Short.MAX_VALUE) {
                        //System.out.println(inputline);
                        System.out.println("Skip this site, chr "+String.valueOf(this.chromosome) + ", position " + String.valueOf(this.position) +". Allele read count > "+String.valueOf(Short.MAX_VALUE));
                        allele = null;
                        return;
                    }
                }
                for (int i = 0; i < allele.length; i++) alleleCount[i] = (short)tempCnt[i];
                if (allele.length == 1) {

                }
                else if (allele.length == 2) {
                    if (alleleCount[1] > alleleCount[0]) {
                        String temp = allele[0]; allele[0] = allele[1]; allele[1] = temp;
                        short tcnt = alleleCount[0]; alleleCount[0] = alleleCount[1]; alleleCount[1] = tcnt;
                    }
                }
                else {
                    Integer[] idx = new Integer[alleleCount.length];
                    for( int i = 0 ; i < idx.length; i++ ) idx[i] = i;
                    Arrays.sort(idx, new Comparator<Integer>() {
                        public int compare(Integer i1, Integer i2) {                        
                            return -Integer.compare(alleleCount[i1], alleleCount[i2]);
                        }                   
                    });
                    String[] nAllele = new String[allele.length];
                    short[] nCount = new short[allele.length];
                    for (int i = 0; i < nAllele.length; i++) {
                        nAllele[i] = allele[idx[i]];
                        nCount[i] = alleleCount[idx[i]];
                    }
                    this.allele = nAllele;
                    this.alleleCount = nCount;
                }
                this.alleleBaseQ = new byte[this.allele.length][];              
                for (int i = 0; i < allele.length; i++) {
                    this.alleleBaseQ[i] = new byte[alleleCount[i]];
                    int count = 0;
                    for (int j = 0; j < alleleArray.length; j++) {
                        if (allele[i].equals(alleleArray[j])) {
                            alleleBaseQ[i][count] = (byte)(qualC[j]-33);
                            count++;
                        }
                    }
                }
            }
            catch (Exception e) {
                e.printStackTrace();
                System.out.println(inputline);
                System.out.println(bases+"\t"+bases.length());
                System.out.println(quals+"\t"+quals.length());
                System.out.println(this.totalDepth);
                System.exit(1);
            }
        }

        @Override
        public int compareTo(PileUpRecord o) {
            if (this.chromosome == o.chromosome) {
                if (this.position == o.position) return 0;
                else if (this.position < o.position) return -1;
                else return 1;
            }
            else if (this.chromosome < o.chromosome) {
                return -1;
            }
            else {
                return 1;
            }
        }
    }
}
 