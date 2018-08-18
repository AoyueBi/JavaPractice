/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maf;


import com.google.common.base.Splitter;
import format.Fasta;
import format.Table;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.LongAdder;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import utils.Benchmark;
import utils.FArrayUtils;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class MafAttributeHmp3 {
    int[] chroms = null;
    int[] chromLength = null;
    HashMap<String, String[]> taxaBamPathMap = null;
    HashMap<String, String> bamPathPileupPathMap = null;
    String[] taxaNames = null;
    String[] bamPaths = null;
    Fasta genomeFa = null;
    int binSize = 100000;
    //A, C, D, G, I, T
    byte[] possibleAllele = {65, 67, 68, 71, 73, 84};
    //regular bases of a reference genome, sometimes the base is not A, C, G, or T
    byte[] bases = {65, 67, 71, 84};
    
    public MafAttributeHmp3 () {
        //this.basicInfoPipe();
        //this.localTestMafAttribute();
        this.callMafAttribute();
    }
    
    //**********************************************Local test section*******************************************************
    public void localTestMafAttribute () {
        String referenceFileS = "M:\\Database\\maize\\agpV3\\maize_10.agpV3.fa";
        String pileupResultDirS = "M:\\production\\maf\\annotations\\maf\\pipelineTest\\pileupResult_10000\\";
        String mafAttributeFileS = "M:\\production\\maf\\annotations\\maf\\pipelineTest\\mafAttribute_chr10_10000.txt";
//        String referenceFileS = "/Users/fl262/Documents/research/maf/reference/maize_10.agpV3.fa";
//        String pileupResultDirS = "/Users/fl262/Documents/research/maf/pipelineTest/pileupResult_10000/";
        
        int currentChr = 10;
        int regionStart = 1;
        int regionEnd = 10000;
        this.readReference(referenceFileS);
        int chrIndex = Arrays.binarySearch(chroms, currentChr);
        String chrSeq = this.genomeFa.getSeq(chrIndex);
        this.taxaBamPathMap = this.getTaxaBamPathMap(pileupResultDirS);
        int[][] binBound = this.creatBins(currentChr, binSize, regionStart, regionEnd);
        int count = 0;
        try{
            HashMap<String, BufferedReader> bamPathPileupReaderMap = this.getBamPathPileupReaderMapLocal();
            ConcurrentHashMap<BufferedReader, List<String>> readerRemainderMap = this.getReaderRemainderMap(bamPathPileupReaderMap);
            BufferedWriter bw = IoUtils.getTextWriter(mafAttributeFileS);
            bw.write(this.getMafAttributeHeader());
            bw.newLine();
            for (int i = 0; i < binBound.length; i++) {
                count = i;
                int binStart = binBound[i][0];
                int binEnd = binBound[i][1];
                ConcurrentHashMap<String, List<List<String>>> bamPileupResultMap = this.getBamPileupResultMap(currentChr, binStart, binEnd, bamPathPileupReaderMap, readerRemainderMap);
                StringBuilder[][] baseSb = this.getPopulateBaseBuilder(binStart, binEnd);
                int[][] depth = this.getPopulatedDepthArray(binStart, binEnd);
                this.fillDepthAndBase(bamPileupResultMap, baseSb, depth, binStart);
                String[][] base = this.getBaseMatrix(baseSb);
                ArrayList<Integer> positionList = this.getPositionList(binStart, binEnd);
                ConcurrentHashMap<Integer, String> posMafMap = new ConcurrentHashMap((int)((binEnd-binStart+1)*1.5));
                this.calculateMaf(posMafMap, positionList, currentChr, binStart, chrSeq, depth, base);
                for (int j = 0; j < positionList.size(); j++) {
                    bw.write(posMafMap.get(positionList.get(j)));
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            System.out.println(count);
            e.printStackTrace();
        }
    }
    
    private HashMap<String, BufferedReader> getBamPathPileupReaderMapLocal () {
        HashMap<String, BufferedReader> bamPathPileupReaderMap = new HashMap();
        try {
            for (int i = 0; i < this.bamPaths.length; i++) {
                String pileupFileS = bamPaths[i];
                File pileupF = new File(pileupFileS);
                if (!pileupF.exists()) {
                    BufferedReader br = new BufferedReader(new StringReader(""), 1024);
                    bamPathPileupReaderMap.put(bamPaths[i], br);
                    continue;
                }
                BufferedReader br = new BufferedReader (new FileReader(pileupF), 1024);
                bamPathPileupReaderMap.put(bamPaths[i], br);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return bamPathPileupReaderMap;
    }
    
    private ConcurrentHashMap<String, List<List<String>>> getBamPileupResultMap (String pileupResultDirS) {
        ConcurrentHashMap<String, List<List<String>>> bamPileupResultMap = new ConcurrentHashMap(2000);
        File[] fs = new File(pileupResultDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            ArrayList<List<String>> list = new ArrayList();
            try {
                BufferedReader br = IoUtils.getTextReader(f.getAbsolutePath());
                String temp;
                while ((temp = br.readLine()) != null) {
                    list.add(FStringUtils.fastSplit(temp, "\t"));
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            bamPileupResultMap.put(f.getAbsolutePath(), list);
        });
        return bamPileupResultMap;
    }
    
    /**
     * This is for local test
     * @param base
     * @param size 
     */
    private void outputBaseMatrix (String[][] base, int size) {
        String outfileS = "M:\\production\\maf\\annotations\\maf\\pipelineTest\\baseMatrix.txt";
        StringBuilder sb = new StringBuilder();
        sb.append("position");
        for (int i = 0; i < this.taxaNames.length; i++) {
            sb.append("\t").append(taxaNames[i]);
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < base.length; i++) {
                sb = new StringBuilder();
                for (int j = 0; j < base[i].length; j++) {
                    sb.append(base[i][j]).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(String.valueOf(i+1));
                bw.newLine();
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    //********************************************************************************************
    
    
    public void callMafAttribute () {
        String bamDirS = "/workdir/fl262/bam/";
        String referenceDirS = "/workdir/fl262/reference/";
        String mafAttributeDirS = "/workdir/fl262/mafAttribute/";
        String pileupDirS = "/workdir/fl262/pileup/";
        new File(mafAttributeDirS).mkdir();
        File[] dirs = new File(bamDirS).listFiles();
        Arrays.sort(dirs);
        for (int i = 0; i < dirs.length; i++) {
            String bamByChrDirS = dirs[i].getAbsolutePath();
            String mafAttributeByChrDirS = new File (mafAttributeDirS, dirs[i].getName()).getAbsolutePath();
            int chr = Integer.valueOf(new File(bamByChrDirS).getName().replaceFirst("c", ""));
            String referenceFileS = "chr"+FStringUtils.getNDigitNumber(2, chr)+".AGPV3.fa";
            referenceFileS = new File(referenceDirS, referenceFileS).getAbsolutePath();
            this.readReference(referenceFileS);
            new File(pileupDirS).mkdir();
            this.callMafAttributeByChr(bamByChrDirS, mafAttributeByChrDirS, referenceFileS, pileupDirS, chr);
            try {
                FileUtils.deleteDirectory(new File(pileupDirS));
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
    }
    
   
    
    private void callMafAttributeByChr (String bamByChrDirS, String mafAttributeByChrDirS, String referenceFileS, String pileupDirS, int currentChr) {
        System.out.println("Start calculating MAF attributes of chromosome " + String.valueOf(currentChr));
        new File (mafAttributeByChrDirS).mkdir();
        String outfileS = "chr"+FStringUtils.getNDigitNumber(3, currentChr)+"_mafAttribute.txt";
        outfileS = new File (mafAttributeByChrDirS, outfileS).getAbsolutePath();
        this.taxaBamPathMap = this.getTaxaBamPathMap(bamByChrDirS);
        this.bamPathPileupPathMap = this.getBamPathPileupPathMap(pileupDirS);
        int chrIndex = Arrays.binarySearch(chroms, currentChr);
        String chrSeq = this.genomeFa.getSeq(chrIndex).toUpperCase();
        int regionStart = 1;
        int regionEnd = chrSeq.length();
        this.performPileup(currentChr, regionStart, regionEnd, referenceFileS);
        int[][] binBound = this.creatBins(currentChr, binSize, regionStart, regionEnd);
        try {
            HashMap<String, BufferedReader> bamPathPileupReaderMap = this.getBamPathPileupReaderMap();
            ConcurrentHashMap<BufferedReader, List<String>> readerRemainderMap = this.getReaderRemainderMap(bamPathPileupReaderMap);
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write(this.getMafAttributeHeader());
            bw.newLine();
            for (int i = 0; i < binBound.length; i++) {
                long startTimePoint = System.nanoTime();
                int binStart = binBound[i][0];
                int binEnd = binBound[i][1];
                ConcurrentHashMap<String, List<List<String>>> bamPileupResultMap = this.getBamPileupResultMap(currentChr, binStart, binEnd, bamPathPileupReaderMap, readerRemainderMap);
                StringBuilder[][] baseSb = this.getPopulateBaseBuilder(binStart, binEnd);
                int[][] depth = this.getPopulatedDepthArray(binStart, binEnd);
                this.fillDepthAndBase(bamPileupResultMap, baseSb, depth, binStart);
                String[][] base = this.getBaseMatrix(baseSb);
                ArrayList<Integer> positionList = this.getPositionList(binStart, binEnd);
                ConcurrentHashMap<Integer, String> posMafMap = new ConcurrentHashMap((int)((binEnd - binStart + 1)*1.5));
                this.calculateMaf(posMafMap, positionList, currentChr, binStart, chrSeq, depth, base);
                for (int j = 0; j < positionList.size(); j++) {
                    bw.write(posMafMap.get(positionList.get(j)));
                    bw.newLine();
                }
                StringBuilder sb = new StringBuilder();
                sb.append("Bin from ").append(binStart).append(" to ").append(binEnd).append(" is finished. Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
                System.out.println(sb.toString());
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Chromosome " + String.valueOf(currentChr) + " is finished. File written to " + outfileS + "\n");
    }
    
    private void calculateMaf (ConcurrentHashMap<Integer, String> posMafMap, List<Integer> positionList, int currentChr, int startPos, String chrSeq, int[][] depth, String[][] base) {
        positionList.parallelStream().forEach(position -> {
            int index = position-startPos;
            byte refBase = (byte)(chrSeq.charAt(position-1));
            int baseIndex = Arrays.binarySearch(bases, refBase);
            if (baseIndex < 0) {
                posMafMap.put(position, this.getEmptyMafAttributeString(currentChr, position, refBase));
            }
            else {
                posMafMap.put(position, this.getMafAttributeString(base[index], depth[index], currentChr, position, refBase));
            }
        });
    }
    
    private String getMafAttributeString (String[] base, int[] depth, int currentChr, int position, byte refBase) {
        TByteArrayList bList;
        boolean ifRecordedDeletion = false;
        TIntHashSet insertionLengthSet = new TIntHashSet();
        TIntHashSet deletionLengthSet = new TIntHashSet();
        int[][] alleleCount = new int[base.length][this.possibleAllele.length];
        for (int i = 0; i < base.length; i++) {
            bList = new TByteArrayList();
            byte[] ba = base[i].getBytes();
            for (int j = 0; j < ba.length; j++) {
                if (ba[j] == '.') {
                    bList.add(refBase);
                }
                else if (ba[j] == ',') {
                    bList.add(refBase);
                }
                else if (ba[j] == 'A') {
                    bList.add((byte)65);
                }
                else if (ba[j] == 'a') {
                    bList.add((byte)65);
                }
                else if (ba[j] == 'C') {
                    bList.add((byte)67);
                }
                else if (ba[j] == 'c') {
                    bList.add((byte)67);
                }
                else if (ba[j] == 'G') {
                    bList.add((byte)71);
                }
                else if (ba[j] == 'g') {
                    bList.add((byte)71);
                }
                else if (ba[j] == 'T') {
                    bList.add((byte)84);
                }
                else if (ba[j] == 't') {
                    bList.add((byte)84);
                }
                else if (ba[j] == '+') {
                    int endIndex = j+2;
                    for (int k = j+1; k < ba.length; k++) {
                        if (ba[k] > 57) {
                            endIndex = k;
                            break;
                        }
                    }
                    StringBuilder sb = new StringBuilder();
                    for (int k = j+1; k < endIndex; k++) {
                        sb.append((char)ba[k]);
                    }
                    int length = Integer.valueOf(sb.toString());
                    insertionLengthSet.add(length);
                    j+=sb.length();
                    j+=length;
                    bList.add((byte)73);
                }
                else if (ba[j] == '-') {
                    int endIndex = j+2;
                    for (int k = j+1; k < ba.length; k++) {
                        if (ba[k] > 57) {
                            endIndex = k;
                            break;
                        }
                    }
                    StringBuilder sb = new StringBuilder();
                    for (int k = j+1; k < endIndex; k++) {
                        sb.append((char)ba[k]);
                    }
                    int length = Integer.valueOf(sb.toString());
                    deletionLengthSet.add(length);
                    j+=sb.length();
                    j+=length;
                    bList.add((byte)68);
                }
                else if (ba[j] == '^') {
                    j++;
                }
                else if (ba[j] == '*') {
                    bList.add(refBase);
                    ifRecordedDeletion = true;
                }
                //N, n, $, >, <
                else {
                    //do nothing
                }
            }
            byte[] taxonBase = bList.toArray();
            for (int j = 0; j < taxonBase.length; j++) {
                int index = Arrays.binarySearch(this.possibleAllele, taxonBase[j]);
                alleleCount[i][index]++;
            }
            
        }
        int totalDepth = 0;
        int allelePresence = 0;
        for (int i = 0; i < depth.length; i++) {
            if (depth[i] == 0) continue;
            totalDepth+=depth[i];
            allelePresence++;
        }
        int[] alleleDepth = new int[this.possibleAllele.length];
        int[] alleleOnce = new int[this.possibleAllele.length];
        int[] alleleTwice = new int[this.possibleAllele.length];
        int[] alleleThrice = new int[this.possibleAllele.length];
        int[] alleleHets = new int[this.possibleAllele.length];
        float[] hetMean = new float[this.possibleAllele.length];
        float[] hetSD = new float[this.possibleAllele.length];
        TDoubleArrayList[] ratioList = new TDoubleArrayList[this.possibleAllele.length];
        for (int i = 0; i < this.possibleAllele.length; i++) ratioList[i] = new TDoubleArrayList();
        for (int i = 0; i < alleleCount.length; i++) {
            for (int j = 0; j < alleleDepth.length; j++) {
                alleleDepth[j]+=alleleCount[i][j];
                if (alleleCount[i][j] > 0) {
                    alleleOnce[j]++;
                }
                if (alleleCount[i][j] > 1) {
                    alleleTwice[j]++;
                }
                if (alleleCount[i][j] > 2) {
                    alleleThrice[j]++;
                }
                if (alleleCount[i][j] != 0) {
                    if (alleleCount[i][j] < depth[i]) {
                        alleleHets[j]++;
                        ratioList[j].add((double)alleleCount[i][j]/depth[i]);
                    }
                }
                
            }
        }
        for (int i = 0; i < ratioList.length; i++) {
            double[] value = ratioList[i].toArray();
            DescriptiveStatistics d = new DescriptiveStatistics (value);
            hetMean[i] = (float)d.getMean();
            hetSD[i] = (float)d.getStandardDeviation();
            if (value.length == 1) {
                hetMean[i] = (float)value[0];
            }
        }
        StringBuilder sb = new StringBuilder();
        sb.append(currentChr).append("\t").append(position).append("\t").append((char)refBase).append("\t").append(totalDepth).append("\t")
            .append(allelePresence).append("\t").append(insertionLengthSet.size()).append("\t");
        if (ifRecordedDeletion) sb.append("NA");
        else sb.append(deletionLengthSet.size());
        for (int i = 0; i < this.possibleAllele.length; i++) {
            sb.append("\t").append(alleleDepth[i]);
        }
        for (int i = 0; i < this.possibleAllele.length; i++) {
            sb.append("\t").append(alleleOnce[i]);
        }
        for (int i = 0; i < this.possibleAllele.length; i++) {
            sb.append("\t").append(alleleTwice[i]);
        }
        for (int i = 0; i < this.possibleAllele.length; i++) {
            sb.append("\t").append(alleleThrice[i]);
        }
        for (int i = 0; i < this.possibleAllele.length; i++) {
            sb.append("\t").append(alleleHets[i]);
        }
        for (int i = 0; i < this.possibleAllele.length; i++) {
            sb.append("\t").append(hetMean[i]);
        }
        for (int i = 0; i < this.possibleAllele.length; i++) {
            sb.append("\t").append(hetSD[i]);
        }
        return sb.toString();
    }
    
    private String getMafAttributeHeader () {
        StringBuilder sb = new StringBuilder();
        sb.append("Chromosome\tPosition\tRefAllele\tTotalDepth\tAllelePresence\tInsertionTypeCount\tDeletionTypeCount");
        for (int i = 0; i < this.possibleAllele.length; i++) {
            sb.append("\tAlleleDepth").append((char)this.possibleAllele[i]);
        }
        for (int i = 0; i < this.possibleAllele.length; i++) {
            sb.append("\tAllelePresenceOnce").append((char)this.possibleAllele[i]);
        }
        for (int i = 0; i < this.possibleAllele.length; i++) {
            sb.append("\tAllelePresenceTwice").append((char)this.possibleAllele[i]);
        }
        for (int i = 0; i < this.possibleAllele.length; i++) {
            sb.append("\tAllelePresenceThrice").append((char)this.possibleAllele[i]);
        }
        for (int i = 0; i < this.possibleAllele.length; i++) {
            sb.append("\tAlleleHetPresence").append((char)this.possibleAllele[i]);
        }
        for (int i = 0; i < this.possibleAllele.length; i++) {
            sb.append("\tAlleleHetRatioMean").append((char)this.possibleAllele[i]);
        }
        for (int i = 0; i < this.possibleAllele.length; i++) {
            sb.append("\tAlleleHetRatioSD").append((char)this.possibleAllele[i]);
        }
        return sb.toString();
    }
    
    private String getEmptyMafAttributeString (int currentChr, int position, byte refBase) {
        StringBuilder sb = new StringBuilder();
        sb.append(currentChr).append("\t").append(position).append("\t").append((char)refBase);
        for (int i = 0; i < 46; i++) {
            sb.append("\t").append("NA");
        }
        return sb.toString();
    }
    
    private ArrayList<Integer> getPositionList (int startPos, int endPos) {
        ArrayList<Integer> positionList = new ArrayList();
        for (int i = startPos; i <= endPos; i++) {
            positionList.add(i);
        }
        return positionList;
    }
    
    private String[][] getBaseMatrix (StringBuilder[][] baseSb) {
        String[][] base = new String[baseSb.length][baseSb[0].length];
        for (int i = 0; i < base.length; i++) {
            for (int j = 0; j < base[0].length; j++) base[i][j] = baseSb[i][j].toString();
        }
        return base;
    }
    
    private void fillDepthAndBase (ConcurrentHashMap<String, List<List<String>>> bamPileupResultMap, StringBuilder[][] baseSb, int[][] depth, int startPos) {
        Set<Map.Entry<String, String[]>> entries = this.taxaBamPathMap.entrySet();
        entries.parallelStream().forEach(e -> {
            String taxa = e.getKey();
            int taxaIndex = Arrays.binarySearch(this.taxaNames, taxa);
            String[] bams = e.getValue();
            
            int count = 0;
            String b = null;
            try {
                
            
            for (int i = 0; i < bams.length; i++) {
                List<List<String>> lines = bamPileupResultMap.get(bams[i]);
                count = lines.size();
                b = bams[i];
                for (int j = 0; j < lines.size(); j++) {
                    List<String> split = lines.get(j);                 
                    if (split.get(2).startsWith("N") || split.get(2).startsWith("n")) continue;
                    int siteIndex = Integer.valueOf(split.get(1)) - startPos;
                    depth[siteIndex][taxaIndex]+=Integer.valueOf(split.get(3));
                    baseSb[siteIndex][taxaIndex].append(split.get(4));
                }
            }
            
            }
            catch (Exception ee) {
                System.out.println(b);
                System.out.println(count);
                
                ee.printStackTrace();
                System.exit(1);
            }
        });
    }
    
    private StringBuilder[][] getPopulateBaseBuilder (int startPos, int endPos) {
        StringBuilder[][] sbs = new StringBuilder[endPos-startPos+1][this.taxaNames.length];
        for (int i = 0; i < sbs.length; i++) {
            for (int j = 0; j < sbs[0].length; j++) sbs[i][j] = new StringBuilder();
        }
        return sbs;
    }
    
    private int[][] getPopulatedDepthArray (int startPos, int endPos) {
        int[][] depth = new int[endPos-startPos+1][this.taxaNames.length];
        return depth;
    }
    
    private void performPileup (int currentChr, int startPos, int endPos, String referenceFileS) {
        System.out.println("Pileup is being performed on chromosome "+String.valueOf(currentChr)+" from "+String.valueOf(startPos)+" to "+String.valueOf(endPos));
        long timeStart = System.nanoTime();
        List<String> bamList = Arrays.asList(bamPaths);
        LongAdder counter = new LongAdder();
        bamList.parallelStream().forEach(bamFileS -> {
            String pileupFileS = this.bamPathPileupPathMap.get(bamFileS);
            StringBuilder sb = new StringBuilder();
            sb.append("samtools mpileup -A -B -q 30 -Q 10 -f ").append(referenceFileS).append(" ").append(bamFileS).append(" -r ");
            sb.append(currentChr).append(":").append(startPos).append("-").append(endPos).append(" -o ").append(pileupFileS);
            String command = sb.toString();
            ArrayList<String> lineList = new ArrayList();
            try {
                Runtime rt = Runtime.getRuntime();
                Process p = rt.exec(command);
                p.waitFor();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            counter.increment();
            int cnt = counter.intValue();
            if (cnt%10 == 0) System.out.println("Pileuped " + String.valueOf(cnt) + " bam files. Total: " + String.valueOf(this.bamPaths.length));
        });
        System.out.println("Pileup is finished. Time took " + Benchmark.getTimeSpanMinutes(timeStart) + " mins");
    }
    
    private ConcurrentHashMap<String, List<List<String>>> getBamPileupResultMap (int currentChr, int binStart, int binEnd, HashMap<String, BufferedReader> bamPathPileupReaderMap, ConcurrentHashMap<BufferedReader, List<String>> readerRemainderMap) {
        ArrayList<String> empty = new ArrayList();
        ConcurrentHashMap<String, List<List<String>>> bamPileupMap = new ConcurrentHashMap(2048);
        List<String> bamList = Arrays.asList(bamPaths);
        bamList.parallelStream().forEach(bamFileS -> {
            ArrayList<List<String>> lineList = new ArrayList();
            BufferedReader br = bamPathPileupReaderMap.get(bamFileS);
            List<String> remainder = readerRemainderMap.get(br);
            boolean flag = false;
            if (remainder.size() == 0) {
                String temp = null;
                try {
                    temp = br.readLine();
                }
                catch (Exception e) {}
                if (temp != null) {
                    List<String> split = FStringUtils.fastSplit(temp, "\t");
                    int currentPos = Integer.valueOf(split.get(1));
                    if (currentPos > binEnd) {
                        readerRemainderMap.put(br, split);
                    }
                    else {
                        lineList.add(split);
                        flag = true;
                    }
                }
            }
            else {
                int currentPos = Integer.valueOf(remainder.get(1));
                if (currentPos <= binEnd) {
                    lineList.add(remainder);
                    flag = true;
                    readerRemainderMap.put(br, empty);
                }
            }
            if (flag == true) {
                try {
                    String temp;
                    while ((temp = br.readLine()) != null) {
                        List<String> split = FStringUtils.fastSplit(temp, "\t");
                        int currentPos = Integer.valueOf(split.get(1));
                        if (currentPos < binEnd) {
                            lineList.add(split);
                        }
                        else if (currentPos == binEnd){
                            lineList.add(split);
                            readerRemainderMap.put(br, empty);
                            break;
                        }
                        else {
                            readerRemainderMap.put(br, split);
                            break;
                        }
                    }
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
            bamPileupMap.put(bamFileS, lineList);
        });
        return bamPileupMap;
    }
    
    /**
     * Return regions boundaries for pileup
     * @param currentChr
     * @param binSize
     * @param regionStart
     * @param endPos
     * @return 
     */
    private int[][] creatBins (int currentChr, int binSize, int regionStart, int regionEnd) {
        int[][] binBound = FArrayUtils.getSubsetsIndicesBySubsetSize(regionEnd-regionStart+1, binSize);
        for (int i = 0; i < binBound.length; i++) {
            binBound[i][0]+=regionStart;
            binBound[i][1]+=regionStart;
            binBound[i][1]--;
        }
        System.out.println("Maf attribute calculation will performed on chromosome "+String.valueOf(currentChr)+" from "+String.valueOf(regionStart)+" to "+String.valueOf(regionEnd));
        System.out.println("Chromosome " + String.valueOf(currentChr) +"  is devided into " + String.valueOf(binBound.length) + " bins. Bin size: " + String.valueOf(this.binSize) + " bp");
        return binBound;
    }
    
    private ConcurrentHashMap<BufferedReader, List<String>> getReaderRemainderMap (HashMap<String, BufferedReader> bamPathPileupReaderMap) {
        ArrayList<String> empty = new ArrayList();
        ConcurrentHashMap<BufferedReader, List<String>> readerRemainderMap = new ConcurrentHashMap();
        Set<Map.Entry<String, BufferedReader>> enties = bamPathPileupReaderMap.entrySet();
        enties.stream().forEach(e -> {
            readerRemainderMap.put(e.getValue(), empty);
        });
        return readerRemainderMap;
    }
    
    private HashMap<String, BufferedReader> getBamPathPileupReaderMap () {
        HashMap<String, BufferedReader> bamPathPileupReaderMap = new HashMap();
        try {
            for (int i = 0; i < this.bamPaths.length; i++) {
                String pileupFileS = this.bamPathPileupPathMap.get(bamPaths[i]);
                File pileupF = new File(pileupFileS);
                if (!pileupF.exists()) {
                    BufferedReader br = new BufferedReader(new StringReader(""), 1024);
                    bamPathPileupReaderMap.put(bamPaths[i], br);
                    continue;
                }
                BufferedReader br = new BufferedReader (new FileReader(pileupF), 1024);
                bamPathPileupReaderMap.put(bamPaths[i], br);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return bamPathPileupReaderMap;
    }
    
    private HashMap<String, String> getBamPathPileupPathMap (String pileupDirS) {
        HashMap<String, String> bamPileupMap = new HashMap();
        for (int i = 0; i < this.bamPaths.length; i++) {
            String pileupFileS = new File (bamPaths[i]).getName().replaceFirst(".bam", ".pileup.txt");
            pileupFileS = new File (pileupDirS, pileupFileS).getAbsolutePath();
            bamPileupMap.put(bamPaths[i], pileupFileS);
        }
        return bamPileupMap;
    }
    
    private HashMap<String, String[]> getTaxaBamPathMap (String bamDirS)  {
        File[] fs = IoUtils.listFilesEndsWith(new File (bamDirS).listFiles(), "bam");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        String[] taxa = nameSet.toArray(new String[nameSet.size()]);
        Arrays.sort(taxa);
        ArrayList<String>[] bamList = new ArrayList[taxa.length];
        for (int i = 0; i < bamList.length; i++) bamList[i] = new ArrayList();
        for (int i = 0; i < fs.length; i++) {
            String name = fs[i].getName().split("_")[0];
            bamList[Arrays.binarySearch(taxa, name)].add(fs[i].getAbsolutePath());
        }
        HashMap<String, String[]> taxaBamMap = new HashMap();
        for (int i = 0; i < taxa.length; i++) {
            String[] bams = bamList[i].toArray(new String[bamList[i].size()]);
            Arrays.sort(bams);
            taxaBamMap.put(taxa[i], bams);
        }
        this.taxaNames = taxa;
        this.bamPaths = new String[fs.length];
        for (int i = 0; i < bamPaths.length; i++) {
            bamPaths[i] = fs[i].getAbsolutePath();
        }
        Arrays.sort(bamPaths);
        System.out.println("Taxon number: " + String.valueOf(this.taxaNames.length) + ". Bam file number: " + String.valueOf(this.bamPaths.length));
        return taxaBamMap;
    }
    
    private void readReference (String referenceGenome) {
        genomeFa = new Fasta(referenceGenome);
        genomeFa.sortRecordByNameValue();
        chroms = new int[genomeFa.getSeqNumber()];
        chromLength = new int[genomeFa.getSeqNumber()];
        for (int i = 0; i < genomeFa.getSeqNumber(); i++) {
            chroms[i] = Integer.valueOf(genomeFa.getName(i));
            chromLength[i] = genomeFa.getSeqLength(i);
        }
    }
    
    public void basicInfoPipe () {
        //this.writeTaxaBamMapFromList();
        //this.estimateCoverage();
    }
    
    public void estimateCoverage () {
        String pileupResultDirS = "M:\\production\\maf\\annotations\\maf\\pipelineTest\\test\\";
        String taxaDepthFileS = "M:\\production\\maf\\annotations\\maf\\rawDataInfo\\taxaCoverage.txt";
        this.taxaBamPathMap = this.getTaxaBamPathMap(pileupResultDirS);
        int maxSiteNum = 2000000;
        try {
            ConcurrentHashMap<String, Integer> taxaDepthMap = new ConcurrentHashMap();
            List<String> taxaList = Arrays.asList(taxaNames);
            taxaList.parallelStream().forEach(taxon -> {
                int[] depth = new int[maxSiteNum];
                String[] fs = taxaBamPathMap.get(taxon);
                for (int j = 0; j < fs.length; j++) {
                    try {
                        BufferedReader br = IoUtils.getTextReader(fs[j]);
                        String temp;
                        while ((temp = br.readLine()) != null) {
                            List<String> split = FStringUtils.fastSplit(temp, "\t");
                            depth[Integer.valueOf(split.get(1))-1]+=Integer.valueOf(split.get(3));
                        }
                    }
                    catch (Exception e) {
                        e.printStackTrace();
                    }
                    
                }
                int size = 100;
                int[] d = new int[size];
                for (int j = 0; j < depth.length; j++) {
                    if (depth[j] > size-1) continue;
                    d[depth[j]]++;
                }
                int max = -1;
                int mean = -1;
                for (int j = 1; j < size; j++) {
                    if (d[j] > max) {
                        max = d[j];
                        mean = j;
                    }
                }
                taxaDepthMap.put(taxon, mean);
            });
            BufferedWriter bw = IoUtils.getTextWriter(taxaDepthFileS);
            bw.write("Taxa\tCoverage");
            bw.newLine();
            for (int i = 0; i < this.taxaNames.length; i++) {
                bw.write(taxaNames[i]+"\t"+taxaDepthMap.get(taxaNames[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    public void writeTaxaBamMapFromList () {
        String infileS = "M:\\production\\maf\\annotations\\maf\\rawDataInfo\\source\\bamList.txt";
        String outfileS = "M:\\production\\maf\\annotations\\maf\\rawDataInfo\\taxaBamMap.txt";
        Table t = new Table (infileS);
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < t.getRowNumber(); i++) {
            t.content[i][0] = t.content[i][0].replaceFirst("./bam/c10/", "");
            nameSet.add(t.content[i][0].split("_")[0]);
        }
        String[] taxa = nameSet.toArray(new String[nameSet.size()]);
        Arrays.sort(taxa);
        ArrayList<String>[] bamList = new ArrayList[taxa.length];
        for (int i = 0; i < bamList.length; i++) bamList[i] = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            String name = t.content[i][0].split("_")[0];
            bamList[Arrays.binarySearch(taxa, name)].add(t.content[i][0]);
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            for (int i = 0; i < taxa.length; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(taxa[i]);
                for (int j = 0; j < bamList[i].size(); j++) {
                    sb.append("\t").append(bamList[i].get(j));
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
}
