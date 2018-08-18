/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.cassava.wgs;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import cern.jet.math.Arithmetic;
import static cern.jet.math.Arithmetic.factorial;
import static cern.jet.math.Arithmetic.logFactorial;
import static cern.jet.math.Arithmetic.longFactorial;
import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;
import format.Fasta;
import format.Table;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;
import gnu.trove.set.hash.TByteHashSet;
import gnu.trove.set.hash.TIntHashSet;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.LongAdder;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.apache.commons.math3.util.Combinations;

import utils.FArrayUtils;
import utils.FStringUtils;
import utils.IoUtils;

/**
 * This slow because SNP calling is done for each region, each region needs to read in reference to do pileup in that region.
 * @author Fei Lu
 */
public class SNPDiscovery {
    int[] chroms = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22};
    int[] chromLength = {34959721, 32431396, 29412403, 28749345, 28438989, 27939960, 27069033, 34011518, 29417918, 26335333, 27369018, 31602189, 28119335, 24855542, 26192760, 28979753, 27388923, 25234269, 32050860, 31903080, 161453, 1325823};
    String[] taxaNames = null;
    HashMap<String, String[]> taxaBamMap = new HashMap();
    HashMap<String, String> bamPileupMap = new HashMap();
    HashMap<String, ArrayList<String>> taxaPileupMap = new HashMap();
    ConcurrentHashMap<Integer, Double> factorialMap = new ConcurrentHashMap();
    int maxFactorial = 150;
    Fasta genomeFa = null;
    //A, C, G, T
    byte[] possibleSNP = {65, 67, 71, 84};
    //D, I
    byte[] possibleIndel = {68, 73};
    //A, C, D, G, I, T
    byte[] possibleAllele = {65, 67, 68, 71, 73, 84};
    
    
    int regionSize = 100000;
    double individualDepthRatioThresh = 0.4;
    double individualThirdAlleleRatioThresh = 0.2;
    double segregationPValueThresh = 0.001;
    double sequencingErrorRate = 0.1;
    
    
    public SNPDiscovery () {
        this.creatFactorialMap();
        this.localTest(1);
    }
    
    public SNPDiscovery (String parameterFileS) {
        this.callSNP(parameterFileS);
    }
    
    public void callSNP (String parameterFileS) {
        ArrayList<String> pLineList = new ArrayList();
        try {
            BufferedReader br = IoUtils.getTextReader(parameterFileS);
            String temp = null;
            boolean ifOut = false;
            if (!(temp = br.readLine()).equals("Author: Fei Lu")) ifOut = true;
            if (!(temp = br.readLine()).equals("Email: fl262@cornell.edu")) ifOut = true;
            if (!(temp = br.readLine()).equals("Homepage: https://sites.google.com/site/feilu0000/")) ifOut = true;
            if (ifOut) {
                System.out.println("Please keep the authorship in the parameter file. Program stops.");
                System.exit(0);
            }
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) continue;
                if (temp.isEmpty()) continue;
                pLineList.add(temp);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        String referenceFileS = pLineList.get(0);
        String bamDirS = pLineList.get(1);
        String taxaBamMapFileS = pLineList.get(2);
        int currentChr = Integer.valueOf(pLineList.get(3));
        String vcfDirS = pLineList.get(7);
        genomeFa = new Fasta(referenceFileS);
        genomeFa.sortRecordByNameValue();
        chroms = new int[genomeFa.getSeqNumber()];
        chromLength = new int[genomeFa.getSeqNumber()];
        for (int i = 0; i < genomeFa.getSeqNumber(); i++) {
            chroms[i] = Integer.valueOf(genomeFa.getName(i));
            chromLength[i] = genomeFa.getSeqLength(i);
        }
        int[][] positions = null;
        if (pLineList.get(4).startsWith("-")) {
            positions = new int[1][2];
            positions[0][0] = Integer.valueOf(pLineList.get(5));
            positions[0][1] = Integer.valueOf(pLineList.get(6));
        }
        else {
            Table t = new Table (pLineList.get(4));
            positions = new int[t.getRowNumber()][2];
            for (int i = 0; i < t.getRowNumber(); i++) {
                for (int j = 0; j < 2; j++) positions[i][j] = Integer.valueOf(t.content[i][j]);
            }
        }
 
        String pileupDirS = new File(new File(vcfDirS).getParent(), "pileup").getAbsolutePath();
        String tempVcfDirS = new File(new File(vcfDirS).getParent(), "vcfTemp").getAbsolutePath();
        new File(pileupDirS).mkdir();
        new File(vcfDirS).mkdir();
        new File(tempVcfDirS).mkdir();
        this.getTaxaBamMap(taxaBamMapFileS);
        File[] bams = new File(bamDirS).listFiles();
        this.updateTaxaBamMap(bams);
        this.creatPileupMap(pileupDirS);
        this.creatFactorialMap();
        
        for (int i = 0; i < positions.length; i++) {
            int startPos = positions[i][0];
            int endPos = positions[i][1];
            if (startPos < 1) {
                startPos = 1;
            }
            int chrIndex = Arrays.binarySearch(chroms, currentChr);
            if (endPos < 0 || endPos > this.chromLength[chrIndex]) {
                endPos = this.chromLength[chrIndex];
            }
            int[][] regionBound = this.creatRegions(currentChr, regionSize, startPos, endPos);
            for (int j = 0; j < regionBound.length; j++) {
                String[] commands = this.creatPileupCommand(currentChr, regionBound[j][0], regionBound[j][1], referenceFileS, bamDirS, pileupDirS);
                this.performPileUp(commands);
                StringBuilder sb = new StringBuilder();
                this.outputVCF(currentChr, regionBound[j][0], regionBound[j][1], tempVcfDirS, genomeFa.getSeq(chrIndex));
                sb.append("Finished variant calling on chromosome ").append(currentChr).append(" from ").append(regionBound[j][0]).append(" to ").append(regionBound[j][1]);
                System.out.print(sb.toString());
                System.out.println(". Used Memory:" + String.valueOf((double)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024/1024)+ "Gb");
                try {
                    FileUtils.cleanDirectory(new File(pileupDirS));
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        this.mergeVCF(tempVcfDirS, vcfDirS, currentChr);
        try {
            FileUtils.deleteDirectory(new File(pileupDirS));
            FileUtils.deleteDirectory(new File(tempVcfDirS));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Variant calling completed");
    }
      
    private void mergeVCF (String tempVcfDirS, String vcfDirS, int currentChr) {
        System.out.println("\nMerging sub VCF");
        String header = this.getVCFHeader();
        int cnt = 0;
        try {
            String outfileS = "chr"+FStringUtils.getNDigitNumber(3, currentChr)+".VCF.txt";
            outfileS = new File(vcfDirS, outfileS).getAbsolutePath();
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write(header);
            bw.newLine();
            File[] fs = new File (tempVcfDirS).listFiles();
            HashMap<Integer, File> posFileMap = new HashMap();
            for (int i = 0; i < fs.length; i++) {
                int key = Integer.valueOf(fs[i].getName().split("_")[1]);
                posFileMap.put(key, fs[i]);
            }
            Set<Integer> posSet = posFileMap.keySet();
            Integer[] startPoses = posSet.toArray(new Integer[posSet.size()]);
            Arrays.sort(startPoses);
            for (int i = 0; i < startPoses.length; i++) {
                BufferedReader br = IoUtils.getTextReader(posFileMap.get(startPoses[i]).getAbsolutePath());
                String temp;
                while ((temp = br.readLine()) != null) {
                    bw.write(temp);
                    bw.newLine();
                    cnt++;
                }
                br.close();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println(String.valueOf(cnt)+" SNV called");
    }
    
    private void outputVCF (int currentChr, int startPos, int endPos, String tempVcfDirS, String chrSeq) {
        ArrayList<Integer> sitePosList = new ArrayList();
        HashMap<Integer, String[]> posBaseMap = new HashMap();
        HashMap<Integer, int[]> posDepthMap = new HashMap();
        int[] pos = new int[endPos-startPos+1];
        for (int i = 0; i < pos.length; i++) pos[i] = startPos+i;
        ConcurrentHashMap<String, int[]> pileupDepthMap = new ConcurrentHashMap();
        ConcurrentHashMap<String, String[]> pileupBaseMap = new ConcurrentHashMap();
        List<String> taxaList = Arrays.asList(taxaNames);
        taxaList.parallelStream().forEach(taxon -> {
            ArrayList<String> pileupList = this.taxaPileupMap.get(taxon);
            pileupList.parallelStream().forEach(pileupFileS -> {
               String temp = null;
               try {
                    BufferedReader br = IoUtils.getTextReader(pileupFileS);
                    int[] depth = new int[pos.length];
                    String[] base = new String[pos.length];
                    while ((temp = br.readLine()) != null) {
                        String[] tem = temp.split("\t");
                        int position = Integer.valueOf(tem[1]);
                        int index = Arrays.binarySearch(pos, position);
                        depth[index] = Integer.valueOf(tem[3]);
                        if (depth[index] != 0) base[index] = tem[4];
                    }
                    pileupDepthMap.put(pileupFileS, depth);
                    pileupBaseMap.put(pileupFileS, base);
               }
               catch (Exception e) {
                   e.printStackTrace();
                   //System.out.println(temp);
               }
            });
        });
        HashMap<String, Integer> taxaIndexMap = new HashMap ();
        for (int i = 0; i < this.taxaNames.length; i++) {
            taxaIndexMap.put(taxaNames[i], i);
        }
        AtomicInteger atomicIndex = new AtomicInteger();
        for (int i = 0; i < pos.length; i++) {
            String[]  taxaSiteBase = new String[this.taxaNames.length];
            int[] taxaSiteDepth = new int[this.taxaNames.length];
            atomicIndex.set(i);
            taxaList.parallelStream().forEach(taxon -> {
                ArrayList<String> pileupFileList = taxaPileupMap.get(taxon);
                StringBuilder sb = new StringBuilder();
                int d = 0;
                int currentIndex = atomicIndex.get();
                for (int j = 0; j < pileupFileList.size(); j++) {                  
                    int td = pileupDepthMap.get(pileupFileList.get(j))[currentIndex];
                    d+=td;
                    if (td == 0) continue;
                    sb.append(pileupBaseMap.get(pileupFileList.get(j))[currentIndex]);
                }
                taxaSiteBase[taxaIndexMap.get(taxon)] = sb.toString();
                taxaSiteDepth[taxaIndexMap.get(taxon)] = d;
            });
            for (int j = 0; j < taxaSiteDepth.length; j++) {
                if (taxaSiteDepth[j] != 0) {
                    sitePosList.add(pos[i]);
                    posBaseMap.put(pos[i], taxaSiteBase);
                    posDepthMap.put(pos[i], taxaSiteDepth);
                    break;
                }
            }
        }
        ConcurrentHashMap<Integer,String> posVCFMap = new ConcurrentHashMap();
        sitePosList.parallelStream().forEach(sitePos -> {
            String vcfRecord = this.getVCFRecord(posBaseMap.get(sitePos), posDepthMap.get(sitePos), currentChr, sitePos, chrSeq);
            if (vcfRecord != null) {
                posVCFMap.put(sitePos, vcfRecord);
            }
        });
        String outfileS = this.buildFileName(currentChr, startPos, endPos, tempVcfDirS, ".tVCF.txt");
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            Set<Integer> sitePosSet = posVCFMap.keySet();
            Integer[] sitePos = sitePosSet.toArray(new Integer[sitePosSet.size()]);
            Arrays.sort(sitePos);
            for (int i = 0; i < sitePos.length; i++) {
                bw.write(posVCFMap.get(sitePos[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private String buildFileName (int currentChr, int startPos, int endPos, String dirS, String suffix) {
        StringBuilder sb = new StringBuilder();
        sb.append("chr").append(FStringUtils.getNDigitNumber(5, currentChr)).append("_").append(startPos).append("_").append(endPos).append(suffix);
        return new File (dirS, sb.toString()).getAbsolutePath();
    }
    
    private String getVCFHeader () {
        StringBuilder sb = new StringBuilder();
        sb.append("#CHR	POS	REF	ALT	INFO	FORMAT");
        for (int i = 0; i < taxaNames.length; i++) {
            sb.append("\t").append(taxaNames[i]);
        }
        return sb.toString();
    }
    
    private String getVCFRecord (String[] base, int[] depth, int currentChr, int sitePos, String chrSeq) {
        String refBase = String.valueOf(chrSeq.charAt(sitePos-1));
        if (refBase.equals("N")) return null;
        StringBuilder sb = new StringBuilder();
        this.processBase(base);
        Pattern pIndel = Pattern.compile("\\+\\w+|\\-\\w+");
	Matcher m;
        byte[] baseByte;
        int[][] pSNPDepth = new int[this.possibleSNP.length][base.length];
        int[] refDepth = new int[base.length];
        TIntHashSet indelIndexSet = new TIntHashSet(); 
        HashSet<String> currentIndelSet = new HashSet();
        HashSet<String> allIndelSet = new HashSet();
        for (int i = 0; i < base.length; i++) {
            m=pIndel.matcher(base[i]);
            boolean isThereIndel = false;
            while (m.find()) {
                String s = this.getIndelString(base[i].substring(m.start(), m.end()));
                currentIndelSet.add(s);
                isThereIndel = true;
                indelIndexSet.add(i);
            }
            String currentBase = new String(base[i]);
            if (isThereIndel) {
                String[] currentIndel = currentIndelSet.toArray(new String[currentIndelSet.size()]);
                for (int j = 0; j < currentIndel.length; j++) {
                    currentBase = currentBase.replaceAll("\\"+currentIndel[j], "");
                    allIndelSet.add(currentIndel[j]);
                }
            }
            baseByte = currentBase.getBytes();
            for (int j = 0; j < baseByte.length; j++) {
                if (baseByte[j] == 46) {
                    refDepth[i]++;
                    continue;
                }
                int index = Arrays.binarySearch(this.possibleSNP, baseByte[j]);
                pSNPDepth[index][i]++;
            }
        }
        String[] allIndel = allIndelSet.toArray(new String[allIndelSet.size()]);
        int[][] pIndelDepth = null;
        if (allIndel.length != 0) {
            Arrays.sort(allIndel);
            int[] indelTaxaIndex = indelIndexSet.toArray();
            Arrays.sort(indelTaxaIndex);
            pIndelDepth = new int[allIndel.length][indelTaxaIndex.length];
            for (int i = 0; i < indelTaxaIndex.length; i++) {
                m=pIndel.matcher(base[indelTaxaIndex[i]]);
                while (m.find()) {
                    String query = this.getIndelString(base[indelTaxaIndex[i]].substring(m.start(), m.end()));
                    int index = Arrays.binarySearch(allIndel, query);
                    pIndelDepth[index][i]++;
                }
            }
            int[][] updatedIndelDepth = new int[this.possibleIndel.length][base.length];
            for (int i = 0; i < indelTaxaIndex.length; i++) {
                for (int j = 0; j < allIndel.length; j++) {
                    if (allIndel[j].startsWith("-")) {
                        updatedIndelDepth[0][indelTaxaIndex[i]]+=pIndelDepth[j][i];
                    }
                    else {
                        updatedIndelDepth[1][indelTaxaIndex[i]]+=pIndelDepth[j][i];
                    }
                }
            }
            pIndelDepth = updatedIndelDepth;
        }
        
        int[][] pAlleleDepth = new int[this.possibleAllele.length][];
        for (int i = 0; i < this.possibleSNP.length; i++) {
            int index = Arrays.binarySearch(this.possibleAllele, this.possibleSNP[i]);
            pAlleleDepth[index] = pSNPDepth[i];
        }
        if (allIndel.length == 0) {
            for (int i = 0; i < this.possibleIndel.length; i++) {
                int index = Arrays.binarySearch(this.possibleAllele, this.possibleIndel[i]);
                pAlleleDepth[index] = new int[base.length];
            }
        }
        else {
             for (int i = 0; i < this.possibleIndel.length; i++) {
                int index = Arrays.binarySearch(this.possibleAllele, this.possibleIndel[i]);
                pAlleleDepth[index] = pIndelDepth[i];
            }
        }
        //****************************Filter1 Depth_ratio_test*****************************************
        //In any individual, alt allele show up < 2 times, ignore
        //In any individual, alt allele show up 2 times when 2 < depth < 5. Pick up
        //In any individual, alt allele show up > individualDepthRatioThresh, when depth >= 5. Pick up
        //When depth is low, tend to have assembly errors, LTR
        boolean[][] ifAllele = new boolean[this.possibleAllele.length][base.length];
        TByteHashSet alleleSet = new TByteHashSet();
        for (int i = 0; i < ifAllele.length; i++) {
            for (int j = 0; j < ifAllele[i].length; j++) {
                if (depth[j] < 2) {}
                else if (depth[j] < 5) {
                    if (pAlleleDepth[i][j] >=2) {
                        ifAllele[i][j] = true;
                        alleleSet.add(this.possibleAllele[i]);
                    }
                }
                else {
                    if ((double)pAlleleDepth[i][j]/depth[j] > this.individualDepthRatioThresh) {
                        ifAllele[i][j] = true;
                        alleleSet.add(this.possibleAllele[i]);
                    }
                }
                
            }
        }
        byte[] allele = alleleSet.toArray();
        if (allele.length == 0) return null;
        //=======================================Filter1=====================================================       
        
        Arrays.sort(allele);
        int[] indelTypeCount = new int[2];
        if (Arrays.binarySearch(allele, (byte)68) >= 0) {
            HashSet<String> indelTypeSet = new HashSet();
            for (int i = 0; i < allIndel.length; i++) {
                if (allIndel[i].startsWith("-")) indelTypeSet.add(allIndel[i]);
            }
            indelTypeCount[0] = indelTypeSet.size();
        }
        if (Arrays.binarySearch(allele, (byte)73) >= 0) {
            HashSet<String> indelTypeSet = new HashSet();
            for (int i = 0; i < allIndel.length; i++) {
                if (allIndel[i].startsWith("+")) indelTypeSet.add(allIndel[i]);
            }
            indelTypeCount[1] = indelTypeSet.size();
        }
        int[] allele2PAlleleIndex = new int[allele.length];
        for (int i = 0; i < this.possibleAllele.length; i++) {
            int index = Arrays.binarySearch(allele, this.possibleAllele[i]);
            if (index < 0) continue;
            allele2PAlleleIndex[index] = i;
        }       
        int[] alleleTotalDepth = new int[allele.length];
        for (int i = 0; i < allele.length; i++) {
            for (int j = 0; j < base.length; j++) {
                alleleTotalDepth[i]+=pAlleleDepth[allele2PAlleleIndex[i]][j];
            }
        }
        int[] alleleDepthDesendingIndex = FArrayUtils.getIndexByDescendingValue(alleleTotalDepth);
        int refTotalDepth = 0;
        for (int i = 0; i < refDepth.length; i++) refTotalDepth+=refDepth[i];
        int allDepth = refTotalDepth;
        for (int i = 0; i < alleleTotalDepth.length; i++) allDepth+=alleleTotalDepth[i];
        
        
        //****************************Filter2 third_allele_test************************************************
        //individual should not have the third allele
        if (allele.length > 1) {
            for (int i = 0; i < base.length; i++) {
                int[] tempCnt = new int[allele.length];
                for (int j = 0; j < this.possibleAllele.length; j++) {
                    int index = Arrays.binarySearch(allele, this.possibleAllele[j]);
                    if (index < 0) continue;
                    tempCnt[index]+=pAlleleDepth[j][i];
                }
                int sum = refDepth[i];
                for (int j = 0; j < tempCnt.length; j++) sum+=tempCnt[j];
                double[] v = new double[allele.length+1];
                for (int j = 0; j < tempCnt.length; j++) v[j] = (double)tempCnt[j]/sum;
                v[v.length-1] = (double)refDepth[i]/sum;
                Arrays.sort(v);
                if (v[v.length-3] > individualThirdAlleleRatioThresh) return null;
            }
        }
        //===========================Filter2=========================================================
        
        //****************************Filter3 Segregation_test*****************************************
        long[] observed = new long[base.length];
        double[] expected = new double[base.length];
        double[] segregationP = new double[allele.length];
        ChiSquareTest ct = new ChiSquareTest();
        int cnt = 0;
        for (int i = 0; i < allele.length; i++) {
            double r = (double)alleleTotalDepth[i]/(refTotalDepth+alleleTotalDepth[i]);
            for (int j = 0; j < base.length; j++) {
                observed[j] = pAlleleDepth[allele2PAlleleIndex[i]][j];
                expected[j] =  r;
            }
            segregationP[i] = ct.chiSquareTest(expected, observed);
            if (segregationP[i] > this.segregationPValueThresh) cnt++;
        }
        if (cnt == allele.length) return null;
       //===========================Filter3=========================================================
        
        int nonMissingCnt = 0;
        int[] refAndAllelePresence = new int[allele.length+1];
        for (int i = 0; i < base.length; i++) {
            if (refDepth[i] != 0) {
                nonMissingCnt++;
                refAndAllelePresence[0]++;
            }
            else {
                for (int j = 0; j < allele.length; j++) {
                    if (pAlleleDepth[allele2PAlleleIndex[j]][i] != 0) {
                        nonMissingCnt++;
                        break;
                    }
                }
            }
            for (int j = 0; j < allele.length; j++) {
                if (pAlleleDepth[allele2PAlleleIndex[j]][i] != 0) {
                   refAndAllelePresence[j+1]++;
                }
            }
        }
        
        sb.append(currentChr).append("\t").append(sitePos).append("\t").append(refBase).append("\t");
        for (int i = 0; i < allele.length; i++) sb.append(String.valueOf((char)allele[alleleDepthDesendingIndex[i]])).append(",");
        sb.deleteCharAt(sb.length()-1);
        sb.append("\t").append("DP=").append(allDepth).append(";AD=").append(refTotalDepth);
        for (int i = 0; i < alleleTotalDepth.length; i++) sb.append(",").append(alleleTotalDepth[alleleDepthDesendingIndex[i]]);
        sb.append(";NZ=").append(nonMissingCnt).append(";AP=");
        for (int i = 0; i < refAndAllelePresence.length; i++) sb.append(refAndAllelePresence[i]).append(",");
        sb.deleteCharAt(sb.length()-1);
        sb.append(";PV=");
        for (int i = 0; i < alleleTotalDepth.length; i++) sb.append(segregationP[i]).append(",");
        sb.deleteCharAt(sb.length()-1);
        sb.append(";DI=");
        for (int i = 0; i < indelTypeCount.length; i++) sb.append(indelTypeCount[i]).append(",");
        sb.deleteCharAt(sb.length()-1);
        sb.append("\t").append("GT:AD:PL");
        
        for (int i = 0; i < base.length; i++) {
            int[] dep = new int[allele.length+1];
            dep[0] = refDepth[i];
            for (int j = 0; j < allele.length; j++) {
                dep[j+1] = pAlleleDepth[allele2PAlleleIndex[alleleDepthDesendingIndex[j]]][i];
            }
            sb.append("\t").append(this.getGenotype(dep));
        }
        return sb.toString();
    }
    
    private void creatFactorialMap () {
        for (int i = 0; i < this.maxFactorial+1; i++) {
            this.factorialMap.put(i, factorial(i));
        }
    }
    
    public String getGenotype (int[] cnt) {
        int n = cnt.length*(cnt.length+1)/2;
        int[] likelihood = new int[n];
        int sum = 0;
        for (int i = 0; i < cnt.length; i++) sum+=cnt[i];
        if (sum == 0) return "./.";
        else if (sum > this.maxFactorial) {
            double portion = (double)this.maxFactorial/sum;
            for (int i = 0; i < cnt.length; i++) {
                cnt[i] = (int)(cnt[i]*portion);
            }
            sum = this.maxFactorial;
        }
        double coe = this.factorialMap.get(sum);
        for (int i = 0; i < cnt.length; i++) coe = coe/this.factorialMap.get(cnt[i]);
        double max = Double.MAX_VALUE;
        int a1 = 0;
        int a2 = 0;
        for (int i = 0; i < cnt.length; i++) {
            for (int j = i; j < cnt.length; j++) {
                int index = (j*(j+1)/2)+i;
                double value = Double.MAX_VALUE;
                if (i == j) {
                    value = -Math.log10(coe*Math.pow((1-0.75*this.sequencingErrorRate), cnt[i])*Math.pow(this.sequencingErrorRate/4, (sum-cnt[i])));
                }
                else {
                    value = -Math.log10(coe*Math.pow((0.5-this.sequencingErrorRate/4), cnt[i]+cnt[j])*Math.pow(this.sequencingErrorRate/4, (sum-cnt[i]-cnt[j])));
                }
                if (value < max) {
                    max = value;
                    a1 = i;
                    a2 = j;
                }
                likelihood[index] = (int)Math.round(value);
            }
        }
        StringBuilder sb = new StringBuilder();
        sb.append(a1).append("/").append(a2).append(":");
        for (int i = 0; i < cnt.length; i++) sb.append(cnt[i]).append(",");
        sb.deleteCharAt(sb.length()-1); sb.append(":");
        for (int i = 0; i < likelihood.length; i++) sb.append(likelihood[i]).append(",");
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }
    
    public int[] getGTLikelihood (int[] cnt) {
        int n = cnt.length*(cnt.length+1)/2;
        int[] likelihood = new int[n];
        int sum = 0;
        for (int i = 0; i < cnt.length; i++) sum+=cnt[i];
        double coe = factorial(sum);
        for (int i = 0; i < cnt.length; i++) coe = coe/factorial(cnt[i]);
        for (int i = 0; i < cnt.length; i++) {
            for (int j = i; j < cnt.length; j++) {
                int index = (j*(j+1)/2)+i;
                if (i == j) {
                    likelihood[index] = (int)Math.round(-Math.log10(coe*Math.pow((1-0.75*this.sequencingErrorRate), cnt[i])*Math.pow(this.sequencingErrorRate/4, (sum-cnt[i]))));
                }
                else {
                    likelihood[index] = (int)Math.round(-Math.log10(coe*Math.pow((0.5-this.sequencingErrorRate/4), cnt[i]+cnt[j])*Math.pow(this.sequencingErrorRate/4, (sum-cnt[i]-cnt[j]))));
                }
            }
        }
        return likelihood;
    }
    
    private void processBase (String[] base) {
        for (int i = 0; i < base.length; i++) {
            base[i] = base[i].replaceAll("\\^.|\\$|>|<", "");
            base[i] = base[i].replaceAll("a", "A");
            base[i] = base[i].replaceAll("g", "G");
            base[i] = base[i].replaceAll("t", "T");
            base[i] = base[i].replaceAll("c", "C");
            base[i] = base[i].replaceAll("n", "N");
            base[i] = base[i].replaceAll(",", ".");
            base[i] = base[i].replaceAll("\\*", ".");
        }
    }
    
    private String getIndelString (String raw) {
        byte[] b = raw.getBytes();
        int endIndex = 1;
        for (int i = 1; i < b.length; i++) {
            if (b[i] > 57 || b[i] < 48) {
                endIndex = i;
                break;
            }
        }
        int n = Integer.valueOf(raw.substring(1, endIndex));
        return raw.substring(0, endIndex+n);
    }
    
    private void creatPileupMap (String pileupDirS) {
        Set<Entry<String, String[]>> entries = taxaBamMap.entrySet();
        for (Entry<String, String[]> e : entries) {
            String[] bams = e.getValue();
            ArrayList<String> pileupList = new ArrayList();
            for (String bam : bams) {
                String pileupFileS = new File (pileupDirS, bam.replaceFirst(".bam", ".pileup.txt")).getAbsolutePath();
                bamPileupMap.put(bam, pileupFileS);
                pileupList.add(pileupFileS);
            }
            taxaPileupMap.put(e.getKey(), pileupList);
        }
    }
    
    
    private void performPileUp (String[] commands) {
        int coreNumber = Runtime.getRuntime().availableProcessors();
        int[][] batchIndex = FArrayUtils.getSubsetsIndicesBySubsetSize(commands.length, coreNumber);
        List<String> commandList = Arrays.asList(commands);
        for (int i = 0; i < batchIndex.length; i++) {
            List<String> subList = commandList.subList(batchIndex[i][0], batchIndex[i][1]);
            subList.parallelStream().forEach(element -> {
                try {
                    Runtime run = Runtime.getRuntime();
                    Process p = run.exec(element);
                    p.waitFor();
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            });
        }
    }
    
    /**
     * 
     * @param currentChr
     * @param startPos
     * @param endPos is inclusive
     * @param referenceFileS
     * @param bamDirS
     * @param pileupDirS
     * @return 
     */
    private String[] creatPileupCommand (int currentChr, int startPos, int endPos, String referenceFileS, String bamDirS, String pileupDirS) {
        ArrayList<String> commandList = new ArrayList();
        for (int i = 0; i < taxaNames.length; i++) {
            String[] bams = taxaBamMap.get(taxaNames[i]);
            for (int j = 0; j < bams.length; j++) {
                StringBuilder sb = new StringBuilder();
                String bamFileS = new File(bamDirS, bams[j]).getAbsolutePath();        
                String pileupFileS = this.bamPileupMap.get(bams[j]);
                sb.append("samtools mpileup -A -B -q 30 -Q 10 -f ").append(referenceFileS).append(" ").append(bamFileS).append(" -r ");
                sb.append(currentChr).append(":").append(startPos).append("-").append(endPos).append(" -o ").append(pileupFileS);
                commandList.add(sb.toString());
            }
        }
        return commandList.toArray(new String[commandList.size()]);
    }
    
    private int[][] creatRegions (int currentChr, int regionSize, int startPos, int endPos) {
        int chrIndex = Arrays.binarySearch(chroms, currentChr);
        int[][] regionBound = FArrayUtils.getSubsetsIndicesBySubsetSize(endPos-startPos+1, regionSize);
        for (int i = 0; i < regionBound.length; i++) {
            regionBound[i][0]+=startPos;
            regionBound[i][1]+=startPos;
            regionBound[i][1]--;
        }
        System.out.println("\nGenotyping will performed on chromosome "+String.valueOf(currentChr)+" from "+String.valueOf(startPos)+" to "+String.valueOf(endPos));
        System.out.println("Chromosome " + String.valueOf(currentChr) +"  is devided into " + String.valueOf(regionBound.length) + " regions\n");
        return regionBound;
    }
    
    private void updateTaxaBamMap (File[] bams) {
        String[] existingBam = new String[bams.length];
        for (int i = 0; i < bams.length; i++) existingBam[i] = bams[i].getName();
        Arrays.sort(existingBam);
        HashSet<String> existingTaxaSet = new HashSet();
        HashMap<String, String[]> updatedTaxaBamMap = new HashMap();
        int cnt = 0;
        for (int i = 0; i < taxaNames.length; i++) {
            String[] bamNames = taxaBamMap.get(taxaNames[i]);
            ArrayList<String> bamList = new ArrayList();
            for (int j = 0; j < bamNames.length; j++) {
                int index = Arrays.binarySearch(existingBam, bamNames[j]);
                if (index < 0) continue;
                bamList.add(bamNames[j]);
                existingTaxaSet.add(taxaNames[i]);
            }
            if (bamList.isEmpty()) continue;
            bamNames = bamList.toArray(new String[bamList.size()]);
            Arrays.sort(bamNames);
            updatedTaxaBamMap.put(taxaNames[i], bamNames);
            cnt+=bamNames.length;
        }
        String[] updatedTaxaNames = existingTaxaSet.toArray(new String[existingTaxaSet.size()]);
        Arrays.sort(updatedTaxaNames);
        taxaNames = updatedTaxaNames;
        taxaBamMap = updatedTaxaBamMap;
        System.out.println("Actual taxa number:\t"+String.valueOf(taxaNames.length));
        System.out.println("Actual bam file number:\t"+String.valueOf(cnt));
        System.out.println();
    }
    
    private void getTaxaBamMap (String taxaBamMapFileS) {
        try {
            BufferedReader br = IoUtils.getTextReader(taxaBamMapFileS);
            String temp;
            ArrayList<String> taxaList = new ArrayList();
            int nBam = 0;
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                taxaList.add(tem[0]);
                String[] bams = new String[tem.length-1] ;
                for (int i = 0; i < bams.length; i++) bams[i] = tem[i+1];
                Arrays.sort(bams);
                taxaBamMap.put(tem[0], bams);
                nBam+=bams.length;
            }
            taxaNames = taxaList.toArray(new String[taxaList.size()]);
            Arrays.sort(taxaNames);
            System.out.println("Created TaxaBamMap from" + taxaBamMapFileS);
            System.out.println("Taxa number:\t"+String.valueOf(taxaNames.length));
            System.out.println("Bam file number in TaxaBamMap:\t"+String.valueOf(nBam));
            System.out.println();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public static void main (String[] args) {
        new SNPDiscovery (args[0]);
    }
    
    //********************************Local test section***********************************   
    public void localTest (int currentChr) {
        String taxaBamMapFileS = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\GenotypingPipelineTest\\source\\taxaBamMap.txt";
        String pileupDirS = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\GenotypingPipelineTest\\pileup\\";
        String vcfDirS = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\GenotypingPipelineTest\\vcf\\";
        String bamList = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\GenotypingPipelineTest\\bam.list.txt" ;
        String referenceFileS = "M:\\Database\\cassavaReference\\genome\\Manihot esculenta\\cassavaV6_chr01.fa";
        Fasta f = new Fasta(referenceFileS);
        String chrSeq = f.getSeq(0);
        this.getTaxaBamMap(taxaBamMapFileS);
        this.updateTaxaBamMap(bamList);
        this.creatPileupMap(pileupDirS);
        this.creatFactorialMap();
        int[][] regionBound = this.creatRegions(currentChr, regionSize, 1, 100000);
        for (int i = 0; i < regionBound.length; i++) {
            this.outputVCF(currentChr, regionBound[i][0], regionBound[i][1], vcfDirS, chrSeq);
            break;
        }
    }
    
    private void updateTaxaBamMap (String bamListFileS) {
        Table t = new Table (bamListFileS);
        String[] existingBam = new String[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) existingBam[i] = t.content[i][0];
        Arrays.sort(existingBam);
        HashSet<String> existingTaxaSet = new HashSet();
        HashMap<String, String[]> updatedTaxaBamMap = new HashMap();
        int cnt = 0;
        for (int i = 0; i < taxaNames.length; i++) {
            String[] bamNames = taxaBamMap.get(taxaNames[i]);
            ArrayList<String> bamList = new ArrayList();
            for (int j = 0; j < bamNames.length; j++) {
                int index = Arrays.binarySearch(existingBam, bamNames[j]);
                if (index < 0) continue;
                bamList.add(bamNames[j]);
                existingTaxaSet.add(taxaNames[i]);
            }
            if (bamList.isEmpty()) continue;
            bamNames = bamList.toArray(new String[bamList.size()]);
            Arrays.sort(bamNames);
            updatedTaxaBamMap.put(taxaNames[i], bamNames);
            cnt+=bamNames.length;
        }
        String[] updatedTaxaNames = existingTaxaSet.toArray(new String[existingTaxaSet.size()]);
        Arrays.sort(updatedTaxaNames);
        taxaNames = updatedTaxaNames;
        taxaBamMap = updatedTaxaBamMap;
        System.out.println("Actual taxa number:\t"+String.valueOf(taxaNames.length));
        System.out.println("Actual bam file number:\t"+String.valueOf(cnt));
        System.out.println();
    }  
//************************************************************************************************************  
}
