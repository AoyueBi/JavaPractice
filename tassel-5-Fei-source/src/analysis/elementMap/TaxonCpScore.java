/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.elementMap;

import utils.BaseCoder;
import format.Fasta;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import net.openhft.koloboke.collect.map.hash.HashByteByteMap;
import net.openhft.koloboke.collect.map.hash.HashIntIntMap;
import net.openhft.koloboke.collect.map.hash.HashIntIntMapFactory;
import net.openhft.koloboke.collect.map.hash.HashIntIntMaps;
import net.openhft.koloboke.collect.map.hash.HashLongIntMap;
import net.openhft.koloboke.collect.map.hash.HashLongIntMapFactory;
import net.openhft.koloboke.collect.map.hash.HashLongIntMaps;
import net.openhft.koloboke.collect.map.hash.HashShortIntMap;
import net.openhft.koloboke.collect.map.hash.HashShortIntMapFactory;
import net.openhft.koloboke.collect.map.hash.HashShortIntMaps;
import utils.Benchmark;
import utils.FArrayUtils;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class TaxonCpScore {
    HashShortIntMap kmerCountShortMap = null;
    HashIntIntMap kmerCountIntMap = null;
    HashLongIntMap kmerCountLongMap = null;
    int kmerLength = 0;
    
    public TaxonCpScore (String referenceGenomeFileS, String kmerListDirS, String taxonName, List<String> fastqList, String taxaScoreDirS) {
        System.out.println("Start counting kmer in "+taxonName+" ...");
        this.initializeKmerCountMap(kmerListDirS);
        this.countKmer(fastqList);
        this.writeCopyNumberScore(referenceGenomeFileS, taxonName, taxaScoreDirS);
    }
    
    private void writeCopyNumberScore (String referenceGenomeFileS, String taxonName, String taxaScoreDirS) {
        new File(taxaScoreDirS).mkdir();
        File taxonDir = new File(taxaScoreDirS, taxonName);
        taxonDir.mkdir();
        int fragmentSize = 100000;
        HashByteByteMap ascIIByteMap = BaseCoder.getAscIIByteMap();
        Fasta f = new Fasta(referenceGenomeFileS);
        HashMap<Integer, String> chrSeqMap = new HashMap();
        for (int i = 0; i < f.getSeqNumber(); i++) {
            chrSeqMap.put(Integer.valueOf(f.getName(i)), f.getSeq(i));
        }
        Set<Map.Entry<Integer, String>> chrSeqset = chrSeqMap.entrySet();
        System.out.println("Writing copy number score of "+taxonName+" by chromosomes...");
        long start = System.nanoTime();
        chrSeqset.parallelStream().forEach(entry -> {
            int chr = entry.getKey();
            String seq = entry.getValue();
            byte[] bArray = seq.getBytes();
            for (int i = 0; i < bArray.length; i++) {
                bArray[i] = ascIIByteMap.get(bArray[i]);
            }
            StringBuilder sb  = new StringBuilder();
            sb.append(taxonName).append("_chr").append(FStringUtils.getNDigitNumber(3, chr)).append("_CpScore.txt.gz");
            String outfileS = new File(taxonDir, sb.toString()).getAbsolutePath();
            int[][] bound = FArrayUtils.getSubsetsIndicesBySubsetSize(seq.length(), fragmentSize);
            try {
                BufferedWriter bw = IoUtils.getTextGzipWriter(outfileS);
                bw.write("CpScore");
                bw.newLine();
                for (int i = 0; i < bound.length; i++) {
                    int intervalSize = bound[i][1] - bound[i][0];
                    TIntArrayList[] kmerCountList = new TIntArrayList[intervalSize];
                    for (int j = 0; j < kmerCountList.length; j++) kmerCountList[j] = new TIntArrayList();
                    int startIndex = bound[i][0]-kmerLength+1;
                    if (startIndex < 0) startIndex = 0;
                    int endIndex = bound[i][1];
                    if (endIndex-1+kmerLength > seq.length()) endIndex = seq.length()-kmerLength+1;
                    for (int j = startIndex; j < endIndex; j++) {
                        boolean flag = false;
                        for (int k = j; k < j+kmerLength; k++) {
                            if (bArray[k] >3) {
                                j = k;
                                flag = true;
                                break;
                            }
                        }
                        if (flag) continue;
                        int count = 0;
                        if (kmerCountShortMap != null) {
                            short query = BaseCoder.getShortSeqFromByteArray(Arrays.copyOfRange(bArray, j, j + kmerLength));
                            count = kmerCountShortMap.get(query);
                        }
                        else if (kmerCountIntMap != null) {
                            int query = BaseCoder.getIntSeqFromByteArray(Arrays.copyOfRange(bArray, j, j + kmerLength));
                            count = kmerCountIntMap.get(query);
                        }
                        else if (kmerCountLongMap != null) {
                            long query = BaseCoder.getLongSeqFromByteArray(Arrays.copyOfRange(bArray, j, j + kmerLength));
                            count = kmerCountLongMap.get(query);
                        }
                        int offSet = bound[i][0]-j;
                        if (offSet > 0) {
                            for (int k = j; k < j+kmerLength-offSet; k++) {
                                kmerCountList[offSet+k-bound[i][0]].add(count);
                            }
                        }
                        else {
                            int end = j+kmerLength;
                            if (end > bound[i][1]) end = bound[i][1];
                            for (int k = j; k < end; k++) {
                                kmerCountList[k-bound[i][0]].add(count);
                            }
                        }
                    }
                    for (int j = 0; j < kmerCountList.length; j++) {
                        int[] kmerCount = kmerCountList[j].toArray();
                        sb = new StringBuilder();
                        if (kmerCount.length == 0) {
                            sb.append("NA");
                        }
                        else {
                            float aver = 0;
                            for (int k = 0; k < kmerCount.length; k++) {
                                aver+=(float)kmerCount[k]/kmerCount.length;
                            }
                            sb.append(aver);
                        }
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    if ((bound[i][1])%50000000 == 0) System.out.println(String.valueOf(bound[i][1])+" sites on chr "+String.valueOf(chr)+" are processed for "+taxonName);
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        StringBuilder sb = new StringBuilder();
        sb.append("Writing copy number score of ").append(taxonName).append(" is finished. Time span: ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(sb.toString());
    }
    
    private void countKmer (List<String> fastqList) {
        HashByteByteMap ascIIByteMap = BaseCoder.getAscIIByteMap();
        System.out.println("Start processing fastqs. File number: "+String.valueOf(fastqList.size()));
        long start = System.nanoTime();
        if (kmerCountShortMap != null) {
            for (int i = 0; i < fastqList.size(); i++) {
                System.out.println("#"+String.valueOf(i+1)+": "+fastqList.get(i));
                try {
                    BufferedReader br = IoUtils.getTextGzipReader(fastqList.get(i));
                    String temp = null;
                    int cnt = 0;
                    while ((temp = br.readLine()) != null) {
                        cnt++;
                        if (cnt%5000000 == 0) System.out.println(String.valueOf(cnt)+" reads processed from " + fastqList.get(i));
                        temp = br.readLine();
                        byte[] bArray = temp.getBytes();
                        int endIndex = bArray.length;
                        for (int j = 0; j < bArray.length; j++) {
                            bArray[j] = ascIIByteMap.get(bArray[j]);
                        }
                        int searchStart = 0;
                        for (int j = 0; j < endIndex - kmerLength+1; j++) {
                            boolean flag = false;
                            for (int k = searchStart; k < j+kmerLength; k++) {
                                if (bArray[k] >3) {
                                    j = k;
                                    searchStart = k+1;
                                    flag = true;
                                    break;
                                }
                            }
                            if (flag) continue;
                            else searchStart = j+kmerLength;
                            byte[] subArray = Arrays.copyOfRange(bArray, j, j + kmerLength);
                            short kmerF = BaseCoder.getShortSeqFromByteArray(subArray);
                            short kmerR = BaseCoder.getShortReverseComplement(kmerF, kmerLength);
                            int value = kmerCountShortMap.get(kmerF);
                            if (value != kmerCountShortMap.defaultValue() && value != Integer.MAX_VALUE) {
                                kmerCountShortMap.addValue(kmerF, 1);
                            }
                            value = kmerCountShortMap.get(kmerR);
                            if (value != kmerCountShortMap.defaultValue() && value != Integer.MAX_VALUE) {
                                kmerCountShortMap.addValue(kmerR, 1);
                            }
                        }
                        br.readLine();br.readLine();
                    }
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        else if (kmerCountIntMap != null) {
            for (int i = 0; i < fastqList.size(); i++) {
                System.out.println("#"+String.valueOf(i+1)+": "+fastqList.get(i));
                try {
                    BufferedReader br = IoUtils.getTextGzipReader(fastqList.get(i));
                    String temp = null;
                    int cnt = 0;
                    while ((temp = br.readLine()) != null) {
                        cnt++;
                        if (cnt%5000000 == 0) System.out.println(String.valueOf(cnt)+" reads processed");
                        temp = br.readLine();
                        byte[] bArray = temp.getBytes();
                        int endIndex = bArray.length;
                        for (int j = 0; j < bArray.length; j++) {
                            bArray[j] = ascIIByteMap.get(bArray[j]);
                        }
                        int searchStart = 0;
                        for (int j = 0; j < endIndex - kmerLength+1; j++) {
                            boolean flag = false;
                            for (int k = searchStart; k < j+kmerLength; k++) {
                                if (bArray[k] >3) {
                                    j = k;
                                    searchStart = k+1;
                                    flag = true;
                                    break;
                                }
                            }
                            if (flag) continue;
                            else searchStart = j+kmerLength;
                            byte[] subArray = Arrays.copyOfRange(bArray, j, j + kmerLength);
                            int kmerF = BaseCoder.getIntSeqFromByteArray(subArray);
                            int kmerR = BaseCoder.getIntReverseComplement(kmerF, kmerLength);
                            int value = kmerCountIntMap.get(kmerF);
                            if (value != kmerCountIntMap.defaultValue() && value != Integer.MAX_VALUE) {
                                kmerCountIntMap.addValue(kmerF, 1);
                            }
                            value = kmerCountIntMap.get(kmerR);
                            if (value != kmerCountIntMap.defaultValue() && value != Integer.MAX_VALUE) {
                                kmerCountIntMap.addValue(kmerR, 1);
                            }
                        }
                        br.readLine();br.readLine();
                    }
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        else if (kmerCountLongMap != null) {
            for (int i = 0; i < fastqList.size(); i++) {
                System.out.println("#"+String.valueOf(i+1)+": "+fastqList.get(i));
                try {
                    BufferedReader br = IoUtils.getTextGzipReader(fastqList.get(i));
                    String temp = null;
                    int cnt = 0;
                    while ((temp = br.readLine()) != null) {
                        cnt++;
                        if (cnt%5000000 == 0) System.out.println(String.valueOf(cnt)+" reads processed");
                        temp = br.readLine();
                        byte[] bArray = temp.getBytes();
                        int endIndex = bArray.length;
                        for (int j = 0; j < bArray.length; j++) {
                            bArray[j] = ascIIByteMap.get(bArray[j]);
                        }
                        int searchStart = 0;
                        for (int j = 0; j < endIndex - kmerLength+1; j++) {
                            boolean flag = false;
                            for (int k = searchStart; k < j+kmerLength; k++) {
                                if (bArray[k] >3) {
                                    j = k;
                                    searchStart = k+1;
                                    flag = true;
                                    break;
                                }
                            }
                            if (flag) continue;
                            else searchStart = j+kmerLength;
                            byte[] subArray = Arrays.copyOfRange(bArray, j, j + kmerLength);
                            long kmerF = BaseCoder.getLongSeqFromByteArray(subArray);
                            long kmerR = BaseCoder.getLongReverseComplement(kmerF, kmerLength);
                            int value = kmerCountLongMap.get(kmerF);
                            if (value != kmerCountLongMap.defaultValue() && value != Integer.MAX_VALUE) {
                                kmerCountLongMap.addValue(kmerF, 1);
                            }
                            value = kmerCountLongMap.get(kmerR);
                            if (value != kmerCountLongMap.defaultValue() && value != Integer.MAX_VALUE) {
                                kmerCountLongMap.addValue(kmerR, 1);
                            }
                        }
                        br.readLine();br.readLine();
                    }
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        StringBuilder sb = new StringBuilder();
        sb.append("Fastq processing is finished. Time span: ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(sb.toString());
    }
    
    private void initializeKmerCountMap (String kmerListDirS) {
        System.out.println("Initializing kmerCountMap ...");
        long start = System.nanoTime();
        String kmerInfoFileS = new File(kmerListDirS, "kmerInfo.txt").getAbsolutePath();
        String kmerBlockDirS = new File(kmerListDirS, "kmerBlock").getAbsolutePath();
        File[] fs = new File(kmerBlockDirS).listFiles();
        int kmerNumber = 0;
        try {
            BufferedReader br = IoUtils.getTextReader(kmerInfoFileS);
            br.readLine();
            kmerLength = Integer.valueOf(br.readLine().split("\t")[1]);
            kmerNumber = Integer.valueOf(br.readLine().split("\t")[1]);
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("KmerLength = "+String.valueOf(kmerLength)+" bp. Kmer number = "+String.valueOf(kmerNumber));
        if (kmerLength>0 && kmerLength<=8) {
            short[] keys = new short[kmerNumber];
            int[] values = new int[kmerNumber];
            int cnt = 0;
            for (int i = 0; i < fs.length; i++) {
                try {
                    DataInputStream dis = IoUtils.getBinaryReader(fs[i].getAbsolutePath());
                    int size = dis.readInt();
                    for (int j = 0; j < size; j++) {
                        keys[cnt] = dis.readShort();
                        cnt++;
                    }
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
            HashShortIntMapFactory fac = HashShortIntMaps.getDefaultFactory().withDefaultValue(-1);
            kmerCountShortMap = fac.newMutableMap(keys, values);
        }
        else if (kmerLength <= 16) {
            int[] keys = new int[kmerNumber];
            int[] values = new int[kmerNumber];
            int cnt = 0;
            for (int i = 0; i < fs.length; i++) {
                try {
                    DataInputStream dis = IoUtils.getBinaryReader(fs[i].getAbsolutePath());
                    int size = dis.readInt();
                    for (int j = 0; j < size; j++) {
                        keys[cnt] = dis.readInt();
                        cnt++;
                    } 
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
            HashIntIntMapFactory fac = HashIntIntMaps.getDefaultFactory().withDefaultValue(-1);
            kmerCountIntMap = fac.newMutableMap(keys, values);
        }
        else if (kmerLength <= 32) {
            long[] keys = new long[kmerNumber];
            int[] values = new int[kmerNumber];
            int cnt = 0;
            for (int i = 0; i < fs.length; i++) {
                try {
                    DataInputStream dis = IoUtils.getBinaryReader(fs[i].getAbsolutePath());
                    int size = dis.readInt();
                    for (int j = 0; j < size; j++) {
                        keys[cnt] = dis.readLong();
                        cnt++;
                    }
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
            HashLongIntMapFactory fac = HashLongIntMaps.getDefaultFactory().withDefaultValue(-1);
            kmerCountLongMap = fac.newMutableMap(keys, values);
        }
        StringBuilder sb = new StringBuilder();
        sb.append("Initializing kmerCountMap is finished. Time span: ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(sb.toString());
    }
}
