/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.elementMap;

import utils.BaseCoder;
import format.Fasta;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import net.openhft.koloboke.collect.map.hash.HashByteByteMap;
import net.openhft.koloboke.collect.set.hash.HashIntSet;
import net.openhft.koloboke.collect.set.hash.HashIntSets;
import net.openhft.koloboke.collect.set.hash.HashLongSet;
import net.openhft.koloboke.collect.set.hash.HashLongSets;
import net.openhft.koloboke.collect.set.hash.HashShortSet;
import net.openhft.koloboke.collect.set.hash.HashShortSets;
import utils.Benchmark;
import utils.FArrayUtils;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class ReferenceKmer {
    static int shortKmerLength = 8;
    static int intKmerLength = 16;
    static int longKmerLength = 32;
    
    public ReferenceKmer (int kmerLength, String referenceGenomeFileS, String kmerListDirS) {
        new File(kmerListDirS).mkdir();
        this.buildReferenceKmers(kmerLength, referenceGenomeFileS, kmerListDirS);
    }
    
    private void buildReferenceKmers (int kmerLength, String referenceGenomeFileS, String kmerListDirS) {
        if (kmerLength > 0 && kmerLength <= shortKmerLength) {
            short[] kmers = this.getShortKmerList(kmerLength, referenceGenomeFileS);
            this.writeReferenceShortKmers(kmerListDirS, kmers, kmerLength, referenceGenomeFileS);
        }
        else if (kmerLength <= intKmerLength){
            int[] kmers = this.getIntKmerList(kmerLength, referenceGenomeFileS);
            this.writeReferenceIntKmers(kmerListDirS, kmers, kmerLength, referenceGenomeFileS);
        }
        else if (kmerLength <= longKmerLength){
            long[] kmers = this.getLongKmerList(kmerLength, referenceGenomeFileS);
            this.writeReferenceLongKmers(kmerListDirS, kmers, kmerLength, referenceGenomeFileS);
        }
        else {
            
        }
    }
    
    private void writeReferenceShortKmers (String kmerListDirS, short[] kmers, int kmerLength, String referenceGenomeFileS) {
        int kmerFileNumber = 128;
        String kmerInfoFileS = new File (kmerListDirS, "kmerInfo.txt").getAbsolutePath();
        File kmerBlockDir = new File (kmerListDirS, "kmerBlock");
        kmerBlockDir.mkdir();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(kmerInfoFileS);
            bw.write("Reference genome:\t"+referenceGenomeFileS); bw.newLine();
            bw.write("KmerLength:\t"+String.valueOf(kmerLength)); bw.newLine();
            bw.write("KmerNumber:\t"+String.valueOf(kmers.length)); bw.newLine();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        if (kmerFileNumber > kmers.length) kmerFileNumber = kmers.length;
        int[][] bound = FArrayUtils.getSubsetsIndicesBySubsetNumber(kmers.length, kmerFileNumber);
        ArrayList<Integer> fileIndexList = new ArrayList();
        for (int i = 0; i < bound.length; i++) fileIndexList.add(i);
        long start = System.nanoTime();
        System.out.println("Writing kmer list to hard drive...");
        fileIndexList.parallelStream().forEach(fileIndex -> {
            String outfileS = new File (kmerBlockDir, "kmer"+FStringUtils.getNDigitNumber(4, fileIndex)+".bin").getAbsolutePath();
            try {
                DataOutputStream dos = IoUtils.getBinaryWriter(outfileS);
                dos.writeInt(bound[fileIndex][1] - bound[fileIndex][0]);
                for (int i = bound[fileIndex][0]; i < bound[fileIndex][1]; i++) {
                    dos.writeShort(kmers[i]);
                }
                dos.flush();
                dos.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        System.out.println("Finished writing kmers. Time span: " + String.valueOf(Benchmark.getTimeSpanSeconds(start)) + " seconds. Memory used: "+String.valueOf(Benchmark.getUsedMemoryGb())+" Gb");
    } 
    
    private void writeReferenceIntKmers (String kmerListDirS, int[] kmers, int kmerLength, String referenceGenomeFileS) {
        int kmerFileNumber = 128;
        String kmerInfoFileS = new File (kmerListDirS, "kmerInfo.txt").getAbsolutePath();
        File kmerBlockDir = new File (kmerListDirS, "kmerBlock");
        kmerBlockDir.mkdir();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(kmerInfoFileS);
            bw.write("Reference genome:\t"+referenceGenomeFileS); bw.newLine();
            bw.write("KmerLength:\t"+String.valueOf(kmerLength)); bw.newLine();
            bw.write("KmerNumber:\t"+String.valueOf(kmers.length)); bw.newLine();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        if (kmerFileNumber > kmers.length) kmerFileNumber = kmers.length;
        int[][] bound = FArrayUtils.getSubsetsIndicesBySubsetNumber(kmers.length, kmerFileNumber);
        ArrayList<Integer> fileIndexList = new ArrayList();
        for (int i = 0; i < bound.length; i++) fileIndexList.add(i);
        long start = System.nanoTime();
        System.out.println("Writing kmer list to hard drive...");
        fileIndexList.parallelStream().forEach(fileIndex -> {
            String outfileS = new File (kmerBlockDir, "kmer"+FStringUtils.getNDigitNumber(4, fileIndex)+".bin").getAbsolutePath();
            try {
                DataOutputStream dos = IoUtils.getBinaryWriter(outfileS);
                dos.writeInt(bound[fileIndex][1] - bound[fileIndex][0]);
                for (int i = bound[fileIndex][0]; i < bound[fileIndex][1]; i++) {
                    dos.writeInt(kmers[i]);
                }
                dos.flush();
                dos.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        System.out.println("Finished writing kmers. Time span: " + String.valueOf(Benchmark.getTimeSpanSeconds(start)) + " seconds. Memory used: "+String.valueOf(Benchmark.getUsedMemoryGb())+" Gb");
    } 
     
    private void writeReferenceLongKmers (String kmerListDirS, long[] kmers, int kmerLength, String referenceGenomeFileS) {
        int kmerFileNumber = 128;
        String kmerInfoFileS = new File (kmerListDirS, "kmerInfo.txt").getAbsolutePath();
        File kmerBlockDir = new File (kmerListDirS, "kmerBlock");
        kmerBlockDir.mkdir();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(kmerInfoFileS);
            bw.write("Reference genome:\t"+referenceGenomeFileS); bw.newLine();
            bw.write("KmerLength:\t"+String.valueOf(kmerLength)); bw.newLine();
            bw.write("KmerNumber:\t"+String.valueOf(kmers.length)); bw.newLine();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        if (kmerFileNumber > kmers.length) kmerFileNumber = kmers.length;
        int[][] bound = FArrayUtils.getSubsetsIndicesBySubsetNumber(kmers.length, kmerFileNumber);
        ArrayList<Integer> fileIndexList = new ArrayList();
        for (int i = 0; i < bound.length; i++) fileIndexList.add(i);
        long start = System.nanoTime();
        System.out.println("Writing kmer list to hard drive...");
        fileIndexList.parallelStream().forEach(fileIndex -> {
            String outfileS = new File (kmerBlockDir, "kmer"+FStringUtils.getNDigitNumber(4, fileIndex)+".bin").getAbsolutePath();
            try {
                DataOutputStream dos = IoUtils.getBinaryWriter(outfileS);
                dos.writeInt(bound[fileIndex][1] - bound[fileIndex][0]);
                for (int i = bound[fileIndex][0]; i < bound[fileIndex][1]; i++) {
                    dos.writeLong(kmers[i]);
                }
                dos.flush();
                dos.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        System.out.println("Finished writing kmers. Time span: " + String.valueOf(Benchmark.getTimeSpanSeconds(start)) + " seconds. Memory used: "+String.valueOf(Benchmark.getUsedMemoryGb())+" Gb");
    }
    
    private short[] getShortKmerList (int kmerLength, String referenceGenomeFileS) {
        Fasta f = new Fasta(referenceGenomeFileS);
        f.sortRecordByNameValue();
        HashByteByteMap ascIIByteMap = BaseCoder.getAscIIByteMap();
        System.out.println("Building kmer list from reference...");
        System.out.println("KmerLength = "+String.valueOf(kmerLength)+ " bp");
        int totalLength = (int)f.getTotalSeqLength();
        HashShortSet kmerSet = HashShortSets.newMutableSet(totalLength);
        long start = System.nanoTime();
        for (int k = 0; k < f.getSeqNumber(); k++) {
            String seq = f.getSeq(k);
            byte[] bArray = seq.getBytes();
            for (int i = 0; i < bArray.length; i++) {
                bArray[i] = ascIIByteMap.get(bArray[i]);
            }
            for (int i = 0; i < bArray.length-kmerLength+1; i++) {
                boolean flag = false;
                for (int j = i; j < i+kmerLength; j++) {
                    if (bArray[j] >3) {
                        i = j;
                        flag = true;
                        break;
                    }
                }
                if (flag) continue;
                byte[] subArray = Arrays.copyOfRange(bArray, i, i + kmerLength);
                short kmerL = BaseCoder.getShortSeqFromByteArray(subArray);
                if (!kmerSet.contains(kmerL)) kmerSet.add(kmerL);
                int pos = i+1;
                if (pos%10000000 == 0) {
                    System.out.println("Chromosome: "+f.getName(k)+". Length = "+String.valueOf(bArray.length)+"bp. Position: "+String.valueOf(pos));
                }
            }
        }
        StringBuilder sb = new StringBuilder();
        sb.append("Initializing reference kmer set is finished. Time span: ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(sb.toString());
        short[] kmers = kmerSet.toShortArray();
        start = System.nanoTime();
        sb = new StringBuilder();
        sb.append("Building reference kmer list is finished. Time span: ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(sb.toString());
        System.out.println("Kmer number = "+String.valueOf(kmers.length));
        start = System.nanoTime();
        System.out.println("Sorting kmers");
        Arrays.sort(kmers);
        sb = new StringBuilder();
        sb.append("Sorting kmer list is finished. Time span: ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(sb.toString());
        return kmers;
    }
    
    private int[] getIntKmerList (int kmerLength, String referenceGenomeFileS) {
        Fasta f = new Fasta(referenceGenomeFileS);
        f.sortRecordByNameValue();
        HashByteByteMap ascIIByteMap = BaseCoder.getAscIIByteMap();
        System.out.println("Building kmer list from reference...");
        System.out.println("KmerLength = "+String.valueOf(kmerLength)+ " bp");
        int totalLength = (int)f.getTotalSeqLength();
        HashIntSet kmerSet = HashIntSets.newMutableSet(totalLength);
        long start = System.nanoTime();
        for (int k = 0; k < f.getSeqNumber(); k++) {
            String seq = f.getSeq(k);
            byte[] bArray = seq.getBytes();
            for (int i = 0; i < bArray.length; i++) {
                bArray[i] = ascIIByteMap.get(bArray[i]);
            }
            for (int i = 0; i < bArray.length-kmerLength+1; i++) {
                boolean flag = false;
                for (int j = i; j < i+kmerLength; j++) {
                    if (bArray[j] >3) {
                        i = j;
                        flag = true;
                        break;
                    }
                }
                if (flag) continue;
                byte[] subArray = Arrays.copyOfRange(bArray, i, i + kmerLength);
                int kmerL = BaseCoder.getIntSeqFromByteArray(subArray);
                if (!kmerSet.contains(kmerL)) kmerSet.add(kmerL);
                int pos = i+1;
                if (pos%10000000 == 0) {
                    System.out.println("Chromosome: "+f.getName(k)+". Length = "+String.valueOf(bArray.length)+"bp. Position: "+String.valueOf(pos));
                }
            }
        }
        StringBuilder sb = new StringBuilder();
        sb.append("Initializing reference kmer set is finished. Time span: ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(sb.toString());
        int[] kmers = kmerSet.toIntArray();
        start = System.nanoTime();
        sb = new StringBuilder();
        sb.append("Building reference kmer list is finished. Time span: ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(sb.toString());
        System.out.println("Kmer number = "+String.valueOf(kmers.length));
        start = System.nanoTime();
        System.out.println("Sorting kmers");
        Arrays.sort(kmers);
        sb = new StringBuilder();
        sb.append("Sorting kmer list is finished. Time span: ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(sb.toString());
        return kmers;
    }
    
    private long[] getLongKmerList (int kmerLength, String referenceGenomeFileS) {
        Fasta f = new Fasta(referenceGenomeFileS);
        f.sortRecordByNameValue();
        HashByteByteMap ascIIByteMap = BaseCoder.getAscIIByteMap();
        System.out.println("Building kmer list from reference...");
        System.out.println("KmerLength = "+String.valueOf(kmerLength)+ " bp");
        int totalLength = (int)f.getTotalSeqLength();
        HashLongSet kmerSet = HashLongSets.newMutableSet(totalLength);
        long start = System.nanoTime();
        for (int k = 0; k < f.getSeqNumber(); k++) {
            String seq = f.getSeq(k);
            byte[] bArray = seq.getBytes();
            for (int i = 0; i < bArray.length; i++) {
                bArray[i] = ascIIByteMap.get(bArray[i]);
            }
            for (int i = 0; i < bArray.length-kmerLength+1; i++) {
                boolean flag = false;
                for (int j = i; j < i+kmerLength; j++) {
                    if (bArray[j] >3) {
                        i = j;
                        flag = true;
                        break;
                    }
                }
                if (flag) continue;
                byte[] subArray = Arrays.copyOfRange(bArray, i, i + kmerLength);
                long kmerL = BaseCoder.getLongSeqFromByteArray(subArray);
                if (!kmerSet.contains(kmerL)) kmerSet.add(kmerL);
                int pos = i+1;
                if (pos%10000000 == 0) {
                    System.out.println("Chromosome: "+f.getName(k)+". Length = "+String.valueOf(bArray.length)+"bp. Position: "+String.valueOf(pos));
                }
            }
        }
        StringBuilder sb = new StringBuilder();
        sb.append("Initializing reference kmer set is finished. Time span: ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(sb.toString());
        long[] kmers = kmerSet.toLongArray();
        start = System.nanoTime();
        sb = new StringBuilder();
        sb.append("Building reference kmer list is finished. Time span: ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(sb.toString());
        System.out.println("Kmer number = "+String.valueOf(kmers.length));
        start = System.nanoTime();
        System.out.println("Sorting kmers");
        Arrays.sort(kmers);
        sb = new StringBuilder();
        sb.append("Sorting kmer list is finished. Time span: ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(sb.toString());
        return kmers;
    }
}
