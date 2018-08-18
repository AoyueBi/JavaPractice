/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.maf;

import format.Fasta;
import format.Range;
import format.RangeAttribute;
import format.Ranges;
import format.ShortreadAlignment;
import format.Table;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TLongObjectMap;
import gnu.trove.map.hash.TLongObjectHashMap;
import graphcis.r.DensityPlot;
import graphcis.r.DensityPlotMultiClass;
import graphcis.r.Histogram;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.LongAdder;
import java.util.function.BiPredicate;
import java.util.stream.IntStream;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GenomeSequence;
import net.maizegenetics.dna.map.GenomeSequenceBuilder;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import utils.FArrayUtils;
import utils.FrequentFileUtils;
import utils.IOFileFormat;
import utils.IoUtils;


/**
 *
 * @author Fei Lu
 */
class UniqueKmers {
    
    public UniqueKmers () {
        this.rangePipe();
        //this.depthPipe2();
    }
    
    public void depthPipe2() {
        //this.creatDepthOfUniqueKmersInGenomes();
        //this.creatDepthFigureOfUniqueKmerInGenomes();
        //this.creatDepthOfUniqueGenomeOfB73();
        //creatDepthFigureOfUniqueGenomeOfB73();
    }
    
    public void creatDepthFigureOfUniqueGenomeOfB73 () {
        String infileS = "M:\\production\\maf\\wgs\\uniqueKmer\\depth2\\depthOfB73UniqueGenome.txt";
        //String outfileS = "M:\\production\\maf\\wgs\\uniqueKmer\\depth2\\depthOfB73UniqueGenome.pdf";
        String outfileS = "E:\\Research\\wgs_maf\\lowcopyGenome\\depthOfB73UniqueGenome.pdf";
        Table t = new Table (infileS);
        double[] v = new double[t.getRowNumber()] ;
        for (int i = 0; i < t.getRowNumber(); i++) {
            v[i] = Double.valueOf(t.content[i][0]);
        }
        DensityPlot d = new DensityPlot (v);
        d.setTitle("HM3.1 read depth in unique portion of genome");
        d.setYLab("Density");
        d.setXLab("Read depth");
        d.setXLim(0, 5000);
        d.saveGraph(outfileS);
    }
    
    public void creatDepthOfUniqueGenomeOfB73 () {
        String infileS = "M:\\production\\maf\\wgs\\uniqueKmer\\position\\uniqueGenome.pos.r.txt";
        String depthFileS = "N:\\HapMap3\\hmp31_allele_depths\\thr.all_c10_DEPTH_hmp31_bcq10_q30.gz";
        String outfileS = "M:\\production\\maf\\wgs\\uniqueKmer\\depth2\\depthOfB73UniqueGenome.txt";
        int chromosome = 10;
        int size = 10000;
        Ranges r = new Ranges(infileS, IOFileFormat.Text);
        r = r.getRangesByChromosome(chromosome);
        TIntArrayList dList = new TIntArrayList();
        try {
            BufferedReader br = IoUtils.getTextGzipReader(depthFileS);
            String temp = br.readLine();
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                if (!r.isInRanges(10, Integer.valueOf(tem[1]))) continue;
                int d = Integer.valueOf(tem[3]) + Integer.valueOf(tem[5]) + Integer.valueOf(tem[7]) + Integer.valueOf(tem[9]);
                dList.add(d);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        int[] depth = dList.toArray();
        FArrayUtils.shuffleArray(depth);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("DepthOfB73UniqueGenome");
            bw.newLine();
            for (int i = 0; i < size; i++) {
                bw.write(String.valueOf(depth[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void creatDepthFigureOfUniqueKmerInGenomes () {
        String infileS = "M:\\production\\maf\\wgs\\uniqueKmer\\depth2\\depthByUniqueKmerCount.txt";
        //String outfileS = "M:\\production\\maf\\wgs\\uniqueKmer\\depth2\\depthByUniqueKmerCount.pdf";
        String outfileS = "E:\\Research\\wgs_maf\\lowcopyGenome\\depthByUniqueKmerCount.pdf";
        Table t = new Table (infileS);
        double[][] values = new double[t.getColumnNumber()][t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            for (int j = 0; j < t.getColumnNumber(); j++) {
                values[j][i] = Double.valueOf(t.content[i][j]);
            }
        }
        DensityPlotMultiClass d = new DensityPlotMultiClass(values, t.header);
        d.setXLab("Read depth");
        d.setYLab("Density");
        d.setTitle("Read depth distribution in unique portion of genomes");
        d.setXLim(0, 5000);
        d.saveGraph(outfileS);
    }
    
    public void creatDepthOfUniqueKmersInGenomes () {
        String depthFileS = "N:\\HapMap3\\hmp31_allele_depths\\thr.all_c10_DEPTH_hmp31_bcq10_q30.gz";
        String kmerFileS = "M:\\production\\maf\\wgs\\uniqueKmer\\position\\kmerCountInGenomes.txt";
        String outfileS = "M:\\production\\maf\\wgs\\uniqueKmer\\depth2\\depthByUniqueKmerCount.txt";
        ArrayList<Range>[] rLists = null;
        int rangeCapasity = 200000;
        int depthSize = 10000;
        try {
            BufferedReader br = IoUtils.getTextReader(kmerFileS);
            String temp = br.readLine();
            rLists = new ArrayList[temp.split("\t").length-2];
            for (int i = 0; i < rLists.length; i++) rLists[i] = new ArrayList();
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                if (Integer.valueOf(tem[0])!=10) continue;
                Range r = new Range(Integer.valueOf(tem[0]), Integer.valueOf(tem[1]), Integer.valueOf(tem[1])+32);
                int index = -1;
                for (int i = 0; i < rLists.length; i++) {
                    index += Integer.valueOf(tem[i+2]);
                }
                rLists[index].add(r);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        for (int i = 0; i < rLists.length; i++) {
            System.out.println(rLists[i].size());
        }
        Ranges[] rs = new Ranges[rLists.length];
        TIntArrayList[] depthLists = new TIntArrayList[rLists.length];
        for (int i = 0; i < rs.length; i++) {
            rs[i] = new Ranges(rLists[i].subList(0, rangeCapasity), String.valueOf(i+1)+"occurenceInGenomes");
            rs[i].sortByStartPosition();
            rs[i].collapse();
            depthLists[i] = new TIntArrayList();
        }
        
        try {
            BufferedReader br = IoUtils.getTextGzipReader(depthFileS);
            String temp = br.readLine();
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                int d = Integer.valueOf(tem[3]) + Integer.valueOf(tem[5]) + Integer.valueOf(tem[7]) + Integer.valueOf(tem[9]);
                for (int i = 0; i < rs.length; i++) {
                    if (rs[i].isInRanges(10, Integer.valueOf(tem[1]))) {
                        depthLists[i].add(d);
                    }
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        int[][] depth = new int[depthLists.length][];
        for (int i = 0; i < depthLists.length; i++) {
            int[] temp = depthLists[i].toArray();
            FArrayUtils.shuffleArray(temp);
            depth[i] = new int[depthSize];
            System.arraycopy(temp, 0, depth[i], 0, depthSize);
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("B73\t2Genomes\t3Genomes\t4Genomes");
            bw.newLine();
            for (int i = 0; i < depth[0].length; i++) {
                StringBuilder sb = new StringBuilder();
                for (int j = 0; j < depth.length; j++) {
                    sb.append(depth[j][i]).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
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
    
    public void rangePipe () {
        //this.creatB73UniqueRangeMap();
        //this.creatSubKmerCountInGenomes();
        //this.creatFullKmerCountInGenomes();
        //this.creatB73UniqueRangeSort(); Same result with creatB73UniqueRangeMap()
        //this.creatUniqueGenomeRange();
        //this.breakdownUniqueGenomeRange();
        //this.getUniqueGenomeStatistics();
    }
    
    public void getUniqueGenomeStatistics () {
        String mergeRangeFileS = "M:\\production\\maf\\wgs\\uniqueKmer\\position\\uniqueGenome.pos.r.txt";
        String mergeRangesStaFileS = "M:\\production\\maf\\wgs\\uniqueKmer\\position\\uniqueGenome.sta.txt";
        Table t = new Table (FrequentFileUtils.B73GenomeInfoFileS);
        Ranges r = new Ranges(mergeRangeFileS, IOFileFormat.Text);
        int[] chromosomes = r.getChromosomes();
        int[] nRange = new int[chromosomes.length];
        int[] length = new int[chromosomes.length];
        for (int i = 0; i < r.getRangeNumber(); i++) {
            int index = Arrays.binarySearch(chromosomes, r.getRangeChromosome(i));
            nRange[index]++;
            length[index]+=r.getRangeSize(i);
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(mergeRangesStaFileS);
            bw.write("Chr\tRangeNumber\tUqniueLength\tGenomeLength\tUniquePortion");
            bw.newLine();
            for (int i = 0; i < chromosomes.length; i++) {
                if (i < 1 || i > 10) continue;
                int cLength = Integer.valueOf(t.content[i-1][1]);
                bw.write(String.valueOf(chromosomes[i])+"\t"+String.valueOf(nRange[i])+"\t"+String.valueOf(length[i])+"\t"+String.valueOf(cLength)+"\t"+String.valueOf((double)length[i]/cLength));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void breakdownUniqueGenomeRange () {
        //String infileS = "M:\\production\\maf\\wgs\\uniqueKmer\\position\\kmerCountInGenomes_full.txt.gz";
        //String outDirS = "M:\\production\\maf\\wgs\\uniqueKmer\\position\\uniqueGenomeBreakdown\\";
        String infileS = "/workdir/fl262/kmerCountInGenomes_full.txt.gz";
        String outDirS = "/workdir/fl262/uniqueGenomeBreakdown/";
        ArrayList<Range>[] rList = new ArrayList[4];
        for (int i = 0; i < rList.length; i++) rList[i] = new ArrayList();
        try {
            BufferedReader br = IoUtils.getTextGzipReader(infileS);
            String temp = br.readLine();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                int pos = Integer.valueOf(tem[1]);
                Range r = new Range(Integer.valueOf(tem[0]), pos, pos + 32);
                int index = Integer.valueOf(tem[2]) + Integer.valueOf(tem[3]) + Integer.valueOf(tem[4]) + Integer.valueOf(tem[5])-1;
                rList[index].add(r);
                cnt++;
                if (cnt%1000000 == 0 && cnt != 0) System.out.println(String.valueOf(cnt)+" lines scanned");
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        ArrayList<Ranges> rs = new ArrayList();
        for (int i = 0; i < rList.length; i++) rs.add(new Ranges(rList[i], "uniqueGenome_"+String.valueOf(i+1)+"_Genomes"));
        rs.parallelStream().forEach(e -> {
            e.collapse();
        });
        for (int i = 0; i < rs.size(); i++) {
            String outfileS = outDirS+"/uniqueGenome"+"_"+String.valueOf(i+1)+".r.txt";
            rs.get(i).writeFile(outfileS, IOFileFormat.Text);
        }
        System.out.println("Mission completed");
    }
    
    public void creatUniqueGenomeRange () {
        String uniqueRangeFileS = "M:\\production\\maf\\wgs\\uniqueKmer\\position\\uniqueKmer.pos.r.bin";
        String mergeRangeFileS = "M:\\production\\maf\\wgs\\uniqueKmer\\position\\uniqueGenome.pos.r.txt";
        //String uniqueRangeFileS = "/workdir/mingh/uniqueKmer.pos.r.bin";
        //String mergeRangeFileS = "/workdir/mingh/uniqueGenome.pos.r.txt";
        Ranges r = new Ranges (uniqueRangeFileS, IOFileFormat.Binary);
        r.collapse();
        r.writeFile(mergeRangeFileS, IOFileFormat.Text);
    }
    
    public void creatB73UniqueRangeSort() {
        String uniqueRangeFileS = "M:\\production\\maf\\wgs\\uniqueKmer\\position\\test\\uniqueKmer.pos.r.bin";
        String b73DirS = "M:\\Database\\maize\\agpV3\\byChromosome\\";
        //String uniqueRangeFileS = "/workdir/mingh/uniqueKmer.pos.r.bin";
        //String b73DirS = "/workdir/mingh/B73/";
        File[] f = new File(b73DirS).listFiles();
        Arrays.sort(f);
        GenomeSequence[] gs = new GenomeSequence[f.length];
        for (int i = 0; i < f.length; i++) {
            gs[i] = GenomeSequenceBuilder.instance(f[i].getAbsolutePath());
        }
        KmerInfo[] kInfo = new KmerInfo[1_200_000_000];
        KmerInfo nullInfo = new KmerInfo(Long.MAX_VALUE, (byte)-1, -1, Byte.MAX_VALUE);
        for (int i = 0; i < kInfo.length; i++) kInfo[i] = nullInfo;
        int kmerLength = 32;
        int current = 0;
        for (int i = 0; i < f.length; i++) {
            GenomeSequence genome = gs[i];
            Chromosome chrom = genome.chromosomes().toArray(new Chromosome[genome.chromosomes().size()])[0];
            System.out.println("Processing " + chrom.getName());
            byte chr = Byte.valueOf(chrom.getName());
            byte[] bArray = genome.chromosomeSequence(chrom);
            for (int j = 0; j < bArray.length-kmerLength; j++) {
                if (!rightKmerPrefix.test(bArray, j)) continue;
                boolean flag = false;
                for (int k = j; k < j+kmerLength; k++) {
                    if (bArray[k] > 3) {
                        j = k+1;
                        flag = true;
                        break;
                    }
                }
                if (flag) continue;
                long seqL = BaseEncoder.getLongSeqFromByteArray(Arrays.copyOfRange(bArray, j, j + kmerLength));
                
                KmerInfo info = new KmerInfo(seqL, chr, j+1, (byte)1);
                kInfo[current] = info;
                current++;
                if (current%100000 == 0) {
                    System.out.println("Chromosome position index:\t" + String.valueOf(j));
                    System.out.println("Current size:\t" + String.valueOf(current));
                }
            }
            System.out.println("Start sorting");
            Arrays.sort(kInfo, 0, current);
            int cnt = 0;
            for (int k = 0; k < current-1; k++) {
                if (kInfo[k].kmer == kInfo[k+1].kmer) {
                    int d = kInfo[k+1].depth + kInfo[k].depth;
                    if (d>Byte.MAX_VALUE) d = Byte.MAX_VALUE;
                    kInfo[k+1].depth = (byte)d;
                    kInfo[k] = nullInfo;
                    cnt++;
                }
            }
            Arrays.sort(kInfo, 0, current);
            System.out.println("Finished chromosome "+String.valueOf(chr));
            System.out.println("Current size "+String.valueOf(current));
            current = current-cnt;
            System.out.println("Reduce size to "+String.valueOf(current));
        }

        for (int i = 0; i < current; i++) {
            long reKey = BaseEncoder.getReverseComplement(kInfo[i].kmer);
            int index = Arrays.binarySearch(kInfo, new KmerInfo(reKey));
            if (index < 0) continue;
            kInfo[index].addDepth();
        }
        int cnt = 0;
        ArrayList<Range> rList = new ArrayList();
        for (int i = 0; i < current; i++) {
            if (kInfo[i].depth > 1) {
                cnt++;
            }
            else {
                Range r = new Range(kInfo[i].chr, kInfo[i].pos, kInfo[i].pos+kmerLength);
                rList.add(r);
            }
        }
        System.out.println("Total kmer\t"+current);
        System.out.println("Unique kmer\t"+(current-cnt));
        System.out.println("Nonunique kmer\t"+cnt);
        Ranges r = new Ranges(rList, "B73UniqueKmerPosition");
        r.sortByStartPosition();
        r.writeFile(uniqueRangeFileS, IOFileFormat.Binary);
    }
    
    private class KmerInfo implements Comparable<KmerInfo> {
        long kmer;
        byte chr = -1;
        int pos = -1;
        byte depth = 0;
        
        public KmerInfo (long kmer) {
            this.kmer = kmer;
        }
        
        public KmerInfo (long kmer, byte chr, int pos, byte depth) {
            this.kmer = kmer;
            this.chr = chr;
            this.pos = pos;
            this.depth = depth;
        }
        public void addDepth() {
            if (depth != Byte.MAX_VALUE) depth++;
        }

        @Override
        public int compareTo(KmerInfo o) {
            if (this.kmer == o.kmer) return 0;
            else if (this.kmer < o.kmer) return -1;
            return 1;
        }
    }
    
    public void creatSubKmerCountInGenomes () {
        String b73DirS = "/workdir/fl262/B73/";
        String genomeDirS = "/workdir/fl262/other/";
        String outfileS = "/workdir/fl262/kmerCountInGenomes.txt";
        File[] f = new File(b73DirS).listFiles();
        File[] others = new File(genomeDirS).listFiles();
        ArrayList<GenomeSequence> chromList = new ArrayList();
        for (int i = 0; i < f.length; i++) {
            GenomeSequence g = GenomeSequenceBuilder.instance(f[i].getAbsolutePath());
            chromList.add(g);
        }
        int kmerLength = 32;
        ConcurrentHashMap<Long,PosInfoCount> kmerMap = new ConcurrentHashMap(Integer.MAX_VALUE, (float)0.99);
        chromList.stream().forEach(genome -> {
            Chromosome chrom = genome.chromosomes().toArray(new Chromosome[genome.chromosomes().size()])[0];
            int chr = Integer.valueOf(chrom.getName());
            byte[] bArray = genome.chromosomeSequence(chrom);
            for (int j = 0; j < bArray.length-kmerLength; j++) {
                if (!atKmerPrefix.test(bArray, j)) continue;
                boolean flag = false;
                for (int k = j; k < j+kmerLength; k++) {
                    if (bArray[k] > 3) {
                        j = k+1;
                        flag = true;
                        break;
                    }
                }
                if (flag) continue;
                long seqL = BaseEncoder.getLongSeqFromByteArray(Arrays.copyOfRange(bArray, j, j + kmerLength));
                if (kmerMap.containsKey(seqL)) {
                    kmerMap.get(seqL).increase(0);
                } 
                else {
                    PosInfoCount p = new PosInfoCount (chr, j+1, others.length+1);
                    kmerMap.put(seqL, p);
                }
            }
            bArray = this.getReverseComplementary(bArray);
            for (int j = 0; j < bArray.length-kmerLength; j++) {
                if (!atKmerPrefix.test(bArray, j)) continue;
                boolean flag = false;
                for (int k = j; k < j+kmerLength; k++) {
                    if (bArray[k] > 3) {
                        j = k+1;
                        flag = true;
                        break;
                    }
                }
                if (flag) continue;
                long seqL = BaseEncoder.getLongSeqFromByteArray(Arrays.copyOfRange(bArray, j, j + kmerLength));
                if (kmerMap.containsKey(seqL)) {
                    kmerMap.get(seqL).increase(0);
                } 
                else {
                    PosInfoCount p = new PosInfoCount (chr, bArray.length-j-kmerLength+1, others.length+1);
                    kmerMap.put(seqL, p);
                }
            }
            System.out.println("Finished chromosome "+String.valueOf(chr));
            System.out.println("Current size "+String.valueOf(kmerMap.size()));
        });
        int window = 2_000_000;
        for (int i = 0; i < others.length; i++) {
            GenomeSequence genomeSequence = GenomeSequenceBuilder.instance(others[i].getAbsolutePath());
            for (long start = 0; start < genomeSequence.genomeSize(); start += window) {
                System.out.printf("Subseq start %,d %n", start);
                byte[] subSeq = genomeSequence.genomeSequence(start, Math.min(genomeSequence.genomeSize(), start + window - 1));
                for (int j = 0; j < subSeq.length - kmerLength; j++) {
                    if (!atKmerPrefix.test(subSeq, j)) continue;
                    boolean flag = false;
                    for (int k = j; k < j+kmerLength; k++) {
                        if (subSeq[k] > 3) {
                            j = k+1;
                            flag = true;
                            break;
                        }
                    }
                    if (flag) continue;
                    //System.out.println("i = " + i);
                    long seqL = BaseEncoder.getLongSeqFromByteArray(Arrays.copyOfRange(subSeq, j, j + kmerLength));
                    if (kmerMap.containsKey(seqL)) {
                        kmerMap.get(seqL).increase(i+1);
                    } 
                }
                subSeq = this.getReverseComplementary(subSeq);
                for (int j = 0; j < subSeq.length - kmerLength; j++) {
                    if (!atKmerPrefix.test(subSeq, i)) continue;
                    boolean flag = false;
                    for (int k = j; k < j+kmerLength; k++) {
                        if (subSeq[k] > 3) {
                            j = k+1;
                            flag = true;
                            break;
                        }
                    }
                    if (flag) continue;
                    //System.out.println("i = " + i);
                    long seqL = BaseEncoder.getLongSeqFromByteArray(Arrays.copyOfRange(subSeq, j, j + kmerLength));
                    if (kmerMap.containsKey(seqL)) {
                        kmerMap.get(seqL).increase(i+1);
                    }
                }
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Chr\tPos\tB73");
            for (int i = 0; i < others.length; i++) {
                bw.write("\t"+others[i].getName());
            }
            bw.newLine();
            for (Map.Entry<Long,PosInfoCount> e : kmerMap.entrySet()) {
                byte[] count = e.getValue().depth;
                boolean flag = false;
                for (int i = 0; i < count.length; i++) {
                    if (count[i] > 1) {
                        flag = true;
                        break;
                    }
                }
                if (flag) continue;
                StringBuilder sb = new StringBuilder();
                sb.append(e.getValue().chr).append("\t").append(e.getValue().pos).append("\t");
                for (int i = 0; i < count.length; i++) {
                    sb.append(count[i]).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
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
    
    public void creatFullKmerCountInGenomes () {
        String b73DirS = "/workdir/fl262/B73/";
        String genomeDirS = "/workdir/fl262/other/";
        String outfileS = "/workdir/fl262/kmerCountInGenomes_full.txt";
        File[] f = new File(b73DirS).listFiles();
        File[] others = new File(genomeDirS).listFiles();
        ArrayList<GenomeSequence> chromList = new ArrayList();
        for (int i = 0; i < f.length; i++) {
            GenomeSequence g = GenomeSequenceBuilder.instance(f[i].getAbsolutePath());
            chromList.add(g);
        }
        int kmerLength = 32;
        ConcurrentHashMap<Long,PosInfoCount> kmerMap = new ConcurrentHashMap(Integer.MAX_VALUE, (float)0.99);
        chromList.stream().forEach(genome -> {
            Chromosome chrom = genome.chromosomes().toArray(new Chromosome[genome.chromosomes().size()])[0];
            int chr = Integer.valueOf(chrom.getName());
            byte[] bArray = genome.chromosomeSequence(chrom);
            for (int j = 0; j < bArray.length-kmerLength; j++) {
                //if (!atKmerPrefix.test(bArray, j)) continue;
                boolean flag = false;
                for (int k = j; k < j+kmerLength; k++) {
                    if (bArray[k] > 3) {
                        j = k+1;
                        flag = true;
                        break;
                    }
                }
                if (flag) continue;
                long seqL = BaseEncoder.getLongSeqFromByteArray(Arrays.copyOfRange(bArray, j, j + kmerLength));
                if (kmerMap.containsKey(seqL)) {
                    kmerMap.get(seqL).increase(0);
                } 
                else {
                    PosInfoCount p = new PosInfoCount (chr, j+1, others.length+1);
                    kmerMap.put(seqL, p);
                }
            }
            bArray = this.getReverseComplementary(bArray);
            for (int j = 0; j < bArray.length-kmerLength; j++) {
                //if (!atKmerPrefix.test(bArray, j)) continue;
                boolean flag = false;
                for (int k = j; k < j+kmerLength; k++) {
                    if (bArray[k] > 3) {
                        j = k+1;
                        flag = true;
                        break;
                    }
                }
                if (flag) continue;
                long seqL = BaseEncoder.getLongSeqFromByteArray(Arrays.copyOfRange(bArray, j, j + kmerLength));
                if (kmerMap.containsKey(seqL)) {
                    kmerMap.get(seqL).increase(0);
                } 
                else {
                    PosInfoCount p = new PosInfoCount (chr, bArray.length-j-kmerLength+1, others.length+1);
                    kmerMap.put(seqL, p);
                }
            }
            System.out.println("Finished chromosome "+String.valueOf(chr));
            System.out.println("Current size "+String.valueOf(kmerMap.size()));
        });
        int window = 2_000_000;
        for (int i = 0; i < others.length; i++) {
            GenomeSequence genomeSequence = GenomeSequenceBuilder.instance(others[i].getAbsolutePath());
            for (long start = 0; start < genomeSequence.genomeSize(); start += window) {
                System.out.printf("Subseq start %,d %n", start);
                byte[] subSeq = genomeSequence.genomeSequence(start, Math.min(genomeSequence.genomeSize(), start + window - 1));
                for (int j = 0; j < subSeq.length - kmerLength; j++) {
                    //if (!atKmerPrefix.test(subSeq, j)) continue;
                    boolean flag = false;
                    for (int k = j; k < j+kmerLength; k++) {
                        if (subSeq[k] > 3) {
                            j = k+1;
                            flag = true;
                            break;
                        }
                    }
                    if (flag) continue;
                    //System.out.println("i = " + i);
                    long seqL = BaseEncoder.getLongSeqFromByteArray(Arrays.copyOfRange(subSeq, j, j + kmerLength));
                    if (kmerMap.containsKey(seqL)) {
                        kmerMap.get(seqL).increase(i+1);
                    } 
                }
                subSeq = this.getReverseComplementary(subSeq);
                for (int j = 0; j < subSeq.length - kmerLength; j++) {
                    //if (!atKmerPrefix.test(subSeq, i)) continue;
                    boolean flag = false;
                    for (int k = j; k < j+kmerLength; k++) {
                        if (subSeq[k] > 3) {
                            j = k+1;
                            flag = true;
                            break;
                        }
                    }
                    if (flag) continue;
                    //System.out.println("i = " + i);
                    long seqL = BaseEncoder.getLongSeqFromByteArray(Arrays.copyOfRange(subSeq, j, j + kmerLength));
                    if (kmerMap.containsKey(seqL)) {
                        kmerMap.get(seqL).increase(i+1);
                    }
                }
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Chr\tPos\tB73");
            for (int i = 0; i < others.length; i++) {
                bw.write("\t"+others[i].getName());
            }
            bw.newLine();
            for (Map.Entry<Long,PosInfoCount> e : kmerMap.entrySet()) {
                byte[] count = e.getValue().depth;
                boolean flag = false;
                for (int i = 0; i < count.length; i++) {
                    if (count[i] > 1) {
                        flag = true;
                        break;
                    }
                }
                if (flag) continue;
                StringBuilder sb = new StringBuilder();
                sb.append(e.getValue().chr).append("\t").append(e.getValue().pos).append("\t");
                for (int i = 0; i < count.length; i++) {
                    sb.append(count[i]).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
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
    
    private class PosInfoCount {
        int chr;
        int pos;
        byte[] depth;
        public PosInfoCount(int chr, int pos, int genomeNumber) {
            this.chr = chr;
            this.pos = pos;
            depth = new byte[genomeNumber];
            depth[0] = 1;
        }
        public void increase (int index) {
            int a = depth[index]+1;
            if (a > Byte.MAX_VALUE) a = Byte.MAX_VALUE;
            depth[index] = (byte)a;
        }
    }
    
    //takes too much memory, needs at least 150Gb mem
    public void creatB73UniqueRangeMap () {
        //String uniqueRangeFileS = "M:\\production\\maf\\wgs\\uniqueKmer\\position\\test\\uniqueKmer.pos.r.bin";
        //String b73DirS = "M:\\Database\\maize\\agpV3\\byChromosome\\";
        String uniqueRangeFileS = "/workdir/fl262/uniqueKmer.pos.r.bin";
        String b73DirS = "/workdir/fl262/B73/";
        File[] f = new File(b73DirS).listFiles();
        ArrayList<GenomeSequence> chromList = new ArrayList();
        for (int i = 0; i < f.length; i++) {
            GenomeSequence g = GenomeSequenceBuilder.instance(f[i].getAbsolutePath());
            chromList.add(g);
        }
        int kmerLength = 32;
        ConcurrentHashMap<Long,PosInfo> kmerMap = new ConcurrentHashMap(Integer.MAX_VALUE, (float)0.99);
        chromList.stream().forEach(genome -> {
            Chromosome chrom = genome.chromosomes().toArray(new Chromosome[genome.chromosomes().size()])[0];
            int chr = Integer.valueOf(chrom.getName());
            byte[] bArray = genome.chromosomeSequence(chrom);
            for (int j = 0; j < bArray.length-kmerLength; j++) {
                //if (!rightKmerPrefix.test(bArray, j)) continue;
                boolean flag = false;
                for (int k = j; k < j+kmerLength; k++) {
                    if (bArray[k] > 3) {
                        j = k+1;
                        flag = true;
                        break;
                    }
                }
                if (flag) continue;
                long seqL = BaseEncoder.getLongSeqFromByteArray(Arrays.copyOfRange(bArray, j, j + kmerLength));
                if (kmerMap.containsKey(seqL)) {
                    kmerMap.get(seqL).addDepth();
                } 
                else {
                    PosInfo p = new PosInfo (chr, j+1, 1);
                    kmerMap.put(seqL, p);
                }
                
            }
            System.out.println("Finished chromosome "+String.valueOf(chr));
            System.out.println("Current size "+String.valueOf(kmerMap.size()));
        });
        long[] reKeys = new long[kmerMap.size()];
        int cnt = 0;
        for (Map.Entry<Long,PosInfo> e : kmerMap.entrySet()) {
            reKeys[cnt] = BaseEncoder.getReverseComplement(e.getKey());
            cnt++;
        }
        for (int i = 0; i < reKeys.length; i++) {
            if (kmerMap.containsKey(reKeys[i])) {
                kmerMap.get(reKeys[i]).addDepth();
            }
        }
        ArrayList<Range> rList = new ArrayList();
        LongAdder counter =new LongAdder();
        kmerMap.entrySet().stream().forEach(e -> {
            PosInfo p = e.getValue();
            if (p.depth != 1) {
                counter.increment();
            }
            else {
                Range r = new Range(p.chr,p.pos,p.pos+kmerLength);
                rList.add(r);
            }
        });
        cnt = counter.intValue();
        System.out.println("Total kmer\t"+kmerMap.size());
        System.out.println("Unique kmer\t"+(kmerMap.size()-cnt));
        System.out.println("Nonunique kmer\t"+cnt);
        Ranges r = new Ranges(rList, "B73UniqueKmerPosition");
        r.sortByStartPosition();
        r.writeFile(uniqueRangeFileS, IOFileFormat.Binary);
    }
    
    private class PosInfo {
        int chr;
        int pos;
        int depth = 0;
        public PosInfo (int chr, int pos, int depth) {
            this.chr = chr;
            this.pos = pos;
            this.depth = depth;
        }
        public void addDepth () {
            depth++;
        }
    }
    
    private byte[] getReverseComplementary (byte[] a) {
        byte[] b = new byte[a.length];
        for (int i = 0; i < a.length; i++) {
            byte t;
            if (a[i] == 0) t = 3;
            else if (a[i] == 1) t = 2;
            else if (a[i] == 2) t = 1;
            else if (a[i] == 3) t = 0;
            else t = a[i];
            b[a.length-i-1] = t;
        }
        return b;
    }
    private static final BiPredicate<byte[], Integer> atKmerPrefix = (bb, i) -> {
        if ((bb[i] == 0) && (bb[i + 1] == 3)) return true;
        return false;
    };
    
    private static final BiPredicate<byte[], Integer> rightKmerPrefix = (bb, i) -> {
        if ((bb[i] == 0) && (bb[i + 1] == 2) && (bb[i + 2] == 3)) return true;
        return false;
    };
}
