/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package analysis.deprecated.cassava.wgs;

import format.ShortreadAlignment;
import format.ShortreadPEAlignment;
import graphcis.r.Histogram;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.dna.read.FastqChunk;
import net.maizegenetics.dna.read.ReadUtils;
import net.maizegenetics.dna.read.ReadUtils.ReadFormat;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class CaWGSUtils {
    public CaWGSUtils () {}
    
    public void compareAligners (String samBWAMEM, String samBowtie2) {
        ShortreadAlignment saBWA = new ShortreadAlignment ();
        saBWA.readFromBWAMEM(samBWAMEM);

        ShortreadAlignment saBowtie = new ShortreadAlignment ();
        saBowtie.readFromBowtie2(samBowtie2);
        String[] queries = saBowtie.getQuerys();
        int alnCnt = 0;
        int alnSameCnt = 0;
        
        int alnBWA = 0;
        int alnBowtie = 0;
        ArrayList<Short> bwaScoreList = new ArrayList();
        ArrayList<Short> bowtieScoreList = new ArrayList();
        for (int i = 0; i < queries.length; i++) {
            int bwaIndex = saBWA.getAlignmentStartIndexByQuery(queries[i]);
            int bowtieIndex = saBowtie.getAlignmentStartIndexByQuery(queries[i]);
            if (saBWA.isMatch(bwaIndex) && saBowtie.isMatch(bowtieIndex)) {
                if (saBWA.getHit(bwaIndex).equals(saBowtie.getHit(bowtieIndex))) {
                    if (Math.abs(saBWA.getStartPos(bwaIndex)-(saBowtie.getStartPos(bowtieIndex))) > 10) {
                        alnBWA++; alnBowtie++; alnCnt+=2;
                    }
                    else {
                        alnSameCnt++; alnCnt+=1;
                    }
                }
                else {
                    alnBWA++; alnBowtie++; alnCnt+=2;
                }
                bwaScoreList.add(saBWA.getMappingQuality(i));
                bowtieScoreList.add(saBowtie.getMappingQuality(i));
            }
            else if (!saBWA.isMatch(bwaIndex) && saBowtie.isMatch(bowtieIndex)) {
                alnBowtie++; alnCnt++;
                bowtieScoreList.add(saBowtie.getMappingQuality(i));
            }
            else if (saBWA.isMatch(bwaIndex) && !saBowtie.isMatch(bowtieIndex)) {
                bwaScoreList.add(saBWA.getMappingQuality(i));
                alnBWA++; alnCnt++;
            }
        }
        Short[] bwaScore = bwaScoreList.toArray(new Short[bwaScoreList.size()]);
        Arrays.sort(bwaScore);
        Short[] bowtieScore = bowtieScoreList.toArray(new Short[bowtieScoreList.size()]);
        Arrays.sort(bowtieScore);
        System.out.println("Median bwa mapping quality:\t"+ bwaScore[bwaScore.length/2]);
        System.out.println("Median bowtie2 mapping quality:\t"+ bowtieScore[bowtieScore.length/2]);
        System.out.println("Total alignment:\t" + String.valueOf(alnCnt));
        System.out.println("Shared alignment:\t" + String.valueOf(alnSameCnt) + "\t" +String.valueOf((double)alnSameCnt/alnCnt));
        System.out.println("BWA-MEM unique alignment:\t" + String.valueOf(alnBWA) + "\t" +String.valueOf((double)alnBWA/alnCnt));
        System.out.println("Bowtie2 unique alignment:\t" + String.valueOf(alnBowtie) + "\t" +String.valueOf((double)alnBowtie/alnCnt));
    }
    
    public void mkContigSizeDistribution (String contigFileS, String pdfFileS) {
        FastqChunk fq = new FastqChunk (contigFileS, ReadFormat.FastqGzip);
        int[] sizes = new int[fq.getReadNum()];
        for (int i = 0; i < sizes.length; i++) {
            sizes[i] = fq.getRead(i).getReadLength();
        }
        Histogram h = new Histogram(sizes);
        h.setTitle("Contig size distribution");
        h.setXLab("Contig size (bp)");
        h.saveGraph(pdfFileS);
    }
    
    /**
     * Using flash
     * @param fastqDirS
     * @param mergeDirS
     * @param perlScript 
     */
    public void mergePE (String fastqDirS, String mergeDirS, String perlScript) {
        File fastqF = new File(fastqDirS);
        File[] fs = fastqF.listFiles();
        ArrayList<String> prefixList = new ArrayList();
        for (int i = 0; i < fs.length; i++) {
            String name = fs[i].getName();
            if (!name.endsWith("fastq.gz")) continue;
            int index = name.indexOf("_R1");
            if (index < 0) continue;
            prefixList.add(name.substring(0, index));
        }
        String[] prefix = prefixList.toArray(new String[prefixList.size()]);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(perlScript);
            for (int i = 0; i < prefix.length; i++) {
                String file1 = new File(fastqF, prefix[i]+ "_R1.fastq.gz").getAbsolutePath();
                String file2 = new File(fastqF, prefix[i]+ "_R2.fastq.gz").getAbsolutePath();
                String command = "system\"flash " + file1 + " " + file2 + " -M 180 -x 0.38 -z -d " + mergeDirS + " -o " + prefix[i] + "\";";
                bw.write(command);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    public void alignmentBWAMEMStatistics (String samDirS, String scoreFigureDirS, String sizeFigureDirS, String reportFileS) {
        File[] sams = new File(samDirS).listFiles();
        new File(scoreFigureDirS).mkdir();
        new File(sizeFigureDirS).mkdir();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(reportFileS);
            bw.write("Fastq\tAlignRate\tPEAlignRate\tSEAlignRate\tQual-0\tQual-10\tQual-20\tQual-30\tQual-40\tQual-50\tQual-60\tQual-70\tQual-80\tQual-90\tQual-100\tSize-0\tSize-10\tSize-20\tSize-30\tSize-40\tSize-50\tSize-60\tSize-70\tSize-80\tSize-90\tSize-100");
            bw.newLine();
            for (int i = 0; i < sams.length; i++) {
                int bothMatchCnt = 0;
                int singleMatchCnt = 0;
                int matchCnt = 0;
                ArrayList<Integer> scoreList = new ArrayList();
                ArrayList<Integer> sizeList = new ArrayList();
                ShortreadPEAlignment spa = new ShortreadPEAlignment();
                spa.readFromBWAMEM(sams[i].getAbsolutePath());
                String prefix = sams[i].getName();
                int index = prefix.indexOf("sam");
                prefix = prefix.substring(0, index-1);
                String sizeFigureFileS = new File (sizeFigureDirS, prefix+"_size.pdf").getAbsolutePath().replaceAll("\\\\", "/");
                String scoreFigureFileS = new File (scoreFigureDirS, prefix+"_score.pdf").getAbsolutePath().replaceAll("\\\\", "/");;
                for (int j = 0; j < spa.getAlignmentNumber(); j++) {
                    if (spa.isMatchF(j) && spa.isMatchB(j)) {
                        bothMatchCnt++; matchCnt++;
                        scoreList.add((int)spa.getMappingQualityF(j));
                        scoreList.add((int)spa.getMappingQualityB(j));
                        if (spa.isOnReverseDirectionFB(j)) {
                            int size = spa.getPEFragmentSize(j);
                            if (size == -1) continue;
                            if (size > 1000) continue;
                            sizeList.add(size);
                        }
                    }
                    else if (!spa.isMatchF(index) && spa.isMatchB(index)) {
                        singleMatchCnt++; matchCnt++;
                        scoreList.add((int)spa.getMappingQualityB(j));
                    }
                    else if (spa.isMatchF(index) && !spa.isMatchB(index)) {
                        singleMatchCnt++; matchCnt++;
                        scoreList.add((int)spa.getMappingQualityF(j));
                    }
                }
                int[] scores = new int[scoreList.size()];
                for (int j = 0; j < scores.length; j++) {
                    scores[j] = scoreList.get(j);
                }
                Arrays.sort(scores);
                int[] sizes = new int[sizeList.size()];
                for (int j = 0; j < sizes.length; j++) {
                    sizes[j] = sizeList.get(j);
                }
                Arrays.sort(sizes);
                Histogram scoreQualH = new Histogram(scores);
                scoreQualH.setTitle("Mapping quality distribution of " + prefix);
                scoreQualH.setXLab("Mapping quality");
                scoreQualH.saveGraph(scoreFigureFileS);
                Histogram scoreSizeH = new Histogram(sizes);
                scoreSizeH.setBreakNumber(100);
                scoreSizeH.setTitle("Fragment size distribution of " + prefix);
                scoreSizeH.setXLab("Fragment size (bp)");
                scoreSizeH.saveGraph(sizeFigureFileS);
                bw.write(prefix+"\t"+String.valueOf((double)matchCnt/spa.getAlignmentNumber())+"\t"+String.valueOf((double)bothMatchCnt/spa.getAlignmentNumber())+"\t"+String.valueOf((double)singleMatchCnt/spa.getAlignmentNumber()));
                for (int j = 0; j < 11; j++) {
                    index = (int)((0.1*j)*scores.length);
                    if (index > scores.length -1) index = scores.length -1;
                    bw.write("\t"+String.valueOf(scores[index]));
                }
                for (int j = 0; j < 11; j++) {
                    index = (int)((0.1*j)*sizes.length);
                    if (index > sizes.length -1) index = sizes.length -1;
                    bw.write("\t"+String.valueOf(sizes[index]));
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    public void alignFastqBWAMEM (String fastqDirS, String reference, String alignmentDirS, String perlScript) {
        File fastqF = new File(fastqDirS);
        File[] fs = fastqF.listFiles();
        ArrayList<String> prefixList = new ArrayList();
        for (int i = 0; i < fs.length; i++) {
            String name = fs[i].getName();
            if (!name.endsWith("fastq.gz")) continue;
            int index = name.indexOf("_R1");
            if (index < 0) continue;
            prefixList.add(name.substring(0, index));
        }
        String[] prefix = prefixList.toArray(new String[prefixList.size()]);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(perlScript);
            for (int i = 0; i < prefix.length; i++) {
                String file1 = new File(fastqF, prefix[i]+ "_R1.fastq.gz").getAbsolutePath();
                String file2 = new File(fastqF, prefix[i]+ "_R2.fastq.gz").getAbsolutePath();
                String outfileS = new File (alignmentDirS, prefix[i]+"_sam.txt").getAbsolutePath();
                String command = "system\"bwa mem " + reference + " " + file1 + " " + file2 + " -t 23 > " + outfileS+"\";";
                bw.write(command);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    public void sliceFastq (String fastqDirS, String slicedFastqDirS, int startIndex, int readNum) {
        File[] fAll = IoUtils.listRecursiveFiles(new File(fastqDirS));
        File[] fastqs = IoUtils.listFilesEndsWith(fAll, "fastq.gz");
        for (int i = 0; i < fastqs.length; i++) {
            String path = fastqs[i].getPath();
            if (path.contains("GBS") || path.contains("RNASeq")) continue;
            FastqChunk f = new FastqChunk (fastqs[i].getAbsolutePath(), ReadFormat.FastqGzip, startIndex, readNum);
            File nf = new File (slicedFastqDirS, fastqs[i].getName());
            f.writeFastq(nf.getAbsolutePath(), ReadFormat.FastqGzip);
        }
    }
    
    
}
