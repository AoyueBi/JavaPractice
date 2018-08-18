/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.cassava.wgs;

import format.Bins;
import format.Range;
import format.RangeAttribute;
import format.Ranges;
import format.ShortreadAlignment;
import format.ShortreadPEAlignment;
import format.Table;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import graphcis.GenomeTrendPlot;
import graphcis.GenomeTrendPlotTwoYAxis;
import graphcis.r.ScatterPlot;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import net.maizegenetics.analysis.gbs.v2.GBSUtils;
import net.maizegenetics.dna.read.FastqChunk;
import net.maizegenetics.dna.read.ReadUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class SpeciesIntrogression {
    
    public SpeciesIntrogression () {
        //this.getErrorFreeHiSeq();
        //this.getErrorFreeRef();
        //this.editDistance();
        //this.paintHighDivergenceReads();
        //this.paintHighDivergeneReadsWithReference();
        //this.paintHighDivergeneReadsWithReferenceChloroplast();
    }
    
    public void paintHighDivergeneReadsWithReferenceChloroplast () {
        String infoFileS = "M:\\Database\\cassavaReference\\Manihot esculenta\\genomeInfo_Cassava.txt";
        String samDirS = "M:\\production\\cassava\\hmp\\speciesIntrogression\\sam_25_withChloroplast\\";
        String outDirS = "M:\\production\\cassava\\hmp\\speciesIntrogression\\divergentReadsPaintingWithReferenceChloroplast\\";
        String referenceSam = "M:\\production\\cassava\\hmp\\speciesIntrogression\\sam_25_withChloroplast\\refererence_sam.txt.gz";
        int binLength = 500000;
        Table t = new Table (infoFileS);
        int[] chromLength = new int[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) chromLength[i] = t.getIntValue(i, 1);
        ArrayList<Range> rList = new ArrayList();
        ShortreadPEAlignment sa = new ShortreadPEAlignment();
        sa.readFromBWAMEM(referenceSam);
        for (int j = 0; j < sa.getAlignmentNumber(); j++) {
            if (sa.isMatchF(j)) {
                if (Integer.valueOf(sa.getHitF(j))>20) continue;
                double d = (double)sa.getEditDistanceF(j)/sa.getMatchNumberF(j);
                if (d > 0.05) rList.add(new Range(Integer.valueOf(sa.getHitF(j)), sa.getStartPosF(j), sa.getStartPosF(j)+1));
            }
            if (sa.isMatchB(j)) {
                if (Integer.valueOf(sa.getHitB(j))>20) continue;
                double d = (double)sa.getEditDistanceB(j)/sa.getMatchNumberB(j);
                if (d > 0.05) rList.add(new Range(Integer.valueOf(sa.getHitB(j)), sa.getStartPosB(j), sa.getStartPosB(j)+1));
            }
        }
        Ranges r = new Ranges(rList, "highDivergentPosition");
        r.sortByStartPosition();
        int chromNum = r.getChromosomeNumber();
        int[][] cor = new int[chromNum][];
        double[][] refValue = new double[chromNum][];
        int readsNum = r.getRangeNumber();
        for (int j = 0; j < chromNum; j++) {
            Ranges crs = r.getRangesByChromosome(j+1);
            int[] pos = new int[crs.getRangeNumber()];
            double[] value = new double[crs.getRangeNumber()];
            for (int k = 0; k < pos.length; k++) {
                pos[k] = crs.getRangeStart(k);
                value[k] = 1;
            }
            Bins b = new Bins (1, chromLength[j], binLength, pos, value);
            cor[j] = new int[b.getBinNum()];
            refValue[j] = new double[b.getBinNum()];
            for (int k = 0; k < b.getBinNum(); k++) {
                cor[j][k] = (b.getBinStart(k)+b.getBinEnd(k))/2;
                refValue[j][k] = (double)b.getNumValues(k)/readsNum;
            }
        }
        File[] fs = new File(samDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            rList = new ArrayList();
            sa = new ShortreadPEAlignment();
            sa.readFromBWAMEM(fs[i].getAbsolutePath());
            for (int j = 0; j < sa.getAlignmentNumber(); j++) {
                if (sa.isMatchF(j)) {
                    if (Integer.valueOf(sa.getHitF(j))>20) continue;
                    double d = (double)sa.getEditDistanceF(j)/sa.getMatchNumberF(j);
                    if (d > 0.05) rList.add(new Range(Integer.valueOf(sa.getHitF(j)), sa.getStartPosF(j), sa.getStartPosF(j)+1));
                }
                if (sa.isMatchB(j)) {
                    if (Integer.valueOf(sa.getHitB(j))>20) continue;
                    double d = (double)sa.getEditDistanceB(j)/sa.getMatchNumberB(j);
                    if (d > 0.05) rList.add(new Range(Integer.valueOf(sa.getHitB(j)), sa.getStartPosB(j), sa.getStartPosB(j)+1));
                }
            }
            r = new Ranges(rList, "highDivergentPosition");
            r.sortByStartPosition();
            chromNum = r.getChromosomeNumber();
            double[][] vs = new double[chromNum][];
            readsNum = r.getRangeNumber();
            for (int j = 0; j < chromNum; j++) {
                Ranges crs = r.getRangesByChromosome(j+1);
                int[] pos = new int[crs.getRangeNumber()];
                double[] value = new double[crs.getRangeNumber()];
                for (int k = 0; k < pos.length; k++) {
                    pos[k] = crs.getRangeStart(k);
                    value[k] = 1;
                }
                Bins b = new Bins (1, chromLength[j], binLength, pos, value);
                cor[j] = new int[b.getBinNum()];
                vs[j] = new double[b.getBinNum()];
                for (int k = 0; k < b.getBinNum(); k++) {
                    cor[j][k] = (b.getBinStart(k)+b.getBinEnd(k))/2;
                    vs[j][k] = (double)b.getNumValues(k)/readsNum;
                }
            }
            String outfileS = new File (outDirS, fs[i].getName().replaceFirst("_sam.txt.gz", ".pdf")).getAbsolutePath();
            String title = "High divergent reads (>0.05) distribution of "+fs[i].getName()+" Vs reference(blue)";
            GenomeTrendPlot gtp = new GenomeTrendPlot (chromLength, cor, vs, title);
            gtp.addData(cor, refValue);
            gtp.setLowValueLim(0);
            gtp.setHighValueLim(0.02);
            gtp.setIfWithDot(false);
            gtp.saveGragh(outfileS);
        }
    }
    
    public void paintHighDivergeneReadsWithReference () {
        String infoFileS = "M:\\Database\\cassavaReference\\genomeInfo_Cassava.txt";
        String samDirS = "M:\\production\\cassava\\hmp\\speciesIntrogression\\sam_25\\";
        String outDirS = "M:\\production\\cassava\\hmp\\speciesIntrogression\\divergentReadsPaintingWithReference\\";
        String referenceSam = "M:\\production\\cassava\\hmp\\speciesIntrogression\\sam_25\\refererence_sam.txt.gz";
        int binLength = 500000;
        Table t = new Table (infoFileS);
        int[] chromLength = new int[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) chromLength[i] = t.getIntValue(i, 1);
        ArrayList<Range> rList = new ArrayList();
        ShortreadPEAlignment sa = new ShortreadPEAlignment();
        sa.readFromBWAMEM(referenceSam);
        for (int j = 0; j < sa.getAlignmentNumber(); j++) {
            if (sa.isMatchF(j)) {
                double d = (double)sa.getEditDistanceF(j)/sa.getMatchNumberF(j);
                if (d > 0.05) rList.add(new Range(Integer.valueOf(sa.getHitF(j)), sa.getStartPosF(j), sa.getStartPosF(j)+1));
            }
            if (sa.isMatchB(j)) {
                double d = (double)sa.getEditDistanceB(j)/sa.getMatchNumberB(j);
                if (d > 0.05) rList.add(new Range(Integer.valueOf(sa.getHitB(j)), sa.getStartPosB(j), sa.getStartPosB(j)+1));
            }
        }
        Ranges r = new Ranges(rList, "highDivergentPosition");
        r.sortByStartPosition();
        int chromNum = r.getChromosomeNumber();
        int[][] cor = new int[chromNum][];
        double[][] refValue = new double[chromNum][];
        int readsNum = r.getRangeNumber();
        for (int j = 0; j < chromNum; j++) {
            Ranges crs = r.getRangesByChromosome(j+1);
            int[] pos = new int[crs.getRangeNumber()];
            double[] value = new double[crs.getRangeNumber()];
            for (int k = 0; k < pos.length; k++) {
                pos[k] = crs.getRangeStart(k);
                value[k] = 1;
            }
            Bins b = new Bins (1, chromLength[j], binLength, pos, value);
            cor[j] = new int[b.getBinNum()];
            refValue[j] = new double[b.getBinNum()];
            for (int k = 0; k < b.getBinNum(); k++) {
                cor[j][k] = (b.getBinStart(k)+b.getBinEnd(k))/2;
                refValue[j][k] = (double)b.getNumValues(k)/readsNum;
            }
        }
        File[] fs = new File(samDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            rList = new ArrayList();
            sa = new ShortreadPEAlignment();
            sa.readFromBWAMEM(fs[i].getAbsolutePath());
            for (int j = 0; j < sa.getAlignmentNumber(); j++) {
                if (sa.isMatchF(j)) {
                    double d = (double)sa.getEditDistanceF(j)/sa.getMatchNumberF(j);
                if (d > 0.05) rList.add(new Range(Integer.valueOf(sa.getHitF(j)), sa.getStartPosF(j), sa.getStartPosF(j)+1));
                }
                if (sa.isMatchB(j)) {
                    double d = (double)sa.getEditDistanceB(j)/sa.getMatchNumberB(j);
                    if (d > 0.05) rList.add(new Range(Integer.valueOf(sa.getHitB(j)), sa.getStartPosB(j), sa.getStartPosB(j)+1));
                }
            }
            r = new Ranges(rList, "highDivergentPosition");
            r.sortByStartPosition();
            chromNum = r.getChromosomeNumber();
            double[][] vs = new double[chromNum][];
            readsNum = r.getRangeNumber();
            for (int j = 0; j < chromNum; j++) {
                Ranges crs = r.getRangesByChromosome(j+1);
                int[] pos = new int[crs.getRangeNumber()];
                double[] value = new double[crs.getRangeNumber()];
                for (int k = 0; k < pos.length; k++) {
                    pos[k] = crs.getRangeStart(k);
                    value[k] = 1;
                }
                Bins b = new Bins (1, chromLength[j], binLength, pos, value);
                cor[j] = new int[b.getBinNum()];
                vs[j] = new double[b.getBinNum()];
                for (int k = 0; k < b.getBinNum(); k++) {
                    cor[j][k] = (b.getBinStart(k)+b.getBinEnd(k))/2;
                    vs[j][k] = (double)b.getNumValues(k)/readsNum;
                }
            }
            String outfileS = new File (outDirS, fs[i].getName().replaceFirst("_sam.txt.gz", ".pdf")).getAbsolutePath();
            String title = "High divergent reads (>0.05) distribution of "+fs[i].getName()+" Vs reference(blue)";
            GenomeTrendPlot gtp = new GenomeTrendPlot (chromLength, cor, vs, title);
            gtp.addData(cor, refValue);
            gtp.setLowValueLim(0);
            gtp.setHighValueLim(0.02);
            gtp.setIfWithDot(false);
            gtp.saveGragh(outfileS);
        }
    }
    
    public void paintHighDivergenceReads () {
        String infoFileS = "M:\\Database\\cassavaReference\\genomeInfo_Cassava.txt";
        String samDirS = "M:\\production\\cassava\\hmp\\speciesIntrogression\\sam_25\\";
        String outDirS = "M:\\production\\cassava\\hmp\\speciesIntrogression\\divergentReadsPainting\\";
        int binLength = 500000;
        Table t = new Table (infoFileS);
        int[] chromLength = new int[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) chromLength[i] = t.getIntValue(i, 1);
        File[] fs = new File(samDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            ArrayList<Range> rList = new ArrayList();
            ShortreadPEAlignment sa = new ShortreadPEAlignment();
            sa.readFromBWAMEM(fs[i].getAbsolutePath());
            for (int j = 0; j < sa.getAlignmentNumber(); j++) {
                if (sa.isMatchF(j)) {
                    double d = (double)sa.getEditDistanceF(j)/sa.getMatchNumberF(j);
                if (d > 0.05) rList.add(new Range(Integer.valueOf(sa.getHitF(j)), sa.getStartPosF(j), sa.getStartPosF(j)+1));
                }
                if (sa.isMatchB(j)) {
                    double d = (double)sa.getEditDistanceB(j)/sa.getMatchNumberB(j);
                    if (d > 0.05) rList.add(new Range(Integer.valueOf(sa.getHitB(j)), sa.getStartPosB(j), sa.getStartPosB(j)+1));
                }
            }
            Ranges r = new Ranges(rList, "highDivergentPosition");
            r.sortByStartPosition();
            int chromNum = r.getChromosomeNumber();
            int[][] cor = new int[chromNum][];
            double[][] vs = new double[chromNum][];
            for (int j = 0; j < chromNum; j++) {
                Ranges crs = r.getRangesByChromosome(j+1);
                int[] pos = new int[crs.getRangeNumber()];
                double[] value = new double[crs.getRangeNumber()];
                for (int k = 0; k < pos.length; k++) {
                    pos[k] = crs.getRangeStart(k);
                    value[k] = 1;
                }
                Bins b = new Bins (1, chromLength[j], binLength, pos, value);
                cor[j] = new int[b.getBinNum()];
                vs[j] = new double[b.getBinNum()];
                for (int k = 0; k < b.getBinNum(); k++) {
                    cor[j][k] = (b.getBinStart(k)+b.getBinEnd(k))/2;
                    vs[j][k] = b.getNumValues(k);
                }
            }
            String outfileS = new File (outDirS, fs[i].getName().replaceFirst("_sam.txt.gz", ".pdf")).getAbsolutePath();
            String title = "High divergent reads painting from "+fs[i].getName();
            GenomeTrendPlot gtp = new GenomeTrendPlot (chromLength, cor, vs, title);
            gtp.setLowValueLim(0);
            gtp.saveGragh(outfileS);
            
        }
    }
    
    public void editDistance () {
        //use Aligner.editDistanceDist()
    }
    
    public void getErrorFreeRef () {
        String r1FileS = "Q:\\D\\Cassava\\JGI\\FASTQ\\DRJL068_NoIndex_L002_R1_001.fastq.gz";
        String r2FileS = "Q:\\D\\Cassava\\JGI\\FASTQ\\DRJL068_NoIndex_L002_R2_001.fastq.gz";
        String slicedFastqDirS = "M:\\production\\cassava\\hmp\\speciesIntrogression\\errorFreeHiSeq_25\\";
        int readNum = 50000;
        int startIndex = 100000;
        int threshold = 25;
        File fastqR1 = new File (r1FileS);
        File fastqR2 = new File (r2FileS);
        ArrayList<String> fastqListR1 = new ArrayList();
        ArrayList<String> fastqListR2 = new ArrayList();
        System.out.println(fastqR1.getAbsolutePath());
        System.out.println(fastqR2.getAbsolutePath());
        try {
            int phredScale = GBSUtils.determineQualityScoreBase(fastqR1.toPath());
            BufferedReader br1 = IoUtils.getTextGzipReader(fastqR1.getAbsolutePath());
            BufferedReader br2 = IoUtils.getTextGzipReader(fastqR2.getAbsolutePath());
            String temp1 = null;
            String temp2 = null;
            for (int j = 0; j < startIndex; j++) {
                br1.readLine();br1.readLine();br1.readLine();br1.readLine();
                br2.readLine();br2.readLine();br2.readLine();br2.readLine();
            }
            int cnt = 0;
            int count = 0;
            while ((temp1 = br1.readLine())!= null) {
                String header1 = temp1;
                String seq1 = br1.readLine();
                String des1 = br1.readLine();
                String qualS1 = br1.readLine();
                String header2 = br2.readLine();
                String seq2 = br2.readLine();
                String des2 = br2.readLine();
                String qualS2 = br2.readLine();
                byte[] qb1 = qualS1.getBytes();
                byte[] qb2 = qualS2.getBytes();
                if (ifPassThresh(qb1, phredScale, threshold) && ifPassThresh(qb2, phredScale, threshold)) {
                    fastqListR1.add(header1); fastqListR1.add(seq1); fastqListR1.add(des1); fastqListR1.add(qualS1);
                    fastqListR2.add(header2); fastqListR2.add(seq2); fastqListR2.add(des2); fastqListR2.add(qualS2);
                    cnt++;
                }
                if (cnt==readNum) {
                    br1.close();
                    br2.close();
                    break;
                }
                count++;
                if (count%100000 == 0) System.out.println("Read " + String.valueOf(count) + " read");
            }
            String outfile1 = new File (slicedFastqDirS, "refererence_R1.fastq.gz").getAbsolutePath();
            String outfile2 = new File (slicedFastqDirS, "refererence_R2.fastq.gz").getAbsolutePath();
            BufferedWriter bw = IoUtils.getTextGzipWriter(outfile1);
            for (int j = 0; j < fastqListR1.size(); j++) {
                bw.write(fastqListR1.get(j));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            bw = IoUtils.getTextGzipWriter(outfile2);
            for (int j = 0; j < fastqListR2.size(); j++) {
                bw.write(fastqListR2.get(j));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void getErrorFreeHiSeq () {
        //String fastqDirS = "R:\\WGS\\Cassava\\";
        //String slicedFastqDirS = "M:\\production\\cassava\\hmp\\speciesIntrogression\\errorFreeHiSeq_15\\";
        String fastqDirS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\raw\\";
        String slicedFastqDirS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\errorFreeHiSeq_30\\";
        File[] fAll = IoUtils.listRecursiveFiles(new File(fastqDirS));
        File[] fastqR1 = IoUtils.listFilesEndsWith(fAll, "R1.fastq.gz");
        File[] fastqR2 = IoUtils.listFilesEndsWith(fAll, "R2.fastq.gz");
        int readNum = 200000;
        int startIndex = 10000;
        int threshold = 30;
        for (int i = 0; i < fastqR1.length; i++) {
            ArrayList<String> fastqListR1 = new ArrayList();
            ArrayList<String> fastqListR2 = new ArrayList();
            System.out.println(fastqR1[i].getAbsolutePath());
            System.out.println(fastqR2[i].getAbsolutePath());
            
            try {
                int phredScale = GBSUtils.determineQualityScoreBase(fastqR1[i].toPath());
                BufferedReader br1 = IoUtils.getTextGzipReader(fastqR1[i].getAbsolutePath());
                BufferedReader br2 = IoUtils.getTextGzipReader(fastqR2[i].getAbsolutePath());
                String temp1 = br1.readLine(); temp1 = br1.readLine();
                String temp2 = br2.readLine(); temp2 = br2.readLine();
                if (temp1.length() < 170) {  
                    br1.close(); br2.close();
                    System.out.println("This is NextSeq, skipped");
                    continue;
                }
                else {
                    System.out.println("This is HiSeq, looking for error free reads");
                    br1.readLine(); br1.readLine();                    
                    br2.readLine(); br2.readLine();                    
                }
                for (int j = 1; j < startIndex; j++) {
                    br1.readLine();br1.readLine();br1.readLine();br1.readLine();
                    br2.readLine();br2.readLine();br2.readLine();br2.readLine();
                }
                int cnt = 0;
                int count = 0;
                while ((temp1 = br1.readLine())!= null) {
                    String header1 = temp1;
                    String seq1 = br1.readLine();
                    String des1 = br1.readLine();
                    String qualS1 = br1.readLine();
                    String header2 = br2.readLine();
                    String seq2 = br2.readLine();
                    String des2 = br2.readLine();
                    String qualS2 = br2.readLine();
                    byte[] qb1 = qualS1.getBytes();
                    byte[] qb2 = qualS2.getBytes();
                    if (ifPassThresh(qb1, phredScale, threshold) && ifPassThresh(qb2, phredScale, threshold)) {
                        fastqListR1.add(header1); fastqListR1.add(seq1); fastqListR1.add(des1); fastqListR1.add(qualS1);
                        fastqListR2.add(header2); fastqListR2.add(seq2); fastqListR2.add(des2); fastqListR2.add(qualS2);
                        cnt++;
                    }
                    if (cnt==readNum) {
                        br1.close();
                        br2.close();
                        break;
                    }
                    count++;
                    if (count%100000 == 0) {
                        System.out.println("Read through " + String.valueOf(count) + " read, " + String.valueOf(cnt) + " high quality, " + String.format("%.2f", (double)cnt/count));
                    }
                }
                if (cnt < readNum) {
                    br1.close();br2.close();
                    System.out.println("This is a bad lane, doesn't have enough error free reads, skipped");
                    continue;
                } 
                String outfile1 = new File (slicedFastqDirS, fastqR1[i].getName()).getAbsolutePath();
                String outfile2 = new File (slicedFastqDirS, fastqR2[i].getName()).getAbsolutePath();
                BufferedWriter bw = IoUtils.getTextGzipWriter(outfile1);
                for (int j = 0; j < fastqListR1.size(); j++) {
                   bw.write(fastqListR1.get(j));
                   bw.newLine();
                }
                bw.flush();
                bw.close();
                bw = IoUtils.getTextGzipWriter(outfile2);
                for (int j = 0; j < fastqListR2.size(); j++) {
                   bw.write(fastqListR2.get(j));
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
    
    private boolean ifPassThresh (byte[] qb, int phredScale, int threshold) {
        for (int i = qb.length-1; i > 0; i--) {
            if ((qb[i]-phredScale) < threshold) return false;
        }
        return true;
    }
}
