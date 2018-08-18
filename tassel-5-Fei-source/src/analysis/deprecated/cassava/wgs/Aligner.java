/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.cassava.wgs;

import format.ShortreadAlignment;
import format.ShortreadPEAlignment;
import gnu.trove.list.array.TDoubleArrayList;
import graphcis.r.CumulativeDistribution;
import graphcis.r.DensityPlot;
import graphcis.r.Histogram;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import net.maizegenetics.dna.read.FastqChunk;
import net.maizegenetics.dna.read.PEFastqChunk;
import net.maizegenetics.dna.read.ReadUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class Aligner {
    
    public Aligner () {
        this.subsetFastq();
        //this.convertToFasta();
        //this.contigingFastq();
        //this.compareAligner();
        //this.editDistanceDist();
        //this.editDistanceStats();
    }
    
    public void editDistanceStats () {
        String fastqDirS = "M:\\pipelineTest\\cassava\\wgs\\pe\\fastqChunk\\";
        String samDirS = "M:\\pipelineTest\\cassava\\wgs\\pe\\alignmentBWAMEM\\sam\\";
        String outfileS = "M:\\production\\cassava\\hmp\\aligner\\editDistance\\distanceCorrelation.txt";
        File[] sams = new File(samDirS).listFiles();
        String[] sampleNames = new String[sams.length];
        for (int i = 0; i < sampleNames.length; i++) {
            sampleNames[i] = sams[i].getName().replaceFirst("_sam.txt.gz", "");
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("SampleName");
            int n = sampleNames[0].split("_").length;
            for (int i = 0; i < n; i++) {
                bw.write("\tfield_"+String.valueOf(i+1));
            }
            bw.write("\tPEMapRatio\tSEMapRatio\tNoMapRatio\tMQ_mean\tMQ_sd\tMQ_1stQ\tMQ_2ndQ\tMQ_3rdQ\tED_mean\tED_sd\tED_skewness\tED_kurtosis\tED_1stQ\tED_2ndQ\tED_3rdQ");
            bw.newLine();
            for (int i = 0; i < sampleNames.length; i++) {
            //for (int i = 0; i < 2; i++) {
                bw.write(sampleNames[i]);
                String[] temp = sampleNames[i].split("_");
                for (int j = 0; j < temp.length; j++) {
                    bw.write("\t"+temp[j]);
                }
                ShortreadPEAlignment spe = new ShortreadPEAlignment ();
                spe.readFromBWAMEM(sams[i].getAbsolutePath());
                bw.write("\t"+String.valueOf(spe.getPEMappedRatio())+"\t"+String.valueOf(spe.getSEMapptedRatio())+"\t"+String.valueOf(spe.getNoMapptedRatio()));
                TDoubleArrayList vList = new TDoubleArrayList();
                for (int j = 0; j < spe.getAlignmentNumber(); j++) {
                    if (spe.isMatchF(j)) vList.add(spe.getMappingQualityF(j));
                    if (spe.isMatchB(j)) vList.add(spe.getMappingQualityB(j));
                }
                double[] v = vList.toArray();
                Arrays.sort(v);
                DescriptiveStatistics des = new DescriptiveStatistics(v);
                bw.write("\t"+des.getMean()+"\t"+des.getStandardDeviation());
                bw.write("\t"+v[(int)(v.length*0.25)]+"\t"+v[(int)(v.length*0.50)]+"\t"+v[(int)(v.length*0.75)]);
                vList.clear();
                for (int j = 0; j < spe.getAlignmentNumber(); j++) {
                    float va1 = spe.getEditDistanceRatioF(j);
                    float va2 = spe.getEditDistanceRatioB(j);
                    if (!Float.isNaN(va1) && !Float.isNaN(va2)) vList.add((va1+va2)/2);
                }
                v = vList.toArray();
                Arrays.sort(v);
                des = new DescriptiveStatistics(v);
                bw.write("\t"+des.getMean()+"\t"+des.getStandardDeviation()+"\t"+des.getSkewness()+"\t"+des.getKurtosis());
                bw.write("\t"+v[(int)(v.length*0.25)]+"\t"+v[(int)(v.length*0.50)]+"\t"+v[(int)(v.length*0.75)]);
                
                
                
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void editDistanceDist () {
        //String infileDirS = "M:\\pipelineTest\\cassava\\wgs\\pe\\alignmentBWAMEM\\sam\\";
        //String distanceFigureDirS = "M:\\pipelineTest\\cassava\\wgs\\aligner\\editDistance\\pdf";
        
        String infileDirS = "M:\\production\\cassava\\hmp\\speciesIntrogression\\sam_30\\";
        String distanceFigureDirS = "M:\\production\\cassava\\hmp\\speciesIntrogression\\pdf_30\\";
        String nonCassavaReadsReportFileS = "M:\\production\\cassava\\hmp\\speciesIntrogression\\nonCassavaReads.txt";
        
        //String infileDirS = "M:\\production\\cassava\\hmp\\alignMultiGenomes\\sam\\";
        //String distanceFigureDirS = "M:\\production\\cassava\\hmp\\alignMultiGenomes\\pdf\\";
        //String nonCassavaReadsReportFileS = "M:\\production\\cassava\\hmp\\alignMultiGenomes\\nonCassavaReads.txt";
        
        File[] fs = new File(infileDirS).listFiles();
        int fileNumber = 1;
        if (fs.length < fileNumber) fileNumber = fs.length;
        int[] readsNum = new int[fs.length];
        int[] divergent5 = new int[fs.length];
        int[] noMap = new int[fs.length];
        for (int i = 0; i < fs.length; i++) {
            TDoubleArrayList distList = new TDoubleArrayList();
            String fileName = fs[i].getAbsolutePath();
            try {
                BufferedReader br = IoUtils.getTextGzipReader(fileName);
                while(br.readLine().startsWith("@PG")==false) {}
                String temp;
                String[] forwardArray = null;
                String[] backwardArray = null;
                int cnt = 0;
                int cnt5 = 0;
                int cntNoMap = 0;
                ArrayList<String[]> infoL = new ArrayList();
                infoL.add(br.readLine().split("\\s"));
                while((temp = br.readLine())!=null) {
                    String[] tem = temp.split("\\s");
                    if (tem[0].equals(infoL.get(infoL.size()-1)[0])) {
                        infoL.add(tem);
                    }
                    else {
                        Iterator<String[]> it = infoL.iterator();
                        while (it.hasNext()) {
                            String[] te = it.next();
                            int flag = Integer.valueOf(te[1]);
                            //the 11th(index) bit is essentially a bug of bwa-mem, the 11th bit is not defined in SAM protocol
                            if (getBitFromInt(flag, 11) == 1 || getBitFromInt(flag, 8) == 1 || getBitFromInt(flag, 10) == 1) it.remove();
                        }
                        if (infoL.size() != 2) {
                            System.out.println(infoL.get(0)[0]);
                        }
                        if (getBitFromInt(Integer.valueOf(infoL.get(0)[1]), 6) == 1) {
                            forwardArray = infoL.get(0);
                            backwardArray = infoL.get(1);
                        }
                        else {
                            forwardArray = infoL.get(1);
                            backwardArray = infoL.get(0);
                        }
                        if (!forwardArray[2].equals("*") && !backwardArray[2].equals("*")) {
                            double d1 = this.getMatchNumberFromCigar(forwardArray[5]);
                            double d2 = this.getMatchNumberFromCigar(backwardArray[5]);
                            if (d1 != Integer.MIN_VALUE && d2 != Integer.MIN_VALUE) {
                                d1 = (double)this.getDistance(forwardArray)/d1;
                                d2 = (double)this.getDistance(backwardArray)/d2;
                                double d = (d1+d2)/2;
                                //turn on controling mapping quality
                                //if (Integer.valueOf(forwardArray[4])>29 && Integer.valueOf(backwardArray[4])>29) distList.add(d);
                                distList.add(d);
                                if (d > 0.05) cnt5++;
                            }
                        }
                        else {
                            cntNoMap++;
                        }
                        if (cnt%500000 == 0) System.out.println("Read in " + String.valueOf(cnt) + " lines");
                        cnt++;
                        infoL.clear();
                        infoL.add(tem);
                    }
                }
                readsNum[i] = cnt;
                divergent5[i] = cnt5;
                noMap[i] = cntNoMap;
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            double[] distance = distList.toArray();
            String outfileS = fs[i].getName().replaceFirst(".gz", ".pdf");
            outfileS = new File(distanceFigureDirS, outfileS).getAbsolutePath();
            DensityPlot d = new DensityPlot(distance);
            d.setTitle(fs[i].getName());   
//            d.setYLim(0, 40);
//            d.setXLim(0, 0.30);
            d.setXLab("Distance");
            d.setYLab("Density");
            d.saveGraph(outfileS);
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(nonCassavaReadsReportFileS);
            bw.write("Files\tRatio_EDGreaterThan5\tDontMap\tEDGreaterThan5AndDontMap");
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                double d1 = (double)divergent5[i]/readsNum[i];
                double d2 = (double)noMap[i]/readsNum[i];
                double d3 = (double)(divergent5[i]+noMap[i])/readsNum[i];
                bw.write(fs[i].getName()+"\t"+String.valueOf(d1)+"\t"+String.valueOf(d2)+"\t"+String.valueOf(d3));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    private int getMatchNumberFromCigar (String cigar) {
        int cnt = 0;
        if (cigar.contains("M")) {
            String[] temp = cigar.split("M");
            int n;
            if (cigar.endsWith("M")) n = temp.length;
            else n = temp.length-1;
            for (int i = 0; i < n; i++) {
                String[] tem = temp[i].split("\\D");
                cnt+=Integer.valueOf(tem[tem.length-1]);
            }
            return cnt;
        }
        return Integer.MIN_VALUE;
    }
    
    private int getDistance (String temp[]) {
        
        for (int i = 11; i < temp.length; i++) {
            if (temp[i].startsWith("NM")) {
                return Integer.valueOf(temp[i].split(":")[2]);
            }
        }
        return Integer.MIN_VALUE;
    }
    
    private int getBitFromInt (int n, int k) {
        return (n >> k) & 1;
    }
 
    public void compareAligner () {
        String fastaR1FileS = "M:\\production\\cassava\\hmp\\aligner\\sam\\trim\\R1_blast.txt";
        String fastaR2FileS = "M:\\production\\cassava\\hmp\\aligner\\sam\\trim\\R2_blast.txt";
        String bwaSamFileS = "M:\\production\\cassava\\hmp\\aligner\\sam\\trim\\bwa_default.sam";
        String bowtie2SamFileS = "M:\\production\\cassava\\hmp\\aligner\\sam\\trim\\bowtie2_sensitive.sam";
        String peBlastFileS = "M:\\production\\cassava\\hmp\\aligner\\comparison\\blast_pe.txt";
        String pebwaFileS = "M:\\production\\cassava\\hmp\\aligner\\comparison\\bwa_pe.txt";
        String pebowtie2FileS = "M:\\production\\cassava\\hmp\\aligner\\comparison\\bowtie2_pe.txt";
        String proportionFileS = "M:\\production\\cassava\\hmp\\aligner\\comparison\\alignProportion.txt";
        String concordantRatioFileS = "M:\\production\\cassava\\hmp\\aligner\\comparison\\concordant.txt";
        ShortreadPEAlignment blast = new ShortreadPEAlignment();
        blast.readFromBlast(fastaR1FileS, fastaR2FileS, 1e-5, 150);
        ShortreadPEAlignment bwa = new ShortreadPEAlignment();
        bwa.readFromBWAMEM(bwaSamFileS);
        ShortreadPEAlignment bowtie2 = new ShortreadPEAlignment();
        bowtie2.readFromBowtie2(bowtie2SamFileS);
        blast.writeShortreadPEAlignment(peBlastFileS, IOFileFormat.Text);
        bwa.writeShortreadPEAlignment(pebwaFileS, IOFileFormat.Text);
        bowtie2.writeShortreadPEAlignment(pebowtie2FileS, IOFileFormat.Text);
        int[] peAlign = new int[3];
        int[] seAlign = new int[3];
        int[] noAlign = new int[3];
        String[] aligner = new String[3];
        int[] total = new int[3];
        double[] concordantRatio = new double[3];
        int conSize = 20000;
        for (int i = 0; i < bwa.getAlignmentNumber(); i++) {
            aligner[0] = "BWA-MEM";
            total[0] = bwa.getAlignmentNumber();
            if (bwa.isMatchBothEnds(i)) {
                peAlign[0]++;
                if (bwa.getHitF(i).equals(bwa.getHitB(i)) && Math.abs(bwa.getStartPosF(i)-bwa.getEndPosF(i)) < conSize) concordantRatio[0]++;
            }
            else if (bwa.isMatchB(i) || bwa.isMatchF(i)) {
                seAlign[0]++;
            }
            else {
                noAlign[0]++;
            }
        }
        for (int i = 0; i < bowtie2.getAlignmentNumber(); i++) {
            aligner[1] = "Bowtie2";
            total[1] = bowtie2.getAlignmentNumber();
            if (bowtie2.isMatchBothEnds(i)) {
                peAlign[1]++;
                if (bowtie2.getHitF(i).equals(bowtie2.getHitB(i)) && Math.abs(bowtie2.getStartPosF(i)-bowtie2.getEndPosF(i)) < conSize) concordantRatio[1]++;
            }
            else if (bowtie2.isMatchB(i) || bowtie2.isMatchF(i)) {
                seAlign[1]++;
            }
            else {
                noAlign[1]++;
            }
        }
        for (int i = 0; i < blast.getAlignmentNumber(); i++) {
            aligner[2] = "Blast";
            total[2] = blast.getAlignmentNumber();
            if (blast.isMatchBothEnds(i)) {
                peAlign[2]++;
                if (blast.getHitF(i).equals(blast.getHitB(i)) && Math.abs(blast.getStartPosF(i)-blast.getEndPosF(i)) < conSize) concordantRatio[2]++;
            }
            else if (blast.isMatchB(i) || blast.isMatchF(i)) {
                seAlign[2]++;
            }
            else {
                noAlign[2]++;
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(proportionFileS);
            bw.write("Aligner\tPEAligned\tSEAligned\tNoAligned\tTotalNumber\tPEConcordantRatio");
            bw.newLine();
            for (int i = 0; i < aligner.length; i++) {
                concordantRatio[i] = concordantRatio[i]/peAlign[i];
                bw.write(aligner[i]+"\t"+peAlign[i]+"\t"+seAlign[i]+"\t"+noAlign[i]+"\t"+total[i]+"\t"+concordantRatio[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        TDoubleArrayList[] dList = new TDoubleArrayList[2];
        double[][] distance = new double[2][];
        int diChrDistance = 1000000;
        for (int i = 0; i < dList.length; i++) dList[i] = new TDoubleArrayList();
        for (int i = 0; i < blast.getAlignmentNumber(); i++) {
            if (!blast.isMatchF(i)) continue;
            if (bwa.isMatchF(i)) {
                if (bwa.getHitF(i).equals(blast.getHitF(i))) {
                    dList[0].add(Math.abs(bwa.getStartPosF(i)-blast.getStartPosF(i)));
                }
                else {
                    dList[0].add(Math.abs(bwa.getStartPosF(i)-blast.getStartPosF(i))+diChrDistance);
                }
            }
            if (bowtie2.isMatchF(i)) {
                if (bowtie2.getHitF(i).equals(blast.getHitF(i))) {
                    dList[1].add(Math.abs(bowtie2.getStartPosF(i)-blast.getStartPosF(i)));
                }
                else {
                    dList[1].add(Math.abs(bowtie2.getStartPosF(i)-blast.getStartPosF(i))+diChrDistance);
                }
            }
        }
        for (int i = 0; i < dList.length; i++) distance[i] = dList[i].toArray();
        int[] sameCnt = new int[2];
        for (int i = 0; i < sameCnt.length; i++) {
            for (int j = 0; j < distance[i].length; j++)  {
                if (distance[i][j] < 20) sameCnt[i]++;
            }
        }
        System.out.println("Same alignment ratio for bwa-mem is " + String.valueOf((double)sameCnt[0]/distance[0].length));
        System.out.println("Same alignment ratio for bowtie2 is " + String.valueOf((double)sameCnt[1]/distance[1].length));
    }
    
    public void contigingFastq () {
        String infileDirS = "M:\\production\\cassava\\hmp\\aligner\\sampleFastq\\";
        File[] fs = new File(infileDirS).listFiles();
        PEFastqChunk pf = new PEFastqChunk (fs[0].getAbsolutePath(), fs[1].getAbsolutePath(), ReadUtils.ReadFormat.FastqGzip);
        pf.merge(true);
    }
    
    public void convertToFasta () {
        //String infileDirS = "M:\\production\\cassava\\hmp\\aligner\\sampleFastq\\";
        //String outfileDirS = "M:\\production\\cassava\\hmp\\aligner\\sampleFasta\\";
        String infileDirS = "M:\\production\\cassava\\hmp\\aligner\\sampleFastq\\trim\\";
        String outfileDirS = "M:\\production\\cassava\\hmp\\aligner\\sampleFasta\\trim\\";
        File[] fs = new File(infileDirS).listFiles();
        
        for (int i = 0; i < fs.length; i++) {
            //String outfileS = fs[i].getName().replaceFirst("/fastq.gz", ".fa");
            String outfileS = fs[i].getName().replaceFirst(".fq.gz", ".fa");
            FastqChunk fc = new FastqChunk(fs[i].getAbsolutePath(), ReadUtils.ReadFormat.FastqGzip);
            fc.writeFasta(new File(outfileDirS, outfileS).getAbsolutePath());
        }
    }
    
    public void subsetFastq () {
        int size = 100000;
        int start = 100000;
        int fileNumber = 2;
        String infileDirS = "N:\\cassavaWGS\\HF2J5BGXX\\";
        String outfileDirS = "M:\\production\\cassava\\hmp\\aligner\\sampleFastq\\";
        //String infileDirS = "M:\\test\\source\\";
        //String outfileDirS = "M:\\test\\testFastq\\";
        File[] fs = IoUtils.listFilesEndsWith(new File(infileDirS).listFiles(), "fastq.gz");
        for (int i = 0; i < fileNumber; i++) {
            FastqChunk f = new FastqChunk (fs[i].getAbsolutePath(), ReadUtils.ReadFormat.FastqGzip, start, size);
            File nf = new File (outfileDirS, fs[i].getName());
            f.writeFastq(nf.getAbsolutePath(), ReadUtils.ReadFormat.FastqGzip);
        }
    }
    
    
}
