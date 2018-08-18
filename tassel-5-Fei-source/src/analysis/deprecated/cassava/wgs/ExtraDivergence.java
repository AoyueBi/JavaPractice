/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.cassava.wgs;

import format.BlastAlignment;
import format.Fasta;
import format.GeneFeature;
import format.JPileUp;
import format.Range;
import format.Ranges;
import format.Sequence;
import format.ShortreadAlignment;
import format.ShortreadPEAlignment;
import format.Table;
import gnu.trove.list.array.TDoubleArrayList;
import graphcis.r.CumulativeDistribution;
import graphcis.r.DensityPlot;
import graphcis.r.DensityPlotMultiClass;
import graphcis.r.Histogram;
import graphcis.r.ScatterPlot;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import net.maizegenetics.dna.read.FastqChunk;
import net.maizegenetics.dna.read.Read;
import net.maizegenetics.dna.read.ReadUtils;
import net.maizegenetics.dna.read.ReadUtils.ReadFormat;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import utils.FArrayUtils;
import utils.FStringUtils;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class ExtraDivergence {
    
    public ExtraDivergence () {
        /**check chloroplast*/
        //this.chloroplastPipe();
        /**check if there is any contaminaton from other species*/
        //this.contaminationPipe();
        /**check edit distance in high or low depth regions*/
        //this.depthPipe();
        /**look for divergent contig*/
        //this.assemblyPipe();
        /**look for divergent genes*/
        //this.genePipe();
        /**check if there is real diversity*/
        //this.reallDiversityPipe();
        //this.repeatFilterPipe();
    }
    
    public void repeatFilterPipe () {
        //this.creatReadsOnChr1();
        //this.creatReadsFasta();
        //this.creatDivergentReadsByDepth();
        this.creatRepeatFilterStatistics();
    }
    
    private void creatRepeatFilterStatistics () {
        String infileDirS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\repeatFilter\\alignment\\";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\repeatFilter\\removedPortion.txt";
        String[] subDirS =  new File(infileDirS).list();
        System.out.println(subDirS[0]);
        double eThresh = 0.01;
        double[][] removedPortion = new double[subDirS.length][];
        int[] wordSize = null;
        for (int i = 0; i < subDirS.length; i++) {
            String dirS = new File(infileDirS, subDirS[i]).getAbsolutePath();
            File[] fs = new File(dirS).listFiles();
            removedPortion[i] = new double[fs.length];
            wordSize = new int[fs.length];
            for (int j = 0; j < fs.length; j++) {
                BlastAlignment ba = new BlastAlignment();
                ba.readFromBlastTable(fs[j].getAbsolutePath(), eThresh);
                String[] query = ba.getQuerys();
                int cnt = 0;
                for (int k = 0; k < query.length; k++) {
                    if (!ba.isNoMatch(query[k])) cnt++;
                }
                removedPortion[i][j] = (double)cnt/query.length;
                wordSize[j] = Integer.valueOf(fs[j].getName().split("_")[1].replaceFirst("word", ""));
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("WordSize");
            for (int i = 0; i < subDirS.length; i++) {
                bw.write("\t"+subDirS[i]);
            }
            bw.newLine();
            for (int i = 0; i < wordSize.length; i++) {
                bw.write(String.valueOf(wordSize[i]));
                for (int j = 0; j < removedPortion.length; j++) {
                    bw.write("\t"+removedPortion[j][i]);
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
    
    private void creatDivergentReadsByDepth () {
        String jPileupFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\jPileUp\\I000070_chr_01_MQ30.jPileUp.bin";
        String infoFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\repeatFilter\\I000070\\chr1_readList.txt";
        String fastaFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\repeatFilter\\I000070\\chr1_reads.fa";
        String lowDepthDivergentFastaFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\repeatFilter\\I000070\\divergent_reads_0_10.fa";
        String highDepthDivergentFastaFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\repeatFilter\\I000070\\divergent_reads_60_unlimited.fa";
        double divergentThresh = 0.05;
        int lowDepthThresh = 10;
        int highDepthThresh = 60;
        JPileUp jp = new JPileUp (jPileupFileS, IOFileFormat.Binary);
        Fasta fa = new Fasta (fastaFileS);
        boolean[] ifOutLow = new boolean[fa.getSeqNumber()];
        boolean[] ifOutHigh = new boolean[fa.getSeqNumber()];
        Table t = new Table (infoFileS);
        HashMap<String, Double> seqEDMap = new HashMap();
        HashMap<String, Integer> seqPosStartMap = new HashMap();
        HashMap<String, Integer> seqPosEndMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (!t.content[i][1].equals("1")) continue;
            seqEDMap.put(t.content[i][0], t.getDoubleValue(i, 4));
            seqPosStartMap.put(t.content[i][0], Integer.valueOf(t.content[i][2]));
            seqPosEndMap.put(t.content[i][0], Integer.valueOf(t.content[i][3]));
        }
        for (int i = 0; i < fa.getSeqNumber(); i++) {
            if (seqEDMap.get(fa.getName(i)) < divergentThresh) continue;
            int start = seqPosStartMap.get(fa.getName(i));
            int end = seqPosEndMap.get(fa.getName(i));
            int cnt = 0;
            double depth = 0;
            for (int j = 0; j < end-start+1; j++) {
                int pos = j+start;
                int index = jp.getSiteIndex(1, pos);
                if (index < 0) continue;
                cnt++;
                depth+=jp.getSiteDepth(index);
            }
            if (cnt == 0) continue;
            depth = depth/cnt;
            if (depth < lowDepthThresh) ifOutLow[i] = true;
            if (depth > highDepthThresh) ifOutHigh[i] = true;
        }
        fa.writeFasta(lowDepthDivergentFastaFileS, ifOutLow);
        fa.writeFasta(highDepthDivergentFastaFileS, ifOutHigh);
    }
    
    private void creatReadsFasta () {
        String infoFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\repeatFilter\\I000070\\chr1_readList.txt";
        String fastqFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\errorFreeHiSeq_30\\3423_7332_10699_HAK9RADXX_I000070_500_GTCCGC_R1.fastq.gz";
        String fastaFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\repeatFilter\\I000070\\chr1_reads.fa";
        Table t = new Table (infoFileS);
        ArrayList<String> nameList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (!t.content[i][1].equals("1")) continue;
            nameList.add(t.content[i][0]);
        }
        String[] name = nameList.toArray(new String[nameList.size()]);
        Arrays.sort(name);
        FastqChunk fq = new FastqChunk(fastqFileS, ReadFormat.FastqGzip);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(fastaFileS);
            for (int i = 0; i < fq.getReadNum(); i++) {
                Read r = fq.getRead(i);
                String query = r.getID().split(" ")[0].replaceFirst("@", "");
                int index = Arrays.binarySearch(name, query);
                if (index < 0) continue;
                bw.write(">"+query);
                bw.newLine();
                bw.write(r.getSeq());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    private void creatReadsOnChr1 () {
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\sam_30\\I000070.sam.gz";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\repeatFilter\\I000070\\chr1_readList.txt";
        ShortreadPEAlignment sa = new ShortreadPEAlignment();
        sa.readFromBWAMEM(infileS);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Name\tR1Chr\tR1PosStart\tR1PosEnd\tR1ED");
            bw.newLine();
            for (int i = 0; i < sa.getAlignmentNumber(); i++) {
                StringBuilder sb = new StringBuilder();
                if (sa.getHitF(i).equals("1")||sa.getHitB(i).equals("1")) {
                    sb.append(sa.getQuery(i)).append("\t").append(sa.getHitF(i)).append("\t").append(sa.getStartPosF(i)).append("\t").append(sa.getEndPosF(i)).append("\t");
                    sb.append(sa.getEditDistanceRatioF(i));
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void reallDiversityPipe () {
        //this.generateKu50Fastq();
        //this.plotKu50Divergenece();
        //this.generateAM560Fastq();
        //this.plotAm560Divergenece();
        //this.getLowDepthSequence();
        //this.getLowDepthSequenceIllunimaMapped();
        //this.getWordCount();
        //this.getWordCount2();
        //this.getTomatoLikeRef();
        //this.estimateProportionOfContamination();
    }
    
    public void estimateProportionOfContamination () {
        String infileS1 = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\readDepth_max100.txt";
        String infileS2 = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\AM560_diversity\\AM560_max100.depth.txt";
        //String infileS1 = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\readDepth.txt";
        //String infileS2 = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\AM560_diversity\\AM560.depth.txt";
        String pdfFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\AM560_diversity\\AM560.depth.compare.pdf";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\AM560_diversity\\AM560.depth.compare.txt";
        int size = 8000;
        double[] d1 = new double[size];
        double[] d2 = new double[size];
        Table t = new Table (infileS1);
        double[] v = new double[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) v[i] = t.getDoubleValue(i, 0);
        d1 = FArrayUtils.getRandomSubset(v, size);
        t = new Table (infileS2);
        v = new double[t.getRowNumber()];
        for (int i = 0; i < size; i++) v[i] = t.getDoubleValue(i, 0);
        d2 = FArrayUtils.getRandomSubset(v, size);
        double[][] d = new double[2][];
        d[0] = d1;
        d[1] = d2;
        String[] variableNames = {"Depth of I000070", "Depth of I000070 where reference fastqs align"};
        DensityPlotMultiClass den = new DensityPlotMultiClass(d, variableNames);
        den.setXLim(0, 60);
        den.setXLab("Read depth");
        den.setYLab("Density");
        den.saveGraph(pdfFileS);
        int intervalSize = 100;
        double[] depthInterval = new double[intervalSize];
        for (int i = 0; i < depthInterval.length; i++) {
            depthInterval[i] = i;
        }
        double[] dp1 = new double[intervalSize];
        double[] dp2 = new double[intervalSize];
        for (int i = 0; i < d1.length; i++) {
            int index = Arrays.binarySearch(depthInterval, d1[i]);
            if (index < 0) index = -index-1;
            if (index > intervalSize - 1) index = intervalSize - 1;
            dp1[index]++;
        }
        for (int i = 0; i < d2.length; i++) {
            int index = Arrays.binarySearch(depthInterval, d2[i]);
            if (index < 0) index = -index-1;
            if (index > intervalSize - 1) index = intervalSize - 1;
            dp2[index]++;
        }
        for (int i = 0; i < intervalSize; i++) {
            dp1[i] = dp1[i]/size;
            dp2[i] = dp2[i]/size;
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("DepthInterval\tDepthOfI000070\tDepthOfI000070ReferenceAlign");
            bw.newLine();
            for (int i = 0; i < intervalSize; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(depthInterval[i]).append("\t").append(dp1[i]).append("\t").append(dp2[i]);
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
    
    public void getTomatoLikeRef () {
        String fastaFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\lowDepthRefSeq\\lowDepthRef.fa";
        String alignmentFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\lowDepthRefSeq\\lowDepthRef_Alignment.txt";
        String outFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\lowDepthRefSeq\\tomatoLikeRef.fa";
        Fasta f = new Fasta(fastaFileS);
        boolean[] ifOut = new boolean[f.getSeqNumber()];
        f.sortRecordByName();
        BlastAlignment ba = new BlastAlignment ();
        ba.readFromBlast(alignmentFileS, 1e-10);
        String[] querys = ba.getQuerys();
        for (int i = 0; i < querys.length; i++) {
            int index = ba.getAlignmentStartIndexByQuery(querys[i]);
            if (ba.getHit(index).contains("Solanum")) {
                index = f.getIndex(querys[i]);
                ifOut[index] = true;
            }
        }
        f.writeFasta(outFileS, ifOut);
    }
            
    public void getWordCount2 () {
        //String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\lowDepthRefSeq\\lowDepthRef_Alignment.txt";
        //String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\lowDepthRefSeq\\lowDepthRef.wordCnt2.txt";
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\lowDepthRefSeq\\lowDepthRef_Illumina_Alignment.txt";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\lowDepthRefSeq\\lowDepthRef_Illumina.wordCnt2.txt";
        double eThresh = 1e-10;
        ArrayList<String> wordList = new ArrayList();
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            String temp;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith(">")) {
                    String[] tem = temp.split("\\ +");
                    while (!(temp = br.readLine()).startsWith(" Score =")) {}
                    String[] te = temp.split("\\ +");
                    if (Double.valueOf(te[te.length-1]) <= eThresh) {
                        StringBuilder sb = new StringBuilder();
                        for (int i = 1; i < tem.length; i++) {
                            sb.append(tem[i]).append("\t");
                        }
                        wordList.add(sb.toString());
                    }
                }
            }
            br.close();
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            for (int i = 0; i < wordList.size(); i++) {
                bw.write(wordList.get(i));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void getWordCount () {
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\lowDepthRefSeq\\lowDepthRef_Alignment.txt";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\lowDepthRefSeq\\lowDepthRef.wordCnt.txt";
        //String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\lowDepthRefSeq\\lowDepthRef_Illumina_Alignment.txt";
        //String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\lowDepthRefSeq\\lowDepthRef_Illumina.wordCnt.txt";
        double eThresh = 1e-10;
        ArrayList<String> wordList = new ArrayList();
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            String temp;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("Sequences producing")) {
                    temp = br.readLine();
                    while (!(temp = br.readLine()).isEmpty()) {
                       String[] tem = temp.split("\\ +");
                       if (Double.valueOf(tem[tem.length-1])>eThresh) continue;
                       wordList.add(tem[1]);wordList.add(tem[2]);
                    }
                }
            }
            br.close();
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            for (int i = 0; i < wordList.size(); i+=2) {
                bw.write(wordList.get(i)+"\t"+wordList.get(i+1));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void getLowDepthSequenceIllunimaMapped () {
        String jPileupFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\jPileUp\\I000070_chr_01_MQ30.jPileUp.bin";
        String genomeFileS = "M:\\Database\\cassavaReference\\genome\\Manihot esculenta\\cassavaV6_chrAndScaffoldsCombined_numeric.fa";
        String alignmentFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\AM560_diversity\\AM560.sam";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\lowDepthRefSeq\\lowDepthRef_Illumina.fa";
        int length = 200;
        int depthThresh = 5;
        int size = 1000;
        JPileUp jp = new JPileUp (jPileupFileS, IOFileFormat.Binary);
        Fasta f = new Fasta (genomeFileS);
        f.sortRecordByName();
        String chrSeq = f.getSeq(f.getIndex("1"));
        ShortreadAlignment sa = new ShortreadAlignment();
        sa.readFromBWAMEM(alignmentFileS);
        String[] querys = sa.getQuerys();
        ArrayList<String> seqList = new ArrayList();
        for (int i = 0; i < querys.length; i++) {
            int index = sa.getAlignmentStartIndexByQuery(querys[i]);
            int startPos = sa.getStartPos(index);
            int endPos = sa.getEndPos(index);
            if ((endPos-startPos) < length/2) continue;
            int startIndex = jp.getSiteIndex(1, startPos);
            int endIndex = jp.getSiteIndex(1, startPos+length);
            if (startIndex < 0) continue;
            if (endIndex < 0) continue;
            int depthCnt = 0;
            for (int j = startIndex; j < endIndex; j++) {
                depthCnt+=jp.getSiteDepth(j);
            }
            if ((double)depthCnt/(endIndex-startIndex)>depthThresh) continue;
            seqList.add(chrSeq.substring(startPos, startPos+length));
            if (seqList.size() == size) break;
        }
        String[] seqs = seqList.toArray(new String[seqList.size()]);
        String[] names = new String[seqs.length];
        int[] ids = new int[seqs.length];
        for (int i = 0; i < names.length; i++) {
            names[i] = String.valueOf(i);
            ids[i] = i;
        }
        Fasta ff = new Fasta(names, seqs, ids);
        ff.writeFasta(outfileS);
    }
    
    public void getLowDepthSequence () {
        String jPileupFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\jPileUp\\I000070_chr_01_MQ30.jPileUp.bin";
        String genomeFileS = "M:\\Database\\cassavaReference\\genome\\Manihot esculenta\\cassavaV6_chrAndScaffoldsCombined_numeric.fa";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\lowDepthRefSeq\\lowDepthRef.fa";
        int size = 1000;
        int length = 200;
        int depthThresh = 5;
        JPileUp jp = new JPileUp (jPileupFileS, IOFileFormat.Binary);
        Fasta f = new Fasta (genomeFileS);
        f.sortRecordByName();
        String chrSeq = f.getSeq(f.getIndex("1"));
        int cnt = 0;
        ArrayList<String> seqList = new ArrayList();
        while (cnt != size) {
            int index = (int)(jp.getSiteNumber()*Math.random());
            if ((index-200) < 0) continue;
            if (index>jp.getSiteNumber()) continue;
            if (jp.getSiteDepth(index) > depthThresh) continue;
            int pos = jp.getPosition(index);
            int prePos = pos-length/2;
            int postPos = pos+length/2;
            int preIndex = jp.getSiteIndex(1, prePos);
            int postIndex = jp.getSiteIndex(1, postPos);
            if (preIndex < 0) continue;
            if (postIndex < 0) continue;
            int depthCnt = 0;
            for (int j = preIndex; j < postIndex; j++) {
                depthCnt+=jp.getSiteDepth(j);
            }
            if ((double)depthCnt/(postIndex-preIndex) > depthThresh) continue;
            seqList.add(chrSeq.substring(prePos-1, postPos-1));
            cnt++;
        }
        String[] seqs = seqList.toArray(new String[seqList.size()]);
        String[] names = new String[seqs.length];
        int[] ids = new int[seqs.length];
        for (int i = 0; i < names.length; i++) {
            names[i] = String.valueOf(i);
            ids[i] = i;
        }
        Fasta ff = new Fasta(names, seqs, ids);
        ff.writeFasta(outfileS);
    }
    
    public void plotAm560Divergenece () {
        String alignmentFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\AM560_diversity\\AM560.sam";
        String jPileupFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\jPileUp\\I000070_chr_01_MQ30.jPileUp.bin";
        String editDistanceFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\AM560_diversity\\AM560.ed.pdf";
        String depthFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\AM560_diversity\\AM560.depth.pdf";
        String depthDataFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\AM560_diversity\\AM560.depth.txt";
        int size = 10000;
        ShortreadAlignment sa = new ShortreadAlignment();
        sa.readFromBWAMEM(alignmentFileS);
        String[] querys = sa.getQuerys();
        TDoubleArrayList edList =new TDoubleArrayList();
        TDoubleArrayList depthList = new TDoubleArrayList();
        JPileUp jp = new JPileUp (jPileupFileS, IOFileFormat.Binary);
        for (int i = 0; i < querys.length; i++) {
            if (sa.isNoMatch(querys[i])) continue;
            int index = sa.getAlignmentStartIndexByQuery(querys[i]);
            edList.add(sa.getEditDistanceRatio(index));
            if (sa.getHit(index).equals("1")) {
                int siteCnt = 0;
                int startIndex = jp.getSiteIndex(1, sa.getStartPos(index));
                if (startIndex < 0) continue;
                int endIndex = jp.getSiteIndex(1, sa.getEndPos(index));
                if (endIndex < 0) continue;
                endIndex++;
                double d  = 0;
                for (int j = startIndex; j < endIndex; j++) {
                    d+=jp.getSiteDepth(j);
                    siteCnt++;
                }
                if (siteCnt == 0) continue;
                depthList.add(d/siteCnt);
            }
        }
        double[] edArray;
        if (edList.size()> size) edArray = edList.toArray(0, size);
        else edArray =edList.toArray();
        double[] depthArray;
        if (depthList.size() > size) depthArray = depthList.toArray(0, size);
        else depthArray = depthList.toArray();
        DensityPlot d =new DensityPlot(depthArray);
        d.setXLab("Depth");
        d.setYLab("Density");
        d.setTitle("");
        d.setXLim(0, 60);
        d.saveGraph(depthFileS);
        d =new DensityPlot(edArray);
        d.setXLab("Edit distance");
        d.setYLab("Density");
        d.setTitle("Edit distance of reference");
        d.saveGraph(editDistanceFileS);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(depthDataFileS);
            bw.write("Depth");
            bw.newLine();
            for (int i  = 0; i < depthArray.length; i++) {
                bw.write(String.valueOf(depthArray[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void generateAM560Fastq () {
        String infileS = "Q:\\D\\Cassava\\JGI\\FASTQ\\DRJL068_NoIndex_L002_R1_001.fastq.gz";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\AM560_diversity\\AM560.fq";
        int size = 200000;
        FastqChunk fq = new FastqChunk(infileS, ReadFormat.FastqGzip, 100000, size);
        fq.writeFastq(outfileS, ReadFormat.FastqText);
    }
    
    public void plotKu50Divergenece () {
        String alignmentFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\ku50_diversity\\ku50.sam";
        String jPileupFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\jPileUp\\I000070_chr_01_MQ30.jPileUp.bin";
        String editDistanceFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\ku50_diversity\\ku50.ed.pdf";
        String depthFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\ku50_diversity\\ku50.depth.pdf";
        String depthDataFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\ku50_diversity\\ku50.depth.txt";
        int size = 10000;
        ShortreadAlignment sa = new ShortreadAlignment();
        sa.readFromBWAMEM(alignmentFileS);
        String[] querys = sa.getQuerys();
        TDoubleArrayList edList =new TDoubleArrayList();
        TDoubleArrayList depthList = new TDoubleArrayList();
        JPileUp jp = new JPileUp (jPileupFileS, IOFileFormat.Binary);
        for (int i = 0; i < querys.length; i++) {
            if (sa.isNoMatch(querys[i])) continue;
            int index = sa.getAlignmentStartIndexByQuery(querys[i]);
            edList.add(sa.getEditDistanceRatio(index));
            if (sa.getHit(index).equals("1")) {
                int siteCnt = 0;
                int startIndex = jp.getSiteIndex(1, sa.getStartPos(index));
                if (startIndex < 0) continue;
                int endIndex = jp.getSiteIndex(1, sa.getEndPos(index));
                if (endIndex < 0) continue;
                endIndex++;
                double d  = 0;
                for (int j = startIndex; j < endIndex; j++) {
                    d+=jp.getSiteDepth(j);
                    siteCnt++;
                }
                if (siteCnt == 0) continue;
                depthList.add(d/siteCnt);
            }
        }
        double[] edArray;
        if (edList.size()> size) edArray = edList.toArray(0, size);
        else edArray =edList.toArray();
        double[] depthArray;
        if (depthList.size() > size) depthArray = depthList.toArray(0, size);
        else depthArray = depthList.toArray();
        DensityPlot d =new DensityPlot(depthArray);
        d.setXLab("Depth");
        d.setYLab("Density");
        d.setTitle("");
        d.setXLim(0, 60);
        d.saveGraph(depthFileS);
        d =new DensityPlot(edArray);
        d.setXLab("Edit distance");
        d.setYLab("Density");
        d.setTitle("Edit distance of ku50");
        d.saveGraph(editDistanceFileS);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(depthDataFileS);
            bw.write("Depth");
            bw.newLine();
            for (int i  = 0; i < depthArray.length; i++) {
                bw.write(String.valueOf(depthArray[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void generateKu50Fastq () {
        String infileS = "M:\\Database\\cassavaReference\\genome\\otherGenomes\\ku50\\ku50.fa";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\realDiversity\\ku50_diversity\\ku50.fq";
        Fasta f = new Fasta (infileS);
        int length = 200;
        int size = 200000;
        int cnt = 0;
        String[] seqs = new String[size];
        String[] names = new String[size];
        int[] ids =  new int[size];
        boolean flag = true;
        for (int i = 0; i < f.getSeqNumber(); i++) {
            if (flag == false) break;
            int n = f.getSeqLength(i)/length;
            for (int j = 0; j < n; j++) {
                seqs[cnt] = f.getSeq(i).substring(j*length, j*length+length);
                names[cnt] = String.valueOf(cnt);
                ids[cnt] = cnt;
                cnt++;
                if (cnt == size) {
                    flag = false;
                    break;
                }
            }
        }
        Fasta fout = new Fasta(names, seqs, ids);
        fout.writeFastq(outfileS);
    }
    
    public void genePipe () {
        //this.getGeneSequence();
        
        //this.getConservedOrthologs();
        //this.getConvervedCassavaGene();
        //this.getDivergenceOfGene();
        
        //this.getConservedGene();
        this.getConservedGeneEditDistance();
        
    }
    
    public void getConservedGeneEditDistance() {
        String alignmentFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\geneDivergence\\align_cassava6OnAssembly_oout.txt";
        String geneFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\geneDivergence\\cassava_conservedGene.fa";
        String geneFeatureFileS = "M:\\Database\\cassavaReference\\Gene\\cassava6.gf.txt";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\geneDivergence\\cassava_conservedGene_report.txt";
        BlastAlignment ba = new BlastAlignment();
        ba.readFromBlast(alignmentFileS, 1);
        ba.sortByQuery();
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        gf.sortGeneByName();
        Fasta f = new Fasta(geneFileS);
        f.sortRecordByName();
        String[] query = ba.getQuerys();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Gene\tAllLength\tExonLength\tIntronLength\tHitContigName\tAllAlignedLength\tExonAlignedLength\tIntronAlignedLength\tAllED\tExonED\tIntronED");
            bw.newLine();
            for (int i = 0; i < query.length; i++) {
                if (ba.isNoMatch(query[i])) continue;
                int length = f.getSeqLength(f.getIndex(query[i]));
                int index = gf.getGeneIndex(query[i].split("_")[3]);
                int geneStart = gf.getGeneStart(index);
                ArrayList<Range> cdsList = gf.getCDSList(index, 0);
                ArrayList<Range> intronList = gf.getIntronList(index, 0);
                int exonLength = 0;
                for (int j = 0; j < cdsList.size(); j++) {
                    Range r = cdsList.get(j);
                    r.setRangeStart(r.getRangeStart()-geneStart+1);
                    r.setRangeEnd(r.getRangeEnd()-geneStart+1);
                    exonLength+=(r.getRangeEnd()-r.getRangeStart());
                }
                int intronLength = 0;
                if (intronList.isEmpty()) {
                    
                }
                else {
                    for (int j = 0; j < intronList.size(); j++) {
                        Range r = intronList.get(j);
                        r.setRangeStart(r.getRangeStart()-geneStart+1);
                        r.setRangeEnd(r.getRangeEnd()-geneStart+1);
                        intronLength+=(r.getRangeEnd()-r.getRangeStart());
                    }
                }
                double[] value = new double[length];
                for (int j = 0; j < value.length; j++) {
                    value[j] = -1;
                }
                int startIndex = ba.getAlignmentStartIndexByQuery(query[i]);
                String contigName = ba.getHit(startIndex);
                int endIndex = ba.getAlignmentEndIndexByQuery(query[i]);
                
                for (int j = startIndex; j < endIndex; j++) {
                    if (!ba.getHit(j).equals(contigName)) {
                        endIndex = j;
                        break;
                    }
                }
                
                for (int j = startIndex; j < endIndex; j++) {
                    double d = ba.getEditDistanceRatio(j);
                    int qStartIndex = ba.getQueryStartPos(j)-1;
                    int qEndIndex = ba.getQueryEndPos(j)-1;
                    for (int k = qStartIndex; k < qEndIndex; k++) {
                        if (value[k] == -1) value[k] = d;
                    }
                }
                int exonAlignedLength = 0;
                double exonED = 0;
                for (int j = 0; j < cdsList.size(); j++) {
                    Range r = cdsList.get(j);
                    for (int k = r.getRangeStart()-1; k < r.getRangeEnd()-1; k++) {
                        if (value[k] == -1) continue;
                        exonAlignedLength++;
                        exonED+=value[k];
                    }
                }
                if (exonAlignedLength == 0) {
                    exonED = Double.NaN;
                }
                else exonED = exonED/exonAlignedLength;
                int intronAlignedLength = 0;
                double intronED = 0;
                if (intronList.isEmpty()) {
                    intronED = Double.NaN;
                }
                else {
                    for (int j = 0; j < intronList.size(); j++) {
                        Range r = intronList.get(j);
                        for (int k = r.getRangeStart()-1; k < r.getRangeEnd()-1; k++) {
                            if (value[k] == -1) continue;
                            intronAlignedLength++;
                            intronED+=value[k];
                        }
                    }
                    if (intronAlignedLength == 0) {
                        intronED = Double.NaN;
                    }
                    else intronED = intronED/intronAlignedLength;
                }
                int allAlignedLength = 0;
                double allED = 0;
                for (int j = 0; j < value.length; j++) {
                    if (value[j] == -1) continue;
                    allAlignedLength++;
                    allED+=value[j];
                }
                if (allAlignedLength == 0) {
                    allED = Double.NaN;
                }
                else allED = allED/allAlignedLength;
                StringBuilder sb = new StringBuilder();
                sb.append(query[i]).append("\t").append(length).append("\t").append(exonLength).append("\t").append(intronLength).append("\t").append(contigName).append("\t");
                sb.append(allAlignedLength).append("\t").append(exonAlignedLength).append("\t").append(intronAlignedLength).append("\t");
                sb.append(allED).append("\t").append(exonED).append("\t").append(intronED);
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
    public void getConservedGene () {
        String geneFileS = "M:\\Database\\cassavaReference\\Gene\\cassava6_gene.fa";
        String alignmentFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\geneDivergence\\align_cassavaOnAra.txt";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\geneDivergence\\cassava_conservedGene.fa";
        ShortreadAlignment sa = new ShortreadAlignment();
        sa.readFromBlast(alignmentFileS, 1, 1);
        String[] querys = sa.getQuerys();
        ArrayList<String> geneList = new ArrayList();
        for (int i = 0; i < querys.length; i++) {
            if (sa.isNoMatch(querys[i])) continue;
            int startIndex = sa.getAlignmentStartIndexByQuery(querys[i]);
            int endIndex = sa.getAlignmentEndIndexByQuery(querys[i]);
            HashSet<String> s = new HashSet();
            for (int j = startIndex; j < endIndex; j++) {
                s.add(sa.getHit(j).split("\\.")[1]);
            }
            if (s.size() == 1) geneList.add(querys[i]);
        }
        Fasta f = new Fasta (geneFileS);
        f.sortRecordByName();
        boolean[] ifOut = new boolean[f.getSeqNumber()];
        for (int i = 0; i < geneList.size(); i++) {
            ifOut[f.getIndex(geneList.get(i))] = true;
        }
        f.writeFasta(outfileS, ifOut);
    }
    
    public void getDivergenceOfGene () {
        String alignFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\geneDivergence\\conserved_discovar_fm7.txt";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\geneDivergence\\conserved_discovar_divergent.txt";
        ShortreadAlignment sra = new ShortreadAlignment();
        sra.readFromBlast(alignFileS, 1, 100);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            String[] queries = sra.getQuerys();
            for (int i = 0; i < queries.length; i++) {
                int index = sra.getAlignmentStartIndexByQuery(queries[i]);
                double d = (double)sra.getEditDistance(index)/sra.getMatchNumber(index);
                bw.write(queries[i]+"\t"+String.valueOf(d));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void getConvervedCassavaGene () {
        String alignFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\geneDivergence\\align_out_fm7.txt";
        String cassavaGeneFileS = "M:\\Database\\cassavaReference\\Gene\\cassava6_gene.fa";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\geneDivergence\\cassava_convervedGene.fa";
        ShortreadAlignment sra = new ShortreadAlignment();
        sra.readFromBlast(alignFileS, 1, 100);
        String[] hits = sra.getHits();
        Fasta f = new Fasta(cassavaGeneFileS);
        f.sortRecordByName();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            for (int i = 0; i < hits.length; i++) {
                int index = f.getIndex(hits[i]);
                bw.write(">"+f.getName(index));
                bw.newLine();
                bw.write(FStringUtils.getMultiplelineString(60, f.getSeq(index)));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void getConservedOrthologs () {
        String infileDirS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\geneDivergence\\cosii\\ortholog_sequences\\";
        String cosFileS = "M:\\Database\\conservedGene\\cos.cds.fa";
        File[] fs = IoUtils.listFilesEndsWith(new File(infileDirS).listFiles(), ".cds.fasta");
        try {
            BufferedWriter bw = IoUtils.getTextWriter(cosFileS);
            for (int i = 0; i < fs.length; i++) {
                Fasta f = new Fasta(fs[i].getAbsolutePath());
                for (int j = 0; j < f.getSeqNumber(); j++) {
                    if (f.getName(j).startsWith("At")) {
                        String seq = f.getSeq(j).replaceAll("-", "");
                        String[] subSeqs = FStringUtils.getMultilineString(60, seq);
                        bw.write(">"+f.getName(j));
                        bw.newLine();
                        for (int k = 0; k < subSeqs.length; k++) {
                            bw.write(subSeqs[k]);
                            bw.newLine();
                        }
                    }
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println(String.valueOf(fs.length)+" genes processed");
    }
    
    public void getGeneSequence () {
        String gffFileS = "M:\\Database\\cassavaReference\\Gene\\MesculentaAnnotationV6.1_pr2.gff3.txt";
        String cassavaGenomeFileS = "M:\\Database\\cassavaReference\\genome\\Manihot esculenta\\cassavaV6_chrAndScaffoldsCombined_numeric.fa";
        String geneSequenceFileS = "M:\\Database\\cassavaReference\\Gene\\cassava6_gene.fa";
        String cdsSequenceFileS = "M:\\Database\\cassavaReference\\Gene\\cassava6_cds.fa";
        Fasta f = new Fasta(cassavaGenomeFileS);
        GeneFeature gf = new GeneFeature();
        gf.readFromCassavaGFF(gffFileS);
        gf.writeGeneSequence(f, geneSequenceFileS);
        gf.writeCDSSequence(f, cdsSequenceFileS);
    }
    
    public void assemblyPipe () {
        //this.renameAssembly();
        //this.getAssemlyStatistics();
        //this.mkFastq();
        //this.getEditDistanceOfContig();
        //this.contigDistanceDensity();
        //this.getOneContigAndAlign();
    }
    
    public void getOneContigAndAlign () {
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\assembly\\genome\\discovar_I0000070_GTCCGC.fa";
        String oufileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\assembly\\contig7.fa";
        Fasta f = new Fasta (infileS);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(oufileS);
            for (int i = 0; i < f.getSeqNumber(); i++) {
                if (f.getName(i).equals("7")) {
                    bw.write(f.getName(i));
                    bw.newLine();
                    bw.write(f.getSeq(i));
                    bw.newLine();
                    break;
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void getEditDistanceOfContig () {
        String samFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\assembly\\out.sam.gz";
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\assembly\\genome\\discovar_I0000070_GTCCGC.fa";
        Fasta f = new Fasta(infileS);
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\assembly\\contig_distance.txt";
        String outfileS2 = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\assembly\\contig_distance2.txt";
        ShortreadAlignment sa = new ShortreadAlignment();
        sa.readFromBWAMEM(samFileS);
        String[] query = sa.getQuerys();
        HashSet<Integer> contigSet = new HashSet();
        for (int i = 0; i < query.length; i++) {
            contigSet.add(Integer.valueOf(query[i].split("_")[0]));
        }
        TDoubleArrayList[] dLists = new TDoubleArrayList[contigSet.size()];
        for (int i = 0; i < dLists.length; i++) {
            dLists[i] = new TDoubleArrayList();
        }
        for (int i = 0; i < query.length; i++) {
            int startIndex = sa.getAlignmentStartIndexByQuery(query[i]);
            double d;
            if (sa.isMatch(startIndex)) {
                d = (double)sa.getEditDistance(startIndex)/sa.getMatchNumber(startIndex);
            }
            else {
                d = 1;
            }
            int index = Integer.valueOf(query[i].split("_")[0]);
            dLists[index].add(d);
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            BufferedWriter bw2 = IoUtils.getTextWriter(outfileS2);
            bw.write("Contig\tLength\tEditDistance\tEditDistanceWithoutNoMap");
            bw.newLine();
            for (int i = 0; i < dLists.length; i++) {
                double[] value = dLists[i].toArray();
                double d = 0;
                double d2 = 0;
                int cnt = 0;
                for (int j = 0; j < value.length; j++) {
                    d+=(value[j]/value.length);
                    bw2.write(String.valueOf(value[j])+"\t");
                    if (value[j] != 1) {
                        d2+=value[j];
                        cnt++;
                    }
                }
                bw2.newLine();
                d2 = d2/cnt;
                bw.write(String.valueOf(i)+"\t"+String.valueOf(f.getSeqLength(i))+"\t"+String.valueOf(d)+"\t"+String.valueOf(d2));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            bw2.flush();
            bw2.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkFastq () {
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\assembly\\genome\\discovar_I0000070_GTCCGC.fa";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\assembly\\query.fq";
        int minLength = 3000;
        int readLength = 300;
        Fasta f = new Fasta(infileS);
        Fasta ff = new Fasta();
        for (int i = 0; i < f.getSeqNumber(); i++) {
            if (f.getSeqLength(i) < minLength) continue;
            Sequence s = new Sequence(f.getSeq(i));
            String[] fragments = s.getFragments(readLength);
            int n = fragments.length;
            if (fragments[fragments.length-1].length() < readLength) n--;
            String[] names = new String[n];
            String[] seqs = new String[n];
            int[] ids = new int[n];
            for (int j = 0; j < names.length; j++) {
                names[j] = f.getName(i)+"_"+String.valueOf(j);
                seqs[j] = fragments[j];
                ids[j] = j;
            }
            Fasta nf = new Fasta(names, seqs, ids);
            ff = ff.getMergedFasta(nf);
            if (i!= 0 && i%1000 == 0) System.out.println(i+"\t"+f.getSeqLength(i));
        }
        ff.writeFastq(outfileS);
    }
    
    public void getAssemlyStatistics () {
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\assembly\\genome\\discovar_I0000070_GTCCGC.fa";
        Fasta f = new Fasta(infileS);
        
        System.out.println("Total length: "+ f.getTotalSeqLength());
        System.out.println("Sequence number: "+f.getSeqNumber());
        System.out.println("L50: "+ f.getL50());
        System.out.println("N50: " + f.getN50());
    }
    
    public void renameAssembly () {
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\assembly\\genome\\discovar_I0000070_GTCCGC.fa";
        Fasta f = new Fasta(infileS);
        f.sortRecordByLengthDescending();
        for (int i = 0; i < f.getSeqNumber(); i++) {
            f.setName(String.valueOf(i), i);
        }
        f.writeFasta(infileS);
    }
    
    public void depthPipe() {
        //this.extractDepthFromSam();
        //this.splitDepthByChromosome();
        //this.convertToJPileUp();
        //this.depthDensity();
        //this.depthAndEditDistance();
        //this.depthAndEditDistancePlot();
        //this.readDepthLength();
        //this.readSegregation();
        //this.depthDistributionOnRatio();
        //this.editDistanceInDepthRange();
        //this.convertJPileUP10Clone();
        //this.depthCorrelationBetweenClones();
        //this.depthDensityReference();
        this.depthDensityReference2();
    }
    
    public void depthCorrelationBetweenClones () {
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileupTenClone\\jPileup\\I000070_chr_01_MQ30.jPileUp.bin";
        String infileDirS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileupTenClone\\jPileup\\";
        String outDirS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileupTenClone\\pdf_depthCorrelation\\";
        JPileUp mjp = new JPileUp(infileS, IOFileFormat.Binary);
        int size = 1000;
        int depthLowThresh = 1;
        int depthHighThresh  = 100;
        
        int yLimLow = 0;
        int yLimHigh = 500;
        int[] position = new int[size];
        int cnt = 0;
        while (cnt < size) {
            int ind = (int)(mjp.getSiteNumber()*Math.random());
            if (mjp.getSiteDepth(ind)< depthLowThresh || mjp.getSiteDepth(ind) > depthHighThresh) continue;
            position[cnt] = mjp.getPosition(ind);
            cnt++;
        }
        File[] fs = new File (infileDirS).listFiles();
        double[][] depths = new double[fs.length][size];
        for (int i = 0; i < fs.length; i++) {
            JPileUp jp = new JPileUp(fs[i].getAbsolutePath(), IOFileFormat.Binary);
            TDoubleArrayList xList = new TDoubleArrayList();
            TDoubleArrayList yList = new TDoubleArrayList();
            
            for (int j = 0; j < size; j++) {
                int index = jp.getSiteIndex(1, position[j]);
                if (index < 0) continue;
                yList.add(jp.getSiteDepth(index));
                index = mjp.getSiteIndex(1, position[j]);
                xList.add(mjp.getSiteDepth(index));
            }
            double[] x = xList.toArray();
            double[] y = yList.toArray();
            String taxaName = fs[i].getName().split("_")[0];
            String outfileS = new File(outDirS, taxaName+".pdf").getAbsolutePath();
            double r = new PearsonsCorrelation().correlation(x, y);
            ScatterPlot s = new ScatterPlot(x, y);
            
            s.setXLim(depthLowThresh, depthHighThresh);
            s.setYLim(yLimLow, yLimHigh);
            s.setPlottingCharacter(19);
            s.setColor(255, 0, 0, 125);
            s.setXLab("Depth of I000070");
            s.setYLab("Depth of " + taxaName);
            s.setTitle("Depth correlation of two clones, range from " +String.valueOf(depthLowThresh)+" to " + String.valueOf(depthHighThresh) +", r = " + String.valueOf(r));
            s.saveGraph(outfileS);
        }
    }
    
    public void convertJPileUP10Clone () {
        String infileDirS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileupTenClone\\pileup\\";
        String outDirS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileupTenClone\\jPileup\\";
        File[] fs = new File(infileDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String outfileS = f.getName().replaceFirst(".txt", ".jPileUp.bin");
            outfileS = new File(outDirS, outfileS).getAbsolutePath();
            JPileUp jp = new JPileUp();
            jp.readFromPileup(f.getAbsolutePath());
            jp.writeFile(outfileS, IOFileFormat.Binary);
        });
    }
    
    public void editDistanceInDepthRange () {
        int lowThresh = 40;
        int highThresh = 5000;
        String pileupFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\jPileUp\\I000070_chr_01_MQ30.jPileUp.bin";
        String samFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\sam_30\\I000070.sam.gz";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\pdf\\edit_distance_40_5000.pdf";
        JPileUp jp = new JPileUp(pileupFileS, IOFileFormat.Binary);
        TDoubleArrayList dList = new TDoubleArrayList();
        ShortreadPEAlignment spe = new ShortreadPEAlignment();
        spe.readFromBWAMEM(samFileS);
        int chr = jp.getChromosome(0);
        for (int i = 0; i < spe.getAlignmentNumber(); i++) {
            if (!spe.isMatchF(i)) continue;
            if (!spe.isMatchB(i)) continue;
            if (Integer.valueOf(spe.getHitF(i)) != chr || Integer.valueOf(spe.getHitB(i)) != chr) continue;
            int startF = spe.getStartPosF(i);
            int startB = spe.getStartPosB(i);
            int endF = spe.getEndPosF(i);
            int endB = spe.getEndPosB(i);
            double depth = 0;
            double siteCnt = 0;
            for (int j = startF; j < endF; j++) {
                int index = jp.getSiteIndex(chr, j);
                if (index < 0) continue;
                depth+=jp.getSiteDepth(index);
                siteCnt++;
            }
            for (int j = startB; j < endB; j++) {
                int index = jp.getSiteIndex(chr, j);
                if (index < 0) continue;
                depth+=jp.getSiteDepth(index);
                siteCnt++;
            }
            if (siteCnt == 0) continue;
            depth= depth/siteCnt;
            if (depth < lowThresh) continue;
            if (depth > highThresh) continue;
            double dF = (double)spe.getEditDistanceF(i)/spe.getMatchNumberF(i);
            double dB = (double)spe.getEditDistanceB(i)/spe.getMatchNumberB(i);
            dList.add((dF+dB)/2);
            if (dList.size() > 5000) break;
        }
        double[] dArray = dList.toArray();
        DensityPlot d = new DensityPlot(dArray);
        d.setXLab("Edit distance");
        d.setYLab("Density");
        d.setTitle("Edit distance distribution from I000070, read depth from " +String.valueOf(lowThresh) + " to " +String.valueOf(highThresh));
        d.saveGraph(outfileS);
        
    }
    
    public void depthDistributionOnRatio () {
        int size = 50000;
        int lowThresh = 10;
        int highThresh = 40;
        int qThresh = 20;
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\jPileUp\\I000070_chr_01_MQ30.jPileUp.bin";
        String outfileS1 = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\pdf\\depthDistributionOnRatio_01.pdf";
        String outfileS2 = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\pdf\\depthDistributionOnRatio_02.pdf";
        String outfileS3 = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\pdf\\depthDistributionOnRatio_03.pdf";
        TDoubleArrayList[] lists = new TDoubleArrayList[3];
        for (int i = 0; i < lists.length; i++) lists[i] = new TDoubleArrayList();
        JPileUp jp = new JPileUp(infileS, IOFileFormat.Binary);
        int cnt = 0;
        while (cnt < size) {
            int index = (int)(Math.random()*jp.getSiteNumber());
            if (jp.getAlleleNumber(index) == 0) continue;
            if (jp.getSiteDepth(index) < lowThresh) continue;
            if (jp.getSiteDepth(index) > highThresh) continue;
            int alleleNum = jp.getAlleleDepth(index, 0);
            int alleleQNum = 0;
            for (int j = 0; j < jp.getAlleleDepth(index, 0); j++) {
                if (jp.getAlleleBaseQuality(index, 0, j) < qThresh) continue;
                alleleQNum++;
            }
            int refNum = jp.getRefDepth(index);
            int refQNum = 0;
            for (int j = 0; j < jp.getRefDepth(index); j++) {
                if (jp.getRefBaseQuality(index, j) < qThresh) continue;
                refQNum++;
            }
            int sum = jp.getRefDepth(index)+ jp.getAlleleDepth(index, 0);
            int sumQ = alleleQNum+refQNum;
            double ratio = (double)alleleQNum/(sumQ);
            if (ratio < 0.2) {
                if (sumQ>=lowThresh) lists[0].add(sumQ);
            }
            else if (ratio  < 0.8) {
                if (sumQ>=lowThresh) lists[1].add(sumQ);
            }
            else {
                if (sumQ>=lowThresh) lists[2].add(sumQ);
            }
            cnt++;
            if (cnt%100 == 0) System.out.println(cnt);
        }
        String[] titleNames = {"Read depth of sites with segregation ratio 0-0.2", "Read depth of sites with segregation ratio 0.2-0.8", "Read depth of sites with segregation ratio 0.8-0.1"};
        double[][] dis = new double[3][];
        for (int i = 0; i < dis.length; i++) {
            dis[i] = lists[i].toArray();
        }
        DensityPlot d = new DensityPlot(dis[0]);
        d.setTitle(titleNames[0]);
        d.setXLab("Read depth");
        d.setYLab("Frequency");
        d.saveGraph(outfileS1);
        d = new DensityPlot(dis[1]);
        d.setTitle(titleNames[1]);
        d.setXLab("Read depth");
        d.setYLab("Frequency");
        d.saveGraph(outfileS2);
        d = new DensityPlot(dis[2]);
        d.setTitle(titleNames[2]);
        d.setXLab("Read depth");
        d.setYLab("Frequency");
        d.saveGraph(outfileS3);
    }
    
    public void readSegregation () {
        int size = 10000;
        int lowThresh = 4;
        int highThresh = 10;
        int qThresh = 20;
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\jPileUp\\I000070_chr_01_MQ30.jPileUp.bin";
        String outfileS1 = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\pdf\\alleleReadSegregaton_4_10.pdf";
        String outfileS2 = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\pdf\\alleleReadSegregaton_4_10_Q20.pdf";
        JPileUp jp = new JPileUp(infileS, IOFileFormat.Binary);
        double[] ratio = new double[size];
        TDoubleArrayList ratioQList = new TDoubleArrayList();
        int cnt = 0;
        while (cnt < size) {
            int index = (int)(Math.random()*jp.getSiteNumber());
            if (jp.getAlleleNumber(index) == 0) continue;
            if (jp.getSiteDepth(index) < lowThresh) continue;
            if (jp.getSiteDepth(index) > highThresh) continue;
            int alleleNum = jp.getAlleleDepth(index, 0);
            int alleleQNum = 0;
            for (int j = 0; j < jp.getAlleleDepth(index, 0); j++) {
                if (jp.getAlleleBaseQuality(index, 0, j) < qThresh) continue;
                alleleQNum++;
            }
            int refNum = jp.getRefDepth(index);
            int refQNum = 0;
            for (int j = 0; j < jp.getRefDepth(index); j++) {
                if (jp.getRefBaseQuality(index, j) < qThresh) continue;
                refQNum++;
            }
            if (alleleQNum != 0) {
                ratioQList.add((double)alleleQNum/(alleleQNum+refQNum));
            }
            ratio[cnt] = (double)alleleNum/(alleleNum+refNum);
            
            cnt++;
            if (cnt%100 == 0) System.out.println(cnt);
        }
        Histogram d = new Histogram(ratio);
        d.setXLab("Allele read segregation ratio");
        d.setYLab("Frequency");
        d.setTitle("Allele read segregation ratio in one WGS clone (I000070)");
        d.setBreakNumber(50);
        d.saveGraph(outfileS1);
        double[] ratioQ = ratioQList.toArray();
        d = new Histogram(ratioQ);
        d.setXLab("Allele read segregation ratio");
        d.setYLab("Frequency");
        d.setTitle("Allele read segregation ratio in one WGS clone (I000070), baseQ > " +String.valueOf(qThresh));
        d.setBreakNumber(50);
        d.saveGraph(outfileS2);
    }
    
    public void readDepthLength() {
        int highThresh = 45;
        int lowThresh = 5;
        int stepSize = 500;
        int stepNumber = 1000;
        int size = 200;
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\jPileUp\\I000070_chr_01_MQ30.jPileUp.bin";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\readDepthChange.txt";
        JPileUp jp = new JPileUp(infileS, IOFileFormat.Binary);
        
        ArrayList<double[]> highList = new ArrayList();
        ArrayList<double[]> lowList = new ArrayList();
        int[] highIndex = new int[size];
        int[] lowIndex = new int[size];
        int upLimit = stepSize*stepNumber;
        int downLimit = jp.getPosition(jp.getSiteNumber()-1);
        int cnt = 0;
        while (cnt < size) {
            int index = (int)(Math.random()*jp.getSiteNumber());
            if (jp.getSiteDepth(index)> lowThresh) continue;
            if (jp.getPosition(index) < upLimit) continue;
            if (jp.getPosition(index) > downLimit) continue;
            lowIndex[cnt] = index;
            cnt++;
        }
        cnt = 0;
        while (cnt < size) {
            int index = (int)(Math.random()*jp.getSiteNumber());
            if (jp.getSiteDepth(index) < highThresh) continue;
            if (jp.getPosition(index) < upLimit) continue;
            if (jp.getPosition(index) > downLimit) continue;
            highIndex[cnt] = index;
            cnt++;
        }
        for (int i = 0; i < lowIndex.length; i++) {
            double[] value = new double[stepNumber];
            for (int j = 0; j < stepNumber; j++) {
                int lastPos = jp.getPosition(lowIndex[i])-(j+1)*stepSize/2;
                int nextPos = jp.getPosition(lowIndex[i])+(j+1)*stepSize/2;
                int lastIndex = jp.getNearbySiteIndex(1, lastPos);
                int nextIndex = jp.getNearbySiteIndex(1, nextPos);
                int totalDepth = 0;
                for (int k = lastIndex; k < nextIndex; k++) {
                    totalDepth+=jp.getSiteDepth(k);
                }
                value[j] = (double)totalDepth/(nextIndex-lastIndex);
            }
            lowList.add(value);
            System.out.println(i);
        }
        for (int i = 0; i < highIndex.length; i++) {
            double[] value = new double[stepNumber];
            for (int j = 0; j < stepNumber; j++) {
                int lastPos = jp.getPosition(highIndex[i])-(j+1)*stepSize/2;
                int nextPos = jp.getPosition(highIndex[i])+(j+1)*stepSize/2;
                int lastIndex = jp.getNearbySiteIndex(1, lastPos);
                int nextIndex = jp.getNearbySiteIndex(1, nextPos);
                int totalDepth = 0;
                for (int k = lastIndex; k < nextIndex; k++) {
                    totalDepth+=jp.getSiteDepth(k);
                }
                value[j] = (double)totalDepth/(nextIndex-lastIndex);
            }
            highList.add(value);
            System.out.println(i);
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Size\tLowDepthTrend\tHighDepthTrend");
            bw.newLine();
            for (int i = 0; i < stepNumber; i++) {
                double lowMean = 0;
                double highMean = 0;
                for (int j = 0; j < lowList.size(); j++) lowMean+=lowList.get(j)[i]/lowList.size();
                for (int j = 0; j < highList.size(); j++) highMean+=highList.get(j)[i]/highList.size();
                bw.write(String.valueOf(stepSize*(i+1))+"\t"+String.valueOf(lowMean)+"\t"+String.valueOf(highMean));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void depthAndEditDistancePlot () {
        //String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\readDepth_editDistance.txt";
        //String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\readDepth_editDistance.pdf";
        
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileupTenClone\\readDepth_editDistance_reference.txt";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileupTenClone\\readDepth_editDistance_reference.pdf";
        Table t = new Table (infileS);
        int size = 1000;
        if (size > t.getRowNumber()) size = t.getRowNumber();
        double[] depth = new double[size];
        double[] dis = new double[size];
        for (int i = 0; i < size; i++) {
            depth[i] = t.getDoubleValue(i, 0);
            dis[i] = t.getDoubleValue(i, 1);
        }
        DensityPlot d = new DensityPlot(dis);
        d.setXLim(0, 0.3);
        //d.showGraph();
        PearsonsCorrelation p = new PearsonsCorrelation();
        double r = p.correlation(depth, dis);
        ScatterPlot s = new ScatterPlot(depth, dis);
        s.setXLim(0, 100);
        s.setXLab("Read depth");
        s.setYLab("Edit distance");
        s.setTitle("r = "+ String.valueOf(r));
        s.saveGraph(outfileS);
    }
    
    public void depthAndEditDistance () {
        //String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\jPileUp\\I000070_chr_01_MQ30.jPileUp.bin";
        //String outfile1S = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\readDepth_editDistance.txt";
        
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileupTenClone\\jPileup\\reference_chr_01_MO30.jPileUp.bin";
        String outfile1S = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileupTenClone\\readDepth_editDistance_reference.txt";
        JPileUp jp = new JPileUp(infileS, IOFileFormat.Binary);
        int size = 20000;
        int length = 500;
        int qThresh = 0;
        int[] index = new int[size];
        int cnt = 0;
        int lastPos = jp.getPosition(jp.getSiteNumber()-1)-length;
        while (cnt < size) {
            int v = (int)(Math.random()*jp.getSiteNumber());
            if (jp.getPosition(v)> length && jp.getPosition(v) < lastPos) {
                index[cnt] = v;
            }
            else {
                cnt--;
                if (cnt < 0) cnt = 0;
            }
            cnt++;
        }
        double[] depth = new double[size];
        double[] distance = new double[size];
        for (int i = 0; i < index.length; i++) {
            int startIndex = index[i]-length/2;
            int endIndex = index[i]+length/2;
            int totalDepth = 0;
            int d = 0;
            for (int j = startIndex; j < endIndex; j++) {
                totalDepth+=jp.getSiteDepth(j);
                for (int k = 0; k < jp.getAlleleNumber(j); k++) {
                    byte[] quals = jp.getAlleleBaseQuality(j, k);
                    for (int u = 0; u < quals.length; u++) {
                        if (quals[u] < qThresh) continue;
                        d++;
                    }
                }
            }
            depth[i] = totalDepth/length;
            if (totalDepth == 0) {
                distance[i] = Double.NaN;
            }
            else {
                distance[i] = (double)d/totalDepth;
            }
            
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfile1S);
            bw.write("DepthInRegion\tEditDisance");
            bw.newLine();
            for (int i = 0; i < depth.length; i++) {
                if (Double.isNaN(distance[i])) continue;
                bw.write(String.valueOf(depth[i])+"\t"+String.valueOf(distance[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void depthDensityReference () {
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileupTenClone\\jPileup\\reference_chr_01_MQ10.jPileUp.bin";
        String outfile1S = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileupTenClone\\read_depth.hist.pdf";
        String outfile2S = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileupTenClone\\read_depth.cdf.pdf";
        String depthFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileupTenClone\\readDepth.txt";
        JPileUp jp = new JPileUp(infileS, IOFileFormat.Binary);
        int size = 20000;
        double[] depth = new double[size];
        for (int i = 0; i < size; i++) {
            int index = (int)(Math.random()*jp.getSiteNumber());
            depth[i] = jp.getSiteDepth(index);
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(depthFileS);
            bw.write("readDepth");
            bw.newLine();
            for (int i = 0; i < depth.length; i++) {
                bw.write(String.valueOf(depth[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch(Exception e) {
            e.printStackTrace();
        }
        Histogram h  =new Histogram(depth);
        h.setXLim(0, 300);
        h.setBreakNumber(40000);
        h.setXLab("Read depth");
        h.setTitle("Read depth distribution of I000070");
        h.saveGraph(outfile1S);
        CumulativeDistribution cd = new CumulativeDistribution(depth);
        cd.setXLab("Read depth");
        cd.setYLab("Frequency");
        cd.setTitle("Cumulative read depth distribution of I000070");
        cd.setXLim(0, 100);
        cd.saveGraph(outfile2S);
    }
    
    public void depthDensityReference2 () {
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileupTenClone\\jPileup\\reference_chr_01_MQ10.jPileUp.bin";
        String referenceFileS = "M:\\Database\\cassavaReference\\genome\\Manihot esculenta\\cassavaV6_chrAndScaffoldsCombined_numeric.fa";
        String outfile1S = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileupTenClone\\refDepth\\read_depth.hist.pdf";
        String depthFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileupTenClone\\refDepth\\readDepth.txt";
        JPileUp jp = new JPileUp(infileS, IOFileFormat.Binary);
        int size = 20000;
        Fasta f = new Fasta(referenceFileS);
        String seq = f.getSeq(0);
        double[] depth = new double[size];
        for (int i = 0; i < size; i++) {
            int index = (int)(Math.random()*seq.length());
            char c = seq.charAt(index);
            if (c == 'N' || c == 'n') {
                i--;
                continue;
            }
            index = jp.getSiteIndex(1, index+1);
            if (index < 0) {
                depth[i] = 0;
            }
            else {
                depth[i] = jp.getSiteDepth(index);
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(depthFileS);
            bw.write("readDepth");
            bw.newLine();
            for (int i = 0; i < depth.length; i++) {
                bw.write(String.valueOf(depth[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch(Exception e) {
            e.printStackTrace();
        }
        Histogram h  =new Histogram(depth);
        h.setXLim(0, 300);
        h.setBreakNumber(40000);
        h.setXLab("Read depth");
        h.setTitle("Read depth distribution of I000070");
        h.saveGraph(outfile1S);
    }
    
    public void depthDensity () {
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\jPileUp\\I000070_chr_01_MQ30.jPileUp.bin";
        String outfile1S = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\pdf\\read_depth.hist.pdf";
        String outfile2S = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\pdf\\read_depth.cdf.pdf";
        String depthFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\readDepth.txt";
        JPileUp jp = new JPileUp(infileS, IOFileFormat.Binary);
        int size = 20000;
        double[] depth = new double[size];
        for (int i = 0; i < size; i++) {
            int index = (int)(Math.random()*jp.getSiteNumber());
            depth[i] = jp.getSiteDepth(index);
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(depthFileS);
            bw.write("readDepth");
            bw.newLine();
            for (int i = 0; i < depth.length; i++) {
                bw.write(String.valueOf(depth[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch(Exception e) {
            e.printStackTrace();
        }
        Histogram h  =new Histogram(depth);
        h.setXLim(0, 60);
        h.setBreakNumber(5000);
        h.setXLab("Read depth");
        h.setTitle("Read depth distribution of I000070");
        h.saveGraph(outfile1S);
        CumulativeDistribution cd = new CumulativeDistribution(depth);
        cd.setXLab("Read depth");
        cd.setYLab("Frequency");
        cd.setTitle("Cumulative read depth distribution of I000070");
        cd.setXLim(0, 100);
        cd.saveGraph(outfile2S);
    }
    
    public void convertToJPileUp () {
        //String inDirS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\bychr\\";
        //String outDirS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\mpileup\\jPileUp\\";
        String inDirS = "/workdir/fl262/pileup/";
        String outDirS = "/workdir/fl262/jPileup/";
        File[] fs = new File(inDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String outfileS = f.getName().replaceFirst(".pileup.txt", ".jPileUp.bin");
            outfileS = new File(outDirS, outfileS).getAbsolutePath();
            JPileUp jp = new JPileUp();
            jp.readFromPileup(f.getAbsolutePath());
            jp.writeFile(outfileS, IOFileFormat.Binary);
        });
    }
    
    public void splitDepthByChromosome () {
        String depthFileS = "/workdir/mingh/I000070.MQ30.pileup.txt";
        String depthDirS = "/workdir/mingh/byChr/";
        try {
            BufferedReader br = IoUtils.getTextReader(depthFileS);
            BufferedWriter bw = null;
            String temp;
            String[] tem;
            int current = -1;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                tem = temp.substring(0, 20).split("\t");
                int chr = Integer.valueOf(tem[0]);
                if (chr == current) {
                    bw.write(temp);
                    bw.newLine();
                }
                else {
                    if (bw != null) {
                        bw.flush();
                        bw.close();
                    }
                    current = chr;
                    String outfileS = "I000070_chr_"+FStringUtils.getNDigitNumber(2, chr)+"_MQ30.pileup.txt";
                    outfileS = new File(depthDirS, outfileS).getAbsolutePath();
                    bw = IoUtils.getTextWriter(outfileS);
                    bw.write(temp);
                    bw.newLine();
                }
                cnt++;
                if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt)+" sites");
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void extractDepthFromSam () {
        //samtools mpileup -q 30 -B
    }
    
    public void contaminationPipe () {
        this.convertToFasta();
    }
    
    public void convertToFasta () {
        String infileS = "M:\\production\\cassava\\hmp\\speciesIntrogression\\errorFreeHiSeq_30\\3422_7332_10683_H9WVUADXX_I30572_ATCACG_R1.fastq.gz";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\contamination\\3422_7332_10683_H9WVUADXX_I30572_ATCACG_R1.fa";
        FastqChunk fc = new FastqChunk(infileS, ReadFormat.FastqGzip);
        fc.writeFasta(outfileS);
    }
    
    public void chloroplastPipe () {
        this.chloroplastContentAndEditDistance();
    }
    
    public void chloroplastContentAndEditDistance () {
        String infileDirS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\chloroplast\\sam\\";
        String distanceFigureDirS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\chloroplast\\pdf\\";
        String chloroReadsReportFileS = "M:\\pipelineTest\\cassava\\wgs\\extraDivergence\\chloroplast\\chloroplastaReads.txt";
        
        File[] fs = new File(infileDirS).listFiles();
        int fileNumber = 1;
        if (fs.length < fileNumber) fileNumber = fs.length;
        int[] readsNum = new int[fs.length];
        int[] chloroNum = new int[fs.length];
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
                                
                                //distList.add(d1);
                                if (Integer.valueOf(forwardArray[2]) == 21 || Integer.valueOf(backwardArray[2]) == 21) {
                                    chloroNum[i]++;
                                }
                                else {
                                    distList.add(d);
                                }
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
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            double[] distance = distList.toArray();
            String outfileS = fs[i].getName().replaceFirst(".gz", ".pdf");
            outfileS = new File(distanceFigureDirS, outfileS).getAbsolutePath();
            CumulativeDistribution d = new CumulativeDistribution(distance);
            d.setTitle(fs[i].getName());   
//            d.setYLim(0, 40);
//            d.setXLim(0, 0.30);
            d.setXLab("Distance");
            d.setYLab("Density");
            d.saveGraph(outfileS);
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(chloroReadsReportFileS);
            bw.write("Files\tChloroplastRatio");
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                double d1 = (double)chloroNum[i]/readsNum[i];
                bw.write(fs[i].getName()+"\t"+String.valueOf(d1));
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
}
