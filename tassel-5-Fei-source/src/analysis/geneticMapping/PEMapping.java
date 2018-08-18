/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package analysis.geneticMapping;

import format.Bins;
import format.Fasta;
import format.ShortreadAlignment;
import format.Table;
import graphcis.r.DensityPlot;
import graphcis.r.DensityPlotMultiClass;
import graphcis.r.Histogram;
import graphcis.r.HistogramRF;
import graphcis.r.ScatterPlot;
import graphcis.r.ScatterPlotMultiClass;
import java.io.BufferedWriter;
import java.util.ArrayList;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.map.TagsOnGeneticMap;
import net.maizegenetics.dna.tag.PETagCounts;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author Fei Lu <fl262@cornell.edu>
 */
public class PEMapping {
    
    public PEMapping () {
        this.PEstatistics();
        //this.outputContigFasta();
        //this.convertAlignment();
        //this.checkCorrectAlignment();
        //this.contigLengthFigure();
        //this.mkAccuracyVsLength();
    }
    
    public void mkAccuracyVsLength () {
        String checkFileS = "E:\\Research\\geneticMapping\\pe\\contigGeneticCheck_k10.txt";
        String accuracyLengthFileS = "E:\\Research\\geneticMapping\\pe\\acuLen.txt";
        String plotFile = "E:\\Research\\geneticMapping\\pe\\acuLen.pdf";
        Table t = new Table (checkFileS);
        int[] cor = new int[t.getRowNumber()];
        double[] d = new double[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            cor[i] = Integer.valueOf(t.content[i][1]);
            d[i] = Double.valueOf(t.content[i][6]);
        }
        Bins b =new Bins (0,500,10,cor,d);
        ArrayList<Double> xList = new ArrayList();
        ArrayList<Double> yList = new ArrayList();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(accuracyLengthFileS);
            bw.write("Length\tAlignmentNum\tAccuracy");
            bw.newLine();
            for (int i = 0; i < b.getBinNum(); i++) {
                double aver = b.getBinAverage(i);
                bw.write(String.valueOf(i*10)+"\t"+String.valueOf(b.getBinValues(i).length)+"\t"+String.valueOf(aver));
                bw.newLine();
                if (Double.isNaN(aver)) continue;
                xList.add((double)i*10);
                yList.add(aver);
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        double[] x = new double[xList.size()];
        double[] y = new double[yList.size()];
        for (int i = 0; i < x.length; i++) {
            x[i] = xList.get(i);
            y[i] = yList.get(i);
        }
        ScatterPlot s = new ScatterPlot(x,y);
        s.setYLim(0, 1);
        s.setXLim(0, 500);
        s.setPlottingCharacter(16);
        s.setSlideMode();
        s.setXLab("Sequence length");
        s.setYLab("Alignment accuracy");
        s.setTitle("");
        s.saveGraph(plotFile);
    }
    
    public void contigLengthFigure () {
        String peTagCountFileS = "M:\\production\\pe\\mergePETagCounts\\merge.con.pe.cnt";
        String contigLengthPdf = "E:\\Research\\geneticMapping\\pe\\contigLengthHisto.pdf";
        PETagCounts ptc = new PETagCounts (peTagCountFileS, FilePacking.Byte);
        int size = 500000;
        System.out.println(String.valueOf(size)+" contigs");
        int cnt = 0; 
        int[] lengths = new int[size];
        for (int i = 0; i < ptc.getTagCount(); i++) {
            if (ptc.getContigLength(i) == 0) continue;
            lengths[cnt] = ptc.getContigLength(i);
            cnt++;
            if (cnt == size) break;
        }
        HistogramRF h = new HistogramRF (lengths);
        h.setXLab("Contig length");
        h.setTitle("");
        h.setSlideMode();
        h.saveGraph(contigLengthPdf);
    }
    
    public void checkCorrectAlignment () {
        //String simpleAlignment = "E:\\Research\\geneticMapping\\pe\\contigAlign_M.sa";
        String simpleAlignment = "E:\\Research\\geneticMapping\\pe\\contigAlign_k10.sa";
        String contigFileS = "E:\\Research\\geneticMapping\\pe\\contig.fa";
        String togmFileS = "M:\\pav\\PhyGenMapping\\v1.togm.txt";
        String outputFileS = "E:\\Research\\geneticMapping\\pe\\contigGeneticCheck.txt";
        TagsOnGeneticMap togm = new TagsOnGeneticMap(togmFileS, FilePacking.Text);
        ShortreadAlignment sa = new ShortreadAlignment(simpleAlignment, IOFileFormat.Binary);
        sa.sortByQuery();
        Fasta f = new Fasta(contigFileS);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outputFileS);
            bw.write("FastqName\tLength\tGChr\tGPos\tPChr\tPPos\t10MCorrent");
            bw.newLine();
            for (int i = 0; i < f.getSeqNumber(); i++) {
                String name = f.getName(i);
                String seq = (f.getSeq(i)+"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").substring(0, 64);
                long[] t = BaseEncoder.getLongArrayFromSeq(seq);
                int index1 = togm.getTagIndex(t);
                if (index1 < 0) continue;
                int index2 = sa.getAlignmentStartIndexByQuery(name);
                if (sa.getStartPos(index2) < 0) continue;
                bw.write(name+"\t"+String.valueOf(f.getSeq(i).length())+"\t"+String.valueOf(togm.getGChr(index1))+"\t"+String.valueOf(togm.getGPos(index1))+"\t");
                bw.write(sa.getHit(index2)+"\t"+String.valueOf(sa.getStartPos(index2))+"\t");
                if (togm.getGChr(index1) != Integer.valueOf(sa.getHit(index2))) {
                    bw.write("0");
                }
                else {
                    if (Math.abs(sa.getStartPos(index2)-togm.getGPos(index1)) < 10000000) {
                        bw.write("1");
                    }
                    else {
                        bw.write("0");
                    }
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }
    
    public void convertAlignment () {
        //String samFileS = "E:\\Research\\geneticMapping\\pe\\contig.sam";
        //String simpleAlignment = "E:\\Research\\geneticMapping\\pe\\contigAlign.sa";
        String samFileS = "E:\\Research\\geneticMapping\\pe\\contig_k10.sam";
        String simpleAlignment = "E:\\Research\\geneticMapping\\pe\\contigAlign_k10.sa";
        ShortreadAlignment sa = new ShortreadAlignment();
        sa.readFromBowtie2(samFileS);
        sa.writeSimpleAlignment(simpleAlignment, IOFileFormat.Binary);
    }
    
    public void outputContigFasta () {
        String peTagCountFileS = "M:\\production\\pe\\mergePETagCounts\\merge.con.pe.cnt";
        String contigFileS = "E:\\Research\\geneticMapping\\pe\\contig.fa";
        PETagCounts ptc = new PETagCounts (peTagCountFileS, FilePacking.Byte);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(contigFileS);
            int cnt = 0;
            for (int i = 0; i < ptc.getTagCount(); i++) {
                if (ptc.getContigLength(i) == 0) continue;
                bw.write(">"+String.valueOf(cnt));
                bw.newLine();
                bw.write(BaseEncoder.getSequenceFromLong(ptc.getContig(i)).substring(0, ptc.getContigLength(i)));
                bw.newLine();
                cnt++;
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void PEstatistics () {
        String peTagCountFileS = "M:\\production\\pe\\mergePETagCounts\\merge.con.pe.cnt";
        String staFileS = "E:\\Research\\geneticMapping\\pe\\sta.txt";
        PETagCounts ptc = new PETagCounts (peTagCountFileS, FilePacking.Byte);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(staFileS);
            bw.write("Tag count is " + ptc.getTagCount());
            bw.newLine();
            bw.write("Read count is " + ptc.getTotalReadCount());
            bw.newLine();
            bw.write("Contig count is " + ptc.getContigCount());
            bw.newLine();
            int min = 500;
            int max = 0;
            int contigMin = 1000;
            int contigMax = 0;

            for (int i = 0; i < ptc.getTagCount(); i++) {
                if (ptc.getTagFLength(i) > max) max = ptc.getTagFLength(i);
                if (ptc.getTagBLength(i) > max) max = ptc.getTagBLength(i);
                if (ptc.getTagFLength(i) < min) min = ptc.getTagFLength(i);
                if (ptc.getTagBLength(i) < min) min = ptc.getTagBLength(i);
                if (ptc.getContigLength(i) == 0) continue;
                if (ptc.getContigLength(i) > contigMax) contigMax = ptc.getContigLength(i);
                if (ptc.getContigLength(i) < contigMin) contigMin = ptc.getContigLength(i);
            }
            bw.write("Min tag length is " + String.valueOf(min));
            bw.newLine();
            bw.write("Max tag length is " + String.valueOf(max));
            bw.newLine();
            bw.write("Min contig length is " + String.valueOf(contigMin));
            bw.newLine();
            bw.write("Max contig length is " + String.valueOf(contigMax));
            bw.newLine();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
