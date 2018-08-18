/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.cml247;

import format.Fasta;
import format.Table;
import graphcis.r.DensityPlot;
import graphcis.r.Histogram;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.HashSet;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class QualityControl {
    
    public QualityControl () {
        //this.checkRatioOnMainChr();
        this.mkReport();
    }
    
    public void mkReport () {
        int minAnchorNum = 0;
        String fastaFileS = "E:\\Research\\cml247\\fasta\\cml247.fas";
        String infileS = "E:\\Research\\cml247\\assemblyEvaluation\\agpv3\\result\\scaffoldQuality_0bp_2anchors.txt";
        String reportFileS = "E:\\Research\\cml247\\assemblyEvaluation\\agpv3\\report\\report.txt";
        
//        String fastaFileS = "N:\\panGenome\\W22\\W22_contigs.fasta";
//        String infileS = "M:\\production\\panGenome\\w22\\evaluation\\contig\\result\\scaffoldQuality_0bp_2anchors.txt";
//        String reportFileS = "M:\\production\\panGenome\\w22\\evaluation\\contig\\report\\report.txt";
        
//        String fastaFileS = "N:\\panGenome\\W22\\W22_scaffolds.fasta";
//        String infileS = "M:\\production\\panGenome\\w22\\evaluation\\scaffold\\result\\scaffoldQuality_0bp_2anchors.txt";
//        String reportFileS = "M:\\production\\panGenome\\w22\\evaluation\\scaffold\\report\\report.txt";
        
        Fasta f = new Fasta(fastaFileS);
        Table t = new Table(infileS);
        int actualMinAnchorNum = Integer.MAX_VALUE;
        int actualMaxAnchorNum = 0;
        int evaluatedScaffoldNum = 0;
        int evaluatedLength = 0;
        ArrayList<Double> dList = new ArrayList();
        int multiRegionScaffoldNum = 0;
        int multiRegionScaffoldLength = 0;
        int multiChromScaffoldNum = 0;
        int multiChromScaffoldLength = 0;
        for (int i = 0; i < t.getRowNumber(); i++) {
            int n = Integer.valueOf(t.content[i][3]);
            if (n < minAnchorNum) continue;
            if (n < actualMinAnchorNum) actualMinAnchorNum = n;
            if (n > actualMaxAnchorNum) actualMaxAnchorNum = n;
            if (Integer.valueOf(t.content[i][7]) > 1) {
                String[] temp = t.content[i][8].split(";");
                HashSet<Integer> hs = new HashSet();
                for (int j = 0; j < temp.length; j++) {
                    hs.add(Integer.valueOf(temp[j].split("-")[0]));
                }
                if (hs.size()>1) {
                    multiChromScaffoldNum++;
                    multiChromScaffoldLength+=Integer.valueOf(t.content[i][2]);
                }
                else {
                    
                }
                multiRegionScaffoldNum++;
                multiRegionScaffoldLength+=Integer.valueOf(t.content[i][2]);
            }
            dList.add(Double.valueOf(t.content[i][6]));
            evaluatedScaffoldNum++;
            evaluatedLength+=Integer.valueOf(t.content[i][2]);
        }
        double[] qvalue = new double[dList.size()];
        for (int i = 0; i < qvalue.length; i++) {
            qvalue[i] = dList.get(i);
        }
        DescriptiveStatistics d = new DescriptiveStatistics(qvalue);
        long tLength = f.getTotalSeqLength();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(reportFileS);
            bw.write(fastaFileS);bw.newLine();
            bw.write("Scaffold number:\t " + String.valueOf(f.getSeqNumber()));bw.newLine();
            bw.write("Total length:\t "+ String.valueOf(tLength));bw.newLine();
            bw.write("L50 length:\t"+String.valueOf(f.getL50()));bw.newLine();
            bw.write("Min anchor number:\t " + String.valueOf(actualMinAnchorNum));bw.newLine();
            bw.write("Max anchor number:\t " + String.valueOf(actualMaxAnchorNum));bw.newLine();
            bw.write("Quality value mean:\t" + String.valueOf(d.getMean()));bw.newLine();
            bw.write("Quality value min:\t" + String.valueOf(d.getMin()));bw.newLine();
            bw.write("Quality value max:\t" + String.valueOf(d.getMax()));bw.newLine();
            bw.write("Quality value standard deviation:\t" + String.valueOf(d.getStandardDeviation()));bw.newLine();
            bw.write("Number of scaffold from multiple region:\t" + String.valueOf(multiRegionScaffoldNum));bw.newLine();
            bw.write("Length of scaffold from multiple region:\t" + String.valueOf(multiRegionScaffoldLength)+"\t"+String.valueOf((double)multiRegionScaffoldLength/tLength));bw.newLine();
            bw.write("Number of scaffold from multiple chromosome:\t" + String.valueOf(multiChromScaffoldNum));bw.newLine();
            bw.write("Length of scaffold from multiple chromosome:\t" + String.valueOf(multiChromScaffoldLength)+"\t"+String.valueOf((double)multiChromScaffoldLength/tLength));bw.newLine();
            bw.write("Evaluated scaffold number:\t"+String.valueOf(qvalue.length));bw.newLine();
            bw.write("Evaluated scaffold length:\t"+String.valueOf(evaluatedLength));bw.newLine();
            bw.write("Evaluated scaffold ratio(eva/all):\t"+String.valueOf((double)evaluatedLength/tLength));bw.newLine();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            
            e.printStackTrace();
        }
    }
    
    public void checkRatioOnMainChr () {
        int minAnchorNum = 20;
        String infileS = "E:\\Research\\cml247\\assemblyEvaluation\\agpv3\\result\\scaffoldQuality_0bp_2anchors.txt";
        String pdfFileS = "E:\\Research\\cml247\\assemblyEvaluation\\agpv3\\report\\20Anchors.pdf";
        
//        String infileS = "M:\\production\\panGenome\\w22\\evaluation\\contig\\result\\scaffoldQuality_0bp_2anchors.txt";
//        String pdfFileS = "M:\\production\\panGenome\\w22\\evaluation\\contig\\report\\20Anchors.pdf";
        
//        String infileS = "M:\\production\\panGenome\\w22\\evaluation\\scaffold\\result\\scaffoldQuality_0bp_2anchors.txt";
//        String pdfFileS = "M:\\production\\panGenome\\w22\\evaluation\\scaffold\\report\\20Anchors.pdf";
        
        Table t = new Table (infileS);
        ArrayList<Double> d = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (Integer.valueOf(t.content[i][3]) < minAnchorNum) continue;
            //if (Integer.valueOf(t.content[i][2]) < 20000) continue;
            d.add(Double.valueOf(t.content[i][6]));
        }
        double[] value = new double[d.size()];
        for (int i = 0; i < value.length; i++) value[i] = d.get(i);
        DensityPlot p = new DensityPlot(value);
        Histogram h = new Histogram(value);
        h.setBreakNumber(20);
        h.setTitle("Assembly quality evaluation using genetic anchors");
        h.setXLab("Proportion of anchors on the main chromosome");
        h.setYLab("Number of CML247 contigs (with > 20 anchors)");
        //h.setYLim(0, 35000);
        h.setLargerSize();
        h.setLargerSize();
        h.saveGraph(pdfFileS);
    }
}
