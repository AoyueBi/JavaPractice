/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.maf;


import format.Fasta;
import format.GeneFeature;
import format.Range;
import format.Ranges;
import format.Table;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import graphcis.r.BoxPlot;
import graphcis.r.DensityPlot;
import graphcis.r.DensityPlotMultiClass;
import graphcis.r.LineChart;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import org.apache.commons.math3.stat.inference.TTest;
import utils.FArrayUtils;
import utils.IoUtils;

/**
 *
 * @author Fei Lu
 */
class InvariantSiteValidation {
    
    public InvariantSiteValidation () {
        //this.bacteriaConservedSite();
        //this.codonBaseUnderMafCut();
        //this.codonBaseProportionByMAF();
        this.codonBaseproportionByMAFFigure();
    }
    
    
    public void codonBaseproportionByMAFFigure () {
        String infileS = "E:\\Research\\wgs_maf\\invariantSite_validation\\codonMAF.txt";
        String infileS2 = "E:\\Research\\wgs_maf\\invariantSite_validation\\codonMAF_density.txt";
        String pdfAllFileS = "E:\\Research\\wgs_maf\\invariantSite_validation\\codonMAF_all.pdf";
        String pdfRareFileS = "E:\\Research\\wgs_maf\\invariantSite_validation\\codonMAF_rare.pdf";
        String pdfBoxFileS = "E:\\Research\\wgs_maf\\invariantSite_validation\\codonMAF_box.pdf";
        String pdfDensityFileS = "E:\\Research\\wgs_maf\\invariantSite_validation\\codonMAF_density.pdf";
        String tTestFileS = "E:\\Research\\wgs_maf\\invariantSite_validation\\codonMAF_TTest.txt";
        Table t = new Table (infileS);
        double[] averageMaf = new double[t.getRowNumber()];
        double[][] mafRatio = new double[3][t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            averageMaf[i] = t.getDoubleValue(i, 1);
            for (int j = 0; j < mafRatio.length; j++) {
                mafRatio[j][i] = t.getDoubleValue(i, j+5);
            }
        }
        String[] variableNames = {"1st Base", "2nd Base", "3rd Base"};
        LineChart l = new LineChart(averageMaf, mafRatio, variableNames);
        l.setTitle("Minor allele read frequency proportion of codon bases");
        l.setXLab("Minor allele read frequency");
        l.setYLab("Proportion of 3 codon bases");
        l.setLegendPosition(0.35, 1.0);
        l.saveGraph(pdfAllFileS);
        l.setXLim(0, 0.005);
        l.setLegendPosition(0.004, 1.0);
        l.saveGraph(pdfRareFileS);
        
        t = new Table (infileS2);
        double[][] rareMaf = new double[3][t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            rareMaf[0][i] = t.getDoubleValue(i, 1);
            rareMaf[1][i] = t.getDoubleValue(i, 0);
            rareMaf[2][i] = t.getDoubleValue(i, 2);
        }
        variableNames[0] = "2nd Base";
        variableNames[1] = "1st Base";
        BoxPlot b = new BoxPlot(rareMaf, variableNames);
        b.setTitle("Rare allele distribution in codon bases");
        b.setYLab("Minor allele read frequency");
        b.saveGraph(pdfBoxFileS);
        DensityPlotMultiClass d = new DensityPlotMultiClass(rareMaf, variableNames);
        d.setTitle("Rare allele distribution in codon bases");
        d.setXLab("Minor allele read frequency");
        d.setYLab("Density");
        d.saveGraph(pdfDensityFileS);
        try {
            TTest test = new TTest();
            BufferedWriter bw  = IoUtils.getTextWriter(tTestFileS);
            bw.write("TTest");
            for (int i = 0; i < variableNames.length; i++) bw.write("\t"+variableNames[i]);
            bw.newLine();
            for (int i = 0; i < variableNames.length; i++) {
                bw.write(variableNames[i]);
                for (int j = 0; j < variableNames.length; j++) {
                    bw.write("\t"+String.valueOf(test.pairedTTest(rareMaf[i], rareMaf[j])));
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
                    
        }
        catch (Exception e) {
            e.toString();
        }
    }
    
    public void codonBaseProportionByMAF () {
        int groupNumber = 500;
        String geneFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
        String outfileS = "E:\\Research\\wgs_maf\\invariantSite_validation\\codonMAF.txt";
        String outfileS2 = "E:\\Research\\wgs_maf\\invariantSite_validation\\codonMAF_density.txt";
        int sampleSize = 30000;
        double mafCutoff = 0.01;
        GeneFeature gf = new GeneFeature(geneFileS);
        ArrayList<Range>[] lists = new ArrayList[3];
        for (int i = 0; i < lists.length; i++) lists[i] = new ArrayList();
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            ArrayList<Range> cds = gf.getCDSList(i, 0);
            int sum = 0;
            for (int j = 0; j < cds.size(); j++) {
                sum+=(cds.get(j).getRangeEnd()-cds.get(j).getRangeStart());
            }
            if (sum%3!=0) continue;
            int cnt = 0;
            if (gf.getTranscriptStrand(i, 0) == 0) {
                for (int j = 0; j < cds.size(); j++) {
                    for (int k = cds.get(j).getRangeStart(); k < cds.get(j).getRangeEnd(); k++) {
                        int remainder = cnt%3;
                        Range r = new Range(gf.getTranscriptChromosome(i, 0), k, k+1);
                        if (remainder == 0) {
                            lists[2].add(r);
                        }
                        else if (remainder == 1) {
                            lists[1].add(r);
                        }
                        else {
                            lists[0].add(r);
                        }
                        cnt++;
                    }
                }
            }
            else if (gf.getTranscriptStrand(i, 0) == 1){
                for (int j = 0; j < cds.size(); j++) {
                    for (int k = cds.get(j).getRangeStart(); k < cds.get(j).getRangeEnd(); k++) {
                        int remainder = cnt%3;
                        Range r = new Range(gf.getTranscriptChromosome(i, 0), k, k+1);
                        if (remainder == 0) {
                            lists[0].add(r);
                        }
                        else if (remainder == 1) {
                            lists[1].add(r);
                        }
                        else {
                            lists[2].add(r);
                        }
                        cnt++;
                    }
                }
            }
        }
        Ranges[] codonPos = new Ranges[lists.length];
        for (int i = 0; i < lists.length; i++) {
            codonPos[i] = new Ranges(lists[i], "CodonPos"+String.valueOf(i+1));
        }
        System.out.println("There are "+String.valueOf(codonPos[0].getRangeNumber())+ " AA");
        TDoubleArrayList[] codonMAFList = new TDoubleArrayList[3];
        double[][] codonMAF = new double[3][];
        for (int i = 0; i < codonMAFList.length; i++) codonMAFList[i] = new TDoubleArrayList();
        AlleleDepth[] ads = AlleleDepth.getAlleleDepthByChrome();
        int avaCnt = 0;
        for (int i = 0; i < ads.length; i++) {
            int currentChr = ads[i].getChromosome(0);
            for (int j = 0; j < codonPos.length; j++) {
                for (int k = 0; k < codonPos[j].getRangeNumber(); k++) {
                    if (codonPos[j].getRangeChromosome(k) != currentChr) continue;
                    int index = ads[i].getSiteIndex(codonPos[j].getRangeChromosome(k), codonPos[j].getRangeStart(k));
                    if (index < 0) continue;
                    double ma = ads[i].getCorrectedMinorAlleleFrequency(index);
                    if (ma == 0) continue;
                    avaCnt++;
                    codonMAFList[j].add(ma);
                }
            }
        }
        System.out.println("Available codon number:\t " + String.valueOf(avaCnt));
        for (int i = 0; i < codonMAF.length; i++) {
            codonMAF[i] = codonMAFList[i].toArray();
            Arrays.sort(codonMAF[i]);
        }
        double[][] mafMean = new double[codonMAF.length][groupNumber];
        double[][] mafRatio = new double[codonMAF.length][groupNumber];
        for (int i = 0; i < codonMAF.length; i++) {
            int[][] index = FArrayUtils.getSubsetsIndicesBySubsetNumber(codonMAF[i].length, groupNumber);
            for (int j = 0; j < groupNumber; j++) {
                double v = 0;
                for (int k = index[j][0]; k < index[j][1]; k++) {
                    v+=(codonMAF[i][k]/(index[j][1]-index[j][0]));
                }
                mafMean[i][j] = v;
            }  
        }
        for (int i = 0; i < groupNumber; i++) {
            double sum = 0;
            for (int j = 0; j < mafMean.length; j++) {
                sum+=mafMean[j][i];
            }
            if (sum == 0) {
                for (int j = 0; j < mafRatio.length; j++) {
                    mafRatio[j][i] = (double)1/mafRatio.length;
                }
            }
            else {
                for (int j = 0; j < mafRatio.length; j++) {
                    mafRatio[j][i] = mafMean[j][i]/sum;
                }
            }
            
        }
        double[] averageMaf = new double[mafMean[0].length];
        for (int i = 0; i < averageMaf.length; i++) {
            double v = 0;
            for (int j = 0; j < mafMean.length; j++) {
                v+=mafMean[j][i];
            }
            averageMaf[i] = v/mafMean.length;
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("MAF_Quantile\tMAF_Mean\tMAF_1stBase\tMAF_2ndBase\tMAF_3rdBase\tRatio_1stBase\tRatio_2ndBase\tRatio_3rdBase");
            bw.newLine();
            for (int i = 0; i < mafMean[0].length; i++) {
                bw.write(String.valueOf(i+1)+"\t");
                bw.write(String.valueOf(averageMaf[i]));
                for (int j = 0; j < mafMean.length; j++) {
                    bw.write("\t"+mafMean[j][i]);
                }
                for (int j = 0; j < mafRatio.length; j++) {
                    bw.write("\t"+mafRatio[j][i]);
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        TDoubleArrayList[] rareCodonMafList = new TDoubleArrayList[3];
        double[][] rareCodonMaf = new double[3][];
        for (int i = 0; i < rareCodonMafList.length; i++) rareCodonMafList[i] = new TDoubleArrayList();
        for (int i = 0; i < codonMAF.length; i++) {
            for (int j = 0; j < codonMAF[i].length; j++) {
                if (codonMAF[i][j] > mafCutoff) continue;
                rareCodonMafList[i].add(codonMAF[i][j]);
            }
            rareCodonMaf[i] = rareCodonMafList[i].toArray();
            FArrayUtils.shuffleArray(rareCodonMaf[i]);
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS2);
            bw.write("1stBase\t2ndBase\t3rdBase");
            bw.newLine();
            for (int i = 0; i < sampleSize; i++) {
                StringBuilder sb = new StringBuilder();
                for (int j = 0; j < rareCodonMaf.length; j++) {
                    sb.append(rareCodonMaf[j][i]+"\t");
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
    /**
     * @deprecated 
     */
    public void codonBaseUnderMafCut () {
        double mafCut = 0.001;
        String geneFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
        GeneFeature gf = new GeneFeature(geneFileS);
        ArrayList<Range>[] lists = new ArrayList[3];
        for (int i = 0; i < lists.length; i++) lists[i] = new ArrayList();
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            ArrayList<Range> cds = gf.getCDSList(i, 0);
            int sum = 0;
            for (int j = 0; j < cds.size(); j++) {
                sum+=(cds.get(j).getRangeEnd()-cds.get(j).getRangeStart());
            }
            if (sum%3!=0) continue;
            int cnt = 0;
            if (gf.getTranscriptStrand(i, 0) == 0) {
                for (int j = 0; j < cds.size(); j++) {
                    for (int k = cds.get(j).getRangeStart(); k < cds.get(j).getRangeEnd(); k++) {
                        int remainder = cnt%3;
                        Range r = new Range(gf.getTranscriptChromosome(i, 0), k, k+1);
                        if (remainder == 0) {
                            lists[2].add(r);
                        }
                        else if (remainder == 1) {
                            lists[1].add(r);
                        }
                        else {
                            lists[0].add(r);
                        }
                        cnt++;
                    }
                }
            }
            else if (gf.getTranscriptStrand(i, 0) == 1){
                for (int j = 0; j < cds.size(); j++) {
                    for (int k = cds.get(j).getRangeStart(); k < cds.get(j).getRangeEnd(); k++) {
                        int remainder = cnt%3;
                        Range r = new Range(gf.getTranscriptChromosome(i, 0), k, k+1);
                        if (remainder == 0) {
                            lists[0].add(r);
                        }
                        else if (remainder == 1) {
                            lists[1].add(r);
                        }
                        else {
                            lists[2].add(r);
                        }
                        cnt++;
                    }
                }
            }
        }
        Ranges[] codonPos = new Ranges[lists.length];
        for (int i = 0; i < lists.length; i++) {
            codonPos[i] = new Ranges(lists[i], "CodonPos"+String.valueOf(i+1));
        }
        System.out.println("There are "+String.valueOf(codonPos[0].getRangeNumber())+ " AA");
        AlleleDepth[] ads = AlleleDepth.getAlleleDepthByChrome();
        int[] existCnt = new int[3];
        int[] invariantCnt = new int[3];
        int cntOnChr = 0;
        for (int i = 0; i < ads.length; i++) {
            int currentChr = ads[i].getChromosome(0);
            for (int j = 0; j < codonPos.length; j++) {
                cntOnChr = 0;
                for (int k = 0; k < codonPos[j].getRangeNumber(); k++) {
                    if (codonPos[j].getRangeChromosome(k) != currentChr) continue;
                    cntOnChr++;
                    int index = ads[i].getSiteIndex(codonPos[j].getRangeChromosome(k), codonPos[j].getRangeStart(k));
                    if (index < 0) continue;
                    existCnt[j]++;
                    if (ads[i].getCorrectedMinorAlleleFrequency(index) > mafCut) continue;
                    //if (ads[i].getMinorAlleleFrequency(index) < mafCut) continue;
                    invariantCnt[j]++;
                }
            }
        }
        System.out.println("codon number on chromosomes:\t " + String.valueOf(cntOnChr));
        for (int i = 0; i < 3; i++) {
            System.out.println(String.valueOf(i+1)+"st in codon:\tExistNumber: "+String.valueOf(existCnt[i])+"\tInvariantsNumber: "+String.valueOf(invariantCnt[i])+"\tRatio: "+String.valueOf((double)invariantCnt[i]/existCnt[i]));
        }
    }
    
    public void bacteriaConservedSite () {
        String conservedSiteFileS = "E:\\Research\\wgs_maf\\invariantSite_validation\\ConservedSitesAgPv3.txt";
        AlleleDepth[] ads = InvariantSiteDistribution.getInvariantSiteByChrome();
        ArrayList<Range> rList = new ArrayList();
        String temp = null;
        try {
            BufferedReader br = IoUtils.getTextReader(conservedSiteFileS);
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split(",");
                for (int i = 0; i < tem.length; i++) {
                    tem[i] = tem[i].replaceAll("\\D", "");
                }
                if (tem[0].isEmpty()) continue;
                int chr = Integer.valueOf(tem[0]);
                if (chr < 1 || chr > 10) continue;
                for (int i = 0; i < tem.length-3; i++) {
                    int pos = Integer.valueOf(tem[i+3]);
                    Range r = new Range(chr,pos,pos+1);
                    rList.add(r);
                }
            }
        }
        catch (Exception e) {
            System.out.println(temp);
            e.printStackTrace();
        }
        Ranges r = new Ranges (rList, "bacteriaConservedSite");
        int totalSite = r.getRangeNumber();
        int cnt = 0;
        for (int i = 0; i < r.getRangeNumber(); i++) {
            AlleleDepth cad = ads[r.getRangeChromosome(i)-1];
            int index = cad.getSiteIndex(r.getRangeChromosome(i), r.getRangeStart(i));
            if (index < 0) continue;
            cnt++;
        }
        System.out.println("Total conserved sites:\t" +String.valueOf(totalSite));
        System.out.println("Total invariant sites:\t" +String.valueOf(cnt));
        System.out.println("Ratio:\t" + String.valueOf((double)cnt/totalSite));
    }
    
}
