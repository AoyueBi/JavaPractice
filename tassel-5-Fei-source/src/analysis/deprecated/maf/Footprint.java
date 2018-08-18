/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.maf;

import format.GeneFeature;
import format.Range;
import format.RangeAttribute;
import format.Table;
import gnu.trove.list.array.TDoubleArrayList;
import graphcis.r.LineChart;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class Footprint {
    
    public Footprint () {
        //this.GeneFootprint2();
        this.GeneFootprintFigure2();
    }
    
    public void GeneFootprintFigure2() {
        String GeneFootprintFileS = "M:\\production\\maf\\wgs\\footprint\\GeneFootprintByMAF_10chrs.txt";
        String GeneFootprintFigure = "E:\\Research\\wgs_maf\\footprint\\GeneFootprintByMAF.pdf";
        TDoubleArrayList[][] mafList = null;
        String[] mafLevel = null;
        try {
            BufferedReader br = IoUtils.getTextReader(GeneFootprintFileS);
            String temp = br.readLine();
            String[] tem = temp.split("\t");
            int geneNum = Integer.valueOf(tem[0]);
            int resolution = Integer.valueOf(tem[1]);
            tem = br.readLine().split("\t");
            mafLevel = new String[tem.length-1];
            for (int i = 0; i < mafLevel.length; i++) {
                mafLevel[i] = tem[i+1];
            }
            mafList = new TDoubleArrayList[mafLevel.length][5*resolution]; //[mafClass][resolution]
            for (int i = 0; i < mafList.length; i++) {
                for (int j = 0; j < mafList[0].length; j++) mafList[i][j] = new TDoubleArrayList();
            }
            for (int i = 0; i < geneNum; i++) {
                for (int j = 0; j < mafList[0].length; j++) {
                    tem = br.readLine().split("\t");
                    for (int k = 0; k < mafList.length; k++) {
                        double d = Double.valueOf(tem[k+1]);
                        if (Double.isNaN(d)) continue;
                        mafList[k][j].add(d);
                    }
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        double[][] averMaf = new double[mafList.length][mafList[0].length];
        double[][] relaAverMaf = new double[mafList.length][mafList[0].length];
        double[] sum = new double[mafList[0].length];
        for (int i = 0; i < averMaf.length; i++) {
            for (int j = 0; j < averMaf[0].length; j++) {
                DescriptiveStatistics d = new DescriptiveStatistics(mafList[i][j].toArray());
                averMaf[i][j] = d.getMean();
                sum[i]+=averMaf[i][j];
                //sd[i] = d.getStandardDeviation();
            }
        }
        for (int i = 0; i < relaAverMaf.length; i++) {
            for (int j = 0; j < relaAverMaf[0].length; j++) {
                relaAverMaf[i][j] = averMaf[i][j]/sum[i];
            }
        }
        double[] x = new double[averMaf[0].length];
        for (int i = 0; i < x.length; i++) x[i] = i+1;
        //LineChart l = new LineChart(x, relaAverMaf, mafLevel);
        LineChart l = new LineChart(x, relaAverMaf[3], mafLevel[3]);
        l.setIfSmoothLine(false);
        //l.showOnlyLine();
        l.setXLab("Site");
        l.setYLab("Variation density along genes");
        l.setTitle("Variation distribution in genes by MARF interval");
        l.setLegendPosition(250, 0.0028);
        //l.showGraph();
        l.saveGraph(GeneFootprintFigure);
        
    }
    
    public void GeneFootprint2 () {
        int windowSize = 100;
        int streamSize = 500;
        int resolution = 100;
        String gffFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
        String GeneFootprintFileS = "M:\\production\\maf\\wgs\\footprint\\GeneFootprintByMAF.txt";
        double[] mafLow = {0, 0.001, 0.005, 0.01};
        double[] mafUp = {0.001, 0.005, 0.01, 0.5};
        GeneFeature gf = new GeneFeature(gffFileS); 
        File[] ads = AlleleDepth.getAlleleDepthFileByChrome();
        Arrays.sort(ads);
        ArrayList<float[][]> upstream = new ArrayList();
        ArrayList<float[][]> utr5 = new ArrayList();
        ArrayList<float[][]> geneBody = new ArrayList();
        ArrayList<float[][]> utr3 = new ArrayList();
        ArrayList<float[][]> downstream = new ArrayList();
        int geneNumber = 0;
        AlleleDepth cad;
        for (int k = 0; k < ads.length; k++) {
            cad = new AlleleDepth (ads[k].getAbsolutePath(), IOFileFormat.Binary);
            int currentChr = cad.getChromosome(0);
            int startGeneIndex = gf.getStartIndexOfChromosome(currentChr);
            int endGeneIndex = gf.getEndIndexOfChromosome(currentChr);
            for (int i = startGeneIndex; i< endGeneIndex; i++) {
                if (!gf.isThere5UTR(i, 0)) continue;
                if (!gf.isThere3UTR(i, 0)) continue;
                byte strand = gf.getTranscriptStrand(i, 0);
                float[][] upA = new float[mafLow.length][];
                float[][] utr5A = new float[mafLow.length][];
                float[][] geneBodyA = new float[mafLow.length][];
                float[][] utr3A = new float[mafLow.length][];
                float[][] downA = new float[mafLow.length][];
                for (int j = 0; j < mafLow.length; j++) {
                    //utr5A[j] = this.getMAFDensityAlongSequence(gf.get5UTRStart(i, 0), gf.get5UTREnd(i, 0), strand, windowSize, resolution, cad, mafLow[j], mafUp[j]);
                    //utr3A[j] = this.getMAFDensityAlongSequence(gf.get3UTRStart(i, 0), gf.get3UTREnd(i, 0), strand, windowSize, resolution, cad, mafLow[j], mafUp[j]);
                    ArrayList<Range> cdsList = gf.getCDSList(i, 0);
                    //geneBodyA[j] = this.getMAFDensityAlongSequence(cdsList.get(0).getRangeStart(), cdsList.get(cdsList.size()-1).getRangeEnd(), strand, windowSize, resolution, cad, mafLow[j], mafUp[j]);
                    geneBodyA[j] = this.getMAFDensityAlongSequence(cdsList, strand, windowSize, resolution, cad, mafLow[j], mafUp[j]);
                    //geneBodyA[j] = this.getMAFDensityAlongCDS(cdsList, strand, windowSize, resolution, cad, mafLow[j], mafUp[j]);
                    if (gf.getTranscriptStrand(i, 0) == 1) {
                        //upA[j] = this.getMAFDensityAlongSequence(gf.get5UTRStart(i, 0)-streamSize-1, gf.get5UTRStart(i, 0), strand, windowSize, resolution, cad, mafLow[j], mafUp[j]);
                        //downA[j] = this.getMAFDensityAlongSequence(gf.get3UTREnd(i, 0), gf.get3UTREnd(i, 0)+streamSize, strand, windowSize, resolution, cad, mafLow[j], mafUp[j]);
                    }
                    else {
                        //upA[j] = this.getMAFDensityAlongSequence(gf.get5UTREnd(i, 0), gf.get5UTREnd(i, 0)+streamSize, strand, windowSize, resolution, cad, mafLow[j], mafUp[j]);
                        //downA[j] = this.getMAFDensityAlongSequence(gf.get3UTRStart(i, 0)-streamSize-1, gf.get3UTRStart(i, 0), strand, windowSize, resolution, cad, mafLow[j], mafUp[j]);
                    }
                }
                upstream.add(upA);
                utr5.add(utr5A);
                geneBody.add(geneBodyA);
                utr3.add(utr3A);
                downstream.add(downA);
                geneNumber++;
                if (geneNumber%100 == 0) System.out.println(String.valueOf(geneNumber) + " genes are processed");
            }
            System.gc();
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(GeneFootprintFileS);
            bw.write(String.valueOf(geneNumber)+"\t"+String.valueOf(resolution));
            bw.newLine();
            bw.write("GeneStructure");
            for (int i = 0; i < mafLow.length; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(mafLow[i]).append("-").append(mafUp[i]);
                bw.write("\t"+sb.toString());
            }
            bw.newLine();
            for (int i = 0; i < geneNumber; i++) {
                for (int j = 0; j < resolution; j++) {
                    StringBuilder sb = new StringBuilder("Gene_"+String.valueOf(i)+"_Upstream_"+String.valueOf(j));
                    for (int k = 0; k < mafLow.length; k++) sb.append("\t").append(upstream.get(i)[k][j]);
                    bw.write(sb.toString()); bw.newLine();
                }
                for (int j = 0; j < resolution; j++) {
                    StringBuilder sb = new StringBuilder("Gene_"+String.valueOf(i)+"_5UTR_"+String.valueOf(j));
                    for (int k = 0; k < mafLow.length; k++) sb.append("\t").append(utr5.get(i)[k][j]);
                    bw.write(sb.toString()); bw.newLine();
                }
                for (int j = 0; j < resolution; j++) {
                    StringBuilder sb = new StringBuilder("Gene_"+String.valueOf(i)+"_GeneBody_"+String.valueOf(j));
                    for (int k = 0; k < mafLow.length; k++) sb.append("\t").append(geneBody.get(i)[k][j]);
                    bw.write(sb.toString()); bw.newLine();
                }
                for (int j = 0; j < resolution; j++) {
                    StringBuilder sb = new StringBuilder("Gene_"+String.valueOf(i)+"_3UTR_"+String.valueOf(j));
                    for (int k = 0; k < mafLow.length; k++) sb.append("\t").append(utr3.get(i)[k][j]);
                    bw.write(sb.toString()); bw.newLine();
                }
                for (int j = 0; j < resolution; j++) {
                    StringBuilder sb = new StringBuilder("Gene_"+String.valueOf(i)+"_Downstream_"+String.valueOf(j));
                    for (int k = 0; k < mafLow.length; k++) sb.append("\t").append(downstream.get(i)[k][j]);
                    bw.write(sb.toString()); bw.newLine();
                }
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
    public void GeneFootprintFigure () {
        String GeneFootprintFileS = "M:\\production\\maf\\wgs\\footprint\\GeneFootprint.txt";
        //String GeneFootprintFigure = "M:\\production\\maf\\wgs\\footprint\\GeneFootprint.pdf";
        String GeneFootprintFigure = "E:\\Research\\wgs_maf\\footprint\\GeneFootprint.pdf";
        String GeneFootprintSummary = "M:\\production\\maf\\wgs\\footprint\\GeneFootprint.sum.txt";
        Table t = new Table (GeneFootprintFileS);
        double[] mean = new double[t.getColumnNumber()-1];
        double[] sd = new double[t.getColumnNumber()-1];
        for (int i = 0; i < t.getColumnNumber()-1; i++) {
            double[] value = new double[t.getRowNumber()];
            for (int j = 0; j < value.length; j++) {
                value[j] = Double.valueOf(t.content[j][i+1]);
            }
            DescriptiveStatistics d = new DescriptiveStatistics(value);
            mean[i] = d.getMean();
            sd[i] = d.getStandardDeviation();
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(GeneFootprintSummary);
            bw.write("TSSSite\tMean\tSD");
            bw.newLine();
            for (int i = 0; i < mean.length; i++) {
                bw.write(t.header[i+1]+"\t"+String.valueOf(mean[i])+"\t"+String.valueOf(sd[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        double[] x = new double[mean.length];
        double[][] ys = new double[3][mean.length];
        for (int i = 0; i < x.length; i++) x[i] = Double.valueOf(t.header[i+1]);
        ys[0] = mean;
        for (int i = 0; i < x.length; i++) {
            ys[1][i] = mean[i] - sd[i] ;
            ys[2][i] = mean[i] + sd[i] ;
        }
        //String[] names = {"Mean", "Mean-SD", "Mean+SD"};
        //LineChart l = new LineChart(x, ys, names);
        String[] names = {"Mean"};
        LineChart l = new LineChart(x, mean, "Mean");
        l.setIfSmoothLine(false);
        l.showOnlyLine();
        l.setXLab("Site");
        l.setYLab("Conversation (invariant sites / bp)");
        l.setTitle("Gene conservation footprint");
        l.saveGraph(GeneFootprintFigure);
    }
    
    /**
     * @deprecated 
     */
    public void GeneFootprint () {
        int windowSize = 100;
        int streamSize = 500;
        int resolution = 100;
        String gffFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
        String GeneFootprintFileS = "M:\\production\\maf\\wgs\\footprint\\GeneFootprint.txt";
        GeneFeature gf = new GeneFeature(gffFileS); 
        int geneNumber = 0;
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            if (!gf.isThere5UTR(i, 0)) continue;
            if (!gf.isThere3UTR(i, 0)) continue;
            geneNumber++;
        }
        float[][] upstream =new float[geneNumber][];
        float[][] utr5 =new float[geneNumber][];
        float[][] cdsIntron =new float[geneNumber][];
        float[][] utr3 =new float[geneNumber][];
        float[][] downstream =new float[geneNumber][];
        AlleleDepth[] invariant = this.getInvariantSiteByChrome();
        Arrays.sort(invariant);
        int cnt = 0;
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            if (!gf.isThere5UTR(i, 0)) continue;
            if (!gf.isThere3UTR(i, 0)) continue;
            AlleleDepth cad = invariant[gf.getGeneChromosome(i)-1];
            byte strand = gf.getTranscriptStrand(i, 0);
            //utr5[cnt] = this.getInvariantFreqeuency(gf.get5UTRStart(i, 0), gf.get5UTREnd(i, 0), strand, windowSize, resolution, cad);
            //utr3[cnt] = this.getInvariantFreqeuency(gf.get3UTRStart(i, 0), gf.get3UTREnd(i, 0), strand, windowSize, resolution, cad);
            ArrayList<Range> cdsList = gf.getCDSList(i, 0);
            cdsIntron[cnt] = this.getInvariantFreqeuency(cdsList.get(0).getRangeStart(), cdsList.get(cdsList.size()-1).getRangeEnd(), strand, windowSize, resolution, cad);
            if (gf.getTranscriptStrand(i, 0) == 1) {
                //upstream[cnt] = this.getInvariantFreqeuency(gf.get5UTRStart(i, 0)-streamSize-1, gf.get5UTRStart(i, 0), strand, windowSize, resolution, cad);
                //downstream[cnt] = this.getInvariantFreqeuency(gf.get3UTREnd(i, 0), gf.get3UTREnd(i, 0)+streamSize, strand, windowSize, resolution, cad);
            }
            else {
                //upstream[cnt] = this.getInvariantFreqeuency(gf.get5UTREnd(i, 0), gf.get5UTREnd(i, 0)+streamSize, strand, windowSize, resolution, cad);
                //downstream[cnt] = this.getInvariantFreqeuency(gf.get3UTRStart(i, 0)-streamSize-1, gf.get3UTRStart(i, 0), strand, windowSize, resolution, cad);
            }
            cnt++;
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(GeneFootprintFileS);
            bw.write("Gene\\site");
            for (int i = 0; i < resolution*5; i++) {
                bw.write("\t"+String.valueOf(i+1));
            }
            bw.newLine();
            for (int i = 0; i < upstream.length; i++) {
                bw.write("Gene"+String.valueOf(i+1));
                for (int j = 0; j < upstream[0].length; j++) bw.write("\t"+String.valueOf(upstream[i][j]));
                for (int j = 0; j < utr5[0].length; j++) bw.write("\t"+String.valueOf(utr5[i][j]));
                for (int j = 0; j < cdsIntron[0].length; j++) bw.write("\t"+String.valueOf(cdsIntron[i][j]));
                for (int j = 0; j < utr3[0].length; j++) bw.write("\t"+String.valueOf(utr3[i][j]));
                for (int j = 0; j < downstream[0].length; j++) bw.write("\t"+String.valueOf(downstream[i][j]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    private float[] getMAFDensityAlongCDS (List<Range> intervals, byte strand, int windowSize, int resolution, AlleleDepth cad, double mafLow, double mafUp) {
        int size = 0;
        for (int i = 0; i < intervals.size(); i++) {
            Range r = intervals.get(i);
            size+=r.getRangeEnd()-r.getRangeStart();
        }
        float[] value = new float[size];
        int[] cdsPos = new int[size];
        int cnt = 0;
        for (int i = 0; i < intervals.size(); i++) {
            Range r = intervals.get(i);
            for (int j = r.getRangeStart(); j < r.getRangeEnd(); j++) {
                cdsPos[cnt] = j;
                cnt++;
            }
        }
        cnt = 0;
        for (int i = 0; i < intervals.size(); i++) {
            Range r = intervals.get(i);
            for (int j = r.getRangeStart(); j < r.getRangeEnd(); j++) {
                int hit = Arrays.binarySearch(cdsPos, j);
                int startIndex = hit-windowSize/2;
                if (startIndex<0) startIndex = 0;
                int endIndex = hit+windowSize/2;
                if (endIndex > cdsPos.length) endIndex = cdsPos.length;
                int actualSize = 0;
                int vCnt = 0;
                for (int k = startIndex; k < endIndex; k++) {
                    hit = cad.getSiteIndex(r.getRangeChromosome(), cdsPos[k]);
                    if (hit < 0) continue;
                    double maf = cad.getCorrectedMinorAlleleFrequency(hit);
                    if (maf>mafLow && maf<mafUp) vCnt++;
                    actualSize++;
                }
                value[cnt] = (float)((double)vCnt/actualSize);
                cnt++;
            }
        }
        if (strand == 0) {
            ArrayUtils.reverse(value);
        }
        float[] frequency = new float[resolution];
        for (int i = 0; i < frequency.length; i++) {
            int index = (int)((double)i/resolution*value.length);
            frequency[i] = value[index];
        }
        return frequency;
    }
    
    private float[] getMAFDensityAlongSequence (List<Range> intervals, byte strand, int windowSize, int resolution, AlleleDepth cad, double mafLow, double mafUp) {
        int size = 0;
        for (int i = 0; i < intervals.size(); i++) {
            Range r = intervals.get(i);
            size+=r.getRangeEnd()-r.getRangeStart();
        }
        float[] value = new float[size];
        int cnt = 0;
        for (int i = 0; i < intervals.size(); i++) {
            Range r = intervals.get(i);
            for (int j = r.getRangeStart(); j < r.getRangeEnd(); j++) {
                int actualSize = cad.getSiteNumberOfNbpInterval(cad.getChromosome(0), j, windowSize);
                value[cnt] = (float)((double)cad.getSiteNumberOfNbpIntervalByMAF(cad.getChromosome(0), j, windowSize, mafLow, mafUp)/actualSize);
                cnt++;
            }
        }
        if (strand == 0) {
            ArrayUtils.reverse(value);
        }
        float[] frequency = new float[resolution];
        for (int i = 0; i < frequency.length; i++) {
            int index = (int)((double)i/resolution*value.length);
            frequency[i] = value[index];
        }
        return frequency;
    }
    
    private float[] getMAFDensityAlongSequence (int seqStart, int seqEnd, byte strand, int windowSize, int resolution, AlleleDepth cad, double mafLow, double mafUp) {
        float[] value = new float[seqEnd-seqStart];
        if (strand == 1) {
            for (int i = 0; i < value.length; i++) {
                int actualSize = cad.getSiteNumberOfNbpInterval(cad.getChromosome(0), i+seqStart, windowSize);
                value[i] = (float)((double)cad.getSiteNumberOfNbpIntervalByMAF(cad.getChromosome(0), i+seqStart, windowSize, mafLow, mafUp)/actualSize);
            }
        }
        else {
            for (int i = 0; i < value.length; i++) {
                int actualSize = cad.getSiteNumberOfNbpInterval(cad.getChromosome(0), seqEnd-1-i, windowSize);
                value[i] = (float)((double)cad.getSiteNumberOfNbpIntervalByMAF(cad.getChromosome(0), seqEnd-1-i, windowSize, mafLow, mafUp)/actualSize);
            }
        }
        float[] frequency = new float[resolution];
        for (int i = 0; i < frequency.length; i++) {
            int index = (int)((double)i/resolution*value.length);
            frequency[i] = value[index];
        }
        return frequency;
    }
    
    private float[] getInvariantFreqeuency (int start, int end, byte strand, int windowSize, int resolution, AlleleDepth cad) {
        float[] value = new float[end-start];
        if (strand == 1) {
            for (int i = 0; i < value.length; i++) {
                value[i] = (float)((double)cad.getSiteNumberOfNbpInterval(cad.getChromosome(0), i+start, windowSize)/windowSize);
            }
        }
        else {
            for (int i = 0; i < value.length; i++) {
                value[i] = (float)((double)cad.getSiteNumberOfNbpInterval(cad.getChromosome(0), end-1-i, windowSize)/windowSize);
            }
        }
        float[] frequency = new float[resolution];
        for (int i = 0; i < frequency.length; i++) {
            int index = (int)((double)i/resolution*value.length);
            frequency[i] = value[index];
        }
        return frequency;
    }
    
    /**
     * @deprecated 
     */
    public void TSSFootprintFigure () {
        String TSSFootprintFileS = "M:\\production\\maf\\wgs\\footprint\\TSSFootprint.txt";
        //String TSSFootprintFigure = "M:\\production\\maf\\wgs\\footprint\\TSSFootprint.pdf";
        String TSSFootprintFigure = "E:\\Research\\wgs_maf\\footprint\\TSSFootprint.pdf";;
        String TSSFootprintSummary = "M:\\production\\maf\\wgs\\footprint\\TSSFootprint.sum.txt";
        Table t = new Table (TSSFootprintFileS);
        double[] mean = new double[t.getColumnNumber()-1];
        double[] sd = new double[t.getColumnNumber()-1];
        for (int i = 0; i < t.getColumnNumber()-1; i++) {
            double[] value = new double[t.getRowNumber()];
            for (int j = 0; j < value.length; j++) {
                value[j] = Double.valueOf(t.content[j][i+1]);
            }
            DescriptiveStatistics d = new DescriptiveStatistics(value);
            mean[i] = d.getMean();
            sd[i] = d.getStandardDeviation();
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(TSSFootprintSummary);
            bw.write("TSSSite\tMean\tSD");
            bw.newLine();
            for (int i = 0; i < mean.length; i++) {
                bw.write(t.header[i+1]+"\t"+String.valueOf(mean[i])+"\t"+String.valueOf(sd[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        double[] x = new double[mean.length];
        double[][] ys = new double[3][mean.length];
        for (int i = 0; i < x.length; i++) x[i] = Double.valueOf(t.header[i+1]);
        ys[0] = mean;
        for (int i = 0; i < x.length; i++) {
            ys[1][i] = mean[i] - sd[i] ;
            ys[2][i] = mean[i] + sd[i] ;
        }
        //String[] names = {"Mean", "Mean-SD", "Mean+SD"};
        //LineChart l = new LineChart(x, ys, names);
        String[] names = {"Mean"};
        LineChart l = new LineChart(x, mean, "Mean");
        l.setIfSmoothLine(false);
        l.showOnlyLine();
        l.setXLab("Site");
        l.setYLab("Conversation (invariant sites / bp)");
        l.setTitle("TSS conservation footprint");
        l.saveGraph(TSSFootprintFigure);
    }
    
    /**
     * @deprecated 
     */
    public void TSSFoorprint () {
        int halfintervalSize = 500;
        int windowSize = 200;
        String utr5FileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\5UTR.ra.txt";
        String TSSFootprintFileS = "M:\\production\\maf\\wgs\\footprint\\TSSFootprint.txt";
        RangeAttribute utr = new RangeAttribute(utr5FileS, IOFileFormat.Text);
        AlleleDepth[] invariant = this.getInvariantSiteByChrome();
        Arrays.sort(invariant);
        float[][] frequency = new float[utr.getRangeNumber()][2*halfintervalSize+1];
        for (int i = 0; i < utr.getRangeNumber(); i++) {
            AlleleDepth cad = invariant[utr.getRangeChromosome(i)-1];
            int cnt = 0;
            if (utr.getStrand(i) == 1) {
                for (int j = utr.getRangeStart(i)-halfintervalSize-1; j < utr.getRangeStart(i)+halfintervalSize; j++) {
                    frequency[i][cnt] = (float)((double)(cad.getSiteNumberOfNbpInterval(utr.getRangeChromosome(i), j, windowSize))/windowSize);
                    cnt++;
                }
            }
            else {
                for (int j = utr.getRangeEnd(i)-1+halfintervalSize; j > utr.getRangeEnd(i)-2-halfintervalSize; j--) {
                    frequency[i][cnt] = (float)((double)(cad.getSiteNumberOfNbpInterval(utr.getRangeChromosome(i), j, windowSize))/windowSize);
                    cnt++;
                }
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(TSSFootprintFileS);
            bw.write("Gene\\site");
            for (int i = 0; i < frequency[0].length; i++) {
                bw.write("\t"+String.valueOf(-halfintervalSize+i));
            }
            bw.newLine();
            for (int i = 0; i < frequency.length; i++) {
                bw.write("Gene"+String.valueOf(i+1));
                for (int j=  0; j < frequency[i].length; j++) {
                    bw.write("\t"+String.valueOf(frequency[i][j]));
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

    private AlleleDepth[] getInvariantSiteByChrome () {
        String infileDirS = "M:\\production\\maf\\wgs\\invariantSite";
        File[] fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        ArrayList<AlleleDepth> adList = new ArrayList();
        fsList.parallelStream().forEach(e -> {
            adList.add(new AlleleDepth(e.getAbsolutePath(), IOFileFormat.Binary));
        });
        AlleleDepth[] ads = adList.toArray(new AlleleDepth[adList.size()]);
        Arrays.sort(ads);
        return ads;
    }
}
