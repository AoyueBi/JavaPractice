/*
* To change this license header, choose License Headers in Project Properties.
* To change this template file, choose Tools | Templates
* and open the template in the editor.
*/
package analysis.maf;

import format.GeneFeature;
import format.Range;
import format.Ranges;
import format.Table;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import utils.FArrayUtils;
import utils.FStringUtils;
import utils.IoUtils;
/**
 *
 * @author fl262
 */
class Expression {
    public Expression () {
        //this.processRawPopBase();
        //this.mergeConfidenceGene(); //chose polbase or b73base
        //this.binData(); //deprecated
        //this.scatterData();
        //this.deleteriousAndGeneByTissue();
        //this.transitionAbility1(); //deprecated
        //this.transitionAbility2();
        //this.mkBarData();
        
        //this.pollenSpecific();
        //this.mkPollenTrainingSet();
        //this.selectOnlyPollenDeleterious();
    }
    
    void scatterData () {
        String tranFileS = "M:\\production\\maf\\003_expression\\delAndExpression\\transcriptome.txt";
        String proFileS = "M:\\production\\maf\\003_expression\\delAndExpression\\proteome.txt";
        String phoFileS = "M:\\production\\maf\\003_expression\\delAndExpression\\phosphoproteome.txt";
        String tranBinFileS = "M:\\production\\maf\\003_expression\\delAndExpression\\scatter\\transcriptome.bin.txt";
        String proBinFileS = "M:\\production\\maf\\003_expression\\delAndExpression\\scatter\\proteome.bin.txt";
        String phoBinFileS = "M:\\production\\maf\\003_expression\\delAndExpression\\scatter\\phosphoproteome.bin.txt";
        int binSize = 200;
        this.subBinData(tranFileS, tranBinFileS, binSize);
        this.subBinData(proFileS, proBinFileS, binSize);
        this.subBinData(phoFileS, phoBinFileS, binSize);
//        this.subBinDataSortByExpressionValue(tranFileS, tranBinFileS, binSize);
//        this.subBinDataSortByExpressionValue(proFileS, proBinFileS, binSize);
//        this.subBinDataSortByExpressionValue(phoFileS, phoBinFileS, binSize);
    }
    
     void subBinDataSortByExpressionValue (String inputFileS, String outputFileS, int binSize) {
        Table t = new Table (inputFileS);
        int[][] bound = FArrayUtils.getSubsetsIndicesBySubsetSize(t.getRowNumber(), binSize);
        double[][][] exAndError = new double[bound.length][4][2];
        double[][][] tsAndError = new double[bound.length][4][2];
        t.sortByValueOnColumnName("Mean");
        double[][] muValues = new double[4][];
        muValues[0] = t.getDoubleArrayByColumn(t.getColumnIndex("Syn"));
        muValues[1] = t.getDoubleArrayByColumn(t.getColumnIndex("Del"));
        muValues[2] = t.getDoubleArrayByColumn(t.getColumnIndex("RatioDelVsSyn"));
        muValues[3] = t.getDoubleArrayByColumn(t.getColumnIndex("Mean"));
        for (int i = 0; i <  bound.length; i++) {
            TDoubleArrayList[] muList = new TDoubleArrayList[4];
            for (int j = 0; j < muList.length; j++) muList[j] = new TDoubleArrayList();
            for (int j = bound[i][0]; j < bound[i][1]; j++) {
                for (int k = 0; k < muValues.length; k++) {
                    muList[k].add(muValues[k][j]);
                }
            }
            for (int j = 0; j < muList.length; j++) {
                double[] muV = muList[j].toArray();
                DescriptiveStatistics d = new DescriptiveStatistics(muV);
                double mean = d.getMean();
                double sd = d.getStandardDeviation();
                double se = sd/Math.sqrt(muV.length);
                exAndError[i][j][0] = mean;
                exAndError[i][j][1] = se;
            }
        }
        t.sortByValueOnColumnName("RSD");
        muValues = new double[4][];
        muValues[0] = t.getDoubleArrayByColumn(t.getColumnIndex("Syn"));
        muValues[1] = t.getDoubleArrayByColumn(t.getColumnIndex("Del"));
        muValues[2] = t.getDoubleArrayByColumn(t.getColumnIndex("RatioDelVsSyn"));
        muValues[3] = t.getDoubleArrayByColumn(t.getColumnIndex("RSD"));
        for (int i = 0; i <  bound.length; i++) {
            TDoubleArrayList[] muList = new TDoubleArrayList[4];
            for (int j = 0; j < muList.length; j++) muList[j] = new TDoubleArrayList();
            for (int j = bound[i][0]; j < bound[i][1]; j++) {
                for (int k = 0; k < muValues.length; k++) {
                    muList[k].add(muValues[k][j]);
                }
            }
            for (int j = 0; j < muList.length; j++) {
                double[] muV = muList[j].toArray();
                DescriptiveStatistics d = new DescriptiveStatistics(muV);
                double mean = d.getMean();
                double sd = d.getStandardDeviation();
                double se = sd/Math.sqrt(muV.length);
                tsAndError[i][j][0] = mean;
                tsAndError[i][j][1] = se;
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outputFileS);
            bw.write("Group rank\tExpressSynMean\tExpressSynError\tExpressDelMean\tExpressDelError\tExpressRatioDVsSMean\tExpressRatioDVsSError\tExpressionMean\tExpressionError\tTsSynMean\tTsSynError\tTsDelMean\tTsDelError\tTsRatioDVsSMean\tTsRatioDVsSError\tTsMean\tTsError");
            bw.newLine();
            for (int i = 0; i < bound.length; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(i+1);
                for (int j = 0; j < exAndError[i].length; j++) {
                    for (int k = 0; k < exAndError[i][j].length; k++) {
                        sb.append("\t").append(exAndError[i][j][k]);
                    }
                }
                for (int j = 0; j < tsAndError[i].length; j++) {
                    for (int k = 0; k < tsAndError[i][j].length; k++) {
                        sb.append("\t").append(tsAndError[i][j][k]);
                    }
                }
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
    
    void binData () {
        String tranFileS = "M:\\production\\maf\\003_expression\\delAndExpression\\transcriptome.txt";
        String proFileS = "M:\\production\\maf\\003_expression\\delAndExpression\\proteome.txt";
        String phoFileS = "M:\\production\\maf\\003_expression\\delAndExpression\\phosphoproteome.txt";
        String tranBinFileS = "M:\\production\\maf\\003_expression\\delAndExpression\\bin\\transcriptome.bin.txt";
        String proBinFileS = "M:\\production\\maf\\003_expression\\delAndExpression\\bin\\proteome.bin.txt";
        String phoBinFileS = "M:\\production\\maf\\003_expression\\delAndExpression\\bin\\phosphoproteome.bin.txt";
        int binSize = 1000;
        this.subBinData(tranFileS, tranBinFileS, binSize);
        this.subBinData(proFileS, proBinFileS, binSize);
        this.subBinData(phoFileS, phoBinFileS, binSize);
    }
    
    void subBinData (String inputFileS, String outputFileS, int binSize) {
        Table t = new Table (inputFileS);
        int[][] bound = FArrayUtils.getSubsetsIndicesBySubsetSize(t.getRowNumber(), binSize);
        double[][][] exAndError = new double[bound.length][3][2];
        double[][][] tsAndError = new double[bound.length][3][2];
        double[][][] vsAndError = new double[bound.length][3][2];
        t.sortByValueOnColumnName("Syn");
        System.out.println(t.content[t.getRowNumber()-1][0]);
        double[] ex = t.getDoubleArrayByColumn(1);
        double[] ts = t.getDoubleArrayByColumn(2);
        double[] vs = t.getDoubleArrayByColumn(t.getColumnIndex("Syn"));
        for (int i = 0; i <  bound.length; i++) {
            TDoubleArrayList exList = new TDoubleArrayList();
            TDoubleArrayList tsList = new TDoubleArrayList();
            TDoubleArrayList vsList = new TDoubleArrayList();
            for (int j = bound[i][0]; j < bound[i][1]; j++) {
                exList.add(ex[j]);
                tsList.add(ts[j]);
                vsList.add(vs[j]);
            }
            double[] exSub = exList.toArray();
            double[] tsSub = tsList.toArray();
            double[] vsSub = vsList.toArray();
            DescriptiveStatistics d = new DescriptiveStatistics(exSub);
            double mean = d.getMean();
            double sd = d.getStandardDeviation();
            double se = sd/Math.sqrt(exSub.length);
            exAndError[i][0][0] = mean;
            exAndError[i][0][1] = se;
            d = new DescriptiveStatistics(tsSub);
            mean = d.getMean();
            sd = d.getStandardDeviation();
            se = sd/Math.sqrt(exSub.length);
            tsAndError[i][0][0] = mean;
            tsAndError[i][0][1] = se;
            d = new DescriptiveStatistics(vsSub);
            mean = d.getMean();
            sd = d.getStandardDeviation();
            se = sd/Math.sqrt(exSub.length);
            vsAndError[i][0][0] = mean;
            vsAndError[i][0][1] = se;
        }
        t.sortByValueOnColumnName("Del");
        ex = t.getDoubleArrayByColumn(1);
        ts = t.getDoubleArrayByColumn(2);
        vs = t.getDoubleArrayByColumn(t.getColumnIndex("Del"));
        for (int i = 0; i <  bound.length; i++) {
            TDoubleArrayList exList = new TDoubleArrayList();
            TDoubleArrayList tsList = new TDoubleArrayList();
            TDoubleArrayList vsList = new TDoubleArrayList();
            for (int j = bound[i][0]; j < bound[i][1]; j++) {
                exList.add(ex[j]);
                tsList.add(ts[j]);
                vsList.add(vs[j]);
            }
            double[] exSub = exList.toArray();
            double[] tsSub = tsList.toArray();
            double[] vsSub = vsList.toArray();
            DescriptiveStatistics d = new DescriptiveStatistics(exSub);
            double mean = d.getMean();
            double sd = d.getStandardDeviation();
            double se = sd/Math.sqrt(exSub.length);
            exAndError[i][1][0] = mean;
            exAndError[i][1][1] = se;
            d = new DescriptiveStatistics(tsSub);
            mean = d.getMean();
            sd = d.getStandardDeviation();
            se = sd/Math.sqrt(exSub.length);
            tsAndError[i][1][0] = mean;
            tsAndError[i][1][1] = se;
            d = new DescriptiveStatistics(vsSub);
            mean = d.getMean();
            sd = d.getStandardDeviation();
            se = sd/Math.sqrt(exSub.length);
            vsAndError[i][1][0] = mean;
            vsAndError[i][1][1] = se;
        }
        t.sortByValueOnColumnName("RatioDelVsSyn");
        ex = t.getDoubleArrayByColumn(1);
        ts = t.getDoubleArrayByColumn(2);
        vs = t.getDoubleArrayByColumn(t.getColumnIndex("RatioDelVsSyn"));
        for (int i = 0; i <  bound.length; i++) {
            TDoubleArrayList exList = new TDoubleArrayList();
            TDoubleArrayList tsList = new TDoubleArrayList();
            TDoubleArrayList vsList = new TDoubleArrayList();
            for (int j = bound[i][0]; j < bound[i][1]; j++) {
                exList.add(ex[j]);
                tsList.add(ts[j]);
                vsList.add(vs[j]);
            }
            double[] exSub = exList.toArray();
            double[] tsSub = tsList.toArray();
            double[] vsSub = vsList.toArray();
            DescriptiveStatistics d = new DescriptiveStatistics(exSub);
            double mean = d.getMean();
            double sd = d.getStandardDeviation();
            double se = sd/Math.sqrt(exSub.length);
            exAndError[i][2][0] = mean;
            exAndError[i][2][1] = se;
            d = new DescriptiveStatistics(tsSub);
            mean = d.getMean();
            sd = d.getStandardDeviation();
            se = sd/Math.sqrt(exSub.length);
            tsAndError[i][2][0] = mean;
            tsAndError[i][2][1] = se;
            d = new DescriptiveStatistics(vsSub);
            mean = d.getMean();
            sd = d.getStandardDeviation();
            se = sd/Math.sqrt(exSub.length);
            vsAndError[i][2][0] = mean;
            vsAndError[i][2][1] = se;
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outputFileS);
            bw.write("Group rank\tExpressSynMean\tExpressSynError\tExpressDelMean\tExpressDelError\tExpressRatioDVsSMean\tExpressRatioDVsSError\tTsSynMean\tTsSynError\tTsDelMean\tTsDelError\tTsRatioDVsSMean\tTsRatioDVsSError\tSynMean\tSynError\tDelMean\tDelError\tRatioDVsSMean\tRatioDVsSError");
            bw.newLine();
            for (int i = 0; i < bound.length-1; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(i+1);
                for (int j = 0; j < exAndError[i].length; j++) {
                    for (int k = 0; k < exAndError[i][j].length; k++) {
                        sb.append("\t").append(exAndError[i][j][k]);
                    }
                }
                for (int j = 0; j < tsAndError[i].length; j++) {
                    for (int k = 0; k < tsAndError[i][j].length; k++) {
                        sb.append("\t").append(tsAndError[i][j][k]);
                    }
                }
                for (int j = 0; j < vsAndError[i].length; j++) {
                    for (int k = 0; k < vsAndError[i][j].length; k++) {
                        sb.append("\t").append(vsAndError[i][j][k]);
                    }
                }
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
    
    
    void mkBarData () {
        String infileS = "M:\\production\\maf\\003_expression\\transition\\delTransition2.txt";
        String outfileS = "M:\\production\\maf\\003_expression\\transition\\barError.txt";
        Table t = new Table (infileS);
        double step = 1;
        int minSize = 50;
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;
        for (int i = 0; i < t.getRowNumber(); i++) {
            double v = t.getDoubleValue(i, 6);
            if (Double.isNaN(v)) continue;
            if (v == 0) continue;
            if (v >max ) max = v;
            if (v < min) min = v;
        }
        System.out.println(min);
        System.out.println(max);
        double logMin = (int)(Math.log(min) / Math.log(2));
        double logMax = (int)(Math.log(max) / Math.log(2));
        int groupNum = (int)((logMax-logMin)/step);
        double[] bound = new double[groupNum+1];
        for (int i = 0; i < bound.length; i++) {
            bound[i] = logMin+i*step;
        }
        TDoubleArrayList[][] valueList = new TDoubleArrayList[groupNum][6];
        for (int i = 0; i < valueList.length; i++) {
            for (int j = 0; j < valueList[i].length; j++) {
                valueList[i][j] = new TDoubleArrayList();
            }
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            double v = t.getDoubleValue(i, 6);
            if (Double.isNaN(v)) continue;
            if (v == 0) continue;
            double logV = (int)(Math.log(v) / Math.log(2));
            if (logV < logMin) continue;
            if (logV >= logMax) continue;
            logV = (Math.log(v) / Math.log(2));
            int index = Arrays.binarySearch(bound, logV);
            if (index < 0) index = -index-1;
            if (index >= groupNum) index = groupNum-1;
            valueList[index][0].add(t.getDoubleValue(i, 7));
            valueList[index][1].add(t.getDoubleValue(i, 8));
            valueList[index][2].add(t.getDoubleValue(i, 9));
            valueList[index][3].add(t.getDoubleValue(i, 10));
            valueList[index][4].add(t.getDoubleValue(i, 2));
            valueList[index][5].add(t.getDoubleValue(i, 3));
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Range\tSize\tSynMean\tSynError\tDelMean\tDelError\tRatioMean\tRatioError\tGERPMean\tGERPError\tExpressionMean\tExpressionError\tTranslationMean\tTranslationError");
            bw.newLine();
            for (int i = 0; i < groupNum; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(bound[i]).append("\t").append(valueList[i][0].size());
                for (int j = 0; j < valueList[0].length; j++) {
                    double[] value = valueList[i][j].toArray();
                    DescriptiveStatistics d = new DescriptiveStatistics(value);
                    double mean = d.getMean();
                    double sd = d.getStandardDeviation();
                    double se = sd/Math.sqrt(value.length);
                    sb.append("\t").append(mean).append("\t").append(se);
                }
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
    
    void transitionAbility2 () {
        String expFileS = "M:\\production\\maf\\003_expression\\delAndExpression\\transcriptome.txt";
        String proFileS = "M:\\production\\maf\\003_expression\\delAndExpression\\proteome.txt";
        String outfileS = "M:\\production\\maf\\003_expression\\transition\\delTransition2.txt";
        Table tE = new Table (expFileS);
        Table tP = new Table (proFileS);
        String[] geneE = tE.getStringArrayByColumn(0);
        String[] geneP = tP.getStringArrayByColumn(0);
        String[] commonGene = FArrayUtils.getIntersection(geneE, geneP);
        String[] tissueE = new String[tE.getColumnNumber()-8];
        String[] tissueP = new String[tP.getColumnNumber()-8];
        for (int i = 0; i < tissueE.length; i++) {
            tissueE[i] = tE.header[i+4];
        }
        for (int i = 0; i < tissueP.length; i++) {
            tissueP[i] = tP.header[i+4];
        }
        String[] commonTissue = FArrayUtils.getIntersection(tissueE, tissueP);
        Arrays.sort(commonTissue);
        System.out.println(commonTissue.length);
        HashMap<String, String> geneDelMap = new HashMap();
        for (int i = 0; i < tE.getRowNumber(); i++) {
            StringBuilder sb = new StringBuilder();
            sb.append(tE.content[i][tE.getColumnNumber()-4]).append("\t").append(tE.content[i][tE.getColumnNumber()-3]).append("\t");
            sb.append(tE.content[i][tE.getColumnNumber()-2]).append("\t").append(tE.content[i][tE.getColumnNumber()-1]);
            geneDelMap.put(tE.content[i][0], sb.toString());
        }
        HashMap<String, Double>[] geneEMap = new HashMap[commonTissue.length];
        HashMap<String, Double>[] genePMap = new HashMap[commonTissue.length];
        for (int i = 0; i < commonTissue.length; i++) {
            geneEMap[i] = new HashMap();
            genePMap[i] = new HashMap();
        }
        double[] sumE = new double[commonTissue.length];
        int[] indices = new int[commonTissue.length];
        for (int i = 0; i < commonTissue.length; i++) {
            for (int j = 0; j < tE.header.length; j++) {
                if (commonTissue[i].equals(tE.header[j])) {
                    indices[i] = j;
                    break;
                }
            }
        }
        for (int i = 0; i < tE.getRowNumber(); i++) {
            if (Arrays.binarySearch(commonGene, tE.content[i][0]) < 0) continue;
            for (int j = 0; j < commonTissue.length; j++) {
                geneEMap[j].put(tE.content[i][0], Double.valueOf(tE.content[i][indices[j]]));
                sumE[j]+=Double.valueOf(tE.content[i][indices[j]]);
            }
        }
        double[] sumP = new double[commonTissue.length];
        indices = new int[commonTissue.length];
        for (int i = 0; i < commonTissue.length; i++) {
            for (int j = 0; j < tP.header.length; j++) {
                if (commonTissue[i].equals(tP.header[j])) {
                    indices[i] = j;
                    break;
                }
            }
        }
        for (int i = 0; i < tP.getRowNumber(); i++) {
            if (Arrays.binarySearch(commonGene, tP.content[i][0]) < 0) continue;
            for (int j = 0; j < commonTissue.length; j++) {
                genePMap[j].put(tP.content[i][0], Double.valueOf(tP.content[i][indices[j]]));
                sumP[j]+=Double.valueOf(tP.content[i][indices[j]]);
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Gene\tTissue\tExpresson\tTranslation\tRelativeExpression\tRelativeTranslation\tRatioTVsE\tSyn\tDel\tRatioDelVsSyn\tMeanGerp");
            bw.newLine();
            for (int i = 0; i < commonGene.length; i++) {
                for (int j = 0; j < commonTissue.length; j++) {
                    StringBuilder sb = new StringBuilder(commonGene[i]);
                    sb.append("\t").append(commonTissue[j]).append("\t").append(geneEMap[j].get(commonGene[i])).append("\t").append(genePMap[j].get(commonGene[i])).append("\t");
                    double abE = geneEMap[j].get(commonGene[i])/sumE[j];
                    double abP = genePMap[j].get(commonGene[i])/sumP[j];
                    double ratio = Double.NaN;
                    if (abE != 0) {
                        ratio = abP/abE; 
                    }
                    sb.append(abE).append("\t").append(abP).append("\t").append(ratio).append("\t").append(geneDelMap.get(commonGene[i]));
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
    
    void transitionAbility1 () {
        String expFileS = "M:\\production\\maf\\003_expression\\delAndExpression\\transcriptome.txt";
        String proFileS = "M:\\production\\maf\\003_expression\\delAndExpression\\proteome.txt";
        String outfileS = "M:\\production\\maf\\003_expression\\transition\\delTransition1.txt";
        Table tE = new Table (expFileS);
        Table tP = new Table (proFileS);
        String[] geneE = tE.getStringArrayByColumn(0);
        String[] geneP = tP.getStringArrayByColumn(0);
        String[] commonGene = FArrayUtils.getIntersection(geneE, geneP);
        HashMap<String, Double> geneExpMap = new HashMap();
        HashMap<String, Double> geneProMap = new HashMap();
        HashMap<String, String> geneDelMap = new HashMap();
        double sumE = 0;
        for (int i = 0; i < tE.getRowNumber(); i++) {
            if (Arrays.binarySearch(commonGene, tE.content[i][0]) < 0) continue;
            geneExpMap.put(tE.content[i][0], tE.getDoubleValue(i, 1));
            sumE+=tE.getDoubleValue(i, 1);
            StringBuilder sb = new StringBuilder();
            sb.append(tE.content[i][tE.getColumnNumber()-4]).append("\t").append(tE.content[i][tE.getColumnNumber()-3]).append("\t");
            sb.append(tE.content[i][tE.getColumnNumber()-2]).append("\t").append(tE.content[i][tE.getColumnNumber()-1]);
            geneDelMap.put(tE.content[i][0], sb.toString());
        }
        double sumP = 0;
        for (int i = 0; i < tP.getRowNumber(); i++) {
            if (Arrays.binarySearch(commonGene, tP.content[i][0]) < 0) continue;
            geneProMap.put(tP.content[i][0], tP.getDoubleValue(i, 1));
            sumP+=tP.getDoubleValue(i, 1);
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Gene\tExpresson\tTranslation\tRatioTVsE\tSyn\tDel\tRatioDelVsSyn\tMeanGerp");
            bw.newLine();
            for (int i = 0; i < commonGene.length; i++) {
                StringBuilder sb = new StringBuilder(commonGene[i]);
                double abE = geneExpMap.get(commonGene[i])/sumE;
                double abP = geneProMap.get(commonGene[i])/sumP;
                sb.append("\t").append(abE).append("\t").append(abP).append("\t").append(abP/abE).append("\t").append(geneDelMap.get(commonGene[i]));
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
    
    
    void selectOnlyPollenDeleterious () {
        String infileS = "M:\\production\\maf\\expression\\pollenVariantCalling\\trainingSet\\pollenTraining.txt";
        String deleteriousFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\class\\Non_Synonymous_Deleterious.txt";
        String outfileS = "M:\\production\\maf\\expression\\pollenVariantCalling\\trainingSet\\pollenTrainingDel.txt";
        int chrNum = 10;
        Table t = new Table (deleteriousFileS);
        TIntArrayList[] posList = new TIntArrayList[chrNum];
        for (int i = 0; i < posList.length; i++) {
            posList[i] = new TIntArrayList();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            int pos = t.getIntValue(i, 1);
            int chrIndex = t.getIntValue(i, 0) - 1;
            posList[chrIndex].add(pos);
        }
        int[][] poss = new int[chrNum][];
        for (int i = 0; i < posList.length; i++) {
            poss[i] = posList[i].toArray();
            Arrays.sort(poss[i]);
        }
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write(br.readLine());
            bw.newLine();
            String temp = null;
            while ((temp = br.readLine()) != null) {
                List<String> l = FStringUtils.fastSplit(temp);
                if (l.get(2).equals("0")) {
                    int chrIndex = Integer.valueOf(l.get(0))-1;
                    int pos = Integer.valueOf(l.get(1));
                    if (Arrays.binarySearch(poss[chrIndex], pos) < 0) continue;
                }
                bw.write(temp);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    void mkPollenTrainingSet () {
        String pollenSpecificFileS = "M:\\production\\maf\\expression\\pollenVariantCalling\\geneClass\\pollenSpecificGene.txt";
        String pollenFileS = "M:\\production\\maf\\expression\\pollenVariantCalling\\geneClass\\pollenGene.txt";
        String hmpDirS = "M:\\production\\maf\\annotations\\siftScore\\001_hmp321Info\\gerpAncestral\\";
        String geneFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
        String outFileS = "M:\\production\\maf\\expression\\pollenVariantCalling\\trainingSet\\pollenTraining.txt";
        String outFileS2 = "M:\\production\\maf\\expression\\pollenVariantCalling\\trainingSet\\pollenSpeTraining.txt";
        Table t = new Table (pollenSpecificFileS);
        String[] pollenSpeGene = t.getStringArrayByColumn(0);
        Arrays.sort(pollenSpeGene);
        t = new Table (pollenFileS);
        String[] pollenGene = t.getStringArrayByColumn(0);
        Arrays.sort(pollenGene);
        ArrayList<Range> pollenCDSList = new ArrayList();
        ArrayList<Range> pollenSpeCDSList = new ArrayList();
        ArrayList<Range> nonPollenCDSList = new ArrayList();
        GeneFeature gf = new GeneFeature(geneFileS);
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            String query = gf.getGeneName(i);
            if (Arrays.binarySearch(pollenGene, query) < 0) {
                nonPollenCDSList.addAll(gf.getCDSList(i, 0));
            }
            else if (Arrays.binarySearch(pollenGene, query) >= 0) {
                pollenCDSList.addAll(gf.getCDSList(i, 0));
                if (Arrays.binarySearch(pollenSpeGene, query) >= 0) {
                    pollenSpeCDSList.addAll(gf.getCDSList(i, 0));
                }
            }
        }
        Ranges nonPollenR = new Ranges(nonPollenCDSList, "nonPollen");
        Ranges pollenR = new Ranges(pollenCDSList, "pollen");
        Ranges pollenSpeR = new Ranges(pollenSpeCDSList, "pollenSpe");
        File[] fs = new File (hmpDirS).listFiles();
        ArrayList<String>[] nonPollen = new ArrayList[fs.length];
        ArrayList<String>[] pollen = new ArrayList[fs.length];
        ArrayList<String>[] pollenSpe = new ArrayList[fs.length];
        for (int i = 0; i < fs.length; i++) {
            nonPollen[i] = new ArrayList();
            pollen[i] = new ArrayList();
            pollenSpe[i] = new ArrayList();
        }
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            int chr = Integer.valueOf(f.getName().replaceFirst("hmp32Info_chr", "").replaceFirst(".txt.gz", ""));
            int chrIndex = chr-1;
            try {
                BufferedReader br = IoUtils.getTextGzipReader(f.getAbsolutePath());
                String temp = br.readLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    List<String> l = FStringUtils.fastSplit(temp);
                    int pos = Integer.valueOf(l.get(1));
                    if (nonPollenR.isInRanges(chr, pos)) {
                        StringBuilder sb = new StringBuilder();
                        String[] tem = l.get(8).split(",");
                        double maf = Double.valueOf(tem[tem.length-1]);
                        
                        sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append("1").append("\t").append(tem.length).append("\t");
                        sb.append(maf).append("\t").append(l.get(9)).append("\t").append(l.get(10));
                        nonPollen[chrIndex].add(sb.toString());
                    }
                    else if (pollenR.isInRanges(chr, pos)){
                        StringBuilder sb = new StringBuilder();
                        String[] tem = l.get(8).split(",");
                        double maf = Double.valueOf(tem[tem.length-1]);
                        sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append("0").append("\t").append(tem.length).append("\t");
                        sb.append(maf).append("\t").append(l.get(9)).append("\t").append(l.get(10));
                        pollen[chrIndex].add(sb.toString());
                        if (pollenSpeR.isInRanges(chr, pos)){
                            pollenSpe[chrIndex].add(sb.toString());
                        }
                    }
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt)+"\t"+f.getName());
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outFileS);
            bw.write("Chr\tPos\tIfPollen\tNumAlleles\tMinorAlleleFrequency\tSiteDepth\tHetCount");
            bw.newLine();
            for (int i = 0; i < pollen.length; i++) {
                for (int j = 0; j < pollen[i].size(); j++) {
                    bw.write(pollen[i].get(j));
                    bw.newLine();
                }
            }
            for (int i = 0; i < nonPollen.length; i++) {
                for (int j = 0; j < nonPollen[i].size(); j++) {
                    bw.write(nonPollen[i].get(j));
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
            bw = IoUtils.getTextWriter(outFileS2);
            bw.write("Chr\tPos\tIfPollen\tNumAlleles\tMinorAlleleFrequency\tSiteDepth\tHetCount");
            bw.newLine();
            for (int i = 0; i < pollenSpe.length; i++) {
                for (int j = 0; j < pollenSpe[i].size(); j++) {
                    bw.write(pollenSpe[i].get(j));
                    bw.newLine();
                }
            }
            for (int i = 0; i < nonPollen.length; i++) {
                for (int j = 0; j < nonPollen[i].size(); j++) {
                    bw.write(nonPollen[i].get(j));
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
    
    
    
    void pollenSpecific () {
        String proteomeFileS = "M:\\production\\maf\\expression\\source\\flat\\proteome.txt";
        String outfileDirS = "M:\\production\\maf\\expression\\pollenVariantCalling\\geneClass\\";
        Table t = new Table (proteomeFileS);
        HashSet<String> pollenSet = new HashSet();
        HashSet<String> nonPollenSet = new HashSet();
        TIntArrayList indexList = new TIntArrayList();
        for (int i = 0; i < t.header.length; i++) {
            if (t.header[i].contains("Pollen")) indexList.add(i);
        }
        int[] pollenIndex = indexList.toArray();
        Arrays.sort(pollenIndex);
        for (int i = 0; i < t.getRowNumber(); i++) {
            for (int j = 1; j < t.getColumnNumber(); j++) {
                if (t.getDoubleValue(i, j) == 0) continue;
                int index = Arrays.binarySearch(pollenIndex, j);
                if (index < 0) {
                    nonPollenSet.add(t.content[i][0]);
                }
                else {
                    pollenSet.add(t.content[i][0]);
                }
            }
        }
        String[] pollenGene = pollenSet.toArray(new String[pollenSet.size()]);
        String[] nonPollenGene = nonPollenSet.toArray(new String[nonPollenSet.size()]);
        Arrays.sort(pollenGene);
        Arrays.sort(nonPollenGene);
        ArrayList<String> pollenSpecific = new ArrayList();
        for (int i = 0; i < pollenGene.length; i++) {
            if (Arrays.binarySearch(nonPollenGene, pollenGene[i]) < 0) pollenSpecific.add(pollenGene[i]);
        }
        String[] polleneSpecificGene = pollenSpecific.toArray(new String[pollenSpecific.size()]);
        Arrays.sort(polleneSpecificGene);
        String pollenFileS = new File (outfileDirS, "pollenGene.txt").getAbsolutePath();
        String nonPollenFileS = new File (outfileDirS, "nonPollenGene.txt").getAbsolutePath();
        String pollenSpecificFileS = new File (outfileDirS, "pollenSpecificGene.txt").getAbsolutePath();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(pollenFileS);
            bw.write("Gene");
            bw.newLine();
            for (int i = 0; i < pollenGene.length; i++) {
                bw.write(pollenGene[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
            bw = IoUtils.getTextWriter(nonPollenFileS);
            bw.write("Gene");
            bw.newLine();
            for (int i = 0; i < nonPollenGene.length; i++) {
                bw.write(nonPollenGene[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
            bw = IoUtils.getTextWriter(pollenSpecificFileS);
            bw.write("Gene");
            bw.newLine();
            for (int i = 0; i < polleneSpecificGene.length; i++) {
                bw.write(polleneSpecificGene[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    void deleteriousAndGeneByTissue () {
        String sourceDirS = "M:\\production\\maf\\003_expression\\delAndExpression\\";
        String outputDirS = "M:\\production\\maf\\003_expression\\delAndGeneByTissue\\";
        File[] fs = new File (sourceDirS).listFiles();
        fs = IoUtils.listFilesEndsWith(fs, ".txt");
        for (int i = 0; i < fs.length; i++) {
            HashMap<String, Double> synMap = new HashMap();
            HashMap<String, Double> delMap = new HashMap();
            HashMap<String, Double> ratioMap = new HashMap();
            Table t = new Table (fs[i].getAbsolutePath());
            for (int j = 0; j < t.getRowNumber(); j++) {
                synMap.put(t.content[j][0], Double.valueOf(t.content[j][t.getColumnNumber()-4]));
                delMap.put(t.content[j][0], Double.valueOf(t.content[j][t.getColumnNumber()-3]));
                ratioMap.put(t.content[j][0], Double.valueOf(t.content[j][t.getColumnNumber()-2]));
            }
            String[] tissue = new String[t.getColumnNumber()-8];
            TDoubleArrayList[] synList = new TDoubleArrayList[tissue.length];
            TDoubleArrayList[] delList = new TDoubleArrayList[tissue.length];
            TDoubleArrayList[] ratioList = new TDoubleArrayList[tissue.length];
            ArrayList<String>[] geneList = new ArrayList[tissue.length];
            double[][] syn  = new double[tissue.length][];
            double[][] del = new double[tissue.length][];
            double[][] ratio = new double[tissue.length][];
            for (int j = 0; j < tissue.length; j++) {
                synList[j] = new TDoubleArrayList();
                delList[j] = new TDoubleArrayList();
                ratioList[j] = new TDoubleArrayList();
                tissue[j] = t.header[j+4];
                geneList[j] = new ArrayList();
            }
            for (int j = 0; j < t.getRowNumber(); j++) {
                for (int k = 0; k < tissue.length; k++) {
                    if (Double.valueOf(t.content[j][k+4])==0) continue;
                    synList[k].add(synMap.get(t.content[j][0]));
                    delList[k].add(delMap.get(t.content[j][0]));
                    ratioList[k].add(ratioMap.get(t.content[j][0]));
                    geneList[k].add(t.content[j][0]);
                }
            }
            double[] meanSyn = new double[tissue.length];
            double[] meanDel = new double[tissue.length];
            double[] meanRatio = new double[tissue.length];
            double[] seSyn = new double[tissue.length];
            double[] seDel = new double[tissue.length];
            double[] seRatio = new double[tissue.length];
            for (int j = 0; j < tissue.length; j++) {
                syn[j] = synList[j].toArray();
                del[j] = delList[j].toArray();
                ratio[j] = ratioList[j].toArray();
                DescriptiveStatistics d = new DescriptiveStatistics(syn[j]);
                meanSyn[j] = d.getMean();
                seSyn[j] = d.getStandardDeviation()/Math.sqrt(syn[j].length);
                d = new DescriptiveStatistics(del[j]);
                meanDel[j] = d.getMean();
                seDel[j] = d.getStandardDeviation()/Math.sqrt(del[j].length);
                d = new DescriptiveStatistics(ratio[j]);
                meanRatio[j] = d.getMean();
                seRatio[j] = d.getStandardDeviation()/Math.sqrt(ratio[j].length);
            }
            int[] indices = FArrayUtils.getIndexByAscendingValue(meanRatio);
            for (int j = 0; j < indices.length; j++) {
                System.out.println(String.valueOf(meanSyn[indices[j]])+"\t"+String.valueOf(tissue[indices[j]]));
            }
            String outfileS = new File (outputDirS, fs[i].getName()).getAbsolutePath();
            try {
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                bw.write("Tissue\tMeanSyn\tSESyn\tMeanDel\tSEDel\tMeanRatio\tSERatio");
                bw.newLine();
                for (int j = 0; j < tissue.length; j++) {
                    StringBuilder sb = new StringBuilder(tissue[indices[j]]);
                    sb.append("\t").append(meanSyn[indices[j]]).append("\t").append(seSyn[indices[j]]);
                    sb.append("\t").append(meanDel[indices[j]]).append("\t").append(seDel[indices[j]]);
                    sb.append("\t").append(meanRatio[indices[j]]).append("\t").append(seRatio[indices[j]]);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
//            try {
//                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
//                bw.write("Gene\tSyn\tDel\tRatio\tTissueID\tTissue");
//                bw.newLine();
//                for (int j = 0; j < tissue.length; j++) {
//                    for (int k = 0; k < geneList[indices[j]].size(); k++) {
//                        StringBuilder sb = new StringBuilder(geneList[indices[j]].get(k));
//                        sb.append("\t").append(syn[indices[j]][k]).append("\t").append(del[indices[j]][k]).append("\t").append(ratio[indices[j]][k]);
//                        sb.append("\t").append(j+1).append("\t").append(tissue[indices[j]]);
//                        bw.write(sb.toString());
//                        bw.newLine();
//                    }
//                }
//                bw.flush();
//                bw.close();
//            }
//            catch (Exception e) {
//                e.printStackTrace();
//            }
            
            
        }
    }
    
    void mergeConfidenceGene () {
        String inputDirS = "M:\\production\\maf\\003_expression\\processedExpression\\";
        String highConfidenceFileS = "M:\\production\\maf\\001_variantDiscovery\\005_snpClass\\highConfidence_transcript.txt";
        String outputDirS = "M:\\production\\maf\\003_expression\\delAndExpression\\";
        new File (outputDirS).mkdir();
        Table t = new Table (highConfidenceFileS);
        HashMap<String, Double> geneSynMap = new HashMap();
        HashMap<String, Double> geneDelMap = new HashMap();
        HashMap<String, Double> geneRioMap = new HashMap();
        HashMap<String, Double> geneGerpMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.content[i][4].equals("0")) continue;
            String name = t.content[i][0];
            if (name.startsWith("GRM")) {
                name = name.split("_")[0];
            }
            double del = Double.valueOf(t.content[i][11]);
            double v = Double.NaN;
            if (Double.valueOf(t.content[i][6]) != 0) {
                v = del/Double.valueOf(t.content[i][6]);
            }
            geneSynMap.put(name, Double.valueOf(t.content[i][6]));
            geneDelMap.put(name, del);
            geneRioMap.put(name, v);
            geneGerpMap.put(name, Double.valueOf(t.content[i][16]));
        }
//        for (int i = 0; i < t.getRowNumber(); i++) {
//            if (t.content[i][4].equals("0")) continue;
//            String name = t.content[i][0];
//            if (name.startsWith("GRM")) {
//                name = name.split("_")[0];
//            }
//            double del = Double.valueOf(t.content[i][27]);
//            double v = Double.NaN;
//            if (Double.valueOf(t.content[i][23]) != 0) {
//                v = del/Double.valueOf(t.content[i][23]);
//            }
//            geneSynMap.put(name, Double.valueOf(t.content[i][25]));
//            geneDelMap.put(name, del);
//            geneRioMap.put(name, v);
//            geneGerpMap.put(name, Double.valueOf(t.content[i][16]));
//        }
        Set<String> geneSet = geneSynMap.keySet();
        File[] fs = new File (inputDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            try {
                BufferedReader br =IoUtils.getTextReader(fs[i].getAbsolutePath());
                BufferedWriter bw = IoUtils.getTextWriter(new File (outputDirS, fs[i].getName()).getAbsolutePath());
                bw.write(br.readLine()+"\tSyn\tDel\tRatioDelVsSyn\tMeanGerp");
                bw.newLine();
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    List<String> l = FStringUtils.fastSplit(temp);
                    if (!geneSet.contains(l.get(0))) continue;
                    if (Double.isNaN(geneRioMap.get(l.get(0)))) continue;
                    if (Double.isNaN(Double.valueOf(l.get(2)))) continue;
                    StringBuilder sb = new StringBuilder(temp);
                    sb.append("\t").append(geneSynMap.get(l.get(0))).append("\t").append(geneDelMap.get(l.get(0))).append("\t").append(geneRioMap.get(l.get(0)));
                    sb.append("\t").append(geneGerpMap.get(l.get(0)));
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
    }
    
    void processRawB73Base () {
        String inputDirS = "M:\\production\\maf\\003_expression\\source\\flat\\";
        String outputDirS = "M:\\production\\maf\\003_expression\\processedExpression\\";
        File[] fs = new File(inputDirS).listFiles();
        for (int i = 0; i <  fs.length; i++) {
            TIntArrayList skipList = new TIntArrayList();
            try {
                Table t = new Table (fs[i].getAbsolutePath());
                for (int j = 0; j < t.header.length; j++) {
                    if (t.header[j].toLowerCase().contains("w23")) skipList.add(j);
                    if (t.header[j].toLowerCase().contains("mo17")) skipList.add(j);
                }
                int[] skipIndex = skipList.toArray();
                Arrays.sort(skipIndex);
                String outfileS = new File (outputDirS, fs[i].getName()).getAbsolutePath();
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                StringBuilder sb = new StringBuilder("Gene\tMean\tRSD\tZeroCount");
                for (int k = 1; k < t.header.length; k++) {
                    if (Arrays.binarySearch(skipIndex, k) >= 0) continue;
                    sb.append("\t").append(t.header[k]);
                }
                bw.write(sb.toString());
                bw.newLine();
                HashSet<String> geneSet = new HashSet();
                
                for (int j = 0; j < t.getRowNumber(); j++) {
                    String name = t.content[j][0];
                    if (name.startsWith("GRM")) {
                        name = name.split("_")[0];
                    }
                    geneSet.add(name);
                }
                String[] genes = geneSet.toArray(new String[geneSet.size()]);
                Arrays.sort(genes);
                boolean[] isthere = new boolean[genes.length];
                int[] indices = new int[genes.length];
                for (int j = 0; j < t.getRowNumber(); j++) {
                    String name = t.content[j][0];
                    if (name.startsWith("GRM")) {
                        name = name.split("_")[0];
                    }
                    int index = Arrays.binarySearch(genes, name);
                    if (isthere[index] == true) continue;
                    isthere[index] = true;
                    indices[index] = j;
                }
                for (int j = 0; j < skipIndex.length; j++) {
                    skipIndex[j]--;
                }
                for (int j = 0; j < genes.length; j++) {
                    double[] value = new double[t.getColumnNumber()-1];
                    for (int k = 0; k < value.length; k++) {
                        value[k] = Double.valueOf(t.content[indices[j]][k+1]);
                    }
                    bw.write(this.getValueString(genes[j], value, skipIndex));
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
    
    String getValueString (String gene, double[] value, int[] skipIndex) {
        StringBuilder sb = new StringBuilder(gene);
        DescriptiveStatistics ds = new DescriptiveStatistics(value);
        double mean = ds.getMean();
        double sd = ds.getStandardDeviation();
        int cnt = 0;
        for (int i = 0; i < value.length; i++) {
            if (Arrays.binarySearch(skipIndex, i) >= 0) continue;
            if (value[i] == 0) {
                cnt++;
            }
        }
        sb.append("\t").append(mean).append("\t").append(sd/mean).append("\t").append(cnt);
        for (int i = 0; i < value.length; i++) {
            if (Arrays.binarySearch(skipIndex, i) >= 0) continue;
            sb.append("\t").append(value[i]);
        }
        return sb.toString();
    }
    
    void processRawPopBase () {
        String inputDirS = "M:\\production\\maf\\003_expression\\source\\flat\\";
        String outputDirS = "M:\\production\\maf\\003_expression\\processedExpression\\";
        File[] fs = new File(inputDirS).listFiles();
        for (int i = 0; i <  fs.length; i++) {
            try {
                Table t = new Table (fs[i].getAbsolutePath());
                String outfileS = new File (outputDirS, fs[i].getName()).getAbsolutePath();
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                StringBuilder sb = new StringBuilder("Gene\tMean\tRSD\tZeroCount");
                for (int k = 1; k < t.header.length; k++) {
                    sb.append("\t").append(t.header[k]);
                }
                bw.write(sb.toString());
                bw.newLine();
                HashSet<String> geneSet = new HashSet();
                
                for (int j = 0; j < t.getRowNumber(); j++) {
                    String name = t.content[j][0];
                    if (name.startsWith("GRM")) {
                        name = name.split("_")[0];
                    }
                    geneSet.add(name);
                }
                String[] genes = geneSet.toArray(new String[geneSet.size()]);
                Arrays.sort(genes);
                boolean[] isthere = new boolean[genes.length];
                int[] indices = new int[genes.length];
                for (int j = 0; j < t.getRowNumber(); j++) {
                    String name = t.content[j][0];
                    if (name.startsWith("GRM")) {
                        name = name.split("_")[0];
                    }
                    int index = Arrays.binarySearch(genes, name);
                    if (isthere[index] == true) continue;
                    isthere[index] = true;
                    indices[index] = j;
                }
                for (int j = 0; j < genes.length; j++) {
                    double[] value = new double[t.getColumnNumber()-1];
                    for (int k = 0; k < value.length; k++) {
                        value[k] = Double.valueOf(t.content[indices[j]][k+1]);
                    }
                    bw.write(this.getValueString(genes[j], value));
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
    
    String getValueString (String gene, double[] value) {
        StringBuilder sb = new StringBuilder(gene);
        DescriptiveStatistics ds = new DescriptiveStatistics(value);
        double mean = ds.getMean();
        double sd = ds.getStandardDeviation();
        int cnt = 0;
        for (int i = 0; i < value.length; i++) {
            if (value[i] == 0) {
                cnt++;
            }
        }
        sb.append("\t").append(mean).append("\t").append(sd/mean).append("\t").append(cnt);
        for (int i = 0; i < value.length; i++) {
            sb.append("\t").append(value[i]);
        }
        return sb.toString();
    }
}


