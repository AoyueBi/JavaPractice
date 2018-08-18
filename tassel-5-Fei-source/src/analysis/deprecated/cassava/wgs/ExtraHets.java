/*
* To change this license header, choose License Headers in Project Properties.
* To change this template file, choose Tools | Templates
* and open the template in the editor.
*/
package analysis.deprecated.cassava.wgs;

import format.Table;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import graphcis.r.DensityPlot;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.apache.commons.math3.stat.inference.TTest;
import utils.FArrayUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class ExtraHets {
    
    public ExtraHets () {
        this.extractUnimputatedVCF();
        //this.extractCultivatedClone();
        //this.summarizeDeleSites();
        //this.matchSNPsByDAF();
        //this.matchSNPsByDAF2();
        //this.matchSNPsByDAF3();
        //this.simulateHWE();
        //this.test();
    }
    
    void simulateHWE () {
        String infoFileS = "E:\\Research\\cassava\\revision\\extraHets\\source\\DeleteriousAllelesInformation.txt";
        String infileS = "E:\\Research\\cassava\\revision\\extraHets\\dele_impute_cultivated_LAC.txt";
        String outfileS = "E:\\Research\\cassava\\revision\\extraHets\\simulationHWE_LAC.txt";
        int nRep = 10000;
        int sampleSize = 100;
        HashMap<String, String> posDeleMap = new HashMap();
        Table t = new Table (infoFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            posDeleMap.put(t.content[i][0]+"_"+t.content[i][1], t.content[i][7]);
        }
        t = new Table (infileS);
        boolean[][][] geno = new boolean[t.getColumnNumber()-9][t.getRowNumber()][2];
        for (int i = 0; i < t.getRowNumber(); i++) {
            String derived = posDeleMap.get(t.content[i][0]+"_"+t.content[i][1]);
            int dInt = 1;
            if (t.content[i][3].equals(derived)) dInt = 0;
            for (int j = 0; j < geno.length; j++) {
                String[] tem = t.content[i][j+9].split(":")[0].split("\\|");
                for (int k = 0; k < tem.length; k++) {
                    if (Integer.valueOf(tem[k]) == dInt) {
                        geno[j][i][k] = true;
                    }
                }
                
            }
        }
        double[] oHets = new double[nRep];
        double[] oHomos = new double[nRep];
        double[] eHets = new double[nRep];
        double[] eHomos = new double[nRep];
        for (int i = 0; i < nRep; i++) {
            //int[] indices = FArrayUtils.getNonredundantRandomIntArray(geno.length, sampleSize);
            sampleSize = geno.length;
            int[] indices = FArrayUtils.getRandomIntArray(geno.length, sampleSize);
            Arrays.sort(indices);
            double[] daf = new double[geno[0].length];
            double[] oHet = new double[sampleSize];
            double[] oHomo = new double[sampleSize];
            for (int j = 0; j < indices.length; j++) {
                for (int k = 0; k < geno[0].length; k++) {
                    if (geno[indices[j]][k][0] == true) {
                        if (geno[indices[j]][k][1] == true) {
                            oHomo[j]++;
                        }
                        else {
                            oHet[j]++;
                        }
                    }
                    else {
                        if (geno[indices[j]][k][1] == true) {
                            oHet[j]++;
                        }
                    }
                }
            }
            for (int j = 0; j < geno[0].length; j++) {
                for (int k = 0; k < indices.length; k++) {
                    for (int l = 0; l < 2; l++) {
                        if (geno[indices[k]][j][l] == true) daf[j]++;
                    }
                }
                daf[j] = daf[j]/sampleSize/2;
            }
            DescriptiveStatistics d = new DescriptiveStatistics(oHet);
            double meanOHet = d.getMean();
            d = new DescriptiveStatistics(oHomo);
            double meanOHomo = d.getMean();
            double eHet = 0;
            double eHomo = 0;
            for (int j = 0; j < geno[0].length; j++) {
                eHet += 2*(daf[j])*(1-daf[j]);
                eHomo += daf[j]*daf[j];
            }
            oHets[i] = meanOHet;
            oHomos[i] = meanOHomo;
            eHets[i] = eHet;
            eHomos[i] = eHomo;
            if (i%100 == 0) System.out.println(String.valueOf(i));
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("OHomo\tEHomo\tOHet\tEHet");
            bw.newLine();
            for (int i = 0; i < nRep; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(oHomos[i]).append("\t").append(eHomos[i]).append("\t");
                sb.append(oHets[i]).append("\t").append(eHets[i]);
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
    
    void matchSNPsByDAF3 () {
        String infileS = "E:\\Research\\cassava\\revision\\extraHets\\dele_raw_cultivated_summary.txt";
        String outfileS = "E:\\Research\\cassava\\revision\\extraHets\\extraHet_ByDAF.txt";
        double maxDaf = 0.5;
        double frequencyInterval = 0.100;
        Table t = new Table (infileS);
        DeleSNP[] snps = new DeleSNP[t.getRowNumber()];
        for (int i = 0; i <  t.getRowNumber(); i++) {
            snps[i] = new DeleSNP(t.content[i]);
        }
        TDoubleArrayList boundList = new TDoubleArrayList();
        double current = 0;
        while (current < maxDaf) {
            boundList.add(current);
            current+=frequencyInterval;
        }
        boundList.add(current);
        double[] bound = boundList.toArray();
        ArrayList<DeleSNP>[] snpGroupList = new ArrayList[bound.length-1];
        DeleSNP[][] snpGroup = new DeleSNP[bound.length-1][];
        for (int i = 0; i < snpGroupList.length; i++) {
            snpGroupList[i] = new ArrayList();
        }
        for (int i = 0; i < snps.length; i++) {
            int index = Arrays.binarySearch(bound, snps[i].daf);
            if (index < 0) index = -index - 2;
            if (index >= snpGroup.length) continue;
            snpGroupList[index].add(snps[i]);
        }
        for (int i = 0; i < snpGroup.length; i++) {
            snpGroup[i] = snpGroupList[i].toArray(new DeleSNP[snpGroupList[i].size()]);
            Arrays.sort(snpGroup[i]);
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder();
            sb.append("DAFGroup\tSNPCount");
            String[] header = {"SampleSize", "RatioMean", "RatioError","SiftMean","SiftError","GerpMean", "GerpError", "P_ttest_SIFT", "P_ttest_GERP"};
            for (int i = 0; i < header.length; i++) {
                sb.append("\t").append("High"+header[i]);
            }
            for (int i = 0; i < header.length; i++) {
                sb.append("\t").append("Low"+header[i]);
            }
            for (int i = 0; i < header.length; i++) {
                sb.append("\t").append("HWE"+header[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < snpGroup.length; i++) {
                System.out.println(i);
                sb = new StringBuilder();
                sb.append(bound[i]).append("\t").append(snpGroup[i].length);

                double[] values = this.getGroupEstimates(snpGroup[i]);
                for (int j = 0; j < values.length; j++) {
                    sb.append("\t").append(values[j]);
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
    
    public double[] getGroupEstimates (DeleSNP[] snpGroup) {
        double[] values = new double[27];
        TIntArrayList lowList = new TIntArrayList();
        TIntArrayList highList = new TIntArrayList();
        TIntArrayList hweList = new TIntArrayList();
        for (int i = 0; i < snpGroup.length; i++) {
            if (snpGroup[i].p < 0.05) {
                if (snpGroup[i].hetRatio > 1) highList.add(i);
                else lowList.add(i);
            }
            else if (snpGroup[i].p > 0.99) hweList.add(i);
        }
        //{"SampleSize", "RatioMean", "RatioError","SiftMean","SiftError","GerpMean", "GerpError", "P_ttest"};
        int[] index = highList.toArray();
        double[] ratio = new double[index.length];
        double[] sift = new double[index.length];
        double[] gerp = new double[index.length];
        for (int i = 0; i < index.length; i++) {
            ratio[i] = snpGroup[index[i]].hetRatio;
            sift[i] = snpGroup[index[i]].sift;
            gerp[i] = snpGroup[index[i]].gerp;
        }
        values[0] = index.length;
        DescriptiveStatistics ds = new DescriptiveStatistics(ratio);
        double mean = ds.getMean();
        double sd = ds.getStandardDeviation();
        double error = sd/Math.sqrt(snpGroup.length);
        values[1] = mean; values[2] = error;
        ds = new DescriptiveStatistics(sift);
        mean = ds.getMean();
        sd = ds.getStandardDeviation();
        error = sd/Math.sqrt(snpGroup.length);
        values[3] = mean; values[4] = error;
        ds = new DescriptiveStatistics(gerp);
        mean = ds.getMean();
        sd = ds.getStandardDeviation();
        error = sd/Math.sqrt(snpGroup.length);
        values[5] = mean; values[6] = error;
        values[7] = Double.NaN; values[8] = Double.NaN;
        double[] highSift = sift;
        double[] highGerp = gerp;
        
        index = lowList.toArray();
        ratio = new double[index.length];
        sift = new double[index.length];
        gerp = new double[index.length];
        for (int i = 0; i < index.length; i++) {
            ratio[i] = snpGroup[index[i]].hetRatio;
            sift[i] = snpGroup[index[i]].sift;
            gerp[i] = snpGroup[index[i]].gerp;
        }
        values[9] = index.length;
        ds = new DescriptiveStatistics(ratio);
        mean = ds.getMean();
        sd = ds.getStandardDeviation();
        error = sd/Math.sqrt(snpGroup.length);
        values[10] = mean; values[11] = error;
        ds = new DescriptiveStatistics(sift);
        mean = ds.getMean();
        sd = ds.getStandardDeviation();
        error = sd/Math.sqrt(snpGroup.length);
        values[12] = mean; values[13] = error;
        ds = new DescriptiveStatistics(gerp);
        mean = ds.getMean();
        sd = ds.getStandardDeviation();
        error = sd/Math.sqrt(snpGroup.length);
        values[14] = mean; values[15] = error;
        double[] lowSift = sift;
        double[] lowGerp = gerp;
        TTest st = new TTest();
        double p = Double.NaN;
        if (lowSift.length > 1 && highSift.length > 1) {
            p = st.homoscedasticTTest(lowSift, highSift);
        }
        values[16] = p;
        p = Double.NaN;
        if (lowGerp.length > 1 && highGerp.length > 1) {
            p = st.homoscedasticTTest(lowGerp, highGerp);
        }
        values[17] = p;
        
        index = hweList.toArray();
        ratio = new double[index.length];
        sift = new double[index.length];
        gerp = new double[index.length];
        for (int i = 0; i < index.length; i++) {
            ratio[i] = snpGroup[index[i]].hetRatio;
            sift[i] = snpGroup[index[i]].sift;
            gerp[i] = snpGroup[index[i]].gerp;
        }
        values[18] = index.length;
        ds = new DescriptiveStatistics(ratio);
        mean = ds.getMean();
        sd = ds.getStandardDeviation();
        error = sd/Math.sqrt(snpGroup.length);
        values[19] = mean; values[20] = error;
        ds = new DescriptiveStatistics(sift);
        mean = ds.getMean();
        sd = ds.getStandardDeviation();
        error = sd/Math.sqrt(snpGroup.length);
        values[21] = mean; values[22] = error;
        ds = new DescriptiveStatistics(gerp);
        mean = ds.getMean();
        sd = ds.getStandardDeviation();
        error = sd/Math.sqrt(snpGroup.length);
        values[23] = mean; values[24] = error;
        double[] hweSift = sift;
        double[] hweGerp = gerp;
        st = new TTest();
        p = Double.NaN;
        if (hweSift.length > 1 && highSift.length > 1) {
            p = st.homoscedasticTTest(hweSift, highSift);
        }
        values[25] = p;
        p = Double.NaN;
        if (hweGerp.length > 1 && highGerp.length > 1) {
            p = st.homoscedasticTTest(hweGerp, highGerp);
        }
        values[26] = p;
        return values;
    }
    
    void matchSNPsByDAF2 () {
        String infileS = "E:\\Research\\cassava\\revision\\extraHets\\dele_raw_cultivated_summary.txt";
        String outfileS = "E:\\Research\\cassava\\revision\\extraHets\\extraHet_ByDAF.txt";
        double maxDaf = 0.3;
        double frequencyInterval = 0.050;
        double sampleProportion = 0.1;
        Table t = new Table (infileS);
        DeleSNP[] snps = new DeleSNP[t.getRowNumber()];
        for (int i = 0; i <  t.getRowNumber(); i++) {
            snps[i] = new DeleSNP(t.content[i]);
        }
        TDoubleArrayList boundList = new TDoubleArrayList();
        double current = 0;
        while (current < maxDaf) {
            boundList.add(current);
            current+=frequencyInterval;
        }
        boundList.add(current);
        double[] bound = boundList.toArray();
        ArrayList<DeleSNP>[] snpGroupList = new ArrayList[bound.length-1];
        DeleSNP[][] snpGroup = new DeleSNP[bound.length-1][];
        for (int i = 0; i < snpGroupList.length; i++) {
            snpGroupList[i] = new ArrayList();
        }
        for (int i = 0; i < snps.length; i++) {
            int index = Arrays.binarySearch(bound, snps[i].daf);
            if (index < 0) index = -index - 2;
            if (index >= snpGroup.length) continue;
            snpGroupList[index].add(snps[i]);
        }
        for (int i = 0; i < snpGroup.length; i++) {
            snpGroup[i] = snpGroupList[i].toArray(new DeleSNP[snpGroupList[i].size()]);
            Arrays.sort(snpGroup[i]);
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder();
            sb.append("DAFGroup\tSNPCount");
            String[] header = {"RatioMean", "RatioError","SiftMean","SiftError","GerpMean", "GerpError"};
            for (int i = 0; i < header.length; i++) {
                sb.append("\t").append("Low"+header[i]);
            }
            for (int i = 0; i < header.length; i++) {
                sb.append("\t").append("High"+header[i]);
            }
            for (int i = 0; i < header.length; i++) {
                sb.append("\t").append("HWE"+header[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < snpGroup.length; i++) {
                System.out.println(i);
                sb = new StringBuilder();
                sb.append(bound[i]).append("\t").append(snpGroup[i].length);
                int size = (int)(snpGroup[i].length*sampleProportion);
                int startIndex = 0; 
                int endIndex = startIndex+size;
                double[] values = this.getGroupEstimates(snpGroup[i], startIndex, endIndex);
                for (int j = 0; j < values.length; j++) {
                    sb.append("\t").append(values[j]);
                }
                endIndex = snpGroup[i].length;
                startIndex = endIndex-size; 
                values = this.getGroupEstimates(snpGroup[i], startIndex, endIndex);
                for (int j = 0; j < values.length; j++) {
                    sb.append("\t").append(values[j]);
                }
                Arrays.sort(snpGroup[i], new SortByDHWE());
                startIndex = 0; 
                endIndex = startIndex+size;
                values = this.getGroupEstimates(snpGroup[i], startIndex, endIndex);
                for (int j = 0; j < values.length; j++) {
                    sb.append("\t").append(values[j]);
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
    
    void matchSNPsByDAF () {
        String infileS = "E:\\Research\\cassava\\revision\\extraHets\\dele_raw_cultivated_summary.txt";
        String outfileS = "E:\\Research\\cassava\\revision\\extraHets\\extraHet_ByDAF.txt";
        int dafMax = 30;
        double sampleP = 0.1;
        Table t = new Table (infileS);
        DeleSNP[] snps = new DeleSNP[t.getRowNumber()];
        TIntHashSet freSet = new TIntHashSet();
        for (int i = 0; i <  t.getRowNumber(); i++) {
            snps[i] = new DeleSNP(t.content[i]);
            int fre = (int)(snps[i].daf*100);
            freSet.add(fre);
        }
        int[] dafInt = freSet.toArray();
        Arrays.sort(dafInt);
        ArrayList<DeleSNP>[] snpGroupList = new ArrayList[dafInt.length];
        DeleSNP[][] snpGroup = new DeleSNP[dafInt.length][];
        for (int i = 0; i < snpGroupList.length; i++) {
            snpGroupList[i] = new ArrayList();
        }
        for (int i = 0; i < snps.length; i++) {
            int fre = (int)(snps[i].daf*100);
            int index = Arrays.binarySearch(dafInt, fre);
            snpGroupList[index].add(snps[i]);
        }
        for (int i = 0; i < dafInt.length; i++) {
            snpGroup[i] = snpGroupList[i].toArray(new DeleSNP[snpGroupList[i].size()]);
            Arrays.sort(snpGroup[i]);
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder();
            sb.append("DAFGroup\tSNPCount");
            String[] header = {"RatioMean", "RatioError","SiftMean","SiftError","GerpMean", "GerpError"};
            for (int i = 0; i < header.length; i++) {
                sb.append("\t").append("Low"+header[i]);
            }
            for (int i = 0; i < header.length; i++) {
                sb.append("\t").append("High"+header[i]);
            }
            for (int i = 0; i < header.length; i++) {
                sb.append("\t").append("HWE"+header[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < dafMax; i++) {
                System.out.println(i);
                sb = new StringBuilder();
                sb.append(i).append("\t").append(snpGroup[i].length);
                int size = (int)(snpGroup[i].length*sampleP);
                int startIndex = 0; 
                int endIndex = startIndex+size;
                double[] values = this.getGroupEstimates(snpGroup[i], startIndex, endIndex);
                for (int j = 0; j < values.length; j++) {
                    sb.append("\t").append(values[j]);
                }
                endIndex = snpGroup[i].length;
                startIndex = endIndex-size; 
                values = this.getGroupEstimates(snpGroup[i], startIndex, endIndex);
                for (int j = 0; j < values.length; j++) {
                    sb.append("\t").append(values[j]);
                }
                Arrays.sort(snpGroup[i], new SortByDHWE());
                startIndex = 0; 
                endIndex = startIndex+size;
                values = this.getGroupEstimates(snpGroup[i], startIndex, endIndex);
                for (int j = 0; j < values.length; j++) {
                    sb.append("\t").append(values[j]);
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
    
    public double[] getGroupEstimates (DeleSNP[] snpGroup, int startIndex, int endIndex) {
        double[] values = new double[6];
        int size = endIndex - startIndex;
        double[] vs = new double[size];
        for (int i = 0; i < size; i++) {
            vs[i] = snpGroup[i+startIndex].hetRatio;
        }
        DescriptiveStatistics ds = new DescriptiveStatistics(vs);
        double mean = ds.getMean();
        double sd = ds.getStandardDeviation();
        double error = sd/Math.sqrt(snpGroup.length);
        values[0] = mean;values[1] = error;
        for (int i = 0; i < size; i++) {
            vs[i] = snpGroup[i+startIndex].sift;
        }
        ds = new DescriptiveStatistics(vs);
        mean = ds.getMean();
        sd = ds.getStandardDeviation();
        error = sd/Math.sqrt(snpGroup.length);
        values[2] = mean;values[3] = error;
        for (int i = 0; i < size; i++) {
            vs[i] = snpGroup[i+startIndex].gerp;
        }
        ds = new DescriptiveStatistics(vs);
        mean = ds.getMean();
        sd = ds.getStandardDeviation();
        error = sd/Math.sqrt(snpGroup.length);
        values[4] = mean;values[5] = error;
        return values;
    }
    
    class SortByDHWE implements Comparator <DeleSNP> {
	public int compare (DeleSNP o1, DeleSNP o2) {
            if (o1.dFromHWE - o2.dFromHWE < 0) return -1;
            else if (o1.dFromHWE - o2.dFromHWE > 0) return 1;
            return 0;
	}
    }
    
    class DeleSNP implements Comparable<DeleSNP> {
        public int chr;
        public int pos;
        public double sift;
        public double gerp;
        public double daf;
        public double hetRatio;
        public double p;
        public double dFromHWE;
        public int freInt;
        
        public DeleSNP (String[] content) {
            this.chr = Integer.valueOf(content[0]);
            this.pos = Integer.valueOf(content[1]);
            this.sift = Double.valueOf(content[4]);
            this.gerp = Double.valueOf(content[5]);
            this.daf = Double.valueOf(content[8]);
            this.hetRatio = Double.valueOf(content[15]);
            this.p = Double.valueOf(content[16]);
            this.dFromHWE = Math.abs(hetRatio-1);
            this.freInt = (int)(daf*100);
        }

        @Override
        public int compareTo(DeleSNP o) {
            if (this.hetRatio < o.hetRatio) return -1;
            else if (this.hetRatio > o.hetRatio) return 1;
            return 0;
        }
    }
    
    void summarizeDeleSites () {
        String siteInfoFileS = "E:\\Research\\cassava\\revision\\extraHets\\source\\DeleteriousAllelesInformation.txt";
        String imputeFileS = "E:\\Research\\cassava\\revision\\extraHets\\dele_impute_cultivated.txt";
        String rawFileS = "E:\\Research\\cassava\\revision\\extraHets\\dele_raw_cultivated.txt";
        String sumImputeFileS = "E:\\Research\\cassava\\revision\\extraHets\\dele_impute_cultivated_summary.txt";
        String sumRawFileS = "E:\\Research\\cassava\\revision\\extraHets\\dele_raw_cultivated_summary.txt";
        Table t = new Table (siteInfoFileS);
        HashMap<String, String> siteRecordMap = new HashMap();
        HashMap<String, String> siteDerivedMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            String key = t.content[i][0]+"_"+t.content[i][1];
            siteDerivedMap.put(key, t.content[i][7]);
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < 8; j++) {
                sb.append(t.content[i][j]).append("\t");
            }
            sb.deleteCharAt(sb.length()-1);
            siteRecordMap.put(key, sb.toString());
        }
        this.summarize(imputeFileS, sumImputeFileS, siteRecordMap, siteDerivedMap);
        this.summarize(rawFileS, sumRawFileS, siteRecordMap, siteDerivedMap);
    }
    
    void summarize (String inputFileS, String outputFileS, HashMap<String, String> siteRecordMap, HashMap<String, String> siteDerivedMap) {
        Table t = new Table (inputFileS);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outputFileS);
            bw.write("Chr\tPos\tREF\tALT\tSIFTscore\tGERP\tGerpTree\tDerivedAleRubber\tDAF\tOHomoDele\tOHet\tOHomoAnce\tEHomoDele\tEHet\tEHomoAnce\tRatioOVsEHets\tPValue");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                String key = t.content[i][0]+"_"+t.content[i][1];
                String derivedA = siteDerivedMap.get(key);
                int dValue = 1;
                if (t.content[i][3].equals(derivedA)) dValue = 0;
                int cnt = 0;
                double oHomoDele = 0;
                double oHet = 0;
                double oHomoAnce = 0;
                for (int j = 9; j < t.getColumnNumber(); j++) {
                    String temp = t.content[i][j].split(":")[0];
                    if (temp.startsWith(".")) continue;
                    cnt++;
                    String[] tem = null;
                    if (temp.contains("|")) {
                        tem = temp.split("\\|");
                    }
                    else {
                        tem = temp.split("\\/");
                    }
                    int[] value = new int[2];
                    for (int k = 0; k < tem.length; k++) {
                        value[k] = Integer.valueOf(tem[k]);
                    }
                    if (value[0] == dValue) {
                        if (value[1] == dValue) {
                            oHomoDele++;
                        }
                        else {
                            oHet++;
                        }
                    }
                    else {
                        if (value[1] == dValue) {
                            oHet++;
                        }
                        else {
                            oHomoAnce++;
                        }
                    }
                }
                double dFre = (oHomoDele*2+oHet)/cnt/2;
                StringBuilder sb = new StringBuilder(siteRecordMap.get(key));
                double eHomoDele = dFre*dFre;
                double eHet = 2*dFre*(1-dFre);
                double eHomoAnce = Math.pow(1-dFre, 2);
                double hetRatio = oHet/cnt/eHet;
                ChiSquareTest c = new ChiSquareTest();
                double[] expected = {eHomoDele, eHet, eHomoAnce};
                long[] observed = {(long)oHomoDele, (long)oHet, (long)oHomoAnce};
                double p = Double.NaN;
                try {
                    p = c.chiSquareTest(expected, observed);
                }
                catch (Exception e) {
                    p = Double.NaN;
                    continue;
                }
                sb.append("\t").append(dFre).append("\t").append(oHomoDele/cnt).append("\t").append(oHet/cnt).append("\t").append(oHomoAnce/cnt);
                sb.append("\t").append(eHomoDele).append("\t").append(eHet).append("\t").append(eHomoAnce).append("\t").append(hetRatio).append("\t").append(p);
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
    
    void extractCultivatedClone () {
        String infoFileS = "E:\\Research\\cassava\\revision\\extraHets\\source\\cloneInfo.txt";
        String sourceImputedFileS = "E:\\Research\\cassava\\revision\\extraHets\\source\\phasedData_delMut_GERPandTree2_filtered.vcf";
        String sourceRawFileS = "E:\\Research\\cassava\\revision\\extraHets\\source\\dele_raw.vcf";
        String outImputedFileS = "E:\\Research\\cassava\\revision\\extraHets\\dele_impute_cultivated_LAC.txt";
        String outRawFileS = "E:\\Research\\cassava\\revision\\extraHets\\dele_raw_cultivated_LAC.txt";
        ArrayList<String> cultiList = new ArrayList();
        Table t = new Table (infoFileS);
//        for (int i = 0; i < t.getRowNumber(); i++) {
//            if (!t.content[i][13].equals("Culti")) continue;
//            cultiList.add(t.content[i][1]);
//        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (!t.content[i][12].equals("LAC")) continue;
            cultiList.add(t.content[i][1]);
        }
        String[] culti = cultiList.toArray(new String[cultiList.size()]);
        Arrays.sort(culti);
        this.exportCulti(sourceImputedFileS, outImputedFileS, culti);
        this.exportCulti(sourceRawFileS, outRawFileS, culti);
    }
    
    void exportCulti (String infileS, String outfileS, String[] cloneName) {
        Table t = new Table(infileS);
        TIntArrayList indices = new TIntArrayList();
        for (int i = 0; i < 9; i++) indices.add(i);
        for (int i = 9; i < t.header.length; i++) {
            if (Arrays.binarySearch(cloneName, t.header[i]) < 0) continue;
            indices.add(i);
        }
        int[] index = indices.toArray();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < index.length; i++) {
                sb.append(t.header[index[i]]).append("\t");
            }
            sb.deleteCharAt(sb.length()-1);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                sb = new StringBuilder();
                for (int j = 0; j < index.length; j++) {
                    sb.append(t.content[i][index[j]]).append("\t");
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
    
    void extractUnimputatedVCF () {
        String infileDirS = "/home/rp444/WGS/preFilterBAMs/hapmap_fei/rawData_filt_noIndels/";
        String siteInfoFileS = "/local/workdir/mingh/phasedData_delMut_GERPandTree2_filtered.vcf";
        String outfileS = "/local/workdir/mingh/dele_raw.vcf";
        int chrNum = 18;
        TIntArrayList[] siteList = new TIntArrayList[chrNum];
        int[][] sites = new int[chrNum][];
        for (int i = 0; i < siteList.length; i++) {
            siteList[i] = new TIntArrayList();
        }
        try {
            BufferedReader br = IoUtils.getTextReader(siteInfoFileS);
            String temp = br.readLine();
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.substring(0, 100).split("\t");
                siteList[Integer.valueOf(tem[0])-1].add(Integer.valueOf(tem[1]));
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        for (int i = 0; i < sites.length; i++) {
            sites[i] = siteList[i].toArray();
            Arrays.sort(sites[i]);
        }
        File[] fs = new File(infileDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        ArrayList<String>[] outputList = new ArrayList[chrNum];
        fList.parallelStream().forEach(f -> {
            int chr = Integer.valueOf(f.getName().replaceFirst("_noIndels_st01.vcf.gz", "").replaceFirst("chr", ""));
            int chrIndex = chr-1;
            outputList[chrIndex] = new ArrayList();
            int cnt = 0;
            try {
                BufferedReader br = IoUtils.getTextGzipReader(f.getAbsolutePath());
                String temp = br.readLine();
                outputList[chrIndex].add(temp);
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%100000 == 0) System.out.println(String.valueOf(cnt)+"\t"+f.getName());
                    String[] tem = temp.substring(0, 100).split("\t");
                    int index = Arrays.binarySearch(sites[chrIndex], Integer.valueOf(tem[1]));
                    if (index < 0) continue;
                    outputList[chrIndex].add(temp);
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write(outputList[0].get(0));
            bw.newLine();
            for (int i = 0; i < outputList.length; i++) {
                for (int j = 1; j < outputList[i].size(); j++) {
                    bw.write(outputList[i].get(j));
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
}
