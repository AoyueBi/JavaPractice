/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maf;

import format.Table;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.set.hash.TIntHashSet;
import graphcis.r.DensityPlot;
import graphcis.r.ScatterPlot;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import utils.FArrayUtils;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class MafAnnotation {
    //A, C, D, G, I, T
    byte[] possibleAllele = {65, 67, 68, 71, 73, 84};
    byte[] allele = {65, 67, 71, 84};
    
    public MafAnnotation () {
        //this.chr10Test();
        this.MafPipe();
    }
    
    public void MafPipe () {
        this.attributeToAnnotation();
    }
    
    public void attributeToAnnotation() {
        String infileDirS = "M:\\production\\maf\\annotations\\maf\\mafAttribute\\";
        String outfileDirS = "M:\\production\\maf\\annotations\\maf\\mafAnnotation\\data\\";
        File[] fs = new File (infileDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String outfileS = f.getName().split("_")[0]+"_maf.txt.gz";
            outfileS = new File (outfileDirS, outfileS).getAbsolutePath();
            try {
                BufferedReader br = null;
                if (f.getName().endsWith("gz")) br = IoUtils.getTextGzipReader(f.getAbsolutePath());
                else br = IoUtils.getTextReader(f.getAbsolutePath());
                BufferedWriter bw = IoUtils.getTextGzipWriter(outfileS);
                //BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                bw.write("Chr\tPos\tRefAllele\tTotalDepth\tAllelePresence\tMinorAllele\tMAF");
                bw.newLine();
                int cnt = 0;
                String temp = br.readLine();
                while ((temp = br.readLine()) != null) {
                    bw.write(this.getMafString(temp));
                    bw.newLine();
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt)+" sites processed from " + f.getName());
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
    
    public void chr10Test () {
        //this.getRandomSetChr10();
        //this.attributeDistribution();
        //this.extractMaf();
        //this.checkDensity();
    }
    
    private void checkDensity () {
        String infileS = "M:\\production\\maf\\annotations\\maf\\mafAnnotation\\randomChr10\\chr10_sub.maf.txt";
        TDoubleArrayList l = new TDoubleArrayList();
        Table t = new Table (infileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.content[i][6].startsWith("N")) continue;
            l.add(Double.valueOf(t.content[i][6]));
        }
        double[] value = l.toArray();
        DensityPlot d = new DensityPlot(value);
        d.setSmoothN(100000);
        d.setXLim(0, 0.005);
        d.showGraph();
    }
    
    
    
    public void extractMaf () {
        String infileS = "M:\\production\\maf\\annotations\\maf\\mafAnnotation\\randomChr10\\chr10_sub.mafAttribute.txt";
        String outfileS = "M:\\production\\maf\\annotations\\maf\\mafAnnotation\\randomChr10\\chr10_sub.maf.txt";
        String temp = null;
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Chr\tPos\tRefAllele\tTotalDepth\tAllelePresence\tMinorAllele\tMAF");
            bw.newLine();
            temp = br.readLine();
            while ((temp = br.readLine()) != null) {
                bw.write(this.getMafString(temp));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            System.out.println(temp);
            e.printStackTrace();
        }
    }
    
    private String getMafString (String line) {
        List<String> split = FStringUtils.fastSplit(line, "\t");
        StringBuilder sb = new StringBuilder();
        sb.append(split.get(0)).append("\t").append(split.get(1)).append("\t").append(split.get(2)).append("\t");
        if (split.get(3).startsWith("N")) {
            sb.append("NA\tNA\tNA\tNA");
            return sb.toString();
        }
        int allelePresence = Integer.valueOf(split.get(4));
        sb.append(split.get(3)).append("\t").append(split.get(4)).append("\t");
        int[] alleleDepth = new int[4];
        int[] alleleOnce = new int[4];
        int[] alleleTwice = new int[4];
        int[] alleleHets = new int[4];
        double[] hetMean = new double[4];
        int cnt = 0;
        for (int i = 0; i < +this.possibleAllele.length; i++) {
            if (i == 2 || i == 4) continue;
            alleleDepth[cnt] = Integer.valueOf(split.get(i+7));
            cnt++;
        }
        cnt = 0;
        for (int i = 0; i < +this.possibleAllele.length; i++) {
            if (i == 2 || i == 4) continue;
            alleleOnce[cnt] = Integer.valueOf(split.get(i+13));
            cnt++;
        }
        cnt = 0;
        for (int i = 0; i < +this.possibleAllele.length; i++) {
            if (i == 2 || i == 4) continue;
            alleleTwice[cnt] = Integer.valueOf(split.get(i+19));
            cnt++;
        }
        cnt = 0;
        for (int i = 0; i < +this.possibleAllele.length; i++) {
            if (i == 2 || i == 4) continue;
            alleleHets[cnt] = Integer.valueOf(split.get(i+31));
            cnt++;
        }
        cnt = 0;
        for (int i = 0; i < +this.possibleAllele.length; i++) {
            if (i == 2 || i == 4) continue;
            hetMean[cnt] = Double.valueOf(split.get(i+37));
            cnt++;
        }
        int[] pIndex = FArrayUtils.getIndexByDescendingValue(alleleOnce);
        if (alleleOnce[pIndex[1]] == alleleOnce[pIndex[2]]) {
            if (alleleDepth[pIndex[1]] > alleleDepth[pIndex[2]]) {
                int tempInt = pIndex[1];
                pIndex[1] = pIndex[2];
                pIndex[2] = tempInt;
            }
        }
        int majorIndex  = pIndex[0];
        int minorIndex = pIndex[1];
        double maf = -1;
        double hetFrequency = -1;
        byte minorAllele = Byte.MIN_VALUE;
        if (alleleTwice[minorIndex] > 0) {
            double hetCount = alleleHets[minorIndex];
            double errorCount = (double)(alleleOnce[pIndex[2]]+alleleOnce[pIndex[3]])/2;
            double errorRatio = (double)(errorCount)/alleleOnce[minorIndex];
            double minorCount = (double)alleleOnce[minorIndex]* (1-errorRatio);
            hetCount = hetCount*(1-errorRatio);
            double homoMinorCount = minorCount-hetCount;
            double homoMajorCount = (double)alleleOnce[majorIndex]-hetCount; 
            double minor = 2*homoMinorCount+hetCount;
            double major = 2*homoMajorCount+hetCount;
            maf = minor/(major+minor);
            minorAllele = this.allele[minorIndex];
        }
        else {
            if (allelePresence == alleleOnce[majorIndex]) {
                maf = 0;
            }
            else {
                minorAllele = this.allele[minorIndex];
                maf = (double)(allelePresence-alleleOnce[majorIndex])/allelePresence;
            }
        }
        if (minorAllele == Byte.MIN_VALUE) {
            sb.append("NA").append("\t");
            if (allelePresence == 0) sb.append("NaN");
            else sb.append((float)maf);
        }
        else {
            sb.append((char)minorAllele).append("\t");
            if (allelePresence == 0) sb.append("NaN");
            else sb.append((float)maf);
            
        }
        return sb.toString();
    }
    
    public void attributeDistribution () {
        String infileS = "M:\\production\\maf\\annotations\\maf\\mafAnnotation\\randomChr10\\chr10_sub.mafAttribute.txt";
        String outfileDirS = "M:\\production\\maf\\annotations\\maf\\mafAnnotation\\randomChr10\\distribution\\";
        Table t = new Table (infileS);
        for (int i = 0; i < t.getColumnNumber(); i++) {
            TDoubleArrayList l = new TDoubleArrayList();
            for (int j = 0; j < t.getRowNumber(); j++) {
                char c = t.content[j][i].charAt(0);
                if (c < 48 || c > 57) continue;
                l.add(t.getDoubleValue(j, i));
            }
            double[] value = l.toArray();
            String outfileS = new File(outfileDirS, t.header[i]+".pdf").getAbsolutePath();
            DensityPlot d =new DensityPlot(value);
            if (t.header[i].startsWith("TotalDepth")) {
                d.setSmoothN(100000);
                d.setXLim(0, 5000);
            }
            else if (t.header[i].startsWith("AlleleDepth")){
                d.setSmoothN(1000000);
                d.setXLim(0, 100);
            }
            else if (t.header[i].startsWith("AlleleHetPresence")){
                d.setSmoothN(10000);
                d.setXLim(0, 100);
            }
            else if (t.header[i].startsWith("AllelePresence")){
                d.setSmoothN(10000);
                if (t.header[i].equals("AllelePresence")) {
                    
                }
                else {
                    d.setXLim(0, 400);
                }
                
            }
            d.setXLab(t.header[i]);
            d.saveGraph(outfileS);
        }
    }
    
    public void getRandomSetChr10 () {
        int size = 10000;
        int chrLength = 149627545;
        String infileS = "M:\\production\\maf\\annotations\\maf\\mafAttribute\\chr010_mafAttribute.txt";
        String oufileS = "M:\\production\\maf\\annotations\\maf\\mafAnnotation\\randomChr10\\chr10_sub.mafAttribute.txt";
        int[] index = FArrayUtils.getRandomIntArray(chrLength, size);
        Arrays.sort(index);
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            BufferedWriter bw = IoUtils.getTextWriter(oufileS);
            int cnt = -1;
            String temp = br.readLine();
            bw.write(temp);
            bw.newLine();
            while ((temp = br.readLine()) != null) {
                cnt++;
                if (cnt%1000000 == 0) System.out.println(cnt);
                if (Arrays.binarySearch(index, cnt) < 0) continue;
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
    
}
