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
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.TTest;
import utils.FArrayUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class Pollen {
    
    public Pollen () {
        this.countDelAndSynEachClass();
        this.mkMeanSEAndPValue();
        
    }
    
    private void mkMeanSEAndPValue () {
        String infileDirS = "M:\\production\\maf\\pollen\\classCount\\";
        String outDirS = "M:\\production\\maf\\pollen\\meanSEP\\";
        new File(outDirS).mkdir();
        File[] fs = new File (infileDirS).listFiles();
        Table[] ts = new Table[fs.length];
        String[] tissue = new String[fs.length];
        for (int i = 0; i < ts.length; i++) {
            tissue[i] = fs[i].getName().replaceFirst(".txt", "");
            ts[i] = new Table(fs[i].getAbsolutePath());
        }
        for (int i = 0; i < ts[0].getColumnNumber()-1; i++) {
            String oufileS1 = new File (outDirS, ts[0].header[i+1]+".meanSE.txt").getAbsolutePath();
            String oufileS2 = new File (outDirS, ts[0].header[i+1]+".pValue.txt").getAbsolutePath();
            try {
                BufferedWriter bw = IoUtils.getTextWriter(oufileS1);
                StringBuilder sb = new StringBuilder("Statistics");
                double[][] value = new double[tissue.length][2];
                for (int j = 0; j < tissue.length; j++) {
                    sb.append("\t").append(tissue[j]);
                    double[] vs = new double[ts[j].getRowNumber()];
                    for (int k = 0; k < vs.length; k++) {
                        vs[k] = ts[j].getDoubleValue(k, i+1);  
                    }
                    vs = FArrayUtils.removeNaN(vs);
                    DescriptiveStatistics d = new DescriptiveStatistics(vs);
                    value[j][0] = d.getMean();
                    value[j][1] = d.getStandardDeviation()/Math.sqrt(vs.length);
                }
                bw.write(sb.toString());
                bw.newLine();
                for (int j = 0; j < value[0].length; j++) {
                    sb = new StringBuilder();
                    if (j%2 == 0) {
                        sb.append("Mean");
                    }
                    else {
                        sb.append("SE");
                    }
                    for (int k = 0; k < value.length; k++) {
                        sb.append("\t").append(value[k][j]);
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                bw = IoUtils.getTextWriter(oufileS2);
                double[][] values = new double[tissue.length][];
                for (int j = 0; j < tissue.length; j++) {
                    sb.append("\t").append(tissue[j]);
                    double[] vs = new double[ts[j].getRowNumber()];
                    for (int k = 0; k < vs.length; k++) {
                        vs[k] = ts[j].getDoubleValue(k, i+1);  
                    }
                    vs = FArrayUtils.removeNaN(vs);
                    values[j] = vs;
                }
                sb = new StringBuilder("Class");
                for (int j = 0; j < tissue.length; j++) {
                    sb.append("\t").append(tissue[j]);
                }
                bw.write(sb.toString());
                bw.newLine();
                double[][] matrix = new double[tissue.length][tissue.length];
                for (int j = 0; j < tissue.length; j++) {
                    for (int k = 0; k < tissue.length; k++) {
                        matrix[j][k] = new TTest().tTest(values[j], values[k]);
                    }
                }
                for (int j = 0; j < matrix.length; j++) {
                    sb = new StringBuilder(tissue[j]);
                    for (int k = 0; k < tissue.length; k++) {
                        sb.append("\t").append(matrix[j][k]);
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
    }
    
    private void countDelAndSynEachClass () {
        String geneModelS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\highConfidence_transcript.txt";
        String expressionDirS = "M:\\production\\maf\\pollen\\source\\expression\\";
        String delSNPFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\class\\Non_Synonymous_Deleterious.txt";
        String synSNPFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\class\\Synonymous.txt";
        String geneFeatureFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
        String outfileDirS = "M:\\production\\maf\\pollen\\classCount\\";
        new File (outfileDirS).mkdir();
        Table t = new Table (geneModelS);
        ArrayList<String> tranList = new ArrayList();
        HashMap<String, Double> geneLengthMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.content[i][4].equals("0")) continue;
            String gene = t.content[i][0];
            if (gene.startsWith("GRM")) {
                gene = gene.split("_")[0];
            }
            tranList.add(gene);
            geneLengthMap.put(gene, Double.valueOf(t.content[i][1]));
        }
        String[] confiGene = tranList.toArray(new String[tranList.size()]);
        Arrays.sort(confiGene);
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        t = new Table (delSNPFileS);
        TIntArrayList cList = new TIntArrayList();
        TIntArrayList pList = new TIntArrayList();
        TDoubleArrayList freList= new TDoubleArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            cList.add(t.getIntValue(i, 0));
            pList.add(t.getIntValue(i, 1));
            freList.add(t.getDoubleValue(i, 3));
        }
        int[] delChr = cList.toArray();
        int[] delPos = pList.toArray();
        double[] delFre = freList.toArray();
        t = new Table (synSNPFileS);
        cList = new TIntArrayList();
        pList = new TIntArrayList();
        freList= new TDoubleArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            cList.add(t.getIntValue(i, 0));
            pList.add(t.getIntValue(i, 1));
            freList.add(t.getDoubleValue(i, 3));
        }
        int[] synChr = cList.toArray();
        int[] synPos = pList.toArray();
        double[] synFre = freList.toArray();
        File[] fs = new File(expressionDirS).listFiles();
        
        for (int i = 0; i < fs.length; i++) {
            String tissue = fs[i].getName().replace(".txt", "");
            t = new Table (fs[i].getAbsolutePath());
            ArrayList<String> geneList = new ArrayList();
            for (int j = 0; j < t.getRowNumber(); j++) {
                int index = Arrays.binarySearch(confiGene, t.content[j][0]);
                if (index < 0) continue;
                geneList.add(t.content[j][0]);
            }
            String[] genes = geneList.toArray(new String[geneList.size()]);
            Arrays.sort(genes);
            double[] synCount = new double[t.getRowNumber()];
            double[] delCount = new double[t.getRowNumber()];
            double[] synFreCount = new double[t.getRowNumber()];
            double[] delFreCount = new double[t.getRowNumber()];
            for (int j = 0; j < delChr.length; j++) {
                int geneIndex = gf.getGeneIndex(delChr[j], delPos[j]);
                if (geneIndex < 0) continue;
                int cdsIndex = gf.getCDSIndex(geneIndex, 0, delChr[j], delPos[j]);
                if (cdsIndex < 0) continue;
                int index = Arrays.binarySearch(genes, gf.getGeneName(geneIndex));
                if (index < 0) continue;
                delCount[index]+=1/geneLengthMap.get(gf.getGeneName(geneIndex));
                delFreCount[index]+=delFre[j]/geneLengthMap.get(gf.getGeneName(geneIndex));
            }
            for (int j = 0; j < synChr.length; j++) {
                int geneIndex = gf.getGeneIndex(synChr[j], synPos[j]);
                if (geneIndex < 0) continue;
                int cdsIndex = gf.getCDSIndex(geneIndex, 0, synChr[j], synPos[j]);
                if (cdsIndex < 0) continue;
                int index = Arrays.binarySearch(genes, gf.getGeneName(geneIndex));
                if (index < 0) continue;
                synCount[index]+=1/geneLengthMap.get(gf.getGeneName(geneIndex));
                synFreCount[index]+=synFre[j]/geneLengthMap.get(gf.getGeneName(geneIndex));
            }
            String outfileS = new File (outfileDirS, tissue+".txt").getAbsolutePath();
            try {
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                bw.write("Genes\tSynPerSite\tDeleteriousPerSite\tDeleteriousRatioPerSite\tSynPerSitePerLine\tDeleteriousPerSitePerLine\tDeleteriousRatioPerSitePerLine");
                bw.newLine();
                for (int j = 0; j < genes.length; j++) {
                    StringBuilder sb = new StringBuilder(genes[j]);
                    sb.append("\t");
                    if (synCount[j] != 0) {
                        sb.append(synCount[j]).append("\t").append(delCount[j]).append("\t").append(delCount[j]/synCount[j]).append("\t");
                        sb.append(synFreCount[j]).append("\t").append(delFreCount[j]).append("\t").append(delFreCount[j]/synFreCount[j]);
                    }
                    else {
                        sb.append(synCount[j]).append("\t").append(delCount[j]).append("\t").append(Double.NaN).append("\t");
                        sb.append(synFreCount[j]).append("\t").append(delFreCount[j]).append("\t").append(Double.NaN);
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
    }
    
    /**
     * @deprecated 
     */
    private void countDelAndSynEachGene () {
        String geneModelS = "M:\\production\\maf\\annotations\\siftScore\\005_hmp32SNPClassMAF\\highConfidence_transcript.txt";
        String expressionDirS = "M:\\production\\maf\\pollen\\source\\expression\\";
        String delSNPFileS = "M:\\production\\maf\\annotations\\siftScore\\005_hmp32SNPClassMAF\\class\\Non_Synonymous_Deleterious.txt";
        String synSNPFileS = "M:\\production\\maf\\annotations\\siftScore\\005_hmp32SNPClassMAF\\class\\Synonymous.txt";
        String geneFeatureFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
        String outfileS = "M:\\production\\maf\\pollen\\delSynPollenEachGene.txt";
        Table t = new Table (geneModelS);
        String[] confiGene = new String[t.getRowNumber()];
        HashMap<String, Double> geneLengthMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            String gene = t.content[i][0];
            if (gene.startsWith("GRM")) {
                gene = gene.split("_")[0];
            }
            confiGene[i] = gene;
            geneLengthMap.put(gene, Double.valueOf(t.content[i][1]));
        }
        Arrays.sort(confiGene);
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        t = new Table (delSNPFileS);
        TIntArrayList cList = new TIntArrayList();
        TIntArrayList pList = new TIntArrayList();
        TDoubleArrayList freList= new TDoubleArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            cList.add(t.getIntValue(i, 0));
            pList.add(t.getIntValue(i, 1));
            freList.add(t.getDoubleValue(i, 3));
        }
        int[] delChr = cList.toArray();
        int[] delPos = pList.toArray();
        double[] delFre = freList.toArray();
        t = new Table (synSNPFileS);
        cList = new TIntArrayList();
        pList = new TIntArrayList();
        freList= new TDoubleArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            cList.add(t.getIntValue(i, 0));
            pList.add(t.getIntValue(i, 1));
            freList.add(t.getDoubleValue(i, 3));
        }
        int[] synChr = cList.toArray();
        int[] synPos = pList.toArray();
        double[] synFre = freList.toArray();
        File[] fs = new File(expressionDirS).listFiles();
        String[] tissue = new String[fs.length];
        
        double[][] synCount = new double[fs.length][];
        double[][] delCount = new double[fs.length][];
        double[][] ratioCount = new double[fs.length][];
        double[][] synFreCount = new double[fs.length][];
        double[][] delFreCount = new double[fs.length][];
        double[][] ratioFreCount = new double[fs.length][];
        
        for (int i = 0; i < fs.length; i++) {
            tissue[i] = fs[i].getName().replace(".txt", "");
            t = new Table (fs[i].getAbsolutePath());
            ArrayList<String> geneList = new ArrayList();
            for (int j = 0; j < t.getRowNumber(); j++) {
                int index = Arrays.binarySearch(confiGene, t.content[j][0]);
                if (index < 0) continue;
                geneList.add(t.content[j][0]);
            }
            String[] genes = geneList.toArray(new String[geneList.size()]);
            Arrays.sort(genes);
            synCount[i] = new double[t.getRowNumber()];
            delCount[i] = new double[t.getRowNumber()];
            synFreCount[i] = new double[t.getRowNumber()];
            delFreCount[i] = new double[t.getRowNumber()];
            for (int j = 0; j < delChr.length; j++) {
                int index = gf.getGeneIndex(delChr[j], delPos[j]);
                if (index < 0) continue;
                int geneIndex = index;
                index = gf.getCDSIndex(index, 0, delChr[j], delPos[j]);
                if (index < 0) continue;
                index = Arrays.binarySearch(genes, gf.getGeneName(geneIndex));
                if (index < 0) continue;
                
                delCount[i][index]+=1/geneLengthMap.get(gf.getGeneName(geneIndex));
                delFreCount[i][index]+=delFre[j]/geneLengthMap.get(gf.getGeneName(geneIndex));
            }
            for (int j = 0; j < synChr.length; j++) {
                int index = gf.getGeneIndex(synChr[j], synPos[j]);
                int geneIndex = index;
                if (index < 0) continue;
                index = gf.getCDSIndex(index, 0, synChr[j], synPos[j]);
                if (index < 0) continue;
                index = Arrays.binarySearch(genes, gf.getGeneName(geneIndex));
                if (index < 0) continue;
                synCount[i][index]+=1/geneLengthMap.get(gf.getGeneName(geneIndex));
                synFreCount[i][index]+=synFre[j]/geneLengthMap.get(gf.getGeneName(geneIndex));
            }
            TDoubleArrayList ratioCountList = new TDoubleArrayList();
            TDoubleArrayList ratioFreCountList = new TDoubleArrayList();
            for (int j = 0; j < t.getRowNumber(); j++) {
                if (synCount[i][j] == 0) continue;
                ratioCountList.add(delCount[i][j]/synCount[i][j]);
                ratioFreCountList.add(delFreCount[i][j]/synFreCount[i][j]);
            }
            ratioCount[i] = ratioCountList.toArray();
            ratioFreCount[i] = ratioFreCountList.toArray();
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder("Class");
            for (int i = 0; i < tissue.length; i++) {
                sb.append("\t").append(tissue[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            
            sb = new StringBuilder("Synonomous_Mean");
            for (int i = 0; i < tissue.length; i++) {
                DescriptiveStatistics d = new DescriptiveStatistics(synCount[i]);
                for (int j = 0; j < synCount[i].length; j++) {
                }
                sb.append("\t").append(d.getMean());
            }
            bw.write(sb.toString());
            bw.newLine();
            sb = new StringBuilder("Synonomous_SD");
            for (int i = 0; i < tissue.length; i++) {
                DescriptiveStatistics d = new DescriptiveStatistics(synCount[i]);
                sb.append("\t").append(d.getStandardDeviation());
            }
            bw.write(sb.toString());
            bw.newLine();
            sb = new StringBuilder("Deleterious_Mean");
            for (int i = 0; i < tissue.length; i++) {
                DescriptiveStatistics d = new DescriptiveStatistics(delCount[i]);
                sb.append("\t").append(d.getMean());
            }
            bw.write(sb.toString());
            bw.newLine();
            sb = new StringBuilder("Deleterious_SD");
            for (int i = 0; i < tissue.length; i++) {
                DescriptiveStatistics d = new DescriptiveStatistics(delCount[i]);
                sb.append("\t").append(d.getStandardDeviation());
            }
            bw.write(sb.toString());
            bw.newLine();
            sb = new StringBuilder("DelRatio_Mean");
            for (int i = 0; i < tissue.length; i++) {
                DescriptiveStatistics d = new DescriptiveStatistics(ratioCount[i]);
                sb.append("\t").append(d.getMean());
            }
            bw.write(sb.toString());
            bw.newLine();
            sb = new StringBuilder("DelRatio_SD");
            for (int i = 0; i < tissue.length; i++) {
                DescriptiveStatistics d = new DescriptiveStatistics(ratioCount[i]);
                sb.append("\t").append(d.getStandardDeviation());
            }
            bw.write(sb.toString());
            bw.newLine();
            
           sb = new StringBuilder("SynonomousPerLine_Mean");
            for (int i = 0; i < tissue.length; i++) {
                DescriptiveStatistics d = new DescriptiveStatistics(synFreCount[i]);
                sb.append("\t").append(d.getMean());
            }
            bw.write(sb.toString());
            bw.newLine();
            sb = new StringBuilder("SynonomousPerLine_SD");
            for (int i = 0; i < tissue.length; i++) {
                DescriptiveStatistics d = new DescriptiveStatistics(synFreCount[i]);
                sb.append("\t").append(d.getStandardDeviation());
            }
            bw.write(sb.toString());
            bw.newLine();
            sb = new StringBuilder("DeleteriousPerLine_Mean");
            for (int i = 0; i < tissue.length; i++) {
                DescriptiveStatistics d = new DescriptiveStatistics(delFreCount[i]);
                sb.append("\t").append(d.getMean());
            }
            bw.write(sb.toString());
            bw.newLine();
            sb = new StringBuilder("DeleteriousPerLine_SD");
            for (int i = 0; i < tissue.length; i++) {
                DescriptiveStatistics d = new DescriptiveStatistics(delFreCount[i]);
                sb.append("\t").append(d.getStandardDeviation());
            }
            bw.write(sb.toString());
            bw.newLine();
            sb = new StringBuilder("DelRatioPerLine_Mean");
            for (int i = 0; i < tissue.length; i++) {
                DescriptiveStatistics d = new DescriptiveStatistics(ratioFreCount[i]);
                sb.append("\t").append(d.getMean());
            }
            bw.write(sb.toString());
            bw.newLine();
            sb = new StringBuilder("DelRatioPerLine_SD");
            for (int i = 0; i < tissue.length; i++) {
                DescriptiveStatistics d = new DescriptiveStatistics(ratioFreCount[i]);
                sb.append("\t").append(d.getStandardDeviation());
            }
            bw.write(sb.toString());
            bw.newLine();
            
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
    private void countDeleterious () {
        String infileS = "M:\\production\\maf\\annotations\\siftScore\\005_hmp32SNPClassMAF\\highConfidence_transcript.txt";
        String expressionDirS = "M:\\production\\maf\\pollen\\source\\expression\\";
        String outfileS = "M:\\production\\maf\\pollen\\delCountPollen.txt";
        Table t = new Table(infileS);
        HashMap<String, Integer> geneSynMap = new HashMap();
        HashMap<String, Integer> geneDeleDMap = new HashMap();
        HashMap<String, Integer> geneLengthMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            String gene = t.content[i][0];
            if (gene.startsWith("GRM")) {
                gene = gene.split("_")[0];
            }
            geneLengthMap.put(gene, Integer.valueOf(t.content[i][1]));
            geneSynMap.put(gene, Integer.valueOf(t.content[i][2]));
            geneDeleDMap.put(gene, Integer.valueOf(t.content[i][7]));
        }
        File[] fs = new File(expressionDirS).listFiles();
        double[] syn = new double[fs.length];
        double[] del = new double[fs.length];
        double[] ratio = new double[fs.length];
        for (int i = 0; i < fs.length; i++) {
            t = new Table(fs[i].getAbsolutePath());
            double sumS = 0;
            double sumD = 0;
            double sumL = 0;
            int cnt = 0;
            for (int j = 0; j < t.getRowNumber(); j++) {
                Integer s = geneSynMap.get(t.content[j][0]);
                Integer l = geneLengthMap.get(t.content[j][0]);
                Integer c = geneDeleDMap.get(t.content[j][0]);
                if (s == null) continue;
                sumS+=s;
                sumD+=c;
                sumL+=l;
                cnt++;
            }
            syn[i] = (double)sumS/sumL;
            del[i] = (double)sumD/sumL;
            ratio[i] = del[i]/syn[i];
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Tissue\tSynonomous\tDeleterious\tDelRatio");
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                StringBuilder sb = new StringBuilder(fs[i].getName().replace(".txt", ""));
                sb.append("\t").append(syn[i]).append("\t").append(del[i]).append("\t").append(ratio[i]);
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
