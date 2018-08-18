/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maf;

import format.Table;
import java.io.BufferedReader;
import java.io.BufferedWriter;
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
class Pathway {
    
    public Pathway () {
        this.delAndPathway();
    }
    
    void delAndPathway () {
        String pathWayFileS = "M:\\production\\maf\\pathway\\source\\MAPMAM_Anotation_Corn.txt";
        String highConfidenceGeneFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\highConfidence_transcript.txt";
        String outfileS = "M:\\production\\maf\\pathway\\delPathway.txt";
        int level = 2;
        int minSize = 10;
        Table t = new Table (highConfidenceGeneFileS);
        HashMap<String, String> geneDelMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getDoubleValue(i, 4) == 0) continue;
            String name = t.content[i][0];
            if (name.startsWith("GRM")) {
                name = name.split("_")[0];
            }
            StringBuilder sb = new StringBuilder();
            double syn = t.getDoubleValue(i, 6);
            double delS = t.getDoubleValue(i, 11);
            double delSG = t.getDoubleValue(i, 13);
            double ratioS = Double.NaN;
            double ratioSG = Double.NaN;
            if (syn != 0) {
                ratioS = delS/syn;
                ratioSG = delSG/syn;
            }
            sb.append(syn).append("\t").append(delS).append("\t").append(ratioS).append("\t").append(delSG).append("\t").append(ratioSG);
            geneDelMap.put(name, sb.toString());
        }
        Set<String> delGeneSet = geneDelMap.keySet();
        String[] classes = null;
        HashSet<String>[] geneSet = null;
        String[][] geneArray = null;
        try {
            HashSet<String> classSet = new HashSet();
            BufferedReader br = IoUtils.getTextReader(pathWayFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                List<String> l = FStringUtils.fastSplit(temp);
                String[] tem = l.get(1).split("\\.");
                if (tem.length < level) continue;
                tem = l.get(2).split("\\.");
                String className = tem[0];
                for (int i = 1; i < level; i++) {
                    className = className+"."+tem[i];
                }
                classSet.add(className);
            }
            br.close();
            classes = classSet.toArray(new String[classSet.size()]);
            Arrays.sort(classes);
            geneSet = new HashSet[classes.length];
            geneArray = new String[classes.length][];
            for (int i = 0; i < geneSet.length; i++) {
                geneSet[i] = new HashSet();
            }
            br = IoUtils.getTextReader(pathWayFileS);
            while ((temp = br.readLine()) != null) {
                List<String> l = FStringUtils.fastSplit(temp);
                String[] tem = l.get(1).split("\\.");
                if (tem.length < level) continue;
                tem = l.get(2).split("\\.");
                String className = tem[0];
                for (int i = 1; i < level; i++) {
                    className = className+"."+tem[i];
                }
                int index = Arrays.binarySearch(classes, className);
                geneSet[index].add(l.get(0));
            }
            br.close();
            for (int i = 0; i < geneArray.length; i++) {
                geneArray[i] = geneSet[i].toArray(new String[geneSet[i].size()]);
                Arrays.sort(geneArray[i]);
            }
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Pathway\tNumOfGenes\tSyn\tSynError\tDelS\tDelSError\tRatioDelSVsSyn\tRatioDelSVsSynError\tDelSG\tDelSGError\tRatioDelSGVsSyn\tRatioDelSGVsSynError");
            bw.newLine();
            for (int i = 0; i < geneArray.length; i++) {
                HashSet<String> set = new HashSet();
                for (int j = 0; j < geneArray[i].length; j++) {
                    if (delGeneSet.contains(geneArray[i][j])) {
                        set.add(geneArray[i][j]);
                    }
                }
                if (set.size() < minSize) continue;
                String[] genes = set.toArray(new String[set.size()]);
                Arrays.sort(genes);
                StringBuilder sb = new StringBuilder(classes[i]);
                sb.append("\t").append(geneArray[i].length);
                double[][] values = new double[5][genes.length];
                for (int j = 0; j < genes.length; j++) {
                    String[] tem = geneDelMap.get(genes[j]).split("\t");
                    for (int k = 0; k < values.length; k++) {
                        values[k][j] = Double.valueOf(tem[k]);
                    }
                }
                for (int j = 0; j < values.length; j++) {
                    DescriptiveStatistics d = new DescriptiveStatistics(FArrayUtils.removeNaN(values[j]));
                    double mean = d.getMean();
                    double sd = d.getStandardDeviation();
                    double se = sd/Math.sqrt(values[j].length);
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
}
