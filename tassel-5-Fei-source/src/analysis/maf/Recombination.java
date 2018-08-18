/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maf;

import format.GeneFeature;
import format.Range;
import format.Table;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import utils.FArrayUtils;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class Recombination {
    
    Recombination () {
        //this.makeRecombinationPointTable();
        //this.mkBinTable();
    }
    
    void mkBinTable () {
//        String transcriptFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\highConfidence_transcript.txt";
//        String crossOverFileS = "M:\\production\\maf\\recombination\\recombinationPoint.txt";
//        String deleteriousSNPFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\class\\Non_Synonymous_Deleterious_High_Gerp.txt";
//        String synSNPFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\class\\Synonymous.txt";
//        String geneFeatureFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
//        String infoFileS = "E:\\Database\\InfoFile\\ChrLenCentPosi_agpV3.txt";
//        String outfileS = "M:\\production\\maf\\recombination\\binTable.txt";
        String transcriptFileS = "M:\\production\\maf\\001_variantDiscovery\\005_snpClass\\highConfidence_transcript.txt";
        String crossOverFileS = "M:\\production\\maf\\002_recombination\\recombinationPoint.txt";
        String deleteriousSNPFileS = "M:\\production\\maf\\001_variantDiscovery\\005_snpClass\\class\\Non_Synonymous_Deleterious_High_Gerp.txt";
        String synSNPFileS = "M:\\production\\maf\\001_variantDiscovery\\005_snpClass\\class\\Synonymous.txt";
        String geneFeatureFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
        String infoFileS = "E:\\Database\\InfoFile\\ChrLenCentPosi_agpV3.txt";
        String outfileS = "M:\\production\\maf\\002_recombination\\binTable.txt";
        int chrNum = 10;
        int binSize = 10000000;
        int[][] bounds = new int[chrNum][];
        int[][] crossover = new int[chrNum][];
        int[][] cdsLength = new int[chrNum][];
        double[][] del = new double[chrNum][];
        double[][] syn = new double[chrNum][];
        Table t = new Table (infoFileS);
        for (int i = 0; i < chrNum; i++) {
            int[][] bound = FArrayUtils.getSubsetsIndicesBySubsetSize(Integer.valueOf(t.content[i][1]), binSize);
            bounds[i] = new int[bound.length];
            cdsLength[i] = new int[bound.length];
            crossover[i] = new int[bound.length];
            del[i] = new double[bound.length];
            syn[i] = new double[bound.length];
            for (int j = 0; j < bound.length; j++) {
                bounds[i][j] = bound[j][0];
            }
        }
        
        t = new Table (transcriptFileS);
        ArrayList<String> tranList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.content[i][0].equals("0")) continue;
            tranList.add(t.content[i][0]);
        }
        String[] transName = tranList.toArray(new String[tranList.size()]);
        Arrays.sort(transName);
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            String name = gf.getTranscriptName(i, 0);
            if (Arrays.binarySearch(transName, name) < 0) continue;
            int chrIndex = gf.getGeneChromosome(i) - 1;
            List<Range> cds = gf.getCDSList(i, 0);
            for (int j = 0; j < cds.size(); j++) {
                int index = Arrays.binarySearch(bounds[chrIndex], cds.get(j).start);
                if (index < 0) index = -index-2;
                if (index == bounds[chrIndex].length-1) {
                    cdsLength[chrIndex][index]+=(cds.get(j).end-cds.get(j).start);
                }
                else {
                    if (cds.get(j).end <= bounds[chrIndex][index+1]) {
                        cdsLength[chrIndex][index]+=(cds.get(j).end-cds.get(j).start);
                    }
                    else {
                        cdsLength[chrIndex][index]+=(bounds[chrIndex][index+1]-cds.get(j).start);
                        cdsLength[chrIndex][index+1]+=(cds.get(j).end-bounds[chrIndex][index+1]);
                    }
                }
            }
        }
        t = new Table (crossOverFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            int chrIndex = Integer.valueOf(t.content[i][0])-1;
            int index = Arrays.binarySearch(bounds[chrIndex], Integer.valueOf(t.content[i][1]));
            if (index < 0) index = -index-2;
            crossover[chrIndex][index]++;
        }
        t = new Table (deleteriousSNPFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            int chrIndex = Integer.valueOf(t.content[i][0])-1;
            int index = Arrays.binarySearch(bounds[chrIndex], Integer.valueOf(t.content[i][1]));
            if (index < 0) index = -index-2;
            del[chrIndex][index]++;
        }
        t = new Table (synSNPFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            int chrIndex = Integer.valueOf(t.content[i][0])-1;
            int index = Arrays.binarySearch(bounds[chrIndex], Integer.valueOf(t.content[i][1]));
            if (index < 0) index = -index-2;
            syn[chrIndex][index]++;
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Chr\tBinStart\tCrossoverCount\tDelCount\tSynCount\tDelCountPerSite\tSynCountPerSite\tDelSynRatio");
            bw.newLine();
            for (int i = 0; i < chrNum; i++) {
                for (int j = 0; j < bounds[i].length; j++) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(i+1).append("\t").append(bounds[i][j]).append("\t").append(crossover[i][j]).append("\t").append(del[i][j]).append("\t").append(syn[i][j]);
                    sb.append("\t").append((double)del[i][j]/cdsLength[i][j]).append("\t").append((double)syn[i][j]/cdsLength[i][j]);
                    if (syn[i][j] == 0) {
                        sb.append("\t").append(Double.NaN);
                    }
                    else {
                        sb.append("\t").append((double)del[i][j]/syn[i][j]);
                    }
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
    
    void makeRecombinationPointTable () {
        String inDirS = "M:\\production\\maf\\002_recombination\\source\\";
        String outfileS = "M:\\production\\maf\\002_recombination\\recombinationPoint.txt";
        int chrNum = 10;
        File[] fs = new File(inDirS).listFiles();
        TIntArrayList[] posList = new TIntArrayList[chrNum];
        for (int i = 0; i < posList.length; i++) {
            posList[i] = new TIntArrayList();
        }
        try {
            for (int i = 0; i < fs.length; i++) {
                BufferedReader br = IoUtils.getTextReader(fs[i].getAbsolutePath());
                String temp = br.readLine();
                while ((temp = br.readLine()) != null) {
                    List<String> l = FStringUtils.fastSplit(temp);
                    if (l.get(3).startsWith("het")) continue;
                    int chrIndex = Integer.valueOf(l.get(0))-1;
                    int pos = (Integer.valueOf(l.get(4))+Integer.valueOf(l.get(5)))/2;
                    posList[chrIndex].add(pos);
                }
                br.close();
            }
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Chr\tPos");
            bw.newLine();
            for (int i = 0; i < posList.length; i++) {
                int[] pos = posList[i].toArray();
                Arrays.sort(pos);
                for (int j = 0; j < pos.length; j++) {
                    bw.write(String.valueOf(i+1)+"\t"+String.valueOf(pos[j]));
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
