/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maf;

import format.Table;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class PopGenGroup {
    
    public PopGenGroup () {
        //this.mkSubsetOfHapMap();
        //this.convertToFasta();//deprecated by using Mega7
        //this.convertCSTToCSV();
        //this.summarizeClustering();
        //this.getDepthOfHmp321();
        //this.getHighDepthTaxaList();
    }
    
    void getHighDepthTaxaList () {
        String depthInfoFileS = "M:\\production\\maf\\popgen\\highDepthGroup\\hmp321.depth.txt";
        String hmp321GeneticGroup = "M:\\production\\maf\\popgen\\group\\geneticGroup.txt";
        String outfileS = "M:\\production\\maf\\popgen\\highDepthGroup\\hmp321.highDepth.group.txt";
        double minDepth = 5;
        Table t = new Table (depthInfoFileS);
        HashMap<String, Double> taxaDepthMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.content[i][2].startsWith("null")) continue;
            double depth = Double.valueOf(t.content[i][2]);
            if (depth < minDepth) continue;
            taxaDepthMap.put(t.content[i][0], depth);
        }
        Set<String> keys = taxaDepthMap.keySet();
        t = new Table (hmp321GeneticGroup);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Taxa\tPC1\tPC2\tGroupIndex\tGeneticGroup\tDepth");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                String taxon = t.content[i][0];
                if (!keys.contains(taxon)) continue;
                StringBuilder sb = new StringBuilder();
                sb.append(taxon).append("\t").append(t.content[i][1]).append("\t").append(t.content[i][2]).append("\t").append(t.content[i][3]).append("\t").append(t.content[i][4]).append("\t").append(taxaDepthMap.get(taxon));
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
    
    void getDepthOfHmp321 () {
        String robertDepthFileS = "M:\\production\\maf\\popgen\\highDepthGroup\\source\\depthByRobert.txt";
        String siteDepthFileS = "M:\\production\\maf\\annotations\\siftScore\\009_hmp321Depth\\taxaDepth.summary.txt";
        String outfileS = "M:\\production\\maf\\popgen\\highDepthGroup\\hmp321.depth.txt";
        Table t = new Table (robertDepthFileS);
        HashMap<String, Double> robertMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            robertMap.put(t.content[i][0], Double.valueOf(t.content[i][1]));
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Taxa\tSNPSiteDepth\tSequencingDepth");
            bw.newLine();
            t = new Table (siteDepthFileS);
            for (int i = 0; i < t.getRowNumber(); i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(t.content[i][0]).append("\t").append(t.content[i][2]).append("\t").append(robertMap.get(t.content[i][0]));
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
    
    void summarizeClustering () {
        String infileS = "M:\\production\\maf\\popgen\\group\\kmeans_model_cluster.txt.arff";
        String sourceFileS = "M:\\production\\maf\\popgen\\group\\source\\source.txt";
        String teocinte = "M:\\production\\maf\\popgen\\group\\teosinte.txt";
        String outfileS = "M:\\production\\maf\\popgen\\group\\geneticGroup.txt";
        Table t = new Table(teocinte);
        String[] teo = new String[t.getRowNumber()];
        for (int i = 0; i < teo.length; i++) {
            teo[i] = t.content[i][0];
        }
        Arrays.sort(teo);
        t = new Table (sourceFileS);
        HashMap<String, String> cauMap = new HashMap();
        HashMap<String, String> ciMap = new HashMap();
        HashMap<String, String> gerMap = new HashMap();
        HashMap<String, String> hp2Map = new HashMap();
        HashMap<String, String> hp2eMap = new HashMap();
        HashMap<String, String> p2822Map = new HashMap();
        HashMap<String, String> p2824Map = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.content[i][1].isEmpty()) cauMap.put(t.content[i][0], "0");
            else cauMap.put(t.content[i][0], "1");
            if (t.content[i][2].isEmpty()) ciMap.put(t.content[i][0], "0");
            else ciMap.put(t.content[i][0], "1");
            if (t.content[i][3].isEmpty()) gerMap.put(t.content[i][0], "0");
            else gerMap.put(t.content[i][0], "1");
            if (t.content[i][4].isEmpty()) hp2Map.put(t.content[i][0], "0");
            else hp2Map.put(t.content[i][0], "1");
            if (t.content[i][5].isEmpty()) hp2eMap.put(t.content[i][0], "0");
            else hp2eMap.put(t.content[i][0], "1");
            if (t.content[i][6].isEmpty()) p2822Map.put(t.content[i][0], "0");
            else p2822Map.put(t.content[i][0], "1");
            if (t.content[i][7].isEmpty()) p2824Map.put(t.content[i][0], "0");
            else p2824Map.put(t.content[i][0], "1");
        }
        HashMap<Integer, String> groupMap = new HashMap();
        groupMap.put(0, "Teosinte");
        groupMap.put(1, "Tropical-subtropical");
        groupMap.put(2, "Non-stiff stalk");
        groupMap.put(3, "Stiff stalk");
        groupMap.put(4, "Mixed");
        groupMap.put(5, "China specific");
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            String temp = br.readLine();
            while (!((temp = br.readLine()).startsWith("@data"))) {}
            bw.write("Taxa\tPC1\tPC2\tGroupIndex\tGeneticGroup\tCAU	CIMMYT-BGI	German	HapMap2	HapMap2_extra	282-2x	282-4x");
            bw.newLine();
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split(",");
                StringBuilder sb = new StringBuilder();
                String taxon = tem[1].replaceFirst("\\.0", "");
                int groupIndex = Integer.valueOf(tem[7].replaceFirst("cluster", ""))+1;
                if (Arrays.binarySearch(teo, taxon) >= 0) groupIndex = 0;
                sb.append(taxon).append("\t").append(tem[2]).append("\t").append(tem[3]).append("\t").append(groupIndex).append("\t").append(groupMap.get(groupIndex));
                sb.append("\t").append(cauMap.get(taxon));
                sb.append("\t").append(ciMap.get(taxon));
                sb.append("\t").append(gerMap.get(taxon));
                sb.append("\t").append(hp2Map.get(taxon));
                sb.append("\t").append(hp2eMap.get(taxon));
                sb.append("\t").append(p2822Map.get(taxon));
                sb.append("\t").append(p2824Map.get(taxon));
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
    
    void convertCSTToCSV () {
        String infileS = "M:\\production\\maf\\popgen\\group\\mds_values.txt";
        String outfileS = "M:\\production\\maf\\popgen\\group\\mds_values.csv";
        Table t = new Table (infileS);
        t.writeTable(outfileS, ",");
    }
    
    void convertToFasta () {
        String infileS = "M:\\production\\maf\\popgen\\group\\fasta_table.txt";
        String outfileS = "M:\\production\\maf\\popgen\\group\\sequence.fa";
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            String temp = br.readLine();
            String[] tem = null;
            while ((temp = br.readLine()) != null) {
                tem = temp.split("\t");
                String name = tem[0];
                StringBuilder sb = new StringBuilder();
                for (int i = 1; i < tem.length; i++) {
                    if (!(tem[i].startsWith("A") || tem[i].startsWith("T") || tem[i].startsWith("G") || tem[i].startsWith("C"))) tem[i] = "N";
                    sb.append(tem[i]);
                }
                bw.write(">"+name);
                bw.newLine();
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
    
    void mkSubsetOfHapMap () {
        String infileS = "O:\\Zea\\Genotypes\\WGS\\HapMap\\v3\\v321\\unimp\\vcf\\merged_flt_c10.vcf.gz";
        String outfileS = "M:\\production\\maf\\popgen\\group\\random_chr10_hmp321.vcf.gz";
        try {
            BufferedReader br = IoUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IoUtils.getTextGzipWriter(outfileS);
            String temp = null;
            while ((temp = br.readLine()).startsWith("##")){}
            bw.write(temp);;
            bw.newLine();
            int cnt = 0;
            int cnt2 = 0;
            while ((temp = br.readLine()) != null) {
                cnt++;
                if (cnt%300 != 0) continue;
                bw.write(temp);
                bw.newLine();
                cnt2++;
                if (cnt2%100 == 0) System.out.println(String.valueOf(cnt2)+" SNPs written");
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
}
