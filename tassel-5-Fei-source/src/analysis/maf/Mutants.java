/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maf;

import format.GeneFeature;
import format.Table;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class Mutants {
    
    public Mutants () {
        this.muStatistics();
        //this.mutatorAndDeleterious();
    }
    
    void acStatistics () {
        //only 2457 insertions or deletion, no statistical power
    }
    
    void mutatorAndDeleterious () {
        String geneFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
        String muInfoFileS = "M:\\production\\maf\\mutants\\statistics\\mutatorInfo.txt";
        String geneModelS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\highConfidence_transcript.txt";
        String geneWithMutatorFileS = "M:\\production\\maf\\mutants\\geneWithMutator.txt";
        int upstream = 2000;
        int downstream = 2000;
        Table t = new Table (geneModelS);
        ArrayList<String> geneList = new ArrayList();
        HashMap<String, Double> geneDeleMap = new HashMap();
        HashMap<String, Double> geneGerpMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.content[i][4].equals("0")) continue;
            String gene = t.content[i][0];
            if (gene.startsWith("GRM")) {
                gene = gene.split("_")[0];
            }
            else {
                gene = gene.replaceFirst("_FGT", "_FG");
            }
            geneList.add(gene);
            geneDeleMap.put(gene, Double.valueOf(t.content[i][11]));
            geneGerpMap.put(gene, Double.valueOf(t.content[i][16]));
        }
        String[] gene = geneList.toArray(new String[geneList.size()]);
        Arrays.sort(gene);
        GeneFeature gf = new GeneFeature(geneFileS);
        gf.sortGeneByName();
        t = new Table (muInfoFileS);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(geneWithMutatorFileS);
            bw.write("Gene\tChr\tStart\tEnd\tDeleteriousPerSite\tMeanGerpScore\tNumOfMutator\tTotalFrequency\tMutatorFrequencys");
            bw.newLine();
            for (int i = 0; i < gene.length; i++) {
                StringBuilder sb = new StringBuilder(gene[i]);
                int index = gf.getGeneIndex(gene[i]);
                int chr = gf.getGeneChromosome(index);
                int start = gf.getTranscriptStart(index, 0);
                int end = gf.getTranscriptEnd(index, 0);
                sb.append("\t").append(chr).append("\t").append(start).append("\t").append(end).append("\t").append(geneDeleMap.get(gene[i])).append("\t").append(geneGerpMap.get(gene[i])).append("\t");
                int cnt = 0;
                ArrayList<String> valueList = new ArrayList();
                for (int j = 0; j < t.getRowNumber(); j++) {
                    if (t.content[j][1].startsWith("s")) continue;
                    if (Integer.valueOf(t.content[j][1]) != chr) continue;
                    int pos = Integer.valueOf(t.content[j][5]);
                    if (pos < start-upstream) continue;
                    if (pos > end+downstream) continue;
                    valueList.add(t.content[j][7]);
                    cnt++;
                }
                String fre = "0";
                double f = 0;
                if (cnt != 0) {
                    StringBuilder ssb = new StringBuilder();
                    for (int j = 0; j < valueList.size(); j++) {
                        ssb.append(valueList.get(j)).append(";");
                        f+=Double.valueOf(valueList.get(j));
                    }
                    ssb.deleteCharAt(ssb.length()-1);
                    fre = ssb.toString();
                }
                sb.append(cnt).append("\t").append(f).append("\t").append(fre);
                bw.write(sb.toString());
                bw.newLine();
                if (i % 500 == 0) System.out.println(i);
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    void muStatistics () {
        String infileS = "M:\\production\\maf\\004_mutants\\source\\UniformMU.gff3";
        String stockMutatorCountFileS = "M:\\production\\maf\\004_mutants\\statistics\\stockMuCount.txt";
        String muInfoFileS = "M:\\production\\maf\\004_mutants\\statistics\\mutatorInfo.txt";
        ArrayList<String> infoList = new ArrayList();
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            br.readLine();br.readLine();br.readLine();
            String temp = null;
            while ((temp = br.readLine()) != null) {
                infoList.add(temp);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        String[] infoArray = infoList.toArray(new String[infoList.size()]);
        HashSet<String> stockSet = new HashSet();
        for (int i = 0; i < infoArray.length; i++) {
            String[] temp = infoArray[i].split("Stock=")[1].split(",");
            for (int j = 0; j < temp.length; j++) {
                stockSet.add(temp[j]);
            }
        }
        String[] stocks = stockSet.toArray(new String[stockSet.size()]);
        Arrays.sort(stocks);
        int[] stockMuCount = new int[stocks.length];
        int totalMuCount = 0;
        for (int i = 0; i < infoArray.length; i++) {
            String[] temp = infoArray[i].split("Stock=")[1].split(",");
            for (int j = 0; j < temp.length; j++) {
                int index = Arrays.binarySearch(stocks, temp[j]);
                stockMuCount[index]++;
                totalMuCount++;
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(stockMutatorCountFileS);
            bw.write("Stock\tNumOfMutator\tProportionOfMutator");
            bw.newLine();
            for (int i = 0; i < stocks.length; i++) {
                StringBuilder sb = new StringBuilder(stocks[i]);
                sb.append("\t").append(stockMuCount[i]).append("\t").append((double)stockMuCount[i]/totalMuCount);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            
            bw = IoUtils.getTextWriter(muInfoFileS);
            bw.write("ID\tChr\tStart\tEnd\tLength\tPosition\tNumOfStocks\tFrequency");
            bw.newLine();
            for (int i = 0; i < infoArray.length; i++) {
                StringBuilder sb = new StringBuilder();
                String[] temp = infoArray[i].split("\t");
                String[] tem = temp[8].split(";");
                sb.append(tem[1].replaceFirst("ID=", "")).append("\t").append(temp[0].replaceFirst("Chr", "")).append("\t");
                int start = Integer.valueOf(temp[3]);
                int end = Integer.valueOf(temp[4]);
                int length = end-start+1;
                int position = (start+end)/2;
                int numOfStocks = tem[2].split(",").length;
                sb.append(start).append("\t").append(end).append("\t").append(length).append("\t").append(position).append("\t").append(numOfStocks).append("\t").append((double)numOfStocks/stocks.length);
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
