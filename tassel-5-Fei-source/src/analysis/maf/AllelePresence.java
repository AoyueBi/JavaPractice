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
import gnu.trove.list.array.TIntArrayList;
import graphcis.r.Histogram;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class AllelePresence {
    
    public AllelePresence () {
        //this.extractBuscoGeneInMaize();
        //this.mkBuscoGeneCDSRange();
        //this.getBuscoGeneAllelePresence();
        //this.mkHistogramOfAllelePresence();
    }
    
    private void mkHistogramOfAllelePresence () {
        String infileS = "M:\\production\\maf\\annotations\\allelePresence\\buscoPlant\\buscoInMaize.cds.chr10.ap.txt";
        String outfileS = "M:\\production\\maf\\annotations\\allelePresence\\buscoPlant\\buscoInMaize.cds.chr10.ap.histo.pdf";
        Table t = new Table (infileS);
        int size = 100000;
        double[] depth = new double[size];
        for (int i = 0; i < size; i++) {
            int idx = (int)(Math.random()*t.getRowNumber());
            depth[i] = t.getDoubleValue(idx, 0);
        }
        Histogram h = new Histogram(depth);
        h.setBreakNumber(100);
        h.setTitle("BUSCO genes presence in maize HapMap3 lines");
        h.setXLab("Number of lines in HapMap3");
        h.saveGraph(outfileS);
        
    }
    
    private void getBuscoGeneAllelePresence () {
        String maizeGeneRangeFileS = "M:\\production\\maf\\annotations\\allelePresence\\buscoPlant\\buscoInMaize.cds.range.txt";
        String infileS = "M:\\production\\maf\\annotations\\maf\\mafAnnotation\\data\\chr010_maf.txt.gz";
        String buscoAPFileS = "M:\\production\\maf\\annotations\\allelePresence\\buscoPlant\\buscoInMaize.cds.chr10.ap.txt";
        Ranges rs = new Ranges(maizeGeneRangeFileS, IOFileFormat.Text);
        rs.sortByStartPosition();
        TIntArrayList apList = new TIntArrayList();
        try {
            BufferedReader br = IoUtils.getTextGzipReader(infileS);
            String temp = br.readLine();
            int cnt = -1;
            while ((temp = br.readLine()) != null) {
                cnt++;
                if(cnt%10000000 == 0) System.out.println(cnt);
                String[] tem = temp.split("\t");
                if (tem[4].startsWith("NA")) continue;
                int chr = Integer.valueOf(tem[0]);
                int pos = Integer.valueOf(tem[1]);
                int ap = Integer.valueOf(tem[4]);
                if (rs.isInRanges(chr, pos)) {
                    apList.add(ap);
                }
            }
            br.close();
            BufferedWriter bw = IoUtils.getTextWriter(buscoAPFileS);
            bw.write("AllelePresence");
            bw.newLine();
            for (int i = 0; i < apList.size(); i++) {
                bw.write(String.valueOf(apList.get(i)));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void mkBuscoGeneCDSRange () {
        String geneFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
        String maizeGeneFileS = "M:\\production\\maf\\annotations\\allelePresence\\buscoPlant\\buscoInMaize.txt";
        String maizeGeneRangeFileS = "M:\\production\\maf\\annotations\\allelePresence\\buscoPlant\\buscoInMaize.cds.range.txt";
        Table t = new Table (maizeGeneFileS);
        String[] name = new String[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) name[i] = t.content[i][0];
        GeneFeature gf = new GeneFeature(geneFileS);
        ArrayList<Range> rList = new ArrayList();
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            List<Range> cdsList = gf.getCDSList(i, 0);
            for (int j = 0; j < cdsList.size(); j++) {
                rList.add(cdsList.get(j));
            }
        }
        Ranges r = new Ranges (rList, "BuscoGeneInMaize_cds");
        r.writeFile(maizeGeneRangeFileS, IOFileFormat.Text);
    }
    
    private void extractBuscoGeneInMaize () {
        String alignmentFileS = "M:\\production\\maf\\annotations\\allelePresence\\buscoPlant\\buscoOnMaizeCDS_1e50.alignment.txt";
        String maizeGeneFileS = "M:\\production\\maf\\annotations\\allelePresence\\buscoPlant\\buscoInMaize.txt";
        try {
            BufferedReader br = IoUtils.getTextReader(alignmentFileS);
            String temp = null;
            int bcnt = 0;
            int cnt = 0;
            ArrayList<String> nameList = new ArrayList();
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("Query=")) {
                    bcnt++;
                    br.readLine();br.readLine();br.readLine();
                    temp = br.readLine();
                    if (temp.contains("No hits found")) continue;
                    String[] tem = br.readLine().split(" +");
                    nameList.add(tem[1]);
                    cnt++;
                }
            }
            br.close();
            BufferedWriter bw = IoUtils.getTextWriter(maizeGeneFileS);
            bw.write("Gene");
            bw.newLine();
            String[] name = nameList.toArray(new String[nameList.size()]);
            Arrays.sort(name);
            for (int i = 0; i < nameList.size(); i++) {
                bw.write(name[i]);
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
