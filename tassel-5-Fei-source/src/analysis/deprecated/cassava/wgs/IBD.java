/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.cassava.wgs;

import format.Table;
import gnu.trove.list.array.TIntArrayList;
import graphcis.r.CumulativeDistribution;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import net.maizegenetics.analysis.distance.IBSDistanceMatrix;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import utils.FArrayUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class IBD {
    
    public IBD () {
        this.findIBDRegions();
        //this.getIBSDistribution();
        //this.mkDistribution();
    }
    
    public void mkDistribution () {
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\IBD\\ibs_cassava_list.txt";
        Table t = new Table (infileS);
        int size = 20000;
        double[] v = new double[size];
        for (int i = 0; i < size; i++) {
            v[i] = Double.valueOf(t.content[(int)(Math.random()*t.getRowNumber())][5]);
        }
        CumulativeDistribution d =new CumulativeDistribution(v);
        d.showGraph();
    }
    
    public void getIBSDistribution () {
        double missingCut = 0.3;
        double mafCut = 0.05;
        double maxHetCut = 0.2;
        int regionNumber = 20;
        double IBSCut = 0.02;
        String infileDirS = "Q:\\D\\Cassava\\GBSdataForHMPfilter\\";
        String taxaFileS = "M:\\pipelineTest\\cassava\\wgs\\IBD\\source\\WGSsamples.txt";
        String ibsFileS = "M:\\pipelineTest\\cassava\\wgs\\IBD\\ibs_cassava_list.txt";
        File[] fs = new File(infileDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        ArrayList<String>[] IBDLists = new ArrayList[fs.length];
        fList.parallelStream().forEach(f -> {
            String[] temp = f.getName().split("_");
            int chr = Integer.valueOf(temp[4].replaceFirst(".h5", "").replaceFirst("chr", ""));
            IBDLists[chr-1] = this.getIBSList(f.getAbsolutePath(), taxaFileS, missingCut, regionNumber);
        });
        try {
            BufferedWriter bw = IoUtils.getTextWriter(ibsFileS);
            bw.write("Taxa_1\tTaxa_2\tChr\tStart\tEnd\tIBSValue");
            bw.newLine();
            for (int i = 0; i < IBDLists.length; i++) {
                for (int j = 0; j < IBDLists[i].size(); j++) {
                    bw.write(IBDLists[i].get(j));
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
    
    public ArrayList<String> getIBSList (String infileS, String taxaFileS, double missingCut, int regionNumber) {
        GenotypeTable gt = ImportUtils.read(infileS);
        Table t = new Table (taxaFileS);
        String[] hapMapTaxa= new String[t.getRowNumber()];
        for (int i = 0; i < hapMapTaxa.length; i++) hapMapTaxa[i] = t.content[i][0];
        TaxaList tList = this.getSharedTaxa(gt, hapMapTaxa);
        GenotypeTable tgt = FilterGenotypeTable.getInstance(gt, tList);
        TIntArrayList siteList = new TIntArrayList();
        for (int i = 0; i < tgt.numberOfSites(); i++) {
            if ((double) tgt.totalNonMissingForSite(i)/tgt.numberOfTaxa() < missingCut) continue;
            //if (tgt.minorAlleleFrequency(i) < mafCut) continue;
            double heteo = 0;
            for (int j = 0; j < tgt.numberOfTaxa(); j++) {
                if (tgt.isHeterozygous(j, i)) heteo+=1;
            }
            //if (heteo/tgt.totalNonMissingForSite(i) > maxHetCut) continue;
            siteList.add(i);
        }
        //FilterGenotypeTable sgt = FilterGenotypeTable.getInstance(tgt, siteList.toArray());
        FilterGenotypeTable sgt = null;
        //int[][] bound = FArrayUtils.getSubsetsIndicesBySubsetNumber(sgt.numberOfSites(), regionNumber);
        int[][] bound = null;
        ArrayList<String> IBSValueList = new ArrayList();
        StringBuilder sb;
        for (int i = 0; i < bound.length; i++) {
            //GenotypeTable rgt = FilterGenotypeTable.getInstance(sgt, bound[i][0], bound[i][1]-1);
            GenotypeTable rgt = null;
            double[][] distances=IBSDistanceMatrix.getInstance(rgt).getDistances();
            for (int j = 0; j < distances.length-1; j++) {
                for (int k = j+1; k < distances[j].length; k++) {
                    //if (distances[j][k] > IBSCut) continue;
                    sb = new StringBuilder();
                    sb.append(rgt.taxaName(j)).append("\t").append(rgt.taxaName(k)).append("\t");
                    sb.append(rgt.chromosome(0).getName()).append("\t").append(rgt.physicalPositions()[0]).append("\t").append(rgt.positions().get(bound[i][1]-bound[i][0]-1).getPosition()).append("\t");
                    sb.append(distances[j][k]);
                    IBSValueList.add(sb.toString());
                    
                }
            }
        }
        return IBSValueList;
    }
    
    public void findIBDRegions () {
        double missingCut = 0.3;
        double mafCut = 0.05;
        double maxHetCut = 0.2;
        int regionNumber = 5;
        double IBSCut = 0.02;
        String infileDirS = "Q:\\D\\Cassava\\GBSdataForHMPfilter\\";
        String outfileDirS = "M:\\pipelineTest\\cassava\\wgs\\IBD\\Genotype\\";
        String taxaFileS = "M:\\pipelineTest\\cassava\\wgs\\IBD\\source\\WGSsamples.txt";
        String ibdFileS = "M:\\pipelineTest\\cassava\\wgs\\IBD\\ibd_cassava_raw.txt";
        File[] fs = new File(infileDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        ArrayList<String>[] IBDLists = new ArrayList[fs.length];
        fList.parallelStream().forEach(f -> {
            String[] temp = f.getName().split("_");
            int chr = Integer.valueOf(temp[4].replaceFirst(".h5", "").replaceFirst("chr", ""));
            String outfileS = new File(outfileDirS, temp[4].replaceAll(".h5", ".hmp.txt")).getAbsolutePath();
            IBDLists[chr-1] = this.findIBDByChr(f.getAbsolutePath(), outfileS, taxaFileS, missingCut, mafCut, maxHetCut, IBSCut, regionNumber);
        });
        try {
            BufferedWriter bw = IoUtils.getTextWriter(ibdFileS);
            bw.write("Taxa_1\tTaxa_2\tChr\tStart\tEnd\tIBSValue");
            bw.newLine();
            for (int i = 0; i < IBDLists.length; i++) {
                for (int j = 0; j < IBDLists[i].size(); j++) {
                    bw.write(IBDLists[i].get(j));
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
    
    public ArrayList<String> findIBDByChr (String infileS, String outGenotypeFileS, String taxaFileS, double missingCut, double mafCut, double maxHetCut, double IBSCut, int regionNumber) {
        GenotypeTable gt = ImportUtils.read(infileS);
        Table t = new Table (taxaFileS);
        String[] hapMapTaxa= new String[t.getRowNumber()];
        for (int i = 0; i < hapMapTaxa.length; i++) hapMapTaxa[i] = t.content[i][0];
        TaxaList tList = this.getSharedTaxa(gt, hapMapTaxa);
        GenotypeTable tgt = FilterGenotypeTable.getInstance(gt, tList);
        TIntArrayList siteList = new TIntArrayList();
        for (int i = 0; i < tgt.numberOfSites(); i++) {
            if ((double) tgt.totalNonMissingForSite(i)/tgt.numberOfTaxa() < missingCut) continue;
            if (tgt.minorAlleleFrequency(i) < mafCut) continue;
            double heteo = 0;
            for (int j = 0; j < tgt.numberOfTaxa(); j++) {
                if (tgt.isHeterozygous(j, i)) heteo+=1;
            }
            if (heteo/tgt.totalNonMissingForSite(i) > maxHetCut) continue;
            siteList.add(i);
        }
        //FilterGenotypeTable sgt = FilterGenotypeTable.getInstance(tgt, siteList.toArray());
        FilterGenotypeTable sgt = null;
        //int[][] bound = FArrayUtils.getSubsetsIndicesBySubsetNumber(sgt.numberOfSites(), regionNumber);
        int[][] bound = null;
        ArrayList<String> IBSValueList = new ArrayList();
        StringBuilder sb;
        for (int i = 0; i < bound.length; i++) {
            //GenotypeTable rgt = FilterGenotypeTable.getInstance(sgt, bound[i][0], bound[i][1]-1);
            GenotypeTable rgt = null;
            double[][] distances=IBSDistanceMatrix.getInstance(rgt).getDistances();
            for (int j = 0; j < distances.length-1; j++) {
                for (int k = j+1; k < distances[j].length; k++) {
                    if (distances[j][k] > IBSCut) continue;
                    sb = new StringBuilder();
                    sb.append(rgt.taxaName(j)).append("\t").append(rgt.taxaName(k)).append("\t");
                    sb.append(rgt.chromosome(0).getName()).append("\t").append(rgt.physicalPositions()[0]).append("\t").append(rgt.positions().get(bound[i][1]-bound[i][0]-1).getPosition()).append("\t");
                    sb.append(distances[j][k]);
                    IBSValueList.add(sb.toString());
                    
                }
            }
        }
        //ExportUtils.writeToHapmap(sgt, outGenotypeFileS);
        return IBSValueList;
    }
    
    public TaxaList getSharedTaxa (GenotypeTable gt,  String[] hapMapTaxa) {
        boolean[] ifExist = new boolean[hapMapTaxa.length];
        Arrays.sort(hapMapTaxa);
        ArrayList<String> l = new ArrayList();
        TaxaListBuilder tlb=new TaxaListBuilder();
        ArrayList<String> a = new ArrayList();
        for (int i = 0; i < gt.numberOfTaxa(); i++) {
            //System.out.println(gt.taxaName(i));
            int index = Arrays.binarySearch(hapMapTaxa, gt.taxaName(i).split(":")[0]);
            if (index < 0) continue;
            if (ifExist[index]) continue;
            a.add(gt.taxaName(i));
            ifExist[index] = true;
        }
        String[] sharedTaxa = a.toArray(new String[a.size()]);
        Arrays.sort(sharedTaxa);
        for (int i = 0; i < sharedTaxa.length; i++) {
            Taxon t = new Taxon (sharedTaxa[i]);
            tlb.add(t);
        }
        TaxaList tList = tlb.build();
        return tList;
    }
    
}
