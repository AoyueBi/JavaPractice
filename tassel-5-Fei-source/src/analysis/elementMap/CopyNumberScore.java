/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.elementMap;

import format.Fasta;
import format.Table;
import gnu.trove.list.array.TDoubleArrayList;
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
import utils.Benchmark;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class CopyNumberScore {
    
    public CopyNumberScore () {
        
    }
    
    public void writeCombinedScore (String taxaDepthFileS, String taxaScoreDirS, String referenceGenomeFileS, String combinedScoreDirS, double minDepth) {
        new File(combinedScoreDirS).mkdir();
        Table t = new Table (taxaDepthFileS);
        ArrayList<String> taxaList = new ArrayList();
        TDoubleArrayList depthList = new TDoubleArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            double d = Double.valueOf(t.content[i][1]);
            if (d < minDepth) continue;
            taxaList.add(t.content[i][0]);
            depthList.add(Double.valueOf(t.content[i][1]));
        }
        String[] taxa = taxaList.toArray(new String[taxaList.size()]);
        double[] depth = depthList.toArray();
        System.out.println(String.valueOf(taxa.length)+" taxa are selected based on min depth of "+ String.valueOf(minDepth));
        Fasta f = new Fasta(referenceGenomeFileS);
        f.sortRecordByNameValue();
        ArrayList<Integer> chrIndexList = new ArrayList();
        for (int i = 0; i < f.getSeqNumber(); i++) chrIndexList.add(i);
        System.out.println(String.valueOf(f.getSeqNumber())+" chromosomes. Start combining scores...");
        chrIndexList.parallelStream().forEach(chrIndex -> {
            int chr = Integer.valueOf(f.getName(chrIndex));
            int chrLength = f.getSeqLength(chrIndex);
            String[] path = new String[taxa.length];
            for (int i = 0; i < path.length; i++) {
                StringBuilder sb  = new StringBuilder();
                sb.append(taxa[i]).append("_chr").append(FStringUtils.getNDigitNumber(3, chr)).append("_CpScore.txt.gz");
                path[i] = new File(new File(taxaScoreDirS, taxa[i]).getAbsolutePath(), sb.toString()).getAbsolutePath();
            }
            StringBuilder sb  = new StringBuilder();
            sb.append("chr").append(FStringUtils.getNDigitNumber(3, chr)).append("_combinedScore.txt.gz");
            String outfileS = new File(combinedScoreDirS, sb.toString()).getAbsolutePath();
            BufferedReader[] brs = new BufferedReader[path.length];
            try{
                BufferedWriter bw = IoUtils.getTextGzipWriter(outfileS);
                bw.write("Chr\tPos\tCpScore\tCpScoreSD");
                bw.newLine();
                String[] temps = new String[brs.length];
                for (int i = 0; i < brs.length; i++) {
                    brs[i] = IoUtils.getTextGzipReader(path[i], 1024);
                    temps[i] = brs[i].readLine();
                }
                String nullS = "\tNA\tNA\tNA";
                int cnt = 0;
                for (int i = 0; i < chrLength; i++) {
                    sb = new StringBuilder();
                    for (int j = 0; j < brs.length; j++) {
                        temps[j] = brs[j].readLine();
                    }
                    if (temps[0].startsWith("N")) {
                        sb.append(chr).append("\t").append(i+1).append(nullS);
                    }
                    else {
                        double[] cs = new double[brs.length];
                        for (int j = 0; j < cs.length; j++) {
                            cs[j] = Double.valueOf(temps[j])/depth[j];
                        }
                        DescriptiveStatistics sta = new DescriptiveStatistics(cs);
                        double cp = sta.getMean();
                        double cpsd = sta.getStandardDeviation();
                        sb.append(chr).append("\t").append(i+1).append("\t").append((float)cp).append("\t").append((float)cpsd);
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                    cnt++;
                    if (cnt%10000000 == 0) System.out.println(String.valueOf(cnt)+" sites processed on chromosome " +String.valueOf(chr));
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        System.out.println("Combining score is finished. Stored at "+combinedScoreDirS);
    }
    
    public void calibrateTaxonDepth (String singleCopyFileS, String taxaScoreDirS, String taxaDepthDirS) {
        double maxDepth = 200;
        double interval = 0.1;
        int binNum = (int)(maxDepth/interval);
        File taxaDepthDir = new File(taxaDepthDirS);
        taxaDepthDir.mkdir();
        File depthDir = new File(taxaDepthDir, "depth");
        depthDir.mkdir();
//        File pdfDir = new File(taxaDepthDir, "pdf");
//        pdfDir.mkdir();
        String taxaDepthFileS = new File(taxaDepthDir, "taxaDepth.txt").getAbsolutePath();
        HashSet<Integer> siteSet = new HashSet();
        Table t = new Table(singleCopyFileS);
        int chr = Integer.valueOf(t.content[0][0]);
        for (int i = 0; i < t.getRowNumber(); i++) siteSet.add(Integer.valueOf(t.content[i][1]));
        File[] dir = new File(taxaScoreDirS).listFiles();
        String[] taxaName = new String[dir.length];
        double[] taxaDepth = new double[dir.length];
        for (int i = 0; i < taxaName.length; i++) {
            taxaName[i] = dir[i].getName();
        }
        Arrays.sort(taxaName);
        List<File> dirList = Arrays.asList(dir);
        dirList.parallelStream().forEach(f -> {
            String taxonName = f.getName();
            String depthFileS = new File (depthDir, taxonName+"_depth.txt").getAbsolutePath();
//            String pdfFileS = new File (pdfDir, taxonName+"_depth.txt").getAbsolutePath();
            StringBuilder sb  = new StringBuilder();
            sb.append(taxonName).append("_chr").append(FStringUtils.getNDigitNumber(3, chr)).append("_CpScore.txt.gz");
            String scoreFileS = new File(f, sb.toString()).getAbsolutePath();
            double[] siteDepth = new double[siteSet.size()];
            try {
                BufferedReader br = IoUtils.getTextGzipReader(scoreFileS);
                String temp = br.readLine();
                int cnt = 0;
                int counter = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%50000000==0) System.out.println(String.valueOf(cnt)+" sites processed from " + scoreFileS);
                    if (!siteSet.contains(cnt)) continue;
                    siteDepth[counter] = Double.valueOf(temp);
                    counter++;
                }
                br.close();
                double[] bins = new double[binNum];
                int[] binCount = new int[binNum];
                for (int i = 0; i < bins.length; i++) {
                    bins[i] = interval*i;
                }
                BufferedWriter bw = IoUtils.getTextWriter(depthFileS);
                bw.write("CpScore_SingleCopySite");
                bw.newLine();
                for (int i = 0; i < siteDepth.length; i++) {
                    bw.write(String.valueOf(siteDepth[i]));
                    bw.newLine();
                    int index = Arrays.binarySearch(bins, siteDepth[i]);
                    if (index < 0) index = -index-2;
                    if (index < bins.length-1) {
                        binCount[index]++;
                    }  
                }
                bw.flush();
                bw.close();
                int maxIndex = 0;
                int maxValue = 0;
                for (int i = 0; i < binCount.length; i++) {
                    if (binCount[i] > maxValue) {
                        maxValue = binCount[i];
                        maxIndex = i;
                    }
                }
                int index = Arrays.binarySearch(taxaName, taxonName);
                taxaDepth[index] = bins[maxIndex];
                System.out.println("Depth of "+ taxonName +" is "+String.valueOf(taxaDepth[index]));
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        try {
            BufferedWriter bw = IoUtils.getTextWriter(taxaDepthFileS);
            bw.write("Taxa\tDepth");
            bw.newLine();
            for (int i = 0; i < taxaName.length; i++) {
                bw.write(taxaName[i]+"\t"+String.valueOf(taxaDepth[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void writeCpScoreOfTaxa (String taxaFqMapFileS, String fastqDirS, String kmerListDirS, String referenceGenomeFileS, String taxaScoreDirS, int concurrencyLevel) {
        System.out.println("Calculating CpScore for all taxa...");
        long start = 0;
        HashMap<String, String[]> taxaFqPathMap = this.getTaxaFastqPathMap(taxaFqMapFileS, fastqDirS);
        Set<String> taxaSet = taxaFqPathMap.keySet();
//        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", String.valueOf(concurrencyLevel));
//        taxaSet.parallelStream().forEach(taxon -> {
//            String[] paths = taxaFqPathMap.get(taxon);
//            List<String> pathList = Arrays.asList(paths);
//            System.out.println("Start calculating CpScore for "+taxon);
//            new TaxonCpScore(referenceGenomeFileS, kmerListDirS, taxon, pathList, taxaScoreDirS);
//        });
        String[] taxa = taxaSet.toArray(new String[taxaSet.size()]);
        Arrays.sort(taxa);
        int cnt = 0;
        for (int i = 0; i < taxa.length; i+=concurrencyLevel) {
            System.out.println("Batch ID: "+ String.valueOf(cnt+1)+". ("+String.valueOf(concurrencyLevel)+" taxa / batch)");
            int endIndex = i + concurrencyLevel;
            if (endIndex >=taxa.length) endIndex = taxa.length;
            ArrayList<String> subTaxaList = new ArrayList();
            for (int j = i; j < endIndex; j++) {
                subTaxaList.add(taxa[j]);
            }
            subTaxaList.parallelStream().forEach(taxon -> {
                String[] paths = taxaFqPathMap.get(taxon);
                List<String> pathList = Arrays.asList(paths);
                System.out.println("Start calculating CpScore for "+taxon);
                TaxonCpScore tcs = new TaxonCpScore(referenceGenomeFileS, kmerListDirS, taxon, pathList, taxaScoreDirS);
                tcs = null;
            });
            cnt++;
            System.gc();
        }
        StringBuilder sb = new StringBuilder();
        sb.append("Calculating copy number score of all taxa is finished. Time span: ").append(Benchmark.getTimeSpanHours(start)).append(" hours. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(sb.toString());
    }
    
    private HashMap<String, String[]> getTaxaFastqPathMap (String taxaFqMapFileS, String fastqDirS) {
        HashSet<String> actualFqSet = new HashSet();
        File[] fs = new File(fastqDirS).listFiles();
        for (int i = 0; i < fs.length; i++) actualFqSet.add(fs[i].getName());
        String[] actualFq = actualFqSet.toArray(new String[actualFqSet.size()]);
        Arrays.sort(actualFq);
        HashMap<String, String[]> taxaFqPathMap = new HashMap(); 
        int fqCnt = 0;
        int cnt = 0;
        int taxaCnt = 0;
        try {
            BufferedReader br = IoUtils.getTextReader(taxaFqMapFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                taxaCnt++;
                List<String> tempList = FStringUtils.fastSplit(temp);
                ArrayList<String> pathList = new ArrayList();
                for (int i = 1; i < tempList.size(); i++) {
                    String q = tempList.get(i);
                    fqCnt++;
                    int index = Arrays.binarySearch(actualFq, q);
                    if (index < 0) continue;
                    String path = new File(fastqDirS, q).getAbsolutePath();
                    pathList.add(path);
                    cnt++;
                }
                if (pathList.isEmpty()) continue;
                String[] path = pathList.toArray(new String[pathList.size()]);
                Arrays.sort(path);
                taxaFqPathMap.put(tempList.get(0), path);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Number of taxa in taxaFqMapFile: "+String.valueOf(taxaCnt)+". Number of existing taxa: " +String.valueOf(taxaFqPathMap.size()));
        System.out.println("Number of fastq in taxaFqMapFile: "+String.valueOf(fqCnt)+". Number of existing fastq: " +String.valueOf(cnt));
        return taxaFqPathMap;
    }
    
    public void writeReferenceKmers (int kmerLength, String referenceGenomeFileS, String kmerListDirS) {
        new ReferenceKmer (kmerLength,referenceGenomeFileS,kmerListDirS);
    }
}
