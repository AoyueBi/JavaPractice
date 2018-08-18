/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maf;

import format.GeneFeature;
import format.Range;
import format.RangeAttribute;
import format.Ranges;
import format.Table;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TCharArrayList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import graphcis.r.BarPlot;
import graphcis.r.DensityPlotMultiClass;
import graphcis.r.ScatterPlot;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import utils.FArrayUtils;
import utils.FStringUtils;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class SiftScore {
    
    public SiftScore () {
        //this.extractHmp321Info();
        //this.updateHmp321InfoWithAncestralAllele();
        //this.summarizeTranscript();
        //this.combineHmpSiftGerpUScore();
        //this.mkHmp32MafPlot();
        //this.siftAndGerp();
        //this.classifySNPs3();
        this.mkBarplotOfSNPs();
        //this.deleteriousCountOnChr();
        
        //this.countDeleteriousHmp321();
        //this.mkDepthOfHmp321();
        //this.mkDepthSummary();
        //this.mkTaxaGroupFile();
        //this.countDeleteriousHmp32HighDepth();
        
        
        //this.getGeneticDistanceToB73();
        //this.mergeDeleteriousTraitDisctance();
    }
    
    private void summarizeTranscript () {
        String geneFeatureFileS = "/workdir/mingh/Zea_mays.AGPv3.26.gf.txt";
        String hmpInfoDirS = "/local/workdir/mingh/hmp321Info/";
        String gerpDirS = "/local/workdir/mingh/gerp/";
        String siftDirS = "/local/workdir/mingh/sift/";
        String outfileS = "/local/workdir/mingh/transcriptSummary.txt";
        double gerpCut = 2;
        File[] fs = new File(hmpInfoDirS).listFiles();
        int chrNum = fs.length;
        HashMap<Integer, String>[] posGeneMap = new HashMap[chrNum];
        int[][] snpPos = new int[chrNum][];
        byte[][] snps = new byte[chrNum][];
        for (int i = 0; i < chrNum; i++) {
            posGeneMap[i] = new HashMap();
        }
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        HashMap<String, Integer> geneCDSLengthMap = new HashMap();
        
        String[] genes = new String[gf.getGeneNumber()];
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            String geneName = gf.getTranscriptName(i, 0);
            genes[i] = geneName;
            List<Range> cdsList = gf.getCDSList(i, 0);
            int cnt = 0;
            int chrIndex = gf.getGeneChromosome(i)-1;
            for (int j = 0; j < cdsList.size(); j++) {
                int rStart = cdsList.get(j).start;
                int rEnd = cdsList.get(j).end;
                for (int k = rStart; k < rEnd; k++) {
                    posGeneMap[chrIndex].put(k, geneName);
                    cnt++;
                }
            }
            geneCDSLengthMap.put(geneName, cnt);
        }
        Arrays.sort(genes);
        int[] snpCount = new int[genes.length];
        List<File> hmpList = Arrays.asList(fs);
        hmpList.parallelStream().forEach(f -> {
            int chrIndex = Integer.valueOf(f.getName().split("_chr")[1].replaceFirst(".txt.gz", ""))-1;
            TIntArrayList snpPosList = new TIntArrayList();
            TByteArrayList snpList = new TByteArrayList();
            try {
                BufferedReader br = IoUtils.getTextGzipReader(f.getAbsolutePath());
                String temp = br.readLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println("Hmp\tchr"+String.valueOf(chrIndex+1)+"\t"+String.valueOf(cnt));
                    List<String> l = FStringUtils.fastSplit(temp);
                    if (l.get(3).contains("<") || l.get(3).contains(",")) continue;
                    int pos = Integer.valueOf(l.get(1));
                    String hit = posGeneMap[chrIndex].get(pos);
                    if (hit == null) continue;
                    int index = Arrays.binarySearch(genes, hit);
                    snpCount[index]++;
                    snpPosList.add(pos);
                    snpList.add(l.get(3).getBytes()[0]);
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            snpPos[chrIndex] = snpPosList.toArray();
            snps[chrIndex] = snpList.toArray();
        });
        int[] synCount = new int[genes.length];
        int[] nonCount = new int[genes.length];
        int[] delCount = new int[genes.length];
        int[] delHGCount = new int[genes.length];
        int[] naCount = new int[genes.length];
        TIntArrayList[] delPosList = new TIntArrayList[chrNum];
        fs = new File (siftDirS).listFiles();
        List<File> siftList = Arrays.asList(fs);
        siftList.parallelStream().forEach(f -> {
            int chrIndex = Integer.valueOf(f.getName().replaceFirst(".gz", ""))-1;
            delPosList[chrIndex] = new TIntArrayList();
            try {
                BufferedReader br = IoUtils.getTextGzipReader(f.getAbsolutePath());
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println("Sift\tchr"+String.valueOf(chrIndex+1)+"\t"+String.valueOf(cnt));
                    List<String> l = FStringUtils.fastSplit(temp, "\t");
                    if (!l.get(6).startsWith("CDS")) continue;
                    if (l.get(3).startsWith("GRM") && (!l.get(3).endsWith("01"))) continue;
                    String gene = l.get(3);
                    int geneIndex = Arrays.binarySearch(genes, gene);
                    if (geneIndex < 0) continue;
                    int pos = Integer.valueOf(l.get(0));
                    int index = Arrays.binarySearch(snpPos[chrIndex], pos);
                    if (index < 0) continue;
                    if (snps[chrIndex][index] != l.get(2).getBytes()[0]) continue;
                    
                    String type = null;
                    if (l.get(7).equals("*")) {
                        if (l.get(8).equals("*")) {
                            type = "StopStop";
                        }
                        else {
                            type = "StopLoss";
                        }
                    }
                    else {
                        if (l.get(8).equals("*")) {
                            type = "StopGain";
                        }
                        else {
                            if (l.get(7).equals(l.get(8))) {
                                type = "Syn";
                                synCount[geneIndex]++;
                            }
                            else {
                                type = "Non";
                                nonCount[geneIndex]++;
                                if (l.get(10).startsWith("N")) {
                                    naCount[geneIndex]++;
                                }
                                else {
                                    if (Double.valueOf(l.get(10)) < 0.05) {
                                        delCount[geneIndex]++;
                                        delPosList[chrIndex].add(pos);
                                    }
                                }
                            }
                        }
                    }
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        int[][] delPos = new int[chrNum][];
        for (int i = 0; i < chrNum; i++) {
            delPos[i] = delPosList[i].toArray();
            Arrays.sort(delPos[i]);
        }
        int[] gerpAlignCount = new int[genes.length];
        int[] snpGerpAlignCount = new int[genes.length];
        double[] gerpTree = new double[genes.length];
        double[] gerpScore = new double[genes.length];
        double[] snpGerpTree = new double[genes.length];
        double[] snpGerpScore = new double[genes.length];
        fs = new File(gerpDirS).listFiles();
        List<File> gerpList = Arrays.asList(fs);
        gerpList.parallelStream().forEach(f -> {
            int chrIndex = Integer.valueOf(f.getName().replaceFirst("roast.chrom.", "").replaceFirst(".msa.in.rates.full", ""))-1;
            try {
                BufferedReader br = IoUtils.getTextReader(f.getAbsolutePath());
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println("Gerp\tchr"+String.valueOf(chrIndex+1)+"\t"+String.valueOf(cnt));
                    String gene = posGeneMap[chrIndex].get(cnt);
                    if (gene == null) continue;
                    int geneIndex = Arrays.binarySearch(genes, gene);
                    List<String> l = FStringUtils.fastSplit(temp);
                    double treeValue = Double.valueOf(l.get(0));
                    double scoreValue = Double.valueOf(l.get(1));
                    if (treeValue == 0) continue;
                    gerpAlignCount[geneIndex]++;
                    gerpTree[geneIndex]+=treeValue;
                    gerpScore[geneIndex]+=scoreValue;
                    int index = Arrays.binarySearch(snpPos[chrIndex], cnt);
                    if (index < 0) continue;
                    snpGerpAlignCount[geneIndex]++;
                    snpGerpTree[geneIndex]+=treeValue;
                    snpGerpScore[geneIndex]+=scoreValue;
                    index = Arrays.binarySearch(delPos[chrIndex], cnt);
                    if (index < 0) continue;
                    if (scoreValue <= gerpCut) continue;
                    delHGCount[geneIndex]++;
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        
        try {
            String header = "Transcript\tCDSLength\tSNPNumber\tSNPPercentage\tIfSiftAligned\tNumberOfSyn\tPercentageSyn\tNumberOfNon\tPercentageNon\tNonVsSynRatio\tNumberOfDeleterious\tPercentageDeleterious\tNumberOfHGDeleterious\tPercentageHGDeleterious\tIfGerpAligned\tGerpAlignedCount\tPercentageGerpAlignedCount\tMeanGerpTreeLength\tMeanGerpScore\tSNPMeanGerpTreeLength\tSNPMeanGerpScore";
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write(header);
            bw.newLine();
            for (int i = 0; i < genes.length; i++) {
                StringBuilder sb =new StringBuilder(genes[i]);
                int cdsLength = geneCDSLengthMap.get(genes[i]);
                sb.append("\t").append(cdsLength).append("\t").append(snpCount[i]).append("\t").append((double)snpCount[i]/cdsLength).append("\t");
                int ifSiftAligned = 1;
                if (naCount[i] == nonCount[i]) ifSiftAligned = 0;
                sb.append(ifSiftAligned).append("\t").append(synCount[i]).append("\t").append((double)synCount[i]/cdsLength).append("\t");
                sb.append(nonCount[i]).append("\t").append((double)nonCount[i]/cdsLength).append("\t");
                double ratio = 0;
                if (synCount[i] == 0) ratio = Double.NaN;
                else ratio = (double)nonCount[i]/synCount[i];
                sb.append(ratio).append("\t").append(delCount[i]).append("\t").append((double)delCount[i]/cdsLength).append("\t");
                sb.append(delHGCount[i]).append("\t").append((double)delHGCount[i]/cdsLength).append("\t");
                int ifGerpAligned = 1;
                if (gerpAlignCount[i] == 0) ifGerpAligned = 0;
                sb.append(ifGerpAligned).append("\t").append(gerpAlignCount[i]).append("\t").append((double)gerpAlignCount[i]/cdsLength).append("\t");
                if (gerpAlignCount[i] == 0) sb.append(Double.NaN).append("\t").append(Double.NaN).append("\t");
                else sb.append((double)gerpTree[i]/gerpAlignCount[i]).append("\t").append((double)gerpScore[i]/gerpAlignCount[i]).append("\t");
                if (snpGerpAlignCount[i] == 0) sb.append(Double.NaN).append("\t").append(Double.NaN);
                else sb.append((double)snpGerpTree[i]/snpGerpAlignCount[i]).append("\t").append((double)snpGerpScore[i]/snpGerpAlignCount[i]);
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
    
    public void mergeDeleteriousTraitDisctance () {
        String distanceFileS = "M:\\production\\maf\\annotations\\siftScore\\hmp32TaxaGroup\\distanceToB73.txt";
        String deleteriousFileS = "M:\\production\\maf\\annotations\\siftScore\\hmp32DeleCount\\reccesiveDeleterious_hmp32_highDepth.txt";
        String traitFileS = "M:\\production\\maf\\annotations\\siftScore\\hmp32TaxaGroup\\Fei_load.txt";
        String mergedFileS = "M:\\production\\maf\\annotations\\siftScore\\hmp32DeleCount\\reccesiveDeleterious_merge.txt";
        HashMap<String, Double> taxaHeightMap = new HashMap ();
        HashMap<String, Double> taxaDTAMap = new HashMap();
        Table t = new Table (traitFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaHeightMap.put(t.content[i][0], Double.NaN);
            taxaDTAMap.put(t.content[i][0], Double.NaN);
            if (!t.content[i][4].equals("")) taxaHeightMap.put(t.content[i][0], Double.valueOf(t.content[i][4]));
            if (!t.content[i][5].equals("")) taxaDTAMap.put(t.content[i][0], Double.valueOf(t.content[i][5]));
        }
        t = new Table (distanceFileS);
        HashMap<String, Double> taxaDisMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaDisMap.put(t.content[i][1], Double.valueOf(t.content[i][2]));
        }
        try {
            BufferedReader br = IoUtils.getTextReader(deleteriousFileS);
            BufferedWriter bw = IoUtils.getTextWriter(mergedFileS);
            String temp = br.readLine();
            bw.write(temp+"\tDistance\tPH_BLUP_peiffer\tDTA_BLUP_peiffer");
            bw.newLine();
            while ((temp = br.readLine()) != null) {
                StringBuilder sb = new StringBuilder();
                String taxon = temp.split("\t")[0];
                sb.append(temp).append("\t").append(taxaDisMap.get(taxon)).append("\t").append(taxaHeightMap.get(taxon)).append("\t").append(taxaDTAMap.get(taxon));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void getGeneticDistanceToB73 () {
        String infileS = "Q:\\Zea\\Genotypes\\HapMap\\v3\\LDKNNi\\imp_VCF_hmp321\\merged_flt_c10.imputed.vcf.gz";
        String outfileS = "M:\\production\\maf\\annotations\\siftScore\\hmp32TaxaGroup\\sample_merged_flt_c10.imputed.vcf.txt";
        int siteNumber = 6027602;
        int size = 10000;
        int[][] bound = FArrayUtils.getSubsetsIndicesBySubsetNumber(siteNumber, size);
        int[] indices = new int[size];
        for (int i = 0; i < size; i++) {
            indices[i] = bound[i][0];
        }
        Arrays.sort(indices);
        try {
            BufferedReader br = IoUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            String temp;
            while ((temp = br.readLine()).startsWith("##")) {}
            String header = temp;
            bw.write(header);
            bw.newLine();
            int cnt = -1;
            while ((temp = br.readLine()) != null) {
                cnt++;
                if (cnt%100000 == 0) System.out.println(cnt);
                if (Arrays.binarySearch(indices, cnt) < 0) continue;
                bw.write(temp);
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }
    
    private void deleteriousCountOnChr () {
        String transcriptFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\highConfidence_transcript.txt";
        String deleteriousSNPFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\class\\Non_Synonymous_Deleterious.txt";
        String synSNPFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\class\\Synonymous.txt";
        String geneFeatureFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
        String infoFileS = "E:\\Database\\InfoFile\\ChrLenCentPosi_agpV3.txt";
        String outfileDirS = "M:\\production\\maf\\annotations\\siftScore\\007_hmp321DeleDisOnChr\\";
        int chr = 1;
        int regionSize = 10000000;
        Table t = new Table (deleteriousSNPFileS);
        TIntArrayList pList = new TIntArrayList();
        TDoubleArrayList freList= new TDoubleArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            int c = t.getIntValue(i, 0);
            if (c != chr) continue;
            pList.add(t.getIntValue(i, 1));
            freList.add(t.getDoubleValue(i, 3));
        }
        int[] delePos = pList.toArray();
        double[] deleFre = freList.toArray();
        t = new Table (synSNPFileS);
        pList = new TIntArrayList();
        freList= new TDoubleArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            int c = t.getIntValue(i, 0);
            if (c != chr) continue;
            pList.add(t.getIntValue(i, 1));
            freList.add(t.getDoubleValue(i, 3));
        }
        int[] synPos = pList.toArray();
        double[] synFre = freList.toArray();
        t = new Table (infoFileS);
        int chrLength = 0;
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getIntValue(i, 0) == chr) {
                chrLength = t.getIntValue(i, 1);
                break;
            }
        }
        int[][] bound = FArrayUtils.getSubsetsIndicesBySubsetSize(chrLength, regionSize);
        int[] bondary = new int[bound.length];
        for (int i = 0; i < bound.length; i++) bondary[i] = bound[i][0]+1;
        int[] delSnpNumber = new int[bound.length];
        double[] delSnpNumberPerLine = new double[bound.length];
        for (int i = 0; i <  delePos.length; i++) {
            int index = Arrays.binarySearch(bondary, delePos[i]);
            if (index < 0) index = -index-2;
            delSnpNumber[index]++;
            delSnpNumberPerLine[index] += deleFre[i];
        }
        int[] synSnpNumber = new int[bound.length];
        double[] synSnpNumberPerLine = new double[bound.length];
        for (int i = 0; i <  synPos.length; i++) {
            int index = Arrays.binarySearch(bondary, synPos[i]);
            if (index < 0) index = -index-2;
            synSnpNumber[index]++;
            synSnpNumberPerLine[index] += synFre[i];
        }
        t = new Table (transcriptFileS);
        ArrayList<String> tranList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            tranList.add(t.content[i][0]);
        }
        String[] transName = tranList.toArray(new String[tranList.size()]);
        Arrays.sort(transName);
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        int[] cdsSiteNumber = new int[bound.length];
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            if (gf.getGeneChromosome(i) != chr) continue;
            String name = gf.getTranscriptName(i, 0);
            if (Arrays.binarySearch(transName, name) < 0) continue;
            List<Range> cds = gf.getCDSList(i, 0);
            int index1 = Arrays.binarySearch(bondary, gf.getTranscriptStart(i, 0));
            if (index1 < 0) index1 = -index1-2;
//            int index2 = Arrays.binarySearch(bound, gf.getTranscriptEnd(i, 0));
//            if (index2 < 0) index2 = -index2-2;
//            if (index1 == index2) {
//                
//            }
            for (int j = 0; j < cds.size(); j++) {
                cdsSiteNumber[index1]+=cds.get(j).getRangeSize();
            }
        }
        try {
            String outfileS = new File (outfileDirS, "deleteriousOnChr1.txt").getAbsolutePath();
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Position\tDeleteriousCount\tDeleteriousPerCDSSite\tSynonymousPerCDSSite\tDeleteriousSynonymousRatio");
            bw.newLine();
            for (int i = 0; i < bondary.length; i++) {
                StringBuilder sb = new StringBuilder();
                double del = (double)delSnpNumber[i]/cdsSiteNumber[i];
                double syn = (double)synSnpNumber[i]/cdsSiteNumber[i];
                sb.append(bondary[i]).append("\t").append(delSnpNumber[i]).append("\t").append(del).append("\t").append(syn).append("\t").append(del/syn);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            outfileS = new File (outfileDirS, "deleteriousPerLineOnChr1.txt").getAbsolutePath();
            bw = IoUtils.getTextWriter(outfileS);
            bw.write("Position\tDeleteriousCount\tDeleteriousPerCDSSitePerLine\tSynonymousPerCDSSitePerLine\tDeleteriousSynonymousRatio");
            bw.newLine();
            for (int i = 0; i < bondary.length; i++) {
                StringBuilder sb = new StringBuilder();
                double del = (double)delSnpNumberPerLine[i]/cdsSiteNumber[i];
                double syn = (double)synSnpNumberPerLine[i]/cdsSiteNumber[i];
                sb.append(bondary[i]).append("\t").append(delSnpNumber[i]).append("\t").append(del).append("\t").append(syn).append("\t").append(del/syn);
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
    
    private void countDeleteriousHmp32HighDepth () {
        String taxaSummaryFileS = "M:\\production\\maf\\annotations\\siftScore\\009_hmp321Depth\\taxaDepth.summary.txt";
        String addInfileS = "M:\\production\\maf\\annotations\\siftScore\\008_hmp321DeleCount\\additiveDeleterious_hmp32.txt";
        String addOutfileS = "M:\\production\\maf\\annotations\\siftScore\\008_hmp321DeleCount\\additiveDeleterious_hmp32_highDepth.txt";
        String recInfileS = "M:\\production\\maf\\annotations\\siftScore\\008_hmp321DeleCount\\reccesiveDeleterious_hmp32.txt";
        String recOutfileS = "M:\\production\\maf\\annotations\\siftScore\\008_hmp321DeleCount\\reccesiveDeleterious_hmp32_highDepth.txt";
        String taxaGroupFileS = "M:\\production\\maf\\annotations\\siftScore\\010_hmp321TaxaGroup\\hmp32.taxaGroup.txt";
        //1.54=3x, 2.75 = 5x
        double depthCut = 1.54;
        Table t = new Table (taxaSummaryFileS);
        ArrayList<String> taxaList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getDoubleValue(i, 2) < depthCut) continue;
            taxaList.add(t.content[i][0]);
        }
        String[] taxa = taxaList.toArray(new String[taxaList.size()]);
        Arrays.sort(taxa);
        t = new Table (taxaGroupFileS);
        HashMap<String, String> taxaGroupMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaGroupMap.put(t.content[i][0], t.content[i][1]);
        }
        t = new Table (addInfileS);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(addOutfileS);
            bw.write("Taxa\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth\tGroup\tRatio");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                int index = Arrays.binarySearch(taxa, t.content[i][0]);
                if (index < 0) continue;
                double ratio = Double.valueOf(t.content[i][1])/Double.valueOf(t.content[i][2]);
                bw.write(taxa[index]+"\t"+t.content[i][1]+"\t"+t.content[i][2]+"\t"+taxaGroupMap.get(taxa[index])+"\t"+String.valueOf(ratio));
                bw.newLine();
            }
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        t = new Table (recInfileS);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(recOutfileS);
            bw.write("Taxa\tDeleteriousCountPerLine\tSiteCountWithMinDepth\tGroup\tRatio");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                int index = Arrays.binarySearch(taxa, t.content[i][0]);
                if (index < 0) continue;
                double ratio = Double.valueOf(t.content[i][1])/Double.valueOf(t.content[i][2]);
                bw.write(taxa[index]+"\t"+t.content[i][1]+"\t"+t.content[i][2]+"\t"+taxaGroupMap.get(taxa[index])+"\t"+String.valueOf(ratio));
                bw.newLine();
            }
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void mkTaxaGroupFile () {
        String infileS = "M:\\production\\maf\\annotations\\siftScore\\010_hmp321TaxaGroup\\Fei_load.txt";
        String outfileS = "M:\\production\\maf\\annotations\\siftScore\\010_hmp321TaxaGroup\\hmp32.taxaGroup.txt";
        Table t = new Table (infileS);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Taxa\tGroup");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                String group = t.content[i][2];
                if (group.equals("")) group = "unknown";
                bw.write(t.content[i][0]+"\t"+group);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    private void mkDepthSummary () {
        String taxaDepthDirS = "M:\\production\\maf\\annotations\\siftScore\\009_hmp321Depth\\taxa\\";
        String taxaSummaryFileS = "M:\\production\\maf\\annotations\\siftScore\\009_hmp321Depth\\taxaDepth.summary.txt";
        File[] fs = new File (taxaDepthDirS).listFiles();
        Arrays.sort(fs);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(taxaSummaryFileS);
            bw.write("Taxa\tID\tMeanDepth");
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                Table t = new Table (fs[i].getAbsolutePath());
                String taxaName = t.header[0].replaceFirst("_siteDepth", "");
                double value = 0;
                for (int j = 0; j < t.getRowNumber(); j++) {
                    value+=t.getDoubleValue(j, 0);
                }
                bw.write(taxaName+"\t"+String.valueOf(i+1)+"\t"+String.valueOf((double)value/t.getRowNumber()));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void mkDepthOfHmp321 () {
        String hmpFileS = "O:\\Zea\\Genotypes\\WGS\\HapMap\\v3\\v321\\unimp\\vcf\\merged_flt_c10.vcf.gz";
        String hmpInfoFileS = "M:\\production\\maf\\annotations\\siftScore\\001_hmp321Info\\gerpAncestral\\hmp32Info_chr010.txt.gz";
        String taxaDepthDirS = "M:\\production\\maf\\annotations\\siftScore\\009_hmp321Depth\\taxa\\";
        int snpNum = 0;
        int size = 10000;
        try {
            BufferedReader br = IoUtils.getTextGzipReader(hmpInfoFileS);
            String temp = br.readLine();
            int cnt = 0; 
            while ((temp = br.readLine()) != null) {
                cnt++;
            }
            snpNum = cnt;
            int[] indices = new int[size];
            for (int i = 0; i < size; i++) {
                indices[i] = (int)(Math.random()*snpNum);
            }
            Arrays.sort(indices);
            br = IoUtils.getTextGzipReader(hmpFileS);
            while ((temp = br.readLine()).startsWith("##")) {}
            List<String> l = FStringUtils.fastSplit(temp, "\t");
            String[] taxa = new String[l.size()-9];
            for (int i = 0; i < taxa.length; i++) {
                taxa[i] = l.get(i+9);
            }
            TIntArrayList[] depthList = new TIntArrayList[taxa.length];
            for (int i = 0; i < taxa.length; i++) depthList[i] = new TIntArrayList();
            cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnt++;
                if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt)+" lines");
                int idx = Arrays.binarySearch(indices, cnt-1);
                if (idx < 0) continue;
                l = FStringUtils.fastSplit(temp, "\t");
                for (int i = 0; i < taxa.length; i++) {
                    String genoS = l.get(i+9);
                    if (genoS.startsWith(".")) {
                        depthList[i].add(0);
                        continue;
                    }
                    List<String> ll = FStringUtils.fastSplit(genoS, ":");
                    List<String> lll = FStringUtils.fastSplit(ll.get(1), ",");
                    int depth = Integer.valueOf(lll.get(0))+Integer.valueOf(lll.get(1));
                    depthList[i].add(depth);
                }
            }
            for (int i = 0; i < taxa.length; i++) {
                String outfileS = new File (taxaDepthDirS, FStringUtils.getNDigitNumber(5, i+1)+"depth.txt").getAbsolutePath();
                try {
                    BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                    bw.write(taxa[i]+"_siteDepth");
                    bw.newLine();
                    int[] depth = depthList[i].toArray();
                    for (int j = 0; j < depth.length; j++) {
                        bw.write(String.valueOf(depth[j]));
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
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    
    private void countDeleteriousHmp321 () {
        String hmpDirS = "O:\\Zea\\Genotypes\\WGS\\HapMap\\v3\\v321\\unimp\\vcf\\";
        String infoFileS = "E:\\Database\\InfoFile\\ChrLenCentPosi_agpV3.txt";
        String deleFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\class\\Non_Synonymous_Deleterious.txt";
        String addCountFileS = "M:\\production\\maf\\annotations\\siftScore\\008_hmp321DeleCount\\additiveDeleterious_hmp32.txt";
        String recCountFileS = "M:\\production\\maf\\annotations\\siftScore\\008_hmp321DeleCount\\reccesiveDeleterious_hmp32.txt";
        int minDepth = 2;//inclusive
        Table t = new Table (infoFileS);
        int chrNum = t.getRowNumber();
        int[] chrLength = new int[chrNum];
        ArrayList<Integer> chrList = new ArrayList();
        for (int i = 0; i < chrNum; i++) {
            chrLength[i] = t.getIntValue(i, 1);
            chrList.add(i+1);
        }
        t = new Table (deleFileS);
        TIntArrayList[] posList = new TIntArrayList[chrNum];
        TCharArrayList[] charList = new TCharArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) {
            posList[i] = new TIntArrayList();
            charList[i] = new TCharArrayList();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = t.getIntValue(i, 0)-1;
            if (t.content[i][4].length()>1) continue;
            posList[index].add(t.getIntValue(i, 1));
            charList[index].add(t.content[i][4].charAt(0));
        }
        int[][] delePos = new int[chrNum][];
        char[][] deleChar = new char[chrNum][];
        for (int i = 0; i < chrNum; i++) {
            delePos[i] = posList[i].toArray();
            deleChar[i] = charList[i].toArray();
        }
        String hmpChr10FileS = "merged_flt_c"+String.valueOf(10)+".vcf.gz";
        hmpChr10FileS = new File(hmpDirS, hmpChr10FileS).getAbsolutePath();
        String[] taxa = null;
        try {
            BufferedReader br = IoUtils.getTextGzipReader(hmpChr10FileS);
            String temp = br.readLine();
            while ((temp = br.readLine()).startsWith("##")) {}
            List<String> l = FStringUtils.fastSplit(temp, "\t");
            taxa = new String[l.size()-9];
            for (int i = 9; i < l.size(); i++) {
                taxa[i-9] = l.get(i);
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        int taxaNum = taxa.length;
        double[] addCount = new double[taxa.length];
        int[] recCount = new int[taxa.length];
        int[] siteWithMinDepthCount = new int[taxa.length];
        chrList.parallelStream().forEach(chr -> {
            String hmpFileS = "merged_flt_c"+String.valueOf(chr)+".vcf.gz";
            hmpFileS = new File(hmpDirS, hmpFileS).getAbsolutePath();
            BufferedReader br = IoUtils.getTextGzipReader(hmpFileS);
            int chrIndex = chr-1;
            try {
                String temp = br.readLine();
                while ((temp = br.readLine()).startsWith("##")) {}
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt)+" lines on chr "+String.valueOf(chr));
                    List<String> l = FStringUtils.fastSplit(temp.substring(0, 50), "\t");
                    int pos = Integer.valueOf(l.get(1));
                    int index = Arrays.binarySearch(delePos[chrIndex], pos);
                    if (index < 0) continue;
                    l = FStringUtils.fastSplit(temp, "\t");
                    int[] idx = new int[2];
                    if (l.get(3).charAt(0) == deleChar[chrIndex][index]) {
                        idx[0] = 0; idx[1] = 1;
                    }
                    else idx[0] = 1; idx[1] = 0;
                    for (int i = 0; i < taxaNum; i++) {
                        String genoS = l.get(i+9);
                        if (genoS.startsWith(".")) continue;
                        List<String> ll = FStringUtils.fastSplit(genoS, ":");
                        List<String> lll = FStringUtils.fastSplit(ll.get(1), ",");
                        int depth = Integer.valueOf(lll.get(0))+Integer.valueOf(lll.get(1));
                        if (depth < minDepth) continue;
                        lll = FStringUtils.fastSplit(ll.get(0), "/");
                        int v1 = Integer.valueOf(lll.get(0));
                        int v2 = Integer.valueOf(lll.get(1));
                        int sum = 0;
                        if (v1 == idx[0]) sum++;
                        if (v2 == idx[0]) sum++;
                        if (sum == 0) {}
                        else if (sum == 1) {
                            addCount[i] += 0.5;
                        }
                        else {
                            addCount[i] += 1;
                            recCount[i] += 1;
                        }
                        siteWithMinDepthCount[i]++;
                    }
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }); 
        try {
            BufferedWriter bw = IoUtils.getTextWriter(addCountFileS);
            bw.write("Taxa\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth");
            bw.newLine();
            for (int i = 0; i < addCount.length; i++) {
                bw.write(taxa[i]+"\t"+String.valueOf(addCount[i])+"\t"+String.valueOf(siteWithMinDepthCount[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            bw = IoUtils.getTextWriter(recCountFileS);
            bw.write("Taxa\tDeleteriousCountPerLine\tSiteCountWithMinDepth");
            bw.newLine();
            for (int i = 0; i < recCount.length; i++) {
                bw.write(taxa[i]+"\t"+String.valueOf(recCount[i])+"\t"+String.valueOf(siteWithMinDepthCount[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void mkBarplotOfSNPs () {
        String infileDirS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\class\\";
        String countFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\classCount.txt";
        String mafDistrubutionFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\mafSFS.txt";
        String dafDistrubutionFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\dafSFS.txt";
        int sampleSize = 10000;
        File[] fs = new File(infileDirS).listFiles();
        int size = 100;
        double[] bound = new double[size];
        for (int i = 1; i < bound.length; i++) {
            bound[i] = (double)1/size*i;
        }
        double[][] mafFrequency = new double[fs.length][size];
        double[][] dafFrequency = new double[fs.length][size];
        int[] count = new int[fs.length];
        int[] dafCount = new int[fs.length];
        TDoubleArrayList[] mafList = new TDoubleArrayList[fs.length];
        TDoubleArrayList[] dafList = new TDoubleArrayList[fs.length];
        for (int i = 0; i < fs.length; i++) {
            mafList[i] = new TDoubleArrayList();
            dafList[i] = new TDoubleArrayList();
            String infileS = fs[i].getAbsolutePath();
            Table t = new Table(infileS);
            count[i] = t.getRowNumber();
            for (int j = 0; j < t.getRowNumber(); j++) {
                double value = t.getDoubleValue(j, 3);
                mafList[i].add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) index = -index -2;
                mafFrequency[i][index]++;
                if (!t.content[j][5].startsWith("N")) {
                    value = t.getDoubleValue(j, 5);
                    dafList[i].add(value);
                    index = Arrays.binarySearch(bound, value);
                    if (index < 0) index = -index -2;
                    dafFrequency[i][index]++;
                    dafCount[i]++;
                }
            }
            for (int j = 0; j < mafFrequency[i].length; j++) {
                mafFrequency[i][j] = mafFrequency[i][j]/count[i];
                dafFrequency[i][j] = dafFrequency[i][j]/dafCount[i];
            }
        }
        try {
            BufferedWriter bw =IoUtils.getTextWriter(countFileS);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < fs.length; i++) {
                sb.append(fs[i].getName().replaceFirst(".txt", "")).append("\t");
            }
            sb.deleteCharAt(sb.length()-1);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < count.length-1; i++) {
                bw.write(String.valueOf(count[i])+"\t");
            }
            bw.write(String.valueOf(count[count.length-1]));
            bw.newLine();
            bw.flush();
            bw.close();
            bw =IoUtils.getTextWriter(mafDistrubutionFileS);
            bw.write("MAF");
            for (int i = 0; i < fs.length; i++) {
                bw.write("\t"+fs[i].getName().replaceFirst(".txt", ""));
            }
            bw.newLine();
            for (int i = 0; i < mafFrequency[0].length; i++) {
                bw.write(String.valueOf(bound[i]));
                for (int j = 0; j < mafFrequency.length; j++) {
                    bw.write("\t"+mafFrequency[j][i]);
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
            bw =IoUtils.getTextWriter(dafDistrubutionFileS);
            bw.write("DAF");
            for (int i = 0; i < fs.length; i++) {
                bw.write("\t"+fs[i].getName().replaceFirst(".txt", ""));
            }
            bw.newLine();
            for (int i = 0; i < dafFrequency[0].length; i++) {
                bw.write(String.valueOf(bound[i]));
                for (int j = 0; j < dafFrequency.length; j++) {
                    bw.write("\t"+dafFrequency[j][i]);
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Step 1: filtered out low quality gene models. Non/syn ratio less than 2.5
     * Step 2: classify SNPs into different classes based on only high confidence gene models. Deleterious mutations (SIFT less than 0.05, GERP greater than 2)
     */
    private void classifySNPs3 () {
        String geneFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
        String geneSummaryFileS = "M:\\production\\maf\\annotations\\siftScore\\002_transcriptSummary\\transcriptSummary.txt";
        String infileDirS = "M:\\production\\maf\\annotations\\siftScore\\003_hmp321SiftGerpUScore\\";
        String outfileDirS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\class\\";
        String transcriptFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\highConfidence_transcript.txt";
        String[] types = {"Synonymous", "Non_Synonymous_Tolerent", "StopLoss", "StopGain", "Non_Synonymous_Deleterious", "Non_Synonymous_Deleterious_High_GERP"};
        
        double gerpCut = 2;
        double nonSynRatioCut = 2.5;//filter incorrect gene model
        double nonSynGeneCut = 0.02;
        Table t = new Table (geneSummaryFileS);
        boolean[] ifOut = new boolean[t.getRowNumber()];
        ArrayList<String> tranNameList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            //if (!t.content[i][4].equals("1")) continue;
            if (!t.content[i][14].equals("1")) continue; //filter incorrect gene model by picking up gerp aligned genes
            double ratio = t.getDoubleValue(i, 9);
            double nonSyn = t.getDoubleValue(i, 8);
            if (Double.isNaN(ratio)) {
                if (nonSyn > nonSynGeneCut) continue;
                tranNameList.add(t.content[i][0]);
                ifOut[i] = true;
            }
            else {
                if (ratio > nonSynRatioCut) continue;
                tranNameList.add(t.content[i][0]);
                ifOut[i] = true;
            }
        }
        t.writeTable(transcriptFileS, ifOut);
        String[] tranName = tranNameList.toArray(new String[tranNameList.size()]);
        System.out.println(String.valueOf(tranName.length) + "genes kept");
        Arrays.sort(tranName);
        GeneFeature gf = new GeneFeature(geneFileS);
        gf.sortGeneByStartPosition();
        ArrayList<Range> cdsList = new ArrayList();
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            String query = gf.getTranscriptName(i, 0);
            int index = Arrays.binarySearch(tranName, query);
            if (index < 0) continue;
            List<Range> cdsl = gf.getCDSList(i, 0);
            for (int j = 0; j < cdsl.size(); j++) {
                cdsList.add(cdsl.get(j));
            }
        }
        Ranges rs = new Ranges(cdsList, "cds");
        rs.sortByStartPosition();
        File[] fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        try {
            BufferedWriter[] bw = new BufferedWriter[types.length];
            for (int i = 0; i < types.length; i++) {
                String outfileS = new File (outfileDirS, types[i]+".txt").getAbsolutePath();
                bw[i] = IoUtils.getTextWriter(outfileS);
                bw[i].write("Chr\tPos\tMinorAllele\tMAF\tDerivedAllele\tDAF");
                bw[i].newLine();
            }
            for (int i = 0; i < fs.length; i++) {
                t = new Table (fs[i].getAbsolutePath());
                int chr = t.getIntValue(0, 0);
                for (int j = 0; j < t.getRowNumber(); j++) {
                    int pos = t.getIntValue(j, 1);
                    if (!rs.isInRanges(chr, pos)) continue;
                    double maf = t.getDoubleValue(j, 8);
                    String da = "NA";
                    String daf = "NA";
                    if (t.content[j][4].length() == 1) {
                        String an = t.content[j][4];
                        if (an.equals(t.content[j][5])) {
                            da = t.content[j][6];
                            daf = t.content[j][8];
                        }
                        else if (an.equals(t.content[j][6])) {
                            da = t.content[j][5];
                            daf = t.content[j][7];
                        }
                    }
                    if (t.content[j][11].equals("Syn")) {
                        if (t.content[j][12].startsWith("N")) continue;
                        double sift = t.getDoubleValue(j, 12);
                        bw[0].write(String.valueOf(chr)+"\t"+String.valueOf(pos)+"\t"+t.content[j][6]+"\t"+String.valueOf(maf)+"\t"+da+"\t"+daf);
                        bw[0].newLine();
                    }
                    if (t.content[j][11].equals("Non")) {
                        if (t.content[j][12].startsWith("N")) continue;
                        double sift = t.getDoubleValue(j, 12);
                        if (sift < 0.05) {
                            if (t.getDoubleValue(j, 14) > gerpCut) {
                                bw[5].write(String.valueOf(chr)+"\t"+String.valueOf(pos)+"\t"+t.content[j][6]+"\t"+String.valueOf(maf)+"\t"+da+"\t"+daf);
                                bw[5].newLine();
                            }
                            bw[4].write(String.valueOf(chr)+"\t"+String.valueOf(pos)+"\t"+t.content[j][6]+"\t"+String.valueOf(maf)+"\t"+da+"\t"+daf);
                            bw[4].newLine();
                        }
                        else {
                            bw[1].write(String.valueOf(chr)+"\t"+String.valueOf(pos)+"\t"+t.content[j][6]+"\t"+String.valueOf(maf)+"\t"+da+"\t"+daf);
                            bw[1].newLine();
                        }
                        
                    }
                    if (t.content[j][11].equals("StopLoss")) {
                        bw[2].write(String.valueOf(chr)+"\t"+String.valueOf(pos)+"\t"+t.content[j][6]+"\t"+String.valueOf(maf)+"\t"+da+"\t"+daf);
                        bw[2].newLine();
                    }
                    if (t.content[j][11].equals("StopGain")) {
                        bw[3].write(String.valueOf(chr)+"\t"+String.valueOf(pos)+"\t"+t.content[j][6]+"\t"+String.valueOf(maf)+"\t"+da+"\t"+daf);
                        bw[3].newLine();
                    }
                }
            }
            for (int i = 0; i < types.length; i++) {
                bw[i].flush();
                bw[i].close();
            }
            
        }
        catch (Exception e) {
            e.printStackTrace();
        }

    }
    
    
    /**
     * Show relationship between sift and gerp
     */
    private void siftAndGerp () {
        String infileS = "M:\\production\\maf\\annotations\\siftScore\\002_hmp32SiftGerpUScore\\chr010.hmpSiftGerpUScore.txt";
        String outfileS = "M:\\production\\maf\\annotations\\siftScore\\004_siftAndGerp\\siftAndGerp.pdf";
        Table t = new Table (infileS);
        int size = 5000;
        int[] index = FArrayUtils.getRandomIntArray(t.getRowNumber(), 10000);
        TDoubleArrayList siftList = new TDoubleArrayList();
        TDoubleArrayList gerpList = new TDoubleArrayList();
        for (int i = 0; i < index.length; i++) {
            if (t.content[index[i]][11].startsWith("N")) continue;
            if (!t.content[index[i]][10].startsWith("Non")) continue;
            siftList.add(Double.valueOf(t.content[index[i]][11]));
            gerpList.add(Double.valueOf(t.content[index[i]][14]));
        }
        ScatterPlot s = new ScatterPlot (siftList.toArray(), gerpList.toArray());
        s.setColor(255, 0, 0, 20);
        s.setPlottingCharacter(16);
        s.setXLab("SIFT score");
        s.setYLab("GERP score");
        s.saveGraph(outfileS);
    }
    
    /**
     * Combine SIFT GERP into one table in coding variants of hmp321
     * Note: sift annotation using a different version of gene annotation, containing RNA genes, GFF3 do not have RNA, I use GFF3
     */
    private void combineHmpSiftGerpUScore () {
        String geneFileS = "/workdir/mingh/Zea_mays.AGPv3.26.gf.txt";
        String hmpInfoDirS = "/workdir/mingh/001_hmp321Info/";
        String gerpDirS = "/workdir/mingh/gerp/";
        String siftDirS = "/workdir/mingh/sift/";
        String uscoreDirS = "/workdir/mingh/combined_6Genomes/";
        String combineDirS = "/workdir/mingh/combined/";
        new File(combineDirS).mkdir();
        GeneFeature gf = new GeneFeature(geneFileS);
        HashSet<String> transSet = new HashSet();
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            transSet.add(gf.getTranscriptName(i, 0));
        }
        File[] fs =  new File(hmpInfoDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(hmpf -> {
            int chr = Integer.valueOf(hmpf.getName().split("chr")[1].replaceFirst(".txt.gz", ""));
            String siftFileS = new File(siftDirS, String.valueOf(chr)+".gz").getAbsolutePath();
            String gerpFileS = new File (gerpDirS, "roast.chrom."+String.valueOf(chr)+".msa.in.rates.full").getAbsolutePath();
            String hmpFileS = hmpf.getAbsolutePath();
            String uscoreFileS = new File (uscoreDirS, "chr"+FStringUtils.getNDigitNumber(3, chr)+"_uniqueness.txt.gz").getAbsolutePath();
            String combineFileS = new File(combineDirS, "chr"+FStringUtils.getNDigitNumber(3, chr)+".hmpSiftGerpUScore.txt").getAbsolutePath();
            HashMap<String, String> posSiftMap = new HashMap();
            HashMap<Integer, Set<String>> posTransMap = new HashMap();
            try {
                BufferedReader br = IoUtils.getTextGzipReader(siftFileS);
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    List<String> l = FStringUtils.fastSplit(temp, "\t");
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt)+"\t"+siftFileS);
                    if (!l.get(6).startsWith("CDS")) continue;
                    if (!transSet.contains(l.get(3))) continue;
                    if (l.get(3).startsWith("GRM") && (!l.get(3).endsWith("01"))) continue;
                    int pos = Integer.valueOf(l.get(0));
                    if (posTransMap.get(pos) == null) {
                        TreeSet<String> tl = new TreeSet();
                        tl.add(l.get(3));
                        posTransMap.put(pos, tl);
                    }
                    else {
                        Set<String> tl = posTransMap.get(pos);
                        tl.add(l.get(3));
                        posTransMap.put(pos, tl);
                    }
                    String type = null;
                    if (l.get(7).equals("*")) {
                        if (l.get(8).equals("*")) {
                            type = "StopStop";
                        }
                        else {
                            type = "StopLoss";
                        }
                    }
                    else {
                        if (l.get(8).equals("*")) {
                            type = "StopGain";
                        }
                        else {
                            if (l.get(7).equals(l.get(8))) {
                                type = "Syn";
                                //type = l.get(7)+"/"+l.get(8)+"/"+"Syn";
                            }
                            else {
                                type = "Non";
                                //type = l.get(7)+"/"+l.get(8)+"/"+"Non";
                                }
                            }
                        }
                    StringBuilder sb = new StringBuilder();
                    sb.append(l.get(0)).append("\t").append(l.get(2));
                    posSiftMap.put(sb.toString(), type+"\t"+l.get(10));
                }
                br.close();
                
                
                br = IoUtils.getTextGzipReader(hmpFileS);
                BufferedReader brG = IoUtils.getTextReader(gerpFileS);
                BufferedReader brU = IoUtils.getTextGzipReader(uscoreFileS);
                BufferedWriter bw = IoUtils.getTextWriter(combineFileS);
                bw.write("Chr	Pos	Ref	Alt	Ancestral(Gerp)	Major	Minor	MajorAlleleFrequency	MinorAlleleFrequency	SiteDepth	HetCount\tMutationClass\tSift\tGerpTreeLength\tGerp\tUScore\tTranscripts");
                bw.newLine();
                temp = br.readLine();
                String tempG = null;
                String tempU = null;
                cnt = 0;
                int cntG = 0;
                int cntU = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt)+"\tchr"+String.valueOf(chr));
                    List<String> l = FStringUtils.fastSplit(temp, "\t");
                    if (l.get(3).contains("<") || l.get(3).contains(",")) continue;
                    int currentPos = Integer.valueOf(l.get(1));
                    for (int i = cntG; i < currentPos; i++) {
                        tempG = brG.readLine();
                    }
                    cntG=currentPos;
                    List<String> ul = null;
                    for (int i = cntU; i < currentPos; i++) {
                        tempU = brU.readLine();
                        ul = FStringUtils.fastSplit(tempU, "\t");
                    }
                    cntU = currentPos;
                    
                    StringBuilder sb = new StringBuilder();
                    sb.append(l.get(1)).append("\t").append(l.get(3));
                    String key = sb.toString();
                    String hmpS = posSiftMap.get(key);
                    if (hmpS == null) continue;
                    sb = new StringBuilder(); 
                    if (tempG==null) tempG="0\t0";
                    sb.append(temp).append("\t").append(hmpS).append("\t").append(tempG).append("\t").append(ul.get(2)).append("\t");
                    Set<String> tranSet = posTransMap.get(currentPos);
                    String[] trans = tranSet.toArray(new String[tranSet.size()]);
                    for (int i = 0; i < trans.length; i++) {
                        sb.append(trans[i]).append(";");
                    }
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                brG.close();
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        
    }
    
    private void mkHmp32MafPlot () {
        String infileS = "M:\\production\\maf\\annotations\\siftScore\\hmp32Info\\hmp32Info_chr010.txt.gz";
        String outfileS = "M:\\production\\maf\\annotations\\siftScore\\hmp32MafPlot\\hmp32MafFre.txt";
        TDoubleArrayList list = new TDoubleArrayList();
        try {
            BufferedReader br = IoUtils.getTextGzipReader(infileS);
            String temp = br.readLine();
            int sum = 0;
            while ((temp = br.readLine()) != null) {
                List<String> l = FStringUtils.fastSplit(temp);
                if (l.get(7).contains(",")) continue;
                list.add(Double.valueOf(l.get(7)));
                sum++;
                if (sum%100000 == 0) System.out.println(sum);
            }
            int size = 20;
            double[] bound = new double[size];
            int[] count = new int[size];
            double[] ratio = new double[size];
            for (int i = 1; i < bound.length; i++) {
                bound[i] = (double)0.5/size*i;
            }
            double[] value = list.toArray();
            for (int i = 0; i < value.length; i++) {
                int index = Arrays.binarySearch(bound, value[i]);
                if (index < 0) index = -index-2;
                count[index]++;
            }
            for (int i = 0; i < ratio.length; i++) {
                ratio[i] = (double)count[i]/sum;
            }
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("MAF\tFrequency");
            bw.newLine();
            for (int i = 0; i < ratio.length; i++) {
                bw.write(String.valueOf(bound[i])+"\t"+String.valueOf(ratio[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void updateHmp321InfoWithAncestralAllele () {
        String inDirS = "M:\\production\\maf\\annotations\\siftScore\\001_hmp321Info\\tripsacumAncestral\\";
        String ancestralDirS = "M:\\production\\maf\\annotations\\gerp\\ancestralAllele\\";
        String outDirS = "M:\\production\\maf\\annotations\\siftScore\\001_hmp321Info\\gerpAncestral\\";
        File[] fs = new File(inDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            int chr = Integer.valueOf(f.getName().replaceFirst("hmp32Info_chr", "").replaceFirst(".txt.gz", ""));
            String outfileS = new File (outDirS, f.getName()).getAbsolutePath();
            String infileS = new File (ancestralDirS, "chr"+FStringUtils.getNDigitNumber(3, chr)+"ancestral.txt").getAbsolutePath();
            try {
                TIntArrayList posList = new TIntArrayList();
                ArrayList<String> baseList = new ArrayList();
                BufferedReader br = IoUtils.getTextReader(infileS);
                String temp = br.readLine();
                while ((temp = br.readLine()) != null) {
                    List<String> l = FStringUtils.fastSplit(temp);
                    if (l.get(1).length()>1) continue;
                    posList.add(Integer.valueOf(l.get(0)));
                    baseList.add(l.get(1));
                }
                br.close();
                System.out.println("Read ancestral allele from chr"+String.valueOf(chr));
                int[] pos = posList.toArray();
                BufferedWriter bw = IoUtils.getTextGzipWriter(outfileS);
                bw.write("Chr	Pos	Ref	Alt	Ancestral(Gerp)	Major	Minor	MajorAlleleFrequency	MinorAlleleFrequency	SiteDepth	HetCount");
                bw.newLine();
                br = IoUtils.getTextGzipReader(f.getAbsolutePath());
                br.readLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt)+" in chr"+String.valueOf(chr));
                    List<String> l = FStringUtils.fastSplit(temp);
                    int index = Arrays.binarySearch(pos, Integer.valueOf(l.get(1)));
                    if (index < 0) {
                        l.set(4, "NA");
                    }
                    else {
                        l.set(4, baseList.get(index));
                    }
                    StringBuilder sb = new StringBuilder(l.get(0));
                    for (int i = 1; i  < l.size(); i++) {
                        sb.append("\t").append(l.get(i));
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
        });
    }
    
    /**
     * Extract Hmp32 information, including chr, pos, major, minor, MAF, site depth, site count
     */
    private void extractHmp321Info () {
        String infileDirS = "/workdir/mingh/hmp";
        String outfileDirS = "/workdir/mingh/hmp32Info/";
//        String infileDirS = "M:\\production\\maf\\annotations\\siftScore\\test\\hmp";
//        String outfileDirS = "M:\\production\\maf\\annotations\\siftScore\\test\\hmpInfo";
        new File(outfileDirS).mkdir();
        File[] fs = new File(infileDirS).listFiles();
        List<File> fl = Arrays.asList(fs);
        fl.parallelStream().forEach(f -> {
            int chr = Integer.valueOf(f.getName().split("_c")[1].replaceFirst(".fixed.vcf.gz", ""));
            String outfileS = new File(outfileDirS, "hmp32Info_chr"+FStringUtils.getNDigitNumber(3, chr)+".txt.gz").getAbsolutePath();
            String temp = null;
            try {
                BufferedReader br = IoUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedWriter bw = IoUtils.getTextGzipWriter(outfileS);
                bw.write("Chr\tPos\tRef\tAlt\tAncestral(Tripsacum)\tMajor\tMinor\tMajorAlleleFrequency\tMinorAlleleFrequency\tSiteDepth\tHetCount");
                bw.newLine();
                while ((temp = br.readLine()).startsWith("##")) {}
                List<String> l = FStringUtils.fastSplit(temp);
                int index = 0; 
                for (int i = 0; i < l.size(); i++) {
                    if (l.get(i).equals("TDD39103")) {
                        index = i;
                        break;
                    }
                }
                int cnt = 1;
                while ((temp = br.readLine()) != null) {
                    bw.write(this.getHmpInfoString(temp, index));
                    bw.newLine();
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt)+" sites process from "+f.getName());
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }    
        });
        
    }
    
    private String getHmpInfoString (String temp, int triIndex) {
        List<String> l = FStringUtils.fastSplit(temp, "\t");
        StringBuilder sb = new StringBuilder();
        sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(3)).append("\t").append(l.get(4)).append("\t");
        List<String> al = FStringUtils.fastSplit(l.get(4), ",");
        int alleleNum = 1+ al.size();
        int[] alleleCount = new int[alleleNum];
        String[] alleles = new String[alleleNum];
        alleles[0] = l.get(3);
        for (int i = 0; i < al.size(); i++) {
            alleles[i+1] = al.get(i);
        }
        String tripGeno = l.get(triIndex);
        String triS = "";
        if (tripGeno.startsWith(".")) {
            triS = "NA";
        }
        else {
            String[] tem = tripGeno.split(":")[0].split("/");
            if (tem[0].equals(tem[1])) {
                triS = alleles[Integer.valueOf(tem[0])];
            }
            else {
                triS = alleles[Integer.valueOf(tem[0])]+","+alleles[Integer.valueOf(tem[1])];;
            }
        }
        sb.append(triS).append("\t");
        al = FStringUtils.fastSplit(l.get(7), ";");
        int siteDepth = Integer.valueOf(al.get(0).split("DP=")[1]);
        int hetCount = 0;
        for (int i = 9; i < l.size(); i++) {
            if (l.get(i).startsWith(".")) continue;
            List<String> ll = FStringUtils.fastSplit(l.get(i), ":");
            List<String> lll = FStringUtils.fastSplit(ll.get(0), "/");
            alleleCount[Integer.valueOf(lll.get(0))]++;
            alleleCount[Integer.valueOf(lll.get(1))]++;
            if (lll.get(0).equals(lll.get(1))) continue;
            hetCount++;
        }
        int sum = 0;
        for (int i = 0; i < alleleCount.length; i++) {
            sum+=alleleCount[i];
        }
        double[] fre = new double[alleleCount.length];
        for (int i = 0; i < fre.length; i++) fre[i] = (double)alleleCount[i]/sum;
        int[] index = FArrayUtils.getIndexByDescendingValue(alleleCount);
        sb.append(alleles[index[0]]).append("\t");
        StringBuilder minorSB = new StringBuilder();
        for (int i = 1; i < alleles.length; i++) {
            minorSB.append(alleles[index[i]]).append(",");
        }
        minorSB.deleteCharAt(minorSB.length()-1);
        sb.append(minorSB.toString()).append("\t").append((float)fre[index[0]]).append("\t");
        minorSB = new StringBuilder();
        for (int i = 1; i < alleles.length; i++) {
            minorSB.append((float)fre[index[i]]).append(",");
        }
        minorSB.deleteCharAt(minorSB.length()-1);
        sb.append(minorSB.toString());
        sb.append("\t").append(siteDepth).append("\t").append(hetCount);
        return sb.toString();
    }
}
