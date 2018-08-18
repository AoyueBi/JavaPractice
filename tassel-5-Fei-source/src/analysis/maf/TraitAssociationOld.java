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
import graphcis.r.Rgraphics;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import rcaller.RCaller;
import rcaller.RCode;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class TraitAssociationOld {
    
    public TraitAssociationOld () {
        //this.getFilteredGene();
        //this.getSNPDataSet();
        //this.mktestData(); //optional
        //this.countDeleteriousHmp32();
        this.mergeDeleteriousAndTrait();
        this.mkR2s();
    }
    
    private void mkR2s () {
        String inDirS = "M:\\production\\maf\\traitAssociation\\inGene\\traitAndDele\\reccesive\\";
        String outfileS = "M:\\production\\maf\\traitAssociation\\inGene\\deleCount\\deleHeight.txt";
        File[] fs = new File (inDirS).listFiles();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Gerp\tMAF\tSift\tR2WithHeight\tRWithHeight");
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                double[] results = this.getR2(fs[i].getAbsolutePath());
                String temp = fs[i].getName().replaceAll("_", "\t");
                bw.write(temp+"\t"+String.valueOf(results[0])+"\t"+String.valueOf(results[1]));
                bw.newLine();
                System.out.println(String.valueOf(i)+"\t"+String.valueOf(results[0])+"\t"+String.valueOf(results[1]));
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private double[] getR2 (String infileS) {
        infileS = infileS.replaceAll("\\\\", "/");
        RCaller caller = new RCaller();
        caller.setRscriptExecutable(Rgraphics.RPath);
        RCode rCode = new RCode();
        rCode.addRCode("dataFile <- (\""+infileS+"\")");
        rCode.addRCode("data <- read.table(dataFile, header=TRUE, sep=\"\\t\")");
        rCode.addRCode("model <- lm (data$Burden ~ data$Distance)");
        rCode.addRCode("model2 <- lm(data$PH_BLUP_peiffer ~ model$residuals)");
        rCode.addRCode("r2 <- summary(model2)$r.squared");
        rCode.addRCode ("r <- cor (data$PH_BLUP_peiffer, model$residuals)");
        rCode.addRCode("results <- list(m1 = r2, m2 = r)");
        //System.out.println(rCode.getCode().toString());
        caller.setRCode(rCode);
        caller.runAndReturnResult("results");
        double[] result1 = caller.getParser().getAsDoubleArray("m1");
        double[] result2 = caller.getParser().getAsDoubleArray("m2");
        double[] results = new double[2];
        results[0] = result1[0];
        results[1] = result2[0];
        return results;
    }
    
    private void mergeDeleteriousAndTrait () {
        String sourceFileS = "M:\\production\\maf\\annotations\\siftScore\\hmp32DeleCount\\reccesiveDeleterious_merge_sub_noOutlier.txt";
        String inDirS = "M:\\production\\maf\\traitAssociation\\inGene\\deleCount\\reccesive\\";
        String outDirS = "M:\\production\\maf\\traitAssociation\\inGene\\traitAndDele\\reccesive\\";
        Table t = new Table (sourceFileS);
        HashMap<String, String> taxaDistanceMap = new HashMap();
        HashMap<String, String> taxaHeightMap = new HashMap();
        HashMap<String, String> taxaDTAMap = new HashMap();
        HashMap<String, String> taxaGroupMap = new HashMap ();
        String[] taxa = new String[t.getRowNumber()];
        for (int i = 0; i < taxa.length; i++) {
            taxa[i] = t.content[i][0];
            taxaGroupMap.put(taxa[i], t.content[i][3]);
            taxaDistanceMap.put(taxa[i], t.content[i][5]);
            taxaHeightMap.put(taxa[i], t.content[i][6]);
            taxaDTAMap.put(taxa[i], t.content[i][7]);
        }
        Arrays.sort(taxa);
        File[] fs = new File(inDirS).listFiles();
        
        for (int i = 0; i < fs.length; i++) {
            t = new Table (fs[i].getAbsolutePath());
            try {
                String outfileS = new File (outDirS, fs[i].getName()).getAbsolutePath();
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                bw.write("Taxa\tDeleteriousCountPerLine\tSiteCountWithMinDepth\tBurden\tGroup\tDistance\tPH_BLUP_peiffer\tDTA_BLUP_peiffer");
                bw.newLine();
                for (int j = 0; j < t.getRowNumber(); j++) {
                    if (Arrays.binarySearch(taxa, t.content[j][0]) < 0) continue;
                    StringBuilder sb = new StringBuilder();
                    sb.append(t.content[j][0]).append("\t").append(t.content[j][1]).append("\t").append(t.content[j][2]).append("\t").append(t.content[j][3]).append("\t").append(taxaGroupMap.get(t.content[j][0]));
                    sb.append("\t").append(taxaDistanceMap.get(t.content[j][0])).append("\t").append(taxaHeightMap.get(t.content[j][0])).append("\t").append(taxaDTAMap.get(t.content[j][0]));
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
    
    private void countDeleteriousHmp32 () {
        String hmpDirS = "Q:\\Zea\\Genotypes\\HapMap\\v3\\LDKNNi\\unimp_hmp321_vcf\\";
        //String hmpDirS = "M:\\production\\maf\\traitAssociation\\inGene\\test\\hmp\\";
        String snpSetDirS = "M:\\production\\maf\\traitAssociation\\inGene\\snpSubsets\\gerpMafSift\\";
        String stopSetDirS = "M:\\production\\maf\\traitAssociation\\inGene\\snpSubsets\\stopGain\\";
        String mafDirS = "M:\\production\\maf\\annotations\\siftScore\\hmp32SiftGerpUScore\\";
        String addDirS = "M:\\production\\maf\\traitAssociation\\inGene\\deleCount\\additive\\";
        String recDirS = "M:\\production\\maf\\traitAssociation\\inGene\\deleCount\\reccesive\\";
        int minDepth = 2;
        int chrNum = 10;
        ArrayList<Integer> chrList = new ArrayList();
        for (int i = 0; i < chrNum; i++) {
            chrList.add(i+1);
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
        int[][][] position = new int[chrNum][][];
        String[][] type = new String[chrNum][];
        HashMap<Integer, Character>[] posMinorMap = new HashMap[chrNum];
        HashSet<Integer>[] existPosSet = new HashSet[chrNum];
        chrList.parallelStream().forEach(chr -> {
            String snpFileS = "chr"+FStringUtils.getNDigitNumber(3, chr)+".hmp32sub.txt";
            snpFileS = new File (snpSetDirS, snpFileS).getAbsolutePath();
            String stopFileS = "chr"+FStringUtils.getNDigitNumber(3, chr)+".stop.txt";
            stopFileS = new File (stopSetDirS, stopFileS).getAbsolutePath();
            String mafFileS = "chr"+FStringUtils.getNDigitNumber(3, chr)+".hmpSiftGerpUScore.txt";
            mafFileS = new File (mafDirS, mafFileS).getAbsolutePath();
            int chrIndex = chr -1;
            posMinorMap[chrIndex] = new HashMap();
            existPosSet[chrIndex] = new HashSet();
            Table t = new Table (mafFileS);
            for (int i = 0; i < t.getRowNumber(); i++) {
                posMinorMap[chrIndex].put(Integer.valueOf(t.content[i][1]), t.content[i][5].charAt(0));
            }
            t = new Table (stopFileS);
            int stopNum = t.getRowNumber();
            try {
                BufferedReader br = IoUtils.getTextReader(snpFileS);
                int size = Integer.valueOf(br.readLine().split("\t")[0]);
      
                position[chrIndex] = new int[size][];
                type[chrIndex] = new String[size];
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith(">")) {
                        List<String> l = FStringUtils.fastSplit(temp, "\t");
                        type[chrIndex][cnt] = l.get(1)+"_"+l.get(2)+"_"+l.get(3);
                        int total = Integer.valueOf(l.get(4));
                        position[chrIndex][cnt] = new int[total];
                        for (int i = 0; i < total; i++) {
                            int pos = Integer.valueOf(br.readLine());
                            position[chrIndex][cnt][i] = pos;
                            existPosSet[chrIndex].add(pos);
                        }
                        
                        Arrays.sort(position[chrIndex][cnt]);
                        cnt++;
                    }
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        double[][] addCount = new double[type[0].length][taxa.length];
        int[][] recCount = new int[type[0].length][taxa.length];
        int[][] siteWithMinDepthCount = new int[type[0].length][taxa.length];
        System.out.println("Start processing Hmp32");
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
                    if (cnt%100000 == 0) {
                        System.out.println(String.valueOf(cnt)+" lines on chr "+String.valueOf(chr));
                    }
                    List<String> l = FStringUtils.fastSplit(temp.substring(0, 50), "\t");
                    int pos = Integer.valueOf(l.get(1));
                    if (!existPosSet[chrIndex].contains(pos)) continue;
                    l = FStringUtils.fastSplit(temp, "\t");
                    for (int i = 0; i < position[chrIndex].length; i++) {
                        int index = Arrays.binarySearch(position[chrIndex][i], pos);
                        if (index < 0) continue;
                        int[] idx = new int[2];
                        if (l.get(3).charAt(0) == posMinorMap[chrIndex].get(Integer.valueOf(l.get(1)))) {
                            idx[0] = 0; idx[1] = 1;
                        }
                        else idx[0] = 1; idx[1] = 0;
                        
                        for (int j = 0; j < taxaNum; j++) {
                            String genoS = l.get(j+9);
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
                                addCount[i][j] += 0.5;
                            }
                            else {
                                addCount[i][j] += 1;
                                recCount[i][j] += 1;
                            }
                            siteWithMinDepthCount[i][j]++;
                        }
                    }
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        try {
            for (int i = 0; i < type[0].length; i++) {
                String addOutfileS = new File (addDirS, type[0][i]).getAbsolutePath();
                String recOutfileS = new File (recDirS, type[0][i]).getAbsolutePath();
                BufferedWriter bw = IoUtils.getTextWriter(addOutfileS);
                bw.write("Taxa\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth\tBurden");
                bw.newLine();
                for (int j = 0; j < addCount[i].length; j++) {
                    StringBuilder sb = new StringBuilder();
                    double burden = (double)addCount[i][j]/siteWithMinDepthCount[i][j];
                    sb.append(taxa[j]).append("\t").append(addCount[i][j]).append("\t").append(siteWithMinDepthCount[i][j]).append("\t").append(burden);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                bw = IoUtils.getTextWriter(recOutfileS);
                bw.write("Taxa\tDeleteriousCountPerLine\tSiteCountWithMinDepth\tBurden");
                bw.newLine();
                for (int j = 0; j < recCount[i].length; j++) {
                    StringBuilder sb = new StringBuilder();
                    double burden = (double)recCount[i][j]/siteWithMinDepthCount[i][j];
                    sb.append(taxa[j]).append("\t").append(recCount[i][j]).append("\t").append(siteWithMinDepthCount[i][j]).append("\t").append(burden);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    private void mktestData () {
        String hmpDirS = "Q:\\Zea\\Genotypes\\HapMap\\v3\\LDKNNi\\unimp_hmp321_vcf\\";
        String outDirS = "M:\\production\\maf\\traitAssociation\\inGene\\test\\hmp\\";
        int size = 100000;
        ArrayList<Integer> chrList = new ArrayList();
        for (int i = 0; i < 10; i++) {
            chrList.add(i+1);
        }
        chrList.parallelStream().forEach(chr -> {
            String infileS = new File (hmpDirS,"merged_flt_c"+String.valueOf(chr)+".vcf.gz").getAbsolutePath();
            String outfileS = new File (outDirS,"merged_flt_c"+String.valueOf(chr)+".vcf.gz").getAbsolutePath();
            try {
                BufferedReader br = IoUtils.getTextGzipReader(infileS);
                BufferedWriter bw = IoUtils.getTextGzipWriter(outfileS);
                for (int j = 0; j < size; j++) {
                    bw.write(br.readLine());
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
    
    private void getSNPDataSet () {
        String infileDirS = "M:\\production\\maf\\annotations\\siftScore\\hmp32SiftGerpUScore\\";
        String geneFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
        String filetedGenefileS = "M:\\production\\maf\\traitAssociation\\inGene\\filteredGene.txt";
        String outDirS = "M:\\production\\maf\\traitAssociation\\inGene\\snpSubsets\\gerpMafSift\\";
        String stopOutDirS = "M:\\production\\maf\\traitAssociation\\inGene\\snpSubsets\\stopGain\\";
        double gerpStart = 0;
        double gerpStep = 0.4;
        double mafStart = 0.05;
        double mafStep = 0.05;
        double siftStart = 0.005;
        double siftStep = 0.005;
        int size = 10;
        double[] gerpCut = new double[size];
        double[] mafCut = new double[size];
        double[] siftCut = new double[size];
        for (int i = 0; i < size; i++) {
            gerpCut[i] = gerpStart+i*gerpStep;
            mafCut[i] = mafStart+i*mafStep;
            siftCut[i] = siftStart+i*siftStep;
        }
        Table t = new Table (filetedGenefileS);
        String[] gene = new String[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            gene[i] = t.content[i][0];
        }
        Arrays.sort(gene);
        GeneFeature gf = new GeneFeature(geneFileS);
        gf.sortGeneByStartPosition();
        ArrayList<Range> cdsList = new ArrayList();
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            String query = gf.getTranscriptName(i, 0);
            int index = Arrays.binarySearch(gene, query);
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
            for (int i = 0; i < fs.length; i++) {
                TIntArrayList stopGainList = new TIntArrayList();
                TIntArrayList[][][] list =new TIntArrayList[size][size][size];
                for (int j = 0; j < gerpCut.length; j++) {
                    for (int k = 0; k < mafCut.length; k++) {
                        for (int l = 0; l < siftCut.length; l++) {
                            list[j][k][l] = new TIntArrayList();
                        }
                    }
                }
                
                t = new Table (fs[i].getAbsolutePath());
                for (int j = 0; j < t.getRowNumber(); j++) {
                    int pos = t.getIntValue(j, 1);
                    if (t.content[j][11].startsWith("N")) {
                        if (t.content[j][10].equals("StopGain")) {
                            stopGainList.add(pos);
                        }
                        continue;
                    }
                    double maf = t.getDoubleValue(j, 7);
                    double sift = t.getDoubleValue(j, 11);
                    double gerp = t.getDoubleValue(j, 14);
                    for (int m = 0; m < gerpCut.length; m++) {
                        for (int k = 0; k < mafCut.length; k++) {
                            for (int l = 0; l < siftCut.length; l++) {
                                if (gerp> gerpCut[m] && maf < mafCut[k] && sift < siftCut[l]) {
                                    list[m][k][l].add(pos);
                                }
                            }
                        }
                    }
                }
                String outfileS = fs[i].getName().replaceFirst("hmpSiftGerpUScore.txt", "hmp32sub.txt");
                outfileS = new File (outDirS, outfileS).getAbsolutePath();
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                bw.write(String.valueOf(size*size*size)+"\tGerp\tMaf\tSift");
                bw.newLine();
                for (int j = 0; j < gerpCut.length; j++) {
                    for (int k = 0; k < mafCut.length; k++) {
                        for (int l = 0; l < siftCut.length; l++) {
                            int[] pos = list[j][k][l].toArray();
                            bw.write(">\t"+String.valueOf(gerpCut[j]) + "\t" + String.valueOf(mafCut[k]) +"\t"+String.valueOf(siftCut[l]) + "\t"+pos.length);
                            bw.newLine();
                            for (int m = 0; m < pos.length; m++) {
                                bw.write(String.valueOf(pos[m]));
                                bw.newLine();
                            }
                        }
                    }
                }
                bw.flush();
                bw.close();
                outfileS = fs[i].getName().replaceFirst("hmpSiftGerpUScore.txt", "stop.txt");
                outfileS = new File (stopOutDirS, outfileS).getAbsolutePath();
                bw = IoUtils.getTextWriter(outfileS);
                bw.write("StopCodon");
                bw.newLine();
                int[] pos = stopGainList.toArray();
                for (int j = 0; j < pos.length; j++) {
                    bw.write(String.valueOf(pos[j]));
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            
            
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
        
        TIntArrayList[][][] list =new TIntArrayList[size][size][size];
        for (int i = 0; i < gerpCut.length; i++) {
            for (int j = 0; j < mafCut.length; j++) {
                for (int k = 0; k < siftCut.length; k++) {
                    
                }
            }
        }
    }
    
    private void getFilteredGene () {
        String geneMutationFileS = "M:\\production\\maf\\annotations\\siftScore\\hmp32SNPClassMAF\\mutationPerGene.txt";
        String outfileS = "M:\\production\\maf\\traitAssociation\\inGene\\filteredGene.txt";
        double nonSynGeneCut = 0.04;
        double delGeneCut = 0.015;
        Table t = new Table (geneMutationFileS);
        ArrayList<String> tranNameList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getDoubleValue(i, 5) > nonSynGeneCut) continue;
            if (t.getDoubleValue(i, 7) > delGeneCut) continue;
            tranNameList.add(t.content[i][0]);
        }
        String[] trans = tranNameList.toArray(new String[tranNameList.size()]);
        Arrays.sort(trans);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("FilteredGene");
            bw.newLine();
            for (int i = 0; i < trans.length; i++) {
                bw.write(trans[i]);
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
