/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maf;

import format.Table;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntByteHashMap;
import gnu.trove.map.hash.TIntDoubleHashMap;
import gnu.trove.set.hash.TIntHashSet;
import graphcis.r.ScatterPlot;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.taxa.TaxaList;
import net.openhft.koloboke.collect.map.hash.HashByteByteMap;
import net.openhft.koloboke.collect.map.hash.HashIntByteMap;
import net.openhft.koloboke.collect.map.hash.HashIntByteMaps;
import net.openhft.koloboke.collect.map.hash.HashIntDoubleMap;
import net.openhft.koloboke.collect.map.hash.HashIntDoubleMaps;
import net.openhft.koloboke.collect.map.hash.HashShortIntMapFactory;
import net.openhft.koloboke.collect.map.hash.HashShortIntMaps;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import utils.BaseCoder;
import utils.FArrayUtils;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class NAMTraitAssociation {
    
    public NAMTraitAssociation () {
        //this.mkCDSGenotype();
        //this.regularPipe();
        this.effectPipe();
    }
    
    void effectPipe () {
        //this.namSNPSegregation();
        //this.reformatTrait();
        //this.mergeTrait();
        //this.mergeGenenotype();
        //this.ramdomSampleSNPForPCA();
        //this.splitByTrait();
        //this.mergeGLMByMutationType();
        //this.mkSigGLMByMutationType();
    }
    
    void mkSigGLMByMutationType () {
        String glmDirS = "M:\\production\\maf\\traitAssociation\\namEffect\\glmResult\\byTraits\\";
        String delDirS = "M:\\production\\maf\\traitAssociation\\namEffect\\deleEffect\\";
        String synDirS = "M:\\production\\maf\\traitAssociation\\namEffect\\synEffect\\";
        String delSigDirS = "M:\\production\\maf\\traitAssociation\\namEffect\\deleEffectSig\\";
        String synSigDirS = "M:\\production\\maf\\traitAssociation\\namEffect\\synEffectSig\\";
        double thresh = 0.05;
        File[] fs = new File(glmDirS).listFiles();
        String[] traits = new String[fs.length];
        double[] pThresh = new double[fs.length];
        double[] meanE = new double[fs.length];
        double[] sdE = new double[fs.length];
        for (int i = 0; i < fs.length; i++) {
            traits[i] = fs[i].getName().replaceFirst(".gwas.txt", "");
            Table t = new Table (fs[i].getAbsolutePath());
            TDoubleArrayList pList = new TDoubleArrayList();
            for (int j = 0; j < t.getRowNumber(); j++) {
                pList.add(t.getDoubleValue(j, 6));
            }
            double[] ps = pList.toArray();
            Arrays.sort(ps);
            pThresh[i] = ps[(int)(ps.length*thresh)];
            TDoubleArrayList sigEList = new TDoubleArrayList();
            for (int j = 0; j < t.getRowNumber(); j++) {
                if (t.getDoubleValue(j, 6) > pThresh[i]) continue;
                sigEList.add(t.getDoubleValue(j, 5));
            }
            double[] sigE = sigEList.toArray();
            DescriptiveStatistics d = new DescriptiveStatistics(sigE);
            meanE[i] = d.getMean();
            sdE[i] = d.getStandardDeviation();
        }
        this.mkSigGLMByMutationTypeSub(traits, pThresh, meanE, sdE, delDirS, delSigDirS);
        this.mkSigGLMByMutationTypeSub(traits, pThresh, meanE, sdE, synDirS, synSigDirS);
    }
    
    void mkSigGLMByMutationTypeSub (String[] traits, double[] pThresh, double[] meanE, double[] sdE, String inputDirS, String outputDirS) {
        File outDir = new File(outputDirS);
        outDir.mkdir();
        for (int i = 0; i < traits.length; i++) {
            String inputFileS = new File (inputDirS, traits[i]+".effect.txt").getAbsolutePath();
            String outputFileS = new File (outputDirS, traits[i]+".effect.sig.txt").getAbsolutePath();
            try {
                BufferedReader br = IoUtils.getTextReader(inputFileS);
                BufferedWriter bw = IoUtils.getTextWriter(outputFileS);
                bw.write(br.readLine()+"\tZEstimate");
                bw.newLine();
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    List<String> l = FStringUtils.fastSplit(temp);
                    if (l.get(7).startsWith("N")) continue;
                    if (l.get(3).startsWith("N")) continue;
                    double p = Double.valueOf(l.get(7));
                    if (p > pThresh[i]) continue;
                    double z = (Double.valueOf(l.get(6))-meanE[i])/sdE[i];
                    bw.write(temp+"\t"+String.valueOf(z));
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
    
    void mergeGLMByMutationType () {
        String traitDirS = "M:\\production\\maf\\traitAssociation\\namEffect\\glmResult\\byTraits";
        String siftGerpDirS = "M:\\production\\maf\\annotations\\siftScore\\003_hmp321SiftGerpUScore\\";
        String delClassFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\class\\Non_Synonymous_Deleterious_High_GERP.txt";
        String synClassFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\class\\Synonymous.txt";
        String delDirS = "M:\\production\\maf\\traitAssociation\\namEffect\\deleEffect";
        String synDirS = "M:\\production\\maf\\traitAssociation\\namEffect\\synEffect";
        this.mergeGLM(delClassFileS, siftGerpDirS, traitDirS, delDirS);
        this.mergeGLM(synClassFileS, siftGerpDirS, traitDirS, synDirS);
    }
    
    void mergeGLM (String mutationFileS, String siftGerpDirS, String traitDirS, String outputDirS) {
        int chrNum = 10;
        File outputDir = new File (outputDirS);
        outputDir.mkdir();
        Table t = new Table (mutationFileS);
        TIntHashSet[] posSet = new TIntHashSet[chrNum];
        int[][] poss = new int[chrNum][];
        TDoubleArrayList[] mafList = new TDoubleArrayList[chrNum];
        ArrayList<String>[] dafList = new ArrayList[chrNum];
        String[][] sift = new String[chrNum][];
        String[][] gerp = new String[chrNum][];
        String[][] estimate = new String[chrNum][];
        String[][] pvalue = new String[chrNum][];
        for (int i = 0; i < chrNum; i++) {
            posSet[i] = new TIntHashSet();
            mafList[i] = new TDoubleArrayList();
            dafList[i] = new ArrayList();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            int chrIndex = t.getIntValue(i, 0)-1;
            posSet[chrIndex].add(t.getIntValue(i, 1));
            mafList[chrIndex].add(t.getDoubleValue(i, 3));
            dafList[chrIndex].add(t.content[i][5]);
        }
        File[] fs = new File(siftGerpDirS).listFiles();
        for (int i = 0; i < chrNum; i++) {
            sift[i] = new String[posSet[i].size()];
            gerp[i] = new String[posSet[i].size()];
            poss[i] = posSet[i].toArray();
            Arrays.sort(poss[i]);
            t = new Table (fs[i].getAbsolutePath());
            int cnt = 0;
            for (int j = 0; j < t.getRowNumber(); j++) {
                int pos = t.getIntValue(j, 1);
                if (!posSet[i].contains(pos)) continue;
                sift[i][cnt] = t.content[j][12];
                gerp[i][cnt] = t.content[j][14];
                cnt++;
            }
        }
        fs = new File (traitDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            String trait = fs[i].getName().replaceFirst(".gwas.txt", "");
            String outfileS = new File (outputDir, trait+".effect.txt").getAbsolutePath();
            for (int j = 0; j < chrNum; j++) {
                estimate[j] = new String[posSet[j].size()];
                pvalue[j] = new String[posSet[j].size()];
                for (int k = 0; k < estimate[j].length; k++) {
                    estimate[j][k] = "NA";
                    pvalue[j][k] = "NA";
                }
            }
            try {
                BufferedReader br = IoUtils.getTextReader(fs[i].getAbsolutePath());
                String temp = br.readLine();
                while ((temp = br.readLine()) != null) {
                    List<String> l = FStringUtils.fastSplit(temp);
                    int chrIndex = Integer.valueOf(l.get(1)) - 1;
                    int pos = Integer.valueOf(l.get(2));
                    int index = Arrays.binarySearch(poss[chrIndex], pos);
                    if (index < 0) continue;
                    estimate[chrIndex][index] = l.get(5);
                    pvalue[chrIndex][index] = l.get(6);
                }
                br.close();
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                bw.write("Chr\tPos\tMAF\tDAF\tSIFT\tGERP\tEstimate\tPValue");
                bw.newLine();
                for (int j = 0; j < chrNum; j++)    {
                    for (int k = 0; k < poss[j].length; k++) {
                        StringBuilder sb = new StringBuilder();
                        sb.append(j+1).append("\t").append(poss[j][k]).append("\t").append(mafList[j].get(k)).append("\t").append(dafList[j].get(k)).append("\t");
                        sb.append(sift[j][k]).append("\t").append(gerp[j][k]).append("\t").append(estimate[j][k]).append("\t").append(pvalue[j][k]);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
                bw.flush();
                bw.close();
                System.out.println(trait);
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
    }
    
    void splitByTrait () {
        String effectFileS = "M:\\production\\maf\\traitAssociation\\namEffect\\glmResult\\genotype_traits_pcs.txt";
        String pFileS = "M:\\production\\maf\\traitAssociation\\namEffect\\glmResult\\statistics_traits_pcs.txt";
        String outDirS = "M:\\production\\maf\\traitAssociation\\namEffect\\glmResult\\byTraits\\";
        ArrayList<String> traitList = new ArrayList();
        HashMap<String, ArrayList<StringBuilder>> traitResultMap = new HashMap();
        HashMap<String, HashMap<String, Double>> traitPMap = new HashMap();
        HashMap<String, Double>[] markerPMap = null;
        String[] base = {"A", "T", "G", "C"};
        Arrays.sort(base);
        try {
            BufferedReader br = IoUtils.getTextReader(effectFileS);
            String temp = br.readLine();
            String currentTrait = "";
            String currentMarker = "";
            StringBuilder sb = new StringBuilder();
            int cnt1 = -1;
            int cnt2 = -1;
            double v1 = -1;
            double v2 = -1;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnt++;
                if (cnt%5000000 == 0) System.out.println(cnt+"\t estimate");
                List<String> l = FStringUtils.fastSplit(temp);
                if (Arrays.binarySearch(base, l.get(5)) < 0) continue;
                if (currentTrait.equals(l.get(0))) {
                    ArrayList<StringBuilder> currentSbList = traitResultMap.get(currentTrait);
                    if (l.get(1).equals(currentMarker)) {
                        if (cnt1 == -1) {
                            //not possible
                        }
                        else  {
                            if (cnt2 == -1) {
                                cnt2 = Integer.valueOf(l.get(4));
                                if (cnt1 < cnt2) {
                                    int mid = cnt1;
                                    cnt1 = cnt2; cnt2 = mid;
                                }
                                v2 = Double.valueOf(l.get(6));
                                sb.append("\t").append(cnt1).append("\t").append(cnt2).append("\t").append(Math.abs(v1-v2));
                                currentSbList.add(sb);
                                cnt1 = -1; cnt2 = -1;
                            }
                            else {
                                //not possible
                            }
                        }
                    }
                    else {
                        currentMarker = l.get(1);
                        sb = new StringBuilder();
                        sb.append(l.get(1)).append("\t").append(l.get(2)).append("\t").append(l.get(3));
                        cnt1 = Integer.valueOf(l.get(4));
                        v1 = Double.valueOf(l.get(6));
                    }                
                }
                else {
                    currentTrait = l.get(0);
                    traitList.add(l.get(0));
                    ArrayList<StringBuilder> sbList = new ArrayList();
                    traitResultMap.put(l.get(0), sbList);
                    
                    if (l.get(1).equals(currentMarker)) {
                        //not possible
                    }
                    else {
                        currentMarker = l.get(1);
                        sb = new StringBuilder();
                        sb.append(l.get(1)).append("\t").append(l.get(2)).append("\t").append(l.get(3));
                        cnt1 = Integer.valueOf(l.get(4));
                        v1 = Double.valueOf(l.get(6));
                    }
                }
            }
            br.close();
            br = IoUtils.getTextReader(pFileS);
            br.readLine();
            markerPMap = new HashMap[traitList.size()];
            for (int i = 0; i < markerPMap.length; i++) {
                markerPMap[i] = new HashMap();
                traitPMap.put(traitList.get(i), markerPMap[i]);
            }
            cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnt++;
                if (cnt%5000000 == 0) System.out.println(cnt+"\t p");
                List<String> l = FStringUtils.fastSplit(temp);
                if (l.get(4).startsWith("N")) continue;
                HashMap<String, Double> currentMap = traitPMap.get(l.get(0));
                currentMap.put(l.get(1), Double.valueOf(l.get(5)));
            }
            for (int i = 0; i < traitList.size(); i++) {
                String outfileS = new File (outDirS, traitList.get(i)+".gwas.txt").getAbsolutePath();
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                bw.write("Marker\tChr\tPos\tNAM_Major\tNAM_Minor\tEstimate\tPValue");
                bw.newLine();
                ArrayList<StringBuilder> currentSbList = traitResultMap.get(traitList.get(i));
                HashMap<String, Double> currentMap = traitPMap.get(traitList.get(i));
                for (int j = 0; j < currentSbList.size(); j++) {
                    sb = currentSbList.get(j);
                    List<String> l = FStringUtils.fastSplit(sb.toString());
                    double p = 1;
                    Double pd = currentMap.get(l.get(0));
                    if (pd != null) p = pd;
                    sb.append("\t").append(p);
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
    
    void ramdomSampleSNPForPCA () {
        String positionFileS = "M:\\production\\maf\\annotations\\siftScore\\001_hmp321Info\\gerpAncestral\\hmp32Info_chr001.txt.gz";
        String hmpFileS = "M:\\production\\maf\\traitAssociation\\source\\genotype\\NAM_HM321_Unimp.hmp.txt.gz";
        String ramdomHmpFileS = "M:\\production\\maf\\traitAssociation\\namEffect\\pca\\random.hmp.txt.gz";
        int size = 20000;
        int snpNum = 0;
        try {
            BufferedReader br = IoUtils.getTextGzipReader(positionFileS);
            String temp = br.readLine();
            while ((temp = br.readLine()) != null) {
                snpNum++;
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        int[] indices = FArrayUtils.getNonredundantRandomIntArray(snpNum, size);
        Arrays.sort(indices);
        try {
            BufferedWriter bw = IoUtils.getTextGzipWriter(ramdomHmpFileS);
            BufferedReader br = IoUtils.getTextGzipReader(hmpFileS);
            String temp = br.readLine();
            bw.write(temp);
            bw.newLine();
            int cnt = 0;
            int count = 0;
            while ((temp = br.readLine()) != null) {
                if (Arrays.binarySearch(indices, cnt) >= 0) {
                    bw.write(temp);
                    bw.newLine();
                    count++;
                    if (count >= size) break;
                }
                cnt++;
                if (cnt%1000000==0) System.out.println(cnt);
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    void mergeGenenotype () {
        String inDirS = "M:\\production\\maf\\traitAssociation\\genotype\\hmp321CDSSNP\\";
        String mergedGenoFileS = "M:\\production\\maf\\traitAssociation\\genotype\\hmp321CDS.hmp.txt.gz";
        File[] fs = new File(inDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        try {
            BufferedWriter bw = IoUtils.getTextGzipWriter(mergedGenoFileS);
            for (int i = 0; i < fs.length; i++) {
                System.out.println(fs[i].getName());
                BufferedReader br = IoUtils.getTextReader(fs[i].getAbsolutePath());
                String temp = br.readLine();
                if (i == 0) {
                    bw.write(temp); bw.newLine();
                }
                while ((temp = br.readLine()) != null) {
                    bw.write(temp);
                    bw.newLine();
                }
                br.close();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    void mergeTrait() {
        String inDirS = "M:\\production\\maf\\traitAssociation\\phenotype\\namTasselFormat\\";
        String outfileS = "M:\\production\\maf\\traitAssociation\\phenotype\\namTrait_merged.txt";
        File[] fs = new File(inDirS).listFiles();
        HashSet<String> taxaSet = new HashSet();
        HashMap<String, String>[] valueMaps = new HashMap[fs.length];
        String[] traits = new String[fs.length];
        for (int i = 0; i < fs.length; i++) {
            valueMaps[i] = new HashMap();
            Table t = new Table (fs[i].getAbsolutePath());
            traits[i] = t.header[1];
            for (int j = 0; j < t.getRowNumber(); j++) {
                taxaSet.add(t.content[j][0]);
                valueMaps[i].put(t.content[j][0], t.content[j][1]);
            }
        }
        String[] taxa = taxaSet.toArray(new String[taxaSet.size()]);
        Arrays.sort(taxa);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder("<Trait>");
            for (int i = 0; i < traits.length; i++) {
                sb.append("\t").append(traits[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < taxa.length; i++) {
                sb = new StringBuilder(taxa[i]);
                for (int j = 0; j < traits.length; j++) {
                    String value = valueMaps[j].get(taxa[i]);
                    if (value == null) value = "NA";
                    sb.append("\t").append(value);
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
    
    void reformatTrait () {
       String inDirS = "M:\\production\\maf\\traitAssociation\\source\\multiblup_formatted\\";
       String outDirS = "M:\\production\\maf\\traitAssociation\\phenotype\\namTasselFormat\\";
       File[] fs = new File(inDirS).listFiles();
       List<File> fList = Arrays.asList(fs);
       fList.parallelStream().forEach(f -> {
           String outfileS = new File (outDirS, f.getName().replaceFirst("multiblup", "blup_Tassel")).getAbsolutePath();
           try {
               BufferedReader br = IoUtils.getTextReader(f.getAbsolutePath());
               BufferedWriter bw = IoUtils.getTextWriter(outfileS);
               bw.write("<Trait>\t"+f.getName().replaceFirst("_multiblup.txt", ""));
               bw.newLine();
               String temp = null;
               while ((temp = br.readLine()) != null) {
                   String[] tem = temp.split(" ");
                   bw.write(tem[0]+"\t"+tem[2]);
                   bw.newLine();
               }
               bw.flush();
               bw.close();
               br.close();
           }
           catch (Exception e) {
               e.printStackTrace();
           }
       });
    }
    
    void namSNPSegregation () {
        String hmpDirS = "M:\\production\\maf\\traitAssociation\\genotype\\hmp321CDSSNP\\";
        String cdsSNPInfoDirS = "M:\\production\\maf\\annotations\\siftScore\\003_hmp321SiftGerpUScore\\";
        String outDirS = "M:\\production\\maf\\traitAssociation\\namEffect\\hmp321CDSNPSegregation\\";
        double segregationThresh = 0.3;
        File[] fs = new File(hmpDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        new File (outDirS).mkdir();
        fList.parallelStream().forEach(f -> {
            int chr = Integer.valueOf(f.getName().replaceFirst("NAM_hm321_unimp_chr", "").replaceFirst(".hmp.txt", ""));
            String infileS = new File (cdsSNPInfoDirS, "chr"+FStringUtils.getNDigitNumber(3, chr)+".hmpSiftGerpUScore.txt").getAbsolutePath();
            String outfileS = new File (outDirS, "SNPSegregation_chr"+FStringUtils.getNDigitNumber(3, chr)+".txt").getAbsolutePath();
            try {
                BufferedReader br = IoUtils.getTextReader(f.getAbsolutePath());
                BufferedReader bri = IoUtils.getTextReader(infileS);
                String temp = br.readLine();
                bri.readLine();
                List<String> l = FStringUtils.fastSplit(temp);
                HashSet<String> familySet = new HashSet();
                for (int i = 0; i < l.size(); i++) {
                    String taxon = l.get(i);
                    if (taxon.startsWith("Z0")) {
                        familySet.add(taxon.substring(0, 4));
                    }
                }
                String[] family = familySet.toArray(new String[familySet.size()]);
                TIntArrayList[] familyIndexList = new TIntArrayList[family.length];
                int[][] familyIndex = new int[family.length][];
                for (int i = 0; i < family.length; i++) {
                    familyIndexList[i] = new TIntArrayList();
                }
                Arrays.sort(family);
                for (int i = 0; i < l.size(); i++) {
                    String taxon = l.get(i);
                    if (!taxon.startsWith("Z0")) continue;
                    int familyInd = Arrays.binarySearch(family, taxon.substring(0, 4));
                    familyIndexList[familyInd].add(i);
                }
                for (int i = 0; i < family.length; i++) {
                    familyIndex[i] = familyIndexList[i].toArray();
                    Arrays.sort(familyIndex[i]);
                }
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                StringBuilder sb = new StringBuilder();
                sb.append("Chr\tPos\tRef\tAlt\tNumOfSegregationFamily");
                for (int i = 0; i < family.length; i++) sb.append("\t").append(family[i]);
                bw.write(sb.toString());
                bw.newLine();
                while ((temp = br.readLine()) != null) {
                    l = FStringUtils.fastSplit(temp);
                    String tempi = bri.readLine();
//                    if (l.get(1).length() == 1) continue;
//                    if (l.get(1).startsWith("N")) continue;
                    List<String> li = FStringUtils.fastSplit(tempi);
                    String alt = li.get(3);
                    double[] frequency = new double[family.length];
                    int[] nonMissing = new int[family.length];
                    int nSegre = 0;
                    for (int i = 0; i < family.length; i++) {
                        for (int j = 0; j < familyIndex[i].length; j++) {
                            if (l.get(familyIndex[i][j]).equals("N")) continue;
                            nonMissing[i]++;
                            if (l.get(familyIndex[i][j]).equals(alt)) frequency[i]++;
                        }
                        if (nonMissing[i] != 0) {
                            frequency[i] = frequency[i]/nonMissing[i];
                        }
                        else {
                            frequency[i] = Double.NaN;
                            //System.out.print("frequency NaN occurs");
                        }
                        if (frequency[i] > segregationThresh && frequency[i] < (1-segregationThresh)) nSegre++;
                    }
                    sb = new StringBuilder();
                    sb.append(chr).append("\t").append(li.get(1)).append("\t").append(li.get(2)).append("\t").append(3).append("\t").append(nSegre);
                    for (int i = 0; i < family.length; i++) {
                        sb.append("\t").append(frequency[i]);
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
    
    void regularPipe () {
        this.mkLoad();
        this.alltaxaLoadAndPhenotype();
        this.familyLoadAndPhenotype();
    }
    
    public void familyLoadAndPhenotype () {
        String phenotypeDirS = "O:\\Zea\\Phenotypes\\BLUPs\\USNAM\\Wallace41Traits\\multiblup_formatted\\";
        String loadFileS = "M:\\production\\maf\\traitAssociation\\nam\\load.txt";
        String familyFileS = "M:\\production\\maf\\traitAssociation\\nam\\namFamily.txt";
        String outFileS = "M:\\production\\maf\\traitAssociation\\nam\\familyLoad\\familyLoadTrait.txt";
        String pdfDirS = "M:\\production\\maf\\traitAssociation\\nam\\familyLoad\\pdf\\";
        new File (pdfDirS).mkdir();
        Table t = new Table (loadFileS);
        HashSet<String> familySet = new HashSet();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (!t.content[i][0].startsWith("Z0")) continue;
            familySet.add(t.content[i][0].substring(0, 4));
        }
        String[] family = familySet.toArray(new String[familySet.size()]);
        Arrays.sort(family);
        HashMap<String, Double>[] taxaLoadMap = new HashMap[family.length];
        for (int i = 0; i < taxaLoadMap.length; i++) taxaLoadMap[i] = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (!t.content[i][0].startsWith("Z0")) continue;
            String query = t.content[i][0].substring(0, 4);
            int index = Arrays.binarySearch(family, query);
            taxaLoadMap[index].put(t.content[i][0], t.getDoubleValue(i, 10));
        }
        HashMap<String, String> familyTaxaMap = new HashMap();
        t = new Table (familyFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            familyTaxaMap.put(t.content[i][0], t.content[i][1]);
        }
        File[] fs = new File(phenotypeDirS).listFiles();
        String[] traits = new String[fs.length];
        double[][] rs = new double[traits.length][family.length];
        for (int i = 0; i < fs.length; i++) {
            traits[i] = fs[i].getName();
            try {
                ArrayList<String> taxaList = new ArrayList();
                ArrayList<Double> pheValueList = new ArrayList();
                BufferedReader br = IoUtils.getTextReader(fs[i].getAbsolutePath());
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    String[] tem = temp.split(" ");
                    if (tem[2].startsWith("N")) continue;
                    taxaList.add(tem[0]);
                    pheValueList.add(Double.valueOf(tem[2]));
                }
                br.close();
                String[] taxa = taxaList.toArray(new String[taxaList.size()]);
                Double[] value = pheValueList.toArray(new Double[pheValueList.size()]);
                for (int j = 0; j < family.length; j++) {
                    taxaList = new ArrayList();
                    pheValueList = new ArrayList();
                    ArrayList<Double> loadList = new ArrayList();
                    for (int k = 0; k < taxa.length; k++) {
                        if (!taxa[k].startsWith(family[j])) continue;
                        Double l = taxaLoadMap[j].get(taxa[k]);
                        if (l == null) continue;
                        taxaList.add(taxa[k]);
                        pheValueList.add(value[k]);
                        loadList.add(l);
                    }
                    double[] pheno = new double[pheValueList.size()];
                    double[] load = new double[loadList.size()];
                    for (int k = 0; k < pheno.length; k++) {
                        pheno[k] = pheValueList.get(k);
                        load[k] = loadList.get(k);
                    }
                    double[][] matrix = new double[pheno.length][2];
                    for (int k = 0; k < pheno.length; k++) {
                        matrix[k][0] = pheno[k];
                        matrix[k][1] = load[k];
                    }
                    double r = Double.NaN;
                    if (pheno.length != 0) {
                        PearsonsCorrelation cor = new PearsonsCorrelation(matrix);
                        r = cor.getCorrelationMatrix().getEntry(0, 1);
                        double pv = cor.getCorrelationPValues().getEntry(0, 1);
                    }
                    rs[i][j] = r;
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outFileS);
            StringBuilder sb = new StringBuilder("Class");
            for (int i = 0; i < family.length; i++) {
                sb.append("\t").append(familyTaxaMap.get(family[i]));
                //sb.append("\t").append(family[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < traits.length; i++) {
                sb = new StringBuilder(traits[i]);
                for (int j = 0; j < family.length; j++) {
                    sb.append("\t").append(rs[i][j]);
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
    
    public void alltaxaLoadAndPhenotype () {
        String phenotypeDirS = "O:\\Zea\\Phenotypes\\BLUPs\\USNAM\\Wallace41Traits\\multiblup_formatted\\";
        String loadFileS = "M:\\production\\maf\\traitAssociation\\nam\\load.txt";
        String loadPhenoDirS = "M:\\production\\maf\\traitAssociation\\nam\\allTaxaLoad\\loadPheno\\";
        String loadPhenoCorrFileS = "M:\\production\\maf\\traitAssociation\\nam\\allTaxaLoad\\loadPheno_r2.txt\\";
        Table t = new Table (loadFileS);
        String[] types = {"all", "centro", "distal"};
        HashMap<String, Double>[] taxaLoadMap = new HashMap[types.length];
        for (int i = 0; i < taxaLoadMap.length; i++) {
            taxaLoadMap[i] = new HashMap();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
//            taxaLoadMap[0].put(t.content[i][0], t.getDoubleValue(i, 7));
//            taxaLoadMap[1].put(t.content[i][0], t.getDoubleValue(i, 8));
//            taxaLoadMap[2].put(t.content[i][0], t.getDoubleValue(i, 9));
            taxaLoadMap[0].put(t.content[i][0], t.getDoubleValue(i, 10));
            taxaLoadMap[1].put(t.content[i][0], t.getDoubleValue(i, 11));
            taxaLoadMap[2].put(t.content[i][0], t.getDoubleValue(i, 12));
        }
        File loadPhenoDir = new File (loadPhenoDirS);
        loadPhenoDir.mkdir();
        File[] subDir = new File[types.length];
        for (int i = 0; i < subDir.length; i++) {
            subDir[i] = new File(loadPhenoDirS, types[i]);
            subDir[i].mkdir();
        }
        File[] fs = new File(phenotypeDirS).listFiles();
        String[] traits = new String[fs.length];
        double[][] r2 = new double[fs.length][types.length];
        double[][] ps = new double[fs.length][types.length];
        for (int i = 0; i < fs.length; i++) {
            traits[i] = fs[i].getName();
            try {
                ArrayList<String> taxaList = new ArrayList();
                ArrayList<Double> pheValueList = new ArrayList();
                BufferedReader br = IoUtils.getTextReader(fs[i].getAbsolutePath());
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    String[] tem = temp.split(" ");
                    if (tem[2].startsWith("N")) continue;
                    taxaList.add(tem[0]);
                    pheValueList.add(Double.valueOf(tem[2]));
                }
                br.close();
                String[] taxa = taxaList.toArray(new String[taxaList.size()]);
                Double[] value = pheValueList.toArray(new Double[pheValueList.size()]);
                for (int j = 0; j < types.length; j++) {
                    taxaList = new ArrayList();
                    pheValueList = new ArrayList();
                    ArrayList<Double> loadList = new ArrayList();
                    for (int k = 0; k < taxa.length; k++) {
                        Double l = taxaLoadMap[j].get(taxa[k]);
                        if (l == null) continue;
                        taxaList.add(taxa[k]);
                        pheValueList.add(value[k]);
                        loadList.add(l);
                    }
                    double[] pheno = new double[pheValueList.size()];
                    double[] load = new double[loadList.size()];
                    for (int k = 0; k < pheno.length; k++) {
                        pheno[k] = pheValueList.get(k);
                        load[k] = loadList.get(k);
                    }
                    String outfileS = new File (subDir[j], traits[i]).getAbsolutePath();
                    String pdfFileS = outfileS.replaceFirst(".txt", ".pdf");
                    BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                    bw.write("Taxa\tLoad\tPhenotype");
                    bw.newLine();
                    for (int k = 0; k < taxaList.size(); k++) {
                        StringBuilder sb = new StringBuilder();
                        sb.append(taxaList.get(k)).append("\t").append(load[k]).append("\t").append(pheno[k]);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    bw.flush();
                    bw.close();
                    double[][] matrix = new double[pheno.length][2];
                    for (int k = 0; k < pheno.length; k++) {
                        matrix[k][0] = pheno[k];
                        matrix[k][1] = load[k];
                    }
                    PearsonsCorrelation cor = new PearsonsCorrelation(matrix);
                    double r = cor.getCorrelationMatrix().getEntry(0, 1);
                    double pv = cor.getCorrelationPValues().getEntry(0, 1);
                    r2[i][j] = r;
                    ps[i][j] = pv;
                    ScatterPlot p = new ScatterPlot(load, pheno);
                    p.setXLab("Genetic load estimate from CDS");
                    p.setYLab("Phenotype");
                    p.setTitle(new File (pdfFileS).getName().replace(".pdf", "")+", r="+String.valueOf((float)r));
                    p.setColor(255, 0, 0, 30);
                    p.setPlottingCharacter(19);
                    p.saveGraph(pdfFileS);
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(loadPhenoCorrFileS);
            bw.write("Traits\tAllLoad\tCentromereLoad\tDistalLoad\tAllLoadPValue\tCentromereLoadPValue\tDistalLoadPValue");
            bw.newLine();
            for (int i = 0; i < traits.length; i++) {
                StringBuilder sb = new StringBuilder(traits[i]);
                for (int j = 0; j < r2[0].length; j++) {
                    sb.append("\t").append(r2[i][j]);
                }
                for (int j = 0; j < ps[0].length; j++) {
                    sb.append("\t").append(ps[i][j]);
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
    
    public void mkLoad () {
        String delSNPFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\class\\Non_Synonymous_Deleterious.txt";
        String NAMGenotypeDirS = "M:\\production\\maf\\traitAssociation\\genotype\\hmp321CDSSNP\\";
        String infoFileS = "E:\\Database\\InfoFile\\ChrLenCentPosi_agpV3.txt";
        String outfileS = "M:\\production\\maf\\traitAssociation\\nam\\load.txt";
        double mafCut = 1;
        int centRadius = 25_000_000;
        Table t = new Table (infoFileS);
        int chrNum = t.getRowNumber();
        int[][] centroPos = new int[chrNum][2];
        for (int i = 0; i < t.getRowNumber(); i++) {
            int mid = (t.getIntValue(i, 2)+t.getIntValue(i, 3))/2;
            centroPos[i][0] = mid-centRadius;
            centroPos[i][1] = mid+centRadius;
        }
        t = new Table (delSNPFileS);
        HashIntByteMap[] posMaMap = new HashIntByteMap[chrNum];
        HashIntDoubleMap[]  posMafMap = new HashIntDoubleMap[chrNum];
        for (int i = 0; i < chrNum; i++) {
            posMaMap[i] = HashIntByteMaps.getDefaultFactory().withDefaultValue((byte)-1).newMutableMap();
            posMafMap[i] = HashIntDoubleMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
        }
        HashByteByteMap ascIIByteMap = BaseCoder.getAscIIByteMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            int chrIndex = t.getIntValue(i, 0)-1;
            int pos = t.getIntValue(i, 1);
//            byte ma = ascIIByteMap.get(t.content[i][2].getBytes()[0]);
//            double maf = t.getDoubleValue(i, 3);
            if (t.content[i][4].length() > 1) continue;
            byte ma = ascIIByteMap.get(t.content[i][4].getBytes()[0]);
            double maf = t.getDoubleValue(i, 5);
            if (maf > mafCut) continue;
            posMaMap[chrIndex].put(pos, ma);
            posMafMap[chrIndex].put(pos, maf);
        }
        File[] fs = new File(NAMGenotypeDirS).listFiles();
        GenotypeTable gt = ImportUtils.readFromHapmap(fs[0].getAbsolutePath());
        TaxaList tl = gt.taxa();
        String[] taxa = new String[tl.size()];
        for (int i = 0; i < taxa.length; i++) {
            taxa[i] = tl.get(i).getName();
        }
        int[] allDel = new int[taxa.length];
        int[] centroDel = new int[taxa.length];
        int[] distalDel = new int[taxa.length];
        double[] allDelMaf = new double[taxa.length];
        double[] centroDelMaf = new double[taxa.length];
        double[] distalDelMaf = new double[taxa.length];
        int[] allNonmissing = new int[taxa.length];
        int[] centroNonmissing = new int[taxa.length];
        int[] distalNonmissing = new int[taxa.length];
        for (int i = 0; i < fs.length; i++) {
            gt = ImportUtils.readFromHapmap(fs[i].getAbsolutePath());
            int chrIndex = gt.chromosome(0).getChromosomeNumber()-1;
            for (int j = 0; j < gt.numberOfSites(); j++) {
                if (gt.alleles(j).length != 2) continue;
                int pos = gt.chromosomalPosition(j);
                byte ma = posMaMap[chrIndex].get(pos);
                if (ma == posMaMap[chrIndex].defaultValue()) continue;
                double maf = posMafMap[chrIndex].get(pos);
                if (pos > centroPos[chrIndex][0] && pos < centroPos[chrIndex][1]) {
                    for (int k = 0; k < taxa.length; k++) {
                        byte[] geno = gt.genotypeArray(k, j);
                        if (geno[0] == GenotypeTable.UNKNOWN_ALLELE || geno[1] == GenotypeTable.UNKNOWN_ALLELE) continue;
                        centroNonmissing[k]++;
                        allNonmissing[k]++;
                        if (geno[0] == ma || geno[1] == ma) {
                            double value = 1/maf;
                            centroDel[k]++;
                            allDel[k]++;
                            centroDelMaf[k]+=value;
                            allDelMaf[k]+=value;
                        }
                    }
                }
                else {
                    for (int k = 0; k < taxa.length; k++) {
                        byte[] geno = gt.genotypeArray(k, j);
                        if (geno[0] == GenotypeTable.UNKNOWN_ALLELE || geno[1] == GenotypeTable.UNKNOWN_ALLELE) continue;
                        distalNonmissing[k]++;
                        allNonmissing[k]++;
                        if (geno[0] == ma || geno[1] == ma) {
                            double value = 1/maf;
                            distalDel[k]++;
                            allDel[k]++;
                            distalDelMaf[k]+=value;
                            allDelMaf[k]+=value;
                        }
                    }
                }  
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Taxa\tAllDel\tCentroDel\tDistalDel\tAllNonMissing\tCentroNonMissing\tDistalNonMissing\tAllLoad\tCentroLoad\tDistalLoad\tAllLoadByMaf\tCentroLoadByMaf\tDistalLoadByMaf");
            bw.newLine();
            for (int i = 0; i < taxa.length; i++) {
                StringBuilder sb = new StringBuilder(taxa[i]);
                sb.append("\t").append(allDel[i]).append("\t").append(centroDel[i]).append("\t").append(distalDel[i]);
                sb.append("\t").append(allNonmissing[i]).append("\t").append(centroNonmissing[i]).append("\t").append(distalNonmissing[i]);
                sb.append("\t").append((double)allDel[i]/allNonmissing[i]).append("\t").append((double)centroDel[i]/centroNonmissing[i]).append("\t").append((double)distalDel[i]/distalNonmissing[i]);
                sb.append("\t").append((double)allDelMaf[i]/allNonmissing[i]).append("\t").append((double)centroDelMaf[i]/centroNonmissing[i]).append("\t").append((double)distalDelMaf[i]/distalNonmissing[i]);
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
    
    
    public void mkCDSGenotype () {
        String genotypeFileS = "/local/workdir/mingh/genotype/NAM_HM321_Unimp.hmp.txt.gz";
        String inDirS = "/workdir/mingh/hmp321SiftGerpUScore";
        String outDirS = "/workdir/mingh/hmp321CDSSNP";
        new File (outDirS).mkdir();
        File[] fs = new File(inDirS).listFiles();
        int[][] positions = new int[fs.length][];
        for (int i = 0; i < fs.length; i++) {
            Table t = new Table(fs[i].getAbsolutePath());
            int chrIndex = Integer.valueOf(t.content[0][0])-1;
            positions[chrIndex] = new int[t.getRowNumber()];
            for (int j = 0; j < t.getRowNumber(); j++) {
                positions[chrIndex][j] = t.getIntValue(j, 1);
            }
            Arrays.sort(positions[chrIndex]);
        }
        try {
            BufferedReader br = IoUtils.getTextGzipReader(genotypeFileS);
            String header = br.readLine();
            BufferedWriter[] bws = new BufferedWriter[positions.length];
            for (int i = 0; i < bws.length; i++) {
                String outfileS = "NAM_hm321_unimp_chr"+FStringUtils.getNDigitNumber(3, i+1)+".hmp.txt";
                outfileS = new File(outDirS, outfileS).getAbsolutePath();
                bws[i] = IoUtils.getTextWriter(outfileS);
                bws[i].write(header);
                bws[i].newLine();
            }
            String temp;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnt++;
                if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt)+" SNPs processed");
                List<String> l = FStringUtils.fastSplit(temp.substring(0, 100));
                int chrIndex = Integer.valueOf(l.get(2))-1;
                int pos = Integer.valueOf(l.get(3));
                int index = Arrays.binarySearch(positions[chrIndex], pos);
                if (index < 0) continue;
                bws[chrIndex].write(temp);
                bws[chrIndex].newLine();
            }
            for (int i = 0; i < bws.length; i++) {
                bws[i].flush();
                bws[i].close();
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
