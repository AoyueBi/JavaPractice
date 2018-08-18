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
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import graphcis.r.DensityPlot;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import utils.FArrayUtils;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class VariantDiscovery {
    
    public VariantDiscovery () {
        //this.extractHmp321Info();
        //this.updateHmp321InfoWithAncestralAllele();
        //this.subSet4Test();
        //this.densityTest();
        //this.filterHmp321Info();
        //this.countSite();
        //this.summarizeTranscript();
        //this.combineHmpSiftGerpUScore();
        //this.classifySNPs();
        //this.mkBarplotOfSNPs();
    }
    
    private void mkBarplotOfSNPs () {
        String infileDirS = "M:\\production\\maf\\001_variantDiscovery\\005_snpClass\\class\\";
        String countFileS = "M:\\production\\maf\\001_variantDiscovery\\005_snpClass\\classCount.txt";
        String mafDistrubutionFileS = "M:\\production\\maf\\001_variantDiscovery\\005_snpClass\\mafSFS.txt";
        String dafDistrubutionFileS = "M:\\production\\maf\\001_variantDiscovery\\005_snpClass\\dafSFS.txt";
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
    private void classifySNPs () {
        String geneFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
        String geneSummaryFileS = "M:\\production\\maf\\001_variantDiscovery\\003_transcriptSummary\\transcriptSummary.txt";
        String infileDirS = "M:\\production\\maf\\001_variantDiscovery\\004_hmp321SiftGerpUScore\\";
        String outfileDirS = "M:\\production\\maf\\001_variantDiscovery\\005_snpClass\\class\\";
        String transcriptFileS = "M:\\production\\maf\\001_variantDiscovery\\005_snpClass\\highConfidence_transcript.txt";
        String[] types = {"Synonymous", "Non_Synonymous_Tolerent", "StopLoss", "StopGain", "Non_Synonymous_Deleterious", "Non_Synonymous_Deleterious_High_GERP"};
        
        double gerpCut = 0;
        double nonSynRatioCut = 2.5;//filter incorrect gene model
        double nonSynGeneCut = 0.02;
        Table t = new Table (geneSummaryFileS);
        boolean[] ifOut = new boolean[t.getRowNumber()];
        ArrayList<String> tranNameList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (!t.content[i][4].equals("1")) continue;
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
                    double maf = t.getDoubleValue(j, 7);
                    String da = "NA";
                    String daf = "NA";
                    if (t.content[j][4].length() == 1) {
                        String an = t.content[j][4];
                        if (an.equals(t.content[j][5])) {
                            da = t.content[j][6];
                            daf = t.content[j][7];
                        }
                        else if (an.equals(t.content[j][6])) {
                            da = t.content[j][5];
                            daf = String.valueOf(1-Double.valueOf(t.content[j][7]));
                        }
                    }
                    if (t.content[j][14].equals("Syn")) {
                        if (t.content[j][15].startsWith("N")) continue;
                        double sift = t.getDoubleValue(j, 15);
                        bw[0].write(String.valueOf(chr)+"\t"+String.valueOf(pos)+"\t"+t.content[j][6]+"\t"+String.valueOf(maf)+"\t"+da+"\t"+daf);
                        bw[0].newLine();
                    }
                    if (t.content[j][14].equals("Non")) {
                        if (t.content[j][15].startsWith("N")) continue;
                        double sift = t.getDoubleValue(j, 15);
                        if (sift < 0.05) {
                            if (t.getDoubleValue(j, 17) > gerpCut) {
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
                    if (t.content[j][14].equals("StopLoss")) {
                        bw[2].write(String.valueOf(chr)+"\t"+String.valueOf(pos)+"\t"+t.content[j][6]+"\t"+String.valueOf(maf)+"\t"+da+"\t"+daf);
                        bw[2].newLine();
                    }
                    if (t.content[j][14].equals("StopGain")) {
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
     * Combine SIFT GERP into one table in coding variants of hmp321
     * Note: sift annotation using a different version of gene annotation, containing RNA genes, GFF3 do not have RNA, I use GFF3
     */
    private void combineHmpSiftGerpUScore () {
        String geneFileS = "/workdir/mingh/Zea_mays.AGPv3.26.gf.txt";
        String hmpInfoDirS = "/workdir/mingh/hmp321Info_filter/";
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
                bw.write("Chr\tPos\tRef\tAlt\tAncestral(Gerp)\tMajor\tMinor\tMinorAlleleFrequency\tMinorPresence\tMinorTotalDepth\tMaxMinorDepth\tSiteDepth\tSiteCount\tHetCount\tMutationClass\tSift\tGerpTreeLength\tGerp\tUScore\tTranscripts");
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
    
    private void summarizeTranscript () {
        String geneFeatureFileS = "/workdir/mingh/Zea_mays.AGPv3.26.gf.txt";
        String hmpInfoDirS = "/workdir/mingh/hmp321Info_filter/";
        String gerpDirS = "/workdir/mingh/gerp/";
        String siftDirS = "/workdir/mingh/sift/";
        String outfileS = "/workdir/mingh/transcriptSummary.txt";
        double gerpCut = 0;
        File[] fs = new File(hmpInfoDirS).listFiles();
        int chrNum = fs.length;
        HashMap<Integer, ArrayList<String>>[] posGeneMap = new HashMap[chrNum];
        int[][] snpPos = new int[chrNum][];
        byte[][] snps = new byte[chrNum][];
        byte[][] snpAnc = new byte[chrNum][];
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
                    ArrayList<String> geneNameList = posGeneMap[chrIndex].get(k);
                    if (geneNameList == null) {
                        geneNameList = new ArrayList();
                        geneNameList.add(geneName);
                        posGeneMap[chrIndex].put(k, geneNameList);
                    }
                    else {
                        geneNameList.add(geneName);
                        posGeneMap[chrIndex].put(k, geneNameList);
                    }
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
            TByteArrayList snpAncList = new TByteArrayList();
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
                    ArrayList<String> geneNameList = posGeneMap[chrIndex].get(pos);
                    if (geneNameList == null) continue;
                    for (int i = 0; i < geneNameList.size(); i++) {
                        int index = Arrays.binarySearch(genes, geneNameList.get(i));
                        snpCount[index]++;
                    }
                    snpPosList.add(pos);
                    snpList.add(l.get(3).getBytes()[0]);
                    snpAncList.add(l.get(4).getBytes()[0]);
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            snpPos[chrIndex] = snpPosList.toArray();
            snps[chrIndex] = snpList.toArray();
            snpAnc[chrIndex] = snpAncList.toArray();
        });
        int[] synCount = new int[genes.length];
        int[] nonCount = new int[genes.length];
        int[] delCount = new int[genes.length];
        int[] delHGCount = new int[genes.length];
        int[] naCount = new int[genes.length];
        
        int[] b73SynCount = new int[genes.length];
        int[] b73NonCount = new int[genes.length];
        int[] b73DelCount = new int[genes.length];
        int[] b73DelHGCount = new int[genes.length];
        int[] noAncCount = new int[genes.length];
        
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
                    byte ref = l.get(1).getBytes()[0];
                    byte derivedState = -1; //mean ancestral allele is not defined
                    if (snpAnc[chrIndex][index] == snps[chrIndex][index]) {
                        derivedState = 1; //mean b73 carries derived allele
                    }
                    else if (snpAnc[chrIndex][index] == ref) {
                        derivedState = 0;
                    }
                    
                    if (derivedState == -1) noAncCount[geneIndex]++;
                    
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
                                if (derivedState == 1) {
                                    b73SynCount[geneIndex]++;
                                }
                            }
                            else {
                                type = "Non";
                                nonCount[geneIndex]++;
                                if (derivedState == 1) {
                                    b73NonCount[geneIndex]++;
                                }
                                if (l.get(10).startsWith("N")) {
                                    naCount[geneIndex]++;
                                }
                                else {
                                    if (Double.valueOf(l.get(10)) < 0.05) {
                                        delCount[geneIndex]++;
                                        delPosList[chrIndex].add(pos);
                                        if (derivedState == 1) {
                                            b73DelCount[geneIndex]++;
                                        }
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
                    ArrayList<String> geneNameList = posGeneMap[chrIndex].get(cnt);
                    if (geneNameList == null) continue;
                    for (int i = 0; i < geneNameList.size(); i++) {
                        String gene = geneNameList.get(i);
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

                        byte derivedState = 0; //mean ancestral allele is not defined or ancestral allele is alt
                        if (snpAnc[chrIndex][index] == snps[chrIndex][index]) {
                            derivedState = 1; //mean b73 carries derived allele
                        }
                        snpGerpAlignCount[geneIndex]++;
                        snpGerpTree[geneIndex]+=treeValue;
                        snpGerpScore[geneIndex]+=scoreValue;
                        index = Arrays.binarySearch(delPos[chrIndex], cnt);
                        if (index < 0) continue;
                        if (scoreValue <= gerpCut) continue;
                        delHGCount[geneIndex]++;
                        if (derivedState == 1) {
                            b73DelHGCount[geneIndex]++;
                        }
                    }
                    
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        
        try {
            String header = "Transcript\tCDSLength\tSNPNumber\tSNPPercentage\tIfSiftAligned\tNumberOfSyn\tPercentageSyn\tNumberOfNon\tPercentageNon\tNonVsSynRatio\tNumberOfDeleterious\tPercentageDeleterious\tNumberOfHGDeleterious\tPercentageHGDeleterious\tIfGerpAligned\tGerpAlignedCount\tPercentageGerpAlignedCount\tMeanGerpTreeLength\tMeanGerpScore\tSNPMeanGerpTreeLength\tSNPMeanGerpScore";
            header = header +"\tNumAmbigousAnc\tB73NumberOfSyn\tB73PercentageSyn\tB73NumberOfNon\tB73PercentageNon\tB73NumberOfDeleterious\tB73PercentageDeleterious\tB73NumberOfHGDeleterious\tB73PercentageHGDeleterious";
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
                
                double cdsL = (double)(snpCount[i]-noAncCount[i])/snpCount[i]*cdsLength;
                sb.append("\t").append(noAncCount[i]).append("\t").append(b73SynCount[i]).append("\t").append((double)b73SynCount[i]/cdsL).append("\t");
                sb.append(b73NonCount[i]).append("\t").append((double)b73NonCount[i]/cdsL).append("\t");
              
                sb.append(b73DelCount[i]).append("\t").append((double)b73DelCount[i]/cdsL).append("\t");
                sb.append(b73DelHGCount[i]).append("\t").append((double)b73DelHGCount[i]/cdsL);
                
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
    
    private void countSite () {
        String inDirS = "M:\\production\\maf\\001_variantDiscovery\\002_hmp321Info_filter\\";
        File[] fs = new File(inDirS).listFiles();
        int sum = 0;
        for (int i = 0; i < fs.length; i++) {
            int cnt = 0;
            try {
                BufferedReader br = IoUtils.getTextGzipReader(fs[i].getAbsolutePath());
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println(String.valueOf(cnt)+"\t"+fs[i].getName());
            sum+=cnt;
        }
        System.out.println(String.valueOf(sum));
       
    }
    
    private void filterHmp321Info () {
        String inDirS = "M:\\production\\maf\\001_variantDiscovery\\001_hmp321Info\\";
        String outDirS = "M:\\production\\maf\\001_variantDiscovery\\002_hmp321Info_filter\\";
        int minMinorDepth = 3;
        double maxIndiDepth = 7;
        int siteDepthCut = 5000;
        double siteHeterozygousCut = 0.15;
        File[] fs = new File(inDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String outfileS = new File(outDirS, f.getName().replaceFirst("hmp321Info", "hmp321Info_filter")).getAbsolutePath();
            try {
                BufferedReader br = IoUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedWriter bw = IoUtils.getTextGzipWriter(outfileS);
                bw.write(br.readLine());
                bw.newLine();
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    List<String> l = FStringUtils.fastSplit(temp);
                    if (l.get(3).contains(",")) continue;
                    if (Integer.valueOf(l.get(10)) < minMinorDepth) continue;
                    if (Integer.valueOf(l.get(11)) > siteDepthCut) continue;
                    double siteHeterozygous = Double.valueOf(l.get(13))/Double.valueOf(l.get(12));
                    if (siteHeterozygous > siteHeterozygousCut) continue;
                    double indiDepth = Double.valueOf(l.get(11))/Double.valueOf(l.get(12));
                    if (indiDepth>maxIndiDepth) continue;
                    bw.write(temp);
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
    
    private void densityTest () {
        String infileS = "M:\\production\\maf\\001_variantDiscovery\\testFilter\\test.txt";
        Table t = new Table (infileS);
        double[] depth = t.getDoubleArrayByColumn(11);
        double[] siteCount = t.getDoubleArrayByColumn(12);
        double[] d = new double[depth.length];
        for (int i = 0; i < d.length; i++) {
            d[i] = depth[i]/siteCount[i];
        }
        DensityPlot dd = new DensityPlot(d);
        dd.showGraph();
    }
    
    private void subSet4Test () {
        String infileS = "M:\\production\\maf\\001_variantDiscovery\\001_hmp321Info\\hmp321Info_chr010.txt.gz";
        String outfileS = "M:\\production\\maf\\001_variantDiscovery\\testFilter\\test.txt";
        try {
            BufferedReader br = IoUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write(br.readLine());
            bw.newLine();
            String temp = null;
            while ((temp = br.readLine()) != null) {
                double r = Math.random();
                if (r > 0.001) continue;
                List<String> l = FStringUtils.fastSplit(temp);
                if (l.get(3).contains(",")) continue;
                bw.write(temp);
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
        String inDirS = "/workdir/mingh/hmp321Info/";
        String ancestralDirS = "/workdir/mingh/ancestralAllele/";
        String outDirS = "/workdir/mingh/hmp321InfoAnce/";
        File[] fs = new File(inDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            int chr = Integer.valueOf(f.getName().replaceFirst("hmp321Info_chr", "").replaceFirst(".txt.gz", ""));
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
                bw.write("Chr\tPos\tRef\tAlt\tAncestral(Gerp)\tMajor\tMinor\tMinorAlleleFrequency\tMinorPresence\tMinorTotalDepth\tMaxMinorDepth\tSiteDepth\tSiteCount\tHetCount");
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
    
    private void extractHmp321Info () {
        String infileDirS = "/workdir/mingh/hmp";
        String outfileDirS = "/workdir/mingh/hmp321Info/";
//        String infileDirS = "M:\\production\\maf\\annotations\\siftScore\\test\\hmp";
//        String outfileDirS = "M:\\production\\maf\\annotations\\siftScore\\test\\hmpInfo";
        new File(outfileDirS).mkdir();
        File[] fs = new File(infileDirS).listFiles();
        List<File> fl = Arrays.asList(fs);
        fl.parallelStream().forEach(f -> {
            int chr = Integer.valueOf(f.getName().split("_c")[1].replaceFirst(".vcf.gz", ""));
            String outfileS = new File(outfileDirS, "hmp321Info_chr"+FStringUtils.getNDigitNumber(3, chr)+".txt.gz").getAbsolutePath();
            String temp = null;
            System.out.println("Processing " + f.getName());
            try {
                BufferedReader br = IoUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedWriter bw = IoUtils.getTextGzipWriter(outfileS);
                bw.write("Chr\tPos\tRef\tAlt\tAncestral(Tripsacum)\tMajor\tMinor\tMinorAlleleFrequency\tMinorPresence\tMinorTotalDepth\tMaxMinorDepth\tSiteDepth\tSiteCount\tHetCount");
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
        int siteCount = Integer.valueOf(al.get(1).split("NZ=")[1]);
        String[] tem = al.get(3).replaceFirst("AN=", "").split(",");
        for (int i = 0; i < alleleCount.length; i++) {
            int cnt = 0;
            if (!tem[i].isEmpty()) cnt =Integer.valueOf(tem[i]);
            alleleCount[i] = Integer.valueOf(cnt);
        }
        int[] index = FArrayUtils.getIndexByDescendingValue(alleleCount);
        sb.append(alleles[index[0]]).append("\t");
        StringBuilder minorSB = new StringBuilder();
        for (int i = 1; i < alleles.length; i++) {
            minorSB.append(alleles[index[i]]).append(",");
        }
        minorSB.deleteCharAt(minorSB.length()-1).append("\t");
        int sum = 0;
        for (int i = 0; i < alleleCount.length; i++) {
            sum+=alleleCount[i];
        }
        double[] fre = new double[alleleCount.length];
        for (int i = 0; i < fre.length; i++) fre[i] = (double)alleleCount[i]/sum;
        for (int i = 1; i < alleles.length; i++) {
            minorSB.append((float)fre[index[i]]).append(",");
        }
        minorSB.deleteCharAt(minorSB.length()-1).append("\t");
        for (int i = 1; i < alleles.length; i++) {
            minorSB.append(alleleCount[index[i]]).append(",");
        }
        minorSB.deleteCharAt(minorSB.length()-1).append("\t");
        tem = al.get(2).replaceFirst("AD=", "").split(",");
        for (int i = 1; i < alleles.length; i++) {
            minorSB.append(tem[index[i]]).append(",");
        }
        minorSB.deleteCharAt(minorSB.length()-1).append("\t");
        
        int[] maxDepth = new int[alleleNum];
        for (int i = 0; i < maxDepth.length; i++) {
            maxDepth[i] = Integer.MIN_VALUE;
        }
        int hetCount = 0;
        for (int i = 9; i < l.size(); i++) {
            if (l.get(i).startsWith(".")) continue;
            List<String> ll = FStringUtils.fastSplit(l.get(i), ":");
            List<String> lll0 = FStringUtils.fastSplit(ll.get(0), "/");
            List<String> lll1 = FStringUtils.fastSplit(ll.get(1), ",");
            for (int j = 0; j < maxDepth.length; j++) {
                if (Integer.valueOf(lll0.get(0)) == j){
                    int cnt = Integer.valueOf(lll1.get(j));
                    if (cnt > maxDepth[j]) maxDepth[j] = cnt;
                }
                else if (Integer.valueOf(lll0.get(1)) == j) {
                    int cnt = Integer.valueOf(lll1.get(j));
                    if (cnt > maxDepth[j]) maxDepth[j] = cnt;
                }
            }
            if (lll0.get(0).equals(lll0.get(1))) continue;
            hetCount++;
        }
        for (int i = 1; i < alleles.length; i++) {
            minorSB.append(maxDepth[index[i]]).append(",");
        }
        minorSB.deleteCharAt(minorSB.length()-1);
        sb.append(minorSB.toString());
        sb.append("\t").append(siteDepth).append("\t").append(siteCount).append("\t").append(hetCount);
        return sb.toString();
    }
}
