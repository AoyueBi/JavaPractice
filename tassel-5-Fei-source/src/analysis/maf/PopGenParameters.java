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
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class PopGenParameters {
    
    public PopGenParameters () {
        //this.sliceVCFTest();
        //this.getGroups();
        //this.mkCommands();
        //this.fstTable();
        this.siteFst();
        //this.averagePi();
        //this.averagePiWindow();
        //this.avergePiFromGene();
    }
    
    
    void avergePiFromGene () {
        String geneFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
        String infileDirS = "M:\\production\\maf\\popgen\\paremeters\\pi\\piByGroups\\";
        String outfileS = "M:\\production\\maf\\popgen\\paremeters\\pi\\average_pi_gene.txt";
        GeneFeature gf = new GeneFeature(geneFileS);
        int startIndex = gf.getStartIndexOfChromosome(10);
        int endIndex = gf.getEndIndexOfChromosome(10);
        int length = 0;
        ArrayList<Range> rList = new ArrayList();
        for (int i = startIndex; i < endIndex; i++) {
            Range r = new Range (gf.getGeneChromosome(i), gf.getGeneStart(i), gf.getGeneEnd(i));
            rList.add(r);
            length += gf.getGeneEnd(i)-gf.getGeneStart(i);
        }
        Ranges gr = new Ranges(rList, "chr10");
        File[] fs = new File(infileDirS).listFiles();
        String[] groups = new String[fs.length];
        double[] pi = new double[fs.length];
        for (int i = 0; i < fs.length; i++) {
            TDoubleArrayList piList = new TDoubleArrayList();
            try {
                BufferedReader br = IoUtils.getTextGzipReader(fs[i].getAbsolutePath());
                groups[i] = fs[i].getName().replaceFirst(".sites.pi.gz", "");
                String temp = br.readLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    
                    String[] tem = temp.split("\t");
                    int pos = Integer.valueOf(tem[1]);
                    if (!gr.isInRanges(10, pos)) continue;
                    cnt++;
                    if (tem[2].startsWith("-n") || tem[2].startsWith("n")) continue;
                    piList.add(Double.valueOf(tem[2]));
                }
                br.close();
                DescriptiveStatistics d = new DescriptiveStatistics(piList.toArray());
                pi[i] = d.getMean()*cnt/length;
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Groups\tPi");
            bw.newLine();
            for (int i = 0; i < groups.length; i++) {
                bw.write(groups[i]+"\t"+String.valueOf(pi[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    void averagePiWindow () {
        String infileDirS = "M:\\production\\maf\\popgen\\paremeters\\pi\\piByGroups_window10000\\";
        String outfileS = "M:\\production\\maf\\popgen\\paremeters\\pi\\average_pi_window.txt";
        File[] fs = new File(infileDirS).listFiles();
        String[] groups = new String[fs.length];
        double[] pi = new double[fs.length];
        for (int i = 0; i < fs.length; i++) {
            TDoubleArrayList piList = new TDoubleArrayList();
            try {
                BufferedReader br = IoUtils.getTextGzipReader(fs[i].getAbsolutePath());
                groups[i] = fs[i].getName().replaceFirst(".window.sites.pi.gz", "");
                String temp = br.readLine();
                while ((temp = br.readLine()) != null) {
                    String[] tem = temp.split("\t");
                    if (tem[4].startsWith("-n") || tem[2].startsWith("n")) continue;
                    piList.add(Double.valueOf(tem[4]));
                }
                br.close();
                DescriptiveStatistics d = new DescriptiveStatistics(piList.toArray());
                pi[i] = d.getMean();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Groups\tPi");
            bw.newLine();
            for (int i = 0; i < groups.length; i++) {
                bw.write(groups[i]+"\t"+String.valueOf(pi[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    void averagePi () {
        String infileDirS = "M:\\production\\maf\\popgen\\paremeters\\pi\\piByGroups\\";
        String outfileS = "M:\\production\\maf\\popgen\\paremeters\\pi\\average_pi.txt";
        File[] fs = new File(infileDirS).listFiles();
        String[] groups = new String[fs.length];
        double[] pi = new double[fs.length];
        for (int i = 0; i < fs.length; i++) {
            TDoubleArrayList piList = new TDoubleArrayList();
            try {
                BufferedReader br = IoUtils.getTextGzipReader(fs[i].getAbsolutePath());
                groups[i] = fs[i].getName().replaceFirst(".sites.pi.gz", "");
                String temp = br.readLine();
                while ((temp = br.readLine()) != null) {
                    String[] tem = temp.split("\t");
                    if (tem[2].startsWith("-n") || tem[2].startsWith("n")) continue;
                    piList.add(Double.valueOf(tem[2]));
                }
                br.close();
                DescriptiveStatistics d = new DescriptiveStatistics(piList.toArray());
                pi[i] = d.getMean();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Groups\tPi");
            bw.newLine();
            for (int i = 0; i < groups.length; i++) {
                bw.write(groups[i]+"\t"+String.valueOf(pi[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    void siteFst () {
        String groupDirS = "M:\\production\\maf\\popgen\\paremeters\\fst\\groups\\";
        String fstDirS = "M:\\production\\maf\\popgen\\paremeters\\fst\\fst\\";
        String delSFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\class\\Non_Synonymous_Deleterious.txt";
        String delSGFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\class\\Non_Synonymous_Deleterious_High_GERP.txt";
        String nonSynFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\class\\Non_Synonymous_Tolerent.txt";
        String synFileS = "M:\\production\\maf\\annotations\\siftScore\\006_hmp321SNPClassMAF\\class\\Synonymous.txt";
        String delSDirS = "M:\\production\\maf\\popgen\\paremeters\\fst\\fstDeleteriousS\\";
        String delSGDirS = "M:\\production\\maf\\popgen\\paremeters\\fst\\fstDeleteriousSG\\";
        String nonSynDirS = "M:\\production\\maf\\popgen\\paremeters\\fst\\fstNonSynonymousTolerent\\";
        String synDirS = "M:\\production\\maf\\popgen\\paremeters\\fst\\fstSynonymous\\";
        int chrNum = 10;
        File[] fs = new File (groupDirS).listFiles();
        String[] groups = new String[fs.length];
        for (int i = 0; i < groups.length; i++) {
            groups[i] = fs[i].getName().replaceFirst(".txt", "");
        }
        Arrays.sort(groups);
//        this.outputSiteFst(chrNum, groups, delSFileS, fstDirS, delSDirS);
//        this.outputSiteFst(chrNum, groups, delSGFileS, fstDirS, delSGDirS);
//        this.outputSiteFst(chrNum, groups, synFileS, fstDirS, synDirS);
        this.outputSiteFst(chrNum, groups, nonSynFileS, fstDirS, nonSynDirS);
    }
    
    void outputSiteFst (int chrNum,String[] groups, String siteFileS, String fstDirS, String outputDirS) {
        TIntArrayList[] posLists = new TIntArrayList[chrNum];
        int[][] poss = new int[chrNum][];
        for (int i = 0; i < posLists.length; i++) {
            posLists[i] = new TIntArrayList();
        }
        Table t = new Table (siteFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            int chrIndex = Integer.valueOf(t.content[i][0])-1;
            int pos = Integer.valueOf(t.content[i][1]);
            posLists[chrIndex].add(pos);
        }
        
        for (int i = 0; i < posLists.length; i++) {
            poss[i] = posLists[i].toArray();
        }
        String[][][] fstS = new String[chrNum][groups.length*groups.length/2][];
        for (int i = 0; i < chrNum; i++) {
            for (int j = 0; j < fstS[0].length; j++) {
                fstS[i][j] = new String[poss[i].length];
            }
        }
        int cnt = 0;
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < groups.length-1; i++) {
            for (int j = i + 1; j < groups.length; j++) {
                for (int k = 0; k < chrNum; k++) {
                    String infileS = groups[i]+"VS"+groups[j]+"_chr"+FStringUtils.getNDigitNumber(3, k+1)+"_fst.txt.weir.fst.gz";
                    infileS = new File (fstDirS, infileS).getAbsolutePath();
                    if (!(new File (infileS).exists())) {
                        infileS = groups[j]+"VS"+groups[i]+"_chr"+FStringUtils.getNDigitNumber(3, k+1)+"_fst.txt.weir.fst.gz";
                        infileS = new File (fstDirS, infileS).getAbsolutePath();
                    }
                    try {
                        BufferedReader br = IoUtils.getTextGzipReader(infileS);
                        String temp = br.readLine();
                        while ((temp = br.readLine()) != null) {
                            List<String> l = FStringUtils.fastSplit(temp);
                            int pos = Integer.valueOf(l.get(1));
                            int index = Arrays.binarySearch(poss[k], pos);
                            if (index < 0) continue;
                            double v = Double.NaN;
                            if (!l.get(2).startsWith("-")) {
                                v = Double.valueOf(l.get(2));
                                if (v < 0) v = 0;
                            }
                            fstS[k][cnt][index] = String.valueOf(v);
                        }
                        br.close();
                    }
                    catch (Exception e) {
                        e.printStackTrace();
                    }
                }
                System.out.println(cnt);
                cnt++;
                sb.append("\t").append(groups[i]).append("_VS_").append(groups[j]);
            }
        }
        String outfileS = new File (outputDirS, "siteFst.txt").getAbsolutePath();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Chr\tPos"+sb.toString());
            bw.newLine();
            for (int i = 0; i < chrNum; i++) {
                for (int j = 0; j < poss[i].length; j++) {
                    sb = new StringBuilder();
                    sb.append(i+1).append("\t").append(poss[i][j]);
                    for (int k = 0; k < cnt; k++) {
                        sb.append("\t").append(fstS[i][k][j]);
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
    
    void fstTable () {
        String groupDirS = "M:\\production\\maf\\popgen\\paremeters\\fst\\groups\\";
        String fstDirS = "M:\\production\\maf\\popgen\\paremeters\\fst\\fst\\";
        String outfileS = "M:\\production\\maf\\popgen\\paremeters\\fst\\fstTable.txt";
        File[] fs = new File (groupDirS).listFiles();
        String[] groups = new String[fs.length];
        for (int i = 0; i < groups.length; i++) {
            groups[i] = fs[i].getName().replaceFirst(".txt", "");
        }
        Arrays.sort(groups);
        double[][] fstValues = new double[groups.length][groups.length];
        for (int i = 0; i < groups.length-1; i++) {
            for (int j = i + 1; j < groups.length; j++) {
                String infileS = groups[i]+"VS"+groups[j]+"_chr"+FStringUtils.getNDigitNumber(3, 10)+"_fst.txt.weir.fst.gz";
                infileS = new File (fstDirS, infileS).getAbsolutePath();
                if (!(new File (infileS).exists())) {
                    infileS = groups[j]+"VS"+groups[i]+"_chr"+FStringUtils.getNDigitNumber(3, 10)+"_fst.txt.weir.fst.gz";
                    infileS = new File (fstDirS, infileS).getAbsolutePath();
                }
                try {
                    BufferedReader br = IoUtils.getTextGzipReader(infileS);
                    String temp = br.readLine();
                    TDoubleArrayList vList = new TDoubleArrayList();
                    while ((temp = br.readLine()) != null) {
                        List<String> l = FStringUtils.fastSplit(temp);
                        String t = l.get(2);
                        if (t.startsWith("-n") || t.startsWith("n")) continue;
                        double v = Double.valueOf(t);
                        if (v < 0) v = 0;
                        vList.add(v);
                    }
                    br.close();
                    double[] v = vList.toArray();
                    DescriptiveStatistics d = new DescriptiveStatistics(v);
                    fstValues[i][j] = d.getMean();
                    System.out.println(infileS);
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder("Pupoluation");
            for (int i = 0; i < groups.length; i++) {
                sb.append("\t").append(groups[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < groups.length; i++) {
                sb = new StringBuilder(groups[i]);
                for (int j = 0; j < groups.length; j++) {
                    if (fstValues[i][j] == 0) {
                        if (i == j) {
                            sb.append("\tNA");
                        }
                        else {
                            sb.append("\t");
                        }
                    }
                    else {
                        sb.append("\t").append(fstValues[i][j]);
                    }
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
    
    void mkCommands () {
        String genotypeDirS = "/workdir/mingh/geno/";
        String groupDirS = "/workdir/mingh/groups/";
        String outputDirS = "/workdir/mingh/fst/";
        String perlScriptDirS = "/workdir/mingh/perl_fst/";
        File[] genotypeFiles = new File(genotypeDirS).listFiles();
        File[] groupFiles = new File(groupDirS).listFiles();
        new File(outputDirS).mkdir();
        new File (perlScriptDirS).mkdir();
        ArrayList<String> perlList = new ArrayList();
        for (int i = 0; i < groupFiles.length-1; i++) {
            String group1 = groupFiles[i].getName().replace(".txt", "");
            for (int j = i+1; j < groupFiles.length; j++) {
                String group2 = groupFiles[j].getName().replace(".txt", "");
                try {
                    String outfileS = new File (perlScriptDirS, group1+"VS"+group2+".pl").getAbsolutePath();
                    BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                    //vcftools --gzvcf test.vcf.gz --weir-fst-pop ../groups/Teosinte.txt --weir-fst-pop ../groups/Stiff_stalk.txt --weir-fst-pop ../groups/Non_stiff_stalk.txt --out out.txt
                    for (int k = 0; k < genotypeFiles.length; k++) {
                        int chr = Integer.valueOf(genotypeFiles[k].getName().replaceFirst("merged_flt_c", "").replaceFirst(".vcf.gz", ""));
                        String output = new File (outputDirS, group1+"VS"+group2+"_chr"+FStringUtils.getNDigitNumber(3, chr)+"_fst.txt").getAbsolutePath();
                        StringBuilder sb = new StringBuilder("vcftools --gzvcf ");
                        sb.append(genotypeFiles[k].getAbsolutePath()).append(" --weir-fst-pop ").append(groupFiles[i]).append(" --weir-fst-pop ").append(groupFiles[j]);
                        sb.append(" --out ").append(output);
                        bw.write("system (\""+sb.toString()+"\");");
                        bw.newLine();
                    }
                    bw.flush();
                    bw.close();
                    perlList.add(outfileS);
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        try {
            String outfileS = new File (perlScriptDirS, "perl_fst.pl").getAbsolutePath();
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            for (int i = 0; i < perlList.size(); i++) {
                StringBuilder sb = new StringBuilder("system (\"perl ");
                sb.append(perlList.get(i)).append(" &\");");
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
    
    void getGroups () {
        String infileS = "M:\\production\\maf\\popgen\\group\\geneticGroup.manual.txt";
        String outfileDirS = "M:\\production\\maf\\popgen\\paremeters\\fst\\groups\\";
        Table t = new Table (infileS);
        HashSet<String> groupSet = new HashSet();
        for (int i = 0; i < t.getRowNumber(); i++) {
            String group = t.content[i][4];
            groupSet.add(group);
        }
        String[] groups = groupSet.toArray(new String[groupSet.size()]);
        Arrays.sort(groups);
        ArrayList<String>[] taxaLists = new ArrayList[groups.length];
        for (int i = 0; i < taxaLists.length; i++) {
            taxaLists[i] = new ArrayList();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = Arrays.binarySearch(groups, t.content[i][4]);
            taxaLists[index].add(t.content[i][0]);
        }
        for (int i = 0; i < groups.length; i++) {
            String group = groups[i].replaceAll("-", "_").replaceAll(" ", "_");
            String outfileS = new File(outfileDirS, group+".txt").getAbsolutePath();
            try {
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                for (int j = 0; j < taxaLists[i].size(); j++) {
                    bw.write(taxaLists[i].get(j));
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
    
    void sliceVCFTest () {
        String infileS = "O:\\Zea\\Genotypes\\WGS\\HapMap\\v3\\v321\\unimp\\vcf\\merged_flt_c10.vcf.gz";
        String outfileS = "M:\\production\\maf\\popgen\\paremeters\\test\\test.vcf.gz";
        int n = 1000;
        try {
            BufferedReader br = IoUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IoUtils.getTextGzipWriter(outfileS);
            for (int i = 0; i < n; i++) {
                bw.write(br.readLine());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    void fst () {
        
    }
}


