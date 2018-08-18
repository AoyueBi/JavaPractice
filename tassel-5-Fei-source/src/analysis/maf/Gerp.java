/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maf;

import gnu.trove.list.array.TDoubleArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.TreeSet;
import utils.FArrayUtils;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class Gerp {
    
    public Gerp () {
        //this.gerpUscoreCorrelation();
        this.extractAncestralAllele();
        //this.gerpStatistics();
    }
    
    private void gerpStatistics () {
        String infiles = "Q:\\Zea\\Genotypes\\Annotations\\GERP\\dicot_monocot_andro2\\roast.chrom.10.msa.in.rates.full";
        String staFileS = "M:\\production\\maf\\annotations\\gerp\\gerpStatistics.txt";
        try {
            BufferedReader br = IoUtils.getTextReader(infiles);
            String temp = null;
            int total = 0;
            int cntTree = 0;
            double[] thresh = {0, 1, 2, 3, 4, 5};
            int[] cnt = new int[thresh.length];
            while ((temp = br.readLine()) != null) {
                total++;
                List<String> l = FStringUtils.fastSplit(temp, "\t");
                if (Double.valueOf(l.get(0)) == 0) {
                    cntTree++;
                }
                double v = Double.valueOf(l.get(1));
                for (int i = 0; i < thresh.length; i++) {
                    if (v > thresh[i]) cnt[i]++;
                }
            }
            br.close();
            BufferedWriter bw = IoUtils.getTextWriter(staFileS);
            bw.write("TreeLength != 0:\t"+ String.valueOf((double)cntTree/total));
            bw.newLine();
            for (int i = 0; i < thresh.length; i++) {
                bw.write("Score > "+String.valueOf(thresh[i])+" :\t" + String.valueOf((double)cnt[i]/total));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    private void extractAncestralAllele () {
        String alignmentDirS = "/workdir/mingh/gerpAlign/";
        String ancestralAlleleDirS = "/workdir/mingh/ancestralAllele";
        new File(ancestralAlleleDirS).mkdir();
        String[] taxa = {"Zea", "Coelorachis",	"Vossia", "Sorghum", "Setaria",  "Panicum", "Brachypodium", "Oryza"};
        File[] fs = new File(alignmentDirS).listFiles();
        Arrays.sort(fs);
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            int chr = Integer.valueOf(f.getName().replaceFirst("roast.chrom.", "").replaceFirst(".msa.in", ""));
            String outfileS = new File(ancestralAlleleDirS, "chr"+FStringUtils.getNDigitNumber(3, chr)+"ancestral.txt").getAbsolutePath();
            try {
                BufferedReader br = IoUtils.getTextReader(f.getAbsolutePath());
                ArrayList<String> nameList = new ArrayList();
                ArrayList<String> seqList = new ArrayList();
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    nameList.add(temp.replaceFirst(">", ""));
                    seqList.add(br.readLine());
                    System.out.println(temp);
                }
                String[] seq = new String[taxa.length];
                for (int i = 0; i < taxa.length; i++) {
                    for (int j = 0; j < nameList.size(); j++) {
                        if (taxa[i].equals(nameList.get(j))) {
                            seq[i] = seqList.get(j);
                        }
                    }
                }
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                bw.write("Position\tAncestralAllele");
                bw.newLine();
                for (int i = 0; i < seq[0].length(); i++) {
                    char[] chars = new char[taxa.length];
                    for (int j = 0; j < taxa.length; j++) {
                        chars[j] = seq[j].charAt(i);
                    }
                    TreeSet<Character> cSet = new TreeSet();
                    for (int j = 1; j < 4; j++) {
                        if (chars[j] == '-') continue;
                        cSet.add(chars[j]);
                    }
                    Character[] cs = cSet.toArray(new Character[cSet.size()]);
                    int[] cCount = new int[cs.length];
                    int sum = 0;
                    for (int j = 1; j < 4; j++) {
                        if (chars[j] == '-') continue;
                        int index = Arrays.binarySearch(cs, chars[j]);
                        cCount[index]++;
                        sum++;
                    }
                    if (sum == 3) {
                        if (cs.length == 1) {
                            bw.write(String.valueOf(i+1)+"\t"+String.valueOf(cs[0]));
                            bw.newLine();
                            continue;
                        }
                        else if (cs.length == 2) {
                            if (cCount[0]>cCount[1]) {
                                bw.write(String.valueOf(i+1)+"\t"+String.valueOf(cs[0]));
                                bw.newLine();
                                continue;
                            }
                            else if (cCount[0]<cCount[1]) {
                                bw.write(String.valueOf(i+1)+"\t"+String.valueOf(cs[1]));
                                bw.newLine();
                                continue;
                            }
                            else {}
                        }
                        else {
                            bw.write(String.valueOf(i+1)+"\tUnknown");
                            bw.newLine();
                            continue;
                        }
                    }
                    else if (sum == 2) {
                        if (cs.length == 1) {
                            bw.write(String.valueOf(i+1)+"\t"+String.valueOf(cs[0]));
                            bw.newLine();
                            continue;
                        }
                        else {}
                    }
                    else if (sum == 1){
                        if (cs[0] == chars[0]) {
                            bw.write(String.valueOf(i+1)+"\t"+String.valueOf(cs[0]));
                            bw.newLine();
                            continue;
                        }
                        else {}
                    }
                    for (int j = 4; j < chars.length; j++) {
                        if (chars[j] == '-') continue;
                        cSet.add(chars[j]);
                    }
                    if (cSet.size() == 0) continue;
                    cs = cSet.toArray(new Character[cSet.size()]);
                    cCount = new int[cs.length];
                    sum = 0;
                    for (int j = 1; j < chars.length; j++) {
                        if (chars[j] == '-') continue;
                        int index = Arrays.binarySearch(cs, chars[j]);
                        cCount[index]++;
                        sum++;
                    }
                    int[] index = FArrayUtils.getIndexByDescendingValue(cCount);
                    if (cs.length < 1) {}
                    else if (cs.length < 2) {
                        if (cs[0] == chars[0]) {
                            bw.write(String.valueOf(i+1)+"\t"+String.valueOf(cs[0]));
                            bw.newLine();
                            continue;
                        }
                        else {}
                    }
                    else {
                        if (cCount[index[0]]>cCount[index[1]]) {
                            bw.write(String.valueOf(i+1)+"\t"+String.valueOf(cs[index[0]]));
                            bw.newLine();
                            continue;
                        }
                        else {
                            bw.write(String.valueOf(i+1)+"\tUnknown");
                            bw.newLine();
                            continue;
                        }
                    }
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        
    }
    
    public void gerpUscoreCorrelation () {
        String infiles = "Q:\\Zea\\Genotypes\\Annotations\\GERP\\dicot_monocot_andro2p\\roast.chrom.10.msa.in.rates.full";
        String uniqueFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\chr010_uniqueness_6Genomes.txt";
        String outfileS = "M:\\production\\maf\\annotations\\gerp\\gerpAndUScore\\gerp_UScore.txt";
        try {
            BufferedReader br = IoUtils.getTextReader(infiles);
            String temp = null;
            TDoubleArrayList vList = new TDoubleArrayList();
            TDoubleArrayList tList = new TDoubleArrayList();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                tList.add(Double.valueOf(temp.split("\t")[0]));
                vList.add(Double.valueOf(temp.split("\t")[1]));
                cnt++;
                if (cnt%10000000 == 0) System.out.println(cnt);
            }
            br = IoUtils.getTextReader(uniqueFileS);
            temp = null;
            TDoubleArrayList uList = new TDoubleArrayList();
            br.readLine();
            cnt = 0;
            while ((temp = br.readLine()) != null) {
                double v = 0;
                if (temp.split("\t")[2].startsWith("N")) v = -1;
                else v = Double.valueOf(temp.split("\t")[2]);
                uList.add(v);
                cnt++;
                if (cnt%10000000 == 0) System.out.println(cnt);
            }
            double[] tree = tList.toArray();
            double[] gerp = vList.toArray();
            double[] uScore = uList.toArray();
            int size = 100000;
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("GerpScore\tTreeLength\tUScore");
            bw.newLine();
            for (int i = 0; i < size; i++) {
                int index = (int)(Math.random()*gerp.length);
                bw.write(String.valueOf(gerp[index])+"\t"+String.valueOf(tree[index])+"\t"+String.valueOf(uScore[index]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
                
    }
    
}
