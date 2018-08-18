/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.cassava.wgs;

import format.Table;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
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
class BurdenRevisit {
    
    public BurdenRevisit () {
        //this.getDeleSNPList();
        //this.extractHmp();
        this.countByGroup();
    }
    
    public void countByGroup () {
        String vcfFileS = "E:\\Research\\cassava\\revision\\burden\\dele.vcf.txt";
        String deleFileS = "E:\\Research\\cassava\\revision\\burden\\deleList.txt";
        String infoFileS = "E:\\Research\\cassava\\revision\\burden\\source\\cloneInfo.txt";
        String outfileS = "E:\\Research\\cassava\\revision\\burden\\total.txt";
        Table t = new Table (infoFileS);
        ArrayList<String> lacList = new ArrayList();
        ArrayList<String> acList = new ArrayList();
        ArrayList<String> proList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.content[i][12].startsWith("PRO")) {
                proList.add(t.content[i][1]);
            }
            else if (t.content[i][12].startsWith("LAC") && !t.content[i][13].startsWith("Ref")) {
                lacList.add(t.content[i][1]);
            } 
            else if (t.content[i][12].startsWith("Africa") && !t.content[i][13].startsWith("Intro")) {
                acList.add(t.content[i][1]);
            } 
        }
        String[] lacTaxa = lacList.toArray(new String[lacList.size()]);
        Arrays.sort(lacTaxa);
        String[] acTaxa = acList.toArray(new String[acList.size()]);
        Arrays.sort(acTaxa);
        String[] proTaxa = proList.toArray(new String[proList.size()]);
        Arrays.sort(proTaxa);
        t = new Table (deleFileS);
        HashMap<Integer, String> posAncMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            posAncMap.put(t.getIntValue(i, 0)*100_000_000+t.getIntValue(i, 1), t.content[i][6]);
        }
        int[] lacCnt = new int[lacTaxa.length];
        int[] acCnt = new int[acTaxa.length];
        int[] proCnt = new int[proTaxa.length];
        int[] lacIndices = new int[lacTaxa.length];
        int[] acIndices = new int[acTaxa.length];
        int[] proIndices = new int[proTaxa.length];
        try {
            BufferedReader br = IoUtils.getTextReader(vcfFileS);
            String temp = br.readLine();
            String[] tem = temp.split("\t");
            for (int i = 9; i < tem.length; i++) {
                int index = Arrays.binarySearch(lacTaxa, tem[i]);
                if (index < 0) continue;
                lacIndices[index] = i;
            }
            for (int i = 9; i < tem.length; i++) {
                int index = Arrays.binarySearch(acTaxa, tem[i]);
                if (index < 0) continue;
                acIndices[index] = i;
            }
            for (int i = 9; i < tem.length; i++) {
                int index = Arrays.binarySearch(proTaxa, tem[i]);
                if (index < 0) continue;
                proIndices[index] = i;
            }
            while ((temp = br.readLine()) != null) {
                tem = temp.split("\t");
                int pos = Integer.valueOf(tem[0])*100_000_000+Integer.valueOf(tem[1]);
                String anc =  posAncMap.get(pos);
                int v = 0;
                if (anc.equals(tem[3])) v = 1;
//homo                
//                for (int i = 0; i < lacTaxa.length; i++) {
//                    String[] te = tem[lacIndices[i]].split(":")[0].split("\\|");
//                    if (Integer.valueOf(te[0]) == v && Integer.valueOf(te[1]) == v) lacCnt[i]++;
//                }
//                for (int i = 0; i < acTaxa.length; i++) {
//                    String[] te = tem[acIndices[i]].split(":")[0].split("\\|");
//                    try {
//                        if (Integer.valueOf(te[0]) == v && Integer.valueOf(te[1]) == v) acCnt[i]++;
//                    }
//                    catch (Exception e) {
//                        System.out.println(te[0]);
//                        e.printStackTrace();
//                    }
//                }
//                for (int i = 0; i < proTaxa.length; i++) {
//                    String[] te = tem[proIndices[i]].split(":")[0].split("\\|");
//                    if (Integer.valueOf(te[0]) == v && Integer.valueOf(te[1]) == v) proCnt[i]++;
//                }
//hetero                
//                for (int i = 0; i < lacTaxa.length; i++) {
//                    String[] te = tem[lacIndices[i]].split(":")[0].split("\\|");
//                    if (Integer.valueOf(te[0]) != v && Integer.valueOf(te[1]) == v) lacCnt[i]++;
//                }
//                for (int i = 0; i < acTaxa.length; i++) {
//                    String[] te = tem[acIndices[i]].split(":")[0].split("\\|");
//                    if (Integer.valueOf(te[0]) != v && Integer.valueOf(te[1]) == v) acCnt[i]++;
//                }
//                for (int i = 0; i < proTaxa.length; i++) {
//                    String[] te = tem[proIndices[i]].split(":")[0].split("\\|");
//                    if (Integer.valueOf(te[0]) != v && Integer.valueOf(te[1]) == v) proCnt[i]++;
//                }
//total                
                for (int i = 0; i < lacTaxa.length; i++) {
                    String[] te = tem[lacIndices[i]].split(":")[0].split("\\|");
                    for (int j = 0; j < te.length; j++) {
                        if (Integer.valueOf(te[j]) == v) lacCnt[i]++;
                    }
                }
                for (int i = 0; i < acTaxa.length; i++) {
                    String[] te = tem[acIndices[i]].split(":")[0].split("\\|");
                    for (int j = 0; j < te.length; j++) {
                        if (Integer.valueOf(te[j]) == v) acCnt[i]++;
                    }
                }
                for (int i = 0; i < proTaxa.length; i++) {
                    String[] te = tem[proIndices[i]].split(":")[0].split("\\|");
                    for (int j = 0; j < te.length; j++) {
                        if (Integer.valueOf(te[j]) == v) proCnt[i]++;
                    }
                }
            }
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("LAC");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            StringBuilder sbb = new StringBuilder();
            for (int i = 0; i < lacCnt.length; i++) {
                sbb.append(lacTaxa[i]).append("\t");
                sb.append(lacCnt[i]).append("\t");
            }
            sbb.deleteCharAt(sbb.length()-1);
            bw.write(sbb.toString());
            bw.newLine();
            sb.deleteCharAt(sb.length()-1);
            bw.write(sb.toString());
            bw.newLine();
            bw.write("AC");
            bw.newLine();
            sb = new StringBuilder();
            sbb = new StringBuilder();
            for (int i = 0; i < acCnt.length; i++) {
                sbb.append(acTaxa[i]).append("\t");
                sb.append(acCnt[i]).append("\t");
            }
            sbb.deleteCharAt(sbb.length()-1);
            bw.write(sbb.toString());
            bw.newLine();
            sb.deleteCharAt(sb.length()-1);
            bw.write(sb.toString());
            bw.newLine();
            bw.write("Pro");
            bw.newLine();
            sb = new StringBuilder();
            sbb = new StringBuilder();
            for (int i = 0; i < proCnt.length; i++) {
                sbb.append(proTaxa[i]).append("\t");
                sb.append(proCnt[i]).append("\t");
            }
            sbb.deleteCharAt(sbb.length()-1);
            bw.write(sbb.toString());
            bw.newLine();
            sb.deleteCharAt(sb.length()-1);
            bw.write(sb.toString());
            bw.newLine();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void extractHmp () {
        String vcf = "E:\\Research\\cassava\\revision\\burden\\source\\new.txt";
        String deleFileS = "E:\\Research\\cassava\\revision\\burden\\deleList.txt";
        String infoFileS = "E:\\Research\\cassava\\revision\\burden\\source\\cloneInfo.txt";
        String outfileS = "E:\\Research\\cassava\\revision\\burden\\dele.vcf.txt";
        
//        String vcf = "/workdir/mingh/vcf/globalCassava_Phased.vcf";
//        String deleFileS = "/workdir/mingh/deleList.txt";
//        String infoFileS = "/workdir/mingh/cloneInfo.txt";
//        String outfileS = "/workdir/mingh/dele.vcf.txt";
        Table t = new Table (deleFileS);
        int[] tpos = new int[t.getRowNumber()];
        for (int i = 0; i < tpos.length; i++) {
            tpos[i] = t.getIntValue(i, 0)*100_000_000+t.getIntValue(i, 1);
        }
        ArrayList<String> taxaList = new ArrayList();
        t = new Table (infoFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.content[i][13].startsWith("Culti") || t.content[i][13].startsWith("PRO") || t.content[i][13].startsWith("Ref")) taxaList.add(t.content[i][1]);
        }
        String[] taxa = taxaList.toArray(new String[taxaList.size()]);
        Arrays.sort(taxa);
        try {
            BufferedReader br = IoUtils.getTextReader(vcf);
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            String[] tem = br.readLine().split("\t");
            int[] indices = new int[taxa.length];
            for (int i = 0; i < indices.length; i++) {
                indices[i] = -1;
            }
            for (int i = 9; i < tem.length; i++) {
                int index = Arrays.binarySearch(taxa, tem[i]);
                if (index < 0)  {
                    continue;
                }
                indices[index] = i;
            }
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < 9; i++) {
                sb.append(tem[i]).append("\t");
            }
            for (int i = 0; i < indices.length; i++) {
                sb.append(tem[indices[i]]).append("\t");
            }
            sb.deleteCharAt(sb.length()-1);
            bw.write(sb.toString());
            bw.newLine();
            String temp;
            while ((temp = br.readLine()) != null) {
                tem = temp.substring(0, 50).split("\t");
                int pos = Integer.valueOf(tem[0])*100_000_000+Integer.valueOf(tem[1]);
                int index = Arrays.binarySearch(tpos, pos);
                if (index < 0) continue;
                if (tem[4].length() > 1) continue;
                tem = temp.split("\t");
                sb = new StringBuilder();
                for (int i = 0; i < 9; i++) {
                    sb.append(tem[i]).append("\t");
                }
                HashSet<String> s = new HashSet();
                for (int i = 0; i < indices.length; i++) {
                    sb.append(tem[indices[i]]).append("\t");
                    String[] te = tem[indices[i]].split(":")[0].split("\\|");
                    s.add(te[0]);
                    s.add(te[1]);
                }
                sb.deleteCharAt(sb.length()-1);
                if (s.size() == 1) continue;
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
    
    public void getDeleSNPList () {
        //String siftFileS = "/workdir/mingh/sift/cassavaHMPfiltImp_SIFTannotations.txt";
        //String gerpscoreDirS = "/workdir/mingh/gerp/score/";
        //String gerpalignDirS = "/workdir/mingh/gerp/score/alignment";
        String siftFileS = "E:\\Research\\cassava\\revision\\burden\\source\\cassavaHMPfiltImp_SIFTannotations.txt";
        String gerpscoreDirS = "E:\\Research\\cassava\\revision\\burden\\source\\gerp\\score\\";
        String gerpalignDirS = "E:\\Research\\cassava\\revision\\burden\\source\\gerp\\alignment\\";
        String outfileS = "E:\\Research\\cassava\\revision\\burden\\deleList.txt";
        HashSet<DeleRecord> deleSet = new HashSet();
        try {
            BufferedReader br = IoUtils.getTextReader(siftFileS);
            String temp = br.readLine();
            while ((temp = br.readLine()) != null) {
                List<String> l  = FStringUtils.fastSplit(temp);
                if (!l.get(7).startsWith("CDS")) continue;
                if (!l.get(8).startsWith("NONSYNONYMOUS")) continue;
                if (l.get(12).startsWith("NA")) continue;
                double sift = Double.valueOf(l.get(12));
                if (sift >= 0.05) continue;
                DeleRecord d = new DeleRecord(Integer.valueOf(l.get(0)), Integer.valueOf(l.get(1)), sift, l.get(2), l.get(3));
                deleSet.add(d);
            }
            br.close();
            File[] fs = new File(gerpscoreDirS).listFiles();
            DeleRecord[] deleArray = deleSet.toArray(new DeleRecord[deleSet.size()]);
            Arrays.sort(deleArray);
            List<File> fList = Arrays.asList(fs);
            fList.parallelStream().forEach(f -> {
                int chr  = Integer.valueOf(f.getName().replaceFirst("roast.chrom.", "").replaceFirst(".msa.in.rates.full", ""));
                int cnt = 0;
                String tem = null;
                try {
                    BufferedReader brr = IoUtils.getTextReader(f.getAbsolutePath());
                    while ((tem = brr.readLine()) != null) {
                        cnt++;
                        List<String> l  = FStringUtils.fastSplit(tem);
                        double gerp = Double.valueOf(l.get(1));
                        DeleRecord query = new DeleRecord(chr, cnt);
                        if (!deleSet.contains(query)) continue;
                        int index = Arrays.binarySearch(deleArray, query);
                        deleArray[index].addGerp(gerp);
                    }
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
                System.out.println(String.valueOf(cnt)+"\t"+f.getName());
            });
            fs = new File(gerpalignDirS).listFiles();
            fList = Arrays.asList(fs);
            fList.parallelStream().forEach(f -> {
                int chr  = Integer.valueOf(f.getName().replaceFirst("roast.chrom.", "").replaceFirst(".msa.in", ""));
                try {
                    BufferedReader brr = IoUtils.getTextReader(f.getAbsolutePath());
                    brr = IoUtils.getTextReader(f.getAbsolutePath());
                    brr.readLine();
                    String ref = brr.readLine();
                    brr.readLine();
                    String rubber = brr.readLine();
                    for (int j = 0; j < ref.length(); j++) {
                        String re = ref.substring(j, j+1);
                        String ru = rubber.substring(j, j+1);
                        DeleRecord query = new DeleRecord(chr, j+1);
                        if (deleSet.contains(query)) {
                            int index = Arrays.binarySearch(deleArray, query);
                            deleArray[index].addDerived(ru);
                        }
                    }
                    System.out.println(ref.length()+"\t"+f.getName());
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            });
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Chr\tPos\tSift\tGerp\tRef\tAlt\tDerived");
            bw.newLine();
            for (int i = 0; i < deleArray.length; i++) {
                if (deleArray[i].gerp <= 2) continue;
                if (deleArray[i].ref.equals("N") || deleArray[i].ref.equals("-")) continue;
                if (deleArray[i].derived.equals("N") || deleArray[i].derived.equals("-")) continue;
                if (!deleArray[i].derived.equals(deleArray[i].ref) && !deleArray[i].derived.equals(deleArray[i].alt)) continue;
                bw.write(deleArray[i].getOutputString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    class DeleRecord implements Comparable<DeleRecord> {
        int chr;
        int pos;
        double sift;
        double gerp;
        String ref = null;
        String alt = null;
        String derived = null;
        
        public DeleRecord (int chr, int pos, double sift, String ref, String alt) {
            this.chr = chr;
            this.pos = pos;
            this.sift = sift;
            this.ref = ref;
            this.alt = alt;
        }
        
        public DeleRecord (int chr, int pos) {
            this.chr = chr;
            this.pos = pos;
        }
        
        public void addGerp (double gerp) {
            this.gerp = gerp;
        }
        
        public void addDerived (String derived) {
            this.derived = derived;
        }
        
        public String getOutputString () {
            StringBuilder sb = new StringBuilder();
            sb.append(chr).append("\t").append(pos).append("\t").append(sift).append("\t").append(gerp);
            sb.append("\t").append(ref).append("\t").append(alt).append("\t").append(derived);
            return sb.toString();
        }
        
        @Override
        public boolean equals(Object obj) {
            if (obj instanceof DeleRecord) {
                DeleRecord o = (DeleRecord)obj;
                if (this.chr == o.chr && this.pos == o.pos) return true;
            }
            return false;
        }
        
        @Override
        public int hashCode() {
            return this.pos;
        }

        @Override
        public int compareTo(DeleRecord o) {
            if (chr == o.chr) {
                return this.pos - o.pos;
            }
            return this.chr - o.chr;
        }
    }
    
}


