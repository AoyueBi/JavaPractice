/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.cml247;

import format.Fasta;
import format.PairwiseAlignment;
import format.Table;
import graphcis.r.DensityPlot;
import graphcis.r.Histogram;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeSet;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class CompareGene {
    
    public CompareGene () {
        //this.firstTranscriptB73();
        //this.transcriptCML247();
        //this.countAlignedGene();
        //this.alignedProportion();
        //this.drawAlignedProportion();
    }
    
    public void drawAlignedProportion () {
        String alignProportionFileS = "E:\\Research\\cml247\\gene\\alignment\\cml247VsB73.proportion.txt";
        String pdfFileS = "E:\\Research\\cml247\\gene\\alignment\\cml247VsB73.proportion.pdf";
        Table t = new Table (alignProportionFileS);
        double[] value = new double[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            value[i] = Double.valueOf(t.content[i][4]);
        }
        Histogram h = new Histogram(value);
        h.setBreakNumber(25);
        h.setTitle("Shared genes between CML247 and B73");
        h.setXLab("Proportion of CML247 gene sequence covered by B73 genome");
        h.setYLab("Number of CML247 genes");
        h.setYLim(0, 35000);
        h.setLargerSize();
        h.setLargerSize();
        h.saveGraph(pdfFileS);
    }
    
    public void alignedProportion () {
        String lastzDirS = "E:\\Research\\cml247\\gene\\alignment\\lastz\\";
        String geneFileS = "E:\\Research\\cml247\\gene\\cml247_cds.fasta";
        String alignmentFileS = "E:\\Research\\cml247\\gene\\alignment\\cml247VsB73.lastz.pw.txt";
        String alignProportionFileS = "E:\\Research\\cml247\\gene\\alignment\\cml247VsB73.proportion.txt";
        Fasta f = new Fasta(geneFileS);
        String[] names = new String[f.getSeqNumber()];
        HashMap<String, Integer> lMap = new HashMap();
        for (int i = 0; i < names.length; i++) {
            names[i] = f.getName(i);
            lMap.putIfAbsent(names[i], f.getSeqLength(i));
        }
        Arrays.sort(names);
        File[] fs = new File(lastzDirS).listFiles();
        PairwiseAlignment[] ps = new PairwiseAlignment[fs.length];
        for (int i = 0; i < fs.length; i++) {
            ps[i] = new PairwiseAlignment();
            ps[i].readFromLastzAXF(fs[i].getAbsolutePath());
        }
        PairwiseAlignment prime = ps[0];
        for (int i = 1; i < fs.length; i++) {
            prime.mergeAlignment(ps[i]);
        }
        prime.writeTxtPairwiseAlignment(alignmentFileS);
        String[] alignedNames = prime.getQueries();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(alignProportionFileS);
            bw.write("Gene\tLength\tChr\tPos\tCoveredProportion");
            bw.newLine();
            for (int i = 0; i < names.length; i++) {
                int hit = Arrays.binarySearch(alignedNames, names[i]);
                int chr = -1;
                int pos = -1;
                double portion = 0;
                if (hit < 0) {
                    
                }
                else {
                    int firstIndex = prime.getFirstIndexOfQuery(names[i]);
                    int lastIndex = prime.getLastIndexOfQuery(names[i]);
                    int mainChr = Integer.valueOf(prime.getHit(firstIndex));
                    int mainPos = (prime.getHitStart(firstIndex)+prime.getHitEnd(firstIndex))/2;
                    int upBound = mainPos-500000;
                    int downBound = mainPos+500000;
                    boolean[] cover = new boolean[lMap.get(names[i])];
                    for (int j = firstIndex; j < lastIndex; j++) {
                        if (Integer.valueOf(prime.getHit(j))!=mainChr) continue;
                        if (prime.getHitStart(j) < upBound) continue;
                        if (prime.getHitEnd(j) > downBound) continue;
                        for (int k = prime.getQueryStart(j); k < prime.getQueryEnd(j); k++) {
                            cover[k-1] = true;
                        }
                    }
                    int cnt = 0;
                    for (int j = 0; j < cover.length; j++) {
                        if (cover[j]) cnt++;
                    }
                    portion  = (double)cnt/cover.length;
                    chr = mainChr;
                    pos = mainPos;
                }
                bw.write(names[i]+"\t"+String.valueOf(lMap.get(names[i]))+"\t"+String.valueOf(chr)+"\t"+String.valueOf(pos)+"\t"+String.valueOf(portion));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
       
    }
    
    public void countAlignedGene () {
        //String blastFileS = "E:\\Research\\cml247\\gene\\CML247VsB73.txt";
        String blastFileS = "E:\\Research\\cml247\\gene\\CML247VsB73Genome.txt";
        //String blastFileS = "E:\\Research\\cml247\\gene\\B73VsCML247Genome.txt";
        //String blastFileS = "E:\\Research\\cml247\\gene\\B73VsCML247.txt";
        double cut = 1e-20;
        int presentCnt = 0;
        int totalCnt = 0;
        String temp = null;
        try {
            BufferedReader br = IoUtils.getTextReader(blastFileS);
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("# BLAST ")) break;
                if (temp.startsWith("# BLASTN")) {
                    totalCnt++;
                    for (int i = 0; i < 3; i++) temp = br.readLine();
                }
                if (temp.startsWith("# Fields")) temp = br.readLine();
                String[] tem = temp.split(" ");
                int n = Integer.valueOf(tem[1]);
                boolean flag = true;
                for (int i = 0; i < n; i++) {
                    temp = br.readLine();
                    tem = temp.split("\\s+");
                    if (Double.valueOf(tem[10]) < cut) {
                        if (flag) {
                            presentCnt++;
                            flag = false;
                        }
                    }
                }
            }
            
        }
        catch(Exception e) {
            System.out.println(temp);
            e.printStackTrace();
        }
        double presentRatio = (double)presentCnt/totalCnt;
        System.out.println("Absent:\t"+(totalCnt-presentCnt));
        System.out.println("Total queries:\t:\t"+(totalCnt));
        System.out.println("Absent ratio:\t "+(1-presentRatio));
    }
    
    public void transcriptCML247 () {
        String inputFileS = "E:\\Research\\cml247\\gene\\CML247_Gene_Jiao\\3_Maker-p\\CML247_maker.cdna.fasta";
        String highConfidenceList = "E:\\Research\\cml247\\gene\\CML247_Gene_Jiao\\4_Filter\\summary\\high_confidence.list";
        String outputFileS = "E:\\Research\\cml247\\gene\\cml247_cds.fasta";
        Fasta f = new Fasta(inputFileS);
        ArrayList<String> highList = new ArrayList();
        try {
            BufferedReader br = IoUtils.getTextReader(highConfidenceList);
            String temp;
            while ((temp = br.readLine()) != null) {
                highList.add(temp);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        String[] high = highList.toArray(new String[highList.size()]);
        Arrays.sort(high);
        boolean[] ifOut = new boolean[f.getSeqNumber()];
        for (int i = 0; i < f.getSeqNumber(); i++) {
            if (f.getSeq(i).length()<1) {
                continue;
            }
            String[] temp = f.getName(i).split("-mRNA");
            int hit = Arrays.binarySearch(high, temp[0]);
            if (hit < 0) continue;
            ifOut[i] = true;
        }
        f.writeFasta(outputFileS, ifOut);
    }
    
    public void firstTranscriptB73 () {
        String inputFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.25.cds.all.fa";
        String outputFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.firstCDS.fa";
        Fasta f = new Fasta(inputFileS);
        boolean[] ifOut = new boolean[f.getSeqNumber()];
        for (int i = 0; i < f.getSeqNumber(); i++) {
            String[] temp = f.getName(i).split(" +");
            int n = Integer.valueOf(temp[4].substring(temp[4].length()-2, temp[4].length()));
            if (temp[0].startsWith("GRM")){
                if (n>1) ifOut[i] = false;
                else ifOut[i] = true;
            }
            else {
                ifOut[i] = true;
            }
        }
        f.writeFasta(outputFileS, ifOut);
    }
    
    
}
