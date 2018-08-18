/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.geneticMapping;

import format.ShortreadAlignment;
import format.Table;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import net.maizegenetics.dna.map.TagsOnGeneticMap;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author Fei Lu
 */
public class JumpLibAccuracy {
    
    public JumpLibAccuracy () {
        //this.mkFastaOfAnchor();
        //this.mkFastaOfLib();
        //this.cleanSorghum();
        //this.convertSam();
        //this.mkSelectInfo();
        //this.selectReads();
        //this.compareGeneticAndPhysicalPEMode();
        //this.compareGeneticAndPhysicalSEMode();
        //this.getStatistics();
        this.others();
    }
    
   
    public void others () {
        String forwardFileS = "M:\\2220_7332_6928_H9CD2ADXX_CML247_ACAGTG_R1.fastq.gz";
        String reverseFileS = "M:\\2220_7332_6928_H9CD2ADXX_CML247_ACAGTG_R2.fastq.gz";
        String forwardFasta = "M:\\cml247_f.fas";
        String reverseFasta = "M:\\cml247_r.fas";
        int start = 10000;
        int numSeq = 100000;
        this.writeLibFasta(forwardFileS, forwardFasta, start, numSeq, "f_");
        this.writeLibFasta(reverseFileS, reverseFasta, start, numSeq, "r_");
    }
    
    public void getStatistics () {
        String gAndP = "E:\\Research\\geneticMapping\\jumpLibAccurary\\geneticPhysical_PEMode.txt";
        //String gAndP = "E:\\Research\\geneticMapping\\jumpLibAccurary\\geneticPhysical_SEMode.txt";
        Table t = new Table (gAndP);
        int interval = 10000000;
        int nonMap = 0;
        int wrongMap = 0;
        int rightMap = 0;
        int oneMap = 0;
        int cnt = 0;
        for (int i = 0; i < t.getRowNumber(); i++) {
            
            if (t.content[i][2].equals("*") || t.content[i][4].equals("*")) {
                cnt++;
                nonMap++;
                continue;
            }
            boolean f = this.isSame(Integer.valueOf(t.content[i][0]), Integer.valueOf(t.content[i][1]), Integer.valueOf(t.content[i][2]), Integer.valueOf(t.content[i][3]), interval);
            boolean r = this.isSame(Integer.valueOf(t.content[i][0]), Integer.valueOf(t.content[i][1]), Integer.valueOf(t.content[i][4]), Integer.valueOf(t.content[i][5]), interval);
            if (f && r) {
                rightMap++;
            }
            else if (!f && !r) {
                wrongMap++;
            }
            else {
                oneMap++;
            }
        }
        int bothMapNum = t.getRowNumber()- cnt;
        System.out.println("Total:\t" + String.valueOf(t.getRowNumber()));
        System.out.println("nonMap:\t"+ String.valueOf((double)nonMap/t.getRowNumber()) + "\t"+ nonMap);
        System.out.println("wrongMap:\t"+ String.valueOf((double)wrongMap/t.getRowNumber()) + "\t"+ wrongMap);
        System.out.println("rightMap:\t"+ String.valueOf((double)rightMap/t.getRowNumber()) + "\t"+ rightMap);
        System.out.println("oneMap:\t"+ String.valueOf((double)oneMap/t.getRowNumber()) + "\t"+ oneMap);
        System.out.println();
        System.out.println("Total:\t" + String.valueOf(bothMapNum));
        System.out.println("nonMap:\t"+ String.valueOf((double)nonMap/bothMapNum) + "\t"+ nonMap);
        System.out.println("wrongMap:\t"+ String.valueOf((double)wrongMap/bothMapNum) + "\t"+ wrongMap);
        System.out.println("rightMap:\t"+ String.valueOf((double)rightMap/bothMapNum) + "\t"+ rightMap);
        System.out.println("oneMap:\t"+ String.valueOf((double)oneMap/bothMapNum) + "\t"+ oneMap);
    }
    
    private boolean isSame (int gChr, int gPos, int pChr, int pPos, int interval) {
        if (gChr == pChr) {
            if (Math.abs(gPos-pPos) < interval) return true;
        }
        return false;
    }
    
    public void compareGeneticAndPhysicalSEMode () {
        String fSam = "E:\\Research\\geneticMapping\\jumpLibAccurary\\alignment\\outfile_fS.sam";
        String rSam = "E:\\Research\\geneticMapping\\jumpLibAccurary\\alignment\\outfile_rS.sam";
        String togmFileS = "M:\\production\\panAnchorV2\\v1\\v1.panAnchor.rigid.togm.txt";
        String infoFileS = "E:\\Research\\geneticMapping\\jumpLibAccurary\\query\\info.txt";
        String gAndP = "E:\\Research\\geneticMapping\\jumpLibAccurary\\geneticPhysical_SEMode.txt";
        Table t = new Table (infoFileS);
        HashMap<Integer, Integer> readAnchorMap = new HashMap ();
        for (int i = 0; i < t.getRowNumber(); i++) {
            readAnchorMap.put(Integer.valueOf(t.content[i][1]), Integer.valueOf(t.content[i][0]));
        }
        TagsOnGeneticMap togm = new TagsOnGeneticMap (togmFileS, FilePacking.Text);
        try {
            BufferedReader brF = IoUtils.getTextReader(fSam);
            BufferedReader brR = IoUtils.getTextReader(rSam);
            BufferedWriter bw = IoUtils.getTextWriter(gAndP);
            bw.write("GChr\tGPos\tFChr\tFPos\tRChr\tRPos");
            bw.newLine();
            String tempF = null;
            String tempR = null;
            while (!(tempF = brF.readLine()).startsWith("@PG")){}
            while (!(tempR = brR.readLine()).startsWith("@PG")){}
            while ((tempF = brF.readLine()) != null) {
                tempR = brR.readLine();
                int key = Integer.valueOf(tempF.split("\t")[0].split("_")[1]);
                int index = readAnchorMap.get(key);
                bw.write(String.valueOf(togm.getGChr(index))+"\t"+String.valueOf(togm.getGPos(index))+"\t");
                String[] tem = tempF.split("\t");
                bw.write(tem[2]+"\t"+tem[3]+"\t");
                tem = tempR.split("\t");
                bw.write(tem[2]+"\t"+tem[3]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
            brF.close();
            brR.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void compareGeneticAndPhysicalPEMode () {
        String peSam = "E:\\Research\\geneticMapping\\jumpLibAccurary\\alignment\\outfile_pe.sam";
        String togmFileS = "M:\\production\\panAnchorV2\\v1\\v1.panAnchor.rigid.togm.txt";
        String infoFileS = "E:\\Research\\geneticMapping\\jumpLibAccurary\\query\\info.txt";
        String gAndP = "E:\\Research\\geneticMapping\\jumpLibAccurary\\geneticPhysical_PEMode.txt";
        Table t = new Table (infoFileS);
        HashMap<Integer, Integer> readAnchorMap = new HashMap ();
        for (int i = 0; i < t.getRowNumber(); i++) {
            readAnchorMap.put(Integer.valueOf(t.content[i][1]), Integer.valueOf(t.content[i][0]));
        }
        TagsOnGeneticMap togm = new TagsOnGeneticMap (togmFileS, FilePacking.Text);
        try {
            BufferedReader br = IoUtils.getTextReader(peSam);
            BufferedWriter bw = IoUtils.getTextWriter(gAndP);
            bw.write("GChr\tGPos\tFChr\tFPos\tRChr\tRPos");
            bw.newLine();
            String temp1 = null;
            String temp2 = null;
            while (!(temp1 = br.readLine()).startsWith("@PG")){}
            while ((temp1 = br.readLine()) != null) {
                temp2 = br.readLine();
                int key = Integer.valueOf(temp1.split("\t")[0].split("_")[1]);
                int index = readAnchorMap.get(key);
                bw.write(String.valueOf(togm.getGChr(index))+"\t"+String.valueOf(togm.getGPos(index))+"\t");
                String[] tem = temp1.split("\t");
                bw.write(tem[2]+"\t"+tem[3]+"\t");
                tem = temp2.split("\t");
                bw.write(tem[2]+"\t"+tem[3]);
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
    
    public void selectReads () {
        String forwardFasta = "E:\\Research\\geneticMapping\\jumpLibAccurary\\lib\\cf.fas";
        String reverseFasta = "E:\\Research\\geneticMapping\\jumpLibAccurary\\lib\\cr.fas";
        String forwardS = "E:\\Research\\geneticMapping\\jumpLibAccurary\\query\\fS.fas";
        String reverseS = "E:\\Research\\geneticMapping\\jumpLibAccurary\\query\\rS.fas";
        String infoFileS = "E:\\Research\\geneticMapping\\jumpLibAccurary\\query\\info.txt";
        Table t = new Table (infoFileS);
        int[] index = new int[t.getRowNumber()];
        for (int i = 0; i < index.length; i++) {
            index[i] = Integer.valueOf(t.content[i][1]);
        }
        Arrays.sort(index);
        this.selReads(index, forwardFasta, forwardS);
        this.selReads(index, reverseFasta, reverseS);
    }
    
    private void selReads (int[] index, String inputFileS, String outputFileS) {
        try {
            BufferedReader br = IoUtils.getTextReader(inputFileS);
            BufferedWriter bw = IoUtils.getTextWriter(outputFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                int que = Integer.valueOf(temp.split("_")[1]);
                int hit = Arrays.binarySearch(index, que);
                String seq = br.readLine();
                if (hit < 0) {
                    
                }
                else {
                    bw.write(temp);
                    bw.newLine();
                    bw.write(seq);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkSelectInfo () {
        String saFileS = "E:\\Research\\geneticMapping\\jumpLibAccurary\\alignment\\outfile_cF.sa";
        String infoFileS = "E:\\Research\\geneticMapping\\jumpLibAccurary\\query\\info.txt";
        ShortreadAlignment sa = new ShortreadAlignment(saFileS, IOFileFormat.Binary);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(infoFileS);
            bw.write("AnchorIndex\tReadsIndex");
            bw.newLine();
            for (int i = 0; i < sa.getAlignmentNumber(); i++) {
                if (!sa.isPerfectMatch(i)) continue;
                bw.write(sa.getQuery(i)+"\t"+sa.getHit(i).split("_")[1]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void convertSam () {
        String samFileS = "E:\\Research\\geneticMapping\\jumpLibAccurary\\alignment\\outfile_cF.sam";
        String saFileS = "E:\\Research\\geneticMapping\\jumpLibAccurary\\alignment\\outfile_cF.sa";
        ShortreadAlignment sa = new ShortreadAlignment();
        sa.readFromBowtie2(samFileS);
        sa.writeSimpleAlignment(saFileS, IOFileFormat.Binary);
    }
    
    public void cleanSorghum () {
        String fSamFileS = "E:\\Research\\geneticMapping\\jumpLibAccurary\\alignment\\f.maizeSorghum.sam";
        String rSamFileS = "E:\\Research\\geneticMapping\\jumpLibAccurary\\alignment\\r.maizeSorghum.sam";
        String forwardFasta = "E:\\Research\\geneticMapping\\jumpLibAccurary\\lib\\f.fas";
        String reverseFasta = "E:\\Research\\geneticMapping\\jumpLibAccurary\\lib\\r.fas";
        String cleanForwardFasta = "E:\\Research\\geneticMapping\\jumpLibAccurary\\lib\\cf.fas";
        String cleanReverseFasta = "E:\\Research\\geneticMapping\\jumpLibAccurary\\lib\\cr.fas";
        ShortreadAlignment saf = new ShortreadAlignment();
        saf.readFromBowtie2(fSamFileS);
        ShortreadAlignment sar = new ShortreadAlignment();
        sar.readFromBowtie2(rSamFileS);
        HashSet<Integer> sorghumSet = new HashSet();
        for (int i = 0; i < saf.getAlignmentNumber(); i++) {
            if (saf.getHit(i).startsWith("so")) {
                sorghumSet.add(Integer.valueOf(saf.getQuery(i).split("_")[1]));
            }
        }
        for (int i = 0; i < sar.getAlignmentNumber(); i++) {
            if (sar.getHit(i).startsWith("so")) {
                sorghumSet.add(Integer.valueOf(sar.getQuery(i).split("_")[1]));
            }
        }
        Integer[] outIndex = sorghumSet.toArray(new Integer[sorghumSet.size()]);
        Arrays.sort(outIndex);
        this.outputCleanData(forwardFasta, cleanForwardFasta, outIndex, "f_");
        this.outputCleanData(reverseFasta, cleanReverseFasta, outIndex, "r_");
    }
    private void outputCleanData (String fastaFileS, String outFasta, Integer[] outIndex, String pre) {
        try {
            BufferedReader br = IoUtils.getTextReader(fastaFileS);
            BufferedWriter bw = IoUtils.getTextWriter(outFasta);
            int cnt = 0;
            String temp;
            while ((temp = br.readLine()) != null) {
                String seq = br.readLine();
                int hit = Arrays.binarySearch(outIndex, Integer.valueOf(temp.split("_")[1]));
                if (hit < 0) {
                    bw.write(">"+pre+String.valueOf(cnt));
                    bw.newLine();
                    bw.write(seq);
                    bw.newLine();
                    cnt++;
                }
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkFastaOfLib () {
        String forwardFileS = "N:\\panGenome\\CML247\\mp\\hiseq\\H94Y1ADXX\\3311_7332_10258_H94Y1ADXX_CML_3_5kb_ACAGTG_R1.fastq.gz";
        String reverseFileS = "N:\\panGenome\\CML247\\mp\\hiseq\\H94Y1ADXX\\3311_7332_10258_H94Y1ADXX_CML_3_5kb_ACAGTG_R2.fastq.gz";
        String forwardFasta = "E:\\Research\\geneticMapping\\jumpLibAccurary\\lib\\f.fas";
        String reverseFasta = "E:\\Research\\geneticMapping\\jumpLibAccurary\\lib\\r.fas";
        int start = 10000;
        int numSeq = 500000;
        this.writeLibFasta(forwardFileS, forwardFasta, start, numSeq, "f_");
        this.writeLibFasta(reverseFileS, reverseFasta, start, numSeq, "r_");
    }
    
    private void writeLibFasta (String inputFileS, String outputFileS, int start, int numSeq, String pre) {
        try {
            BufferedReader br = IoUtils.getTextGzipReader(inputFileS);
            BufferedWriter bw = IoUtils.getTextWriter(outputFileS);
            for (int i = 0; i < start*4; i++) br.readLine();
            for (int i = 0; i < numSeq; i++) {
                br.readLine();
                bw.write(">"+pre+String.valueOf(i));
                bw.newLine();
                bw.write(br.readLine());
                bw.newLine();
                br.readLine();
                br.readLine();
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkFastaOfAnchor () {
        String togmFileS = "M:\\production\\panAnchorV2\\v1\\v1.panAnchor.rigid.togm.txt";
        String fastaFileS = "E:\\Research\\geneticMapping\\jumpLibAccurary\\anchor\\v1.togm.rigid.fas";
        TagsOnGeneticMap togm = new TagsOnGeneticMap (togmFileS, FilePacking.Text);
        togm.writeFastA(fastaFileS);
    }
    
    
    
}
