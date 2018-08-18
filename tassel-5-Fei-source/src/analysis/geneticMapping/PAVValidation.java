/*
* To change this license header, choose License Headers in Project Properties.
* To change this template file, choose Tools | Templates
* and open the template in the editor.
*/

package analysis.geneticMapping;

import format.Fasta;
import format.Sequence;
import format.ShortreadAlignment;
import format.Table;
import graphcis.GenomeAlignmentPlot;
import graphcis.r.Histogram;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import net.maizegenetics.dna.map.TagsOnGeneticMap;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import utils.FStringUtils;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author Fei Lu <fl262@cornell.edu>
 */
public class PAVValidation {
    
    public PAVValidation () {
        //this.selectBestCML247();
        //this.splitBestCML247();
        //this.getB73OrthologInfo();
        //this.getB73Ortholog();
        //this.mkBlastLibrary();
        //this.mkBlastAlignment();
        //this.drawBlastAlignment();
        //this.getPAVValidationRate();
        //this.mkPaperFigure();
        
        //*****deprecated******
        //this.mkBWALibrary();
        //this.runBWA();
        //this.processSam();
        //this.drawBWAAlignment();
        //**********************
    }
    
    public void mkPaperFigure () {
        String contigPAVPosDirS = "E:\\Research\\geneticMapping\\cml247\\contigPAVPos\\";
        String selectContigFileS = "E:\\Research\\geneticMapping\\cml247\\paperFigure\\selectedContigInfo.txt";
        String alignmentDirS = "E:\\Research\\geneticMapping\\cml247\\blastAlignmentFilter\\";
        String paperFigureDirS = "E:\\Research\\geneticMapping\\cml247\\paperFigure\\";
        Table t = new Table (selectContigFileS);
        String[] selectContig = new String[t.getRowNumber()];
        int[] contigLength = new int[t.getRowNumber()];
        int[] selectStart = new int[t.getRowNumber()];
        int[] selectEnd = new int[t.getRowNumber()];
        int[] selectLength = new int[t.getRowNumber()];
        int[][] pavPos = new int[t.getRowNumber()][];
        for (int i = 0; i < t.getRowNumber(); i++) {
            selectContig[i] = t.content[i][0];
            selectStart[i] = Integer.valueOf(t.content[i][1]);
            selectEnd[i] = Integer.valueOf(t.content[i][2]);
            selectLength[i] = selectEnd[i] - selectStart[i] + 1;
        }
        for (int i = 0; i < selectContig.length; i++) {
            String infileS = contigPAVPosDirS+selectContig[i]+".pavPos.txt";
            t = new Table (infileS);
            pavPos[i] = new int[t.getRowNumber()];
            for (int j = 0; j < t.getRowNumber(); j++) {
                pavPos[i][j] = Integer.valueOf(t.content[j][0]);
            }
        }
        for (int i = 0; i < selectContig.length; i++) {
            String infileS = alignmentDirS+selectContig[i]+".aln.fil.txt";
            String outfileS = paperFigureDirS + selectContig[i]+".aln.pdf";
            t = new Table (infileS);
            contigLength[i] = Integer.valueOf(t.content[0][1]);
            ArrayList<Integer> qStartList = new ArrayList();
            ArrayList<Integer> qEndList = new ArrayList();
            ArrayList<Integer> hStartList = new ArrayList();
            ArrayList<Integer> hEndList = new ArrayList();
            for (int j = 0; j < t.getRowNumber(); j++) {
                int hStart = Integer.valueOf(t.content[j][5]);
                int hEnd = Integer.valueOf(t.content[j][6]);
                if (hStart < selectStart[i]) continue;
                if (hEnd < selectStart[i]) continue;
                if (hStart > selectEnd[i]) continue;
                if (hEnd > selectEnd[i]) continue;
                qStartList.add(Integer.valueOf(t.content[j][3]));
                qEndList.add(Integer.valueOf(t.content[j][4]));
                hStartList.add(hStart - selectStart[i] + 1);
                hEndList.add(hEnd - selectStart[i]);
            }
            GenomeAlignmentPlot gap = new GenomeAlignmentPlot ("CML247", "B73", contigLength[i], selectLength[i], qStartList, qEndList, hStartList, hEndList);
            int[] pavPosEnd = new int[pavPos[i].length];
            for (int j = 0; j < pavPosEnd.length; j++) {
                pavPosEnd[j] = pavPos[i][j]+1;
            }
            gap.addQueryBrick(pavPos[i], pavPosEnd);
            gap.saveGraph(outfileS);
        }
    }
    
    public void getPAVValidationRate () {
        String filterAlignmentDirS = "E:\\Research\\geneticMapping\\cml247\\blastAlignmentFilter\\";
        String anchorOnContigFileS = "M:\\production\\panGenome\\cml247\\nrgene\\anchorOnContig\\anchorOnContig.txt";
        String validationRateFileS = "E:\\Research\\geneticMapping\\cml247\\pav_validation.txt";
        String anchorFileS = "E:\\Research\\geneticMapping\\identifyPAV\\anchor\\v1.togm.Bowtie2Validation.txt";
        String orthologInfo = "E:\\Research\\geneticMapping\\cml247\\b73_ortho_info.txt";
        Table t = new Table (orthologInfo);
        HashMap<String, Integer> contigLengthMap = new HashMap();
        HashMap<String, Integer> b73LengthMap = new HashMap();
        HashMap<String, Integer> contigChrMap = new HashMap();
        HashMap<String, Integer> contigPosMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            contigLengthMap.put(t.content[i][0], Integer.valueOf(t.content[i][1]));
            b73LengthMap.put(t.content[i][0], Integer.valueOf(t.content[i][5]));
            contigChrMap.put(t.content[i][0], Integer.valueOf(t.content[i][2]));
            contigPosMap.put(t.content[i][0], (Integer.valueOf(t.content[i][3])+Integer.valueOf(t.content[i][4]))/2);
        }
        
        TagsOnGeneticMap togm = new TagsOnGeneticMap(anchorFileS, FilePacking.Text);
        File[] fs = new File(filterAlignmentDirS).listFiles();
        int contigNum = 0;
        String[] contigName = null;
        int[] contigID = null;
        int[] contigLength = null;
        int[] contigAnchorNum = null;
        int[][] anchorIndex = null;
        int[][] anchorContigPos = null;
        boolean[][] ifPAV = null;
        int[][] anchorChr = null;
        int[][] anchorChrPos = null;
        try {
            BufferedReader br = new BufferedReader (new FileReader(anchorOnContigFileS), 65536);
            contigNum = Integer.parseInt(br.readLine().split("\t")[1]);
            contigID = new int[contigNum];
            contigName = new String[contigNum];
            contigLength = new int[contigNum];
            contigAnchorNum = new int[contigNum];
            anchorIndex = new int[contigNum][];
            anchorContigPos = new int[contigNum][];
            ifPAV = new boolean[contigNum][];
            anchorChr = new int[contigNum][];
            anchorChrPos = new int[contigNum][];
            for (int i = 0; i < 2; i++) br.readLine();
            for (int i = 0; i < contigNum; i++) {
                String[] temp = br.readLine().substring(1).split("\t");
                contigName[i] = temp[0];
                contigID[i] = i;
                contigLength[i] = Integer.parseInt(temp[1]);
                contigAnchorNum[i] = Integer.parseInt(temp[2]);
                anchorIndex[i] = new int[contigAnchorNum[i]];
                anchorContigPos[i] = new int[contigAnchorNum[i]];
                ifPAV[i] = new boolean[contigAnchorNum[i]];
                anchorChr[i] = new int[contigAnchorNum[i]];
                anchorChrPos[i] = new int[contigAnchorNum[i]];
                for (int j = 0; j < contigAnchorNum[i]; j++) {
                    temp = br.readLine().split("\t");
                    anchorIndex[i][j] = Integer.parseInt(temp[0]);
                    anchorContigPos[i][j] = Integer.parseInt(temp[1]);
                    anchorChr[i][j] = Integer.parseInt(temp[2]);
                    anchorChrPos[i][j] = Integer.parseInt(temp[3]);
                    if (temp[4].equals("1")) ifPAV[i][j] = true;
                    else ifPAV[i][j] = false;
                }
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        int[] validationCnt = new int[fs.length];
        int[] contigPavNum = new int[fs.length];
        
       
        for (int i = 0; i < fs.length; i++) {
            t = new Table (fs[i].getAbsolutePath());
            int index = 0;
            for (int j = 0; j < contigName.length; j++) {
                if (contigName[j].equals(t.content[0][0])) {
                    index = j;
                    break;
                }
            }
            ArrayList<Integer> pavList = new ArrayList ();
            
                    
            for (int j = 0; j < contigAnchorNum[index]; j++) {
                if (!ifPAV[index][j]) continue;
                int tagIndex = anchorIndex[index][j];
                int gChr = togm.getGChr(tagIndex);
                int gPos = togm.getGPos(tagIndex);
                if (gChr != contigChrMap.get(t.content[0][0])) continue;
                if (Math.abs(gPos- contigPosMap.get(t.content[0][0]))> 10000000) continue;
                
                if (ifPAV[index][j]) pavList.add(anchorContigPos[index][j]);
            }
            Integer[] pavPos = pavList.toArray(new Integer[pavList.size()]);
            int[] contigStart = new int[t.getRowNumber()];
            int[] contigEnd = new int[t.getRowNumber()];
            for (int j = 0; j < t.getRowNumber(); j++) {
                contigStart[j] = Integer.valueOf(t.content[j][3]);
                contigEnd[j] = Integer.valueOf(t.content[j][4]);
            }
            int cnt = 0;
            for (int j = 0; j < pavPos.length; j++) {
                for (int k = 0; k < contigStart.length; k++) {
                    if (pavPos[j]>contigStart[k] && pavPos[j] < contigEnd[k]) {
                      cnt--;
                      break;
                    }
                }
                cnt++;
            }
            validationCnt[i] = cnt;
            contigPavNum[i] = pavPos.length;
        }
        int pavSum = 0;
        int validationSum = 0;
        try {
            BufferedWriter bw = IoUtils.getTextWriter(validationRateFileS);
            bw.write("Contig\tTotalPAV\tValidatedPAV\tRate");
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                bw.write(fs[i].getName().split("\\.")[0]+"\t"+String.valueOf(contigPavNum[i])+"\t"+String.valueOf(validationCnt[i])+"\t"+String.valueOf((double)validationCnt[i]/contigPavNum[i]));
                bw.newLine();
                pavSum+=contigPavNum[i];
                validationSum+=validationCnt[i];
            }
            bw.flush();
            bw.close();
        }
        catch(Exception e) {
            e.printStackTrace();
        }
        System.out.println("\nTotal validation rate:\t" + String.valueOf((double)validationSum/pavSum));
    }
    
    
    
    public void drawBlastAlignment () {
        String orthologInfo = "E:\\Research\\geneticMapping\\cml247\\b73_ortho_info.txt";
        String alignmentDirS = "E:\\Research\\geneticMapping\\cml247\\blastAlignment\\";
        String alignmentFigureDirS = "E:\\Research\\geneticMapping\\cml247\\blastAlignmentFigureFilter\\";
        String filterAlignmentDirS = "E:\\Research\\geneticMapping\\cml247\\blastAlignmentFilter\\";
        String anchorOnContigFileS = "M:\\production\\panGenome\\cml247\\nrgene\\anchorOnContig\\anchorOnContig.txt";
        String anchorFileS = "E:\\Research\\geneticMapping\\identifyPAV\\anchor\\v1.togm.Bowtie2Validation.txt";
        String contigPAVPosDirS = "E:\\Research\\geneticMapping\\cml247\\contigPAVPos\\";
        new File (contigPAVPosDirS).mkdir();
        boolean ifBackBone = true;
        double identityCut = 95;
        double pCut = 1e-10;
        double scoreCut = 0;
        String[] reverseList = {"2", "54", "80", "114", "177", "183", "315", "329", "346", "379", "423", "433", "458"};
        Arrays.sort(reverseList);
        int contigNum = 0;
        String[] contigName = null;
        int[] contigID = null;
        int[] contigLength = null;
        int[] contigAnchorNum = null;
        int[][] anchorIndex = null;
        int[][] anchorContigPos = null;
        boolean[][] ifPAV = null;
        int[][] anchorChr = null;
        int[][] anchorChrPos = null;
        TagsOnGeneticMap togm = new TagsOnGeneticMap(anchorFileS, FilePacking.Text);
        try {
            BufferedReader br = new BufferedReader (new FileReader(anchorOnContigFileS), 65536);
            contigNum = Integer.parseInt(br.readLine().split("\t")[1]);
            contigID = new int[contigNum];
            contigName = new String[contigNum];
            contigLength = new int[contigNum];
            contigAnchorNum = new int[contigNum];
            anchorIndex = new int[contigNum][];
            anchorContigPos = new int[contigNum][];
            ifPAV = new boolean[contigNum][];
            anchorChr = new int[contigNum][];
            anchorChrPos = new int[contigNum][];
            for (int i = 0; i < 2; i++) br.readLine();
            for (int i = 0; i < contigNum; i++) {
                String[] temp = br.readLine().substring(1).split("\t");
                contigName[i] = temp[0];
                contigID[i] = i;
                contigLength[i] = Integer.parseInt(temp[1]);
                contigAnchorNum[i] = Integer.parseInt(temp[2]);
                anchorIndex[i] = new int[contigAnchorNum[i]];
                anchorContigPos[i] = new int[contigAnchorNum[i]];
                ifPAV[i] = new boolean[contigAnchorNum[i]];
                anchorChr[i] = new int[contigAnchorNum[i]];
                anchorChrPos[i] = new int[contigAnchorNum[i]];
                for (int j = 0; j < contigAnchorNum[i]; j++) {
                    temp = br.readLine().split("\t");
                    anchorIndex[i][j] = Integer.parseInt(temp[0]);
                    anchorContigPos[i][j] = Integer.parseInt(temp[1]);
                    anchorChr[i][j] = Integer.parseInt(temp[2]);
                    anchorChrPos[i][j] = Integer.parseInt(temp[3]);
                    if (temp[4].equals("1")) ifPAV[i][j] = true;
                    else ifPAV[i][j] = false;
                }
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        
        Table t = new Table (orthologInfo);
        HashMap<String, Integer> contigLengthMap = new HashMap();
        HashMap<String, Integer> b73LengthMap = new HashMap();
        HashMap<String, Integer> contigChrMap = new HashMap();
        HashMap<String, Integer> contigPosMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            contigLengthMap.put(t.content[i][0], Integer.valueOf(t.content[i][1]));
            b73LengthMap.put(t.content[i][0], Integer.valueOf(t.content[i][5]));
            contigChrMap.put(t.content[i][0], Integer.valueOf(t.content[i][2]));
            contigPosMap.put(t.content[i][0], (Integer.valueOf(t.content[i][3])+Integer.valueOf(t.content[i][4]))/2);
        }
        File[] fs = new File(alignmentDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            String outfileS = fs[i].getName().replaceFirst("txt", "pdf");
            outfileS = alignmentFigureDirS+outfileS;
            String pavPosFileS = fs[i].getName().replaceFirst(".aln.txt", ".pavPos.txt");
            pavPosFileS = contigPAVPosDirS+pavPosFileS;
            try {
                ArrayList<Integer> qStartList = new ArrayList();
                ArrayList<Integer> qEndList = new ArrayList();
                ArrayList<Integer> hStartList = new ArrayList();
                ArrayList<Integer> hEndList = new ArrayList();
                BufferedReader br = IoUtils.getTextReader(fs[i].getAbsolutePath());
                String temp;
                String[] tem = null;
                int cnt = 0;
                boolean reverse = false;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) continue;
                    tem = temp.split("\t");
                    if (Integer.valueOf(tem[8]) > Integer.valueOf(tem[9])) cnt++;
                    tem = br.readLine().split("\t");
                    if (Integer.valueOf(tem[8]) > Integer.valueOf(tem[9])) cnt++;
                    tem = br.readLine().split("\t");
                    if (Integer.valueOf(tem[8]) > Integer.valueOf(tem[9])) cnt++;
                    if (cnt == 3) reverse = true;
                    int index = Arrays.binarySearch(reverseList, tem[0]);
                    if (index >=0) reverse = true;
                    break;
                }
//********************filter by backbone*************************** 
                
                if (ifBackBone) {
                    br = IoUtils.getTextReader(fs[i].getAbsolutePath());
                    int backBoneSize = 10;
                    int[] qbStart = new int[backBoneSize];
                    int[] qbEnd = new int[backBoneSize];
                    int[] hbStart = new int[backBoneSize];
                    int[] hbEnd = new int[backBoneSize];
                    cnt = 0;               
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith("#")) continue;
                        tem = temp.split("\t");
                        int qStart = Integer.valueOf(tem[6]);
                        int qEnd = Integer.valueOf(tem[7]);
                        int hStart;
                        int hEnd;
                        if (reverse) {
                            hStart = b73LengthMap.get(tem[0])-Integer.valueOf(tem[8])+1;
                            hEnd = b73LengthMap.get(tem[0])-Integer.valueOf(tem[9])+1;
                        }
                        else {
                            hStart = Integer.valueOf(tem[8]);
                            hEnd = Integer.valueOf(tem[9]);
                        }
                        if (hStart > hEnd) continue;
                        qbStart[cnt] = qStart;
                        qbEnd[cnt] = qEnd;
                        hbStart[cnt] = hStart;
                        hbEnd[cnt] = hEnd;
                        cnt++;
                        qStartList.add(qStart);
                        qEndList.add(qEnd);
                        hStartList.add(hStart);
                        hEndList.add(hEnd);
                        if (cnt == backBoneSize) break;
                    }

                    br = IoUtils.getTextReader(fs[i].getAbsolutePath());
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith("#")) continue;
                        tem = temp.split("\t");
                        if (Double.valueOf(tem[10]) > pCut) continue;
                        if (Double.valueOf(tem[11]) < scoreCut) continue;
                        if (Double.valueOf(tem[2]) < identityCut) continue;
                        int qStart = Integer.valueOf(tem[6]);
                        int qEnd = Integer.valueOf(tem[7]);
                        int hStart;
                        int hEnd;
                        if (reverse) {
                            hStart = b73LengthMap.get(tem[0])-Integer.valueOf(tem[8])+1;
                            hEnd = b73LengthMap.get(tem[0])-Integer.valueOf(tem[9])+1;
                        }
                        else {
                            hStart = Integer.valueOf(tem[8]);
                            hEnd = Integer.valueOf(tem[9]);
                        }
                        boolean flag = false;
                        for (int j = 0; j < qbStart.length; j++) {
                            if (qbStart[j] == qStart) {
                                flag = true;
                                break;
                            }
                        }
                        if (flag) continue;
                        flag = false;
                        for (int j = 0; j < qbStart.length; j++) {
                            if (qStart < qbStart[j] && hStart > hbStart[j]) {
                                flag = true;
                                break;
                            }
                            if (qStart > qbStart[j] && hStart < hbStart[j]) {
                                flag = true;
                                break;
                            }
                        }
                        if (flag) continue;
                        qStartList.add(qStart);
                        qEndList.add(qEnd);
                        hStartList.add(hStart);
                        hEndList.add(hEnd);
                    }
                }
                    
  //**************************************************************************
                
  //***********************simgple filter by e-evalue and score**************** 
                else {
                    br = IoUtils.getTextReader(fs[i].getAbsolutePath());
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith("#")) continue;
                        tem = temp.split("\t");
                        if (Double.valueOf(tem[10]) > pCut) continue;
                        if (Double.valueOf(tem[11]) < scoreCut) continue;
                        if (Double.valueOf(tem[2]) < identityCut) continue;

                        qStartList.add(Integer.valueOf(tem[6]));
                        qEndList.add(Integer.valueOf(tem[7]));
                        if (reverse) {
                            hStartList.add(b73LengthMap.get(tem[0])-Integer.valueOf(tem[8])+1);
                            hEndList.add(b73LengthMap.get(tem[0])-Integer.valueOf(tem[9])+1);
                        }
                        else {
                            hStartList.add(Integer.valueOf(tem[8]));
                            hEndList.add(Integer.valueOf(tem[9]));
                        }
                    }
                }
                
                int index = 0;
                for (int j = 0; j < contigName.length; j++) {
                    if (contigName[j].equals(tem[0])) {
                        index = j;
                        break;
                    }
                }
                ArrayList<Integer> pavPosStartList = new ArrayList();
                ArrayList<Integer> pavPosEndList = new ArrayList();
                for (int j = 0; j < contigAnchorNum[index]; j++) {
                    if (!ifPAV[index][j]) continue;
                    int tagIndex = anchorIndex[index][j];
                    int gChr = togm.getGChr(tagIndex);
                    int gPos = togm.getGPos(tagIndex);
                    if (gChr != contigChrMap.get(tem[0])) continue;
                    if (Math.abs(gPos- contigPosMap.get(tem[0]))> 10000000) continue;
                    pavPosStartList.add(anchorContigPos[index][j]);
                    pavPosEndList.add(anchorContigPos[index][j]+1);
                }
                BufferedWriter bw = IoUtils.getTextWriter(pavPosFileS);
                bw.write("ContigStartPos");
                bw.newLine();
                for (int j = 0; j < pavPosStartList.size(); j++) {
                    bw.write(String.valueOf(pavPosStartList.get(j)));
                    bw.newLine();
                }
                bw.flush();
                bw.close();
  //*******************************************************************************************             
                GenomeAlignmentPlot gap = new GenomeAlignmentPlot ("CML247", "B73", contigLengthMap.get(tem[0]), b73LengthMap.get(tem[0]), qStartList, qEndList, hStartList, hEndList);
                gap.addQueryBrick(pavPosStartList, pavPosEndList);
                gap.saveGraph(outfileS);
                String alignmentOutfileS = fs[i].getName().replaceFirst("txt", "fil.txt");
                alignmentOutfileS = filterAlignmentDirS + alignmentOutfileS;
                bw = IoUtils.getTextWriter(alignmentOutfileS);
                bw.write("Name\tQueryLength\tHitLength\tQStart\tQEnd\tHStart\tHEnd\tIfReverse");
                bw.newLine();
                for (int j = 0; j < qStartList.size(); j++) {
                    bw.write(tem[0]+"\t"+contigLengthMap.get(tem[0])+"\t"+b73LengthMap.get(tem[0])+"\t"+String.valueOf(qStartList.get(j))+"\t"+String.valueOf(qEndList.get(j))+"\t");
                    bw.write(String.valueOf(hStartList.get(j))+"\t"+String.valueOf(hEndList.get(j))+"\t");
                    if (reverse) bw.write("1");
                    else bw.write("0");
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
    
    public void mkBlastAlignment () {
        String b73OrthoDirS = "E:\\Research\\geneticMapping\\cml247\\b73ortho\\";
        String cml247BestDirS = "E:\\Research\\geneticMapping\\cml247\\contigBest\\";
        String resultDirS = "E:\\Research\\geneticMapping\\cml247\\blastAlignment\\";
        String script = "E:\\Research\\geneticMapping\\cml247\\runBlast.pl";
        
        File[] queries = new File(cml247BestDirS).listFiles();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(script);
            for (int i = 0; i < queries.length; i++) {
                String name = queries[i].getName().replaceAll(".fa", "");
                String db = new File(b73OrthoDirS, "b73_"+name+".fa").getAbsolutePath().replaceAll("\\\\", "/");
                String outfileS = new File(resultDirS, name+".aln.txt").getAbsolutePath().replaceAll("\\\\", "/");
                bw.write("system \"blastn -query " + queries[i].getAbsolutePath().replaceAll("\\\\", "/"));
                bw.write(" -db "+db+" -out " + outfileS + " -evalue 1e-20 -outfmt 7 -num_threads 3\";");
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkBlastLibrary () {
        String b73OrthoDirS = "E:\\Research\\geneticMapping\\cml247\\b73ortho\\";
        String script = "E:\\Research\\geneticMapping\\cml247\\mkLib.pl";
        
        File[] fs = new File(b73OrthoDirS).listFiles();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(script);
            for (int i = 0; i < fs.length; i++) {
                bw.write("system \"makeblastdb -in " + fs[i].getAbsolutePath().replaceAll("\\\\", "/") + " -dbtype nucl\";");
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void getB73Ortholog () {
        String orthologInfo = "E:\\Research\\geneticMapping\\cml247\\b73_ortho_info.txt";
        String b73Genome = "N:\\Zea\\reference_genome\\ZmB73_RefGen_v2.fa";
        String b73OrthoDirS = "E:\\Research\\geneticMapping\\cml247\\b73ortho\\";
        Table t = new Table(orthologInfo);
        Fasta f = new Fasta(b73Genome);
        f.sortRecordByName();
        try {
            for (int i = 0; i < t.getRowNumber(); i++) {
                String outfileS = b73OrthoDirS+"b73_"+t.content[i][0]+".fa";
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                String chr = t.content[i][2];
                int startPos = Integer.valueOf(t.content[i][3]);
                int endPos = Integer.valueOf(t.content[i][4]);
                bw.write(">b73_"+t.content[i][0]);
                bw.newLine();
                int index = f.getIndex(chr);
                bw.write(f.getSeq(index).substring(startPos, endPos));
                bw.newLine();
                bw.flush();
                bw.close();
            }
            
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void getB73OrthologInfo () {
        String cml247BestFastaS = "E:\\Research\\geneticMapping\\cml247\\contig_best.fas";
        String anchorOnContigFileS = "M:\\production\\panGenome\\cml247\\nrgene\\anchorOnContig\\anchorOnContig.txt";
        String orthologInfo = "E:\\Research\\geneticMapping\\cml247\\b73_ortho_info.txt";
        String b73Genome = "N:\\Zea\\reference_genome\\ZmB73_RefGen_v2.fa";
        int contigNum = 0;
        String[] contigName = null;
        int[] contigID = null;
        int[] contigLength = null;
        int[] contigAnchorNum = null;
        int[][] anchorContigPos = null;
        int[][] anchorChr = null;
        int[][] anchorChrPos = null;
        try {
            BufferedReader br = new BufferedReader (new FileReader(anchorOnContigFileS), 65536);
            contigNum = Integer.parseInt(br.readLine().split("\t")[1]);
            contigID = new int[contigNum];
            contigName = new String[contigNum];
            contigLength = new int[contigNum];
            contigAnchorNum = new int[contigNum];
            anchorContigPos = new int[contigNum][];
            anchorChr = new int[contigNum][];
            anchorChrPos = new int[contigNum][];
            for (int i = 0; i < 2; i++) br.readLine();
            for (int i = 0; i < contigNum; i++) {
                String[] temp = br.readLine().substring(1).split("\t");
                contigName[i] = temp[0];
                contigID[i] = i;
                contigLength[i] = Integer.parseInt(temp[1]);
                contigAnchorNum[i] = Integer.parseInt(temp[2]);
                anchorContigPos[i] = new int[contigAnchorNum[i]];
                anchorChr[i] = new int[contigAnchorNum[i]];
                anchorChrPos[i] = new int[contigAnchorNum[i]];
                for (int j = 0; j < contigAnchorNum[i]; j++) {
                    temp = br.readLine().split("\t");
                    anchorContigPos[i][j] = Integer.parseInt(temp[1]);
                    anchorChr[i][j] = Integer.parseInt(temp[2]);
                    anchorChrPos[i][j] = Integer.parseInt(temp[3]);
                }
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        Fasta fs = new Fasta(cml247BestFastaS);
        String[] bestContig = new String[fs.getSeqNumber()];
        for (int i = 0; i < fs.getSeqNumber(); i++) {
            bestContig[i] = fs.getName(i);
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(orthologInfo);
            bw.write("ContigName\tContigLength\tB73Chr\tB73StartPos\tB73EndPos\tB73Length\tLengthRatio");
            bw.newLine();
            for (int i = 0; i < bestContig.length; i++) {
                int index = 0;
                for (int j = 0; j < contigName.length; j++) {
                    if (bestContig[i].equals(contigName[j])) {
                        index = j;
                        break;
                    }
                }
                int[] temp = new int[11];
                for (int j = 0; j < contigAnchorNum[index]; j++) {
                    temp[anchorChr[index][j]]++;
                }
                int chr = 1;
                int max = 0;
                for (int j = 1; j < temp.length; j++) {
                    if (temp[j] > max) {
                        max = temp[j];
                        chr = j;
                    }
                }
                ArrayList<Integer> posList = new ArrayList();
                for (int j = 0; j < contigAnchorNum[index]; j++) {
                    if (anchorChr[index][j] != chr) continue;
                    posList.add(anchorChrPos[index][j]);
                }
                Integer[] chrPos = posList.toArray(new Integer[posList.size()]);
                Arrays.sort(chrPos);
                int startIndex = (int)(chrPos.length*0.05);
                int endIndex = (int)(chrPos.length*0.95);
                int startPos = chrPos[startIndex] - (int)(fs.getSeqLength(i)*0.25);
                int endPos = chrPos[endIndex] + (int)(fs.getSeqLength(i)*0.25);
                int b73Length = endPos-startPos;
                bw.write(String.valueOf(bestContig[i])+"\t"+String.valueOf(fs.getSeqLength(i))+"\t"+String.valueOf(chr)+"\t"+String.valueOf(startPos)+"\t"+String.valueOf(endPos));
                bw.write("\t"+b73Length+"\t"+(double)fs.getSeqLength(i)/b73Length);
                bw.newLine();
                System.out.println(">"+bestContig[i]);
                for (int j = 0; j < chrPos.length; j++) {
                    System.out.println(chrPos[j]);
                }
                System.out.println("\n\n\n");
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void splitBestCML247 () {
        String cml247BestFastaS = "E:\\Research\\geneticMapping\\cml247\\contig_best.fas";
        String cml247BestDirS = "E:\\Research\\geneticMapping\\cml247\\contigBest\\";
        Fasta f = new Fasta (cml247BestFastaS);
        for (int i = 0; i < f.getSeqNumber(); i++) {
            String outfileS = cml247BestDirS+f.getName(i)+".fa";
            try{
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                bw.write(">"+f.getName(i));
                bw.newLine();
                bw.write(f.getSeq(i));
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
    }
    
    public void selectBestCML247 () {
        String maizeSorghumRatioFileS = "M:\\production\\panGenome\\cml247\\nrgene\\maizeSorghumRatio\\contig_maizeSorghumRatio.txt";
        String qualityFileS = "M:\\production\\panGenome\\cml247\\nrgene\\anchorOnContig\\contig_quality.txt";
        String cml247FastaS = "M:\\production\\panGenome\\cml247\\nrgene\\fasta\\contig_sorted.fas";
        String cml247BestFastaS = "E:\\Research\\geneticMapping\\cml247\\contig_best.fas";
        new GMUtils().selectBestCML247(maizeSorghumRatioFileS, qualityFileS, cml247FastaS, cml247BestFastaS, 200, 0.95, 100, 0.98);
    }
    
    public void drawBWAAlignment () {
        String orthologInfo = "E:\\Research\\geneticMapping\\cml247\\b73_ortho_info.txt";
        String alignmentDirS = "E:\\Research\\geneticMapping\\cml247\\bwaAlignment\\";
        String alignmentFigureDirS = "E:\\Research\\geneticMapping\\cml247\\bwaAlignmentFigure\\";
        int fragmentSize = 200;
       
        Table t = new Table (orthologInfo);
        HashMap<String, Integer> contigLengthMap = new HashMap();
        HashMap<String, Integer> b73LengthMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            contigLengthMap.put(t.content[i][0], Integer.valueOf(t.content[i][1]));
            b73LengthMap.put(t.content[i][0], Integer.valueOf(t.content[i][5]));
        }
        File[] fs = new File(alignmentDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            String outfileS = fs[i].getName().replaceFirst("sa", "pdf");
            outfileS = alignmentFigureDirS+outfileS;
            try {
                ArrayList<Integer> qStartList = new ArrayList();
                ArrayList<Integer> qEndList = new ArrayList();
                ArrayList<Integer> hStartList = new ArrayList();
                ArrayList<Integer> hEndList = new ArrayList();
                ShortreadAlignment sa = new ShortreadAlignment (fs[i].getAbsolutePath(), IOFileFormat.Text);
                for (int j = 0; j < sa.getAlignmentNumber(); j++) {
                    if (!sa.isMatch(j)) continue;
                    int index = Integer.valueOf(sa.getQuery(j));
                    qStartList.add(1+index*fragmentSize);
                    qEndList.add(index*fragmentSize+fragmentSize);
                    hStartList.add(sa.getStartPos(j));
                    hEndList.add(sa.getEndPos(j));
                }
                
                String name = fs[i].getName().replaceFirst(".sa", "");
                GenomeAlignmentPlot gap = new GenomeAlignmentPlot ("CMl247", "B73", contigLengthMap.get(name), b73LengthMap.get(name), qStartList, qEndList, hStartList, hEndList);
                gap.saveGraph(outfileS);
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
    }
    
    public void processSam () {
        String samDirS = "/SSD/mingh/bwaSam/";
        String alignmentDirS = "/SSD/mingh/bwaAlignment/";
        new File(alignmentDirS).mkdir();
        
        File[] fs = new File(samDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            String outfileS = fs[i].getName().replaceFirst("sam", "sa");
            outfileS = alignmentDirS+outfileS;
            ShortreadAlignment sa = new ShortreadAlignment();
            sa.readFromBWAMEM(fs[i].getAbsolutePath());
            sa.writeSimpleAlignment(outfileS, IOFileFormat.Text);
        }
    }
    
    public void runBWA () {
        String b73OrthoDirS = "/SSD/mingh/b73ortho/";
        String bestContigFileS = "/SSD/mingh/contig_best.fas";
        String samDirS = "/SSD/mingh/bwaSam/";
        String fastqDirS = "/SSD/mingh/fastq/";
        String script = "/SSD/mingh/runBWA.pl"; 
        Fasta f = new Fasta(bestContigFileS);
        int fragmentSize = 200;
        String qual  = "";
        for (int i = 0; i < fragmentSize; i++) qual = qual+"a";
        new File(fastqDirS).mkdir();
        for (int i = 0; i < f.getSeqNumber(); i++) {
            String[] subs = new Sequence(f.getSeq(i)).getFragments(fragmentSize);
            try {
                String outfileS = fastqDirS + f.getName(i)+".fq";
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                for (int j = 0; j < subs.length; j++) {
                    bw.write("@"+FStringUtils.getNDigitNumber(5, j));
                    bw.newLine();
                    bw.write(subs[j]);
                    bw.newLine();
                    bw.write("+");
                    bw.newLine();
                    bw.write(qual.substring(0, subs[j].length()));
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        new File(samDirS).mkdir();
        
        try {
            BufferedWriter bw = IoUtils.getTextWriter(script);
            for (int i = 0; i < f.getSeqNumber(); i++) {
                String samFileS = samDirS + f.getName(i)+".sam";
                String dbFileS = b73OrthoDirS+"b73_"+f.getName(i)+".fa";
                String queryFileS = fastqDirS + f.getName(i)+".fq";
                bw.write("system \"" + "bwa mem " + dbFileS + " " + queryFileS + " > " + samFileS + "\";");
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkBWALibrary () {
        String b73OrthoDirS = "/SSD/mingh/b73ortho/";
        String script = "/SSD/mingh/indexBWA.pl";
        
        File[] fs = new File(b73OrthoDirS).listFiles();
        try {
            BufferedWriter bw =  IoUtils.getTextWriter(script);
            for (int i = 0; i < fs.length; i++) {
                bw.write("system \"bwa  index " + fs[i].getAbsolutePath() + "\";");
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
