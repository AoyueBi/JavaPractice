/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package analysis.deprecated.panAnchor;


import com.itextpdf.awt.DefaultFontMapper;
import com.itextpdf.awt.PdfGraphics2D;
import com.itextpdf.text.Document;
import com.itextpdf.text.DocumentException;
import com.itextpdf.text.Rectangle;
import com.itextpdf.text.pdf.PdfContentByte;
import com.itextpdf.text.pdf.PdfWriter;
import format.Fasta;
import format.Sequence;
import format.ShortreadAlignment;
import format.Table;
import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import net.maizegenetics.dna.map.TagsOnGeneticMap;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import utils.IOFileFormat;

/**
 *
 * @author Fei Lu <fl262@cornell.edu>
 */
public class Mo17Utils {
    
    public Mo17Utils() {}
    
    public void sortGenomeByLength (String inputFileS, String outputFileS) {
        Fasta f = new Fasta (inputFileS);
        f.sortRecordByLengthDescending();
        f.writeFasta(outputFileS);
    }
    
    public void mkDisagreementBlastSyntenyGraph (String samFileS, String compareFileS, String chrInfoFileS,String syntenyDirS) {
        Table t = new Table (chrInfoFileS);
        int chrNum = t.getRowNumber();
        int[] chrID = new int[t.getRowNumber()];
        int[] chrLength = new int[t.getRowNumber()];
        int[] centPos = new int[t.getRowNumber()];
        for (int i = 0; i < chrNum; i++) {
            chrID[i] = Integer.parseInt(t.content[i][0]);
            chrLength[i] = Integer.parseInt(t.content[i][1]);
            centPos[i] = Integer.parseInt(t.content[i][2]);
        }
        t = new Table (compareFileS);
        int cnt = 0;
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.content[i][5].equals("0")) cnt++;
        }
        int contigNum = cnt;
        int[] contigID = new int[contigNum];
        int[] contigLength = new int[contigNum];
        int[] contigAnchorNum = new int[contigNum];
        int[][] anchorContigPos = new int[contigNum][];
        int[][] anchorChr = new int[contigNum][];
        int[][] anchorChrPos = new int[contigNum][];
        cnt = 0;
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.content[i][5].equals("1")) continue;
            contigID[cnt] = Integer.valueOf(t.content[i][0]);
            contigLength[cnt] = Integer.valueOf(t.content[i][1]);
            cnt++;
        }
        ShortreadAlignment sa = new ShortreadAlignment();
        sa.readFromBowtie2(samFileS);
        boolean[] keepList = new boolean[sa.getAlignmentNumber()];
        for (int i = 0; i < keepList.length; i++) {
            String hit = sa.getHit(i);
            if (hit == null) {
                keepList[i] = false;
                continue;
            }
            Pattern p = Pattern.compile("\\D");
            Matcher m = p.matcher(sa.getHit(i));
            if (m.find()) {
                keepList[i] = false;
                continue;
            }
            keepList[i] = true;
        }
        sa.toSubset(keepList);
        sa.sortByQuery();
        for (int i = 0; i < contigID.length; i++) {
            String contigName = "scaffold"+String.valueOf(contigID[i])+"_";
            int startIndex = sa.getAlignmentStartIndexByQueryStartWith(contigName);
            int endIndex = sa.getAlignmentEndIndexByQueryStartWith(contigName);
            contigAnchorNum[i] = endIndex-startIndex;
            anchorContigPos[i] = new int[contigAnchorNum[i]];
            anchorChr[i] = new int[contigAnchorNum[i]];
            anchorChrPos[i] = new int[contigAnchorNum[i]];
            for (int j = 0; j < contigAnchorNum[i]; j++) {
                anchorContigPos[i][j] = Integer.valueOf(sa.getQuery(j+startIndex).split("_")[1])*300;
                anchorChr[i][j] = Integer.valueOf(sa.getHit(j+startIndex));
                anchorChrPos[i][j] = Integer.valueOf(sa.getStartPos(j+startIndex));
            }
        }
        
        int width = 1200;
        int height = 400;
        int maxRecLength = 200;
        int xStart = 50;
        int xInterval = 20;
        int recHeight = 6;
        int yStart = 50;
        int yInterval = 150;
        int contigRecLength = maxRecLength*5+xInterval*4;
        double ratio = (double)maxRecLength/chrLength[0];
        for (int i = 0; i < chrNum; i++) {
            chrLength[i] = (int)(ratio*chrLength[i]);
            centPos[i] = (int)(ratio*centPos[i]);
        }
        
        for (int i = 0; i < contigNum; i++) {
            double ratioContig = (double)contigRecLength/contigLength[i];
            for (int j = 0; j < contigAnchorNum[i]; j++) {
                anchorContigPos[i][j] = (int)(anchorContigPos[i][j]*ratioContig)+xStart;
                int n = (anchorChr[i][j]-1)%5;
                anchorChrPos[i][j] = xStart+n*(maxRecLength+xInterval)+(int)(anchorChrPos[i][j]*ratio);
            }
            try {
                String fileNameS = new File(syntenyDirS, "scaffold"+this.getSixFigureString(contigID[i])).getAbsolutePath()+".pdf";
                Document pd = new Document(new Rectangle(width,height));
                PdfWriter pw;
                pw = PdfWriter.getInstance(pd, new FileOutputStream (fileNameS));
                pd.open();
                PdfContentByte canvas = pw.getDirectContent();
                DefaultFontMapper mapper = new DefaultFontMapper();
                PdfGraphics2D g2d = new PdfGraphics2D(canvas, width, height, mapper);
                int x = xStart;
                int y = yStart;
                g2d.drawString("Scaffold"+String.valueOf(contigID[i])+"\t\tLength="+String.valueOf(contigLength[i]), x, y-25);
                for (int j = 0; j < 5; j++) {
                    g2d.setColor(Color.black);
                    g2d.drawString(String.valueOf(j+1), x, y-5);
                    g2d.setColor(Color.ORANGE);
                    g2d.fillRect(x, y, chrLength[j], recHeight);
                    g2d.setColor(Color.black);
                    g2d.fillRect(x+centPos[j], y, 3, recHeight);
                    x+=chrLength[0];
                    x+=xInterval;
                }
                x = xStart;
                y+=recHeight; int y1 = y;
                y+=yInterval; int y2 = y;
                g2d.setColor(Color.cyan);
                g2d.fillRect(x, y, contigRecLength, recHeight);
                y+=recHeight; int y3 = y;
                y+=yInterval; int y4 = y;
                for (int j = 5; j < chrNum; j++) {
                    g2d.setColor(Color.black);
                    g2d.drawString(String.valueOf(j+1), x, y+recHeight+15);
                    g2d.setColor(Color.ORANGE);
                    g2d.fillRect(x, y, chrLength[j], recHeight);
                    g2d.setColor(Color.black);
                    g2d.fillRect(x+centPos[j], y, 3, recHeight);
                    x+=chrLength[0];
                    x+=xInterval;
                }
                g2d.setColor(Color.blue);
                for (int j = 0; j < contigAnchorNum[i]; j++) {
                    if (anchorChr[i][j] <= 5){
                        g2d.drawLine(anchorChrPos[i][j], y1, anchorContigPos[i][j], y2);
                    }
                    else {
                        g2d.drawLine(anchorContigPos[i][j], y3,anchorChrPos[i][j] , y4);
                    }
                }
                g2d.dispose();
                pd.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }
    
    public void mkCompareAnchorAndAlignmentGoodScaffold (String anchorOnContigFileS, String samFileS, String compareFileS) {
        int contigNum = 0;
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
            contigLength = new int[contigNum];
            contigAnchorNum = new int[contigNum];
            anchorContigPos = new int[contigNum][];
            anchorChr = new int[contigNum][];
            anchorChrPos = new int[contigNum][];
            for (int i = 0; i < 2; i++) br.readLine();
            for (int i = 0; i < contigNum; i++) {
                String[] temp = br.readLine().substring(1).split("\t");
                contigID[i] = Integer.parseInt(temp[0]);
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
        ShortreadAlignment sa = new ShortreadAlignment();
        sa.readFromBowtie2(samFileS);
        String[] querys = sa.getQuerys();
        TreeSet<Integer> idSet = new TreeSet();
        for (int i = 0; i < querys.length; i++) {
            idSet.add(Integer.valueOf(querys[i].split("_")[0].replaceAll("\\D", "")));
        }
        Integer[] ids = idSet.toArray(new Integer[idSet.size()]);
        Arrays.sort(ids);
        ArrayList<Integer>[] chrList = new ArrayList[ids.length];
        for (int i = 0; i < chrList.length; i++) chrList[i] = new ArrayList();
        for (int i = 0; i < sa.getAlignmentNumber(); i++) {
            int id = Integer.valueOf(querys[i].split("_")[0].replaceAll("\\D", ""));
            int index = Arrays.binarySearch(ids, id);
            if (sa.getHit(i) == null) continue;
            Pattern p = Pattern.compile("\\D");
            Matcher m = p.matcher(sa.getHit(i));
            if (m.find()) continue;
            chrList[index].add(Integer.valueOf(sa.getHit(i)));
        }
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(compareFileS), 65536);
            bw.write("ScaffoldID\tLength\tAnchorNum\tAnchorChr\tAlignChr\tAgreement");
            bw.newLine();
            for (int i = 0; i < contigNum; i++) {
                int index = Arrays.binarySearch(ids, contigID[i]);
                if (index < 0) continue;
                int anchorMajority = this.getMajority(anchorChr[i]);
                int[] a = new int[chrList[index].size()];
                for (int j = 0; j < a.length; j++) {
                    a[j] = chrList[index].get(j);
                }
                int alignMajority = this.getMajority(a);
                int agreement = anchorMajority == alignMajority?1:0;
                bw.write(String.valueOf(contigID[i])+"\t"+String.valueOf(contigLength[i])+"\t"+String.valueOf(contigAnchorNum[i])+"\t"+String.valueOf(anchorMajority)+"\t"+String.valueOf(alignMajority)+"\t"+String.valueOf(agreement));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private int getMajority (int[] a) {
        TreeSet<Integer> set = new TreeSet();
        for (int i = 0; i < a.length;  i++) {
            set.add(a[i]);
        }
        Integer[] ar = set.toArray(new Integer[set.size()]);
        Arrays.sort(ar);
        int[] count = new int[ar.length];
        for (int i = 0; i < a.length; i++) {
            int index = Arrays.binarySearch(ar, a[i]);
            count[index]++;
        }
        int max = -1;
        int value = -1;
        for (int i = 0; i < ar.length; i++) {
            if (count[i] > max) {
                max = count[i];
                value = ar[i];
            }
        }
        return value;
    }
    
    public void mkGoodScaffoldFragments (String assemblyQualityFileS, String mo17Ref, String fragmentFileS) {
        Fasta f = new Fasta(mo17Ref);
        Table t = new Table (assemblyQualityFileS);
        ArrayList<String> scaffoldList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (Integer.valueOf(t.content[i][2]) < 20) continue;
            if (Double.valueOf(t.content[i][3]) < 0.95) continue;
            scaffoldList.add("scaffold"+t.content[i][0]);
        }
        String[] scaffoldNames = scaffoldList.toArray(new String[scaffoldList.size()]);
        ArrayList<String> fragNameList = new ArrayList();
        ArrayList<String> fragSeqList = new ArrayList();
        for (int i = 0; i < scaffoldNames.length; i++) {
            int index = f.getIndex(scaffoldNames[i]);
            Sequence s = new Sequence(f.getSeq(index));
            String[] frag = s.getFragments(300);
            for (int j = 0; j < frag.length; j++) {
                fragNameList.add(scaffoldNames[i]+"_"+String.valueOf(j));
                fragSeqList.add(frag[j]);
            }
        }
        String[] fragName = fragNameList.toArray(new String[fragNameList.size()]);
        String[] fragSeq = fragSeqList.toArray(new String[fragSeqList.size()]);
        int[] ids = new int[fragName.length];
        for (int i = 0; i < ids.length; i++) ids[i] = i+1;
        Fasta nf = new Fasta(fragName, fragSeq, ids);
        nf.writeFasta(fragmentFileS);
    }
    
    public void mkGoodContigError (String assemblyQualityFileS, String goodContigList, String goodContigError) {
        Table ta = new Table (assemblyQualityFileS);
        Table tg = new Table (goodContigList);
        String[] goodName = new String[tg.getRowNumber()];
        for (int i = 0; i < goodName.length; i++) {
            goodName[i] = tg.content[i][0];
        }
        Arrays.sort(goodName);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(goodContigError), 65536);
            bw.write("Scaffold\tError");
            bw.newLine();
            for (int i = 0; i < ta.getRowNumber(); i++) {
                int index = Arrays.binarySearch(goodName, ta.content[i][0]);
                if (index < 0) continue;
                bw.write(ta.content[i][0]+"\t"+String.valueOf(1-Double.valueOf(ta.content[i][3])));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkFragmentOnChromosome (String alignFileS, String chrInfoFileS, String fragmentOnChrFileS) {
        Table t = new Table(chrInfoFileS);
        ShortreadAlignment sa = new ShortreadAlignment (alignFileS, IOFileFormat.Text);
        boolean[] keepList = new boolean[sa.getAlignmentNumber()];
        for (int i = 0; i < sa.getAlignmentNumber(); i++) {
            keepList[i] = true;
        }
        sa.toSubset(keepList);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(fragmentOnChrFileS), 65536);
            bw.write("Contig number:\t" + String.valueOf(1));
            bw.newLine();
            bw.write(">Contig\tContigLength\tAnchorNumber");
            bw.newLine();
            bw.write("AnchorIndex\tContigPos\tRefChr\tRefPos\tIfPAV\tAnchorReadLength");
            bw.newLine();
            for (int i = 0; i < 1; i++) {
                int alignmentStartIndex = 0;
                int alignmentEndIndex = sa.getAlignmentNumber();
                int alignmentNumber = alignmentEndIndex - alignmentStartIndex;
                bw.write(">"+String.valueOf(162)+"\t"+String.valueOf(363933)+"\t"+String.valueOf(sa.getAlignmentNumber()));
                bw.newLine();
                for (int j = 0; j < alignmentNumber; j++) {
                    String query = sa.getQuery(j+alignmentStartIndex);
                    bw.write(query+"\t"+String.valueOf(j*300)+"\t");
                    bw.write(sa.getHit(j)+"\t"+String.valueOf(sa.getStartPos(j))+"\t");
                    bw.write(String.valueOf("0")+"\t"+String.valueOf(300));
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void mkScaffoldFragments (String mo17Ref, String scaffold, String fragmentFileS) {
        Fasta f = new Fasta(mo17Ref);
        f.sortRecordByName();
        int index = f.getIndex(scaffold);
        Sequence s = new Sequence(f.getSeq(index));
        Fasta sf = s.getFragmentsFasta(300, "");
        sf.writeFasta(fragmentFileS);
    }
    
    public void mkBarQualityAssembly (String qualityFileS, String barPlotFileS, int intervalLength, int max) {
        double accuracyCut = 0.85;
        int intervalNum = max/intervalLength;
        int[] upperBound = new int[intervalNum];
        int[] lowerBound = new int[intervalNum];
        int[] count = new int[intervalNum];
        ArrayList<Double>[] valueList = new ArrayList[intervalNum];
        for (int i = 0; i < intervalNum; i++) {
            lowerBound[i] = i*intervalLength;
            upperBound[i] = i*intervalLength+intervalLength;
            valueList[i] = new ArrayList();
        }
        Table t = new Table (qualityFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = Arrays.binarySearch(lowerBound, Integer.valueOf(t.content[i][1]));
            if (index < 0) index = -index-2;
            valueList[index].add(Double.valueOf(t.content[i][3]));
        }
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(barPlotFileS), 65536);
            bw.write("Interval\tContigCount\tMean\tSD\tAccuracyRate");
            bw.newLine();
            for (int i = 0; i < intervalNum; i++) {
                bw.write(String.valueOf(lowerBound[i])+"-"+String.valueOf(upperBound[i])+"\t");
                bw.write(String.valueOf(valueList[i].size())+"\t");
                double[] value = new double[valueList[i].size()];
                int cnt = 0;
                for (int j = 0; j < valueList[i].size(); j++) {
                    value[j] = valueList[i].get(j);
                    if (value[j] > accuracyCut) cnt++;
                }
                DescriptiveStatistics ds = new DescriptiveStatistics(value);
                bw.write(String.valueOf(ds.getMean())+"\t"+String.valueOf(ds.getStandardDeviation())+"\t"+String.valueOf((double)cnt/value.length));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkStatisticsAssemblyByRegion (String anchorOnContigFileS, String qualityFileS, String chrInfoFileS, int regionNumPerChr) {
        
    }
    
    public void mkStatisticsAssemblyByChr (String anchorOnContigFileS, String qualityFileS) {
        int contigNum = 0;
        String[] contigName = null;
        int[] contigLength = null;
        int[] contigAnchorNum = null;
        int[][] anchorContigPos = null;
        int[][] anchorChr = null;
        int[][] anchorChrPos = null;
        try {
            BufferedReader br = new BufferedReader (new FileReader(anchorOnContigFileS), 65536);
            contigNum = Integer.parseInt(br.readLine().split("\t")[1]);
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
        double[] mainChrRatioOnContig = new double[contigNum];
        for (int i = 0; i < contigNum; i++) {
            int[] c = new int[10];
            for (int j = 0; j < 10; j++) c[j] = 0;
            for (int j = 0; j < contigAnchorNum[i]; j++) {
                c[anchorChr[i][j]-1]++;
            }
            int mainChrID = 0;
            int mainCnt = 0;
            for (int j = 0; j < 10; j++) {
                if (c[j] > mainCnt) {
                    mainCnt=c[j];
                    mainChrID = j+1;
                }
            }
            int cnt = 0;
            for (int j = 0; j < contigAnchorNum[i]; j++) {
                if (anchorChr[i][j] == mainChrID) cnt++;
            }
            mainChrRatioOnContig[i] = (double)cnt/contigAnchorNum[i];
        }
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(qualityFileS), 65536);
            bw.write("Scaffold\tLength\tAnchorNumber\tRatioAnchorOnMainChr");
            bw.newLine();
            for (int i = 0; i < contigNum; i++) {
                bw.write(contigName[i]+"\t"+String.valueOf(contigLength[i])+"\t"+String.valueOf(contigAnchorNum[i])+"\t");
                bw.write(String.valueOf(mainChrRatioOnContig[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void mkSyntenyGraph (String anchorOnContigFileS, String chrInfoFileS,String syntenyDirS, int actualContigNum, boolean ifDrawAll) {
        Table t = new Table (chrInfoFileS);
        int chrNum = t.getRowNumber();
        int[] chrID = new int[t.getRowNumber()];
        int[] chrLength = new int[t.getRowNumber()];
        int[] centPos = new int[t.getRowNumber()];
        for (int i = 0; i < chrNum; i++) {
            chrID[i] = Integer.parseInt(t.content[i][0]);
            chrLength[i] = Integer.parseInt(t.content[i][1]);
            centPos[i] = Integer.parseInt(t.content[i][2]);
        }
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
        
        int width = 1200;
        int height = 400;
        int maxRecLength = 200;
        int xStart = 50;
        int xInterval = 20;
        int recHeight = 6;
        int yStart = 50;
        int yInterval = 150;
        int contigRecLength = maxRecLength*5+xInterval*4;
        double ratio = (double)maxRecLength/chrLength[0];
        for (int i = 0; i < chrNum; i++) {
            chrLength[i] = (int)(ratio*chrLength[i]);
            centPos[i] = (int)(ratio*centPos[i]);
        }
        int drawNum = contigNum;
        if (!ifDrawAll) drawNum = actualContigNum;
        //for (int i = 0; i < contigNum; i++) {
        for (int i = 0; i < drawNum; i++) {
            double ratioContig = (double)contigRecLength/contigLength[i];
            for (int j = 0; j < contigAnchorNum[i]; j++) {
                anchorContigPos[i][j] = (int)(anchorContigPos[i][j]*ratioContig)+xStart;
                int n = (anchorChr[i][j]-1)%5;
                anchorChrPos[i][j] = xStart+n*(maxRecLength+xInterval)+(int)(anchorChrPos[i][j]*ratio);
            }
            try {
                String fileNameS = new File(syntenyDirS, "contig"+this.getSixFigureString(contigID[i])+"_"+contigName[i]).getAbsolutePath()+".pdf";
                Document pd = new Document(new Rectangle(width,height));
                PdfWriter pw;
                pw = PdfWriter.getInstance(pd, new FileOutputStream (fileNameS));
                pd.open();
                PdfContentByte canvas = pw.getDirectContent();
                DefaultFontMapper mapper = new DefaultFontMapper();
                PdfGraphics2D g2d = new PdfGraphics2D(canvas, width, height, mapper);
                int x = xStart;
                int y = yStart;
                g2d.drawString(contigName[i]+"\t\tLength="+String.valueOf(contigLength[i]), x, y-25);
                for (int j = 0; j < 5; j++) {
                    g2d.setColor(Color.black);
                    g2d.drawString(String.valueOf(j+1), x, y-5);
                    g2d.setColor(Color.ORANGE);
                    g2d.fillRect(x, y, chrLength[j], recHeight);
                    g2d.setColor(Color.black);
                    g2d.fillRect(x+centPos[j], y, 3, recHeight);
                    x+=chrLength[0];
                    x+=xInterval;
                }
                x = xStart;
                y+=recHeight; int y1 = y;
                y+=yInterval; int y2 = y;
                g2d.setColor(Color.cyan);
                g2d.fillRect(x, y, contigRecLength, recHeight);
                y+=recHeight; int y3 = y;
                y+=yInterval; int y4 = y;
                for (int j = 5; j < chrNum; j++) {
                    g2d.setColor(Color.black);
                    g2d.drawString(String.valueOf(j+1), x, y+recHeight+15);
                    g2d.setColor(Color.ORANGE);
                    g2d.fillRect(x, y, chrLength[j], recHeight);
                    g2d.setColor(Color.black);
                    g2d.fillRect(x+centPos[j], y, 3, recHeight);
                    x+=chrLength[0];
                    x+=xInterval;
                }
                g2d.setColor(Color.blue);
                for (int j = 0; j < contigAnchorNum[i]; j++) {
                    if (anchorChr[i][j] <= 5){
                        g2d.drawLine(anchorChrPos[i][j], y1, anchorContigPos[i][j], y2);
                    }
                    else {
                        g2d.drawLine(anchorContigPos[i][j], y3,anchorChrPos[i][j] , y4);
                    }
                }
                g2d.dispose();
                pd.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }
    
    private String getSixFigureString (int number) {
        String s = String.valueOf(number);
        int length = s.length();
        for (int i = 0; i < 6-length; i++) {
            s = "0"+s;
        }
        return s;
    }
    
    /**
     * Only count those anchors unique perfect match in all contigs
     * @param alignFileS
     * @param contigLengthFileS
     * @param anchorFileS
     * @param anchorOnContigFileS
     * @param satisticsFileS 
     */
    public void mkAnchorOnContigFile (String alignFileS, String contigLengthFileS, String anchorFileS, String anchorOnContigFileS, String satisticsFileS) {
        Table t = new Table(contigLengthFileS);
        HashMap<String,Integer> contigLengthMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            contigLengthMap.put(t.content[i][0], Integer.valueOf(t.content[i][1]));
        }
        int[] gChr = null;
        int[] gPos = null;
        byte[] anchorLength = null;
        byte[] ifPav = null;
        try {
            BufferedReader br = new BufferedReader (new FileReader(anchorFileS), 65536);
            br.readLine();
            int cnt = Integer.valueOf(br.readLine());
            gChr = new int[cnt]; gPos = new int[cnt]; anchorLength = new byte[cnt]; ifPav = new byte[cnt];
            br.readLine();
            String[] tem;
            for (int i = 0; i < cnt; i++) {
                tem = br.readLine().split("\t");
                gChr[i] = Integer.valueOf(tem[2]);
                gPos[i] = Integer.valueOf(tem[3]);
                anchorLength[i] = Byte.valueOf(tem[1]);
                ifPav[i] = Byte.valueOf(tem[4]);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        ShortreadAlignment sa = new ShortreadAlignment (alignFileS, IOFileFormat.Text);
        boolean[] keepList = new boolean[sa.getAlignmentNumber()];
        for (int i = 0; i < sa.getAlignmentNumber(); i++) {
            String query = sa.getQuery(i);
            if (sa.isOnlyPerfectMatch(query)) {
                if(sa.isPerfectMatch(i)) keepList[i] = true;
            }
        }
        sa.toSubset(keepList);
        sa.sortByHitAndPos();
        String[] chromosomes = sa.getHits();
        ChromLength[] cl = new ChromLength[chromosomes.length];
        for (int i = 0; i < chromosomes.length; i++) {
            cl[i] = new ChromLength(chromosomes[i], contigLengthMap.get(chromosomes[i]));
        }
        Arrays.sort(cl);
        for (int i = 0; i < chromosomes.length; i++) {
            chromosomes[i] = cl[i].chr;
        }
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(anchorOnContigFileS), 65536);
            bw.write("Contig number:\t" + String.valueOf(chromosomes.length));
            bw.newLine();
            bw.write(">Contig\tContigLength\tAnchorNumber");
            bw.newLine();
            bw.write("AnchorIndex\tContigPos\tRefChr\tRefPos\tIfPAV\tAnchorReadLength");
            bw.newLine();
            for (int i = 0; i < chromosomes.length; i++) {
                int alignmentStartIndex = sa.getAlignmentStartIndexByHit(chromosomes[i]);
                int alignmentEndIndex = sa.getAlignmentEndIndexByHit(chromosomes[i]);
                int alignmentNumber = alignmentEndIndex - alignmentStartIndex;
                bw.write(">"+chromosomes[i]+"\t"+String.valueOf(contigLengthMap.get(chromosomes[i]))+"\t"+String.valueOf(alignmentNumber));
                bw.newLine();
                for (int j = 0; j < alignmentNumber; j++) {
                    String query = sa.getQuery(j+alignmentStartIndex);
                    bw.write(query+"\t"+String.valueOf(sa.getStartPos(j+alignmentStartIndex))+"\t");
                    int index = Integer.valueOf(query);
                    bw.write(String.valueOf(gChr[index])+"\t"+String.valueOf(gPos[index])+"\t");
                    bw.write(String.valueOf(ifPav[index])+"\t"+String.valueOf(anchorLength[i]));
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        int contigNumOnAnchor = chromosomes.length;
        int contigNumTotal = t.getRowNumber();
        long contigLengthOnAnchor = 0;
        long contigLengthTotal = 0;
        for (int i = 0; i < chromosomes.length; i++) contigLengthOnAnchor+=contigLengthMap.get(chromosomes[i]);
        for (int i = 0; i < t.getRowNumber(); i++) {
            contigLengthTotal += Long.valueOf(t.content[i][1]);
        }
        int anchorNumTotal = gChr.length;
        int anchorNumOnContig = sa.getQuerys().length;
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(satisticsFileS), 65536);
            bw.write("Total contig number:\t" + String.valueOf(contigNumTotal));
            bw.newLine();
            bw.write("Number of contig with anchor:\t"+String.valueOf(contigNumOnAnchor)+"\t"+String.valueOf((double)contigNumOnAnchor/contigNumTotal)+"%");
            bw.newLine();
            bw.write("Total contig length:\t" + String.valueOf(contigLengthTotal));
            bw.newLine();
            bw.write("Length of contig with anchor:\t"+String.valueOf(contigLengthOnAnchor)+"\t"+String.valueOf((double)contigLengthOnAnchor/contigLengthTotal)+"%");
            bw.newLine();
            bw.write("Total anchor number:\t" + String.valueOf(anchorNumTotal));
            bw.newLine();
            bw.write("Number of anchor on contig:\t"+String.valueOf(anchorNumOnContig)+"\t"+String.valueOf((double)anchorNumOnContig/anchorNumTotal)+"%");
            bw.newLine();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    class ChromLength implements Comparable <ChromLength>{
        String chr;
        int length;
        
        ChromLength (String chr, int length) {
            this.chr = chr;
            this.length = length;
        }
        
        @Override
        public int compareTo(ChromLength o) {
            return o.length-this.length;
        }
        
    }
}
