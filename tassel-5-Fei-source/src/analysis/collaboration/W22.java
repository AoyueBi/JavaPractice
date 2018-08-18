/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.collaboration;

import analysis.panGenome.denovoEvaluation.ScaffoldWithAnchor;
import com.itextpdf.awt.DefaultFontMapper;
import com.itextpdf.awt.PdfGraphics2D;
import com.itextpdf.text.Document;
import com.itextpdf.text.Rectangle;
import com.itextpdf.text.pdf.PdfContentByte;
import com.itextpdf.text.pdf.PdfWriter;
import format.Bins;
import format.Fasta;
import format.Range;
import format.Sequence;
import format.Table;
import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class W22 {
    String scaffoldFileS = "M:\\production\\panGenome\\w22\\source\\201504\\W22.scaffolds.fa";
    
    public W22 () {
        //this.getGenomeStatistics();
        //this.getBz1Scaffold();
        //this.getBz1Region();
        //this.processMaf();
        //this.processMafForDot();
        //this.drawGenomeAlignment();
        this.drawLineChart();
    }
    
    public void drawLineChart () {
        String gbbFileS = "M:\\collaboration\\w22\\source\\EU338354.gb.txt";
        String mafFileS = "M:\\collaboration\\w22\\W22RegionOnW22Bac.maf.txt";
        String alignedRegionFileS = "M:\\collaboration\\w22\\W22RegionOnW22Bac.maf.processed.lineChart.txt";
        String pdfFileS = "M:\\collaboration\\w22\\geneModel.pdf";
        String outfileS = "M:\\collaboration\\w22\\mismatch.linechart.txt";
        int bacLength;
        boolean[] ifMis = null;
        ArrayList<Range> geneL = new ArrayList();
        ArrayList<Range> repeatL = new ArrayList();
        Table t = new Table (alignedRegionFileS);
        int[] hitStarts = new int[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            hitStarts[i] = t.getIntValue(i, 2);
        }
        Arrays.sort(hitStarts);
        try {
            BufferedReader br = IoUtils.getTextReader(gbbFileS);
            String temp = br.readLine();
            String[] tem = temp.split("\\s++");
            bacLength = Integer.valueOf(tem[2]);
            ifMis = new boolean[bacLength];
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("ORIGIN")) break;
                else if (temp.startsWith("     gene")) {
                    temp = temp.replaceFirst("\\s+gene\\s+", "");
                    temp = temp.replaceFirst(".+\\(", "");
                    temp = temp.replaceFirst("\\)", "");
                    temp = temp.replaceAll(">|<", "");
                    tem = temp.split("\\.\\.");
                    temp = br.readLine();
                    Range r = new Range(1, Integer.valueOf(tem[0]), Integer.valueOf(tem[1]));
                    if (temp.contains("gag") || temp.contains("pol") || temp.contains("TNP")) {
                        repeatL.add(r);
                    }
                    else if (temp.contains("hypro")) {
                        
                    }
                    else {
                        System.out.println(temp+"\t"+String.valueOf((double)r.getRangeStart()/bacLength));
                        geneL.add(r);
                    }
                }
                else if (temp.startsWith("     mobile_element")) {
                    temp = temp.replaceFirst("\\s+mobile_element\\s+", "");
                    temp = temp.replaceFirst(".+\\(", "");
                    temp = temp.replaceFirst("\\)", "");
                    tem = temp.split("\\.\\.");
                    Range r = new Range(1, Integer.valueOf(tem[0]), Integer.valueOf(tem[1]));
                    repeatL.add(r);
                }
                else if (temp.startsWith("     LTR")) {
                    temp = temp.replaceFirst("\\s+LTR\\s+", "");
                    temp = temp.replaceFirst(".+\\(", "");
                    temp = temp.replaceFirst("\\)", "");
                    tem = temp.split("\\.\\.");
                    Range r = new Range(1, Integer.valueOf(tem[0]), Integer.valueOf(tem[1]));
                }
            }
            br.close();
            br = IoUtils.getTextReader(mafFileS);
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("a score")) {
                    tem = br.readLine().split("\\s+");
                    int startPoint = Integer.valueOf(tem[2]);
                    if (Arrays.binarySearch(hitStarts, startPoint) < 0) continue;
                    String hitS = tem[6];
                    temp = br.readLine();
                    tem = temp.split("\\s+");
                    int startIndex = Integer.valueOf(tem[2]);
                    String queryS = tem[6];
                    for (int i = 0; i < queryS.length(); i++) {
                        char q = queryS.charAt(i);
                        char h = hitS.charAt(i);
                        if (q == '-') continue;
                        if (h == '-') continue;
                        if (q != h) {
                            ifMis[startIndex] = true;
                        }
                        startIndex++;
                    }
                }
                
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        //int corStart, int corEnd, int rangeLength, int[] coordinate, double[] v
        int[] coordinate = new int[ifMis.length];
        double[] value = new double[ifMis.length];
        for (int i = 0; i < value.length; i++) {
            coordinate[i] = i+1;
            if (ifMis[i]) value[i] = 1;
        }
        Bins bin = new Bins (1, ifMis.length, 1000, coordinate, value);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Position\tMismatch");
            bw.newLine();
            for (int i = 0; i < bin.getBinNum(); i++) {
                bw.write(String.valueOf(bin.getBinStart(i)-1)+"\t"+String.valueOf(bin.getTotalValue(i)));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
        try {
            int width = 1600;
            int height = 1200;
            int xMargin = 100;
            int yMargin = 500;
            int barLength = width -2* xMargin;
            int barHeight = 30;
            Document pd = new Document(new Rectangle(width,height));
            PdfWriter pw;
            pw = PdfWriter.getInstance(pd, new FileOutputStream (pdfFileS));
            pd.open();
            PdfContentByte canvas = pw.getDirectContent();
            DefaultFontMapper mapper = new DefaultFontMapper();
            PdfGraphics2D g2d = new PdfGraphics2D(canvas, width, height, mapper);
            bacLength = ifMis.length;
            int[][] genePos = new int[geneL.size()][2];
            int[][] repeatPos = new int[repeatL.size()][2];
            double ratio = (double)barLength/bacLength;
            for (int i = 0; i <  genePos.length; i++) {
                genePos[i][0] = (int)(geneL.get(i).getRangeStart()*ratio);
                genePos[i][1] = (int)((geneL.get(i).getRangeEnd())*ratio);
            }
            for (int i = 0; i <  repeatPos.length; i++) {
                repeatPos[i][0] = (int)(repeatL.get(i).getRangeStart()*ratio);
                repeatPos[i][1] = (int)((repeatL.get(i).getRangeEnd())*ratio);
            }
            Color c = new Color(0, 0, 255, 30);
            g2d.setColor(c);
            g2d.fillRect(xMargin, yMargin, barLength, barHeight);
            g2d.setColor(Color.RED);
            for (int i = 0; i <  repeatPos.length; i++) {
                g2d.fillRect(xMargin+repeatPos[i][0], yMargin, repeatPos[i][1]-repeatPos[i][0], barHeight);
            }
            g2d.setColor(Color.BLACK);
            for (int i = 0; i <  genePos.length; i++) {
                g2d.fillRect(xMargin+genePos[i][0], yMargin, genePos[i][1]-genePos[i][0], barHeight);
            }
           
            g2d.dispose();
            pd.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    public void drawGenomeAlignment () {
        String infoFileS = "E:\\Database\\InfoFile\\ChrLenCentPosi_agpV3.txt";
        String mappingFileS = "M:\\collaboration\\w22\\geneticValidation\\scaffold242506.txt";
        String pdfFileS = "M:\\collaboration\\w22\\geneticValidation\\bz1GeneticAnchor.pdf";
        Table t = new Table (infoFileS);
        int chrNum = t.getRowNumber();
        int[] chr = new int[t.getRowNumber()];
        int[] chrLength = new int[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            chr[i] = Integer.valueOf(t.content[i][0]);
            chrLength[i] = Integer.valueOf(t.content[i][1]);
        }
        String scaffoldName = null;
        int scaffoldLength = 0;
        int anchorNum = 0;
        int[] anchorScaffoldPos = null;
        int[] anchorGChr = null;
        int[] anchorGPos = null;
        try {
            BufferedReader br = IoUtils.getTextReader(mappingFileS);
            String temp = br.readLine().replaceFirst(">", "");
            String[] tem = temp.split("\t");
            scaffoldName = tem[0];
            scaffoldLength = Integer.valueOf(tem[1]);
            anchorNum = Integer.valueOf(tem[2]);
            anchorScaffoldPos = new int[anchorNum];
            anchorGChr = new int[anchorNum];
            anchorGPos = new int[anchorNum];
            for (int i = 0; i < anchorNum; i++) {
                tem = br.readLine().split("\t");
                anchorScaffoldPos[i] = Integer.valueOf(tem[1]);
                anchorGChr[i] = Integer.valueOf(tem[2]);
                anchorGPos[i] = Integer.valueOf(tem[3]);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
        int width = 1300;
        int height = 400;
        int maxRecLength = 100;
        int xStart = 50;
        int xInterval = 20;
        int recHeight = 6;
        int yStart = 50;
        int yInterval = 150;
        int contigRecLength = maxRecLength*10+xInterval*9;
        double ratio = (double)maxRecLength/chrLength[0];
        for (int i = 0; i < chrNum; i++) {
            chrLength[i] = (int)(ratio*chrLength[i]);
        }
        int markUnit = 1000000;
        int[] mark = new int[scaffoldLength/markUnit];
        int[] bzPos = {3489235, 3721538};
        double ratioScaffold = (double)contigRecLength/scaffoldLength;
        for (int i = 0; i < mark.length; i++) {
            mark[i] = (int)(markUnit*i*ratioScaffold)+xStart;
        }
        for (int i = 0; i < bzPos.length; i++) {
            bzPos[i] = (int)(bzPos[i]*ratioScaffold)+xStart;
        }
        for (int i = 0; i < anchorNum; i++) {
            anchorScaffoldPos[i] = (int)(anchorScaffoldPos[i]*ratioScaffold)+xStart;
            int n = (anchorGChr[i]-1)%10;
            anchorGPos[i] = xStart+n*(maxRecLength+xInterval)+(int)(anchorGPos[i]*ratio);
        }
        
        try {
            Document pd = new Document(new Rectangle(width,height));
            PdfWriter pw;
            pw = PdfWriter.getInstance(pd, new FileOutputStream (pdfFileS));
            pd.open();
            PdfContentByte canvas = pw.getDirectContent();
            DefaultFontMapper mapper = new DefaultFontMapper();
            PdfGraphics2D g2d = new PdfGraphics2D(canvas, width, height, mapper);
            int x = xStart;
            int y = yStart;
            g2d.drawString(scaffoldName+"\t\tLength="+String.valueOf(scaffoldLength), x, y-25);
            for (int j = 0; j < chrNum; j++) {
                g2d.setColor(Color.black);
                g2d.drawString(String.valueOf(j+1), x, y-5);
                g2d.setColor(Color.ORANGE);
                g2d.fillRect(x, y, chrLength[j], recHeight);
                x+=chrLength[0];
                x+=xInterval;
            }
            x = xStart;
            y+=recHeight; int y1 = y;
            y+=1.8*yInterval; int y2 = y;
            g2d.setColor(Color.cyan);
            g2d.fillRect(x, y, contigRecLength, recHeight);
            
            g2d.setColor(Color.BLACK);
            y+=recHeight+2;
            g2d.drawLine(x, y, x+contigRecLength, y);
            for (int i = 0; i < mark.length; i++) {
                g2d.drawLine(mark[i], y, mark[i], y+3);
                g2d.drawString(String.valueOf(i)+" Mb", mark[i], y+16);
            }
            g2d.setColor(Color.cyan);
            for (int i = 0; i < bzPos.length; i++) {
                g2d.drawLine(bzPos[i], y, bzPos[i], y+3);
            }
            
            Color c = new Color (0,0,255,20);
            g2d.setColor(c);
            for (int i = 0; i < anchorNum; i++) {
                g2d.drawLine(anchorGPos[i], y1, anchorScaffoldPos[i], y2);
            }
            g2d.dispose();
            pd.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void processMafForDot () {
        String infileS = "M:\\collaboration\\w22\\W22RegionOnW22Bac.maf.txt";
        String outfileS = "M:\\collaboration\\w22\\W22RegionOnW22Bac.maf.dot.txt";
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("W22_BAC\tW22_Genome_Assembly");
            bw.newLine();
            String temp;
            while ((temp = br.readLine())!= null) {
                if (temp.contains("score=")) {
                    int score = Integer.valueOf(temp.split("=")[1]);
                    String ts = br.readLine();
                    String qs = br.readLine();
                    String[] tem = ts.split("\\s+");
                    int hitStart = Integer.valueOf(tem[2]);
                    int hitLength = Integer.valueOf(tem[3]);
                    String hitStrand = tem[4];
                    tem = qs.split("\\s+");
                    int queryStart = Integer.valueOf(tem[2]);
                    int queryLength = Integer.valueOf(tem[3]);
                    String queryStrand = tem[4];
                    
                    StringBuilder sb = new StringBuilder();
                    if (queryStrand.equals("+")) {
                        sb.append(hitStart+1).append("\t").append(queryStart+1);
                        bw.write(sb.toString());
                        bw.newLine();
                        sb = new StringBuilder();
                        sb.append(hitStart+hitLength).append("\t").append(queryStart+queryLength);
                    }
                    else {
                        queryStart = 238141-queryStart-queryLength;
                        sb.append(hitStart+1).append("\t").append(queryStart+queryLength);
                        bw.write(sb.toString());
                        bw.newLine();
                        sb = new StringBuilder();
                        sb.append(hitStart+hitLength).append("\t").append(queryStart+1);
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                    bw.write("NA\tNA");
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
    
    public void processMaf () {
        String infileS = "M:\\collaboration\\w22\\W22RegionOnW22Bac.maf.txt";
        String outfileS = "M:\\collaboration\\w22\\W22RegionOnW22Bac.maf.processed.txt";
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Score\tHit\tHitStart\tHitLength\tHitStrand\tQuery\tQueryStart\tQueryLength\tQueryStrand");
            bw.newLine();
            String temp;
            while ((temp = br.readLine())!= null) {
                if (temp.contains("score=")) {
                    int score = Integer.valueOf(temp.split("=")[1]);
                    String ts = br.readLine();
                    String qs = br.readLine();
                    String[] tem = ts.split("\\s+");
                    String hit = tem[1];
                    int hitStart = Integer.valueOf(tem[2]);
                    int hitLength = Integer.valueOf(tem[3]);
                    String hitStrand = tem[4];
                    tem = qs.split("\\s+");
                    String query = tem[1];
                    int queryStart = Integer.valueOf(tem[2]);
                    int queryLength = Integer.valueOf(tem[3]);
                    String queryStrand = tem[4];
                    StringBuilder sb = new StringBuilder();
                    sb.append(score).append("\t").append(hit).append("\t").append(hitStart).append("\t").append(hitLength).append("\t").append(hitStrand).append("\t").append(query).append("\t").append(queryStart).append("\t").append(queryLength).append("\t").append(queryStrand);
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
    
    public void getBz1Region () {
        String infileS = "M:\\collaboration\\w22\\bz1W22Scaffold.fa";
        String outfileS = "M:\\collaboration\\w22\\bz1W22Region.fa";
        int[] pos = {3482753,3736302};
        Fasta f = new Fasta (infileS);
        String s = new Sequence(f.getSeq(0).substring(pos[0]-1, pos[1])).getReverseComplementarySeq();
        String[] names = {"bz1W22Region"};
        String[] seqs = {s};
        int[] ids = {1};
        Fasta fr = new Fasta(names, seqs, ids);
        fr.writeFasta(outfileS);
    }
    
    public void getBz1Scaffold () {
        String sName = "scaffold242506";
        String outfileS = "M:\\collaboration\\w22\\bz1W22Scaffold.fa";
        Fasta f = new Fasta (this.scaffoldFileS);
        f.sortRecordByLengthDescending();
        try {
            
            boolean[] ifOut = new boolean[f.getSeqNumber()];
            for (int i = 0; i < f.getSeqNumber(); i++) {
                if (f.getName(i).equals(sName)) {
                    ifOut[i] = true;
                    break;
                }
            }
            f.writeFasta(outfileS, ifOut);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    public void getGenomeStatistics () {
        String outfileS = "M:\\collaboration\\w22\\scaffoldStatistics.txt";
        Fasta f = new Fasta (this.scaffoldFileS);
        f.sortRecordByLengthDescending();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Name\tLength");
            bw.newLine();
            for (int i = 0; i < f.getSeqNumber(); i++) {
                bw.write(f.getName(i)+"\t"+String.valueOf(f.getSeqLength(i)));
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
