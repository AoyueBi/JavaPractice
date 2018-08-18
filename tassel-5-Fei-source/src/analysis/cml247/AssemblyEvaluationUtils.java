/*
* To change this license header, choose License Headers in Project Properties.
* To change this template file, choose Tools | Templates
* and open the template in the editor.
*/

package analysis.cml247;

import analysis.deprecated.panAnchor.Mo17Utils;
import com.itextpdf.awt.DefaultFontMapper;
import com.itextpdf.awt.PdfGraphics2D;
import com.itextpdf.text.Document;
import com.itextpdf.text.Rectangle;
import com.itextpdf.text.pdf.PdfContentByte;
import com.itextpdf.text.pdf.PdfWriter;
import format.Bins;
import format.Fasta;
import format.ShortreadAlignment;
import format.Table;
import format.Vertices;
import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeSet;
import net.maizegenetics.dna.map.TagsOnGeneticMap;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import utils.IOFileFormat;
import utils.IoUtils;
import utils.StringArrayUtils;

/**
 *
 * @author fl262
 */
public class AssemblyEvaluationUtils {
    
    public AssemblyEvaluationUtils () {}
    
    public void mkPathRepeatValue (String vertexRepeatPathFileS, String vertexPathRepeatValueFileS) {
        Table t = new Table (vertexRepeatPathFileS);
        TreeSet<Integer> l = new TreeSet();
        for (int i = 0; i < t.getRowNumber(); i++) {
            l.add(Integer.valueOf(t.content[i][2]));
        }
        Integer[] path = l.toArray(new Integer[l.size()]);
        int[] countA = new int[path.length];
        int[] count1 = new int[path.length];
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = Arrays.binarySearch(path, Integer.valueOf(t.content[i][2]));
            countA[index]++;
            count1[index] += Integer.valueOf(t.content[i][1]);
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(vertexPathRepeatValueFileS);
            bw.write("PathNum\tRepeatValue\tNumRepeatSequence\tNumSequence");
            bw.newLine();
            for (int i = 0; i < path.length; i++) {
                bw.write(String.valueOf(path[i])+"\t"+String.valueOf((double)count1[i]/countA[i])+"\t"+String.valueOf(count1[i])+"\t"+String.valueOf(countA[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkVertexRepeatPathNum (String samFileS, String vertexPathFileS, String vertexRepeatPathFileS) {
        HashMap<String, String> vpMap = new HashMap ();
        Table t = new Table (vertexPathFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            vpMap.put(t.content[i][0], t.content[i][3]);
        }
        try {
            BufferedReader br = IoUtils.getTextReader(samFileS);
            BufferedWriter bw = IoUtils.getTextWriter(vertexRepeatPathFileS);
            bw.write("VID_Index\tRepeat\tPathNum");
            bw.newLine();
            String temp;
            while (!(temp = br.readLine()).startsWith("@PG")) {}
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.substring(0, 30).split("\t");
                String name = tem[0].split("_")[0];
                if (tem[1].equals("4")) {
                    bw.write(tem[0]+"\t0\t"+vpMap.get(name)+"\t");
                }
                else {
                    bw.write(tem[0]+"\t1\t"+vpMap.get(name)+"\t");
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }
    
    public void sliceMaizeContig2 (String sortedFastaFileS, String sliceFastaFileS, String maizeSorghumRatioFileS, int sliceSize, int seqLengthCut, float ratioCut) {
        Table t = new Table (maizeSorghumRatioFileS);
        ArrayList<Integer> maizeList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (Float.valueOf(t.content[i][4]) < ratioCut) continue;
            maizeList.add(Integer.valueOf(t.content[i][0]));
        }
        Integer[] maizeIndices = maizeList.toArray(new Integer[maizeList.size()]);
        Arrays.sort(maizeIndices);
        String temp;
        String name;
        String seq;
        int cnt = 0;
        try {
            BufferedReader br = IoUtils.getTextReader(sortedFastaFileS);
            BufferedWriter bw = new BufferedWriter (new FileWriter(sliceFastaFileS), 65536);
            
            while ((temp = br.readLine()) != null) {
                name = temp.replaceFirst(">", "");
                if (Arrays.binarySearch(maizeIndices, Integer.valueOf(name)) < 0) {
                    br.readLine();
                    continue;
                }
                seq = br.readLine();
                int len = seq.length();
                if (seq.length() < seqLengthCut) continue;
                int left = len%sliceSize;
                int n = len/sliceSize;
                if (left != 0) n++;
                int currentIndex = 0;
                for (int j = 0; j < n; j++) {
                    int actSize = sliceSize;
                    if (j == n-1 && left!=0) continue;
                    bw.write(">"+name+"_"+String.valueOf(j));
                    bw.newLine();
                    bw.write(seq.substring(currentIndex, currentIndex+actSize));
                    bw.newLine();
                    currentIndex+=actSize;
                }
                if (cnt%1000000 == 0) System.out.println("Sliced " + String.valueOf(cnt) + " sequences");
                cnt++;
            }
            bw.flush();
            bw.close();
            System.out.println("Sliced sequences are in " + sliceFastaFileS);
        }
        catch (Exception e) {
            System.out.println(cnt);
            e.printStackTrace();
        }
    }
    
    public void sliceMaizeContig (String sortedFastaFileS, String sliceFastaFileS, String maizeSorghumRatioFileS, int sliceSize, int seqLengthCut, float ratioCut) {
        Table t = new Table (maizeSorghumRatioFileS);
        ArrayList<Integer> maizeList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (Float.valueOf(t.content[i][4]) < ratioCut) continue;
            maizeList.add(Integer.valueOf(t.content[i][0]));
        }
        Integer[] maizeIndices = maizeList.toArray(new Integer[maizeList.size()]);
        Arrays.sort(maizeIndices);
        Fasta f = new Fasta (sortedFastaFileS);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(sliceFastaFileS), 65536);
            for (int i = 0; i < f.getSeqNumber(); i++) {
                int len = f.getSeqLength(i);
                if (len < seqLengthCut) continue;
                if (Arrays.binarySearch(maizeIndices, Integer.valueOf(f.getName(i))) < 0) continue;
                int left = len%sliceSize;
                int n = len/sliceSize;
                if (left != 0) n++;
                int currentIndex = 0;
                for (int j = 0; j < n; j++) {
                    int actSize = sliceSize;
                    if (j == n-1 && left!=0) continue;
                    bw.write(">"+f.getName(i)+"_"+String.valueOf(j));
                    bw.newLine();
                    bw.write(f.getSeq(i).substring(currentIndex, currentIndex+actSize));
                    bw.newLine();
                    currentIndex+=actSize;
                }
                if (i%1000000 == 0) System.out.println("Sliced " + String.valueOf(i) + " sequences");
            }
            bw.flush();
            bw.close();
            System.out.println("Sliced sequences are in " + sliceFastaFileS);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void convertDiscovarSrcFileS (String srcFileS, String edgeSrcFileS, String indexMapFileS) {
        String temp = null;
        String[] tem = null;
        String[] te = null;
        try {
            HashMap<Integer, Integer> edgeMap = new HashMap();
            BufferedReader br = IoUtils.getTextReader(indexMapFileS);
            temp = br.readLine();
            while ((temp = br.readLine()) != null) {
                tem = temp.split("\t");
                edgeMap.put(Integer.valueOf(tem[1].split(" ")[0]), Integer.valueOf(tem[0]));
            }
            br = IoUtils.getTextReader(srcFileS);
            BufferedWriter bw = IoUtils.getTextWriter(edgeSrcFileS);
            while ((temp = br.readLine()) != null) {
                if (temp.indexOf(",") == -1) {
                    bw.write(String.valueOf(edgeMap.get(Integer.valueOf(temp))));
                    bw.newLine();
                    continue;
                }
                tem = temp.split("\\D+");
                te = temp.split("\\d+");
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < tem.length-1; i++) {
                    sb.append(edgeMap.get(Integer.valueOf(tem[i]))).append(te[i+1]);
                }
                sb.append(edgeMap.get(Integer.valueOf(tem[tem.length-1])));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.close();
            br.close();
        }
        catch (Exception e) {
            System.out.println(temp);
            e.printStackTrace();
        }
    }
    
    public void selectVertexPathLength (String vertexPathLengthFileS, int size, String subVertexPathLengthFileS) {
        Table t = new Table (vertexPathLengthFileS);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(subVertexPathLengthFileS);
            bw.write("VID\tPathNum\tLength");
            bw.newLine();
            for (int i = 0; i < size; i++) {
                int index = (int)(Math.random()*t.getRowNumber());
                if (t.content[index][1].equals("null")) continue;
                bw.write(t.content[index][0]+"\t"+t.content[index][1]+"\t"+t.content[index][2]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
        
    }
    
    public void mkVertexPathLength (String vertexPathFileS, String indexMapFileS, String vertexPathLengthFileS) {
        HashMap<String, String> vpMap = new HashMap ();
        Table t = new Table (vertexPathFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            vpMap.put(t.content[i][0], t.content[i][3]);
        }
        t = null;
        t = new Table (indexMapFileS);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(vertexPathLengthFileS);
            bw.write("VID\tPathNum\tLength");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                bw.write(t.content[i][0]+"\t"+vpMap.get(t.content[i][0])+"\t"+t.content[i][2]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkVertexDegreeFileS (String graphFileS, String vertexDegreeFileS, int degreeCut) {
        Vertices vs = new Vertices (graphFileS, IOFileFormat.Text);
        vs.sortByGIndex();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(vertexDegreeFileS);
            bw.write("VID\tInDegree\tOutDegree\tAllDegree\tGIndex");
            bw.newLine();
            int VID;
            int numIn;
            int numOut;
            int numAll;
            int gIndex;
            for (int i = 0; i < vs.getVertexNum(); i++) {
                VID = vs.getVID(i);
                numIn = vs.getInDegreeNum(i);;
                numOut  = vs.getOutDegreeNum(i);
                numAll = numIn+numOut;
                gIndex = vs.getGIndex(i);
                if (numAll < degreeCut) continue;
                bw.write(String.valueOf(VID)+"\t"+String.valueOf(numIn)+"\t"+String.valueOf(numOut)+"\t"+String.valueOf(numAll)+"\t"+String.valueOf(gIndex));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Max in degree: " + String.valueOf(vs.getMaxInDegreeNum()));
        System.out.println("Max out degree: " + String.valueOf(vs.getMaxOutDegreeNum()));
    }
    
    public void mkVertexPathFileS (String srcFileS, String vertexPathFileS) {
        String temp = null;
        String[] tem = null;
        String[] te;
        try {
            BufferedReader br = IoUtils.getTextReader(srcFileS);
            BufferedWriter bw = IoUtils.getTextWriter(vertexPathFileS);
            bw.write("VID\tGIndex\tVNumInG\tNumPath");
            bw.newLine();
            TreeSet<Integer> vertexSet = new TreeSet();
            Integer[] vertexArray = null;
            int[] pathNum = null;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                vertexSet = new TreeSet();
                vertexArray = null;
                temp = temp.replaceAll(",\\{\\{}}", "");
                tem = temp.replaceAll("\\{|}", "").split(",");
                for (int i = 0; i < tem.length; i++) {
                    vertexSet.add(Integer.valueOf(tem[i]));
                }
                vertexArray = vertexSet.toArray(new Integer[vertexSet.size()]);
                pathNum = new int[vertexArray.length];
                for (int i = 0; i < pathNum.length; i++) pathNum[i] = 1;
                int leftBraceIndex = Integer.MIN_VALUE;
                int rightBraceIndex = Integer.MIN_VALUE;
                while ((leftBraceIndex = temp.indexOf(",{{")) != -1) {
                    rightBraceIndex = temp.indexOf("}},");
                    String inBraceStr = temp.substring(leftBraceIndex+3, rightBraceIndex);
                    tem = inBraceStr.split("\\D+");
                    TreeSet<Integer> set = new TreeSet();
                    for (int i = 0; i < tem.length; i++) {
                        set.add(Integer.valueOf(tem[i]));
                    }
                    Integer[] setA = set.toArray(new Integer[set.size()]);
                    int[] setCount = new int[setA.length];
                    for (int i = 0; i < tem.length; i++) {
                        int index = Arrays.binarySearch(setA, Integer.valueOf(tem[i]));
                        setCount[index]++;
                    }
                    int nEdge = inBraceStr.split("},\\{").length;
                    for (int i = 0; i < setA.length; i++) {
                        setCount[i] = nEdge-setCount[i]+1;
                        int index = Arrays.binarySearch(vertexArray, setA[i]);
                        pathNum[index] = setCount[i];
                    }
                    temp = temp.substring(rightBraceIndex+3, temp.length());
                }
                for (int i = 0; i < vertexArray.length; i++) {
                    bw.write(String.valueOf(vertexArray[i])+"\t"+String.valueOf(cnt)+"\t"+String.valueOf(vertexArray.length)+"\t"+String.valueOf(pathNum[i]));
                    bw.newLine();
                }
                cnt++;
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void convertDiscovarToGraph (String discovarFileS, String graphFileS) {
        String temp = null;
        String[] tem = null;
        try {
            BufferedReader br = IoUtils.getTextReader(discovarFileS);
            BufferedWriter bw = IoUtils.getTextWriter(graphFileS);
            
            TreeSet<Integer> vertexSet = new TreeSet();
            Integer[] vertexArray = null;
            TreeSet<Integer>[] inDegrees = null;
            TreeSet<Integer>[] outDegrees = null;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                vertexSet = new TreeSet();
                vertexArray = null;
                inDegrees = null;
                outDegrees = null;
                temp = temp.replaceAll(",\\{\\{}}", "");
                tem = temp.replaceAll("\\{|}", "").split(",");
                for (int i = 0; i < tem.length; i++) {
                    vertexSet.add(Integer.valueOf(tem[i]));
                }
                vertexArray = vertexSet.toArray(new Integer[vertexSet.size()]);
                inDegrees = new TreeSet[vertexArray.length];
                outDegrees = new TreeSet[vertexArray.length];
                for (int i = 0; i < vertexArray.length; i++) {
                    inDegrees[i] = new TreeSet();
                    outDegrees[i] = new TreeSet();
                }
                int leftBraceIndex = Integer.MIN_VALUE;
                int rightBraceIndex = Integer.MIN_VALUE;
                String[] left;
                String[] right;
                String leftV;
                String rightV;
                boolean ifFirst = true;
                String leftS;
                String[] pVArray = null;
                while ((leftBraceIndex = temp.indexOf(",{{")) != -1) {
                    rightBraceIndex = temp.indexOf("}},");
                    String inBraceStr = temp.substring(leftBraceIndex+3, rightBraceIndex);
                    tem = inBraceStr.split("},\\{");
                    int bondaryIndex = leftBraceIndex-25;
                    if (bondaryIndex < 0) bondaryIndex = 0;
                    left = temp.substring(bondaryIndex, leftBraceIndex).split(",");
                    leftV = left[left.length-1].replaceAll("\\{|}", "");
                    bondaryIndex = rightBraceIndex+30;
                    if (bondaryIndex > temp.length()) bondaryIndex = temp.length();
                    right = temp.substring(rightBraceIndex+3, bondaryIndex).split(",");
                    rightV = right[0].replaceAll("\\{|}", "");
                    String[] leftVArray = new String[tem.length];
                    String[] rightVArray = new String[tem.length];
                    for (int i = 0; i < tem.length; i++) {
                        String[] te = tem[i].split(",");
                        leftVArray[i] = te[0];
                        rightVArray[i] = te[te.length-1];
                        this.addVertexFromDiscovar(vertexArray, inDegrees, outDegrees, te);
                        this.addInDegree(vertexArray, inDegrees, leftVArray[i], leftV);
                        this.addOutDegree(vertexArray, outDegrees, rightVArray[i], rightV);
                    }
                    if (ifFirst) {
                        leftS = temp.substring(0, leftBraceIndex);
                        ifFirst = false;
                        String[] te = leftS.split(",");
                        this.addVertexFromDiscovar(vertexArray, inDegrees, outDegrees, te);
                        this.addInDegree(vertexArray, inDegrees, te[0], null);
                        for (int i = 0; i < leftVArray.length; i++) {
                            this.addOutDegree(vertexArray, outDegrees, te[te.length-1], leftVArray[i]);
                        }
                    }
                    else {
                        leftS = temp.substring(0, leftBraceIndex);
                        String[] te = leftS.split(",");
                        this.addVertexFromDiscovar(vertexArray, inDegrees, outDegrees, te);
                        for (int i = 0; i < pVArray.length; i++) {
                            this.addInDegree(vertexArray, inDegrees, te[0], pVArray[i]);
                        }
                        for (int i = 0; i < leftVArray.length; i++) {
                            this.addOutDegree(vertexArray, outDegrees, te[te.length-1], leftVArray[i]);
                        }
                    }
                    pVArray = rightVArray;
                    temp = temp.substring(rightBraceIndex+3, temp.length());
                }
                if (ifFirst) {
                    String[] te = temp.split(",");
                    this.addVertexFromDiscovar(vertexArray, inDegrees, outDegrees, te);
                    this.addInDegree(vertexArray, inDegrees, te[0], null);
                    this.addOutDegree(vertexArray, outDegrees, te[te.length-1], null);
                }
                else {
                    String[] te = temp.split(",");
                    this.addVertexFromDiscovar(vertexArray, inDegrees, outDegrees, te);
                    for (int i = 0; i < pVArray.length; i++) {
                        this.addInDegree(vertexArray, inDegrees, te[0], pVArray[i]);
                    }
                    this.addOutDegree(vertexArray, outDegrees, te[te.length-1], null);
                }
                Integer[] inDegreeArray = null;
                Integer[] outDegreeArray = null;
                for (int i = 0; i < vertexArray.length; i++) {
                    bw.write("GIndex:\t"+String.valueOf(cnt));
                    bw.newLine();
                    bw.write("VID:\t"+String.valueOf(vertexArray[i]));
                    bw.newLine();
                    inDegreeArray = inDegrees[i].toArray(new Integer[inDegrees[i].size()]);
                    outDegreeArray = outDegrees[i].toArray(new Integer[outDegrees[i].size()]);
                    for (int j = 0; j < inDegreeArray.length-1; j++) {
                        bw.write(String.valueOf(inDegreeArray[j])+"\t");
                    }
                    if (inDegreeArray.length > 0)bw.write(String.valueOf(inDegreeArray[inDegreeArray.length-1])+"\t");
                    bw.newLine();
                    for (int j = 0; j < outDegreeArray.length-1; j++) {
                        bw.write(String.valueOf(outDegreeArray[j])+"\t");
                    }
                    if (outDegreeArray.length > 0)bw.write(String.valueOf(outDegreeArray[outDegreeArray.length-1])+"\t");
                    bw.newLine();
                }
                cnt++;
                if (cnt%100000 == 0) System.out.println("Converted " + String.valueOf(cnt) + " graphs");
            }
            
            bw.flush();
            bw.close();
            br.close();
        }
        catch (Exception e) {
            System.out.println(temp);
            e.printStackTrace();
        }
    }
    
    private void addOutDegree (Integer[] vertexArray, TreeSet<Integer>[] outDegrees, String vertexS, String vertexOut) {
        if (vertexOut == null) return;
        int query = Integer.valueOf(vertexS);
        int index = Arrays.binarySearch(vertexArray, query);
        outDegrees[index].add(Integer.valueOf(vertexOut));
    }
    
    private void addInDegree (Integer[] vertexArray, TreeSet<Integer>[] inDegrees, String vertexS, String vertexIn) {
        if (vertexIn == null) return;
        int query = Integer.valueOf(vertexS);
        int index = Arrays.binarySearch(vertexArray, query);
        inDegrees[index].add(Integer.valueOf(vertexIn));
    }
    
    private void addVertexFromDiscovar (Integer[] vertexArray, TreeSet<Integer>[] inDegrees, TreeSet<Integer>[] outDegrees, String[] vertices) {
        for (int i = 0; i < vertices.length-1; i++) {
            int query = Integer.valueOf(vertices[i]);
            int nextQuery = Integer.valueOf(vertices[i+1]);
            int index = Arrays.binarySearch(vertexArray, query);
            int nextIndex = Arrays.binarySearch(vertexArray, nextQuery);
            outDegrees[index].add(nextQuery);
            inDegrees[nextIndex].add(query);
        }
    }
    
    public void mkNodeDegreeFileS (String nodeFileS, String nodeDegreeFileS) {
        try {
            BufferedReader br = IoUtils.getTextReader(nodeFileS);
            BufferedWriter bw = IoUtils.getTextWriter(nodeDegreeFileS);
            String temp;
            bw.write("NodeMultiDegree");
            bw.newLine();
            while ((temp = br.readLine()) != null) {
                ArrayList<Integer> l = new ArrayList();
                l = this.getDegreeList(temp, l, 10);
                for (int i = 0; i < l.size(); i++) {
                    bw.write(String.valueOf(l.get(i)));
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
    
    private ArrayList<Integer> getDegreeList (String record, ArrayList<Integer> l, int nMax) {
        for (int i = 0; i < nMax; i++) {
            String query= "";
            for (int j = 0; j < nMax-i; j++) {
                query = query +"{";
            }
            int startIndex = record.indexOf(query);
            if (startIndex!=-1) {
                String hit = "";
                for (int j = 0; j < query.length(); j++) {
                    hit = hit+"}";
                }
                int endIndex = record.indexOf(hit)+hit.length();
                String subRecord = record.substring(startIndex+1, endIndex-1);
                String restRecord = record.substring(endIndex, record.length());
                int d = this.getDegree(subRecord, query.length()-1);
                if (d>1) l.add(d);
                l = this.getDegreeList(restRecord, l, query.length());
                break;
            }
        }
        return l;
    }
    
    private int getDegree (String subRecord, int n) {
        String query = "";
        String hit = "";
        for (int i = 0; i < n; i++) {
            query = query+"{";
            hit = hit+"}";
        }
        int startIndex = 0;
        int endIndex;
        int cnt = 0;
        while (startIndex > -1) {
            startIndex = subRecord.indexOf(query);
            endIndex = subRecord.indexOf(hit);
            if (startIndex>-1) cnt++;
            else break;
            subRecord = subRecord.substring(endIndex+n, subRecord.length());
        }
        return cnt;
    }
    
    public void mkStatisticsAssemblyByChr (String anchorOnContigFileS, String qualityFileS, String indexMapFileS) {
        Table t = new Table(indexMapFileS);
        HashMap<String,String> indexNameMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            indexNameMap.put(t.content[i][0], t.content[i][1]);
        }
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
            bw.write("Index\tScaffold\tLength\tAnchorNumber\tRatioAnchorOnMainChr");
            bw.newLine();
            for (int i = 0; i < contigNum; i++) {
                bw.write(contigName[i]+"\t"+indexNameMap.get(contigName[i])+"\t"+String.valueOf(contigLength[i])+"\t"+String.valueOf(contigAnchorNum[i])+"\t");
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
     * @param indexMapFileS
     * @param anchorFileS
     * @param anchorOnContigFileS
     * @param satisticsFileS
     */
    public void mkAnchorOnContigFile (String alignFileS, String indexMapFileS, String anchorFileS, String anchorOnContigFileS, String satisticsFileS) {
        Table t = new Table(indexMapFileS);
        HashMap<String,Integer> contigLengthMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            contigLengthMap.put(t.content[i][0], Integer.valueOf(t.content[i][2]));
        }
        TagsOnGeneticMap togm =new TagsOnGeneticMap (anchorFileS, FilePacking.Text);
        int[] gChr = new int[togm.getTagCount()];
        int[] gPos = new int[togm.getTagCount()];
        byte[] anchorLength = new byte[togm.getTagCount()];
        byte[] ifPav = new byte[togm.getTagCount()];
        try {
            for (int i = 0; i < togm.getTagCount(); i++) {
                gChr[i] = togm.getGChr(i);
                gPos[i] = togm.getGPos(i);
                anchorLength[i] = (byte)togm.getTagLength(i);
                ifPav[i] = togm.getIfPAV(i);
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
            contigLengthTotal += Long.valueOf(t.content[i][2]);
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
    
    private class ChromLength implements Comparable <ChromLength>{
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
    
    public void selectLargeContigs (String sourceFastaFileS, String desFastaFileS, String maizeSorghumRatioFileS, int minSize) {
        HashMap<String, Float> ratioMap = new HashMap();
        Table t = new Table (maizeSorghumRatioFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (Integer.valueOf(t.content[i][1]) < minSize) continue;
            ratioMap.put(">"+t.content[i][0], Float.valueOf(t.content[i][4]));
        }
        try {
            BufferedReader br = IoUtils.getTextReader(sourceFastaFileS);
            BufferedWriter bw = IoUtils.getTextWriter(desFastaFileS);
            String name;
            int cnt = 0;
            while ((name = br.readLine()) != null) {
                String seq = br.readLine();
                if (seq.length() < minSize) break;
                Float f = ratioMap.get(name);
                if (f < 0.5) continue;
                if (name == null) continue;
                if (Float.isNaN(f)) continue;
                bw.write(name);
                bw.newLine();
                bw.write(seq);
                bw.newLine();
                if (cnt%100000 == 0) System.out.println("Wrote " + String.valueOf(cnt) + " sequences");
                cnt++;
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkMaizeSorghumRatioFileS (String alignmentFileS, String indexMapFileS, String maizeSorghumRatioFileS) {
        Table t = new Table (indexMapFileS);
        HashMap<Integer,Integer> idLengthMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            idLengthMap.put(Integer.valueOf(t.content[i][0]), Integer.valueOf(t.content[i][2]));
        }
        try {
            DataInputStream dis = IoUtils.getBinaryReader(alignmentFileS);
            BufferedWriter bw = IoUtils.getTextWriter(maizeSorghumRatioFileS);
            bw.write("QueryIndex\tLength\tMaizeCnt\tSorghumCnt\tMaizeRatio");
            bw.newLine();
            dis.readByte();
            int alignNum = dis.readInt();
            String query;
            String hit;
            String[] temp;
            String current = "";
            int mCnt = 0;
            int sCnt = 0;
            for (int i = 0; i < alignNum; i++) {
                query = dis.readUTF();
                hit = dis.readUTF();
                dis.readInt(); dis.readInt(); dis.readByte(); dis.readBoolean();
                temp = query.split("_");
                if (temp[0].equals(current)) {
                    if (hit.startsWith("m")) {
                        mCnt++;
                    }
                    else if (hit.startsWith("s")) {
                        sCnt++;
                    }
                }
                else {
                    if (i==0) {
                        if (hit.startsWith("m")) {
                            mCnt++;
                        }
                        else if (hit.startsWith("s")) {
                            sCnt++;
                        }
                    }
                    else {
                        bw.write(current+"\t"+String.valueOf(idLengthMap.get(Integer.valueOf(current)))+"\t"+String.valueOf(mCnt)+"\t"+String.valueOf(sCnt)+"\t"+String.valueOf((float)mCnt/(mCnt+sCnt)));
                        bw.newLine();
                        mCnt = 0;
                        sCnt = 0;
                    }
                    current = temp[0];
                }
                if (i%5000000 == 0) System.out.println("Checked " + String.valueOf(i) + " alignments for maize sorghum ratio");
            }
            bw.write(current+"\t"+String.valueOf(idLengthMap.get(Integer.valueOf(current)))+"\t"+String.valueOf(mCnt)+"\t"+String.valueOf(sCnt)+"\t"+String.valueOf((float)mCnt/(mCnt+sCnt)));
            bw.newLine();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Maize sorghum ratio file is written at " + maizeSorghumRatioFileS);
    }
    
    public void convertSamToSimpleAlignment (String samFileS, String alignmentFileS) {
        ShortreadAlignment sa = new ShortreadAlignment();
        sa.readFromBowtie2(samFileS);
        sa.sortByQuery();
        sa.writeSimpleAlignment(alignmentFileS, IOFileFormat.Binary);
    }
    
    public void mkAssemblyStatistics (String inputFastaFileS, String statisticsFileS) {
        Fasta f = new Fasta (inputFastaFileS);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(statisticsFileS);
            bw.write(inputFastaFileS);
            bw.newLine();
            bw.write("Scaffold Number:\t" + String.valueOf(f.getSeqNumber()));
            bw.newLine();
            bw.write("Total length:\t" + String.valueOf(f.getTotalSeqLength()));
            bw.newLine();
            bw.write("L50 (bp):\t" + String.valueOf(f.getL50()));
            bw.newLine();
            bw.write("N50 Number:\t" + String.valueOf(f.getN50()));
            bw.newLine();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void selectMaizeAssembly (String sorghumRatioFileS, String inputFastaFileS, String outputFastaFileS) {
        Fasta f = new Fasta (inputFastaFileS);
        Table t = new Table (sorghumRatioFileS);
        HashMap<String, Double> valueMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            valueMap.put(t.content[i][0], Double.valueOf(t.content[i][4]));
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outputFastaFileS);
            for (int i = 0; i < f.getSeqNumber(); i++) {
                double d = valueMap.get(f.getName(i));
                if (Double.isNaN(d)) {
                    
                }
                else {
                    if (d > 0.1) continue;
                }
                bw.write(">"+f.getName(i));
                bw.newLine();
                bw.write(f.getSeq(i));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.toString();
        }
    }
    
    public void mkSorghumRatioFile (ShortreadAlignment sa, String intputFastaFileS, String sorghumRatioFileS) {
        Fasta f = new Fasta(intputFastaFileS);
        f.sortRecordByName();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(sorghumRatioFileS);
            bw.write("Scaffold\tLength\tMaizeNum\tSorghumNum\tSorghumRatio");
            bw.newLine();
            ArrayList<String> nameList = new ArrayList();
            for (int i = 0; i < sa.getAlignmentNumber(); i++) {
                if (!sa.getQuery(i).endsWith("_0")) continue;
                nameList.add(sa.getQuery(i).split("_")[0]);
            }
            String[] name = nameList.toArray(new String[nameList.size()]);
            Arrays.sort(name);
            int[] maizeCount = new int[name.length];
            int[] sorghumCount = new int[name.length];
            for (int i = 0; i < sa.getAlignmentNumber(); i++) {
                if (!sa.isMatch(i)) continue;
                int index = Arrays.binarySearch(name, sa.getQuery(i).split("_")[0]);
                if (sa.getHit(i).startsWith("maize")) {
                    maizeCount[index]++;
                }
                else {
                    sorghumCount[index]++;
                }
            }
            for (int i = 0; i < name.length; i++) {
                int index = f.getIndex(name[i]);
                int sum = sorghumCount[i]+maizeCount[i];
                double ratio;
                if (sum != 0) {
                    ratio = sorghumCount[i]/ (sorghumCount[i]+maizeCount[i]);
                }
                else {
                    ratio = Double.NaN;
                }
                bw.write(name[i]+"\t"+String.valueOf(f.getSeqLength(index))+"\t"+String.valueOf(maizeCount[i])+"\t"+String.valueOf(sorghumCount[i])+"\t"+String.valueOf(ratio));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void sliceFasta (String inputFastaFileS, String outputFastaFileS, int sliceSize) {
        Fasta f = new Fasta (inputFastaFileS);
        f.sortRecordByLengthDescending();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outputFastaFileS);
            for (int i = 0; i < f.getSeqNumber(); i++) {
                int left = f.getSeqLength(i)%sliceSize;
                int sliceNum;
                if (left == 0) {
                    sliceNum = f.getSeqLength(i)/sliceSize;
                    for (int j = 0; j < sliceNum; j++) {
                        bw.write(">"+f.getName(i)+"_"+String.valueOf(j));
                        bw.newLine();
                        bw.write(f.getSeq(i).substring(j*sliceSize, j*sliceSize+sliceSize));
                        bw.newLine();
                    }
                }
                else {
                    sliceNum = f.getSeqLength(i)/sliceSize + 1;
                    for (int j = 0; j < sliceNum-1; j++) {
                        bw.write(">"+f.getName(i)+"_"+String.valueOf(j));
                        bw.newLine();
                        bw.write(f.getSeq(i).substring(j*sliceSize, j*sliceSize+sliceSize));
                        bw.newLine();
                    }
                    bw.write(">"+f.getName(i)+"_"+String.valueOf(sliceNum-1));
                    bw.newLine();
                    bw.write(f.getSeq(i).substring((sliceNum-1)*sliceSize, (sliceNum-1)*sliceSize+left));
                    bw.newLine();
                }
                if (i%10000 == 0) {
                    System.out.println(String.valueOf(i+1)+" sequences processed");
                }
                
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    /**
     * When the leftover length < sliceSliece, descarded
     * @param sortedFastaFileS
     * @param sliceFastaFileS
     * @param sliceSize
     * @param seqLengthCut
     */
    public void sliceFasta (String sortedFastaFileS, String sliceFastaFileS, int sliceSize, int seqLengthCut) {
        Fasta f = new Fasta (sortedFastaFileS);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(sliceFastaFileS), 65536);
            for (int i = 0; i < f.getSeqNumber(); i++) {
                int len = f.getSeqLength(i);
                if (len < seqLengthCut) continue;
                int left = len%sliceSize;
                int n = len/sliceSize;
                if (left != 0) n++;
                int currentIndex = 0;
                for (int j = 0; j < n; j++) {
                    int actSize = sliceSize;
                    if (j == n-1 && left!=0) continue;
                    bw.write(">"+f.getName(i)+"_"+String.valueOf(j));
                    bw.newLine();
                    bw.write(f.getSeq(i).substring(currentIndex, currentIndex+actSize));
                    bw.newLine();
                    currentIndex+=actSize;
                }
                if (i%1000000 == 0) System.out.println("Sliced " + String.valueOf(i) + " sequences");
            }
            bw.flush();
            bw.close();
            System.out.println("Sliced sequences are in " + sliceFastaFileS);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void basicAssemblyStatistics (String indexMapFileS, String statisticFileS) {
        Table t = new Table (indexMapFileS);
        int[] lens = new int[t.getRowNumber()];
        long sum = 0;
        for (int i = 0; i < lens.length; i++) {
            lens[i] = Integer.valueOf(t.content[i][2]);
            sum+=lens[i];
        }
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(statisticFileS), 65536);
            bw.write("Contig number:\t"+String.valueOf(lens.length));
            bw.newLine();
            bw.write("Total length:\t"+String.valueOf(sum)+" bp");
            bw.newLine();
            long halfSum = sum/2;
            long current = 0;
            for (int i = 0; i < lens.length; i++) {
                current+=lens[i];
                if (lens[i] < 13000) System.out.println(String.valueOf(current)+"\t"+String.valueOf((double)current/sum));
                if (current > halfSum) {
                    bw.write("N50: "+String.valueOf(lens[i-1])+" bp");
                    break;
                }
            }
            bw.newLine();
            bw.write("Max:\t"+String.valueOf(lens[0])+" bp");
            bw.newLine();
            bw.write("Min:\t"+String.valueOf(lens[lens.length-1])+" bp");
            bw.newLine();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void renameSortFasta (String oriFastaFileS, String desFastaFileS, String indexMapFileS) {
        Fasta f = new Fasta(oriFastaFileS);
        f.sortRecordByLengthDescending();
        String[] oriNames = new String[f.getSeqNumber()];
        for (int i = 0; i < oriNames.length; i++) {
            oriNames[i] = f.getName(i);
            f.setName(String.valueOf(i), i);
        }
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(indexMapFileS), 65536);
            bw.write("IndexInNewFasta\tOriginalName\tContigLength");
            bw.newLine();
            for (int i = 0; i < oriNames.length; i++) {
                bw.write(String.valueOf(i)+"\t"+oriNames[i]+"\t"+String.valueOf(f.getSeqLength(i)));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println("IndexMap is written to " + indexMapFileS);
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        f.writeFasta(desFastaFileS);
    }
    
}
