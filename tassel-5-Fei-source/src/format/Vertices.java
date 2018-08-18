/*
* To change this license header, choose License Headers in Project Properties.
* To change this template file, choose Tools | Templates
* and open the template in the editor.
*/

package format;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 * Hold vertices information in graphs, including GraphIndex, vertex ID, in-degrees and out-degrees. There can be independent graphs. But the VID of vertex is unique.
 * Support sorting and searching methods.
 * @author fl262
 */
public class Vertices {
    int GraphNum = Integer.MIN_VALUE;
    Vertex[] v;
    int sortType = -1;
    
    public Vertices (String inputFileS, IOFileFormat format) {
        this.readVerticesFile(inputFileS, format);
    }
    
    public Vertices (Vertex[] v) {
        this.v = v;
    }
    
    /**
     * Read vertices file
     * @param inputFileS
     * @param format 
     */
    public void readVerticesFile (String inputFileS, IOFileFormat format) {
        if (format == IOFileFormat.Binary) this.readBinaryFile(inputFileS);
        else if (format == IOFileFormat.Text) this.readTextFile(inputFileS);
        else throw new UnsupportedOperationException("Not supported yet.");
        System.out.println("Vertices file read from " + inputFileS);
    }
    
    private void iniMatrix (int vertexNum) {
        v = new Vertex[vertexNum];
    }
    
    private void readTextFile (String inputFileS) {
        try {
            BufferedReader br = IoUtils.getTextReader(inputFileS);
            String temp = null;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnt++;
            }
            cnt = cnt/4;
            this.iniMatrix(cnt);
            br = IoUtils.getTextReader(inputFileS);
            String[] tem;
            int GIndex;
            int VID;
            int[] inDegree = null;
            int[] outDegree = null;
            for (int i = 0; i < cnt; i++) {
                GIndex = Integer.valueOf(br.readLine().split("\t")[1]);
                VID = Integer.valueOf(br.readLine().split("\t")[1]);
                temp = br.readLine();
                if (temp.equals("")) inDegree = new int[0];
                else {
                    tem = temp.split("\t");
                    inDegree = new int[tem.length];
                    for (int j = 0; j < tem.length; j++) {
                        inDegree[j] = Integer.valueOf(tem[j]);
                    }
                }
                temp = br.readLine();
                if (temp.equals("")) outDegree = new int[0];
                else {
                    tem = temp.split("\t");
                    outDegree = new int[tem.length];
                    for (int j = 0; j < tem.length; j++) {
                        outDegree[j] = Integer.valueOf(tem[j]);
                    }
                }
                v[i] = new Vertex(GIndex, VID, inDegree, outDegree);
                if (i%100000 == 0) System.out.println("Read " + String.valueOf(i) + " vertices");
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void readBinaryFile (String inputFileS) {
        try {
            DataInputStream dis = IoUtils.getBinaryReader(inputFileS);
            this.iniMatrix(dis.readInt());
            int GIndex;
            int VID;
            int[] inDegree = null;
            int[] outDegree = null;
            for (int i = 0; i < this.getVertexNum(); i++) {
                GIndex = dis.readInt();
                VID = dis.readInt();
                inDegree = new int[dis.readInt()];
                for (int j = 0; j < inDegree.length; j++) {
                    inDegree[j] = dis.readInt();
                }
                outDegree = new int[dis.readInt()];
                for (int j = 0; j < outDegree.length; j++) {
                    outDegree[j] = dis.readInt();
                }
                v[i] = new Vertex(GIndex, VID, inDegree, outDegree);
                if (i%100000 == 0) System.out.println("Read " + String.valueOf(i) + " vertices");
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Return the max number of in-degree
     * @return 
     */
    public int getMaxInDegreeNum () {
        int max = 0;
        for (int i = 0; i < this.getVertexNum(); i++) {
            if (this.getInDegreeNum(i) > max) max = this.getInDegreeNum(i);
        }
        return max;
    }
    
    /**
     * Return the max number of out-degree
     * @return 
     */
    public int getMaxOutDegreeNum () {
        int max = 0;
        for (int i = 0; i < this.getVertexNum(); i++) {
            if (this.getOutDegreeNum(i) > max) max = this.getOutDegreeNum(i);
        }
        return max;
    }
    
    /**
     * Return vertices of a graph
     * @param graphIndex
     * @return 
     */
    public Vertices getVerticesOfGraph (int graphIndex) {
        int startIndex = this.getStartIndexOfGraph(graphIndex);
        int endIndex = this.getEndIndexOfGraph(graphIndex);
        int size = endIndex - startIndex;
        Vertex[] nv = new Vertex[size];
        for (int i = 0; i < nv.length; i++) {
            nv[i] = v[startIndex+i];
        }
        return new Vertices(nv);
    }
    
    /**
     * Return start index of a graph, inclusive
     * @param graphIndex
     * @return 
     */
    public int getStartIndexOfGraph (int graphIndex) {
        if (this.sortType != 0) this.sortByGIndex();
        Vertex query = new Vertex (graphIndex, 0, null, null);
        int index = Arrays.binarySearch(v, query);
        if (index == 0) return index;
        else if (index < 0) {
            System.out.println("This graph index is not in the vertices data set. Program quits");
            System.exit(1);
        }
        int nextIndex;
        int step = 1000;
        while (step != 0) {
            nextIndex = index - step;
            if (nextIndex < 0) nextIndex = 0;
            if (this.getGIndex(nextIndex) == graphIndex) {
                index = nextIndex;
                if (index == 0) return index;
            }
            else if (this.getGIndex(nextIndex) > graphIndex){
                
            }
            else {
                step = step/2;
            }
        }
        return index;
    }
    
    /**
     * Return end index of a graph, exclusive
     * @param graphIndex
     * @return 
     */
    public int getEndIndexOfGraph (int graphIndex) {
        if (this.sortType != 0) this.sortByGIndex();
        Vertex query = new Vertex (graphIndex, 0, null, null);
        int index = Arrays.binarySearch(v, query);
        if (index == this.getVertexNum()-1) return index+1;
        else if (index < 0) {
            System.out.println("This graph index is not in the vertices data set. Program quits");
            System.exit(1);
        }
        int nextIndex;
        int step = 1000;
        while (step != 0) {
            nextIndex = index + step;
            if (nextIndex > this.getVertexNum()-1) nextIndex = this.getVertexNum()-1;
            if (this.getGIndex(nextIndex) == graphIndex) {
                index = nextIndex;
                if (index == this.getVertexNum()-1) return index+1;
            }
            else if (this.getGIndex(nextIndex) > graphIndex){
                step = step/2;
            }
            else {
                
            }
        }
        return index+1;
    }
    
    /**
     * Return graph index of a vertex
     * @param index
     * @return 
     */
    public int getGIndex (int index) {
        return v[index].GIndex;
    }
    
    /**
     * Return VID of a vertex
     * @param index
     * @return 
     */
    public int getVID (int index) {
        return v[index].VID;
    }
    
    /**
     * Return number of in-degree to a vertex
     * @param index
     * @return 
     */
    public int getInDegreeNum (int index) {
        return v[index].inDegree.length;
    }
    
    /**
     * Return VIDs of in-degree to a vertex
     * @param index
     * @return 
     */
    public int[] getInDegree (int index) {
        return v[index].inDegree;
    }
    
    /**
     * Return a string of VIDs of in-degree to a vertex, separated by Tab
     * @param index
     * @return 
     */
    public String getInDegreeStr (int index) {
        return v[index].getInDegreeStr();
    }
    
    /**
     * Return number of out-degree from a vertex
     * @param index
     * @return 
     */
    public int getOutDegreeNum (int index)  {
        return v[index].outDegree.length;
    }
    
    /**
     * Return VIDs of out-degree from a vertex
     * @param index
     * @return 
     */
    public int[] getOutDegree (int index) {
        return v[index].outDegree;
    }
    
    /**
     * Return a string of VIDs of out-degree from a vertex, separated by Tab
     * @param index
     * @return 
     */
    public String getOutDegreeStr (int index) {
        return v[index].getOutDegreeStr();
    }
    
    /**
     * Return number of graphs composed by these vertices
     * @return 
     */
    public int getGraphNum () {
        if (this.GraphNum == Integer.MIN_VALUE) {
            this.GraphNum = this.getGraphIndices().length;
        }
        return this.GraphNum;
    }
    
    /**
     * Return a indices list of graphs
     * @return 
     */
    public Integer[] getGraphIndices () {
        TreeSet<Integer> indexSet = new TreeSet();
        for (int i = 0; i < this.getVertexNum(); i++) {
            indexSet.add(this.getGIndex(i));
        }
        return indexSet.toArray(new Integer[indexSet.size()]);
    }
    
    /**
     * Return number of vertices
     * @return 
     */
    public int getVertexNum () {
        return this.v.length;
    }
    
    /**
     * Write vertices file
     * @param outputFileS
     * @param format 
     */
    public void writeVerticesFile (String outputFileS, IOFileFormat format) {
        if (format == IOFileFormat.Binary) this.writeBinaryFile(outputFileS);
        else if (format == IOFileFormat.Text) this.writeTextFile(outputFileS);
        else throw new UnsupportedOperationException("Not supported yet.");
        System.out.println("Vertices file written to " + outputFileS);
    }
    
    private void writeTextFile (String outputFileS) {
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outputFileS);
            for (int i = 0; i < this.getVertexNum(); i++) {
                bw.write("GIndex:\t"+String.valueOf(this.getGIndex(i)));
                bw.newLine();
                bw.write("VID:\t"+String.valueOf(this.getVID(i)));
                bw.newLine();
                bw.write(this.getInDegreeStr(i));
                bw.newLine();
                bw.write(this.getOutDegreeStr(i));
                bw.newLine();
                if (i%100000 == 0) System.out.println("Write " + String.valueOf(i) + " vertices");
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void writeBinaryFile (String outputFileS) {
        try {
            DataOutputStream dos = IoUtils.getBinaryWriter(outputFileS);
            dos.writeInt(this.getVertexNum());
            for (int i = 0; i < this.getVertexNum(); i++) {
                dos.writeInt(this.getGIndex(i));
                dos.writeInt(this.getVID(i));
                dos.writeInt(this.getInDegreeNum(i));
                for (int j = 0; j < this.getInDegreeNum(i); j++) {
                    dos.writeInt(this.getInDegree(i)[j]);
                }
                dos.writeInt(this.getOutDegreeNum(i));
                for (int j = 0; j < this.getOutDegreeNum(i); j++) {
                    dos.writeInt(this.getOutDegree(i)[j]);
                }
                if (i%100000 == 0) System.out.println("Write " + String.valueOf(i) + " vertices");
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Sort vertices by GIndex
     */
    public void sortByGIndex () {
        System.out.println("sorting by graph index");
        sortType = 0;
        Arrays.sort(v);
    }
    
    /**
     * Sort vertices by VID
     */
    public void sortByVID () {
        System.out.println("sorting by vertex ID");
        sortType = 1;
        Arrays.sort(v);
    }
    
    /**
     * Sort vertices by number of in-degree
     */
    public void sortByInDegreeSize () {
        System.out.println("sorting by in-degree size");
        sortType = 2;
        Arrays.sort(v);
    }
    
    /**
     * Sort vertices by number of out-degree
     */
    public void sortByOutDegreeSize () {
        System.out.println("sorting by out-degree size");
        sortType = 3;
        Arrays.sort(v);
    }
    
    private class Vertex implements Comparable<Vertex> {
        public int GIndex;
        public int VID;
        public int[] inDegree = null;
        public int[] outDegree = null;
        
        public Vertex (int GIndex, int ID, int[] inDegree, int[] outDegree) {
            this.GIndex = GIndex;
            this.VID = ID;
            this.inDegree = inDegree;
            this.outDegree = outDegree;
        }
        
        public String getInDegreeStr () {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < inDegree.length; i++) {
                sb.append(inDegree[i]).append("\t");
            }
            if (inDegree.length > 0) sb.deleteCharAt(sb.length()-1);
            return sb.toString();
        }
        
        public String getOutDegreeStr () {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < outDegree.length; i++) {
                sb.append(outDegree[i]).append("\t");
            }
            if (outDegree.length > 0) sb.deleteCharAt(sb.length()-1);
            return sb.toString();
        }
        
        public void sortDegree () {
            Arrays.sort(inDegree);
            Arrays.sort(outDegree);
        }
        
        @Override
        public int compareTo(Vertex o) {
            if (sortType == 0) return GIndex - o.GIndex;
            else if (sortType == 1) return VID - o.VID;
            else if (sortType == 2) return inDegree.length - inDegree.length;
            else if (sortType == 3) return outDegree.length - outDegree.length;
            else throw new UnsupportedOperationException("Not supported yet.");
        }
    }
}
    


