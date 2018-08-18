/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package format;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import utils.IoUtils;

/**
 * Holding info from pairwise alignment result, including tools to convert results from aligners (lastz and blast)
 * @author Fei Lu
 */
public class PairwiseAlignment {
    PWAlignmentInfo[] ai = null;
    /**
     * sortByQueryAndScore is 0, sortByQueryAndScore is 1
     */
    int sortType=0;
    public PairwiseAlignment () {}
    
    /**
     * Return the first index of query, if query does not exist, return -1; if the list is not sorted, return min value
     * @param query
     * @return 
     */
    public int getFirstIndexOfQuery (String query) {
        if (sortType!=0) return Integer.MIN_VALUE;
        PWAlignmentInfo q = new PWAlignmentInfo(query, Integer.MAX_VALUE);
        int index = Arrays.binarySearch(ai, q, new SortByQueryAndScore());
        index = -index -1;
        if (ai[index].query.equals(query)) return index;
        else return -1;
    }
    
    /**
     * Return the last index of query, exclusive, if query does not exist, return -1; if the list is not sorted, return min value
     * @param query
     * @return 
     */
    public int getLastIndexOfQuery (String query) {
        if (sortType!=0) return Integer.MIN_VALUE;
        PWAlignmentInfo q = new PWAlignmentInfo(query, 0);
        int index = Arrays.binarySearch(ai, q, new SortByQueryAndScore());
        index = -index -1;
        if (ai[index-1].query.equals(query)) return index;
        else return -1;
    }
    
    /**
     * Return the number of alignment of a query, if query does not exist, return -1; if the list is not sorted, return min value
     * @param query
     * @return 
     */
    public int getNumberOfAlignmentOfQuery (String query) {
        int first = this.getFirstIndexOfQuery(query);
        int last = this.getLastIndexOfQuery(query);
        if (first == Integer.MIN_VALUE) return Integer.MIN_VALUE;
        else if (first == -1) return -1;
        return last-first;
    }
    
    public String[] getQueries () {
        HashSet<String> hs = new HashSet();
        for (int i = 0; i < this.getNumberOfAlignment(); i++) {
            hs.add(this.getQuery(i));
        }
        String[] qs = hs.toArray(new String[hs.size()]);
        Arrays.sort(qs);
        return qs;
    }
    
    public int getNumberOfQueries () {
        return this.getQueries().length;
    }
    
    public void writeTxtPairwiseAlignment (String alignmentFileS) {
        try {
            BufferedWriter bw = IoUtils.getTextWriter(alignmentFileS);
            bw.write("Query\tHit\tQueryStart\tQueryEnd\tHitStart\tHitEnd\tStrand\tScore\tEvalue");
            bw.newLine();
            for (int i = 0; i < this.getNumberOfAlignment(); i++) {
                bw.write(ai[i].getOutputString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Merge another PairwiseAlignment to this one
     * @param another 
     */
    public void mergeAlignment (PairwiseAlignment another) {
        ArrayList<PWAlignmentInfo> anotherList = new ArrayList<PWAlignmentInfo>(Arrays.asList(another.ai));
        ArrayList<PWAlignmentInfo> thisList = new ArrayList<PWAlignmentInfo>(Arrays.asList(ai));
        thisList.addAll(anotherList);
        ai = thisList.toArray(new PWAlignmentInfo[thisList.size()]);
        this.sortByQueryAndScore();
        System.out.println("Another PairwiseAlignment is merged");
    }
    
    public int getNumberOfAlignment () {
        return this.ai.length;
    }
    
    public String getQuery (int index) {
        return this.ai[index].query;
    }
    
    public String getHit (int index) {
        return this.ai[index].hit;
    }
    
    public int getQueryStart (int index) {
        return this.ai[index].qStart;
    }
    
    public int getQueryEnd (int index) {
        return this.ai[index].qEnd;
    }
    
    public int getHitStart (int index) {
        return this.ai[index].hStart;
    }
    
    public int getHitEnd (int index) {
        return this.ai[index].hEnd;
    }
    
    public int getScore (int index) {
        return this.ai[index].score;
    }
    
    public byte getStrand (int index) {
        return this.ai[index].strand;
    }
    
    public double getEvalue (int index)  {
        return this.ai[index].evalue;
    }
    
    public void readFromLastzAXF (String alignmentFileS) {
        ArrayList<PWAlignmentInfo> aList = new ArrayList();
        System.out.println("Reading SAM format alignment (Bowtie2) from: " + alignmentFileS);
        String temp = null;
        int cnt = 0;
        try {
            BufferedReader br = IoUtils.getTextReader(alignmentFileS);
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) continue;
                aList.add(this.getAlignmentInfoFromLastzAXF(temp));
                br.readLine();br.readLine();br.readLine();
                cnt++;
            }
        }
        catch(Exception e) {
            e.printStackTrace();
        }
        System.out.println(String.valueOf(cnt)+" alignment imported from " + alignmentFileS);
        ai = aList.toArray(new PWAlignmentInfo[aList.size()]);
        this.sortByQueryAndScore();
    }
    
    private PWAlignmentInfo getAlignmentInfoFromLastzAXF (String inputStr) {
        String[] temp = inputStr.split("\\s+");
        String query = temp[4];
        String hit = temp[1];
        int hStart = Integer.valueOf(temp[2]);
        int hEnd = Integer.valueOf(temp[2])+1;
        byte strand;
        if (temp[7].startsWith("+")) strand = 1;
        else strand = -1;
        int qStart = Integer.valueOf(temp[5]);
        int qEnd = Integer.valueOf(temp[6])+1;
        int score = Integer.valueOf(temp[8]);
        return new PWAlignmentInfo (query, hit, qStart, qEnd, hStart, hEnd, score, strand);
    }
    
    public void sortByQuery () {
        Arrays.sort(ai, new SortByQuery());
        this.sortType = 1;
    }
    
    private class SortByQuery implements Comparator <PWAlignmentInfo> {
        @Override
	public int compare (PWAlignmentInfo o1, PWAlignmentInfo o2) {
		return o1.query.compareTo(o2.query);
	}
    }
    
    /**
     * Sort alignment by query and score, score is in descending order
     */
    public void sortByQueryAndScore () {
        System.out.println("Start sortByQueryAndScore");
        Arrays.sort(ai, new SortByQueryAndScore());
        this.sortType=0;
        System.out.println("Finished sortByQueryAndScore");
    }
    
    private class SortByQueryAndScore implements Comparator <PWAlignmentInfo> {
        @Override
	public int compare (PWAlignmentInfo o1, PWAlignmentInfo o2) {
                if (o1.query.equals(o2.query)) {
                    return o2.score-o1.score;
                }
                else {
                    return o1.query.compareTo(o2.query);
                }
		
	}
    }
    
    private class PWAlignmentInfo {
        String query = "";
        String hit = "";
        //strand doesn't matter for qStart and qEnd, no switch
        int qStart = Integer.MIN_VALUE;
        int qEnd = Integer.MIN_VALUE;
        int hStart = Integer.MIN_VALUE;
        int hEnd = Integer.MIN_VALUE;
        int score = Integer.MIN_VALUE;
        //forward is 1, backward is -1
        byte strand = Byte.MIN_VALUE;
        double evalue = 1;
        
        public PWAlignmentInfo (String query, int score) {
            this.query = query;
            this.score = score;
        }
        
        public PWAlignmentInfo (String query, String hit, int qStart, int qEnd, int hStart, int hEnd, int score, byte strand) {
            this.query = query;
            this.hit = hit;
            this.qStart = qStart;
            this.qEnd = qEnd;
            this.hStart = hStart;
            this.hEnd = hEnd;
            this.score = score;
            this.strand = strand;
        }
        
        public PWAlignmentInfo (String query, String hit, int qStart, int qEnd, int hStart, int hEnd, int score, byte strand, double evalue) {
            this(query, hit, qStart,  qEnd,  hStart, hEnd, score,  strand);
            this.evalue = evalue;
        }
        
        public String getOutputString () {
            StringBuilder sb = new StringBuilder();
            sb.append(query).append("\t").append(hit).append("\t").append(qStart).append("\t").append(qEnd).append("\t").append(hStart).append("\t").append(hEnd).append("\t");
            sb.append(strand).append("\t").append(score).append("\t").append(evalue);
            return sb.toString();
        }

    }
}
