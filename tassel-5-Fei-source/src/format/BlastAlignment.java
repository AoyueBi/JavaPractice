/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package format;

import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class BlastAlignment {
    AlignmentInfo[] ai;
    /**0 by query (Default). 1 by hit and pos.*/
    byte sortType = 0; 
    
    /**
     * Read from blast result, outfmt 7
     * @param inputFileS
     * @param eThresh
     * @param readLength read length of fasta or fastq
     */
    public void readFromBlastTable (String inputFileS, double eThresh) {
        ArrayList<AlignmentInfo> aList = new ArrayList();
        System.out.println("Reading blast alignment (outfmt 7) from: " + inputFileS);
        try {
            BufferedReader br = IoUtils.getTextReader(inputFileS);
            String temp = null;
            String query = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("# Query:")) {
                    query = temp.split(" +")[2];
                    temp = br.readLine();
                    temp = br.readLine();
                    if (temp.startsWith("# 0")) {
                        AlignmentInfo ai = new AlignmentInfo (query, "");
                        aList.add(ai);
                    }
                    else {
                        temp = br.readLine();
                        int cnt = 0;
                        while (!(temp = br.readLine()).startsWith("#")) {
                            AlignmentInfo aInfo = this.getAlignmentInfoFromBlast(temp, eThresh);
                            if (aInfo == null) {
                                if (cnt == 0) aList.add(new AlignmentInfo(query, ""));
                            }
                            else {
                                aList.add(aInfo);
                            }
                            cnt++;
                        }
                    }
                }
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        ai = aList.toArray(new AlignmentInfo[aList.size()]);
        this.sortByQuery();
        System.out.println("ShortreadAlignment has "+String.valueOf(ai.length) + " alignments");
    }
    
    public void readFromBlast (String infileS, double eThresh) {
        ArrayList<AlignmentInfo> aList = new ArrayList();
        String temp = null;
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            System.out.println(br.readLine());
            String query = "";
            String hit = "";
            int qStartPos = Integer.MIN_VALUE;
            int qEndPos = Integer.MIN_VALUE;
            int hStartPos = Integer.MIN_VALUE;
            int hEndPos = Integer.MIN_VALUE;
            byte strand = Byte.MIN_VALUE;
            int gaps = Integer.MIN_VALUE;
            double evalue = 1;
            int score = Integer.MIN_VALUE;
            int matchNumber = Integer.MIN_VALUE;
            int editDistance = Integer.MIN_VALUE;
            //String query, String hit, int qStartPos, int qEndPos, int hStartPos, int hEndPos, byte strand, int gaps, double evalue, int score, int matchNumber, int editDistance
            while ((temp = br.readLine())!=null) {
                if (temp.startsWith("Query=")) {
                    query = temp.split("= ")[1];
                    while (!(temp = br.readLine()).startsWith(">")) continue;
                    hit = temp.split(" ")[1];
                }
                else if (temp.startsWith(" Score =")) {
                    String[] tem = temp.split("\\s+");
                    double v = Double.valueOf(tem[3]);
                    score = (int)v;
                    evalue = Double.valueOf(tem[8]);
                    tem = br.readLine().split("\\s+");
                    gaps = Integer.valueOf(tem[7].split("\\/")[0]);
                    matchNumber = Integer.valueOf(tem[3].split("\\/")[1])-gaps;
                    editDistance = matchNumber - Integer.valueOf(tem[3].split("\\/")[0]);
                    if (br.readLine().endsWith("Plus")) strand = 1;
                    else strand = 0;
                    br.readLine();
                    qStartPos = Integer.MIN_VALUE;
                    hStartPos = Integer.MIN_VALUE;
                    TIntArrayList posList = new TIntArrayList();
                    while ((temp = br.readLine()).startsWith("Query ")) {
                        tem = temp.split("\\ +");
                        int startIndex=5;
                        for (int j = startIndex; j <temp.length(); j++) {
                            char c = temp.charAt(j);
                            if (c > 57 || c==45) {
                                startIndex = j;
                                break;
                            }
                        }
                        if (qStartPos == Integer.MIN_VALUE) qStartPos = Integer.valueOf(tem[1]);
                        int currentPosition = Integer.valueOf(tem[1]);
                        qEndPos = Integer.valueOf(tem[3]);
                        String seq1 = tem[2];
                        temp = br.readLine();
                        String match = temp.substring(startIndex, temp.length());
                        temp = br.readLine();
                        tem = temp.split("\\ +");
                        if (hStartPos == Integer.MIN_VALUE) hStartPos = Integer.valueOf(tem[1]);
                        hEndPos = Integer.valueOf(tem[3]);
                        String seq2 = tem[2];
                        for (int i = 0; i < seq1.length(); i++) {
                            if (seq1.charAt(i) == '-') {
                                
                            }
                            else {
                                if (match.charAt(i) == ' ' && seq2.charAt(i) != '-') {
                                    posList.add(currentPosition);
                                }
                                currentPosition++;
                            }
                        }
                        temp = br.readLine();
                    }
                    int[] pos = posList.toArray();
                    AlignmentInfo aInfo = new AlignmentInfo(query, hit, qStartPos, qEndPos, hStartPos, hEndPos, strand, gaps, evalue, score, matchNumber, editDistance);
                    aInfo.addMismatchPosition(pos);
                    aList.add(aInfo);
                } 
            }
        }
        catch (Exception e) {
            System.out.println(temp);
            e.printStackTrace();
        }
        ai = aList.toArray(new AlignmentInfo[aList.size()]);
        this.sortByQuery();
        System.out.println("ShortreadAlignment has "+String.valueOf(ai.length) + " alignments");
    }
    /**
     * Sort alignments by query name
     */
    public void sortByQuery () {
        System.out.println("Start sorting by query");
        sortType = 0;
        Arrays.sort(ai);
        System.out.println("Finished sort");
    }
    
    /**
     * Sort alignments by hit and position on hit
     */
    public void sortByHitAndPos () {
        System.out.println("Start sorting by hit and pos");
        sortType = 1;
        Arrays.sort(ai);
        System.out.println("Finished sort");
    }
    
    AlignmentInfo getAlignmentInfoFromBlast (String inputStr, double eTresh) {
        String[] temp =inputStr.split("\\s+");
        String hit = "";
        byte strand = Byte.MIN_VALUE;
        int qStartPos = Integer.MIN_VALUE;
        int qEndPos = Integer.MIN_VALUE;
        int hStartPos = Integer.MIN_VALUE;
        int hEndPos = Integer.MIN_VALUE;
        int gaps = Integer.MIN_VALUE;
        double evalue = 1;
        int score = Integer.MIN_VALUE;
        int matchNumber = Integer.MIN_VALUE;
        int editDistance = Integer.MIN_VALUE;
        if (Double.valueOf(temp[10]) > eTresh) {
            return null;
        }
        else {
            hit = temp[1];
            gaps = Integer.valueOf(temp[5]);
            qStartPos = Integer.valueOf(temp[6]);
            qEndPos = Integer.valueOf(temp[7]);
            hStartPos = Integer.valueOf(temp[8]);
            hEndPos = Integer.valueOf(temp[9]);
            if (hStartPos <= hEndPos) {
                strand = 1;
            }
            else {
                strand = -1;
                int tempInt = hStartPos;
                hStartPos = hEndPos;
                hEndPos = tempInt;
            }
            matchNumber = (Integer.valueOf(temp[3]) - Integer.valueOf(temp[5]));
            editDistance = Short.valueOf(temp[4]);
            double d = Double.valueOf(temp[11]);
            score = (int)d;
            evalue = Double.valueOf(temp[10]);
            AlignmentInfo aInfo = new AlignmentInfo(temp[0], hit, qStartPos, qEndPos, hStartPos, hEndPos, strand, gaps, evalue, score, matchNumber, editDistance);
            return aInfo;
        }
    }
    
   /**
     * Return name of query
     * @param index
     * @return 
     */
    public String getQuery (int index) {
        return ai[index].query;
    }
    
    /**
     * Return name of hit, return empty String if the query does not match 
     * @param index
     * @return 
     */
    public String getHit (int index) {
        return ai[index].hit;
    }
    
    /**
     * Return start position of the alignment, return Integer.MIN_VALUE if the query does not match
     * @param index
     * @return 
     */
    public int getQueryStartPos (int index) {
        return ai[index].qStartPos;
    }
    
    /**
     * Return end position of the alignment, return Integer.MIN_VALUE if the query does not match
     * @param index
     * @return 
     */
    public int getQueryEndPos (int index) {
        return ai[index].qEndPos;
    }
    
     /**
     * Return start position of the alignment, return Integer.MIN_VALUE if the query does not match
     * @param index
     * @return 
     */
    public int getHitStartPos (int index) {
        return ai[index].hStartPos;
    }
    
    /**
     * Return end position of the alignment, return Integer.MIN_VALUE if the query does not match
     * @param index
     * @return 
     */
    public int getHitEndPos (int index) {
        return ai[index].hEndPos;
    }
    
    public byte getStrand (int index) {
        return ai[index].strand;
    }
    
    public int getGapNumber (int index) {
        return ai[index].gaps;
    }
    
    public double getEvalue (int index) {
        return ai[index].evalue;
    }
    
    public int getScore (int index) {
        return ai[index].score;
    }
    
    public int getMatchNumber (int index) {
        return ai[index].matchNumber;
    }
    
    public int getEditDistance (int index) {
        return ai[index].editDistance;
    }
    
    /**
     * Return number of alignment
     * @return 
     */
    public int getAlignmentNumber () {
        return ai.length;
    }
    
    /**
     * Return non-redundant names of queries
     * @return 
     */
    public String[] getQuerys () {
        TreeSet<String> querySet = new TreeSet();
        for (int i = 0; i < this.getAlignmentNumber(); i++) {
            querySet.add(ai[i].query);
        }
        return querySet.toArray(new String[querySet.size()]);
    }
    
    /**
     * Return non-redundant names of hits
     * @return 
     */
    public String[] getHits () {
        TreeSet<String> hitSet = new TreeSet();
        for (int i = 0; i < this.getAlignmentNumber(); i++) {
            if (ai[i].hit.equals("")) continue;
            hitSet.add(ai[i].hit);
        }
        return hitSet.toArray(new String[hitSet.size()]);
    }
    
    /**
     * Return if a query does not match the database
     * @param query
     * @return 
     */
    public boolean isNoMatch (String query) {
        int index = this.getAlignmentStartIndexByQuery(query);
        if (this.getHit(index).equals("")) return true;
        return false;
    }
    
    /**
     * Note: Need to be sorted by query first, otherwise program quits
     * @param query
     * @return 
     */
    public int getAlignmentNumberByQuery (String query) {
        if (this.sortType != 0) {
            System.out.println("Alignment should be sorted by query first. Program quits");
            System.exit(0);
        }
        return getAlignmentEndIndexByQuery(query) - getAlignmentStartIndexByQuery(query);
    }
    
     /**
     * Note: Need to be sorted by query first, otherwise program quites
     * @param query
     * @return 
     */
    public int getAlignmentStartIndexByQuery (String query) {
        if (this.sortType != 0) {
            System.out.println("Alignment should be sorted by query first. Program quits");
            System.exit(0);
        }
        int index = Arrays.binarySearch(ai, new AlignmentInfo(query, null));
        while (index > 0 && ai[index-1].query.equals(query)) {
            index--;
        }
        return index;
    }
    
    /**
     * Return the index of the last alignment of a query, exclusive
     * Note: Need to be sorted by query first
     * @param query
     * @return 
     */
    public int getAlignmentEndIndexByQuery (String query) {
        if (this.sortType != 0) {
            System.out.println("Alignment should be sorted by query first. Program quits");
            System.exit(0);
        }
        int index = Arrays.binarySearch(ai, new AlignmentInfo(query, null));
        if (index < 0) return index;
        while ((index+1) < this.getAlignmentNumber() && ai[index+1].query.equals(query)) {
            index++;
        }
        return index+1;
    }
    
     /**
     * Return alignment index
     * Note: Need to be sorted by hit and pos first, otherwise program quits
     * @param hit
     * @param pos
     * @return 
     */
    public int getAlignmentIndexByHitPos (String hit, int pos) {
        if (this.sortType != 1) {
            System.out.println("Alignment should be sorted by hit and pos first. Program quits");
            System.exit(0);
        }
        return Arrays.binarySearch(ai, new AlignmentInfo(hit, pos));
    }
    
    /**
     * Note: Need to be sorted by hit and pos first, otherwise program quit
     * @param hit
     * @return 
     */
    public int getAlignmentNumberByHit (String hit) {
        if (this.sortType != 1) {
            System.out.println("Alignment should be sorted by hit and pos first. Program quits");
            System.exit(0);
        }
        return getAlignmentEndIndexByHit(hit) - getAlignmentStartIndexByHit(hit);
    }
    
    /**
     * Return -1 if the hit doesn't exist
     * Note: Need to be sorted by hit and pos first, otherwise program quits
     * @param hit
     * @return 
     */
    public int getAlignmentStartIndexByHit (String hit) {
        if (this.sortType != 1) {
            System.out.println("Alignment should be sorted by hit and pos first. Program quits");
            System.exit(0);
        }
        int index = this.getAlignmentIndexByHitPos(hit, 0);
        if (index < 0) {
            index = -index -1;
            if (ai[index].hit.equals(hit)) return index;
            else return -1;
        }
        else {
            while (index > 0 && ai[index-1].hit.equals(hit)) {
                index--;
            }
            return index;
        }
    }
    
    /**
     * Return the index of the last alignment of a hit, exclusive.
     * Return -1 if the hit doesn't exist
     * Note: Need to be sorted by hit first
     * @param hit
     * @return 
     */
    public int getAlignmentEndIndexByHit (String hit) {
        if (this.sortType != 1) {
            System.out.println("Alignment should be sorted by hit and pos first. Program quits");
            System.exit(0);
        }
        int index = this.getAlignmentIndexByHitPos(hit, Integer.MAX_VALUE);
        if (index < 0) {
            index = - index - 2;
            if (ai[index].hit.equals(hit)) return index+1;
            else return -1;
        }
        else {
            while ((index+1) < this.getAlignmentNumber() && ai[index+1].hit.equals(hit)) {
                index++;
            }
            return index+1;
        }
    }
    
    /**
     * Return the ratio of edit distance vs length of matched range
     * @param index
     * @return NaN when there is no match
     */
    public double getEditDistanceRatio (int index) {
        if (this.getMatchNumber(index) > -1 && this.getEditDistance(index) > -1) return (double)this.getEditDistance(index)/(double)this.getMatchNumber(index);
        return Double.NaN;
    }
    
    class AlignmentInfo implements Comparable <AlignmentInfo>{
        /**Query sequence*/
        String query = null;
        /**Chromosome*/
        String hit = "";
        /**Query starting position of query sequence, inclusive*/
        int qStartPos = Integer.MIN_VALUE;
        /**Query ending position of query sequence, inclusive*/
        int qEndPos = Integer.MIN_VALUE;
        /**Reference starting position of query sequence, inclusive*/
        int hStartPos = Integer.MIN_VALUE;
        /**Reference ending position of query sequence, inclusive*/
        int hEndPos = Integer.MIN_VALUE;
        /**Strand of alignment, + = 1, - = -1*/
        byte strand = Byte.MIN_VALUE;
        /**gap open in alignment*/
        int gaps = Integer.MIN_VALUE;
        /**E value of alignment*/
        double evalue = 1;
        /**Bit score*/
        int score = Integer.MIN_VALUE;
         /**The length of matched range*/
        int matchNumber = Integer.MIN_VALUE;
        /**Mismatch number in the matched range*/
        int editDistance = Integer.MIN_VALUE;
        
        int[] qMismatchPos = null;
        
        
        AlignmentInfo (String query, String hit) {
            this.query = query;
            this.hit = hit;
        }
        
        AlignmentInfo (String hit, int pos) {
            this.hit = hit;
            this.hStartPos = pos;
        }
        
        AlignmentInfo (String query, String hit, int qStartPos, int qEndPos, int hStartPos, int hEndPos, byte strand, int gaps, double evalue, int score, int matchNumber, int editDistance) {
            this.query = query;
            this.hit = hit;
            this.qStartPos = qStartPos;
            this.qEndPos = qEndPos;
            this.hStartPos = hStartPos;
            this.hEndPos = hEndPos;
            this.strand = strand;
            this.gaps = gaps;
            this.evalue = evalue;
            this.score = score;
            this.matchNumber = matchNumber;
            this.editDistance = editDistance;
        }
        
        public void addMismatchPosition (int[] qMismatchPos) {
            this.qMismatchPos = qMismatchPos;
        }
        
        @Override
        public int compareTo(AlignmentInfo o) {
            if (sortType == 0) {
                return query.compareTo(o.query);
            }
            else if (sortType == 1) {
                if (hit.equals(o.hit)) {
                    return hStartPos - o.hStartPos;
                }
                else {
                    return hit.compareTo(o.hit);
                }
            }
            else return 0;
        }
    }
}
