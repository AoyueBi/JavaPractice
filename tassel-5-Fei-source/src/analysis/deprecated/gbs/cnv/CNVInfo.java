/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.gbs.cnv;

/**
 *
 * @author Fei Lu
 */
public class CNVInfo implements Comparable<CNVInfo> {
    String taxon;
    int chr;
    int binStart;
    int length;
    /**Deletion true, duplication false*/
    boolean ifDeletion;
    
    public CNVInfo (String taxon, int chr, int binStart, int length, boolean ifDeletion) {
        this.taxon = taxon;
        this.chr = chr;
        this.binStart = binStart;
        this.length = length;
        this.ifDeletion = ifDeletion;
    }

    @Override
    public int compareTo(CNVInfo o) {
        int v = taxon.compareTo(o.taxon);
        if (v == 0) {
            if (this.chr == o.chr) {
                return this.binStart-o.binStart;
            }
            else if (this.chr > o.chr){
                return 1;
            }
            else {
                return -1;
            }
        }
        else if (v > 0) {
            return 1;
        }
        else {
            return -1;
        }
    }
}
