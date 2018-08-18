/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package analysis.collaboration;

import format.Fasta;
import format.Sequence;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import utils.IoUtils;

/**
 *
 * @author Fei Lu <fl262@cornell.edu>
 */
public class VTE {
    
    public VTE () {
        //this.findOrthologs();
        this.findCML247MinorOrthologs();
        //this.findB73Ortholog();
    }
    
    public void findB73Ortholog () {
        String input = "N:\\Zea\\reference_genome\\ZmB73_RefGen_v2.fa";
        String output = "M:\\collaboration\\vte\\contigs\\B73_contig.fa";
        String chr = "9";
        int startPos = 107400000;
        int endPos = 107700000;
        Fasta f = new Fasta (input);
        f.sortRecordByName();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(output);
            bw.write(">"+chr+"-"+String.valueOf(startPos)+"-"+String.valueOf(endPos));
            bw.newLine();
            bw.write(f.getSeq(f.getIndex(chr)).substring(startPos-1, endPos-1));
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void findCML247MinorOrthologs () {
        String input = "M:\\production\\panGenome\\cml247\\nrgene\\source\\CML247_maize.fasta";
        String contigFileS = "M:\\collaboration\\vte\\contigs\\cml247_contig_minor.fa";
        Fasta f = new Fasta (input);
        f.sortRecordByName();
        String name = null;
        String seq = null;
        try {
            name = "scaffold3676941-add";
            seq = f.getSeq(f.getIndex(name));
            BufferedWriter bw = IoUtils.getTextWriter(contigFileS);
            bw.write(">"+name);
            bw.newLine();
            bw.write(seq);
            bw.newLine();
            name = "scaffold743254-add";
            seq = f.getSeq(f.getIndex(name));
            bw.write(">"+name);
            bw.newLine();
            bw.write(seq);
            bw.newLine();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void findOrthologs () {
        //String input = "M:\\production\\panGenome\\mo17\\genome\\Mo17_scaffolds_201403.fa";
        //String alignmentFileS = "M:\\collaboration\\vte\\alignment\\mo17_alignment_7.txt";
        //String contigFileS = "M:\\collaboration\\vte\\contigs\\mo17_contig.fa";
        
        //String input = "M:\\production\\panGenome\\cml247\\nrgene\\source\\CML247_maize.fasta";
        //String alignmentFileS = "M:\\collaboration\\vte\\alignment\\cml247_alignment_7.txt";
        //String contigFileS = "M:\\collaboration\\vte\\contigs\\cml247_contig.fa";
        
        String input = "M:\\production\\panGenome\\ph207\\scaffold\\PH207.ordered_scaffolds.fa";
        String alignmentFileS = "M:\\collaboration\\vte\\alignment\\PH207_alignment_7.txt";
        String contigFileS = "M:\\collaboration\\vte\\contigs\\PH207_contig.fa";
        
        Fasta f = new Fasta (input);
        f.sortRecordByName();
        String name = null;
        String seq = null;
        try {
            BufferedReader br = IoUtils.getTextReader(alignmentFileS);
            for (int i = 0; i < 5; i++) br.readLine();
            String[] temp = br.readLine().split("\t");
            name = temp[1];
            seq = f.getSeq(f.getIndex(name));
            if (Integer.valueOf(temp[8]) < Integer.valueOf(temp[9])) {
                seq = new Sequence(seq).getReverseComplementarySeq();
                name = name+"-r";
            }
            BufferedWriter bw = IoUtils.getTextWriter(contigFileS);
            bw.write(">"+name);
            bw.newLine();
            bw.write(seq);
            bw.newLine();
            bw.flush();
            bw.close();
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
