/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.collaboration;

import analysis.panGenome.denovoEvaluation.DenovoEvaluation;
import format.Fasta;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class B104 {
    
    public B104 () {
        //this.genomeStatistics();
        this.test();
    }
    
    public void test () {
        
    }
    
    public void genomeStatistics () {
        String infileS = "M:\\collaboration\\nancy\\B104.fa";
        String outfileS = "M:\\collaboration\\nancy\\contigs.txt";
        Fasta f = new Fasta (infileS);
        f.sortRecordByLengthDescending();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            for (int i = 0; i < f.getSeqNumber(); i++) {
                bw.write(f.getName(i)+"\t"+String.valueOf(f.getSeqLength(i)));
                bw.newLine();
            }
            System.out.println("L50: "+String.valueOf(f.getL50()));
            System.out.println("N50: "+String.valueOf(f.getN50()));
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
