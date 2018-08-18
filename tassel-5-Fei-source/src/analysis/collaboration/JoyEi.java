/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.collaboration;

import format.Fasta;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.HashSet;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class JoyEi {
    public JoyEi () {
        //this.getOrthologous1();
        this.getOrthologous2();
    }
    
    public void getOrthologous2 () {
        String infileS = "E:\\Research\\cml247\\fasta\\cml247.fas";
        String alignmentFileS = "M:\\collaboration\\joy-ei\\cds_out.txt";
        String outfileS = "M:\\collaboration\\joy-ei\\cds_seq.fa";
        Fasta f = new Fasta (infileS);
        f.sortRecordByName();
        try {
            HashSet<String> nSet = new HashSet();
            BufferedReader br = IoUtils.getTextReader(alignmentFileS);
            String temp;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) continue;
                nSet.add(temp.split("\t")[1]);
            }
            String[] name = nSet.toArray(new String[nSet.size()]);
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            for (int i = 0; i < name.length; i++) {
                bw.write(">"+name[i]);
                bw.newLine();
                int index = f.getIndex(name[i]);
                bw.write(f.getSeq(index));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void getOrthologous1 () {
        String infileS = "E:\\Research\\cml247\\fasta\\cml247.fas";
        String scaffold = "scaffold386926";
        Fasta f = new Fasta (infileS);
        f.sortRecordByName();
        int index = f.getIndex(scaffold);
        String outfileS = "M:\\collaboration\\joy-ei\\seq.fa";
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write(scaffold);
            bw.newLine();
            bw.write(f.getSeq(index));
            bw.newLine();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
}
