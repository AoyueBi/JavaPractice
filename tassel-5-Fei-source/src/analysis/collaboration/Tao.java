/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.collaboration;

import format.Fasta;

/**
 *
 * @author fl262
 */
class Tao {
    
    Tao () {
        this.reformat();
    }
    
    void reformat () {
        String infileS = "M:\\collaboration\\tao\\Sbicolor_313_v3.0.fa";
        String outfileS = "M:\\collaboration\\tao\\Sorghum_v3.0.fa";
        int seqNum = 11;
        Fasta f = new Fasta (infileS);
        String[] names = new String[seqNum];
        String[] seqs = new String[seqNum];
        int[] ids = new int[seqNum];
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < f.getSeqNumber(); i++) {
            String seqName = f.getName(i);
            if (seqName.startsWith("Chr")) {
                names[i] = String.valueOf(Integer.valueOf(seqName.replaceFirst("Chr", "")));
                seqs[i] = f.getSeq(i);
                ids[i] = i+1;
            }
            else {
                sb.append(f.getSeq(i)).append("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
            }
        }
        sb.delete(sb.length()-50, sb.length());
        names[10] = "11";
        seqs[10] = sb.toString();
        ids[10] = 11;
        Fasta ff = new Fasta (names, seqs, ids);
        ff.writeFasta(outfileS);
    }
}
