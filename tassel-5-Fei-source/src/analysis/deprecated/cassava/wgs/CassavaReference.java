/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.cassava.wgs;

import format.Fasta;
import format.Table;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Set;

/**
 *
 * @author fl262
 */
public class CassavaReference {
    
    public CassavaReference () {
        //this.outputBaitScaffold();
        //this.reformatReference();
    }
    
    public void reformatReference () {
        String infileS = "O:\\Cassava\\GenomeReference\\cassavaV6With23Chromosomes\\cassavaV6_chr1-18And19And20_plastidsAndRepeats.fa";
        String outfileS = "O:\\Cassava\\GenomeReference\\cassavaV6With23Chromosomes\\cassavaV6_23Chr.fa";
        Fasta f = new Fasta(infileS);
        f.writeFasta(outfileS);
    }
    
    public void outputBaitScaffold () {
        String fastaFileS = "M:\\Database\\cassavaReference\\genome\\Manihot esculenta\\source\\MesV6.fa";
        String assemblyScaffoldList = "M:\\Database\\cassavaReference\\genome\\Manihot esculenta\\source\\productionAssembly\\scaffoldList.txt";
        String baitScaffold = "M:\\Database\\cassavaReference\\genome\\Manihot esculenta\\source\\baitScaffold.fa";
        Fasta f = new Fasta (fastaFileS);
        String[] names = f.getNames();
        Table t = new Table (assemblyScaffoldList);
        String[] assSca = new String[t.getRowNumber()] ;
        for (int i = 0; i < t.getRowNumber(); i++) {
            assSca[i] = t.content[i][0].toLowerCase();
        }
        Arrays.sort(assSca);
        boolean[] ifOut = new boolean[f.getSeqNumber()];
        int cnt = 0;
        int cnt2 = 0;
        for (int i = 0; i < f.getSeqNumber(); i++) {
            String head = f.getName(i).split(" ")[0].replaceFirst(">", "");
            if (head.contains("Chromosome")) ifOut[i] = false;
            else if (head.contains("chloroplast")) ifOut[i] = false;
            else {
                int index = Arrays.binarySearch(assSca, head.toLowerCase());
                if (index < 0) {
                    ifOut[i] = true;
                    cnt++;
                }
                else {
                    cnt2++;
                }
            }
        }
        f.writeFasta(baitScaffold, ifOut);
        System.out.println(cnt);
        System.out.println(cnt2);
    }
    
}
