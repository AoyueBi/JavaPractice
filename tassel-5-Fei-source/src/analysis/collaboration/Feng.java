/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.collaboration;

import format.Table;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;

/**
 *
 * @author fl262
 */
class Feng {
    
    public Feng () {
        this.getData();
        
    }
    
    void getData () {
        //String infileS = "M:\\collaboration\\feng\\ZeaGBSv27_20150803_AGPv3_subBy282UniqueHighCov.h5";
        String infileS = "M:\\collaboration\\feng\\ZeaGBSv27impV5_20160627.h5";
        String amesFileS = "E:\\Research\\old\\pav\\gwas\\phenotype\\source\\Ames_282_NAM_NameModified.txt";
        String outfileS = "M:\\collaboration\\feng\\nam.cnnam.ames.GBS27.hmp.txt";
        Table t = new Table (amesFileS);
        ArrayList<String> amesList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.content[i][0].equals("AMES")) {
                amesList.add(t.content[i][2]);
            }
        }
        String[] ames = amesList.toArray(new String[amesList.size()]);
        Arrays.sort(ames);
        GenotypeTable gt = ImportUtils.read(infileS);
        TaxaList tList = gt.taxa();
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < tList.size(); i++) {
            String name = tList.get(i).getName();
            String[] temp = name.split(":");
            //System.out.println(name);
            int index = Arrays.binarySearch(ames, temp[0]);
            if (name.startsWith("CN-") || name.startsWith("Z0") || index >= 0) {
                tlb.add(tList.get(i));  
            }
        }
        tList = tlb.build();
        GenotypeTable sgt = FilterGenotypeTable.getInstance(gt, tList);
        ExportUtils.writeToHapmap(sgt, outfileS);
    }
    
}
