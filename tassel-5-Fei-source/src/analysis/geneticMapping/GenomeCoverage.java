/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package analysis.geneticMapping;

import format.ShortreadAlignment;
import java.util.HashSet;
import net.maizegenetics.dna.tag.TagCounts;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import net.maizegenetics.dna.tag.TagsByTaxaByte;
import utils.IOFileFormat;

/**
 *
 * @author Fei Lu <fl262@cornell.edu>
 */
public class GenomeCoverage {
    
    public GenomeCoverage () {
        //this.convertTCToFasta();
        //this.mkSimpleAlignment();
        //this.countUniquePosition();
    }
    
    public void countUniquePosition () {
        String input = "/SSD/mingh/merged.tc.sa";
        ShortreadAlignment sa = new ShortreadAlignment(input, IOFileFormat.Binary);
        sa.sortByQuery();

        HashSet<String> uniqueSet = new HashSet();
        for (int i = 0; i < sa.getAlignmentNumber(); i++) {
            if (!sa.isPerfectMatch(i)) continue;
            String site = sa.getHit(i)+"_"+String.valueOf(sa.getStartPos(i))+"_"+String.valueOf(sa.getStrand(i));
            uniqueSet.add(site);
        }
        System.out.println(uniqueSet.size()+" unique positions");
        
    }
    
    public void mkSimpleAlignment () {
        String input = "/SSD/mingh/result.sam";
        String output = "/SSD/mingh/merged.tc.sa";
        
        ShortreadAlignment sa = new ShortreadAlignment();
        sa.readFromBowtie2(input);
        sa.writeSimpleAlignment(output, IOFileFormat.Binary);
    }
    
    public void convertTCToFasta () {
        String input = "N:\\Zea\\build20120110\\mergedTagCounts\\merged.cnt";
        String output = "E:\\Research\\geneticMapping\\genomeCoverage\\merged.tc.fa";
        
        
        TagCounts tc = new TagCounts(input, FilePacking.Byte);
        System.out.println(tc.getTagCount());
        int min = 50;
        for (int i = 0; i < tc.getTagCount(); i++) {
            if (tc.getReadCount(i) < min) min = tc.getReadCount(i);
        }
        System.out.println(min);
        tc.toFASTA(output);
    }
}
