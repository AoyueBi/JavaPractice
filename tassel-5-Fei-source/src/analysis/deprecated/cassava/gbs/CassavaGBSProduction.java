/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.cassava.gbs;

import net.maizegenetics.analysis.gbs.DiscoverySNPCallerPlugin;
import net.maizegenetics.analysis.gbs.MergeTagsByTaxaFilesPlugin;
import net.maizegenetics.analysis.gbs.SAMConverterPlugin;
import net.maizegenetics.analysis.gbs.TagCountToFastqPlugin;


/**
 *
 * @author fl262
 * @deprecated 
 */
public class CassavaGBSProduction {
    
    public CassavaGBSProduction () {
        //this.mergeTBT();
        //this.mkTagCount();
        //this.mkFastqFile();
        //this.mkTOPM();
        //this.callSNPOnServer();
        //this.mergeTaxaTBT(); //For Martha (very short term purpose)
    }
    
    public void findHaplotypes () {

        
        String hapMapFileS = "M:/cassava/gbs/hapmap/myGBSGenos.chr1.hmp.txt";
        String haplotypeFileS = "M:/cassava/gbs/haplotype/hap.txt";
        String errorFileS = "M:/cassava/gbs/haplotype/error.txt";
        
        
        String arguments = "-i " + hapMapFileS + " -o " + haplotypeFileS + " -oE " + errorFileS + " 1 1 ";
        
    }
    
    public void callSNPOnServer () {
        String topmFileS = "/workdir/mingh/master.topm";
        //String pedigreeFileS = "/workdir/mingh/AllZeaPedigree2012oct01C.txt";
        String tbtFileS = "/workdir/mingh/master.tbt.byte";
        String referenceFileS = "/workdir/mingh/cassava.fa";
        String hapmapFileS = "/workdir/mingh/hapmap/myGBSGenos.chr+.hmp.txt";
        
        String arguments = "-i " + tbtFileS + " -y -m " + topmFileS + " -o " + hapmapFileS + " -ref " + referenceFileS;
        arguments = arguments + " -mxSites 700000 -mnMAF 0.02 -mnMAC 10 -mnLCov 0.1 -sC 1 -eC 500";
        String[] args = arguments.split(" ");
        
        DiscoverySNPCallerPlugin d = new DiscoverySNPCallerPlugin ();
        d.setParameters(args);
		d.performFunction(null);
    }
    
    public void mkTOPM () {
        String samFileS = "M:/cassava/gbs/alignment/master.sam";
        String topmFileS = "M:/cassava/gbs/topm/master.topm";
        SAMConverterPlugin convert = new SAMConverterPlugin();
        String arguments = "-i " + samFileS + " -o " + topmFileS;
        String[] args = arguments.split(" ");
        convert.setParameters(args);
		convert.performFunction(null);
        
    }
    
    public void mkFastqFile () {
        String tagCountFileS = "M:/cassava/gbs/tagCount/master.cnt";
        String fastqFileS = "M:/cassava/gbs/alignment/master.fq";
        TagCountToFastqPlugin t2f = new TagCountToFastqPlugin();
        String arguments = "-i " + tagCountFileS + " -o " + fastqFileS;
        String[] args = arguments.split(" ");
        t2f.setParameters(args);
		t2f.performFunction(null);
    }
    
    public void mkTagCount () {
        String tbtFileS = "M:/cassava/gbs/tbt/master.tbt.byte";
        String tagCountFileS = "M:/cassava/gbs/tagCount/master.cnt";
        cassavaUtils util = new cassavaUtils();
        util.mkTagCountFromTBT(tbtFileS, tagCountFileS);
    }
    
    public void mergeTaxaTBT () {
        String tbtSourceFileS = "M:/cassava/gbs/tbt/master.tbt.byte";
        String tbtTargetFileS = "M:/cassava/gbs/tbt/master.mTaxa.tbt.byte";
        cassavaUtils util = new cassavaUtils();
        util.mergeTaxaTBT(tbtSourceFileS, tbtTargetFileS);
    }
    
    public void mergeTBT () {
        String tbtSourceDirS = "M:/cassava/gbs/source/tbt/";
        String tbtFileS = "M:/cassava/gbs/tbt/master.tbt.byte";
        MergeTagsByTaxaFilesPlugin mp = new MergeTagsByTaxaFilesPlugin();
        String arguments = "-i " + tbtSourceDirS + " -o " + tbtFileS;
        String[] args = arguments.split(" ");
        mp.setParameters(args);
	mp.performFunction(null);
    }
    
    public static void main (String[] args) {
        new CassavaGBSProduction();
    }
}
