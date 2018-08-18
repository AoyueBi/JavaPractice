/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package analysis.deprecated.cassava.gbs;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;

/**
 *
 * @author Fei Lu
 */
public class CassavaTest {
    
    public CassavaTest () {
        this.filterSite();
        //this.filterSiteEndChr1();
        //this.mkBinFile();
    }
    
    public void filterSiteEndChr1 () {
        new cassavaUtils().filterSiteEndChr1();
        
    }
    
    public void mkBinFile () {
        String statisticFileS = "M:\\pipelineTest\\cassava\\gbs\\siteSummary\\chr1SiteSummary.txt";
        String binFileS = "M:\\pipelineTest\\cassava\\gbs\\siteSummary\\chr1BinSummary.txt";
        new cassavaUtils().mkBinSummary(statisticFileS, binFileS);
    }
    
    public void getChr1Statistics () {
        String inputGenotypeFileS =  "M:\\pipelineTest\\cassava\\gbs\\filtered_genotype\\chr1_filter.h5";
        String statisticFileS = "M:\\pipelineTest\\cassava\\gbs\\siteSummary\\chr1SiteSummary.txt";
        
    }
    
    public void filterSite () {
        //String inputGenotypeFileS = "M:\\pipelineTest\\cassava\\gbs\\original_genotype\\chr1.hmp.txt";
        String inputGenotypeFileS = "M:\\pipelineTest\\cassava\\gbs\\original_genotype\\igdSamples_ch1.h5";
        String outputGenotypeFileS = "M:\\pipelineTest\\cassava\\gbs\\filtered_genotype\\chr1_filter.h5";
        new cassavaUtils().filterSite(inputGenotypeFileS, outputGenotypeFileS, 0.5, 0.1, 0.1, 0.05);
    }
    
    public void mkTaxaNameFile () {
        String inputGenotypeFileS = "M:\\pipelineTest\\cassava\\gbs\\original_genotype\\chr1.hmp.txt";
        String taxaFileS = "M:\\pipelineTest\\cassava\\gbs\\taxaSummary\\taxaFullName.txt";
        new cassavaUtils().mkTaxaNameFileS(inputGenotypeFileS, taxaFileS);
    }
    
    public static void main (String[] args) {
        new CassavaTest();
    }
}
