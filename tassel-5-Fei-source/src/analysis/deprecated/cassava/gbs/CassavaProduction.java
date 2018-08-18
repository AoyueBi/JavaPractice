/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package analysis.deprecated.cassava.gbs;

/**
 *
 * @author Fei Lu 
 */
public class CassavaProduction {
    
    public CassavaProduction () {
        this.filterData();
        //this.splitHapMapToChromosomes();
        this.splitH5Tochromosomes();
    }
    
    public void splitH5Tochromosomes () {
        String inputGenotypeFileS = "M:\\production\\cassava\\gbs\\filtered_genotype\\genotypeFilterByHets.h5";
        String outputDirS = "M:\\production\\cassava\\gbs\\filtered_genotype\\chromosomes\\";
        new cassavaUtils().splictH5ToChromosomes(inputGenotypeFileS, outputDirS);
    }
    
    public void splitHapMapToChromosomes () {
        String inputGenotypeFileS = "M:\\production\\cassava\\gbs\\filtered_genotype\\genotypeFilterByHets.hmp.txt";
        String outputDirS = "M:\\production\\cassava\\gbs\\filtered_genotype\\chromosomes\\";
        new cassavaUtils().splitHapMapToChromosomes(inputGenotypeFileS, outputDirS, 19);
    }
    
    public void filterData () {
        String inputGenotypeFileS = "M:\\production\\cassava\\gbs\\original_genotype\\cassava_mergeProdHapMap_noDepths_20140421.h5";
        String outputGenotypeFileS = "M:\\production\\cassava\\gbs\\filtered_genotype\\genotypeFilterByHets.h5";
        new cassavaUtils().filterSite(inputGenotypeFileS, outputGenotypeFileS, 0.5, 0.1, 0.1, 0.05);
    }
    
    public static void main (String[] args) {
        new CassavaProduction();
    }
}
