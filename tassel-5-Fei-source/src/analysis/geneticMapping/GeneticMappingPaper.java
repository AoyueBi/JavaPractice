/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package analysis.geneticMapping;

import format.Fasta;

/**
 *
 * @author Fei Lu 
 */
public class GeneticMappingPaper {
    
    public GeneticMappingPaper() {
        //this.genomeCoverage();
        //this.PEMapping();
        //this.identifyPAVTag();
        //this.pavValidation();
        //this.MultipleML();
        //this.accuracyFactor();
        //this.accuracyGeneticDistance();
        //this.populationStructure();
        //this.jumpLibArrcuracy();//don't work, ignore
        //this.manyOthers();
        
    }
    
    public void manyOthers () {
        new ManyOthers ();
    }
    
    public void jumpLibArrcuracy () {
        new JumpLibAccuracy();
    }
    
    public void populationStructure () {
        new PopulationStructure();
    }
    
    public void accuracyGeneticDistance () {
        new AccuracyGeneticDistance();
    }
    
    public void accuracyFactor () {
        new AccuracyFactor();
    }
    
    public void MultipleML () {
        new MultipleML();
    }
    
    public void pavValidation () {
        new PAVValidation();
    }
    
    public void identifyPAVTag () {
        new IdentifyPAVTag();
    }
    
     public void PEMapping () {
        new PEMapping ();
    }
    
    public void genomeCoverage () {
        new GenomeCoverage();
    }
    
    public static void main (String[] args) {
        new GeneticMappingPaper();
    }
}
