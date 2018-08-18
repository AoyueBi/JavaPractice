/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.maf;

import graphcis.r.DensityPlot;
import graphcis.r.DensityPlotMultiClass;
import graphcis.r.LineChart;
import graphcis.r.ScatterPlot;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;
import net.maizegenetics.analysis.distance.IBSDistanceMatrix;
import net.maizegenetics.dna.snp.GenotypeTable;

/**
 *
 * @author Fei Lu
 */
public class MafGo {
    
    public MafGo () {
        //this.GBSV2Pipe();
        this.WGSPaperPipe();
    }
    
    public void WGSPaperPipe () {
       new UniqueKmers();
       new WGSDepth();
       //new InvariantSiteValidation ();
       //new Footprint();
       //new RangeAnnotation();
        
       //new InvariantSiteDistribution();
       //new InvariantSiteAnnotation();
       
       
    }
    
    public void GBSV2Pipe () {
        //this.generateAGPV3Reference();
        //this.selectKeyFile();
        //this.seqToTagDB();
        //this.tagExportToFastq();
        //this.SamToTagDB();
        //this.SNPDiscovery();
        //this.SNPQuality();
    }
    
    public void SNPQuality () {
        new MafGBSUtils().SNPQuality();
    }
    
    public void SNPDiscovery () {
        new MafGBSUtils().SNPDiscovery();
    }
    
    public void SamToTagDB () {
        new MafGBSUtils().SamToTagDB();
    }
    public void runBowtie2 () {
        //"bowtie2 -M 4 -p 60 --very-sensitive-local -x ../AGPv3/ZMAGP3_9103 -U tagsForAlign.fa -S tagsForAlignFullvs.sam";
    }
    
    public void tagExportToFastq () {
        new MafGBSUtils().tagExportToDB();
    }
    
    public void seqToTagDB () {
        new MafGBSUtils().seqToTagDB();
    }
    
    public void selectKeyFile () {
        new MafGBSUtils().selectKeyFile();
        new MafGBSUtils().writeUniqueLane();
    }
    
    public void generateAGPV3Reference () {
        new MafGBSUtils().generateAGPV3Reference();
    }
    
    public static void main (String[] args) {
        new MafGo();
    }
}
