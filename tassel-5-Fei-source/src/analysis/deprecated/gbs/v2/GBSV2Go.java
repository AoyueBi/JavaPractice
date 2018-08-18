/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.gbs.v2;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;
import net.maizegenetics.analysis.gbs.v2.DiscoverySNPCallerPluginV2;
import net.maizegenetics.analysis.gbs.v2.GBSSeqToTagDBPlugin;
import net.maizegenetics.analysis.gbs.v2.SAMToGBSdbPlugin;
import net.maizegenetics.analysis.gbs.v2.SNPQualityProfilerPlugin;
import net.maizegenetics.analysis.gbs.v2.TagExportToFastqPlugin;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.tag.TagCounts;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;

/**
 *
 * @author Fei Lu
 */
public class GBSV2Go {
    
    public GBSV2Go () {
        //this.seqToTagDB();
        //this.TagExportToFastq();
        //this.SamToTagDB();
        //this.SNPDiscovery();
        this.SNPQuality();
    }
    
    public void SNPQuality () {
        String tagDBFileS = "M:\\pipelineTest\\gbsV2\\tagDB\\v2Tag.db";
        new SNPQualityProfilerPlugin()
                .dBFile(tagDBFileS)
                .performFunction(null);
    }
    
    public void SNPDiscovery () {
        String tagDBFileS = "M:\\pipelineTest\\gbsV2\\tagDB\\v2Tag.db";
        new DiscoverySNPCallerPluginV2()
                .inputDB(tagDBFileS)
                .minMinorAlleleFreq(0.00001)
                .minLocusCoverage(0.00001)
                .includeGaps(false)
                .startChromosome(new Chromosome("1"))
                .endChromosome(new Chromosome("12"))
                .performFunction(null);
    }
    
    public void SamToTagDB () {
        String tagDBFileS = "M:\\pipelineTest\\gbsV2\\tagDB\\v2Tag.db";
        String samFileS = "M:\\pipelineTest\\gbsV2\\alignment\\tag.sam";
        new SAMToGBSdbPlugin()
                .gBSDBFile(tagDBFileS)
                .sAMInputFile(samFileS)
                .minAlignLength(50)
                .minAlignProportion(0.0)
                .performFunction(null); 
    }
    
    public void TagExportToFastq () {
        String tagDBFileS = "M:\\pipelineTest\\gbsV2\\tagDB\\v2Tag.db";
        String tagFastqFileS = "M:\\pipelineTest\\gbsV2\\tagFastq\\tag.fq";
        new TagExportToFastqPlugin()
                .inputDB(tagDBFileS)
                .outputFile(tagFastqFileS)
                .minCount(1)
                .performFunction(null);
    }
    
    public void seqToTagDB () {
        String fastqDirS = "M:\\pipelineTest\\gbsV2\\Illumina\\smaller\\";
        String keyFileS = "M:\\pipelineTest\\gbsV2\\key\\81FE7ABXX_key.txt";
        
        //String fastqDirS = "M:\\pipelineTest\\gbsV2\\Illumina\\normal\\";
        //String keyFileS = "M:\\pipelineTest\\gbsV2\\key\\normal_key.txt";
        
        String tagDBFileS = "M:\\pipelineTest\\gbsV2\\tagDB\\v2Tag.db";
        new GBSSeqToTagDBPlugin()
                .enzyme("ApeKI")
                .inputDirectory(fastqDirS)
                .outputDatabaseFile(tagDBFileS)
                .keyFile(keyFileS)
                .minimumQualityScore(20)
                .batchSize(6)
                .performFunction(null);
    }
    
    
    public static void main (String[] args) {
        new GBSV2Go ();
    }
    
}
