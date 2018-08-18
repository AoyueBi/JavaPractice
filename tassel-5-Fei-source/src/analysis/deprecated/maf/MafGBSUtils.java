/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.maf;

import format.Fasta;
import format.Table;
import java.io.BufferedWriter;
import java.io.File;
import java.util.TreeSet;
import net.maizegenetics.analysis.gbs.v2.DiscoverySNPCallerPluginV2;
import net.maizegenetics.analysis.gbs.v2.GBSSeqToTagDBPlugin;
import net.maizegenetics.analysis.gbs.v2.SAMToGBSdbPlugin;
import net.maizegenetics.analysis.gbs.v2.SNPQualityProfilerPlugin;
import net.maizegenetics.analysis.gbs.v2.TagExportToFastqPlugin;
import net.maizegenetics.dna.map.Chromosome;
import utils.IoUtils;

/**
 *
 * @author Fei Lu
 */
public class MafGBSUtils {
    
    public MafGBSUtils () {
        
    }
    
    public void SNPQuality () {
        String tagDBFileS = "/workdir/feilu/tagDB/tag.db";
        new SNPQualityProfilerPlugin()
                .dBFile(tagDBFileS)
                .performFunction(null);
    }
    
    public void SNPDiscovery () {
        String tagDBFileS = "/workdir/feilu/tagDB/tag.db";
        new DiscoverySNPCallerPluginV2()
                .inputDB(tagDBFileS)
                .minMinorAlleleFreq(0.0)
                .minLocusCoverage(0.02)
                .includeGaps(false)
                .startChromosome(new Chromosome("0"))
                .endChromosome(new Chromosome("12"))
                .performFunction(null);
    }
    
    public void SamToTagDB () {
        String tagDBFileS = "/workdir/feilu/tagDB/tag.db";
        String samFileS = "/workdir/feilu/alignment/tag.sam";
        new SAMToGBSdbPlugin()
                .gBSDBFile(tagDBFileS)
                .sAMInputFile(samFileS)
                .minAlignLength(50)
                .minAlignProportion(0.0)
                .performFunction(null); 
    }
    
    public void tagExportToDB () {
        String tagDBFileS = "/workdir/feilu/tagDB/tag.db";
        String tagFastqFileS = "/workdir/feilu/alignment/tag.fq";
        new TagExportToFastqPlugin()
                .inputDB(tagDBFileS)
                .outputFile(tagFastqFileS)
                .minCount(1)
                .performFunction(null);
    }
    
    public void seqToTagDB () {
        String fastqDirS = "/workdir/feilu/Illumina/";
        String keyFileS = "/workdir/feilu/key/maf_key_sortByFlowLane.txt";
        String tagDBFileS = "/workdir/feilu/tagDB/tag.db";
        new GBSSeqToTagDBPlugin()
                .enzyme("ApeKI")
                .inputDirectory(fastqDirS)
                .outputDatabaseFile(tagDBFileS)
                .keyFile(keyFileS)
                .minimumQualityScore(30)
                .batchSize(8)
                .performFunction(null);
    }
    
    public void writeUniqueLane () {
        String desKeyFileS = "M:\\production\\maf\\gbs\\key\\maf_key.txt";
        String uniqueLaneFileS = "M:\\production\\maf\\gbs\\key\\uniqueLane.txt";
        Table t = new Table (desKeyFileS);
        TreeSet<String> laneSet = new TreeSet();
        for (int i = 0; i < t.getRowNumber(); i++) {
            laneSet.add(t.content[i][0]+"_"+t.content[i][1]);
        }
        String[] lane = laneSet.toArray(new String[laneSet.size()]);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(uniqueLaneFileS);
            bw.write("Unique lanes");
            for (int i = 0; i < lane.length; i++) {
                bw.write(lane[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }
    
    public void selectKeyFile () {
        String sourceKeyFileS = "M:\\production\\maf\\gbs\\key\\source\\AllZeaGBSv2.7_SeqToGenos_parts01to21_key.txt";
        String[] keyword = {"Ames", "Teo_Div", "TeoDNA", "SEED"};
        String desKeyFileS = "M:\\production\\maf\\gbs\\key\\maf_key.txt";
        Table t = new Table (sourceKeyFileS);
        boolean[] ifOut = new boolean[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            for (int j = 0; j < keyword.length; j++) {
                if (t.content[i][4].startsWith(keyword[j])) {
                    ifOut[i] = true;
                    break;
                }
            }
        }
        t.writeTable(desKeyFileS, ifOut);
    }
    
    public void generateAGPV3Reference () {
        String downloadDirS = "M:\\Database\\maizeReference\\agpV3\\download\\uncompressed\\";
        String outputFileS = "M:\\Database\\maizeReference\\agpV3\\maize.agpV3.fa";
        String[] chromMark = {"1","2","3","4","5","6","7","8","9","10","Mt","Pt","nonchromosomal"};
        String pre = "Zea_mays.AGPv3.25.dna.";
        String post = ".fa";
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outputFileS);
            for (int i = 0; i < chromMark.length-1; i++) {
                String inputFileS = new File(downloadDirS,pre+"chromosome."+chromMark[i]+post).getAbsolutePath();
                Fasta f = new Fasta(inputFileS);
                bw.write(">"+String.valueOf(i+1));
                bw.newLine();
                bw.write(f.getSeq(0));
                bw.newLine();
            }
            String inputFileS = new File(downloadDirS,pre+chromMark[chromMark.length-1]+post).getAbsolutePath();
            Fasta f = new Fasta(inputFileS);
            bw.write(">0");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            String nString = "";
            for (int i = 0; i < 50; i ++) {
                nString+="N";
            }
            for (int i = 0; i < f.getSeqNumber()-1; i++) {
                sb.append(f.getSeq(i)).append(nString);
            }
            sb.append(f.getSeq(f.getSeqNumber()-1));
            bw.write(sb.toString());
            bw.newLine();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
}
