/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.elementMap;


import format.Fasta;
import gnu.trove.map.hash.TByteByteHashMap;

import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import net.maizegenetics.dna.BaseEncoder;
import net.openhft.koloboke.collect.set.hash.HashLongSet;
import net.openhft.koloboke.collect.set.hash.HashLongSets;

import utils.Benchmark;
import utils.FArrayUtils;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class ElementMapGo {
    
    public ElementMapGo () {
        //this.calculationPipe();
        //this.supportPipe();
        this.analysisPipe();
    }
    
    public void analysisPipe () {
        //new ScoreCharacter();
        //new UpstreamExpression();
    }
    
    public void calculationPipe () {
        //this.generateReferenceKmers();
        //this.calculateCpScoreInTaxon(); //most for test
        //this.countKmerInAllTaxa();
        //this.calibrateTaxonDepth();
        //this.combineScore();
    }
    
    private void combineScore () {
        String taxaDepthFileS = "/workdir/mingh/taxaDepth/taxaDepth.txt";
        String taxaScoreDirS = "/workdir/mingh/taxaScore/";
        String combinedScoreDirS = "/workdir/mingh/combinedScore/";
        String referenceGenomeFileS = "/workdir/mingh/referenceGenome/maize.agpV3.fa";
        double minDepth = 5;
        new CopyNumberScore().writeCombinedScore(taxaDepthFileS, taxaScoreDirS, referenceGenomeFileS, combinedScoreDirS, minDepth);
    }
    
    private void calibrateTaxonDepth() {
        String singleCopyFileS = "/workdir/mingh/taxaDepth/singleCopySite.txt";
        String taxaScoreDirS = "/workdir/mingh/taxaScore/";
        String taxaDepthDirS = "/workdir/mingh/taxaDepth/";
        new CopyNumberScore().calibrateTaxonDepth(singleCopyFileS, taxaScoreDirS, taxaDepthDirS);
    }
    
    private void countKmerInAllTaxa () {
//        String taxaFqMapFileS = "/workdir/mingh/taxaFastqMap/taxaFqMap.txt";
//        String fastqDirS = "/workdir/mingh/fastq/";
//        String kmerListDirS = "/workdir/mingh/referenceKmer/32_kmerList/";
//        String referenceGenomeFileS = "/workdir/mingh/referenceGenome/maize.agpV3.fa";
//        String taxaScoreDirS = "/workdir/mingh/taxaScore/";
        String taxaFqMapFileS = "/workdir/fl262/taxaFastqMap/taxaFqMap.txt";
        String fastqDirS = "/workdir/fl262/fastq/";
        String kmerListDirS = "/workdir/fl262/referenceKmer/32_kmerList/";
        String referenceGenomeFileS = "/workdir/fl262/referenceGenome/maize.agpV3.fa";
        String taxaScoreDirS = "/workdir/fl262/taxaScore/";
        new CopyNumberScore().writeCpScoreOfTaxa(taxaFqMapFileS, fastqDirS, kmerListDirS, referenceGenomeFileS, taxaScoreDirS, 9);
    }
    
    private void calculateCpScoreInTaxon () {
        String r1FileS = "M:\\pipelineTest\\elementMap\\fastq\\test_1.fq.gz";
        String r2FileS = "M:\\pipelineTest\\elementMap\\fastq\\test_2.fq.gz";
        String kmerListDirS = "M:\\pipelineTest\\elementMap\\referenceKmer\\kmerList\\";
        String referenceGenomeFileS = "M:\\pipelineTest\\elementMap\\referenceKmer\\referenceGenome\\testReference.txt";
        String taxaScoreDirS = "M:\\pipelineTest\\elementMap\\taxaScores\\";
//        String r1FileS = "/workdir/mingh/fastq/CML11_HJ7HFCCXX_L8_1.clean.fq.gz";
//        String r2FileS = "/workdir/mingh/fastq/CML11_HJ7HFCCXX_L8_2.clean.fq.gz";
//        String kmerListDirS = "/workdir/mingh/referenceKmer/16_kmerList/";
//        String referenceGenomeFileS = "/workdir/mingh/referenceKmer/referenceGenome/maize.agpV3.fa";
//        String taxaScoreDirS = "/workdir/mingh/taxaScore/"; 
        ArrayList<String> fastqList = new ArrayList();
        fastqList.add(r1FileS);
        fastqList.add(r2FileS);
        String taxonName = "CML11";
        new TaxonCpScore(referenceGenomeFileS, kmerListDirS, taxonName, fastqList, taxaScoreDirS);
    }
    
    public void generateReferenceKmers () {
        int kmerLength = 16;
//        String referenceGenomeFileS = "M:\\pipelineTest\\elementMap\\referenceKmer\\referenceGenome\\testReference.txt";
//        String kmerListDirS = "M:\\pipelineTest\\elementMap\\referenceKmer\\kmerList\\";
        String referenceGenomeFileS = "/workdir/mingh/referenceKmer/referenceGenome/maize.agpV3.fa";
        String kmerListDirS = "/workdir/mingh/referenceKmer/kmerList/";
        new CopyNumberScore().writeReferenceKmers(kmerLength, referenceGenomeFileS, kmerListDirS);
    }
    
    public void supportPipe () {
        //new Miscellaneous().getTaxaNamList();
        //new Miscellaneous().checkTaxa();
        //new Miscellaneous().getSingleCopySite();
        new Miscellaneous().getHighCoverageTaxa();
    }
    
    
    public static void main (String[] args) {
        new ElementMapGo ();
    }
}
