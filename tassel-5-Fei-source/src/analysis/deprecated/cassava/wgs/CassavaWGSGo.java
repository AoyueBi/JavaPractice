/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package analysis.deprecated.cassava.wgs;

import format.Fasta;
import format.ShortreadAlignment;
import format.ShortreadPEAlignment;
import graphcis.r.BarPlot;
import graphcis.r.BoxPlot;
import graphcis.r.CumulativeDistribution;
import graphcis.r.Histogram;
import graphcis.r.DensityPlotMultiClass;
import graphcis.r.ScatterPlot;
import java.awt.Desktop;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;
import static net.maizegenetics.analysis.chart.ChartDisplayPlugin.ChartType.Histogram;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author Fei Lu
 */
public class CassavaWGSGo {
    
    public CassavaWGSGo () {
        //this.insertSizePipe();
        //this.QCPipe();
        //this.hapMapPipe();
        this.revisionPipe();
    }
    
    public void revisionPipe () {
        //new ExtraHets();
        //new Location();
        //new Reformat();
        new BurdenRevisit();

    }
    
    public void hapMapPipe () {
        //new IBD();
        //new Heterozygosity();
        //new SpeciesIntrogression();
        //new FastqPrefilter();
        //new FilterHapMap ();
    }
    
    public void QCPipe () {
        //new Aligner();
        //new SequencingQualityPattern();
        //new ExtraDivergence();
        //new CassavaReference ();
        //new SNPDiscovery ();
        //new SNPDiscoveryFast();
        //new SNPFilter();
        //new BurdenCount();
    }
    
    public void insertSizePipe () {
        //this.sliceFastq();
        //this.alignPE();
        //this.peAlignStatistics();
        //this.mergePE();
        //this.mkContigSizeDistribution();
        //this.alignmentOverlap();
    }
    
    public void alignmentOverlap () {
        String seBWAMEM = "M:\\pipelineTest\\cassava\\wgs\\pe\\contiging\\alignment\\r1.bwa.sam";
        String seBowtie2 = "M:\\pipelineTest\\cassava\\wgs\\pe\\contiging\\alignment\\r1.bowtie2.sam";
        new CaWGSUtils().compareAligners(seBWAMEM, seBowtie2);
        String peBWAMEM = "M:\\pipelineTest\\cassava\\wgs\\pe\\contiging\\alignment\\contig.bwa.sam";
        String peBowtie2 = "M:\\pipelineTest\\cassava\\wgs\\pe\\contiging\\alignment\\contig.bowtie2.sam";
        new CaWGSUtils().compareAligners(peBWAMEM, peBowtie2);
    }
    
    public void mkContigSizeDistribution () {
        String contigFileS = "M:\\pipelineTest\\cassava\\wgs\\pe\\contiging\\contig\\contig.extendedFrags.fastq.gz";
        String pdfFileS = "M:\\pipelineTest\\cassava\\wgs\\pe\\contiging\\contig\\contigSize.pdf";
        new CaWGSUtils().mkContigSizeDistribution(contigFileS, pdfFileS);
    }
    
    public void mergePE () {
        String fastqDirS = "/workdir/mingh/fastq/";
        String mergeDirS = "/workdir/mingh/PEcontig/";
        String perlScirpt = "/workdir/mingh/merge.pl";
        new CaWGSUtils().mergePE(fastqDirS, mergeDirS, perlScirpt);
    }
    
    public void peAlignStatistics () {
        String samDirS = "M:\\pipelineTest\\cassava\\wgs\\pe\\alignmentBWAMEM\\sam\\";
        String scoreFigureDirS = "M:\\pipelineTest\\cassava\\wgs\\pe\\alignmentBWAMEM\\scoreFigure\\";
        String sizeFigureDirS = "M:\\pipelineTest\\cassava\\wgs\\pe\\alignmentBWAMEM\\sizeFigure\\";
        String reportFileS = "M:\\pipelineTest\\cassava\\wgs\\pe\\alignmentBWAMEM\\pe_report.txt";
        new CaWGSUtils().alignmentBWAMEMStatistics(samDirS, scoreFigureDirS, sizeFigureDirS, reportFileS);
    }
    
    public void alignPE () {
        String fastqDirS = "/workdir/mingh/fastqChunk";
        //String reference = "/workdir/mingh/bwa/cassavaV6_chrAndScaffoldsCombined_numeric.fa";
        String reference = "/workdir/mingh/bwa/cassavaV6.fa";
        String alignmentDirS = "/workdir/mingh/alignment";
        String perlScript = "/workdir/mingh/runBWAMEM.pl";
        new CaWGSUtils().alignFastqBWAMEM(fastqDirS, reference, alignmentDirS, perlScript);
    }
    
    public void sliceFastq () {
        String fastqDirS = "N:\\cassavaWGS\\";
        String slicedFastqDirS = "M:\\pipelineTest\\cassava\\wgs\\pe\\fastqChunk\\";
        int startIndex = 100000;
        int readNum = 5000;
        new CaWGSUtils().sliceFastq(fastqDirS, slicedFastqDirS, startIndex, readNum);
    }
   
    public static void main (String[] args) {
        new CassavaWGSGo ();
    }
}
