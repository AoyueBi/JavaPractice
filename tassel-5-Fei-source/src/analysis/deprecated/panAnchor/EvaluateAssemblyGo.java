/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package analysis.deprecated.panAnchor;

import format.Fasta;
import format.ShortreadAlignment;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import net.maizegenetics.dna.map.TagsOnGeneticMap;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import utils.IOFileFormat;



/**
 *
 * @author Fei Lu
 */
public class EvaluateAssemblyGo {
    
    public EvaluateAssemblyGo () {
        this.mo17ScaffoldPipe();
        //this.ph207PipeOrdered();
        //this.ph207PipeUnordered();
    }
    
    public void ph207PipeUnordered () {
        //this.statisticsPh207Unordered();
        //this.sortGenomePh207Unordered();
        this.alignByBowtie2();
        //this.readInBowtie2Ph207Unordered();
        //this.mkAnchorOnContigPh207Unordered();
        //this.mkSyntenyPh207Unordered();
        //this.statisticsAssemblyPh207Unordered();
    }
    
    public void statisticsAssemblyPh207Unordered () {
        String anchorOnContigFileS = "M:\\production\\panGenome\\ph207\\anchorOnContig\\anchorOnContig_ph207_unordered.txt";
        String assemblyQualityFileS = "M:\\production\\panGenome\\ph207\\anchorOnContig\\quality_ph207_unordered.txt";
        new Mo17Utils().mkStatisticsAssemblyByChr(anchorOnContigFileS, assemblyQualityFileS);
        
        String barPlotFileS = "M:\\production\\panGenome\\ph207\\anchorOnContig\\bar_ph207_unordered.txt";
        new Mo17Utils().mkBarQualityAssembly(assemblyQualityFileS, barPlotFileS, 50000, 5000000);
    }
    
    public void mkSyntenyPh207Unordered () {
        String anchorOnContigFileS = "M:\\production\\panGenome\\ph207\\anchorOnContig\\anchorOnContig_ph207_unordered.txt";
        String chrInfoFileS = "E:\\Database\\InfoFile\\ChrLenCentPosi.txt";
        String syntenyDirS = "M:\\production\\panGenome\\ph207\\synteny\\unordered\\";
        new Mo17Utils().mkSyntenyGraph(anchorOnContigFileS, chrInfoFileS, syntenyDirS, 5000, true);
    }
     
    public void mkAnchorOnContigPh207Unordered () {
        String alignFileS = "M:\\production\\panGenome\\ph207\\alignment\\togmOnPh207Unordered.saln.txt";
        String ph207OrderedLengthFileS = "M:\\production\\panGenome\\ph207\\scaffold\\contigLength_ph207_unordered.txt";
        String anchorFileS = "M:\\pav\\PhyGenMapping\\v1.paper.txt";
        String anchorOnContigFileS = "M:\\production\\panGenome\\ph207\\anchorOnContig\\anchorOnContig_ph207_unordered.txt";
        String satisticsFileS = "M:\\production\\panGenome\\ph207\\anchorOnContig\\anchorOnContig_ph207_unordered_statistics.txt";
        new Mo17Utils().mkAnchorOnContigFile(alignFileS, ph207OrderedLengthFileS, anchorFileS, anchorOnContigFileS, satisticsFileS);
    }
    
    public void readInBowtie2Ph207Unordered () {
        String bowtie2FileS = "M:\\production\\panGenome\\ph207\\alignment\\togmOnPh207Unordered.sam";
        String alignFileS = "M:\\production\\panGenome\\ph207\\alignment\\togmOnPh207Unordered.saln.txt";
        ShortreadAlignment sa = new ShortreadAlignment();
        sa.readFromBowtie2(bowtie2FileS);
        sa.writeSimpleAlignment(alignFileS, IOFileFormat.Text);
    }
    
    public void statisticsPh207Unordered () {
        String ph207OrderedRef = "M:\\production\\panGenome\\ph207\\scaffold\\PH207.unordered_scaffolds_sortByLength.fa";
        String ph207OrderedLengthFileS = "M:\\production\\panGenome\\ph207\\scaffold\\contigLength_ph207_unordered.txt";
        Fasta f = new Fasta(ph207OrderedRef);
        f.sortRecordByLengthDescending();
        f.writeLength(ph207OrderedLengthFileS);
        System.out.println(f.getSeqNumber()+"\t"+f.getTotalSeqLength());
    }
     
    public void sortGenomePh207Unordered () {
        String inputFileS = "M:\\production\\panGenome\\ph207\\scaffold\\PH207.unordered_scaffolds.fa";
        String outputFileS = "M:\\production\\panGenome\\ph207\\scaffold\\PH207.unordered_scaffolds_sortByLength.fa";
        new Mo17Utils().sortGenomeByLength(inputFileS, outputFileS);
    }
    
    public void ph207PipeOrdered () {
        //this.statisticsPh207Ordered();
        //this.sortGenome();
        this.alignByBowtie2();
        //this.readInBowtie2Ph207Ordered();
        //this.mkAnchorOnContigPh207Ordered();
        //this.mkSyntenyPh207Ordered();
        //this.statisticsAssemblyPh207OrderedByChr();
        this.statisticsAssemblyPh207OrderedByRegion();
    }
    
    public void statisticsAssemblyPh207OrderedByRegion() {
        String chrInfoFileS = "E:\\Database\\InfoFile\\ChrLenCentPosi.txt";
        String anchorOnContigFileS = "M:\\production\\panGenome\\ph207\\anchorOnContig\\anchorOnContig_ph207_ordered.txt";
        String assemblyQualityFileS = "M:\\production\\panGenome\\ph207\\anchorOnContig\\qualityByRegion_ph207_ordered.txt";
        int regionNumPerChr = 4;
        new Mo17Utils().mkStatisticsAssemblyByRegion(anchorOnContigFileS, assemblyQualityFileS, chrInfoFileS, regionNumPerChr);
    }
    
    public void statisticsAssemblyPh207OrderedByChr () {
        String anchorOnContigFileS = "M:\\production\\panGenome\\ph207\\anchorOnContig\\anchorOnContig_ph207_ordered.txt";
        String assemblyQualityFileS = "M:\\production\\panGenome\\ph207\\anchorOnContig\\quality_ph207_ordered.txt";
        new Mo17Utils().mkStatisticsAssemblyByChr(anchorOnContigFileS, assemblyQualityFileS);
        
        String barPlotFileS = "M:\\production\\panGenome\\ph207\\anchorOnContig\\bar_ph207_ordered.txt";
        new Mo17Utils().mkBarQualityAssembly(assemblyQualityFileS, barPlotFileS, 50000, 5000000);
    }
    
    public void mkSyntenyPh207Ordered () {
        String anchorOnContigFileS = "M:\\production\\panGenome\\ph207\\anchorOnContig\\anchorOnContig_ph207_ordered.txt";
        String chrInfoFileS = "E:\\Database\\InfoFile\\ChrLenCentPosi.txt";
        String syntenyDirS = "M:\\production\\panGenome\\ph207\\synteny\\ordered\\";
        new Mo17Utils().mkSyntenyGraph(anchorOnContigFileS, chrInfoFileS, syntenyDirS,5000, true);
    }
    
    public void mkAnchorOnContigPh207Ordered () {
        String alignFileS = "M:\\production\\panGenome\\ph207\\alignment\\togmOnPh207Ordered.saln.txt";
        String ph207OrderedLengthFileS = "M:\\production\\panGenome\\ph207\\scaffold\\contigLength_ph207_order.txt";
        String anchorFileS = "M:\\pav\\PhyGenMapping\\v1.paper.txt";
        String anchorOnContigFileS = "M:\\production\\panGenome\\ph207\\anchorOnContig\\anchorOnContig_ph207_ordered.txt";
        String satisticsFileS = "M:\\production\\panGenome\\ph207\\anchorOnContig\\anchorOnContig_ph207_ordered_statistics.txt";
        new Mo17Utils().mkAnchorOnContigFile(alignFileS, ph207OrderedLengthFileS, anchorFileS, anchorOnContigFileS, satisticsFileS);
    }
     
    public void readInBowtie2Ph207Ordered () {
        String bowtie2FileS = "M:\\production\\panGenome\\ph207\\alignment\\togmOnPh207Ordered.sam";
        String alignFileS = "M:\\production\\panGenome\\ph207\\alignment\\togmOnPh207Ordered.saln.txt";
        ShortreadAlignment sa = new ShortreadAlignment();
        sa.readFromBowtie2(bowtie2FileS);
        sa.writeSimpleAlignment(alignFileS, IOFileFormat.Text);
    }
    
    public void alignByBowtie2 () {
        String command = "-k 5 --very-sensitive-local";
    }
    
    public void statisticsPh207Ordered () {
        String ph207OrderedRef = "M:\\production\\panGenome\\ph207\\scaffold\\PH207.ordered_scaffolds_sortByLength.fa";
        String ph207OrderedLengthFileS = "M:\\production\\panGenome\\ph207\\scaffold\\contigLength_ph207_order.txt";
        Fasta f = new Fasta(ph207OrderedRef);
        f.sortRecordByLengthDescending();
        f.writeLength(ph207OrderedLengthFileS);
        System.out.println(f.getSeqNumber()+"\t"+f.getTotalSeqLength());
    }
    
    public void sortGenome () {
        String inputFileS = "M:\\production\\panGenome\\ph207\\scaffold\\PH207.ordered_scaffolds.fa";
        String outputFileS = "M:\\production\\panGenome\\ph207\\scaffold\\PH207.ordered_scaffolds_sortByLength.fa";
        new Mo17Utils().sortGenomeByLength(inputFileS, outputFileS);
    }
    
    public void mo17ScaffoldPipe () {
        //this.statisticsMo17();
        //this.convertToFastQ();
        //this.readInBowtie2();
        //this.mkAnchorOnContig();
        //this.mkSynteny();
        this.evaluateAnchor1();
        //this.evaluateAnchor2();
        //this.evaluteAnchorByBlast();
        this.statisticsAssembly();
        //this.anchorError();
    }
    
    public void anchorError () {
        String assemblyQualityFileS = "M:\\production\\panGenome\\mo17\\anchorOnContig\\quality.txt";
        //see those figures and come up a list
        String goodContigList = "M:\\production\\panGenome\\mo17\\evaluationAnchor\\goodContigList.txt";
        String goodContigError = "M:\\production\\panGenome\\mo17\\evaluationAnchor\\goodContigError.txt";
        new Mo17Utils().mkGoodContigError(assemblyQualityFileS, goodContigList, goodContigError);
    }
    
    public void statisticsAssembly () {
        String anchorOnContigFileS = "M:\\production\\panGenome\\mo17\\anchorOnContig\\anchorOnContig.txt";
        String assemblyQualityFileS = "M:\\production\\panGenome\\mo17\\anchorOnContig\\quality.txt";
        new Mo17Utils().mkStatisticsAssemblyByChr(anchorOnContigFileS, assemblyQualityFileS);
        
        String barPlotFileS = "M:\\production\\panGenome\\mo17\\anchorOnContig\\bar.txt";
        new Mo17Utils().mkBarQualityAssembly(assemblyQualityFileS, barPlotFileS, 50000, 1000000);
        
    }
    
    public void evaluteAnchorByBlast () {
        String assemblyQualityFileS = "M:\\production\\panGenome\\mo17\\anchorOnContig\\quality.txt";
        String mo17Ref = "M:\\production\\panGenome\\mo17\\genome\\Mo17_scaffolds_201403.fa";
        String fragmentFileS = "M:\\production\\panGenome\\mo17\\evalutetionAnchorByBlast\\goodScaffold_frag.fa";
        //new Mo17Utils().mkGoodScaffoldFragments(assemblyQualityFileS, mo17Ref, fragmentFileS);
        
        String samFileS = "M:\\production\\panGenome\\mo17\\evalutetionAnchorByBlast\\fragOnB73.sam";
        String anchorOnContigFileS = "M:\\production\\panGenome\\mo17\\anchorOnContig\\anchorOnContig.txt";
        String compareFileS = "M:\\production\\panGenome\\mo17\\evalutetionAnchorByBlast\\goodScaffoldAnchorVsAlign.txt";
        //new Mo17Utils().mkCompareAnchorAndAlignmentGoodScaffold(anchorOnContigFileS, samFileS, compareFileS);
        
        String chrInfoFileS = "E:\\Database\\InfoFile\\ChrLenCentPosi.txt";
        String syntenyDirS = "M:\\production\\panGenome\\mo17\\evalutetionAnchorByBlast\\synteny_blastDisagree\\";
        new Mo17Utils().mkDisagreementBlastSyntenyGraph(samFileS, compareFileS, chrInfoFileS, syntenyDirS);
       
    }
    
    public void evaluateAnchor2 () {
        String bowtie2FileS = "M:\\production\\panGenome\\mo17\\evaluationAnchor\\scaffoldOnB73.sam";
        String alignFileS = "M:\\production\\panGenome\\mo17\\evaluationAnchor\\scaffoldOnB73.saln.txt";
        ShortreadAlignment sa = new ShortreadAlignment();
        sa.readFromBowtie2(bowtie2FileS);
        sa.writeSimpleAlignment(alignFileS, IOFileFormat.Text);
        String chrInfoFileS = "E:\\Database\\InfoFile\\ChrLenCentPosi.txt";
        String fragmentOnChrFileS = "M:\\production\\panGenome\\mo17\\evaluationAnchor\\fragmentOnChr.txt";
        new Mo17Utils().mkFragmentOnChromosome(alignFileS, chrInfoFileS, fragmentOnChrFileS);
        
        String syntenyDirS = "M:\\production\\panGenome\\mo17\\evaluationAnchor\\synteny\\";
        new Mo17Utils().mkSyntenyGraph(fragmentOnChrFileS, chrInfoFileS, syntenyDirS,1, false);
    }
    
    public void evaluateAnchor1 () {
        //evalute with scafold162
        String mo17Ref = "M:\\production\\panGenome\\mo17\\genome\\Mo17_scaffolds_201403.fa";
        String scaffold = "scaffold162";
        String fragmentFileS = "M:\\production\\panGenome\\mo17\\evaluationAnchor\\"+scaffold+".fa";
        new Mo17Utils().mkScaffoldFragments(mo17Ref, scaffold, fragmentFileS);  
    }
    
    public void mkSynteny () {
        String anchorOnContigFileS = "M:\\production\\panGenome\\mo17\\anchorOnContig\\anchorOnContig.txt";
        String chrInfoFileS = "E:\\Database\\InfoFile\\ChrLenCentPosi.txt";
        String syntenyDirS = "M:\\production\\panGenome\\mo17\\synteny\\";
        new Mo17Utils().mkSyntenyGraph(anchorOnContigFileS, chrInfoFileS, syntenyDirS,5000, true);
    }
    
    public void mkAnchorOnContig () {
        String alignFileS = "M:\\production\\panGenome\\mo17\\alignment\\togmOnMo17.saln.txt";
        String contigLengthFileS = "M:\\production\\panGenome\\mo17\\genome\\statistics\\contigLength.txt";
        String anchorFileS = "M:\\pav\\PhyGenMapping\\v1.paper.txt";
        String anchorOnContigFileS = "M:\\production\\panGenome\\mo17\\anchorOnContig\\anchorOnContig.txt";
        String satisticsFileS = "M:\\production\\panGenome\\mo17\\anchorOnContig\\anchorOnContig_statistics.txt";
        new Mo17Utils().mkAnchorOnContigFile(alignFileS, contigLengthFileS, anchorFileS, anchorOnContigFileS, satisticsFileS);
    }
    
    public void readInBowtie2 () {
        String bowtie2FileS = "M:\\production\\panGenome\\mo17\\alignment\\togmOnMo17.sam";
        String alignFileS = "M:\\production\\panGenome\\mo17\\alignment\\togmOnMo17.saln.txt";
        ShortreadAlignment sa = new ShortreadAlignment();
        sa.readFromBowtie2(bowtie2FileS);
        sa.writeSimpleAlignment(alignFileS, IOFileFormat.Text);
    }
    
    public void convertToFastQ () {
        String togmFileS = "M:\\pav\\PhyGenMapping\\v1.togm.txt";
        String fastqFileS = "M:\\production\\panGenome\\mo17\\alignment\\v1.togm.fq";
        TagsOnGeneticMap togm = new TagsOnGeneticMap(togmFileS,FilePacking.Text);
        togm.writeFastQ(fastqFileS);
    }
    
    public void statisticsMo17 () {
        String mo17Ref = "M:\\production\\panGenome\\mo17\\genome\\Mo17_scaffolds_201403.fa";
        String mo17LengthFileS = "M:\\production\\panGenome\\mo17\\genome\\statistics\\contigLength.txt";
        Fasta f = new Fasta(mo17Ref);
        f.sortRecordByLengthDescending();
        f.writeLength(mo17LengthFileS);
        System.out.println(f.getSeqNumber()+"\t"+f.getTotalSeqLength());
    }
    
    public static void main (String[] args) {
        new EvaluateAssemblyGo();
    }
}
