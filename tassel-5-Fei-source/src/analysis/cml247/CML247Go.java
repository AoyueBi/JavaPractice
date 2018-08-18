/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package analysis.cml247;

import analysis.deprecated.panAnchor.Mo17Utils;

import format.Vertices;
import format.ShortreadAlignment;
import format.Table;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.TreeSet;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class CML247Go {
    
    
    public CML247Go () {
        //this.oldPipe();
        this.paperPipe();
    }
    
    public void paperPipe () {
        //this.selectMaizeSequence();
        //this.compareGene();
        this.QualityControl();
    }
    
    public void QualityControl () {
        new QualityControl();
    }
    
    public void compareGene () {
        new CompareGene();
    }
    
    public void selectMaizeSequence () {
        //this.sliceScaffold();
        //this.mkSAFile();
        //this.mkSorghumRatioFile();
        //this.selectMaizeAssembly();
        //this.assemblyStatistics();
    }
    
    private void sliceScaffold () {
        String intputFastaFileS = "M:\\production\\panGenome\\cml247\\nrgene\\source\\CML247_maize.fasta";
        String sliceFastaFileS = "E:\\Research\\cml247\\sorghumCleaning\\CML247_maize_slice.fasta";
        int sliceSize = 200;
        new AssemblyEvaluationUtils().sliceFasta(intputFastaFileS, sliceFastaFileS, sliceSize);
    } 
    
    public void mkSAFile () {
        String samFileS = "E:\\Research\\cml247\\sorghumCleaning\\CML247_maize_slice.sam";
        String saFileS = "E:\\Research\\cml247\\sorghumCleaning\\CML247_maize_slice.sa.bin";
        ShortreadAlignment sa = new ShortreadAlignment();
        sa.readFromBowtie2(samFileS);
        sa.sortByQuery();
        sa.writeSimpleAlignment(saFileS, IOFileFormat.Binary);
    }
    
    public void mkSorghumRatioFile () {
        String saFileS = "E:\\Research\\cml247\\sorghumCleaning\\CML247_maize_slice.sa.bin";
        String intputFastaFileS = "M:\\production\\panGenome\\cml247\\nrgene\\source\\CML247_maize.fasta";
        String sorghumRatioFileS = "E:\\Research\\cml247\\sorghumCleaning\\sorghumRatio.txt";
        ShortreadAlignment sa = new ShortreadAlignment(saFileS, IOFileFormat.Binary);
        new AssemblyEvaluationUtils().mkSorghumRatioFile(sa, intputFastaFileS, sorghumRatioFileS);
    }
    
    public void selectMaizeAssembly () {
        String sorghumRatioFileS = "E:\\Research\\cml247\\sorghumCleaning\\sorghumRatio.txt";
        String inputFastaFileS = "M:\\production\\panGenome\\cml247\\nrgene\\source\\CML247_maize.fasta";
        String outputFastaFileS = "E:\\Research\\cml247\\fasta\\cml247.fas";
        new AssemblyEvaluationUtils().selectMaizeAssembly(sorghumRatioFileS, inputFastaFileS, outputFastaFileS);
    }
    
    public void assemblyStatistics () {
        String inputFastaFileS = "E:\\Research\\cml247\\fasta\\cml247.fas";
        String staFileS = "E:\\Research\\cml247\\fasta\\cml247.fas.statistics.txt";
        new AssemblyEvaluationUtils().mkAssemblyStatistics(inputFastaFileS, staFileS);
    }
    
    public void oldPipe () {
        //this.checkLibraryPipe(); //checking if sequence file are contaminated by sorghorm
        this.sorghumControlPipe();
        //this.anchorControlPipe();
        //this.graphPipe(); //specially for DISCOVAR
    }
    
    public void graphPipe () {
        //this.convertDiscovarSrcFileS();
        //this.convertToGraphFile();
        //this.mkVertexPath();
        //this.mkVertexDegree();
        //this.mkVertexPathLength();
        //this.selectVertexPathLength();
        //this.sliceMaizeEdge();
        //this.alignMaizeEdgeSlice();
        //this.mkSimpleAlignEdgeSlice();
        //this.mkSliceRepeatPath();
        //this.mkPathRepeatValue();
    }
    
    public void mkPathRepeatValue () {
        String vertexRepeatPathFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\vertex\\edge.repeat.path.txt";
        String vertexPathRepeatValueFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\vertex\\edge.path.repeatValue.txt";
        new AssemblyEvaluationUtils().mkPathRepeatValue(vertexRepeatPathFileS, vertexPathRepeatValueFileS);
    }
    
    public void mkSliceRepeatPath () {
        String samFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\alignment\\edge_sorted_maizeSeqOnRepeat_slice_200_0.sam";
        String vertexPathFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\graph\\edge.path.txt";
        String vertexRepeatPathFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\vertex\\edge.repeat.path.txt";
        new AssemblyEvaluationUtils().mkVertexRepeatPathNum(samFileS, vertexPathFileS, vertexRepeatPathFileS);
    }
    
    public void mkSimpleAlignEdgeSlice () {
        String samFileS = "/SSD/mingh/edgeSlice.sam";
        String alignmentFileS = "/SSD/mingh/edge_sorted_maizeSeqOnRepeat_slice_200_0.sa";
        new AssemblyEvaluationUtils().convertSamToSimpleAlignment(samFileS, alignmentFileS);
    }
    
    public void alignMaizeEdgeSlice () {
        //should use "-x maizerepeat -f edge_sorted_slice_200_0.fas -M --very-sensitive-local -S edge_sorted_maizeSeqOnRepeat_slice_200_0.sam"
    }
    
    public void sliceMaizeEdge () {
        //String sortedFastaFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\fasta\\edge_sorted.fas";
        //String maizeSorghumRatioFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\maizeSorghumRatio\\edge_maizeSorghumRatio.txt";
        //String sliceFastaFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\fasta\\edge_sorted_maize_slice_200_2000.fas";
        String sortedFastaFileS = "/SSD/mingh/edge_sorted.fas";
        String maizeSorghumRatioFileS = "/SSD/mingh/edge_maizeSorghumRatio.txt";
        String sliceFastaFileS = "/SSD/mingh/edge_sorted_slice_200_0.fas";
        int sliceSize = 200;
        int seqLengthCut = 0;
        float ratioCut = (float)0.97;
        new AssemblyEvaluationUtils().sliceMaizeContig(sortedFastaFileS, sliceFastaFileS, maizeSorghumRatioFileS, sliceSize, seqLengthCut, ratioCut);
    }
    
    public void selectVertexPathLength () {
        String vertexPathLengthFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\vertex\\edge.path.length.txt";
        String subVertexPathLengthFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\vertex\\edge.path.length.sub.txt";
        int size = 100000;
        new AssemblyEvaluationUtils().selectVertexPathLength(vertexPathLengthFileS, size, subVertexPathLengthFileS);
    }
    
    public void mkVertexPathLength () {
        String vertexPathFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\graph\\edge.path.txt";
        String indexMapFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\fasta\\edge_indexMap.txt";
        String vertexPathLengthFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\vertex\\edge.path.length.txt";
        new AssemblyEvaluationUtils().mkVertexPathLength(vertexPathFileS, indexMapFileS, vertexPathLengthFileS);
    }
    
    public void mkVertexDegree () {
        String graphFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\graph\\edge.vtx.txt";
        String vertexDegreeFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\vertex\\vertex.degree.txt";
        int degreeCut = 2;
        new AssemblyEvaluationUtils().mkVertexDegreeFileS(graphFileS, vertexDegreeFileS, degreeCut);
    }
    
    public void mkVertexPath () {
        String srcFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\graph\\edge.src";
        String vertexPathFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\graph\\edge.path.txt";
        new AssemblyEvaluationUtils().mkVertexPathFileS(srcFileS, vertexPathFileS);
    }
     
    public void convertToGraphFile () {
        String discovarFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\graph\\edge.src";
        String graphFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\graph\\edge.vtx.txt";
        new AssemblyEvaluationUtils().convertDiscovarToGraph(discovarFileS, graphFileS);
    }
    
    public void convertDiscovarSrcFileS () {
        String srcFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\source\\a.lines.src";
        String edgeSrcFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\graph\\edge.src";
        String indexMapFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\fasta\\edge_indexMap.txt";
        new AssemblyEvaluationUtils().convertDiscovarSrcFileS(srcFileS, edgeSrcFileS, indexMapFileS);
    }
    
    public void anchorControlPipe () {
        //this.selectLargeContig();
        //this.alignAnchorByBowtie2();
        //this.readInAnchorOnContigSam();
        //this.mkAnchorOnContigFile();
        this.mkSyntenyGraph();
        //this.anchorStatisticsAssembly();

    }
    
    
    public void anchorStatisticsAssembly () {
        String anchorOnContigFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\anchorOnContig\\anchorOnContig.txt";
        String assemblyQualityFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\anchorOnContig\\contig_quality.txt";
        String indexMapFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\fasta\\contig_indexMap.txt";
        new AssemblyEvaluationUtils().mkStatisticsAssemblyByChr(anchorOnContigFileS, assemblyQualityFileS, indexMapFileS);
    }
    
    public void mkSyntenyGraph () {
        String anchorOnContigFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\anchorOnContig\\anchorOnContig.txt";
        String chrInfoFileS = "E:\\Database\\InfoFile\\ChrLenCentPosi.txt";
        String syntenyDirS = "M:\\production\\panGenome\\cml247\\currentApproach\\synteny\\";
        new AssemblyEvaluationUtils().mkSyntenyGraph(anchorOnContigFileS, chrInfoFileS, syntenyDirS, 5000, true);
    }
    
    public void mkAnchorOnContigFile () {
        String alignFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\alignment\\togmOnContig20k.sa.txt";
        String indexMapFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\fasta\\contig_indexMap.txt";
        //String anchorFileS = "M:\\pav\\PhyGenMapping\\v1.paper.txt";
        //String anchorFileS = "E:\\Research\\geneticMapping\\identifyPAV\\anchor\\v1.togm.BlastValidation.txt";
        String anchorFileS = "E:\\Research\\geneticMapping\\identifyPAV\\anchor\\v1.togm.Bowtie2Validation.txt";
        //String anchorFileS = "E:\\Research\\geneticMapping\\identifyPAV\\anchor\\v1.togm.BothValidation.txt";
        String anchorOnContigFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\anchorOnContig\\anchorOnContig.txt";
        String satisticsFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\anchorOnContig\\anchorOnContig_statistics.txt";
        new AssemblyEvaluationUtils().mkAnchorOnContigFile(alignFileS, indexMapFileS, anchorFileS, anchorOnContigFileS, satisticsFileS);
    }
    
    public void readInAnchorOnContigSam () {
        String bowtie2FileS = "M:\\production\\panGenome\\cml247\\currentApproach\\alignment\\anchorOnContig20k.sam";
        String alignFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\alignment\\togmOnContig20k.sa.txt";
        ShortreadAlignment sa = new ShortreadAlignment();
        sa.readFromBowtie2(bowtie2FileS);
        sa.sortByHitAndPos();
        sa.writeSimpleAlignment(alignFileS, IOFileFormat.Text);
    }
     
    public void alignAnchorByBowtie2 () {
        //should be -k 2, then the only best alignment of anchors can be used to evaluate assembly. 
        String command = "-k 2 --very-sensitive-local";
    }
    
    public void selectLargeContig () {
        String sourceFastaFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\fasta\\contig_sorted.fas";
        String desFastaFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\fasta\\contig_sorted_20k.fas";
        String maizeSorghumRatioFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\maizeSorghumRatio\\contig_maizeSorghumRatio.txt";
        int minSize = 20000;
        new AssemblyEvaluationUtils().selectLargeContigs(sourceFastaFileS, desFastaFileS, maizeSorghumRatioFileS, minSize);
    }
    
    public void sorghumControlPipe () {
        //this.renameSortFasta();
        //this.basicAssemblyStatistics();
        //this.sliceContigs();
        //this.mkSimpleAlignment();
        this.mkMaizeSorghumRatio();
    }
    
    public void mkMaizeSorghumRatio () {
        String alignmentFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\alignment\\contig_sorted_slice_100_2000.sa";
        String indexMapFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\fasta\\contig_indexMap.txt";
        String maizeSorghumRatioFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\maizeSorghumRatio\\contig_maizeSorghumRatio.txt";
        new AssemblyEvaluationUtils().mkMaizeSorghumRatioFileS(alignmentFileS, indexMapFileS, maizeSorghumRatioFileS);
    }
    
    public void mkSimpleAlignment () {
        String samFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\alignment\\contig_maizeSorghum.sam";
        String alignmentFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\alignment\\contig_sorted_slice_100_2000.sa";
        //String samFileS = "/workdir/mingh/contig_maizeSorghum.sam";
        //String alignmentFileS = "/workdir/mingh/contig_sorted_slice_50_1000.sa";
        new AssemblyEvaluationUtils().convertSamToSimpleAlignment(samFileS, alignmentFileS);
    }
    
     /**
     * For sorghum contamination detection
     */
    public void sliceContigs () {
        String sortedFastaFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\fasta\\contig_sorted.fas";
        String sliceFastaFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\fasta\\contig_sorted_slice_100_2000.fas";
        //String sortedFastaFileS = "/workdir/mingh/contig_sorted.fas";
        //String sliceFastaFileS = "/workdir/mingh/contig_sorted_slice_50_1000.fas";
        int sliceSize = 100;
        int seqLengthCut = 2000;
        new AssemblyEvaluationUtils().sliceFasta(sortedFastaFileS, sliceFastaFileS, sliceSize, seqLengthCut);
    }
    
    public void basicAssemblyStatistics () {
        String indexMapFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\fasta\\contig_indexMap.txt";
        String statisticFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\fasta\\basicStatistics.txt";
        new AssemblyEvaluationUtils().basicAssemblyStatistics(indexMapFileS, statisticFileS);
    }
    
    public void renameSortFasta () {
        String oriFastaFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\source\\CML247_maize.fasta";
        String desFastaFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\fasta\\contig_sorted.fas";
        String indexMapFileS = "M:\\production\\panGenome\\cml247\\currentApproach\\fasta\\contig_indexMap.txt";
        
        //String oriFastaFileS = "/workdir/mingh/contig.fasta";
        //String desFastaFileS = "/workdir/mingh/contig_sorted.fas";
        //String indexMapFileS = "/workdir/mingh/contig_indexMap.fas";
        new AssemblyEvaluationUtils().renameSortFasta(oriFastaFileS, desFastaFileS, indexMapFileS);
    }
    
    
    public void checkLibraryPipe() {
        //this.mergeMaizeSorghum();
        //this.checkSeqByBlast();
        //this.mkQuerys();
        //this.mkPacBioQuery();
        //this.mergeQuery();
        //this.mkRatioTable();
    }
    
    public void mkRatioTable () {
        String checkSamFileS = "M:\\cml247Check\\queryCheck.sam";
        String checkRatioTable = "M:\\cml247Check\\checkRatioTable.txt";
        new CheckLibraryUtils().mkRatioTable(checkSamFileS, checkRatioTable);
    }
    
    public void mergeQuery () {
        String queryDir = "M:\\cml247Check\\query\\querys\\";
        String queryFileS = "M:\\cml247Check\\query\\query.fa";
        new CheckLibraryUtils().mergeQuery(queryDir, queryFileS);
    }
    
    public void mkPacBioQuery () {
        String pacBioDirS = "M:\\production\\panGenome\\cml247\\data\\pacBio\\Maize_CML247_Library2\\";
        String queryDir = "M:\\cml247Check\\query\\querys\\";
        new CheckLibraryUtils().mkPacBioQuery(pacBioDirS, queryDir);
    }
    
    public void mkQuerys () {
        String fastqDir = "M:\\cml247Check\\otherflowcell\\hapmapV2\\";
        String queryDir = "M:\\cml247Check\\query\\querys\\";
        new CheckLibraryUtils().mkQueries(fastqDir, queryDir);
    }
    
    public void checkSeqByBlast () {
        String fastqDir = "N:\\panGenome\\CML247\\mp\\hiseq\\H94Y1ADXX\\";
        String queryFileS = "M:\\cml247Check\\queryResult\\query.fa";
        String resultFileS = "M:\\cml247Check\\queryResult\\result.txt";
        String alignmentFileS = "M:\\cml247Check\\queryResult\\alignment.txt";
        String dataBase = "M:\\cml247Check\\ref\\maizeSorghum.fa";
        new CheckLibraryUtils().checkSeqByBlast(fastqDir, queryFileS, dataBase, alignmentFileS,resultFileS);
    }
    
    public void mergeMaizeSorghum () {
        String maizeRef = "M:\\cml247Check\\ref\\ZmB73_RefGen_v2.fa";
        String sorghumRef = "M:\\cml247Check\\ref\\sbi1.fasta";
        String mergeRef = "M:\\cml247Check\\ref\\maizeSorghum.fa";
        new CheckLibraryUtils().mergeReference(maizeRef, sorghumRef, mergeRef);
    }
    
    public static void main (String[] args) {
        new CML247Go ();
    }
}
