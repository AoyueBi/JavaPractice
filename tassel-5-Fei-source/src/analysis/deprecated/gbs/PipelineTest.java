/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.gbs;

import format.Bins;
import format.Table;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Arrays;
import net.maizegenetics.analysis.data.MigrateHDF5FromT4T5;
import net.maizegenetics.analysis.gbs.AnnotateTOPM;
import net.maizegenetics.analysis.gbs.ContigPETagCountPlugin;
import net.maizegenetics.analysis.gbs.FastqToPETagCountPlugin;
import net.maizegenetics.analysis.gbs.MergePETagCountPlugin;
import net.maizegenetics.analysis.gbs.pana.PanASplitTBTPlugin;
import net.maizegenetics.analysis.gbs.QseqToPETagCountPlugin;
import net.maizegenetics.analysis.gbs.SimpleGenotypeSBit;
import net.maizegenetics.analysis.gbs.TagAgainstAnchor;
import net.maizegenetics.analysis.gbs.TagAgainstAnchorHypothesis;
import net.maizegenetics.analysis.gbs.TagAgainstAnchorPlugin;
import net.maizegenetics.analysis.gbs.TagBlockPosition;
import net.maizegenetics.analysis.gbs.pana.PanASplitTagBlockPosPlugin;
import net.maizegenetics.dna.map.PETagsOnPhysicalMapV3;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.map.TagGWASMap;
import net.maizegenetics.dna.map.TagMappingInfoV3.Aligner;
import net.maizegenetics.dna.map.TagsOnPhysicalMapV3;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.tag.PETagCounts;
import net.maizegenetics.dna.tag.TagCounts;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import net.maizegenetics.dna.tag.TagsByTaxaByteHDF5TagGroups;


/**
 *
 * @author fl262
 */
public class PipelineTest {
    
    public PipelineTest () {
        geneticMappingPipeline ();
        //PEPipeline ();
        //V3TestPipeline();
    }
    
    public void V3TestPipeline () {
        //this.mkSmallTagCount();
        //this.mkFastq();
        //this.mkFasta();
        //this.initializeTOPMHDF5();
        //this.annotateTOPMHDF5WithAligner();
        
        /*not tested yet*/
        //this.hypothesisGeneticMapping();
        //this.annotateTOPMHDF5WithGMGW();
        //this.predictHypothesis();
        //this.importPrediction();
        //this.callSNPDiscovery();
    }
    
    public void callSNPDiscovery () {
        String topmFileS = "M:/pipelineTest/GBSV3test/topm/ini.topm.h5";
        String pedigreeFileS = "M:/pipelineTest/GBSV3test/keyfile/AllZeaPedigree2012oct01C.txt";
        String tbtFileS = "M:/pipelineTest/GBSV3test/tbt/small4096TBTHDF5_mergedtaxa_pivot_20120921.h5";
        String referenceFileS = "M:/Database/maizeReference/chr10.fasta";
        String hapmapFileS = "M:/pipelineTest/GBSV3test/hapmap/myGBSGenos.chr+.hmp.txt";
        
        String arguments = "-i " + tbtFileS + " -m " + topmFileS + " -o " + hapmapFileS + " -p " + pedigreeFileS + " -ref " + referenceFileS;
        arguments = arguments + " -mxSites 10000 -mnF 0.8 -mnMAF 0.001 -mnMAC 10 -mnLCov 0.1 -cF -sC 10 -eC 10";
        String[] args = arguments.split(" ");
        //DiscoverySNPCallerV3Plugin d = new DiscoverySNPCallerV3Plugin ();
        //d.setParameters(args);
		//d.performFunction(null); 
    }
    
    public void importPrediction () {
        String topmFileS = "M:/pipelineTest/GBSV3test/topm/ini.topm.h5";
        String predictDirS = "M:/pipelineTest/GBSV3test/wekaPrediction/output/";
        String indexDirS = "M:/pipelineTest/GBSV3test/wekaPrediction/index/";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(topmFileS);
        AnnotateTOPM anno = new AnnotateTOPM(topm);
        Aligner priorityAligner = Aligner.Bowtie2;
        anno.annotateBestMappingImport(indexDirS, predictDirS, priorityAligner);

    }
    
    public void predictHypothesis () {
        String topmFileS = "M:/pipelineTest/GBSV3test/topm/ini.topm.h5";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(topmFileS);
        //String wekaFileS = "C:/Program Files/Weka-3-6/weka.jar";
        String modelFileS = "E:/Research/gbsv3/modelDevelop/262144Tags_RandomForest.model";
        String tagCountFileS = "M:/pipelineTest/GBSV3test/tagCount/small.cnt"; // to pull out count attribute
        String inputDirS = "M:/pipelineTest/GBSV3test/wekaPrediction//input/";
        String indexDirS = "M:/pipelineTest/GBSV3test/wekaPrediction/index/";
        String outputDirS = "M:/pipelineTest/GBSV3test/wekaPrediction/output/";
        AnnotateTOPM anno = new AnnotateTOPM(topm);
        anno.annotateBestMappingPredict(modelFileS, tagCountFileS, inputDirS, indexDirS, outputDirS);
    }
    
    public void annotateTOPMHDF5WithGMGW () {
        String inputFileS = "M:/pipelineTest/GBSV3test/topm/ini.topm.h5";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(inputFileS);
        AnnotateTOPM anno = new AnnotateTOPM (topm);
        String TOGMFileS = "M:/pav/PhyGenMapping/v1.togm.txt";
        anno.annotateWithGMGW(TOGMFileS, 1);
    }
    
    public void hypothesisGeneticMapping() {
        //on pc
        String topmFileS = "M:/pipelineTest/GBSV3test/topm/ini.topm.h5";
        String hapMapHDF5 = "M:/pipelineTest/GBSV3test/anchor/GBS27.small1024.sBit.h5"; //1024 sites
        String tbtHDF5 = "M:/pipelineTest/GBSV3test/tbt/small4096TBTHDF5_mergedtaxa_pivot_20120921.h5"; // 4096 tags
        TagAgainstAnchorHypothesis taah = new TagAgainstAnchorHypothesis (hapMapHDF5, tbtHDF5, topmFileS, "Bowtie2", 10, 0.5, 20, 1);
        
        //on server
        //String topmFileS = "/workdir/mingh/small65536.V3.TOPM.h5";
        //String hapMapHDF5 = "/workdir/mingh/GBS27.sBit.h5"; //1024 sites
        //String tbtHDF5 = "/workdir/mingh/mergeTBTHDF5_mergedtaxa_pivot_20120921.h5"; // 4096 tags
        //TagAgainstAnchorHypothesis taah = new TagAgainstAnchorHypothesis (hapMapHDF5, tbtHDF5, topmFileS, "Bowtie2", 1000, 0.000001, 20, -1);
    }
    
    public void annotateTOPMHDF5WithAligner () {
        String inputFileS = "M:/pipelineTest/GBSV3test/topm/ini.topm.h5";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(inputFileS);
        AnnotateTOPM anno = new AnnotateTOPM (topm);
        String bowtie2SamFileS = "M:/pipelineTest/GBSV3test/alignment/small_bowtie2-K5.sam";
        anno.annotateWithBowtie2(bowtie2SamFileS, 5);
        String bwaSamFileS = "M:/pipelineTest/GBSV3test/alignment/small_bwa-N5.sam";
        anno.annotateWithBWA(bwaSamFileS, 5);
        ////String blastFileS = "M:/pipelineTest/GBSV3test/alignment/small_blast.m8.txt";
        ////anno.annotateWithBLAST(blastFileS, 5); //In practice, the file would too large. So just split the fasta to small pieces
        String blastDirS = "M:/pipelineTest/GBSV3test/alignment/blastOut/";
        anno.annotateWithBlastFromDir(blastDirS, 5);
        String bwaMemSamFileS = "M:/pipelineTest/GBSV3test/alignment/small_bwa-mem-a.sam";
        anno.annotateWithBWAMEM(bwaMemSamFileS, 5);
        String PETOPMFileS = "M:/production/pe/ptopm/PE.topm";
        anno.annotateWithPE(PETOPMFileS, 5);
    }
    
    public void initializeTOPMHDF5 () {
        String inputFileS = "M:/pipelineTest/GBSV3test/tagCount/small.cnt";
        String outputFileS = "M:/pipelineTest/GBSV3test/topm/ini.topm.h5";
        TagCounts tc = new TagCounts(inputFileS, FilePacking.Byte);
        TagsOnPhysicalMapV3.createFile(tc, outputFileS);
    }
    
    /**
     * Make FASTA file and do alignment using Blast, -m 8 -e 1e-10
     */
    public void mkFasta () {
        String inputFileS = "M:/pipelineTest/GBSV3test/tagCount/small.cnt";
        String outputFileS = "M:/pipelineTest/GBSV3test/alignment/small.fa";
        TagCounts tc = new TagCounts(inputFileS, FilePacking.Byte);
        tc.toFASTA(outputFileS);
    }
    
    /**
     * Make Fastq file and do alignment using Bowtie2 and BWA
     * Alignment files should be in the alignment folder
     */
    public void mkFastq () {
        String inputFileS = "M:/pipelineTest/GBSV3test/tagCount/small.cnt";
        String outputFileS = "M:/pipelineTest/GBSV3test/alignment/small.fq";
        new V3Utils ().convertTagCount2Fastq(inputFileS, outputFileS);
    }
    
    public void mkSmallTagCount () {
        String inputFileS = "M:/pipelineTest/GBSV3test/tagCount/434GFAAXX_s_4.cnt";
        String outputFileS = "M:/pipelineTest/GBSV3test/tagCount/small.cnt";
        new V3Utils ().mkSmallTagCountsFile(inputFileS, outputFileS, 10001, 500);
    }
    
    public void PEPipeline () {
        //this.parseQseq();
        //this.parseFastq();
        //this.mergePETagCounts();
        //this.contigPETagCounts();
        //this.mkFastaFileFromPE();
        //this.mkPEAlignment();
    }
    
    public void mkPEAlignment () {
        String fFastaFileS = "M:/production/pe/alignment/f.fasta.txt";
        String bFastaFileS = "M:/production/pe/alignment/b.fasta.txt";
        String fSamFileS = "M:/production/pe/alignment/f.k5.sam.gz";
        String bSamFileS = "M:/production/pe/alignment/b.k5.sam.gz";
        String a = null;
        PETagsOnPhysicalMapV3 ptopm = new PETagsOnPhysicalMapV3 (fFastaFileS, bFastaFileS, fSamFileS, bSamFileS);
        String PETOPM = "M:/production/pe/ptopm/PE.topm";
        ptopm.writeBinaryFile(PETOPM);
        ptopm = new PETagsOnPhysicalMapV3(PETOPM);
    }
    
    public void mkFastaFileFromPE () {
        String infileS = "M:/production/pe/mergePETagCounts/merge.con.pe.cnt";
        PETagCounts ptc = new PETagCounts (infileS, FilePacking.Byte);
        String fFastaFileS = "M:/production/pe/alignment/f.fasta.txt";
        String bFastaFileS = "M:/production/pe/alignment/b.fasta.txt";
        ptc.mkFastaFile(fFastaFileS, bFastaFileS);
    }
    
    public void contigPETagCounts () {
        String infileS = "M:/production/pe/mergePETagCounts/merge.pe.cnt";
        String outfileS = "M:/production/pe/mergePETagCounts/merge.con.pe.cnt";
        String arguments = "-i " + infileS + " -o " + outfileS;
        String[] args = arguments.split(" ");
        ContigPETagCountPlugin m = new ContigPETagCountPlugin();
        m.setParameters(args);
		m.performFunction(null);
    }
    
    public void mergePETagCounts () {
        String inputDirS = "M:/production/pe/PETagCounts/";
        String outfileS = "M:/production/pe/mergePETagCounts/merge.pe.cnt";
        String arguments = "-i " + inputDirS + " -o " + outfileS;
        String[] args = arguments.split(" ");
        MergePETagCountPlugin m = new MergePETagCountPlugin();
        m.setParameters(args);
		m.performFunction(null); 
    }
    
    public void parseFastq () {
        String infile1 = "M:/production/pe/Illumina/test/ImputationP15_1_1_fastq.txt";
        String infile2 = "M:/production/pe/Illumina/test/ImputationP15_1_2_fastq.txt";
        
        String keyfile = "M:/production/pe/key/ImputationP15_key.txt";
        String outputDirS = "M:/production/pe/PETagCounts/";
        String arguments = "-iF " + infile1 + " -iB " + infile2 + " -k " + keyfile + " -e ApekI -l 8 -o " + outputDirS;
        String[] args = arguments.split(" ");
        FastqToPETagCountPlugin q = new FastqToPETagCountPlugin();
        q.setParameters(args);
		q.performFunction(null);     
    }
    
    public void parseQseq () {
        String infile1 = "M:/production/pe/Illumina/81546ABXX_8_1_qseq.txt";
        String infile2 = "M:/production/pe/Illumina/81546ABXX_8_2_qseq.txt";
        String keyfile = "M:/production/pe/key/81546ABXX_key.txt";
        String outputDirS = "M:/production/pe/PETagCounts/";
        String arguments = "-iF " + infile1 + " -iB " + infile2 + " -k " + keyfile + " -e ApekI -o " + outputDirS;
        String[] args = arguments.split(" ");
        QseqToPETagCountPlugin q = new QseqToPETagCountPlugin();
        q.setParameters(args);
		q.performFunction(null);     
    }
    
    
    public void geneticMappingPipeline () {
        //**********Anchor*****************//
        //this.migrateHDF5ToTASSEL5();
        //this.mkSmallAnchorHDF5();//alternative
        //this.transformAnchorHDF5();
        //*********************************//
        
        //**********Tag********************//
        //this.mkSmallTBTHDF5();
        //this.splitTBTHDF5();
        //*********************************//
        
        //**********TagBlock***************//
        //this.mkTBTTagBlockFromTOPMV3(); //For SNPs produced from V3 pipeline
        //this.mkTBTTagBlockFromTOPM(); //For GBS build version 1 & 2
        //this.splitTagBlock();
        //*********************************//
        
        //*********Mapping*****************//
        this.geneticMapping();
        //this.geneticMappingPlugin();
        //*********************************//
        
        //*******Model training************//
        //this.mkTagGWASMap();
    }
    
    public void mkTagGWASMap () {
        String gwasMappingFileS = "M:/pipelineTest/geneticMapping/result/outfile.txt";
        String topmFileS = "M:/production/v3gbs/topm/bowtie2.topm.h5";
        String tagCountFileS = "M:/production/v3gbs/tagCount/AllZeaMasterTags_c10_20120606.cnt";
        String tagGWASMapFileS = "M:/pipelineTest/geneticMapping/tagAttributeMap/tagGwasMap.h5";
        new TagGWASMap(gwasMappingFileS, topmFileS, tagCountFileS, tagGWASMapFileS);
    }
    
    public void geneticMappingPlugin () {
        TagAgainstAnchorPlugin taa = new TagAgainstAnchorPlugin();
        String hapMapHDF5 = "M:/pipelineTest/geneticMapping/anchor/GBS27.small1024.sBit.h5";
        String tbtHDF5 = "M:/pipelineTest/geneticMapping/tbt/small4096TBTHDF5_mergedtaxa_pivot_20120921.h5";
        String blockFileS = "M:/pipelineTest/geneticMapping/tagBlock/block.tbp";
        String outfileS = "M:/pipelineTest/geneticMapping/result/outfile.txt";
        
        //String arguments = "-pc " + 1 + " -t " + tbtHDF5;
        String arguments = "-g " + hapMapHDF5 + " -t " + tbtHDF5 + " -b " + blockFileS + " -o " + outfileS + " -m 20 -c max -s 65536 -cs 0 -ce 1";
        
        String[] args = arguments.split(" ");
		taa.setParameters(args);
		taa.performFunction(null);
    }
    
    public void geneticMapping () {
        String hapMapHDF5 = "M:/pipelineTest/geneticMapping/anchor/GBS27.small1024.sBit.h5";
        String tbtHDF5 = "M:/pipelineTest/geneticMapping/tbt/small4096TBTHDF5_mergedtaxa_pivot_20120921.h5";
        String blockFileS = "M:/pipelineTest/geneticMapping/tagBlock/block.tbp";
        String outfileS = "M:/pipelineTest/geneticMapping/result/outfile.txt";

        TagAgainstAnchor taa = new TagAgainstAnchor(hapMapHDF5, tbtHDF5, blockFileS, outfileS, 0.000001, 20, 1, 65536);
    }
    
    public void splitTagBlock () {
        String blockFileS = "M:/pipelineTest/geneticMapping/tagBlock/block.tbp";
        String outputDirS = "M:\\pipelineTest\\geneticMapping\\tagBlock\\subTBP";
        String chunkSize = "48";
        String arguments = "-i " + blockFileS + " -s " + chunkSize + " -o " + outputDirS;
        String[] args = arguments.split(" ");
        PanASplitTagBlockPosPlugin stb = new PanASplitTagBlockPosPlugin();
        stb.setParameters(args);
        stb.performFunction(null);
    }
    
    public void mkTBTTagBlockFromTOPM () {
        String tbtHDF5 = "M:/pipelineTest/geneticMapping/tbt/small4096TBTHDF5_mergedtaxa_pivot_20120921.h5";
        String topmFileS = "M:\\production\\geneticMapping\\tagBlock\\AllZeaGBSv2.6ProdTOPM_20130605.topm.h5";
        String blockFileS = "M:/pipelineTest/geneticMapping/tagBlock/block.tbp";
        int TOPMVersionValue = 1;
        TagBlockPosition tbp = new TagBlockPosition(tbtHDF5, topmFileS, TOPMVersionValue);
        tbp.writeTagBlockPosition(blockFileS);
    }
    
    public void mkTBTTagBlockFromTOPMV3 () {
        String tbtHDF5 = "M:/pipelineTest/geneticMapping/tbt/small4096TBTHDF5_mergedtaxa_pivot_20120921.h5";
        String topmFileS = "M:/pipelineTest/GBSV3test/topm/ini.topm.h5";
        String blockFileS = "M:/pipelineTest/geneticMapping/tagBlock/block.tbp";
        TagBlockPosition tbp = new TagBlockPosition(tbtHDF5, topmFileS, "Bowtie2");
        tbp.writeTagBlockPosition(blockFileS);
    }
    
    public void splitTBTHDF5 () {
        String inputTBTS = "M:/pipelineTest/geneticMapping/tbt/small4096TBTHDF5_mergedtaxa_pivot_20120921.h5";
        String outputDirS = "M:\\pipelineTest\\geneticMapping\\tbt\\subTBT\\";
        String chunkSize = "48";
        String arguments = "-i " + inputTBTS + " -s " + chunkSize + " -o " + outputDirS;
        String[] args = arguments.split(" ");
        PanASplitTBTPlugin sbp = new PanASplitTBTPlugin();
        sbp.setParameters(args);
        sbp.performFunction(null);
    }
    
    public void mkSmallTBTHDF5 () {
        String inputTBTS = "M:/production/geneticMapping/tbt/mergeTBTHDF5_mergedtaxa_pivot_20120921.h5";
        String outputTBTS = "M:/pipelineTest/geneticMapping/tbt/small4096TBTHDF5_mergedtaxa_pivot_20120921.h5";
        TagsByTaxaByteHDF5TagGroups tbt = new TagsByTaxaByteHDF5TagGroups (inputTBTS);
        int[] tagIndex = new int[4096];
        for (int i = 0; i < tagIndex.length; i++) {
            //tagIndex[i] = (int)(tbt.getTagCount()*Math.random()); //random chosen
            tagIndex[i] = i; //first N chosen
        }
        Arrays.sort(tagIndex);
        tbt.writeDistFile(outputTBTS, tagIndex);
    }
    
    public void transformAnchorHDF5 () {
        String hapMapInputFileS = "M:/pipelineTest/geneticMapping/anchor/GBS27.small1024.T5.imp.hmp.h5";
        String sBitFileS = "M:/pipelineTest/geneticMapping/anchor/GBS27.small1024.sBit.h5";
        new SimpleGenotypeSBit(hapMapInputFileS, sBitFileS);
        
    }
    
    public void mkSmallAnchorHDF5 () {//alternative
        String hapMapInputS = "M:/pipelineTest/geneticMapping/anchor/AllZeaGBSv27i3b.T5.imp.hmp.h5";
        String hapMapOutputS = "M:/pipelineTest/geneticMapping/anchor/GBS27.small1024.T5.imp.hmp.h5";
        GenotypeTable gt = ImportUtils.readGuessFormat(hapMapInputS);
        int[] index = new int[1024];
        for (int i = 0; i < index.length; i++) {
            index[i] = i;
        }
        GenotypeTable ngt = FilterGenotypeTable.getInstance(gt, index);
        
        ExportUtils.writeGenotypeHDF5(ngt, hapMapOutputS);
    }
    
    public void migrateHDF5ToTASSEL5 () {
        String hapMapInputS = "M:/production/GBSV3/genotype/GBS27.small1024.imp.hmp.h5";
        String hapMapOutputS = "M:/pipelineTest/geneticMapping/anchor/GBS27.small1024.T5.imp.hmp.h5";
        MigrateHDF5FromT4T5.copyGenotypes(hapMapInputS, hapMapOutputS);
    }
    
    public static void main (String[] args) {
        new PipelineTest ();
    }
    
}
