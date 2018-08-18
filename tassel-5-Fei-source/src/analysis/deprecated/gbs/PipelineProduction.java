/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.gbs;

import format.Table;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import net.maizegenetics.analysis.data.MigrateHDF5FromT4T5;
import net.maizegenetics.analysis.gbs.AnnotateTOPM;
import net.maizegenetics.analysis.gbs.pana.PanASplitTBTPlugin;
import net.maizegenetics.analysis.gbs.SimpleGenotypeSBit;
import net.maizegenetics.analysis.gbs.TagAgainstAnchor;
import net.maizegenetics.analysis.gbs.TagBlockPosition;
import net.maizegenetics.analysis.gbs.TagCountToFastqPlugin;
import net.maizegenetics.analysis.gbs.pana.PanAAddPosToTagMapPlugin;
import net.maizegenetics.analysis.gbs.pana.PanABuildTagGWASMapPlugin;
import net.maizegenetics.analysis.gbs.pana.PanABuildTrainingSetPlugin;
import net.maizegenetics.analysis.gbs.pana.PanAFilteringTagMapPlugin;
import net.maizegenetics.analysis.gbs.pana.PanAMergeMappingResultPlugin;
import net.maizegenetics.analysis.gbs.pana.PanAModelTrainingPlugin;
import net.maizegenetics.analysis.gbs.pana.PanAPredictionPlugin;
import net.maizegenetics.analysis.gbs.pana.PanASamToMultiPositionTOPMPlugin;
import net.maizegenetics.analysis.gbs.pana.PanASplitTagBlockPosPlugin;
import net.maizegenetics.analysis.gbs.pana.PanATagGWASMappingPlugin;
import net.maizegenetics.analysis.gbs.pana.PanATagMapToFastaPlugin;
import net.maizegenetics.dna.map.TagGWASMap;
import net.maizegenetics.dna.map.TagGWASMapInfo;
import net.maizegenetics.dna.map.TagsOnGeneticMap;
import net.maizegenetics.dna.map.TagsOnPhysicalMapV3;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.tag.TagCounts;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import net.maizegenetics.dna.tag.TagsByTaxaByteHDF5TagGroups;
import utils.IoUtils;


/**
 *
 * @author Fei Lu
 */
public class PipelineProduction {
    
    public PipelineProduction () {
        //this.v3Pipeline() //deprecated;
        //this.geneticMappingPipeline();
        this.panAnchorV2Pipeline();
    }
    
    public void panAnchorV2Pipeline () {
        //this.mkRelaxedPanAnchorV2();
        //this.mkRigidPanAnchorV2();
        this.identiryPAVAnchor();
    }
    
    public void identiryPAVAnchor () {
        new IdentifyV2PAVAnchor();
    }
    
    public void mkRigidPanAnchorV2 () {
        String primeFileS = "M:\\production\\geneticMapping\\togm\\anchor.100k.togm.txt";
        String anotherFileS = "M:\\production\\panAnchorV2\\v1\\v1.panAnchor.rigid.togm.txt";
        String mergedFileS = "M:\\production\\panAnchorV2\\rigid\\v2.panAnchor.rigid.togm.txt";
        TagsOnGeneticMap togmPrime = new TagsOnGeneticMap(primeFileS, FilePacking.Text);
        TagsOnGeneticMap togmAnother = new TagsOnGeneticMap(anotherFileS, FilePacking.Text);
        togmPrime.mergeTOGM(togmAnother);
        togmPrime.sort();
        boolean[] ifOut = new boolean[togmPrime.getTagCount()];
        for (int i = 0; i < ifOut.length; i++) {
            if (togmPrime.getGChr(i) < 1 || togmPrime.getGChr(i) > 10) continue;
            ifOut[i] = true;
        }
        togmPrime.writeDistFile(mergedFileS, ifOut, FilePacking.Text);
    }
    
    public void mkRelaxedPanAnchorV2 () {
        String primeFileS = "M:\\production\\geneticMapping\\togm\\anchor.1000k.togm.txt";
        String anotherFileS = "M:\\production\\panAnchorV2\\v1\\v1.panAnchor.relaxed.togm.txt";
        String mergedFileS = "M:\\production\\panAnchorV2\\relaxed\\v2.panAnchor.relaxed.togm.txt";
        TagsOnGeneticMap togmPrime = new TagsOnGeneticMap(primeFileS, FilePacking.Text);
        TagsOnGeneticMap togmAnother = new TagsOnGeneticMap(anotherFileS, FilePacking.Text);
        togmPrime.mergeTOGM(togmAnother);
        togmPrime.sort();
        boolean[] ifOut = new boolean[togmPrime.getTagCount()];
        for (int i = 0; i < ifOut.length; i++) {
            if (togmPrime.getGChr(i) < 1 || togmPrime.getGChr(i) > 10) continue;
            ifOut[i] = true;
        }
        togmPrime.writeDistFile(mergedFileS, ifOut, FilePacking.Text);
    }
    
    public void v3Pipeline() {
        //this.initializeTOPMHDF5();
        //this.mkFastq();
        //this.mkFasta();
        //this.annotateTOPMHDF5WithAligner();
    }
    
    public void annotateTOPMHDF5WithAligner () {
        String inputFileS = "M:/production/v3gbs/topm/ini.topm.h5";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(inputFileS);
        AnnotateTOPM anno = new AnnotateTOPM (topm);
        String bowtie2SamFileS = "M:/production/v3gbs/alignment/tags.bowtie2.sam.gz";
        anno.annotateWithBowtie2(bowtie2SamFileS, 5);
        String bwaSamFileS = "M:/production/v3gbs/alignment/tags.bwa.sam.gz";
        //anno.annotateWithBWA(bwaSamFileS, 5);
        //String blastDirS = "M:/pipelineTest/GBSV3test/alignment/blastOut/";
        //anno.annotateWithBlastFromDir(blastDirS, 5);
        String bwaMemSamFileS = "M:/production/v3gbs/alignment/tags.bwa.sam.gz";
        //anno.annotateWithBWAMEM(bwaMemSamFileS, 5);
        String PETOPMFileS = "M:/production/pe/ptopm/PE.topm";
        //anno.annotateWithPE(PETOPMFileS, 5);
    }
    
    /**
     * Make FASTA file and do alignment using Blast, -m 8 -e 1e-10
     */
    public void mkFasta () {
        String inputFileS = "M:/production/v3gbs/tagCount/AllZeaMasterTags_c10_20120606.cnt";
        String outputFileS = "M:/production/v3gbs/alignment/tags.fa";
        TagCounts tc = new TagCounts(inputFileS, FilePacking.Byte);
        tc.toFASTA(outputFileS);
    }
    
    /**
     * Make Fastq file and do alignment using Bowtie2 and BWA
     * Alignment files should be in the alignment folder
     */
    public void mkFastq () {
        //bowtie2 -x library -q query -k 5 --very-sensitive-local -p 10 -S outfile.sam
        //bwa mem -t 10 -a library query > out.sam
        //bwa aln -t 10 library query > outfile.sai
        //bwa samse -n 4 library outfile.sai query > outfile.sam //-n 4 means 4 secondary alignments
        String inputFileS = "M:/production/v3gbs/tagCount/AllZeaMasterTags_c10_20120606.cnt";
        String outputFileS = "M:/production/v3gbs/alignment/tags.fq";       
        TagCountToFastqPlugin umithm = new TagCountToFastqPlugin();
		String arguments = "-i " + inputFileS + " -o " + outputFileS;
		String[] args = arguments.split(" ");
		umithm.setParameters(args);
		umithm.performFunction(null);
    }
    
    public void initializeTOPMHDF5 () {
        String inputFileS = "M:/production/v3gbs/tagCount/AllZeaMasterTags_c10_20120606.cnt";
        String outputFileS = "M:/production/v3gbs/topm/ini.topm.h5";
        TagCounts tc = new TagCounts(inputFileS, FilePacking.Byte);
        TagsOnPhysicalMapV3.createFile(tc, outputFileS);
    }
    
    public void geneticMappingPipeline () {
        //this.migrateHDF5ToTASSEL5();
        //this.transformAnchorHDF5();
        //this.splitTBTHDF5();
        //this.mkTagBlockPosition();
        //this.splitTagBlockPosition();
        //this.geneticMappingWholeTBT();
        //this.geneticMappingSubTBTCBSU();
        //this.geneticMappingSubTBTInTACC(args); needs to be ran from TACC
        //this.mergeMappingResult();
        //this.buildTagGWASMap();
        //this.tagMapToFasta();
        //this.samToMultiPositionTOPM();
        //this.addPosToTagMap();
        //this.buildTrainingSet();
        //this.modelTraining();
        //this.prediction();
        this.filterTagMap();
    }
    
    public void filterTagMap () {
        String tagGWASMapFileS = "M:\\production\\geneticMapping\\tagMap\\tagGWASMap.h5";
        int distanceCutoff = 1000000;
        String anchorFileS = "M:\\production\\geneticMapping\\togm\\anchor.1000k.togm.txt" ;
        String arguments = "-t " + tagGWASMapFileS +  " -a " + anchorFileS +  " -c " + String.valueOf(distanceCutoff);
        String[] args = arguments.split(" ");
        PanAFilteringTagMapPlugin p = new PanAFilteringTagMapPlugin();
        p.setParameters(args);
        p.performFunction(null);
    }
    
    public void prediction () {
        String wekaPath = "E:\\Database\\Weka-3-6\\weka.jar";
        String tagGWASMapFileS = "M:\\production\\geneticMapping\\tagMap\\tagGWASMap.h5";
        String modelFileS = "M:\\production\\geneticMapping\\training\\m5.mod";
        String boxcoxParemeterFileS = "M:\\production\\geneticMapping\\training\\boxcoxParemeter.txt";
        String arguments = "-t " + tagGWASMapFileS +  " -w " + wekaPath +  " -m " + modelFileS + " -b " + boxcoxParemeterFileS;
        String[] args = arguments.split(" ");
        PanAPredictionPlugin p = new PanAPredictionPlugin();
        p.setParameters(args);
        p.performFunction(null);
    }
    
    public void modelTraining () {
        String trainingSetFileS = "M:\\production\\geneticMapping\\training\\uniqueRefTrain.arff";
        String wekaPath = "E:\\Database\\Weka-3-6\\weka.jar";
        String modelFileS = "M:\\production\\geneticMapping\\training\\m5.mod";
        String trainingReportDirS = "M:\\production\\geneticMapping\\training\\report\\";
        String arguments = "-t " + trainingSetFileS +  " -w " + wekaPath +  " -m " + modelFileS + " -r " + trainingReportDirS;
        String[] args = arguments.split(" ");
        PanAModelTrainingPlugin p = new PanAModelTrainingPlugin();
        p.setParameters(args);
        p.performFunction(null);
    }
    
    public void buildTrainingSet () {
        String tagGWASMapFileS = "M:\\production\\geneticMapping\\tagMap\\tagGWASMap.h5";
        String trainingSetFileS = "M:\\production\\geneticMapping\\training\\uniqueRefTrain.arff";
        String rScriptPath = "C:\\Users\\fl262\\Documents\\R\\R-3.0.2\\bin\\Rscript";
        String boxcoxParemeterFileS = "M:\\production\\geneticMapping\\training\\boxcoxParemeter.txt";
        int maxInstance = 30000;
        String arguments = "-m " + tagGWASMapFileS +  " -t " + trainingSetFileS + " -r " + rScriptPath +  " -b " + boxcoxParemeterFileS + " -i " + String.valueOf(maxInstance);
        String[] args = arguments.split(" ");
        PanABuildTrainingSetPlugin p = new PanABuildTrainingSetPlugin();
        p.setParameters(args);
        p.performFunction(null);
    }
    
    public void addPosToTagMap () {
        String tagGWASMapFileS = "M:\\production\\geneticMapping\\tagMap\\tagGWASMap.h5";
        String topmV3FileS = "M:\\production\\geneticMapping\\alignment\\tagGWASMap.v3.topm.h5";
        String arguments = "-i " + tagGWASMapFileS +  " -t " + topmV3FileS;
        String[] args = arguments.split(" ");
        PanAAddPosToTagMapPlugin  p = new PanAAddPosToTagMapPlugin();
        p.setParameters(args);
        p.performFunction(null);
    }
    
    public void samToMultiPositionTOPM () {
        String samFileS = "M:\\production\\geneticMapping\\alignment\\tagGWASMap.sam.gz";
        String tagGWASMapFileS = "M:\\production\\geneticMapping\\tagMap\\tagGWASMap.h5";
        String topmV3FileS = "M:\\production\\geneticMapping\\alignment\\tagGWASMap.v3.topm.h5";
        String arguments = "-i " + samFileS + " -t " + tagGWASMapFileS +  " -o " + topmV3FileS;
        String[] args = arguments.split(" ");
        PanASamToMultiPositionTOPMPlugin stt = new PanASamToMultiPositionTOPMPlugin();
        stt.setParameters(args);
        stt.performFunction(null);
    }
    
    public void tagMapToFasta () {
        String tagGWASMapFileS = "M:\\production\\geneticMapping\\tagMap\\tagGWASMap.h5";
        String fastaFileS = "M:\\production\\geneticMapping\\alignment\\tagGWASMap.fa";
        String arguments = "-i " + tagGWASMapFileS +  " -o " + fastaFileS;
        String[] args = arguments.split(" ");
        PanATagMapToFastaPlugin tmtf = new PanATagMapToFastaPlugin();
        tmtf.setParameters(args);
        tmtf.performFunction(null);
    }
    
    public void buildTagGWASMap () {
        String mappingResultFileS = "M:\\production\\geneticMapping\\gwasResult\\pivotTBT.gwas.txt";
        String tagCountFileS = "M:/production/v3gbs/tagCount/AllZeaMasterTags_c10_20120606.cnt";
        String tagGWASMapFileS = "M:\\production\\geneticMapping\\tagMap\\tagGWASMap.h5";
        String arguments = "-i " + mappingResultFileS + " -t " + tagCountFileS + " -o " + tagGWASMapFileS;
        String[] args = arguments.split(" ");
        PanABuildTagGWASMapPlugin mrtg = new PanABuildTagGWASMapPlugin();
        mrtg.setParameters(args);
        mrtg.performFunction(null);
    }
    
    public void mergeMappingResult () {
        String subResultDirS = "M:\\production\\geneticMapping\\gwasResult\\subResult\\";
        String mergedResultFileS = "M:\\production\\geneticMapping\\gwasResult\\pivotTBT.gwas.txt";
        String arguments = "-i " + subResultDirS + " -o " + mergedResultFileS;
        String[] args = arguments.split(" ");
        PanAMergeMappingResultPlugin mmr = new PanAMergeMappingResultPlugin();
        mmr.setParameters(args);
        mmr.performFunction(null);
    }
    
    public void geneticMappingSubTBTInTACC (String[] args) {
        PanATagGWASMappingPlugin tgm = new  PanATagGWASMappingPlugin();
        tgm.setParameters(args);
        tgm.performFunction(null);
    }
    
    public void geneticMappingSubTBTCBSU () {
        String sBitGenotypeFileS = "/workdir/mingh/genotype/AllZeaGBSv27i3b.sBit.h5";
        String outDirS = "/workdir/mingh/result/";
        String tbtDirS = "/workdir/mingh/subTBT/";
        String tbpDirS = "/workdir/mingh/subTBP/";
        File[] tbts = new File (tbtDirS).listFiles();
        Arrays.sort(tbts);
        File[] tbps = new File (tbpDirS).listFiles();
        Arrays.sort(tbps);
        for (int i = 0; i < tbts.length; i++) {
            String arguments = "-g " + sBitGenotypeFileS + " -t " + tbts[i].getAbsolutePath() + " -b " + tbps[i].getAbsolutePath() + " -o " + outDirS + " -c max -s 65536 -cs 0 -ce 1";
            
            String[] args = arguments.split(" ");
            PanATagGWASMappingPlugin tgm = new  PanATagGWASMappingPlugin();
            tgm.setParameters(args);
            tgm.performFunction(null);
        }
    }
    
    public void geneticMappingWholeTBT () {
        //String hapMapHDF5 = "M:/production/geneticMapping/genotype/AllZeaGBSv27i3b.sBit.h5";
        //String tbtHDF5 = "M:/production/geneticMapping/tbt/mergeTBTHDF5_mergedtaxa_pivot_20120921.h5";
        //String blockFileS = "M:/production/geneticMapping/tagBlock/GBS27TagBlock.tbp";
        //String outfileS = "M:/production/geneticMapping/result/outfile.txt";
        
        String hapMapHDF5 = "/workdir/mingh/AllZeaGBSv27i3b.sBit.h5";
        String tbtHDF5 = "/workdir/mingh/mergeTBTHDF5_mergedtaxa_pivot_20120921.h5";
        String blockFileS = "/workdir/mingh/GBS27TagBlock.tbp";
        String outfileS = "/workdir/mingh/outfile.txt";

        TagAgainstAnchor taa = new TagAgainstAnchor(hapMapHDF5, tbtHDF5, blockFileS, outfileS, 0.000001, 20, -1, 65536, 0, 10);
        //TagAgainstAnchor.getChunkNum(tbtHDF5, 65536);
    }
    
    public void splitTagBlockPosition () {
        String blockFileS = "M:\\production\\geneticMapping\\tagBlock\\GBS27TagBlock.tbp";
        String outputDirS = "M:\\production\\geneticMapping\\tagBlock\\subTBP\\";
        String chunkSize = "65536";
        String arguments = "-i " + blockFileS + " -s " + chunkSize + " -o " + outputDirS;
        String[] args = arguments.split(" ");
        PanASplitTagBlockPosPlugin stb = new PanASplitTagBlockPosPlugin();
        stb.setParameters(args);
        stb.performFunction(null);
    }
    
    public void mkTagBlockPosition () {
        String tbtHDF5 = "M:/production/geneticMapping/tbt/mergeTBTHDF5_mergedtaxa_pivot_20120921.h5";
        String topmFileS = "M:\\production\\geneticMapping\\tagBlock\\AllZeaGBSv2.6ProdTOPM_20130605.topm.h5";
        String blockFileS = "M:/production/geneticMapping/tagBlock/GBS27TagBlock.tbp";
        int topmVersionValue = 1;
        TagBlockPosition tbp = new TagBlockPosition(tbtHDF5, topmFileS, topmVersionValue);
        tbp.writeTagBlockPosition(blockFileS);
    }
    
    public void splitTBTHDF5 () {
        String inputTBTS = "M:/production/geneticMapping/tbt/mergeTBTHDF5_mergedtaxa_pivot_20120921.h5";
        String outputDirS = "M:\\production\\geneticMapping\\tbt\\subTBT\\";
        String chunkSize = "65536";
        String arguments = "-i " + inputTBTS + " -s " + chunkSize + " -o " + outputDirS;
        String[] args = arguments.split(" ");
        PanASplitTBTPlugin sbp = new PanASplitTBTPlugin();
        sbp.setParameters(args);
        sbp.performFunction(null);
    }
    
    public void transformAnchorHDF5 () {
        //on PC
        //String hapMapInputFileS = "M:/production/geneticMapping/genotype/AllZeaGBSv27i3b.T5.imp.hmp.h5";
        //String sBitFileS = "M:/production/geneticMapping/genotype/AllZeaGBSv27i3b.sBit.h5";
        
        //on Server
        String hapMapInputFileS = "/workdir/mingh/AllZeaGBSv27i3b.T5.imp.hmp.h5";
        String sBitFileS = "/workdir/mingh/geneticMapping/genotype/AllZeaGBSv27i3b.sBit.h5";
        new SimpleGenotypeSBit(hapMapInputFileS, sBitFileS);
        
    }
    
    public void migrateHDF5ToTASSEL5 () {
        String hapMapInputS = "M:/production/geneticMapping/genotype/AllZeaGBSv27i3b.imp.hmp.h5";
        String hapMapOutputS = "M:/production/geneticMapping/genotype/AllZeaGBSv27i3b.T5.imp.hmp.h5";
        MigrateHDF5FromT4T5.copyGenotypes(hapMapInputS, hapMapOutputS);
    }
    
    public static void main (String[] args) {
        //new PipelineProduction(args);
        new PipelineProduction();
    }
}
