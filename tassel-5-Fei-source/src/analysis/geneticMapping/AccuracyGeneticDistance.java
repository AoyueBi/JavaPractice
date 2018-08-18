/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package analysis.geneticMapping;

import format.Table;
import graphcis.r.DensityPlotMultiClass;
import graphcis.r.ScatterPlot;
import graphcis.r.ScatterPlotMultiClass;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.analysis.gbs.TagBlockPosition;
import net.maizegenetics.analysis.gbs.pana.PanAAddPosToTagMapPlugin;
import net.maizegenetics.analysis.gbs.pana.PanABuildTagGWASMapPlugin;
import net.maizegenetics.analysis.gbs.pana.PanABuildTrainingSetPlugin;
import net.maizegenetics.analysis.gbs.pana.PanAH5ToAnchorPlugin;
import net.maizegenetics.analysis.gbs.pana.PanAModelTrainingPlugin;
import net.maizegenetics.analysis.gbs.pana.PanAPredictionPlugin;
import net.maizegenetics.analysis.gbs.pana.PanASamToMultiPositionTOPMPlugin;
import net.maizegenetics.analysis.gbs.pana.PanATagGWASMappingPlugin;
import net.maizegenetics.analysis.gbs.pana.PanATagMapToFastaPlugin;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.map.TagGWASMap;
import net.maizegenetics.dna.map.TagGWASMapInfo;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import net.maizegenetics.dna.tag.TagsByTaxaByte;
import net.maizegenetics.dna.tag.TagsByTaxaByteHDF5TagGroups;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import org.apache.commons.math3.random.RandomDataImpl;
import utils.IoUtils;
import utils.stats.r.PCA;

/**
 *
 * @author Fei Lu <fl262@cornell.edu>
 */
public class AccuracyGeneticDistance {
    
    public AccuracyGeneticDistance () {
        this.getPopulationNames();
        //this.selectGenotype();
        //this.compareGeneticDistance();
        //this.mkTBT();
        //this.convertTBT();
        //this.mkTagBlock();
        //this.convertSBitGenotype();
        //this.mapNAM();
        //this.mapAmes();
        //this.buildTagGWASMap();
        //this.tagMapToFasta();
        //this.samToMultiPositionTOPM();
        //this.addPosToTagMapPlugin();
        //this.buildTrainingSet();
        //this.modelTraining();
        //this.predict();
        //this.compareResolution();
    }
    
    public void compareResolution () {
        String namMap = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\tagMap\\namSelect.tagGWASMap.h5";
        String amesMap = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\tagMap\\amesSelect.tagGWASMap.h5";
        String resolutionComparisonFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\resolutionDensity.pdf";
        String resolutionScatterFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\resolutionScatter.pdf";
        TagGWASMap nam = new TagGWASMap (namMap);
        TagGWASMap ames = new TagGWASMap (amesMap);
        int size = 10000;
        double[] namDis = new double[size];
        double[] amesDis = new double[size];
        int cnt = 0;
        for (int i = 0; i < nam.getTagCount(); i++) {
            TagGWASMapInfo info = nam.getTagGWASMapInfo(i);
            if (!info.isUniqueRef()) continue;
            double d = 0;
            if (info.gChr == info.pChr) {
                d = Math.log10(Math.abs(info.gPos-info.pPos));
            }
            else {
                d = Math.log10(Math.abs(info.gPos-info.pPos) + 1000000000);
            }
            namDis[cnt] = d;
            cnt++;
            if (cnt == size) break;
        }
        cnt = 0;
        for (int i = 0; i < ames.getTagCount(); i++) {
            TagGWASMapInfo info = ames.getTagGWASMapInfo(i);
            if (!info.isUniqueRef()) continue;
            double d = 0;
            if (info.gChr == info.pChr) {
                d = Math.log10(Math.abs(info.gPos-info.pPos));
            }
            else {
                d = Math.log10(Math.abs(info.gPos-info.pPos) + 1000000000);
            }
            amesDis[cnt] = d;
            cnt++;
            if (cnt == size) break;
        }
        double[][] dis = new double[2][];
        String[] names = {"NAM", "Ames"};
        dis[0] = namDis;
        dis[1] = amesDis;
        DensityPlotMultiClass d = new DensityPlotMultiClass(dis, names);
        d.setXLab("Log10 (distance)");
        d.setYLab("Density");
        d.setTitle("GWAS mapping accuracy of unique B73 tags");
        d.saveGraph(resolutionComparisonFileS);
        cnt = 0;
        
        for (int i = 0; i < ames.getTagCount(); i++) {
            long[] t = ames.getTag(i);
            int index = nam.getTagIndex(t);
            if (index < 0) continue;
            TagGWASMapInfo info = ames.getTagGWASMapInfo(i);
            if (!info.isUniqueRef()) continue;
            double dAmes = 0;
            if (info.gChr == info.pChr) {
                dAmes = Math.log10(Math.abs(info.gPos-info.pPos));
            }
            else {
                dAmes = Math.log10(Math.abs(info.gPos-info.pPos) + 1000000000);
            }
            amesDis[cnt] = dAmes;
            info = nam.getTagGWASMapInfo(index);
            double dNam = 0;
            if (info.gChr == info.pChr) {
                dNam = Math.log10(Math.abs(info.gPos-info.pPos));
            }
            else {
                dNam = Math.log10(Math.abs(info.gPos-info.pPos) + 1000000000);
            }
            namDis[cnt] = dNam;
            cnt++;
            if (cnt == size) break;
        }
        ScatterPlot s = new ScatterPlot(namDis, amesDis);
        s.setColor(255, 0, 0, 60);
        s.setPlottingCharacter(16);
        s.setXLab("Log10 (distance in NAM)");
        s.setYLab("Log10 (distance in Ames)");
        s.setTitle("GWAS mapping accuracy in two populations");
        s.setSlideMode();
        //s.showGraph();
        s.saveGraph(resolutionScatterFileS);
    }
    
    public void predict () {
        String wekaPath = "E:\\Database\\Weka-3-6\\weka.jar";
        String tagGWASMapFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\tagMap\\amesSelect.tagGWASMap.h5";
        String modelFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\training\\amesSelect_m5.mod";
        String boxcoxParemeterFileS = "E:/Research/geneticMapping/accuracyGeneticDistance/mapping/training/ames_boxcoxParemeter.txt";
        String arguments = "-t " + tagGWASMapFileS +  " -w " + wekaPath +  " -m " + modelFileS + " -b " + boxcoxParemeterFileS;
        String[] args = arguments.split(" ");
        PanAPredictionPlugin p = new PanAPredictionPlugin();
        p.setParameters(args);
        p.performFunction(null);
    }
     
    public void modelTraining () {
        String trainingSetFileS = "E:/Research/geneticMapping/accuracyGeneticDistance/mapping/training/namSelect_uniqueRefTrain.arff";
        String wekaPath = "E:\\Database\\Weka-3-6\\weka.jar";
        String modelFileS = "E:/Research/geneticMapping/accuracyGeneticDistance/mapping/training/namSelect_m5.mod";
        String trainingReportDirS = "E:/Research/geneticMapping/accuracyGeneticDistance/mapping/training//namReport/";
        String arguments = "-t " + trainingSetFileS +  " -w " + wekaPath +  " -m " + modelFileS + " -r " + trainingReportDirS;
        String[] args = arguments.split(" ");
        PanAModelTrainingPlugin p = new PanAModelTrainingPlugin();
        p.setParameters(args);
        p.performFunction(null);
    }
    
    public void buildTrainingSet () {
        String tagGWASMapFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\tagMap\\amesSelect.tagGWASMap.h5";
        String trainingSetFileS = "E:/Research/geneticMapping/accuracyGeneticDistance/mapping/training/amesSelect_uniqueRefTrain.arff";
        String rScriptPath = "C:\\Users\\fl262\\Documents\\R\\R-3.0.2\\bin\\Rscript";
        String boxcoxParemeterFileS = "E:/Research/geneticMapping/accuracyGeneticDistance/mapping/training/ames_boxcoxParemeter.txt";
        int maxInstance = 30000;
        String arguments = "-m " + tagGWASMapFileS +  " -t " + trainingSetFileS + " -r " + rScriptPath +  " -b " + boxcoxParemeterFileS + " -i " + String.valueOf(maxInstance);
        String[] args = arguments.split(" ");
        PanABuildTrainingSetPlugin p = new PanABuildTrainingSetPlugin();
        p.setParameters(args);
        p.performFunction(null);
    }
    
    public void addPosToTagMapPlugin () {
        String tagGWASMapFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\tagMap\\amesSelect.tagGWASMap.h5";
        String topmV3FileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\alignment\\amesSelect_tagGWASMap.v3.topm.h5";
        String arguments = "-i " + tagGWASMapFileS +  " -t " + topmV3FileS;
        String[] args = arguments.split(" ");
        PanAAddPosToTagMapPlugin  p = new PanAAddPosToTagMapPlugin();
        p.setParameters(args);
        p.performFunction(null);
    }
    
    public void samToMultiPositionTOPM () {
        String samFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\alignment\\namSelect_tagGWASMap.sam";
        String tagGWASMapFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\tagMap\\namSelect.tagGWASMap.h5";
        String topmV3FileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\alignment\\namSelect_tagGWASMap.v3.topm.h5";
        String arguments = "-i " + samFileS + " -t " + tagGWASMapFileS +  " -o " + topmV3FileS;
        String[] args = arguments.split(" ");
        PanASamToMultiPositionTOPMPlugin stt = new PanASamToMultiPositionTOPMPlugin();
        stt.setParameters(args);
        stt.performFunction(null);
    }
    
    public void alignWithBowtie2() {
        //readlized some UB73 is not unique, due to different version of aligner. parameters changed
    }
    
    public void tagMapToFasta () {
        String tagGWASMapFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\tagMap\\namSelect.tagGWASMap.h5";
        String fastaFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\alignment\\namSelect_tagGWASMap.fa";
        String arguments = "-i " + tagGWASMapFileS +  " -o " + fastaFileS;
        String[] args = arguments.split(" ");
        PanATagMapToFastaPlugin tmtf = new PanATagMapToFastaPlugin();
        tmtf.setParameters(args);
        tmtf.performFunction(null);
    }
    
    public void buildTagGWASMap () {
        //String mappingResultFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\mappingResult\\amesSelect.pivotTBT.gwas.txt";
        String mappingResultFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\mappingResult\\namSelect.pivotTBT.gwas.txt";
        String tagCountFileS = "M:\\pav\\tagCount\\merged_20.cnt";
        String tagGWASMapFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\tagMap\\amesSelect.tagGWASMap.h5";
        //String tagGWASMapFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\tagMap\\namSelect.tagGWASMap.h5";
        String arguments = "-i " + mappingResultFileS + " -t " + tagCountFileS + " -o " + tagGWASMapFileS;
        String[] args = arguments.split(" ");
        PanABuildTagGWASMapPlugin mrtg = new PanABuildTagGWASMapPlugin();
        mrtg.setParameters(args);
        mrtg.performFunction(null);
    }
    
    public void mapAmes () {
        String sBitGenotypeFileS = "/SSD/mingh/amesSelect.sBit.h5";
        String outDirS = "/SSD/mingh/sub/";
        String tbtDirS = "/SSD/mingh/subTBT/";
        String tbpDirS = "/SSD/mingh/subTBP/";
        File[] tbts = new File (tbtDirS).listFiles();
        File[] tbps = new File (tbpDirS).listFiles();
        for (int i = 0; i < tbts.length; i++) {
            String arguments = "-g " + sBitGenotypeFileS + " -t " + tbts[i].getAbsolutePath() + " -b " + tbps[i].getAbsolutePath() + " -o " + outDirS + " -c max -s 1000 -cs 0 -ce 1107";
            
            String[] args = arguments.split(" ");
            PanATagGWASMappingPlugin tgm = new  PanATagGWASMappingPlugin();
            tgm.setParameters(args);
            tgm.performFunction(null);
        }
    }
    
    public void mapNAM () {
        String sBitGenotypeFileS = "/SSD/mingh/namSelect.sBit.h5";
        String outDirS = "/SSD/mingh/sub/";
        String tbtDirS = "/SSD/mingh/subTBT/";
        String tbpDirS = "/SSD/mingh/subTBP/";
        File[] tbts = new File (tbtDirS).listFiles();
        File[] tbps = new File (tbpDirS).listFiles();
        for (int i = 0; i < tbts.length; i++) {
            String arguments = "-g " + sBitGenotypeFileS + " -t " + tbts[i].getAbsolutePath() + " -b " + tbps[i].getAbsolutePath() + " -o " + outDirS + " -c max -s 1000 -cs 0 -ce 1107";
            
            String[] args = arguments.split(" ");
            PanATagGWASMappingPlugin tgm = new  PanATagGWASMappingPlugin();
            tgm.setParameters(args);
            tgm.performFunction(null);
        }
    }
    
    public void convertSBitGenotype () {
        String h5GentoypeFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\genotype\\namSelect.hmp.h5";
        String sBitGenotypeFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\genotype\\namSelect.sBit.h5";
        String arguments = "-i " + h5GentoypeFileS +  " -o " + sBitGenotypeFileS;
        String[] args = arguments.split(" ");
        PanAH5ToAnchorPlugin hta = new PanAH5ToAnchorPlugin();
        hta.setParameters(args);
        hta.performFunction(null);
    }
    
    public void mkTagBlock () {
        String topmFileS = "N:\\Zea\\build20120110\\topm\\allZea_20120115.topm";
        String namH5TBT = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\tbt\\namSelect.pivotTBT.h5";
        String tbpFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\tagBlockPosition\\UB73.tbp";
        TagBlockPosition tbp = new TagBlockPosition(namH5TBT, topmFileS, 0);
        tbp.writeTagBlockPosition(tbpFileS);
    }
    
    public void convertTBT () {
        String namTBTByte = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\tbt\\namSelect.tbt.byte";
        String amesTBTByte = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\tbt\\amesSelect.tbt.byte";
        String namH5TBT = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\tbt\\namSelect.pivotTBT.h5";
        String amesH5TBT = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\tbt\\amesSelect.pivotTBT.h5";
        new TagsByTaxaByteHDF5TagGroups(new TagsByTaxaByte(namTBTByte, FilePacking.Byte), namH5TBT);
        new TagsByTaxaByteHDF5TagGroups(new TagsByTaxaByte(amesTBTByte, FilePacking.Byte), amesH5TBT);
    }
    
    public void mkTBT () {
        String tbtSourceFileS = "M:\\pav\\mergedTBT\\UB73.tbt.byte";
        String namTaxaFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\namSelect.txt";
        String amesTaxaFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\amesSelect.txt";
        String namTBTByte = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\tbt\\namSelect.tbt.byte";
        String amesTBTByte = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\mapping\\tbt\\amesSelect.tbt.byte";
        Table t = new Table (namTaxaFileS);
        String[] namTaxa = new String[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) namTaxa[i] = t.content[i][0];
        Arrays.sort(namTaxa);
        t = new Table(amesTaxaFileS);
        String[] amesTaxa = new String[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) amesTaxa[i] = t.content[i][0];
        Arrays.sort(amesTaxa);
        TagsByTaxaByte tbtSource = new TagsByTaxaByte(tbtSourceFileS, FilePacking.Byte);
        this.selectSubsetTBT(tbtSource, namTaxa, namTBTByte);
        this.selectSubsetTBT(tbtSource, amesTaxa, amesTBTByte);
    }
    
    private void selectSubsetTBT (TagsByTaxaByte tbt, String[] selectTaxa, String selectTBTFileS) {
        int[] index = new int[selectTaxa.length];
        String[] sourceTaxa = tbt.getTaxaNames();
        for (int i = 0; i < index.length; i++) {
            for (int j = 0; j < sourceTaxa.length; j++) {
                if (selectTaxa[i].equals(sourceTaxa[j])) {
                    index[i] = j;
                    break;
                }
            }
        }
        try {
            DataOutputStream dos = IoUtils.getBinaryWriter(selectTBTFileS);
            dos.writeInt(tbt.getTagCount());
            dos.writeInt(tbt.getTagSizeInLong());
            dos.writeInt(selectTaxa.length);
            for (int i = 0; i < selectTaxa.length; i++) {
                dos.writeUTF(selectTaxa[i]);
            }
            for (int i = 0; i < tbt.getTagCount(); i++) {
                long[] t = tbt.getTag(i);
                for (int j = 0; j < tbt.getTagSizeInLong(); j++) {
                    dos.writeLong(t[j]);
                }
                dos.writeByte((byte)tbt.getTagLength(i));
                for (int j = 0; j < selectTaxa.length; j++) {
                    dos.writeByte(tbt.getReadCountForTagTaxon(i, index[j]));
                 
                }
                if (i%100000 == 0) System.out.println(String.valueOf(i)+ " tags written");
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void compareGeneticDistance () {
        String combineGenotype = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\combine.hmp.h5";
        String pcaPlotFileS = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\pca.pdf";
        String pcaSummary = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\pca.summary.txt";
        GenotypeTable gt = ImportUtils.readGuessFormat(combineGenotype);
        int siteSize = gt.numberOfTaxa();
        int[] siteIndex = new int[siteSize];
        for (int i = 0; i < siteSize; i++) {
            siteIndex[i] = (int)(gt.numberOfSites()*Math.random());
        }
        double[][] matrix  = new double[siteSize][gt.numberOfTaxa()];
        for (int i = 0; i < siteSize; i++) {
            byte major = gt.majorAllele(i);
            byte minor = gt.minorAllele(i);
            int het = major+minor;
            if (major > minor) het = minor+major;
            int ma = major*2;
            int mi = minor*2;
            for (int j = 0; j < gt.numberOfTaxa(); j++) {
                byte[] genotypes = gt.genotypeArray(j, i);
                int sum = genotypes[0]+genotypes[1];
                if (sum == ma) {
                    matrix[i][j] = 0;
                }
                else if (sum == mi) {
                    matrix[i][j] = 2;
                }
                else if (sum == het) {
                    matrix[i][j] = 1;
                }
                else {
                    matrix[i][j] = 1;
                }
            }
        }
        PCA p = new PCA(matrix);
        String[] classes = new String[gt.numberOfTaxa()];
        for (int i = 0; i < classes.length; i++) {
            if (gt.taxa().get(i).getName().startsWith("Z00")) {
                classes[i] = "NAM";
            }
            else {
                classes[i] = "Ames";
            }
        }
        ScatterPlotMultiClass s = new ScatterPlotMultiClass(p.getScoresOfPC(0), p.getScoresOfPC(1), classes);
        s.setXLab("PC 1");
        s.setYLab("PC 2");
        s.setTitle("PCA of inbreds in NAM and Ames populations");
        s.saveGraph(pcaPlotFileS);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(pcaSummary);
            bw.write(p.getSummary());
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void selectGenotype () {
        String genotypeFileS = "/SSD/mingh/impAll.hmp.h5";
        String namSelect = "/SSD/mingh/namSelect.txt";
        String amesSelect = "/SSD/mingh/amesSelect.txt";
        String namGenotype = "/SSD/mingh/namSelect.hmp.h5";
        String amesGenotype = "/SSD/mingh/amesSelect.hmp.h5";
        String combineGenotype = "/SSD/mingh/combine.hmp.h5";
        
        Table t = new Table(namSelect);
        TaxaListBuilder combineBuilder=new TaxaListBuilder();
        TaxaListBuilder namBuilder=new TaxaListBuilder();
        for (int i = 0; i < t.getRowNumber(); i++) {
            Taxon ta = new Taxon(t.content[i][0]);
            namBuilder.add(ta);
            combineBuilder.add(ta);
        }
        t = new Table(amesSelect);
        TaxaListBuilder amesBuilder = new TaxaListBuilder();
        for (int i = 0; i < t.getRowNumber(); i++) {
            Taxon ta = new Taxon(t.content[i][0]);
            amesBuilder.add(ta);
            combineBuilder.add(ta);
        }
        GenotypeTable gt = ImportUtils.readGuessFormat(genotypeFileS);
        GenotypeTable namGt = FilterGenotypeTable.getInstance(gt, namBuilder.build());
        GenotypeTable amesGt = FilterGenotypeTable.getInstance(gt, amesBuilder.build());
        GenotypeTable combineGt = FilterGenotypeTable.getInstance(gt, combineBuilder.build());
        ExportUtils.writeGenotypeHDF5(namGt, namGenotype);
        ExportUtils.writeGenotypeHDF5(amesGt, amesGenotype);
        ExportUtils.writeGenotypeHDF5(combineGt, combineGenotype);
    }
    
    public void getPopulationNames () {
        String allFullName = "M:\\pav\\cnv2\\taxaName\\Jan.taxa.fullname.txt";
        String namSelect = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\namSelect.txt";
        String amesFullName = "M:\\pav\\cnv2\\taxaName\\282Ames.txt";
        String amesSelect = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\amesSelect_random.txt";
        
        Table t = new Table(allFullName);
        ArrayList<String> l = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.content[i][0].startsWith("Z001") || t.content[i][0].startsWith("Z023")) {
                l.add(t.content[i][0]);
            }
        }
        String[] nam = l.toArray(new String[l.size()]);
        Arrays.sort(nam);
        String[] ames = new String[nam.length];
        t = new Table(amesFullName);
        String[] allAmes = new String[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            allAmes[i] = t.content[i][0];
        }
        int[] index = new RandomDataImpl().nextPermutation(allAmes.length-1, ames.length);
        for (int i = 0; i < index.length; i++) {
            ames[i] = allAmes[index[i]];
        }
        Arrays.sort(ames);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(namSelect);
            bw.write("Taxa");
            bw.newLine();
            for (int i = 0; i < nam.length; i++) {
                bw.write(nam[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
            bw = IoUtils.getTextWriter(amesSelect);
            bw.write("Taxa");
            bw.newLine();
            for (int i = 0; i < index.length; i++) {
                bw.write(ames[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
