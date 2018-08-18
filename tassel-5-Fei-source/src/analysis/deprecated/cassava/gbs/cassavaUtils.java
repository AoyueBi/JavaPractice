/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.cassava.gbs;

import format.Bins;
import format.Table;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import net.maizegenetics.dna.tag.TagsByTaxaByte;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import utils.IoUtils;


/**
 *
 * @author fl262
 */
public class cassavaUtils {
    
    public cassavaUtils () {
        
    }
    
    public void mkBinSummary (String statisticFileS, String binFileS) {
        Table t = new Table (statisticFileS);
        int[] cor = new int[t.getRowNumber()];
        double[] maf = new double[t.getRowNumber()];
        double[] missing = new double[t.getRowNumber()];
        double[] het = new double[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            cor[i] = Integer.valueOf(t.content[i][3]);
            maf[i] = Double.valueOf(t.content[i][4]);
            missing[i] = Double.valueOf(t.content[i][5]);
            het[i] = Double.valueOf(t.content[i][6]);
        }
        Bins mafBin = new Bins(1, cor[cor.length-1],500000, cor, maf);
        Bins missingBin = new Bins(1, cor[cor.length-1],500000, cor, missing);
        Bins hetBin = new Bins(1, cor[cor.length-1],500000, cor, het);
        BufferedWriter bw = IoUtils.getTextWriter(binFileS);
        try {
            bw.write("BinID\tBinStart\tMeanMAF\tMeanMissing\tMeanHet\tSDMAF\tSDMissing\tSDHet");
            bw.newLine();
            for (int i = 0; i < mafBin.getBinNum(); i++) {
                bw.write(String.valueOf(i+1)+"\t"+String.valueOf(mafBin.getBinEnd(i))+"\t"+String.valueOf(mafBin.getBinMedian(i))+"\t"+String.valueOf(missingBin.getBinMedian(i))+"\t"+String.valueOf(hetBin.getBinMedian(i)));
                bw.write("\t"+String.valueOf(mafBin.getBinSD(i))+"\t"+String.valueOf(missingBin.getBinSD(i))+"\t"+String.valueOf(hetBin.getBinSD(i)));
                bw.newLine();
            }
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void splictH5ToChromosomes (String inputGenotypeFileS, String outputDirS) {
        GenotypeTable gt = ImportUtils.readGuessFormat(inputGenotypeFileS);
        Chromosome[] chrs =gt.chromosomes();
        for (int i = 0; i < chrs.length; i++) {
            File outfile = new File(outputDirS, "genotypeFilterByHets.chr"+String.valueOf(chrs[i].getChromosomeNumber())+".h5");
            //FilterGenotypeTable fgt = FilterGenotypeTable.getInstance(gt, chrs[i]);
            FilterGenotypeTable fgt = null;
            //ExportUtils.writeGenotypeHDF5(fgt, outfile.getAbsolutePath());
            //System.out.println(String.valueOf(fgt.numberOfSites())+" sites, " + String.valueOf(fgt.numberOfTaxa()) + " taxa, on chromosome " + chrs[i].getName() + " is written to " + outfile.getAbsolutePath());
        }
    }
    
    public void splitHapMapToChromosomes (String inputHapMapFileS, String outputHapMapDirS, int chrNum) {
        int[] chrID = new int[chrNum];
        File[] outfiles = new File[chrNum];
        for (int i = 0; i < chrNum; i++) {
            chrID[i] = i+1;
            outfiles[i] = new File(outputHapMapDirS, "genotypeFilterByHets.chr"+String.valueOf(chrID[i])+".hmp.txt");
        }
        try {
           BufferedWriter[] bws = new BufferedWriter[chrNum];
           for (int i = 0; i < bws.length; i++) {
               bws[i] = new BufferedWriter(new FileWriter(outfiles[i]), 65536);
           }
           BufferedReader br = new BufferedReader (new FileReader(inputHapMapFileS),65536);
           String temp;
           ArrayList<String> sampleList = new ArrayList();
           while (!(temp = br.readLine()).startsWith("rs")) {
               sampleList.add(temp);
           }
           String header = temp;
           for (int i = 0; i < bws.length; i++) {
               for (int j = 0; j < sampleList.size(); j++) {
                   bws[i].write(sampleList.get(i));
                   bws[i].newLine();
               }
               bws[i].write(header);
               bws[i].newLine();
           }
           String[] tem;
           int hit;
           while ((temp = br.readLine()) != null ) {
               tem = temp.substring(0, 50).split("\t");
               hit = Arrays.binarySearch(chrID, Integer.valueOf(tem[2]));
               bws[hit].write(temp);
               bws[hit].newLine();
           }
           for (int i = 0; i < bws.length; i++) {
               bws[i].flush();
               bws[i].close();
           }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    
    public void filterSiteEndChr1 () {
        String inputGenotypeFileS = "M:\\pipelineTest\\cassava\\gbs\\original_genotype\\chr1.hmp.txt";
        String filterGenotypeFileS = "M:\\pipelineTest\\cassava\\gbs\\filtered_genotype\\chr1_filter.hmp.txt";
        String endChr1GenotypeFileS = "M:\\pipelineTest\\cassava\\gbs\\filtered_genotype\\endChr1.hmp.txt";
        String filterEndChr1GenotypeFileS = "M:\\pipelineTest\\cassava\\gbs\\filtered_genotype\\endChr1_filter.hmp.txt";
        String filterEndChr1GenotypeHomoTaxaFileS = "M:\\pipelineTest\\cassava\\gbs\\filtered_genotype\\endChr1_filter_homoTaxa.hmp.txt";
        int startPos = 19000000;
        int endPos = 20420071;
        //GenotypeTable gt = ImportUtils.readFromHapmap(inputGenotypeFileS);
        //FilterGenotypeTable fgt = FilterGenotypeTable.getInstance(gt, gt.chromosomes()[0], startPos, endPos);
        //ExportUtils.writeToHapmap(fgt, endChr1GenotypeFileS);
        
        GenotypeTable agt = ImportUtils.readFromHapmap(filterGenotypeFileS);
        GenotypeTable afgt = FilterGenotypeTable.getInstance(agt, agt.chromosomes()[0], startPos, endPos);
        ExportUtils.writeToHapmap(afgt, filterEndChr1GenotypeFileS);
       
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < afgt.numberOfTaxa(); i++) {
            if ((double)afgt.heterozygousCountForTaxon(i)/afgt.numberOfSites() < 0.1) tlb.add(afgt.taxa().get(i));
        }
        TaxaList tl = tlb.build();
        afgt = FilterGenotypeTable.getInstance(afgt, tl);
        ExportUtils.writeToHapmap(afgt, filterEndChr1GenotypeHomoTaxaFileS);
        
    }
    
    public void filterSite (String inputGenotype, String outputGenotype, double coverageForSite, double minMaf, double maxHetForTaxon, double maxHetForSite) {
        GenotypeTable gt = ImportUtils.readGuessFormat(inputGenotype);
        Chromosome[] theChr=gt.chromosomes();
        int window=1000;
        ArrayList<Integer> siteKeepList = new ArrayList();
        ArrayList<Integer> hiCoverSiteList;
        for (Chromosome currentChr : theChr) {
            System.out.println("Filtering chromosome "+currentChr.getName());
            hiCoverSiteList = new ArrayList();
            int[] chrStartEnd = gt.firstLastSiteOfChromosome(currentChr);
            for (int i = chrStartEnd[0]; i < chrStartEnd[1]+1; i++) {
                if ((double) gt.totalNonMissingForSite(i)/gt.numberOfTaxa() > coverageForSite && gt.minorAlleleFrequency(i) > minMaf) hiCoverSiteList.add(i);
            }
            Integer[] hiCoverSite = hiCoverSiteList.toArray(new Integer[hiCoverSiteList.size()]);
            int siteNumOnCurrentChr = chrStartEnd[1]-chrStartEnd[0]+1;
            int base = siteNumOnCurrentChr%window;
            int numOfWindow;
            if (base == 0) numOfWindow = siteNumOnCurrentChr/window;
            else numOfWindow = siteNumOnCurrentChr/window + 1;
            int chrKeepCnt = 0;
            for (int i=0; i<numOfWindow; i++) {
                int actualWindow = window;
                if (base != 0 && i == numOfWindow-1) actualWindow = base;
                int actualStart = chrStartEnd[0]+window*i;
                int actualEnd = actualStart+actualWindow;
                ArrayList<Integer> taxaIndexList = new ArrayList();
                int hiCoverCntInWindow = 0;
                for (int j=0; j<gt.numberOfTaxa(); j++) {
                    int nonMissingSiteCnt = 0;
                    hiCoverCntInWindow = 0;
                    int hetCnt = 0;
                    for (int k=actualStart; k<actualEnd; k++) {
                        int hit = Arrays.binarySearch(hiCoverSite, k);
                        if (hit < 0) continue;
                        hiCoverCntInWindow++;
                        if (gt.genotype(j, k) != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                            nonMissingSiteCnt++;
                            if (gt.isHeterozygous(j, k)) hetCnt++;
                        }
                    }
                    if (nonMissingSiteCnt == 0) continue;
                    
                    double hetRatio = (double)hetCnt/nonMissingSiteCnt; // or by hiCoverCnt?
                    //double hetRatio = (double)hetCnt/hiCoverCnt;
                    if (hetRatio < maxHetForTaxon) taxaIndexList.add(j);
                }
                Integer[] selectedTaxaIndex = taxaIndexList.toArray(new Integer[taxaIndexList.size()]);
                System.out.println(hiCoverSite.length+"\t"+hiCoverCntInWindow+"\t");
                System.out.println(gt.positions().get(actualStart).getPosition() + "\t"+ actualStart +" Homo Taxa num is " + selectedTaxaIndex.length);
               
                //Site filters in selected homozygous taxa
                
                for (int j = actualStart; j < actualEnd; j++) {
                    int missingTaxaCnt = 0;
                    int nonMissingTaxaCnt = 0;
                    int hetCnt = 0;
                    for (int k = 0; k < selectedTaxaIndex.length; k++) {
                        if (gt.genotype(selectedTaxaIndex[k], j) != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) nonMissingTaxaCnt++;
                        else missingTaxaCnt++;
                        if (gt.isHeterozygous(selectedTaxaIndex[k], j)) hetCnt++;
                    }
                    if (nonMissingTaxaCnt == 0) continue;
                    double hetRatio = (double) hetCnt/nonMissingTaxaCnt;
                    if (hetRatio < maxHetForSite) {
                        siteKeepList.add(j);
                        chrKeepCnt++;
                    }
 /*                   
                    int hit = Arrays.binarySearch(hiCoverSite, j);
                    if (hit < 0) {
                        System.out.println(j+"\t"+missingTaxaCnt+"\t"+nonMissingTaxaCnt+"\t"+gt.totalNonMissingForSite(j)+"\t"+hetCnt+"\t"+0+"\t"+hetRatio+"\t"+gt.minorAlleleFrequency(j));
                    }
                    else {
                        System.out.println(j+"\t"+missingTaxaCnt+"\t"+nonMissingTaxaCnt+"\t"+gt.totalNonMissingForSite(j)+"\t"+hetCnt+"\t"+1+"\t"+hetRatio+"\t"+gt.minorAlleleFrequency(j));
                    }
 */                   
                }
            }    
            System.out.println("Keep "+String.valueOf(chrKeepCnt)+" sites. "+String.valueOf((double)chrKeepCnt/siteNumOnCurrentChr)+"\n");
        }
        int[] siteKeepArray = new int[siteKeepList.size()];
        for (int i = 0; i < siteKeepList.size(); i++) {
            siteKeepArray[i] = siteKeepList.get(i);
        }
        //FilterGenotypeTable fgt = FilterGenotypeTable.getInstance(gt, siteKeepArray);
        //ExportUtils.writeGenotypeHDF5(fgt, outputGenotype);
        //ExportUtils.writeToHapmap(fgt, "M:\\pipelineTest\\cassava\\gbs\\filtered_genotype\\chr1_filter.hmp.txt");
    }
    
    public void mkTaxaNameFileS (String genotypeFileS, String taxaNameFileS) {
        GenotypeTable gt = ImportUtils.readGuessFormat(genotypeFileS);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(taxaNameFileS), 65536);
            bw.write("Taxa_FullName");
            for (int i = 0; i < gt.numberOfTaxa(); i++) {
                bw.write(gt.taxaName(i));
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    public void mergeTaxaTBT (String sourceTBTFileS, String targetTBTFileS) {
        TagsByTaxaByte tbt = new TagsByTaxaByte(sourceTBTFileS, FilePacking.Byte);
        TreeSet<String> labelSet = new TreeSet();
        String[] temp = null;
        String label;
        for (int i = 0; i < tbt.getTaxaCount(); i++) {
            temp = tbt.getTaxaName(i).split(":");
            label = temp[0]+":"+temp[temp.length-1];
            labelSet.add(label);
        }
        String[] labels = labelSet.toArray(new String[labelSet.size()]);
        Arrays.sort(labels);
        String[][] labelFullName = new String[labels.length][];
        ArrayList<String>[] labelFullNameList = new ArrayList[labelFullName.length];
        for (int i = 0; i < labelFullNameList.length; i++) labelFullNameList[i] = new ArrayList();
        for (int i = 0; i < tbt.getTaxaCount(); i++) {
            temp = tbt.getTaxaName(i).split(":");
            label = temp[0]+":"+temp[temp.length-1];
            int hit = Arrays.binarySearch(labels, label);
            labelFullNameList[hit].add(tbt.getTaxaName(i));
        }
        for (int i = 0; i < labelFullName.length; i++) {
            labelFullName[i] = labelFullNameList[i].toArray(new String[labelFullNameList[i].size()]);
            Arrays.sort(labelFullName[i]);
        }
        int[][] index = new int[labelFullName.length][];
        for (int i = 0; i < index.length; i++) index[i] = new int[labelFullName[i].length];
        for (int i = 0; i < tbt.getTaxaCount(); i++) {
            for (int j = 0; j < labelFullName.length; j++) {
                int hit = Arrays.binarySearch(labelFullName[j], tbt.getTaxaName(i));
                if (hit < 0) continue;
                index[j][hit] = i;
            }
        }
        try {
            DataOutputStream dis = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(targetTBTFileS), 65536));
            dis.writeInt(tbt.getTagCount());
            dis.writeInt(tbt.getTagSizeInLong());
            dis.writeInt(labels.length);
            for (int i = 0; i < labels.length; i++) {
                if (labelFullName[i].length == 1) {
                    dis.writeUTF(labelFullName[i][0]);
                }
                else {
                    temp = labelFullName[i][0].split(":");
                    String newName = temp[0]+":"+temp[1]+":0:"+temp[3];
                    dis.writeUTF(newName);
                }
            }
            for (int i = 0; i < tbt.getTagCount(); i++) {
                long[] tag = tbt.getTag(i);
                for (int j = 0; j < tbt.getTagSizeInLong(); j++) {
                    dis.writeLong(tag[j]);
                }
                dis.writeByte(tbt.getTagLength(i));
                for (int j = 0; j < labels.length; j++) {
                    int count = 0;
                    for (int k = 0; k < labelFullName[j].length; k++) {
                        count += tbt.getReadCountForTagTaxon(i, index[j][k]);
                    }
                    if (count > Byte.MAX_VALUE) count = Byte.MAX_VALUE;
                    dis.writeByte(count);
                }
            }
            dis.flush();
            dis.close();  
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkTagCountFromTBT (String tbtFileS, String tagCountFileS) {
        TagsByTaxaByte tbt = new TagsByTaxaByte(tbtFileS, FilePacking.Byte);
        int minCount = 10000;
        try{
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(tagCountFileS), 65536));
            dos.writeInt(tbt.getTagCount());
            dos.writeInt(tbt.getTagSizeInLong());
            for (int i = 0; i < tbt.getTagCount(); i++) {
                long[] tag = tbt.getTag(i);
                for (int j = 0; j < tag.length; j++) {
                    dos.writeLong(tag[j]);
                }
                dos.writeByte(tbt.getTagLength(i));
                int count = 0;
                for (int j = 0; j < tbt.getTaxaCount(); j++) {
                    count += tbt.getReadCountForTagTaxon(i, j);
                }
                dos.writeInt(count);
                if (count < minCount) minCount = count;
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("tagCountFile is created from TBT at " + tagCountFileS);
        System.out.println("MinCount is " + String.valueOf(minCount));
    }
}
