/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.maf;

import format.Ranges;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import graphcis.r.DensityPlot;
import graphcis.r.Histogram;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import utils.FArrayUtils;
import utils.FStringUtils;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author Fei Lu
 */
class WGSDepth {
    
    public WGSDepth () {
        //this.checkFormat(); //check Robert's HM3.1 alleleDepth format
        
        this.convertAlleleDepthChr10(true);
        //this.extractInvariantSite();
        
        //this.convertAlleleDepth();
        //this.convertAlleleDepthUniqueGenome();
        
    }
    
    public void convertAlleleDepthChr10 (boolean ifUniqueGenome) {
        int minDepth = 1000;
        int maxDepth = 4000;
        String uniqueGenomeFileS = "M:\\production\\maf\\wgs\\uniqueKmer\\position\\uniqueGenome.pos.r.txt";
        String infileS = "M:\\production\\maf\\wgs\\dataCheck\\thr.all_c10_DEPTH_hmp31_bcq10_q30";
        String outfileDirS = "M:\\production\\maf\\wgs\\alleleDepth\\";
        File infile = new File(infileS);
        AlleleDepth ad = new AlleleDepth ();
        ad.readFromHapMap(infile.getAbsolutePath());
        Ranges r;
        TIntArrayList indexList = new TIntArrayList();
        if (ifUniqueGenome) {
            r = new Ranges(uniqueGenomeFileS, IOFileFormat.Text);
            for (int j = 0; j < ad.getSiteNumber(); j++) {
                if (ad.getSiteDepth(j) < minDepth) continue;
                if (ad.getSiteDepth(j) > maxDepth) continue;
                if (!r.isInRanges(ad.getChromosome(j), ad.getPosition(j))) continue;
                indexList.add(j);
            }
        }
        else {
            for (int j = 0; j < ad.getSiteNumber(); j++) {
                if (ad.getSiteDepth(j) < minDepth) continue;
                if (ad.getSiteDepth(j) > maxDepth) continue;
                indexList.add(j);
            }
        }
        AlleleDepth newAD = ad.getSubset(indexList.toArray(new int[indexList.size()]));
        int chr = newAD.getChromosome(0);
        String chrS = FStringUtils.getNDigitNumber(2, chr);
        String outfileS = new File(outfileDirS, "minDepth_hmp31_ad_"+chrS+".bin").getAbsolutePath();
        newAD.writeAlleleDepth(outfileS, IOFileFormat.Binary);
        System.out.println("Depth of chromosome " + String.valueOf(chr)+ " is written");
    }
    
    public void extractInvariantSite () {
        double maxMaf = 0.001;
        String infileDirS = "M:\\production\\maf\\wgs\\alleleDepth\\";
        String outfileDirS = "M:\\production\\maf\\wgs\\invariantSite\\";
        File[] infiles = new File(infileDirS).listFiles();
        for (int i = 0; i < infiles.length; i++) {
            AlleleDepth ad = new AlleleDepth (infiles[i].getAbsolutePath(), IOFileFormat.Binary);
            TIntArrayList indexList = new TIntArrayList();
            for (int j = 0; j < ad.getSiteNumber(); j++) {
                if (ad.getCorrectedMinorAlleleFrequency(j) > maxMaf) continue;
                indexList.add(j);
            }
            AlleleDepth newAD = ad.getSubset(indexList.toArray(new int[indexList.size()]));
            int chr = newAD.getChromosome(0);
            String chrS = FStringUtils.getNDigitNumber(2, chr);
            String outfileS = new File(outfileDirS, "invariant_hmp31_ad_"+chrS+".bin").getAbsolutePath();
            newAD.writeAlleleDepth(outfileS, IOFileFormat.Binary);
            System.out.println("Depth of chromosome " + String.valueOf(chr)+ " is written");
        }
    }
    
    /**
     * Only convert sites in unique portion of genome
     * The problem is really conserved sites have two copies (genome duplication), will be deleted
     */
    public void convertAlleleDepthUniqueGenome () {
        int minDepth = 1000;
        int maxDepth = 4000;
        String uniqueGenomeFileS = "M:\\production\\maf\\wgs\\uniqueKmer\\position\\uniqueGenome.pos.r.txt";
        String infileDirS = "N:\\HapMap3\\hmp31_allele_depths\\";
        String outfileDirS = "M:\\production\\maf\\wgs\\alleleDepth\\";
        File[] infiles = new File(infileDirS).listFiles();
        Ranges r = new Ranges(uniqueGenomeFileS, IOFileFormat.Text);
        for (int i = 0; i < infiles.length; i++) {
            AlleleDepth ad = new AlleleDepth ();
            ad.readFromHapMap(infiles[i].getAbsolutePath());
            TIntArrayList indexList = new TIntArrayList();
            for (int j = 0; j < ad.getSiteNumber(); j++) {
                if (ad.getSiteDepth(j) < minDepth) continue;
                if (ad.getSiteDepth(j) > maxDepth) continue;
                if (!r.isInRanges(ad.getChromosome(j), ad.getPosition(j))) continue;
                indexList.add(j);
            }
            AlleleDepth newAD = ad.getSubset(indexList.toArray(new int[indexList.size()]));
            int chr = newAD.getChromosome(0);
            String chrS = FStringUtils.getNDigitNumber(2, chr);
            String outfileS = new File(outfileDirS, "minDepth_hmp31_ad_"+chrS+".bin").getAbsolutePath();
            newAD.writeAlleleDepth(outfileS, IOFileFormat.Binary);
            System.out.println("Depth of chromosome " + String.valueOf(chr)+ " is written");
        }
    }
    
    /**
     * High coverage can be regions from repetitive sequences
     *
     */
    public void convertAlleleDepth () {
        int minDepth = 1000;
        int maxDepth = 4000;
        String infileDirS = "N:\\HapMap3\\hmp31_allele_depths\\";
        String outfileDirS = "M:\\production\\maf\\wgs\\alleleDepth\\";
        File[] infiles = new File(infileDirS).listFiles();
        for (int i = 0; i < infiles.length; i++) {
            AlleleDepth ad = new AlleleDepth ();
            ad.readFromHapMap(infiles[i].getAbsolutePath());
            TIntArrayList indexList = new TIntArrayList();
            for (int j = 0; j < ad.getSiteNumber(); j++) {
                if (ad.getSiteDepth(j) < minDepth) continue;
                if (ad.getSiteDepth(j) > maxDepth) continue;
                indexList.add(j);
            }
            AlleleDepth newAD = ad.getSubset(indexList.toArray(new int[indexList.size()]));
            int chr = newAD.getChromosome(0);
            String chrS = FStringUtils.getNDigitNumber(2, chr);
            String outfileS = new File(outfileDirS, "minDepth_hmp31_ad_"+chrS+".bin").getAbsolutePath();
            newAD.writeAlleleDepth(outfileS, IOFileFormat.Binary);
            System.out.println("Depth of chromosome " + String.valueOf(chr)+ " is written");
        }
    }
    
    void checkFormat () {
        int siteNum = 10000;
        String depthFileS = "M:\\production\\maf\\wgs\\alleleDepth\\thr.all_c10_DEPTH_hmp31_bcq10_q30.gz";
        String sampleFileS = "M:\\production\\maf\\wgs\\dataCheck\\depthSample.txt";
        try {
            BufferedReader br = IoUtils.getTextGzipReader(depthFileS);
            BufferedWriter bw = IoUtils.getTextWriter(sampleFileS);
            for (int i = 0; i < siteNum; i++) {
                String temp = br.readLine();
                bw.write(temp);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            
        }
    }
}
