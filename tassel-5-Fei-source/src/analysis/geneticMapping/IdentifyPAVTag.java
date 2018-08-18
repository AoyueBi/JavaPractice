/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package analysis.geneticMapping;

import format.Fasta;
import format.Sequence;
import format.ShortreadAlignment;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.map.TagsOnGeneticMap;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import utils.FStringUtils;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author Fei Lu <fl262@cornell.edu>
 */
public class IdentifyPAVTag {
    
    public IdentifyPAVTag () {
        //this.blastnPipe();
        this.bowtie2Pipe();
        //this.categaryPipe();
    }
    
    public void categaryPipe () {
        //this.togmToFastQ();
        //this.convertAlignment();
        this.summary();
    }
    
    public void summary () {
        String togmFileS = "E:\\Research\\geneticMapping\\identifyPAV\\anchor\\v1.panAnchor.togm.txt";
        String alignmentFileS = "E:\\Research\\geneticMapping\\identifyPAV\\anchor\\alignment\\panAnchor.k2.sa";
        String summaryFileS = "E:\\Research\\geneticMapping\\identifyPAV\\anchor\\categarySummary.txt";
        ShortreadAlignment sa = new ShortreadAlignment(alignmentFileS, IOFileFormat.Binary);
        TagsOnGeneticMap togm = new TagsOnGeneticMap(togmFileS, FilePacking.Text);
        sa.sortByQuery();
        String[] queries = sa.getQuerys();
        int noHit = 0;
        int b73Hit = 0;
        int b73PHit = 0;
        int nonHit = 0;
        int nonPHit = 0;
        for (int i = 0;i < queries.length; i++) {
            int index = sa.getAlignmentStartIndexByQuery(queries[i]);
            if (sa.isMatch(index)) {
                if (sa.isPerfectMatch(index)) {
                    if (togm.isPAV(Integer.valueOf(queries[i]))) {
                        b73PHit++;
                    }
                    else {
                        b73Hit++;
                    }
                }
                else {
                    if (togm.isPAV(Integer.valueOf(queries[i]))) {
                        nonPHit++;
                    }
                    else {
                        nonHit++;
                    }
                }
            }
            else {
                noHit++;
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(summaryFileS);
            bw.write("Tags do not align, but genetically mapped\t"+ String.valueOf((double)noHit/togm.getTagCount()));
            bw.newLine();
            bw.write("B73 tags whose physical positions agree with genetic positions\t" + String.valueOf((double)b73Hit/togm.getTagCount()));
            bw.newLine();
            bw.write("Non-B73 tags whose physical positions agree with genetic positions\t" + String.valueOf((double)nonHit/togm.getTagCount()));
            bw.newLine();
            bw.write("B73 tags whose physical positions do not agree with genetic positions\t" + String.valueOf((double)b73PHit/togm.getTagCount()));
            bw.newLine();
            bw.write("Non-B73 tags whose physical positions do not agree with genetic positions\t" + String.valueOf((double)nonPHit/togm.getTagCount()));
            bw.newLine();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void convertAlignment () {
        String samFileS = "E:\\Research\\geneticMapping\\identifyPAV\\anchor\\alignment\\panAnchor.sam";
        String alignmentFileS = "E:\\Research\\geneticMapping\\identifyPAV\\anchor\\alignment\\panAnchor.k2.sa";
        ShortreadAlignment sa = new ShortreadAlignment();
        sa.readFromBowtie2(samFileS);
        sa.writeSimpleAlignment(alignmentFileS, IOFileFormat.Binary);
    }
    
    public void togmToFastQ () {
        String togmFileS = "E:\\Research\\geneticMapping\\identifyPAV\\anchor\\v1.panAnchor.togm.txt";
        String fastqFileS = "E:\\Research\\geneticMapping\\identifyPAV\\anchor\\alignment\\panAnchor.fastq";
        new TagsOnGeneticMap(togmFileS, FilePacking.Text).writeFastQ(fastqFileS);
    }
    
    /**
     * Use this pipeline
     */
    public void bowtie2Pipe () {
        //this.mkBowtie2Lib();
        //this.alignByBowtie2();
        //this.identifyPAVBowtie2();
        this.identifyPAVByBoth();
    }
    
    /**
     * Too stringent
     * @deprecated 
     */
    public void identifyPAVByBoth () {
        String bowtie2AlignmentDirS = "E:\\Research\\geneticMapping\\identifyPAV\\regionAlignmentBowtie2\\";
        String blastAlignmentDirS = "E:\\Research\\geneticMapping\\identifyPAV\\regionAlignmentBlast\\";
        String anchorFileS = "M:\\pav\\PhyGenMapping\\v1.togm.txt";
        String newAnchorFileS = "E:\\Research\\geneticMapping\\identifyPAV\\anchor\\v1.togm.BothValidation.txt";
        
        File[] bowtie2Fs = new File(bowtie2AlignmentDirS).listFiles();
        File[] blastFs = new File(blastAlignmentDirS).listFiles();
        TagsOnGeneticMap togm = new TagsOnGeneticMap(anchorFileS, FilePacking.Text);
        for (int i = 0; i < togm.getTagCount(); i++) {
            togm.setIfPAV(i, 1);
        }
        for (int i = 0; i < bowtie2Fs.length; i++) {
            try {
                BufferedReader br = IoUtils.getTextReader(bowtie2Fs[i].getAbsolutePath());
                String temp;
                String[] tem;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("@")) continue;
                    tem = temp.split("\t");
                    if (tem[2].equals("*")) continue;
                    togm.setIfPAV(Integer.valueOf(tem[0]), 0);
                }
                if (i%100 == 0) System.out.println("Read in " + String.valueOf(i+1) + " files");
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        for (int i = 0; i < blastFs.length; i++) {
            try {
                BufferedReader br = IoUtils.getTextReader(blastFs[i].getAbsolutePath());
                String temp;
                String[] tem;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) continue;
                    tem = temp.split("\t");
                    togm.setIfPAV(Integer.valueOf(tem[0]), 0);
                }
                System.out.println("Read in " + String.valueOf(i+1) + " files");
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        int cnt = 0;
        for (int i = 0; i < togm.getTagCount(); i++) {
            cnt+=togm.getIfPAV(i);
        }
        System.out.println(String.valueOf(cnt)+" PAV identified");
        togm.writeDistFile(newAnchorFileS, FilePacking.Text);
    }
    
    public void identifyPAVBowtie2 () {
        String alignmentDirS = "E:\\Research\\geneticMapping\\identifyPAV\\regionAlignmentBowtie2\\";
        String anchorFileS = "M:\\pav\\PhyGenMapping\\v1.togm.txt";
        String newAnchorFileS = "E:\\Research\\geneticMapping\\identifyPAV\\anchor\\v1.togm.Bowtie2Validation.txt";
        
        File[] fs = new File(alignmentDirS).listFiles();
        TagsOnGeneticMap togm = new TagsOnGeneticMap(anchorFileS, FilePacking.Text);
        for (int i = 0; i < togm.getTagCount(); i++) {
            togm.setIfPAV(i, 1);
        }
        for (int i = 0; i < fs.length; i++) {
            try {
                BufferedReader br = IoUtils.getTextReader(fs[i].getAbsolutePath());
                String temp;
                String[] tem;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("@")) continue;
                    tem = temp.split("\t");
                    if (tem[2].equals("*")) continue;
                    togm.setIfPAV(Integer.valueOf(tem[0]), 0);
                }
                if (i%100 == 0) System.out.println("Read in " + String.valueOf(i+1) + " files");
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        int cnt = 0;
        for (int i = 0; i < togm.getTagCount(); i++) {
            cnt+=togm.getIfPAV(i);
        }
        System.out.println(String.valueOf(cnt)+" PAV identified");
        togm.writeDistFile(newAnchorFileS, FilePacking.Text);
        
    }
    
    public void alignByBowtie2 () {
        String genomeSliceDirS = "/workdir/mingh/bowtie2Lib/";
        String queryDirS = "/workdir/mingh/tagQueries/";
        String alignmentDirS = "/workdir/mingh/regionAlignmentBowtie2/";
        String script = "/workdir/mingh/runBowtie2.pl";
        String coreNumS = "8";
        new File (alignmentDirS).mkdir();
        File[] fs = new File(queryDirS).listFiles();
        File[][] regionQuery = new File[10][];
        int[] regionCnt = new int[10];
        for (int i = 0; i < fs.length; i++) {
            regionCnt[Integer.valueOf(fs[i].getName().split("_")[0])-1]++;
        }
        for (int i = 0; i < regionQuery.length; i++) {
            regionQuery[i] = new File[regionCnt[i]];
        }
        for (int i = 0; i < fs.length; i++) {
            String name = fs[i].getName();
            String[] temp = name.replaceFirst(".fa", "").split("_");
            regionQuery[Integer.valueOf(temp[0])-1][Integer.valueOf(temp[1])] = fs[i];
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(script);
            for (int i = 0; i < regionQuery.length; i++) {
                for (int j = 0; j < regionQuery[i].length; j++) {
                    String[] temp = regionQuery[i][j].getName().replaceFirst(".fa", "").split("_");
                    int chr = Integer.valueOf(temp[0]);
                    int regionIndex = Integer.valueOf(temp[1]);
                    int preRegionIndex = regionIndex-1;
                    int nextRegionIndex = regionIndex+1;
                    
                    if (regionIndex == 0) {
                        String dbFileS = new File(genomeSliceDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)).getAbsolutePath();
                        
                        String nextDbFileS = new File(genomeSliceDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, nextRegionIndex)).getAbsolutePath();
                        String outfileS = new File (alignmentDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+".sam").getAbsolutePath();
                        
                        String nextOutfileS = new File (alignmentDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+"_"+FStringUtils.getNDigitNumber(3, nextRegionIndex)+".sam").getAbsolutePath();
                        bw.write("system \"bowtie2 -x " + dbFileS + " -f " + regionQuery[i][j].getAbsolutePath()+" -M --very-sensitive-local -S "+ outfileS + " -p " + coreNumS + "\";");
                        bw.newLine();
                        
          
                        bw.write("system \"bowtie2 -x " + nextDbFileS + " -f " + regionQuery[i][j].getAbsolutePath()+" -M --very-sensitive-local -S "+ nextOutfileS + " -p " + coreNumS + "\";");
                        bw.newLine();
                    }
                    else if (regionIndex == regionQuery[i].length-1) {
                        String dbFileS = new File(genomeSliceDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)).getAbsolutePath();
                        String preDbFileS = new File(genomeSliceDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, preRegionIndex)).getAbsolutePath();
                        
                        String outfileS = new File (alignmentDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+".sam").getAbsolutePath();
                        String preOutfileS = new File (alignmentDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+"_"+FStringUtils.getNDigitNumber(3, preRegionIndex)+".sam").getAbsolutePath();
                        
                        bw.write("system \"bowtie2 -x " + dbFileS + " -f " + regionQuery[i][j].getAbsolutePath()+" -M --very-sensitive-local -S "+ outfileS + " -p " + coreNumS + "\";");
                        bw.newLine();
                        bw.write("system \"bowtie2 -x " + preDbFileS + " -f " + regionQuery[i][j].getAbsolutePath()+" -M --very-sensitive-local -S "+ preOutfileS + " -p " + coreNumS + "\";");
                        bw.newLine();
                      
                        
                    }
                    else {
                        String dbFileS = new File(genomeSliceDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)).getAbsolutePath();
                        String preDbFileS = new File(genomeSliceDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, preRegionIndex)).getAbsolutePath();
                        String nextDbFileS = new File(genomeSliceDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, nextRegionIndex)).getAbsolutePath();
                        String outfileS = new File (alignmentDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+".sam").getAbsolutePath();
                        String preOutfileS = new File (alignmentDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+"_"+FStringUtils.getNDigitNumber(3, preRegionIndex)+".sam").getAbsolutePath();
                        String nextOutfileS = new File (alignmentDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+"_"+FStringUtils.getNDigitNumber(3, nextRegionIndex)+".sam").getAbsolutePath();
                        bw.write("system \"bowtie2 -x " + dbFileS + " -f " + regionQuery[i][j].getAbsolutePath()+" -M --very-sensitive-local -S "+ outfileS + " -p " + coreNumS + "\";");
                        bw.newLine();
                        bw.write("system \"bowtie2 -x " + preDbFileS + " -f " + regionQuery[i][j].getAbsolutePath()+" -M --very-sensitive-local -S "+ preOutfileS + " -p " + coreNumS + "\";");
                        bw.newLine();
                        bw.write("system \"bowtie2 -x " + nextDbFileS + " -f " + regionQuery[i][j].getAbsolutePath()+" -M --very-sensitive-local -S "+ nextOutfileS + " -p " + coreNumS + "\";");
                        bw.newLine();
                    }
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkBowtie2Lib () {
        String genomeDirS = "/workdir/mingh/bowtie2Lib/";
        String script = "/workdir/mingh/mkLib.pl";
        
        File[] fs = new File(genomeDirS).listFiles();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(script);
            for (int i = 0; i < fs.length; i++) {
                String cmd = "system \"bowtie2-build " + fs[i].getAbsolutePath() + " " + fs[i].getAbsolutePath().replaceFirst(".fa", "")+"\";";
                bw.write(cmd);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * low-complexity sequence are not aligned by blastn,so not using this pipeline
     */
    public void blastnPipe () {
        //this.splitGenomes();
        this.mkTagQueries();
        //this.alignRegions();
        //this.identifyPAV();
    }
    
    public void identifyPAV () {
        String alignmentDirS = "E:\\Research\\geneticMapping\\identifyPAV\\regionAlignment\\";
        String anchorFileS = "M:\\pav\\PhyGenMapping\\v1.togm.txt";
        String newAnchorFileS = "E:\\Research\\geneticMapping\\identifyPAV\\anchor\\v1.togm.BlastValidation.txt";
        double eCut = 1e-5;
        File[] fs = new File(alignmentDirS).listFiles();
        TagsOnGeneticMap togm = new TagsOnGeneticMap(anchorFileS, FilePacking.Text);
        for (int i = 0; i < togm.getTagCount(); i++) {
            togm.setIfPAV(i, 1);
        }
        for (int i = 0; i < fs.length; i++) {
            try {
                BufferedReader br = IoUtils.getTextReader(fs[i].getAbsolutePath());
                String temp;
                String[] tem;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) continue;
                    tem = temp.split("\t");
                    if (Double.valueOf(tem[10]) > eCut) continue;
                    togm.setIfPAV(Integer.valueOf(tem[0]), 0);
                }
                System.out.println("Read in " + String.valueOf(i+1) + " files");
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        int cnt = 0;
        for (int i = 0; i < togm.getTagCount(); i++) {
            cnt+=togm.getIfPAV(i);
        }
        System.out.println(String.valueOf(cnt)+" PAV identified");
        togm.writeDistFile(newAnchorFileS, FilePacking.Text);
    }
    
    public void alignRegions () {
        //String fastaFileS = "E:\\Research\\geneticMapping\\identifyPAV\\tagQueries\\";
        //String genomeSliceDirS = "E:\\Research\\geneticMapping\\identifyPAV\\genomeSlice\\";
        //String alignmentDirS = "E:\\Research\\geneticMapping\\identifyPAV\\regionAlignment\\";
        //String scriptFileS = "E:\\Research\\geneticMapping\\identifyPAV\\runBlast.pl";
        
        String fastaFileS = "/workdir/mingh/tagQueries/";
        String genomeSliceDirS = "/workdir/mingh/genomeSlice/";
        String alignmentDirS = "/workdir/mingh/regionAlignment/";
        String scriptFileS = "/workdir/mingh/runBlast.pl";

        String coreNumS = "8";
        File[] fs = new File(fastaFileS).listFiles();
        File[][] regionQuery = new File[10][];
        int[] regionCnt = new int[10];
        for (int i = 0; i < fs.length; i++) {
            regionCnt[Integer.valueOf(fs[i].getName().split("_")[0])-1]++;
        }
        for (int i = 0; i < regionQuery.length; i++) {
            regionQuery[i] = new File[regionCnt[i]];
        }
        for (int i = 0; i < fs.length; i++) {
            String name = fs[i].getName();
            String[] temp = name.replaceFirst(".fa", "").split("_");
            regionQuery[Integer.valueOf(temp[0])-1][Integer.valueOf(temp[1])] = fs[i];
        }
        new File(alignmentDirS).mkdir();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(scriptFileS);
            for (int i = 0; i < regionQuery.length; i++) {
                for (int j = 0; j < regionQuery[i].length; j++) {
                    String[] temp = regionQuery[i][j].getName().replaceFirst(".fa", "").split("_");
                    int chr = Integer.valueOf(temp[0]);
                    int regionIndex = Integer.valueOf(temp[1]);
                    int preRegionIndex = regionIndex-1;
                    int nextRegionIndex = regionIndex+1;
                    
                    if (regionIndex == 0) {
                        String dbFileS = new File(genomeSliceDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+".fa").getAbsolutePath();
                       
                        String nextDbFileS = new File(genomeSliceDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, nextRegionIndex)+".fa").getAbsolutePath();
                        String outfileS = new File (alignmentDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+".aln.txt").getAbsolutePath();
                        
                        String nextOutfileS = new File (alignmentDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+"_"+FStringUtils.getNDigitNumber(3, nextRegionIndex)+".aln.txt").getAbsolutePath();
                        bw.write(("system \"blastn -query "+regionQuery[i][j].getAbsolutePath()+" -db "+dbFileS+" -out "+outfileS+" -outfmt 6 -num_threads "+coreNumS+"\";").replaceAll("\\\\", "/"));
                        bw.newLine();
                        bw.write(("system \"blastn -query "+regionQuery[i][j].getAbsolutePath()+" -db "+nextDbFileS+" -out "+nextOutfileS+" -outfmt 6 -num_threads "+coreNumS+"\";").replaceAll("\\\\", "/"));
                        bw.newLine();
                    }
                    else if (regionIndex == regionQuery[i].length-1) {
                        String dbFileS = new File(genomeSliceDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+".fa").getAbsolutePath();
                        String preDbFileS = new File(genomeSliceDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, preRegionIndex)+".fa").getAbsolutePath();
                        
                        String outfileS = new File (alignmentDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+".aln.txt").getAbsolutePath();
                        String preOutfileS = new File (alignmentDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+"_"+FStringUtils.getNDigitNumber(3, preRegionIndex)+".aln.txt").getAbsolutePath();
                        
                        bw.write(("system \"blastn -query "+regionQuery[i][j].getAbsolutePath()+" -db "+dbFileS+" -out "+outfileS+" -outfmt 6 -num_threads "+coreNumS+"\";").replaceAll("\\\\", "/"));
                        bw.newLine();
                        bw.write(("system \"blastn -query "+regionQuery[i][j].getAbsolutePath()+" -db "+preDbFileS+" -out "+preOutfileS+" -outfmt 6 -num_threads "+coreNumS+"\";").replaceAll("\\\\", "/"));
                        bw.newLine();
                        
                    }
                    else {
                        String dbFileS = new File(genomeSliceDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+".fa").getAbsolutePath();
                        String preDbFileS = new File(genomeSliceDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, preRegionIndex)+".fa").getAbsolutePath();
                        String nextDbFileS = new File(genomeSliceDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, nextRegionIndex)+".fa").getAbsolutePath();
                        String outfileS = new File (alignmentDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+".aln.txt").getAbsolutePath();
                        String preOutfileS = new File (alignmentDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+"_"+FStringUtils.getNDigitNumber(3, preRegionIndex)+".aln.txt").getAbsolutePath();
                        String nextOutfileS = new File (alignmentDirS, FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, regionIndex)+"_"+FStringUtils.getNDigitNumber(3, nextRegionIndex)+".aln.txt").getAbsolutePath();
                        bw.write(("system \"blastn -query "+regionQuery[i][j].getAbsolutePath()+" -db "+dbFileS+" -out "+outfileS+" -outfmt 6 -num_threads "+coreNumS+"\";").replaceAll("\\\\", "/"));
                        bw.newLine();
                        bw.write(("system \"blastn -query "+regionQuery[i][j].getAbsolutePath()+" -db "+preDbFileS+" -out "+preOutfileS+" -outfmt 6 -num_threads "+coreNumS+"\";").replaceAll("\\\\", "/"));
                        bw.newLine();
                        bw.write(("system \"blastn -query "+regionQuery[i][j].getAbsolutePath()+" -db "+nextDbFileS+" -out "+nextOutfileS+" -outfmt 6 -num_threads "+coreNumS+"\";").replaceAll("\\\\", "/"));
                        bw.newLine();
                    }
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkTagQueries () {
        String referenceGenome = "N:\\Zea\\reference_genome\\ZmB73_RefGen_v2.fa";
        String anchorFileS = "M:\\pav\\PhyGenMapping\\v1.togm.txt";
        String fastaFileS = "E:\\Research\\geneticMapping\\identifyPAV\\tagQueries\\";
        
        Fasta f = new Fasta(referenceGenome);
        BufferedWriter[][] bws = new BufferedWriter[f.getSeqNumber()][];
        int[][] regionStart = new int[f.getSeqNumber()][];
        for (int i = 0; i < f.getSeqNumber(); i++) {
            int chrIndex = Integer.valueOf(f.getName(i))-1;
            System.out.println(chrIndex+"\t"+f.getSeqLength(i));
            int numStream = 0;
            int left = f.getSeqLength(i)%10000000;
            if (left == 0) numStream = f.getSeqLength(i)/10000000;
            else numStream = f.getSeqLength(i)/10000000+1;
            regionStart[chrIndex] = new int[numStream];
            bws[chrIndex] = new BufferedWriter[numStream];
            for (int j = 0; j < numStream; j++) {
                regionStart[chrIndex][j] = j*10000000+1;
                String outfileS = FStringUtils.getNDigitNumber(3, chrIndex+1)+"_"+FStringUtils.getNDigitNumber(3, j)+".fa";
                outfileS = fastaFileS + outfileS;
                bws[chrIndex][j] = IoUtils.getTextWriter(outfileS);
            }
        }
        TagsOnGeneticMap togm = new TagsOnGeneticMap(anchorFileS, FilePacking.Text);
        try {
            for (int i = 0; i < togm.getTagCount(); i++) {
                int chrIndex = togm.getGChr(i)-1;
                int regionIndex = Arrays.binarySearch(regionStart[chrIndex], togm.getGPos(i));
                if (regionIndex < 0) regionIndex = -regionIndex-2;
                bws[chrIndex][regionIndex].write(">"+String.valueOf(i));
                bws[chrIndex][regionIndex].newLine();
                long[] t = togm.getTag(i);
                bws[chrIndex][regionIndex].write(BaseEncoder.getSequenceFromLong(t).substring(0, togm.getTagLength(i)));
                bws[chrIndex][regionIndex].newLine();
            }
            for (int i = 0; i < bws.length; i++) {
                for (int j = 0; j < bws[i].length; j++) {
                    bws[i][j].flush();
                    bws[i][j].close();
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void splitGenomes () {
        String referenceGenome = "N:\\Zea\\reference_genome\\ZmB73_RefGen_v2.fa";
        String genomeSliceDirS = "E:\\Research\\geneticMapping\\identifyPAV\\genomeSlice\\";
        String script = "E:\\Research\\geneticMapping\\identifyPAV\\mkLib.pl";
        Fasta f = new Fasta (referenceGenome);
        ArrayList<String> scList = new ArrayList();
        for (int i = 0; i < f.getSeqNumber(); i++) {
            String[] slices = new Sequence(f.getSeq(i)).getFragments(10000000);
            int chr = Integer.valueOf(f.getName(i));
            for (int j = 0; j < slices.length; j++) {
                String outName = FStringUtils.getNDigitNumber(3, chr)+"_"+FStringUtils.getNDigitNumber(3, j)+".fa";
                String outfileS = new File(genomeSliceDirS, outName).getAbsolutePath();
                try {
                    BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                    bw.write(">"+outName.replaceFirst(".fa", ""));
                    bw.newLine();
                    bw.write(slices[j]);
                    bw.newLine();
                    bw.flush();
                    bw.close();
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
                String command = "system \"makeblastdb -in " + outfileS.replaceAll("\\\\", "/") + " -dbtype nucl\";";
                scList.add(command);
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(script);
            for (int i = 0; i < scList.size(); i++) {
                bw.write(scList.get(i));
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
