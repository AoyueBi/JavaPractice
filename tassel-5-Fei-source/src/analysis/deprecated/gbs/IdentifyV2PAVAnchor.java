/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package analysis.deprecated.gbs;

import analysis.geneticMapping.*;
import format.Fasta;
import format.Sequence;
import format.ShortreadAlignment;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.map.TagGWASMap;
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
public class IdentifyV2PAVAnchor {
    
    public IdentifyV2PAVAnchor () {
        //this.splitGenomes();
        //this.mkTagQueries();
        //this.mkBowtie2Lib();
        //this.alignByBowtie2();
        //this.identifyPAVBowtie2();
        //this.resolutionEstimation();
    }
    
    public void resolutionEstimation () {
        new PanAnchorUtils().resolutionEstimation();
    }
    
    public void identifyPAVBowtie2 () {
        String alignmentDirS = "M:\\production\\panAnchorV2\\rigid\\regionAlignmentBowtie2\\";
        String anchorFileS = "M:\\production\\panAnchorV2\\rigid\\v2.panAnchor.rigid.togm.txt";
        String newAnchorFileS = "M:\\production\\panAnchorV2\\rigid\\v2.panAnchor.rigid.pav.togm.txt";
        
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
    
    
    public void mkTagQueries () {
        String referenceGenome = "N:\\Zea\\reference_genome\\ZmB73_RefGen_v2.fa";
        String anchorFileS = "M:\\production\\panAnchorV2\\rigid\\v2.panAnchor.rigid.togm.txt";
        String fastaFileS = "M:\\production\\panAnchorV2\\rigid\\tagQueries\\";
        
        int binSize = 10000000;
        Fasta f = new Fasta(referenceGenome);
        BufferedWriter[][] bws = new BufferedWriter[f.getSeqNumber()][];
        int[][] regionStart = new int[f.getSeqNumber()][];
        for (int i = 0; i < f.getSeqNumber(); i++) {
            int chrIndex = Integer.valueOf(f.getName(i))-1;
            System.out.println(chrIndex+"\t"+f.getSeqLength(i));
            int numStream = 0;
            int left = f.getSeqLength(i)%binSize;
            if (left == 0) numStream = f.getSeqLength(i)/binSize;
            else numStream = f.getSeqLength(i)/binSize+1;
            regionStart[chrIndex] = new int[numStream];
            bws[chrIndex] = new BufferedWriter[numStream];
            for (int j = 0; j < numStream; j++) {
                regionStart[chrIndex][j] = j*binSize+1;
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
    }
}
