/*
* To change this license header, choose License Headers in Project Properties.
* To change this template file, choose Tools | Templates
* and open the template in the editor.
*/

package analysis.geneticMapping;

import format.Bins;
import format.Table;
import graphcis.r.ScatterPlot;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeSet;
import net.maizegenetics.analysis.gbs.TagAgainstAnchor;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.tag.TagCounts;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import net.maizegenetics.dna.tag.TagsByTaxaByte;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import utils.IoUtils;


/**
 *
 * @author Fei Lu <fl262@cornell.edu>
 */
public class AccuracyFactor {
    
    public AccuracyFactor () {
        //this.compareImputation();
        //this.mkDepthOfTBTs();
        //this.mkDepthFile();
        //this.mkUB73TrainTagTaxa();
        //this.check();
        //this.mkImputationBin();
        this.mkFigures();
        
    }
    
    public void mkFigures () {
        String tagTaxaDepthFileS = "E:\\Research\\geneticMapping\\accuracyFactor\\UB73train.taxaDepth.txt";
        String UB73TrainTransFileS = "M:\\pav\\train\\UB73train_box.csv";
        String UB73TrainFileS = "M:\\pav\\train\\UB73train.mpgmap.txt";
        String imputationBinS = "E:\\Research\\geneticMapping\\accuracyFactor\\imputationBin.txt";
        String accVsTaxaDepthFileS = "E:\\Research\\geneticMapping\\accuracyFactor\\accTaxaDepth.pdf";
        String accVsTaxaCountFileS = "E:\\Research\\geneticMapping\\accuracyFactor\\accTaxaCount.pdf";
        String accVsMissingFileS = "E:\\Research\\geneticMapping\\accuracyFactor\\accMissing.pdf";
        String accVsImputeFileS = "E:\\Research\\geneticMapping\\accuracyFactor\\accImpute.pdf";
        String accFactorFileS = "E:\\Research\\geneticMapping\\accuracyFactor\\accFactor.txt";
        int pointNum = 5000;
        
        Table t = new Table (imputationBinS);
        HashMap<String, Double> missingMap = new HashMap();
        HashMap<String, Double> imputeMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            String s = t.content[i][0]+"_"+t.content[i][1];
            missingMap.put(s, Double.valueOf(t.content[i][2]));
            imputeMap.put(s, Double.valueOf(t.content[i][3]));
        }
        
        t = new Table (UB73TrainTransFileS, ",");
        if (pointNum > t.getRowNumber()) pointNum = t.getRowNumber();
        double ratio = (double)14129/23463;
        double[] distance = new double[pointNum];
        double[] taxaCount = new double[pointNum];
        double[] taxaDepth = new double[pointNum];
        double[] missing = new double[pointNum];
        double[] impute = new double[pointNum];
        for (int i = 0; i < pointNum; i++) {
            distance[i] = Double.valueOf(t.content[i][t.getColumnNumber()-1]);
        }
        t = new Table (UB73TrainFileS);
        for (int i = 0; i < pointNum; i++) {
            taxaCount[i] = Double.valueOf(t.content[i][12]) * ratio;
            String s = t.content[i][7]+"_"+t.content[i][9];
            missing[i] = missingMap.get(s);
            impute[i] = imputeMap.get(s);
        }
        t = new Table (tagTaxaDepthFileS);
        for (int i = 0; i < pointNum; i++) {
            taxaDepth[i] = Math.log10(Double.valueOf(t.content[i][1]));
        }
        
        ScatterPlot s = new ScatterPlot(distance, taxaDepth);
        s.setColor(255, 0, 0, 60);
        s.setPlottingCharacter(16);
        s.setXLab("Log10 (distance)");
        s.setYLab("Log10 (read count of all inbreds having the tag)");
        s.setTitle("");
        s.setSlideMode();
        s.saveGraph(accVsTaxaDepthFileS);
        s = new ScatterPlot(distance, taxaCount);
        s.setColor(255, 0, 0, 60);
        s.setPlottingCharacter(16);
        s.setXLab("Log10 (distance)");
        s.setYLab("Number of inbreds having the tag");
        s.setTitle("");
        s.setSlideMode();
        s.saveGraph(accVsTaxaCountFileS);
        s = new ScatterPlot(distance, missing);
        s.setColor(255, 0, 0, 60);
        s.setPlottingCharacter(16);
        s.setXLab("Log10 (distance)");
        s.setYLab("Proportion of missing genotype");
        s.setTitle("");
        s.setSlideMode();
        s.saveGraph(accVsMissingFileS);
        s = new ScatterPlot(distance, impute);
        s.setColor(255, 0, 0, 60);
        s.setPlottingCharacter(16);
        s.setXLab("Log10 (distance)");
        s.setYLab("Proportion of imputed genotype");
        s.setTitle("");
        s.setSlideMode();
        s.saveGraph(accVsImputeFileS);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(accFactorFileS);
            bw.write("Distance\tTaxaDepth\tTaxaCount\tMissing\tImputation");
            bw.newLine();
            for (int i = 0; i < pointNum; i++) {
                bw.write(String.valueOf(distance[i])+"\t"+String.valueOf(taxaDepth[i])+"\t"+String.valueOf(taxaCount[i])+"\t"+String.valueOf(missing[i])+"\t"+String.valueOf(impute[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        double r;
        double pv;
        try {
            PearsonsCorrelation p = new PearsonsCorrelation();
            double[][] data = new double[distance.length][2];
            for (int i = 0; i < data.length; i++) {
                data[i][0]  = distance[i];
                data[i][1] = taxaDepth[i];
            }
            p = new PearsonsCorrelation(data);
            r = Math.pow(p.getCorrelationMatrix().getEntry(0, 1),2);
            pv = p.getCorrelationPValues().getEntry(0, 1);
            System.out.println("R2 of distance Vs taxaDepth:\t" + r + "\t"+pv);
            for (int i = 0; i < data.length; i++) {
                data[i][0]  = distance[i];
                data[i][1] = taxaCount[i];
            }
            p = new PearsonsCorrelation(data);
            r = Math.pow(p.getCorrelationMatrix().getEntry(0, 1),2);
            pv = p.getCorrelationPValues().getEntry(0, 1);
            System.out.println("R2 of distance Vs taxaCount:\t" + r + "\t"+pv);
             for (int i = 0; i < data.length; i++) {
                data[i][0]  = distance[i];
                data[i][1] = missing[i];
            }
            p = new PearsonsCorrelation(data);
            r = Math.pow(p.getCorrelationMatrix().getEntry(0, 1),2);
            pv = p.getCorrelationPValues().getEntry(0, 1);
            System.out.println("R2 of distance Vs missing:\t" + r + "\t"+pv);
             for (int i = 0; i < data.length; i++) {
                data[i][0]  = distance[i];
                data[i][1] = impute[i];
            }
            p = new PearsonsCorrelation(data);
            r = Math.pow(p.getCorrelationMatrix().getEntry(0, 1),2);
            pv = p.getCorrelationPValues().getEntry(0, 1);
            System.out.println("R2 of distance Vs impute:\t" + r + "\t"+pv);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkImputationBin () {
        String imputationFileS = "E:\\Research\\geneticMapping\\accuracyFactor\\imputationCompare.txt";
        String imputationBinS = "E:\\Research\\geneticMapping\\accuracyFactor\\imputationBin.txt";
        String inforFileS = "E:\\Database\\InfoFile\\ChrLenCentPosi.txt";
        int binSize = 200000;
        Table t = new Table (inforFileS);
        int chrNum = t.getRowNumber();
        int[] chrLength = new int[chrNum];
        for (int i = 0; i < chrNum; i++) {
            chrLength[i] = Integer.valueOf(t.content[i][1]);
        }
        int[][] pos = new int[chrNum][];
        int[] snpCount = new int[chrNum];
        double[][] missing = new double[chrNum][];
        double[][] impute = new double[chrNum][];
        t = new Table(imputationFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            snpCount[Integer.valueOf(t.content[i][0])-1]++;
        }
        for (int i = 0; i < chrNum; i++) {
            pos[i] = new int[snpCount[i]];
            missing[i] = new double[snpCount[i]];
            impute[i] = new double[snpCount[i]];
            snpCount[i] = 0;
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            int chrIndex = Integer.valueOf(t.content[i][0])-1;
            pos[chrIndex][snpCount[chrIndex]] = Integer.valueOf(t.content[i][1]);
            missing[chrIndex][snpCount[chrIndex]] = Double.valueOf(t.content[i][3]);
            impute[chrIndex][snpCount[chrIndex]] = Double.valueOf(t.content[i][7]);
            snpCount[chrIndex]++;
        }
        double[][] averMissing = new double[chrNum][];
        double[][] averImpute = new double[chrNum][];
        int[][] binStart = new int[chrNum][];
        for (int i = 0; i < chrNum; i++) {
            Bins b = new Bins(1, chrLength[i], binSize, pos[i], missing[i]);
            averMissing[i] = new double[b.getBinNum()];
            for (int j = 0; j < averMissing[i].length; j++) {
                averMissing[i][j] = b.getBinAverage(j);
            }
            b = new Bins(1, chrLength[i], binSize, pos[i], impute[i]);
            averImpute[i] = new double[b.getBinNum()];
            binStart[i] = new int[b.getBinNum()];
            for (int j = 0; j < averImpute[i].length; j++) {
                averImpute[i][j] = b.getBinAverage(j);
                binStart[i][j] = b.getBinStart(j);
            }
            
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(imputationBinS);
            bw.write("Chromosome\tSite\tBinAverMissing\tBinAverImpute");
            bw.newLine();
            for (int i = 0; i < chrNum; i++) {
                for (int j = 0; j < missing[i].length; j++) {
                    int index = Arrays.binarySearch(binStart[i], pos[i][j]);
                    if (index < 0) {
                        index = -index-2;
                    }
                    bw.write(String.valueOf(i+1)+"\t"+String.valueOf(pos[i][j])+"\t"+String.valueOf(averMissing[i][index])+"\t"+String.valueOf(averImpute[i][index]));
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void check () {
        String tagTaxaDepthFileS = "E:\\Research\\geneticMapping\\accuracyFactor\\UB73train.taxaDepth.txt";
        String UB73TrainFileS = "M:\\pav\\train\\UB73train.mpgmap.txt";
        Table t = new Table(tagTaxaDepthFileS);
        Table tt = new Table(UB73TrainFileS);
        TagCounts tc = new TagCounts(2, t.getRowNumber());
        for (int i = 0; i < t.getRowNumber(); i++) {
            tc.setTag(BaseEncoder.getLongArrayFromSeq(t.content[i][0]), (byte)64, 1, i);
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            long[] tag = tc.getTag(i);
            if (!BaseEncoder.getSequenceFromLong(tag).equals(tt.content[i][0])) {
                System.out.println("wrong");
                System.out.println(tt.content[i][0]+"\n");
            }
        }
    }
    
    public void mkUB73TrainTagTaxa () {
        String tbtFileS = "M:\\pav\\mergedTBT\\UB73.tbt.byte";
        String UB73TrainFileS = "M:\\pav\\train\\UB73train.mpgmap.txt";
        String taxaDepthFileS = "E:\\Research\\geneticMapping\\accuracyFactor\\taxaDepth.txt";
        String tagTaxaDepthFileS = "E:\\Research\\geneticMapping\\accuracyFactor\\UB73train.taxaDepth.txt";
        
        Table t = new Table (taxaDepthFileS);
        HashMap<String,Integer> depthMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            depthMap.put(t.content[i][0], Integer.valueOf(t.content[i][1]));
        }
        t = new Table (UB73TrainFileS);
        ArrayList<String> l = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) l.add(t.content[i][0]);
        String[] tags = l.toArray(new String[l.size()]);
        Arrays.sort(tags);
        try {
            DataInputStream dis = IoUtils.getBinaryReader(tbtFileS);
            BufferedWriter bw = IoUtils.getTextWriter(tagTaxaDepthFileS);
            bw.write("Tag\tDepth");
            bw.newLine();
            int tagNum = dis.readInt();
            System.out.println("Total " + String.valueOf(tagNum) + " tags");
            int tagLengthInLong = dis.readInt();
            int taxaNum = dis.readInt();
            long[] tag;
            String[] taxa = new String[taxaNum];
            for (int i = 0; i < taxaNum; i++) taxa[i] = dis.readUTF();
            for (int i = 0; i < tagNum; i++) {
                tag = new long[tagLengthInLong];
                for (int j = 0; j < tag.length; j++) {
                    tag[j] = dis.readLong();
                }
                dis.readByte();
                int index = Arrays.binarySearch(tags, BaseEncoder.getSequenceFromLong(tag));
                long dep = 0;
                if (index < 0) {
                    for (int j = 0; j < taxaNum; j++) dis.readByte();
                }
                else {
                    for (int j = 0; j < taxaNum; j++) {
                        int count = dis.readByte();
                        if (count > 0) {
                            dep+=depthMap.get(taxa[j]);
                        }
                    }
                    bw.write(tags[index]+"\t"+String.valueOf(dep));
                    bw.newLine();
                }
                if (i%10000 ==0) {
                    bw.flush();
                    System.out.println(i);
                }
            }
            bw.flush();
            bw.close();
            dis.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkDepthFile () {
        String depthDirS = "E:\\Research\\geneticMapping\\accuracyFactor\\depth\\";
        String depthFileS = "E:\\Research\\geneticMapping\\accuracyFactor\\taxaDepth.txt";
        File[] fs = new File(depthDirS).listFiles();
        TreeSet<String> ts = new TreeSet();
        
        HashMap<String, Integer> map = new HashMap();
        for (int i = 0; i < fs.length; i++) {
            Table t = new Table(fs[i].getAbsolutePath());
            for (int j = 0; j < t.getRowNumber(); j++) {
                ts.add(t.content[j][0]);
                
                map.put(t.content[j][0], Integer.valueOf(t.content[j][1]));
            }
        }
        int[] nameCount = new int[ts.size()];
        String[] uniqueName = ts.toArray(new String[ts.size()]);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(depthFileS);
            bw.write("Taxa\tTagCount");
            bw.newLine();
            for (int i = 0; i < uniqueName.length; i++) {
                bw.write(uniqueName[i]+"\t"+map.get(uniqueName[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkDepthOfTBTs () {
        String tbtDirS = "/SSD/mingh/tbt/";
        String depthDirS = "/SSD/mingh/depth/";
        
        new File(depthDirS).mkdir();
        File[] fs = new File(tbtDirS).listFiles();
        int batchSize = 20;
        int batchNum = 0;
        int[] batchStart = null;
        int[] batchEnd = null;
        int left = fs.length % batchSize;
        if (left == 0) {
            batchNum = fs.length/batchSize;
            batchStart = new int[batchNum];
            batchEnd = new int[batchNum];
            for (int i = 0; i < batchNum; i++) {
                batchStart[i] = i*batchSize;
                batchEnd[i] = batchStart[i] + batchSize;
            }
        }
        else {
            batchNum = fs.length/batchSize + 1;
            batchStart = new int[batchNum];
            batchEnd = new int[batchNum];
            for (int i = 0; i < batchNum; i++) {
                batchStart[i] = i*batchSize;
                batchEnd[i] = batchStart[i] + batchSize;
            }
            batchEnd[batchNum-1] = batchStart[batchNum-1] + left;
        }
        for (int i = 0; i < batchNum; i++) {
            int jobNum = batchEnd[i] - batchStart[i];
            Task[] jobs = new Task[jobNum];
            Thread[] mts = new Thread[jobNum];
            for (int j = 0; j < jobNum; j++) {
                String outfileS = fs[batchStart[i]+j].getName().replaceFirst(".tbt.byte", ".txt");
                outfileS = depthDirS+outfileS;
                jobs[j] = new Task (fs[batchStart[i]+j].getAbsolutePath(), outfileS);
            }
            for (int j = 0; j < jobNum; j++) {
                mts[j] = new Thread(jobs[j]);
                mts[j].start();
            }
            for (int j = 0; j < jobNum; j++) {
                try {
                    mts[j].join();
                }
                catch (Exception e) {
                    System.out.println(e.toString());
                }
            }
        }
    }
    
    class Task implements Runnable {
        String infileS = null;
        String outfileS = null;
        public Task (String infileS, String outfileS) {
            this.infileS = infileS;
            this.outfileS = outfileS;
        }
        @Override
        public void run() {
            TagsByTaxaByte tbt = new TagsByTaxaByte(infileS, FilePacking.Byte);
            String[] taxa = tbt.getTaxaNames();
            int[] taxaTagCnt = new int[taxa.length];
            for (int i = 0; i < tbt.getTagCount(); i++) {
                for (int j = 0; j < taxa.length; j++) {
                    taxaTagCnt[j] += tbt.getReadCountForTagTaxon(i, j);
                }
            }
            try {
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                bw.write("Taxa\tTagCount");
                bw.newLine();
                for (int i = 0; i < taxa.length; i++) {
                    bw.write(taxa[i]+"\t"+String.valueOf(taxaTagCnt[i]));
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println(infileS + " processed");
        }
        
    }
    public void compareImputation () {
        //String unimputatedFileS = "M:\\pav\\anchor\\h5\\bpec.hmp.h5";
        //String imputatedFileS = "M:\\pav\\anchor\\h5\\impAll.hmp.h5";
        //String comparisonFileS = "E:\\Research\\geneticMapping\\imputation\\imputationCompare.txt";
        
        String unimputatedFileS = "/SSD/mingh/h5/bpec.hmp.h5";
        String imputatedFileS = "/SSD/mingh/h5/impAll.hmp.h5";
        String comparisonFileS = "/SSD/mingh/imputationCompare.txt";
        
        GenotypeTable imGt = ImportUtils.readGuessFormat(imputatedFileS);
        GenotypeTable unGt = ImportUtils.readGuessFormat(unimputatedFileS);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(comparisonFileS);
            bw.write("Chromosome\tSite\tUnimputatedMissing\tUnimputatedMissingRate\tImputaedMissing\tImputaedMissingRate\tMissingDecrease\tMissingDecreaseRate");
            bw.newLine();
            for (int i = 0; i < unGt.numberOfSites(); i++) {
                bw.write(unGt.chromosome(i).getName()+"\t"+String.valueOf(unGt.positions().chromosomalPosition(i))+"\t");
                BitSet ma = unGt.allelePresenceForAllTaxa(i, WHICH_ALLELE.Major);
                BitSet mi = unGt.allelePresenceForAllTaxa(i, WHICH_ALLELE.Minor);
                OpenBitSet oma = new OpenBitSet(ma.getBits(),ma.getNumWords());
                OpenBitSet omi = new OpenBitSet(mi.getBits(),mi.getNumWords());
                oma.or(omi);
                int unMissing = unGt.numberOfTaxa()-(int)oma.cardinality();
                bw.write(String.valueOf(unMissing)+"\t"+String.valueOf((double)unMissing/unGt.numberOfTaxa())+"\t");
                ma = imGt.allelePresenceForAllTaxa(i, WHICH_ALLELE.Major);
                mi = imGt.allelePresenceForAllTaxa(i, WHICH_ALLELE.Minor);
                oma = new OpenBitSet(ma.getBits(),ma.getNumWords());
                omi = new OpenBitSet(mi.getBits(),mi.getNumWords());
                oma.or(omi);
                int imMissing = imGt.numberOfTaxa()-(int)oma.cardinality();
                bw.write(String.valueOf(imMissing)+"\t"+String.valueOf((double)imMissing/imGt.numberOfTaxa())+"\t");
                bw.write(String.valueOf(unMissing-imMissing)+"\t"+String.valueOf((double)(unMissing-imMissing)/unGt.numberOfTaxa()));
                bw.newLine();
                if (i%1000 ==0 ) System.out.println(String.valueOf(i)+" sites processed");
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
