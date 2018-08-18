/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.geneticMapping;

import format.Table;
import graphcis.HeatMap;
import graphcis.r.ScatterPlot;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.util.Arrays;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import utils.IoUtils;
import utils.stats.r.PCA;

/**
 *
 * @author fl262
 */
public class PopulationStructure {
    
    public PopulationStructure () {
        //this.mkMatrix();
        //this.mkPCs();
        //this.mkTagAndPCCorrelation();
        this.mkAccuracyPCR2Figure();
        //this.mkPCAndAttributesCorrelation();
        //this.mkHeatMap();
    }
    
    public void mkHeatMap () {
        String pcAttR2FileS = "E:\\Research\\geneticMapping\\populationStructure\\pcAttR2_final.txt";
        String pcHeatMap = "E:\\Research\\geneticMapping\\populationStructure\\pcAttR2_final.pdf";
        Table t = new Table (pcAttR2FileS);
        HeatMap hm = new HeatMap(t);
        hm.saveGraph(pcHeatMap);
    }
    
    public void mkPCAndAttributesCorrelation () {
        String trainArffFileS = "E:\\Research\\geneticMapping\\populationStructure\\box_Num_G.arff.csv";
        String tagPCRFileS = "E:\\Research\\geneticMapping\\populationStructure\\tagPcR.txt";
        String pcAttR2FileS = "E:\\Research\\geneticMapping\\populationStructure\\pcAttR2.txt";
        Table t = new Table(tagPCRFileS);
        double[][] pcr2 = new double[t.getColumnNumber()-1][t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            for (int j = 0; j < pcr2.length; j++) {
                //pcr2[j][i] = Math.pow(Double.valueOf(t.content[i][j+1]),2);
                pcr2[j][i] = Double.valueOf(t.content[i][j+1]);
            }
        }
        t = new Table (trainArffFileS, ",");
        double[][] atts = new double[t.getColumnNumber()][t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            for (int j = 0; j < t.getColumnNumber(); j++) {
                atts[j][i] = Double.valueOf(t.content[i][j]);
            }
        }
        double[][] correlations = new double[atts.length][pcr2.length];
        for (int i = 0; i < atts.length; i++) {
            for (int j = 0; j < pcr2.length; j++) {
                correlations[i][j] = Math.pow(new PearsonsCorrelation().correlation(atts[i], pcr2[j]),2);
                //correlations[i][j] = new PearsonsCorrelation().correlation(atts[i], pcr2[j]);
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(pcAttR2FileS);
            bw.write("Attributes");
            for (int i = 0; i < pcr2.length; i++) bw.write("\tPC"+String.valueOf(i+1));
            bw.newLine();
            for (int i = 0; i < correlations.length; i++) {
                bw.write(t.header[i]);
                for (int j = 0; j < correlations[0].length; j++) {
                    bw.write("\t"+String.valueOf(correlations[i][j]));
                }
                bw.newLine();
            }
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkAccuracyPCR2Figure () {
        String tagPCRFileS = "E:\\Research\\geneticMapping\\populationStructure\\tagPcR.txt";
        String UB73TrainFileS = "E:\\Research\\geneticMapping\\populationStructure\\box_Num_G.arff.csv";
        String figureDirS = "E:\\Research\\geneticMapping\\populationStructure\\";
        Table t = new Table (tagPCRFileS);
        int size = 5000;
        double[][] pcr2 = new double[3][size];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < pcr2.length; j++) {
                //pcr2[j][i] = Math.pow(Double.valueOf(t.content[i][j+1]), 2);
                pcr2[j][i] = Double.valueOf(t.content[i][j+1]);
            }
        }
        t= new Table(UB73TrainFileS, ",");
        double[] d = new double[size];
        for (int i = 0; i < d.length; i++) {
            d[i] = Double.valueOf(t.content[i][9]);
        }
        for (int i = 0; i < pcr2.length; i++) {
            String outfileS = new File(figureDirS, "pc"+String.valueOf(i+1)+".pdf").getAbsolutePath();
            ScatterPlot s = new ScatterPlot(d, pcr2[i]);
            s.setColor(255, 0, 0, 60);
            s.setPlottingCharacter(16);
            s.setXLab("Log10 (distance)");
            s.setYLab("Pearson's r with PC" +String.valueOf(i+1));
            s.setTitle("");
            s.setSlideMode();
            s.saveGraph(outfileS);
            double[][] data = new double[d.length][2];
            for (int j = 0; j < data.length; j++) {
                data[j][0]  = d[j];
                data[j][1] = pcr2[i][j];
            }
            PearsonsCorrelation p = new PearsonsCorrelation(data);
            double r = Math.pow(p.getCorrelationMatrix().getEntry(0, 1),2);
            double pv = 1;
            try {
                pv = p.getCorrelationPValues().getEntry(0, 1);
            }
            catch (Exception e) {e.printStackTrace();}
       
            System.out.println(r + "\t" +pv);
        }
        
    }
    
    public void mkTagAndPCCorrelation () {
        String UB73TrainFileS = "M:\\pav\\train\\train_result\\UB73train.mpgmap.txt";
        String tbtFileS = "M:\\pav\\mergedTBT\\UB73_withIssue.tbt.byte";
        String pcScoreFileS = "E:\\Research\\geneticMapping\\populationStructure\\pcScore.txt";
        String tagPCR2FileS = "E:\\Research\\geneticMapping\\populationStructure\\tagPcR.txt";
        Table t = new Table (UB73TrainFileS);
        String[] seqs = new String[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            seqs[i] = t.content[i][0];
        }
        t = new Table(pcScoreFileS);
        String[] pcTaxa = new String[t.getRowNumber()];
        int[] taxaIndex = new int[pcTaxa.length];
        double[][] pcs = new double[t.getColumnNumber()-1][t.getRowNumber()];
        double[][] r2 = new double[seqs.length][pcs.length];
        for (int i = 0; i < t.getRowNumber(); i++) {
            pcTaxa[i] = t.content[i][0];
            for (int j = 0; j < pcs.length; j++) {
                pcs[j][i] = Double.valueOf(t.content[i][j+1]);
            }
        }
        String[] tbtTaxa = null;
        try {
            DataInputStream dis = IoUtils.getBinaryReader(tbtFileS);
            BufferedWriter bw = IoUtils.getTextWriter(tagPCR2FileS);
            bw.write("Tag\tRPC1\tRPC2\tRPC3");
            bw.newLine();
            int tagNum = dis.readInt();
            int tagLengthInLong = dis.readInt();
            int tbtTaxaNum = dis.readInt();
            tbtTaxa = new String[tbtTaxaNum];
            for (int i = 0; i < tbtTaxa.length; i++) {
                tbtTaxa[i] = dis.readUTF();
            }
            System.out.println("Redicting...");
            for (int i = 0; i < pcTaxa.length; i++) {
                for (int j = 0; j < tbtTaxa.length; j++) {
                    if (pcTaxa[i].equals(tbtTaxa[j])) {
                        taxaIndex[i] = j;
                        break;
                    }
                }
            }
            System.out.println("Redirection done");
            int cnt = 0;
            for (int i = 0; i < tagNum; i++) {
                long[] tag = new long[tagLengthInLong];
                for (int j = 0; j < tagLengthInLong; j++) {
                    tag[j] = dis.readLong();
                }
                dis.readByte();
                String currentSeq = BaseEncoder.getSequenceFromLong(tag);
                if (currentSeq.equals(seqs[cnt])) {
                    double[] value = new double[pcTaxa.length];
                    byte[] geno = new byte[tbtTaxa.length];
                    for (int j = 0; j < geno.length; j++) {
                        geno[j] = dis.readByte();
                    }
                    for (int j = 0; j < value.length; j++) {
                        if (geno[taxaIndex[j]] == 0) value[j] = 0;
                        else value[j] = 1;
                    }
                    for (int j = 0; j < r2[cnt].length; j++) {
                        PearsonsCorrelation p = new PearsonsCorrelation();
                        r2[cnt][j] = p.correlation(pcs[j], value);
                    }
                    bw.write(currentSeq);
                    for (int j = 0; j < r2[cnt].length; j++) {
                        bw.write("\t"+r2[cnt][j]);
                    }
                    bw.newLine();
                    cnt++;
                }
                else {
                    for (int j = 0; j < tbtTaxa.length; j++) dis.readByte();
                }
                if (cnt == seqs.length) break;
            }
            bw.flush();
            bw.close();
            dis.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkPCs() {
        String matrixFileS = "E:\\Research\\geneticMapping\\populationStructure\\matrix.txt";
        String pcFileS = "E:\\Research\\geneticMapping\\populationStructure\\pcScore.txt";
        String pcSummary = "E:\\Research\\geneticMapping\\populationStructure\\pcSummary.txt";
        //String matrixFileS = "E:\\Research\\geneticMapping\\populationStructure\\matrix_missingToMajor.txt";
        //String pcFileS = "E:\\Research\\geneticMapping\\populationStructure\\pcScore_missingToMajor.txt";
        //String pcSummary = "E:\\Research\\geneticMapping\\populationStructure\\pcSummary_missingToMajor.txt";
        Table t =new Table(matrixFileS);
        double[][] matrix = new double[t.getRowNumber()][t.getColumnNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            for (int j = 0; j < t.getColumnNumber(); j++) {
                matrix[i][j] = Double.valueOf(t.content[i][j]);
            }
        }
        PCA p = new PCA(matrix);
        double[][] pcs = new double[3][];
        for (int i = 0; i < pcs.length; i++) {
            pcs[i] = p.getScoresOfPC(i);
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(pcFileS);
            bw.write("Taxa\tPC1\tPC2\tPC3");
            bw.newLine();
            for (int i = 0; i < t.getColumnNumber(); i++) {
                bw.write(t.header[i]+"\t"+String.valueOf(pcs[0][i])+"\t"+String.valueOf(pcs[1][i])+"\t"+String.valueOf(pcs[2][i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            bw = IoUtils.getTextWriter(pcSummary);
            bw.write(p.getSummary());
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkMatrix () {
        String hapmapFileS = "M:\\pav\\anchor\\h5\\impAll.hmp.h5";
        String matrixFileS = "E:\\Research\\geneticMapping\\populationStructure\\matrix.txt";
        //String matrixFileS = "E:\\Research\\geneticMapping\\populationStructure\\matrix_missingToMajor.txt";
        
        int snpNum = 2000;
        int[] snpIndex = new int[snpNum];
        GenotypeTable gt = ImportUtils.readGuessFormat(hapmapFileS);
        for (int i = 0; i < snpNum; i++) {
            snpIndex[i] = (int)(gt.numberOfSites()*Math.random());
        }
        double[][] matrix = new double[snpNum][gt.numberOfTaxa()];
        for (int i = 0; i < snpNum; i++) {
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
        try {
            BufferedWriter bw = IoUtils.getTextWriter(matrixFileS);
            for (int i = 0; i < gt.numberOfTaxa()-1; i++) {
                bw.write(gt.taxaName(i)+"\t");
            }
            bw.write(gt.taxaName(gt.numberOfTaxa()-1));
            bw.newLine();
            for (int i = 0; i < snpNum; i++) {
                for (int j = 0; j < matrix[i].length-1; j++) {
                    bw.write(String.valueOf(matrix[i][j])+"\t");
                }
                bw.write(String.valueOf(matrix[i][matrix[i].length-1])+"\t");
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
