/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.cassava.wgs;

import graphcis.HeatMap;
import java.io.BufferedWriter;
import java.io.File;
import net.maizegenetics.dna.read.FastqChunk;
import net.maizegenetics.dna.read.ReadUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class SequencingQualityPattern {
    
    public SequencingQualityPattern () {
        //this.getBaseMeanAndSD();
        this.mkErrorHeatMap();
    }
    
    public void mkErrorHeatMap () {
        String infileDirS = "M:\\pipelineTest\\cassava\\wgs\\pe\\fastqChunk\\";
        String heatMapDirS = "M:\\pipelineTest\\cassava\\wgs\\sequencingErrorPattern\\heatMap\\";
        File[] fs = new File(infileDirS).listFiles();
        int readNumber = 300;
        for (int i = 0; i < fs.length; i++) {
            FastqChunk fq = new FastqChunk(fs[i].getAbsolutePath(), ReadUtils.ReadFormat.FastqGzip);
            if (readNumber > fq.getReadNum()) readNumber = fq.getReadNum();
            double[][] matrix = new double[readNumber][fq.getRead(0).getReadLength()];
            for (int j = 0; j < matrix.length; j++) {
                for (int k = 0; k < matrix[0].length; k++) {
                    matrix[j][k] = fq.getRead(j).getBaseQuality(k, fq.getPhredScale());
                }
            }
            String outfileS = new File(heatMapDirS, fs[i].getName().replaceFirst("fastq.gz", "pdf")).getAbsolutePath();
            String[] rowName = new String[matrix.length];
            String[] columnName = new String[matrix[0].length];
            for (int j = 0; j < rowName.length; j++) {
                rowName[j] = "Read"+String.valueOf(j+1);
            }
            for (int j = 0; j < columnName.length; j++) {
                columnName[j] = "Base"+String.valueOf(j+1);
            }
            HeatMap hm = new HeatMap(matrix, rowName, columnName);
            hm.setCellSize(10, 10);
            hm.setCapValues(0, 41);
            hm.saveGraph(outfileS);
        }
    }
    
    public void getBaseMeanAndSD () {
        String infileDirS = "M:\\pipelineTest\\cassava\\wgs\\pe\\fastqChunk\\";
        String summaryFileS = "M:\\pipelineTest\\cassava\\wgs\\sequencingErrorPattern\\meanSd.summary.txt";
        File[] fs = new File(infileDirS).listFiles();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(summaryFileS);
            bw.write("Filename\tPair\tReadlength\tPlatform\tStatistics");
            for (int i = 0; i < 201; i++) {
                bw.write("\t"+String.valueOf(i+1)+"th");
            }
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                FastqChunk fq = new FastqChunk(fs[i].getAbsolutePath(), ReadUtils.ReadFormat.FastqGzip);
                int phredScale = fq.getPhredScale();
                int readLength = fq.getRead(0).getReadLength();
                String platform = "NextSeq";
                if (readLength > 170) platform = "HiSeq";
                StringBuilder averS = new StringBuilder();
                StringBuilder sdS = new StringBuilder();
                for (int j = 0; j < readLength; j++) {
                    double[] baseQual = new double[fq.getReadNum()];
                    for (int k = 0; k < baseQual.length; k++) {
                        baseQual[k] = fq.getRead(k).getBaseQuality(j, phredScale);
                    }
                    DescriptiveStatistics ds = new DescriptiveStatistics(baseQual);
                    double aver = ds.getMean();
                    double sd = ds.getStandardDeviation();
                    averS.append("\t").append(aver);
                    sdS.append("\t").append(sd);
                }
                String pair;
                String[] temp = fs[i].getName().split("_");
                if (temp[temp.length-1].split("\\.")[0].equals("R1")) pair = "R1";
                else pair = "R2";
                bw.write(fs[i].getAbsolutePath()+"\t"+pair+"\t"+String.valueOf(readLength)+"\t"+platform+"\tMean"+averS.toString());
                bw.newLine();
                bw.write(fs[i].getAbsolutePath()+"\t"+pair+"\t"+String.valueOf(readLength)+"\t"+platform+"\tSD"+sdS.toString());
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
