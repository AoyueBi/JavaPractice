/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package utils.stats.r;

import graphcis.r.Rgraphics;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import rcaller.RCaller;
import rcaller.RCode;
import utils.IoUtils;

/**
 * PCA analysis using R package. R result output to temp directory, read into java. Delete r output
 * Number of observation should be equal or greater than number of variables
 * @author Fei Lu <fl262@cornell.edu>
 */
public class PCA {
    String rPath = Rgraphics.RPath;
    File tempDir = new File (System.getProperty("java.io.tmpdir"));
    String varianceFileS = new File(tempDir, "var.txt").getAbsolutePath().replaceAll("\\\\", "/");
    String scoreFileS = new File(tempDir, "score.txt").getAbsolutePath().replaceAll("\\\\", "/");
    String valueFileS = new File(tempDir, "value.txt").getAbsolutePath().replaceAll("\\\\", "/");
    //String varianceFileS = "M:/var.txt";
    //String scoreFileS = "M:/score.txt";
    /**
     *          Taxa1   Taxa2   Taxa3   Taxa4
     * Site1
     * Site2
     * Site3
     * 1st dimension is column number (observation), 2nd dimension is row number (variables). 
     * Calculating PCs of 2nd dimension. This transformed matrix is easier to pass to R
     */
    double[][] value;
    double[] exVar;
    double[] var;
    double[] sd;
    double[][] score;
    
    public PCA () {}
    public PCA (double[][] value) {
        this.value = value;
        this.excute();
        this.readIn();
    }
    
    public int getNumOfPC () {
        return sd.length;
    }
    
    public double getSDOfPC (int index) {
        return sd[index];
    }
    
    public double getVarianceExplainedOfPC (int index) {
        return exVar[index];
    }
    
    public double getEigenvalueOfPC (int index) {
        return var[index];
    }
    
    public double[] getScoresOfPC (int index) {
        return score[index];
    }
    
    public double getScoreOfPC (int pcIndex, int observationIndex) {
        return score[pcIndex][observationIndex];
    }
    
    public String getSummary () {
       StringBuilder sb = new StringBuilder();
       sb.append("PCs\tSD\tVariance\tVariaceExplained\n");
       for (int i = 0; i < sd.length; i++) {
           sb.append("PC").append(String.valueOf(i+1)).append("\t").append(sd[i]).append("\t").append(var[i]).append("\t").append(exVar[i]).append("\n");
       }
       return sb.toString();
    }
    
    private void readIn () {
       try {
           BufferedReader br = IoUtils.getTextReader(varianceFileS);
           String temp;
           int cnt = 0;
           while ((temp = br.readLine()) != null) cnt++;
           this.sd = new double[cnt];
           this.var = new double[cnt];
           this.exVar = new double[cnt];
           cnt = 0;
           double sum = 0;
           br = IoUtils.getTextReader(varianceFileS);
           while ((temp = br.readLine()) != null) {
               sd[cnt] = Double.valueOf(temp);
               var[cnt] = Math.pow(sd[cnt], 2);
               sum+=var[cnt];
               cnt++;
           }
           for (int i = 0; i < sd.length; i++) exVar[i] = var[i]/sum;
           br.close();
           score = new double[sd.length][value[0].length];
           br = IoUtils.getTextReader(scoreFileS);
           for (int i = 0; i < sd.length; i++) {
               for (int j = 0; j < value[0].length; j++) {
                   temp = br.readLine();
                   score[i][j] = Double.valueOf(temp);
               }
           }
           br.close();
           new File(this.varianceFileS).delete();
           new File(this.scoreFileS).delete();
           new File(this.valueFileS).delete();
       }
       catch (Exception e) {
           e.printStackTrace();
       }
    }
    
    private void excute () {
        try {
            BufferedWriter bw = IoUtils.getTextWriter(this.valueFileS);
            for (int i = 0; i < value.length-1; i++) {
                bw.write("v"+String.valueOf(i+1)+"\t");
            }
            bw.write("v"+String.valueOf(value.length));
            bw.newLine();
            for (int i = 0; i < value[0].length; i++) {
                for (int j = 0; j < value.length-1; j++) {
                    bw.write(String.valueOf(value[j][i])+"\t");
                }
                bw.write(String.valueOf(value[value.length-1][i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("PCA in progress...\n");
        RCaller caller = new RCaller();
        caller.setRscriptExecutable(this.rPath);
        RCode rCode = new RCode();
        rCode.addRCode("library(psych)");
        rCode.addRCode("df <- read.table(\""+this.valueFileS+"\", header = T, sep = \"\t\")");
        rCode.addRCode("pca<-princomp(df)");
        rCode.addRCode("score <- pca$scores");
        rCode.addRCode("sd <- pca$sdev");
        rCode.addRCode("write (sd, \""+this.varianceFileS+"\", ncolumns = 1)");
        rCode.addRCode("write (score, \""+this.scoreFileS+"\", ncolumns = 1)");
        caller.setRCode(rCode);
        //System.out.println(rCode.getCode());
        caller.runOnly();
        System.out.println("PCA finished");
    }
}
