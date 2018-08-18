/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package analysis.geneticMapping;

import graphcis.r.BoxPlot;
import graphcis.r.DensityPlot;
import graphcis.r.ScatterPlot;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import utils.IoUtils;

/**
 *
 * @author Fei Lu <fl262@cornell.edu>
 */
public class MultipleML {
    
    public MultipleML () {
        this.drawScatterAndBox();
    }
    
    public void drawScatterAndBox () {
        String infileDirS = "E:\\Research\\geneticMapping\\multiML\\trainResult\\";
        String scatterDirS = "E:\\Research\\geneticMapping\\multiML\\scatterPlots\\";
        String boxPlotFileS = "E:\\Research\\geneticMapping\\multiML\\errorBox.pdf";
        String correlationFile = "E:\\Research\\geneticMapping\\multiML\\correlation.txt";
        
        File[] fs = new File(infileDirS).listFiles();
        String[] methodName = new String[fs.length];
        double[][] prediction = new double[fs.length][];
        double[][] observation = new double[fs.length][];
        double[][] error = new double[fs.length][];
        double[] correlation = new double[fs.length];
        double[] meanError = new double[fs.length];
        for (int i = 0; i < fs.length; i++) {
            methodName[i] = fs[i].getName().replaceFirst(".txt", "");
            try {
                BufferedReader br = IoUtils.getTextReader(fs[i].getAbsolutePath());
                String temp;
                while (!(temp = br.readLine()).startsWith("inst#")) {}
                int cnt = 0;
                while (!(temp = br.readLine()).equals("")) {cnt++;}
                prediction[i] = new double[cnt];
                observation[i] = new double[cnt];
                error[i] = new double[cnt];
                br = IoUtils.getTextReader(fs[i].getAbsolutePath());
                while (!(temp = br.readLine()).startsWith("inst#")) {}
                cnt = 0;
                while (!(temp = br.readLine()).equals("")) {
                    String[] tem = temp.split("\\s+");
                    prediction[i][cnt] = Double.valueOf(tem[3]);
                    observation[i][cnt] = Double.valueOf(tem[2]);
                    error[i][cnt] = Math.abs(Double.valueOf(tem[4]));
                    cnt++;
                }
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("Correlation coefficient")) {
                        correlation[i] = Double.valueOf(temp.split("\\s+")[2]);
                    }
                    else if (temp.startsWith("Mean absolute error")) {
                        meanError[i] = Double.valueOf(temp.split("\\s+")[3]);
                    }
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        for (int i = 0; i < fs.length; i++) {
            String outfileS = scatterDirS+methodName[i]+".pdf";
            int sampleSize = 5000;
            double[] obser = new double[sampleSize];
            double[] predi = new double[sampleSize];
            for (int j = 0; j < sampleSize; j++) {
                obser[j] = observation[i][j];
                predi[j] = prediction[i][j];
            }
            ScatterPlot plot = new ScatterPlot(obser, predi);
            plot.setColor(255, 0, 0, 60);
            plot.setPlottingCharacter(16);
            plot.setSlideMode();
            plot.setXLim(1, 10);
            plot.setYLim(1, 10);
            plot.setXLab("Log10 (observed distance)");
            plot.setYLab("Log10 (predicted distance)");
            plot.setTitle(methodName[i]);
            plot.saveGraph(outfileS);
        }
        for (int i = 0; i < methodName.length - 1; i++) {
            for (int j = i; j < methodName.length; j++) {
                if (correlation[i] < correlation[j]) {
                    String temp = methodName[i];
                    methodName[i] = methodName[j];
                    methodName[j] = temp;
                    double td = correlation[i];
                    correlation[i] = correlation[j];
                    correlation[j] = td;
                    td = meanError[i];
                    meanError[i] = meanError[j];
                    meanError[j] = td;
                    double[] te = error[i];
                    error[i] = error[j];
                    error[j] = te;
                }
            }
        }
        BoxPlot b = new BoxPlot(error, methodName);
        b.setColorRainbow();
        b.setXLabVertical();
        b.setYLab("Absolute error of prediction");
        b.setTitle("Error comparison of different models");
        b.saveGraph(boxPlotFileS);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(correlationFile);
            bw.write("MethodName\tCorrelationCoefficient\tMeanError\tP-value");
            bw.newLine();
            for (int i = 0; i < methodName.length; i++) {
                double[][] data = new double[prediction[i].length][2];
                for (int j = 0; j < data.length; j++) {
                    data[j][0] = prediction[i][j];
                    data[j][1] = observation[i][j];
                }
                PearsonsCorrelation co= new PearsonsCorrelation(data);
                double pValue = co.getCorrelationPValues().getEntry(0, 1);
                bw.write(methodName[i]+"\t"+String.valueOf(correlation[i])+"\t"+String.valueOf(meanError[i])+"\t"+String.valueOf(pValue));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch(Exception e) {
            e.printStackTrace();
        }
        
    }
    
}
