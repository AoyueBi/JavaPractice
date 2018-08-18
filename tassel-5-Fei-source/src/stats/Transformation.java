/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package stats;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 *
 * @author fl262
 */
public class Transformation {
    
    public static double[] getZScore (double[] x) {
        DescriptiveStatistics d = new DescriptiveStatistics(x);
        double mean = d.getMean();
        double sd = d.getStandardDeviation();
        double[] z = new double[x.length];
        for (int i = 0; i < z.length; i++) {
            z[i] = (x[i]-mean)/sd;
        }
        return z;
    }
    
}
