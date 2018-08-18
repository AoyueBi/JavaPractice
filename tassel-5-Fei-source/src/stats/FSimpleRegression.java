/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package stats;

import org.apache.commons.math3.stat.regression.SimpleRegression;

/**
 *
 * @author fl262
 */
public class FSimpleRegression {
    SimpleRegression sr = null;
    double[] xs = null;
    double[] ys = null;
    
    public FSimpleRegression (double[] xs, double[] ys) {
        this.initialize(xs, ys);
    }
    
    private void initialize (double[] xs, double[] ys) {
        this.xs = xs;
        this.ys = ys;
        sr = new SimpleRegression();
        for (int i = 0; i < xs.length; i++) {
            sr.addData(xs[i], ys[i]);
        }
    }
    
    public double getSlope () {
        return sr.getSlope();
    }
    
    public double getIntercept () {
        return sr.getIntercept();
    }
    
    public double getPrediction (double x) {
        return sr.predict(x);
    }
    
    public double[] getPrediciton () {
        double[] yba = new double[ys.length];
        for (int i = 0; i < yba.length; i++) {
            yba[i] = this.getPrediction(xs[i]);
        }
        return yba;
    }
    
    public double[] getResidual () {
        double[] res = new double[ys.length];
        for (int i = 0; i < res.length; i++) {
            res[i] = ys[i]-this.getPrediction(xs[i]);
        }
        return res;
    }
}
