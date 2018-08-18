/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.maf;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class WGSDataCheck {
    
    public WGSDataCheck () {
        this.checkFormat();
    }
    
    public void checkFormat () {
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
