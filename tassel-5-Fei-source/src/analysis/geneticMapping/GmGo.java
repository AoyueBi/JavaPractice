/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package analysis.geneticMapping;

import graphcis.r.DensityPlot;
import graphcis.r.Histogram;
import graphcis.r.DensityPlotMultiClass;
import java.util.ArrayList;
import net.maizegenetics.analysis.gbs.pana.PanAH5ToAnchorPlugin;

/**
 *
 * @author Fei Lu 
 */
public class GmGo {
    
    public GmGo () {
        //this.mappingPipe();
        //this.paperPipe();
        
    }
    
    public void paperPipe () {
        new GeneticMappingPaper ();
    }
    
    public void mappingPipe () {
        //this.convertTBTByte();
        //this.mergeHapMap();
        //this.convertHapMapToH5();    
        //this.convertToSBit();
    }
    
    public void convertHapMapToH5 () {
        //String hapMapFileS = "M:\\pav\\anchor\\testData\\impAll.hmp.txt";
        //String h5FileS = "M:\\pav\\anchor\\testData\\impAll.hmp.h5";
        //String hapMapFileS = "/SSD/mingh/impAll.hmp.txt";
        //String h5FileS = "/SSD/mingh/impAll.hmp.h5";
        String hapMapFileS = "/SSD/mingh/bpec.hmp.txt";
        String h5FileS = "/SSD/mingh/bpec.hmp.h5";
        new GMUtils().convertHapMapToH5(hapMapFileS, h5FileS);
    }
    
    public void mergeHapMap () {
        //String inputDirS = "M:\\pav\\anchor\\testData\\impAll\\name\\";
        //String outputHapMap = "M:\\pav\\anchor\\testData\\impAll.hmp.txt";
        //String inputDirS = "/SSD/mingh/impAll/";
        //String outputHapMap = "/SSD/mingh/impAll.hmp.txt";
        String inputDirS = "/SSD/mingh/bpec/";
        String outputHapMap = "/SSD/mingh/bpec.hmp.txt";
        new GMUtils().mergeHapMap(inputDirS, outputHapMap);
    }
    
    public void convertTBTByte () {
        //String sourceTBT = "M:\\pav\\mergedTBT\\sample\\tbt_test_50.byte";
        //String targetTBT = "M:\\pav\\mergedTBT\\h5\\UB73.tagGroup.tbt.h5";
        String sourceTBT = "/SSD/mingh/UB73.tbt.byte";
        String targetTBT = "/SSD/mingh/UB73.tagGroup.tbt.h5";
        new GMUtils().convertTBTByte(sourceTBT, targetTBT);
    }
    
    public static void main (String[] args) {
        new GmGo ();
    }
    
}
