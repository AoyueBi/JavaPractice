/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.gbs.cnv;

import java.util.ArrayList;

/**
 *
 * @author fl262
 */
public class CNVGo {
    
    public CNVGo () {
        this.test();
    }
    
    public void test () {
        String inputFileS = "M:\\pav\\cnv\\bbt\\cnv.ataxa.200k.bbt.txt";
        String ratioFileS = "M:\\production\\cnv\\bbt\\cnv.ataxa.200k.bbt.ratio.txt";
        String pValueFileS = "M:\\production\\cnv\\bbt\\cnv.ataxa.200k.bbt.pValue.txt";
        BinByTaxa bbt = new BinByTaxa(inputFileS, "b73");
        ArrayList<CNVInfo> cnvList = bbt.scanCNV("MO17:81546ABXX:8:B6", 2, -5, 3);
        System.out.print("ok");
        //bbt.writeRatioFile(ratioFileS);
    }
    
    public static void main (String[] args) {
        new CNVGo();
    }
}
