/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maf;

import com.itextpdf.awt.DefaultFontMapper;
import com.itextpdf.awt.PdfGraphics2D;
import com.itextpdf.text.Document;
import com.itextpdf.text.Rectangle;
import com.itextpdf.text.pdf.PdfContentByte;
import com.itextpdf.text.pdf.PdfWriter;
import format.Table;
import graphcis.r.ScatterPlot;
import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.util.Arrays;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class MafGo {
    
    public MafGo () {
        //this.annotationPipe();
        this.analysisPipe();
        
    }
    
    private void analysisPipe () {
        //new VariantDiscovery ();
        //new Recombination();
        //new Expression();
        //new Mutants();
        
        //new Pathway();//deprecated
        //new Pollen();
        //new PopGenGroup();
        //new PopGenParameters ();
        //new PopGenAdaptation();
        //new PopGenSelection();
        //new TraitAssociationOld ();
    
        
        //new NAMTraitAssociation();
    }
    
    public void annotationPipe () {
        //new SiteUniqueness();
        //new MafAttributeHmp3();
        //new MafAnnotation();
        //new VCAP();
        //new AllelePresence ();
        //new SiftScore();
        //new Gerp();
    }
    
    public static void main (String[] args) {
        new MafGo();
    }
}
