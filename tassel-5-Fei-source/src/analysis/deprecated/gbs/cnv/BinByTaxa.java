/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.gbs.cnv;

import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.TreeSet;
import org.apache.commons.math3.stat.inference.TestUtils;
import utils.IoUtils;

/**
 *
 * @author Fei Lu
 */
public class BinByTaxa {
    String[] taxa = null;
    int[] chrID = null;
    int[][] binStart = null;
    int[][][] depth = null;
    int refTaxonIndex = Integer.MIN_VALUE;
    float[][][] ratio = null;
    float[][][] pValue = null;
    
    public BinByTaxa (String binByTaxaFileS, String referenceName) {
        this.readFile(binByTaxaFileS, referenceName);
    }

    /**
     * Reading binary file
     * @param binByTaxaFileS 
     */
    public BinByTaxa (String binByTaxaFileS) {
        this.readBinaryFile(binByTaxaFileS);
    }
    
    private void readBinaryFile (String binByTaxaFileS) {
        try {
            DataInputStream dis = IoUtils.getBinaryReader(binByTaxaFileS);
            int taxaNum = dis.readInt();
            taxa = new String[taxaNum];
            this.refTaxonIndex = dis.readInt();
            int chromNum = dis.readInt();
            chrID = new int[chromNum];
            binStart = new int[chromNum][];
            depth = new int[chromNum][][];
            ratio = new float[chromNum][][];
            pValue = new float[chromNum][][];
            for (int i = 0; i < chrID.length; i++) chrID[i] = dis.readInt();
            for (int i = 0; i < binStart.length; i++) binStart[i] = new int[dis.readInt()];
            for (int i = 0; i < taxa.length; i++) taxa[i] = dis.readUTF();
            for (int i = 0; i < chromNum; i++) {
                for (int j = 0; j < binStart[i].length; j++) {
                    binStart[i][j] = dis.readInt();
                    for (int k = 0; k < taxaNum; k++) {
                        depth[i][j][k] = dis.readInt();
                        ratio[i][j][k] = dis.readFloat();
                        pValue[i][j][k] = dis.readFloat();
                    }
                }
            }
            dis.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
   
    private void readFile (String binByTaxaFileS, String referenceName) {
        int numTaxa = 0;
        int[] numBins = null;
        try {
            BufferedReader br = IoUtils.getTextReader(binByTaxaFileS);
            TIntArrayList chrList = new TIntArrayList();
            TIntArrayList binList = new TIntArrayList();
            String[] tem = br.readLine().split("\t");
            numTaxa = tem.length-2;
            taxa = new String[numTaxa];
            for (int i = 0; i < numTaxa; i++) {
                taxa[i] = tem[i+2];
            }
            this.refTaxonIndex = Arrays.binarySearch(taxa, referenceName);
            if (refTaxonIndex < 0) {
                System.out.println("Please set the reference correctly. Program aborts");
                System.exit(1);
            }
            String temp;
            while ((temp = br.readLine()) != null) {
                tem = temp.substring(0, 50).split("\t");
                chrList.add(Integer.valueOf(tem[0]));
                binList.add(Integer.valueOf(tem[1]));
            }
            int[] allChr = chrList.toArray(new int[chrList.size()]);
            int[] allBin = binList.toArray(new int[binList.size()]);
            TIntHashSet chrSet = new TIntHashSet();
            for (int i = 0; i < allChr.length; i++) {
                chrSet.add(allChr[i]);
            }
            chrID = chrSet.toArray(new int[chrSet.size()]);
            Arrays.sort(chrID);
            numBins = new int[chrID.length];
            for (int i = 0; i < allChr.length; i++) {
                numBins[allChr[i]-1]++;
            }
            binStart = new int[chrID.length][];
            depth = new int[chrID.length][][];
            for (int i = 0; i < chrID.length; i++) {
                binStart[i] = new int[numBins[i]];
                depth[i] = new int[numBins[i]][numTaxa];
            }
            br = IoUtils.getTextReader(binByTaxaFileS);
            br.readLine();
            for (int i = 0; i < this.getNumOfChromsome(); i++) {
                for (int j = 0; j < this.getNumOfBinOnChromosome(i); j++) {
                    tem = br.readLine().split("\t");
                    binStart[i][j] = Integer.valueOf(tem[1]);
                    for (int k = 0; k < this.getNumOfTaxa(); k++) {
                        depth[i][j][k] = Integer.valueOf(tem[k+2]);
                    }
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("BBT read from " + binByTaxaFileS);
        this.calRatioPValue();
    }
    
    private void calRatioPValue () {
        ratio = new float[this.getNumOfChromsome()][][];
        pValue = new float[this.getNumOfChromsome()][][];
        for (int i = 0; i < ratio.length; i++) {
            ratio[i] = new float[this.getNumOfBinOnChromosome(i)][this.getNumOfTaxa()];
            pValue[i] = new float[this.getNumOfBinOnChromosome(i)][this.getNumOfTaxa()];
        }
        long[] taxaDepth = this.getTaxaDepth();
        for (int i = 0; i < this.getNumOfChromsome(); i++) {
            for (int j = 0; j < this.getNumOfBinOnChromosome(i); j++) {
                double refDep = (double)depth[i][j][refTaxonIndex]/taxaDepth[refTaxonIndex];
                for (int k = 0; k < this.getNumOfTaxa(); k++) {
                    double comDep = (double)depth[i][j][k]/taxaDepth[k];
                    double r;
                    double p = 1;
                    if (refDep == 0) {
                        r = Double.NaN;
                    }
                    else {
                        r = comDep/refDep;
                        long[][] count = new long[2][2];
                        count[0][0] = depth[i][j][k];
                        count[0][1] = taxaDepth[k]-count[0][0];
                        count[1][0] = depth[i][j][refTaxonIndex];
                        count[1][1] = taxaDepth[refTaxonIndex]-count[1][0];
                        try {
                            p = TestUtils.chiSquareTest(count);
                        }
                        catch (Exception e) {
                            e.printStackTrace();
                        }
                    }
                    if (Double.isNaN(r)) {
                        ratio[i][j][k] = 0;
                    }
                    else if (r == 0) {
                        ratio[i][j][k] = -15;
                    }
                    else {
                        ratio[i][j][k] = (float)(Math.log(r)/Math.log(2));
                    }
                    pValue[i][j][k] = -(float)Math.log10(p);
                }
            }
            System.out.println("Finished calculation on chromosome " + String.valueOf(i+1));
        }
    }    
    
    public ArrayList<CNVInfo> scanCNV (String taxon, float positiveRatio, float negativeRatio, float pValue) {
        int index = Arrays.binarySearch(taxa, taxon);
        if (index < 0) {
            System.out.println("Please set the reference correctly. Program aborts");
            System.exit(1);
        }
        return this.scanCNV(index, positiveRatio, negativeRatio, pValue);
    }
    
    public ArrayList<CNVInfo> scanCNV (int taxonIndex, float positiveRatio, float negativeRatio, float pValue) {
        ArrayList<CNVInfo> cnvList = new ArrayList();
        int binSize = this.binStart[0][1]-this.binStart[0][0];
        for (int i = 0; i < this.getNumOfChromsome(); i++) {
            boolean ifInitial = true;
            int size = 0;
            int start = Integer.MIN_VALUE;
            for (int j = 0; j < this.getNumOfBinOnChromosome(i); j++) {
                if (this.ratio[i][j][taxonIndex] < negativeRatio && this.pValue[i][j][taxonIndex] > pValue) {
                    if (ifInitial) {
                      start = this.binStart[i][j];
                      size+=binSize;
                      ifInitial = false;
                    }
                    else {
                      size+=binSize;
                    }
                }
                else {
                    CNVInfo inf = new CNVInfo (this.taxa[taxonIndex], this.chrID[i], start, size, true);
                    cnvList.add(inf);
                    size = 0;
                    ifInitial = true;
                }
            }
        }
        for (int i = 0; i < this.getNumOfChromsome(); i++) {
            boolean ifInitial = true;
            int size = 0;
            int start = Integer.MIN_VALUE;
            for (int j = 0; j < this.getNumOfBinOnChromosome(i); j++) {
                if (this.ratio[i][j][taxonIndex] > positiveRatio && this.pValue[i][j][taxonIndex] > pValue) {
                    if (ifInitial) {
                      start = this.binStart[i][j];
                      size+=binSize;
                      ifInitial = false;
                    }
                    else {
                      size+=binSize;
                    }
                }
                else {
                    CNVInfo inf = new CNVInfo (this.taxa[taxonIndex], this.chrID[i], start, size, true);
                    cnvList.add(inf);
                    size = 0;
                    ifInitial = true;
                }
            }
        }
        return cnvList;
    }
    
    public long[] getTaxaDepth () {
        long[] taxaDepth = new long[this.getNumOfTaxa()];
        for (int i = 0; i < this.getNumOfChromsome(); i++) {
            for (int j = 0; j < this.getNumOfBinOnChromosome(i); j++) {
                for (int k = 0; k < this.getNumOfTaxa(); k++) {
                    taxaDepth[k]+= depth[i][j][k];
                }
            }
        }
        return taxaDepth;
    }
    
    
    public void writeRatioFile (String ratioFileS) {
        try {
            BufferedWriter bw = IoUtils.getTextWriter(ratioFileS);
            bw.write("Chr\tBinStart");
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < this.getNumOfTaxa(); i++) {
                sb.append("\t").append(taxa[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < this.getNumOfChromsome(); i++) {
                for (int j = 0; j < this.getNumOfBinOnChromosome(i); j++) {
                    sb = new StringBuilder();
                    sb.append(chrID[i]).append("\t").append(binStart[i][j]);
                    for (int k = 0; k < this.getNumOfTaxa(); k++) {
                        sb.append("\t").append(ratio[i][j][k]);
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void writePValueFile (String pValueFileS) {
        try {
            BufferedWriter bw = IoUtils.getTextWriter(pValueFileS);
            bw.write("Chr\tBinStart");
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < this.getNumOfTaxa(); i++) {
                sb.append("\t").append(taxa[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < this.getNumOfChromsome(); i++) {
                for (int j = 0; j < this.getNumOfBinOnChromosome(i); j++) {
                    sb = new StringBuilder();
                    sb.append(chrID[i]).append("\t").append(binStart[i][j]);
                    for (int k = 0; k < this.getNumOfTaxa(); k++) {
                        sb.append("\t").append(pValue[i][j][k]);
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private String get3DigitDouble (double value) {
        DecimalFormat df = new DecimalFormat("#.000");
        return df.format(value);
    }
    
    public void writeTaxaFile (String taxaFileS) {
        try {
            BufferedWriter bw = IoUtils.getTextWriter(taxaFileS);
            bw.write("Taxa");
            bw.newLine();
            for (int i = 0; i < taxa.length; i++) {
                bw.write(taxa[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void writeTextBinByTaxaFile (String binByTaxaFileS) {
        try {
            BufferedWriter bw = IoUtils.getTextWriter(binByTaxaFileS);
            bw.write("Chr\tBinStart");
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < this.getNumOfTaxa(); i++) {
                sb.append("\t").append(taxa[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < this.getNumOfChromsome(); i++) {
                for (int j = 0; j < this.getNumOfBinOnChromosome(i); j++) {
                    sb = new StringBuilder();
                    sb.append(chrID[i]).append("\t").append(binStart[i][j]);
                    for (int k = 0; k < this.getNumOfTaxa(); k++) {
                        sb.append("\t").append(depth[i][j][k]);
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void writeBinaryBinByTaxaFile (String binByTaxaFileS) {
        try {
            DataOutputStream dos = IoUtils.getBinaryWriter(binByTaxaFileS);
            dos.writeInt(taxa.length);
            dos.writeInt(this.refTaxonIndex);
            dos.writeInt(chrID.length);
            for (int i = 0; i < this.getNumOfChromsome(); i++) dos.writeInt(chrID[i]);
            for (int i = 0; i < this.getNumOfChromsome(); i++) dos.writeInt(binStart[i].length);
            for (int i = 0; i < taxa.length; i++) dos.writeUTF(taxa[i]);
            for (int i = 0; i < this.getNumOfChromsome(); i++) {
                for (int j = 0; j < this.getNumOfBinOnChromosome(i); j++) {
                    dos.writeInt(binStart[i][j]);
                    for (int k = 0; k < this.getNumOfTaxa(); k++) {
                        dos.writeInt(depth[i][j][k]);
                        dos.writeFloat(ratio[i][j][k]);
                        dos.writeFloat(pValue[i][j][k]);
                    }
                }
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public String getTaxon (int taxaIndex) {
        return this.taxa[taxaIndex];
    }
    
    public int getTaxonIndex (String taxonName) {
        return Arrays.binarySearch(taxa, taxonName);
    }
    
    public int getNumOfChromsome () {
        return this.chrID.length;
    }
    
    public int getNumOfBinOnChromosome (int chromIndex) {
        return binStart[chromIndex].length;
    }
    
    public int getNumOfTaxa () {
        return this.taxa.length;
    }
}
