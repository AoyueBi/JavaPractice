/*
* To change this license header, choose License Headers in Project Properties.
* To change this template file, choose Tools | Templates
* and open the template in the editor.
*/

package analysis.geneticMapping;

import format.Fasta;
import format.Table;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeSet;
import net.maizegenetics.analysis.data.MergeGenotypeTablesPlugin;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import net.maizegenetics.dna.tag.TagsByTaxaByte;
import net.maizegenetics.dna.tag.TagsByTaxaByteHDF5TagGroups;
import utils.FArrayUtils;
import utils.IoUtils;

/**
 *
 * @author Fei Lu
 */
public class GMUtils {
    
    public GMUtils () {}
    
    public void selectBestCML247 (String maizeSorghumRatioFileS, String qualityFileS, String cml247FastaS, String cml247BestFastaS, int maxNum, double qualityCut, double anchorNumCut, double ratioCut) {
        HashMap<String, Double> ratioMap = new HashMap();
        Table t = new Table (maizeSorghumRatioFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            ratioMap.put(t.content[i][0], Double.valueOf(t.content[i][4]));
        }
        HashMap<String, Double> qualityMap = new HashMap();
        HashMap<String, Integer> anchorNumMap = new HashMap();
        t = new Table (qualityFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            qualityMap.put(t.content[i][0], Double.valueOf(t.content[i][4]));
            anchorNumMap.put(t.content[i][0], Integer.valueOf(t.content[i][3]));
        }
        Fasta f = new Fasta (cml247FastaS);
        boolean[] ifOut = new boolean[f.getSeqNumber()];
        int cnt = 0;
        ArrayList<Integer> lengthList = new ArrayList();
        for (int i = 0; i < f.getSeqNumber(); i++) {
            String name = f.getName(i);
            if (name.equals("411")) continue;
            if (ratioMap.get(name) == null) continue;
            if (ratioMap.get(name) < ratioCut) continue;
            if (qualityMap.get(name) < qualityCut) continue;
            if (anchorNumMap.get(name) < anchorNumCut) continue;
            ifOut[i] = true;
            lengthList.add(f.getSeqLength(i));
            System.out.println(name+"\t"+f.getSeqLength(i));
            cnt++;
            if (cnt >= maxNum) break;
        }
        f.writeFasta(cml247BestFastaS, ifOut);
        Integer[] lengths = lengthList.toArray(new Integer[lengthList.size()]);
        int sum = 0;
        for (int i = 0; i < lengths.length; i++) sum+=lengths[i];
        System.out.println("Total length:\t" + sum);
        System.out.println("Median length:\t" + lengths[lengths.length/2]);
    }
    
    public void convertHapMapToH5 (String hapMapFileS, String h5FileS) {
        System.out.println("Start converting hapmap to h5");
        GenotypeTable gt = ImportUtils.readGuessFormat(hapMapFileS);
        ExportUtils.writeGenotypeHDF5(gt, h5FileS);
    }
    
    public void mergeHapMap (String inputDirS, String outfileS) {
        File[] fs = new File(inputDirS).listFiles();
        Arrays.sort(fs);
        TreeSet<String> allNameSet = new TreeSet ();
        String[][] chrTaxa = new String[fs.length][];
        String temp = null;
        String[] tem = null;
        for (int i = 0; i < fs.length; i++) {
            try {
                BufferedReader br = IoUtils.getTextReader(fs[i].getAbsolutePath());
                temp = br.readLine();
                tem = temp.split("\t");
                chrTaxa[i] = new String[tem.length-11];
                for (int j = 11; j < tem.length; j++) {
                    chrTaxa[i][j-11] = tem[j];
                    allNameSet.add(tem[j]);
                }
                Arrays.sort(chrTaxa[i]);
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        String[] allName = allNameSet.toArray(new String[allNameSet.size()]);
        ArrayList<String> sharedNameList = new ArrayList();
        for (int i = 0; i < allName.length; i++) {
            int flag = 0;
            for (int j = 0; j < chrTaxa.length; j++) {
                int hit = Arrays.binarySearch(chrTaxa[j], allName[i]);
                if (hit < 1) {
                    flag = 1;
                    break;
                }
            }
            if (flag == 1) continue;
            sharedNameList.add(allName[i]);
        }
        String[] sharedName= sharedNameList.toArray(new String[sharedNameList.size()]);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            for (int i = 0; i < fs.length; i++) { 
                System.out.println("Start processing " + fs[i].getAbsolutePath());
                BufferedReader br = IoUtils.getTextReader(fs[i].getAbsolutePath());
                temp = br.readLine();
                tem = temp.split("\t");
                int[] index = new int[sharedName.length];
                boolean[] ifMatch = new boolean[sharedName.length];
                for (int j = 11; j < tem.length; j++) {
                    int hit = Arrays.binarySearch(sharedName, tem[j]);
                    if (hit < 0) continue;
                    if (ifMatch[hit]) continue;
                    index[hit] = j;
                    ifMatch[hit] = true;
                }
                if (i == 0) {
                    for (int j = 0; j < 11; j++) {
                        bw.write(tem[j]+"\t");
                    }
                    for (int j = 0; j < sharedName.length-1; j++) {
                        bw.write(sharedName[j]+"\t");
                    }
                    bw.write(sharedName[sharedName.length-1]);
                    bw.newLine();
                }
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    tem = temp.split("\t");
                    for (int j = 0; j < 11; j++) {
                        bw.write(tem[j]+"\t");
                    }
                    for (int j = 0; j < index.length-1; j++) {
                        if (tem[index[j]].contains("_")) System.out.println (tem[index[j]]+"\t"+index[j]);
                        bw.write(tem[index[j]]+"\t");
                    }
                    bw.write(tem[index[index.length-1]]);
                    bw.newLine();
                    cnt++;
                    if (cnt%10000 == 0) System.out.println("Processed " + cnt + " SNPs");
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("File written to " + outfileS);
    }
    
    public void convertTBTByte (String sourceTBT, String targetTBT) {
        TagsByTaxaByte tbtb = new TagsByTaxaByte (sourceTBT, FilePacking.Byte);
        TagsByTaxaByteHDF5TagGroups tbth5 = new TagsByTaxaByteHDF5TagGroups (tbtb, targetTBT);
    }
    
}
