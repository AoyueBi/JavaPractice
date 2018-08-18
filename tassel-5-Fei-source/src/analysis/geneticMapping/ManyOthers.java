/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.geneticMapping;

import format.Fasta;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import net.maizegenetics.dna.map.TagsOnGeneticMap;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.tag.TagsByTaxa;

/**
 *
 * @author fl262
 */
public class ManyOthers {
    
    public ManyOthers () {
        //this.mkGWASAccuracyTable();
        //this.skimApproachCoverage();
        this.countSharedAndMissing();
    }
    
    public void countSharedAndMissing () {
        String hapMapH5 = "E:\\Research\\geneticMapping\\accuracyGeneticDistance\\combine.hmp.h5";
        GenotypeTable gt = ImportUtils.readGuessFormat(hapMapH5);
        int taxaNum = gt.taxa().numberOfTaxa();
        int size = 30;
        int[] taxaIndex = new int[size];
        for (int i = 0; i < taxaIndex.length; i++) {
            taxaIndex[i] = (int)(Math.random()*taxaNum);
        }
        ArrayList<Double> shareList = new ArrayList();
        ArrayList<Double> missingList = new ArrayList();
        for (int i = 0; i < taxaIndex.length-1; i++) {
            for (int j = i+1; j < taxaIndex.length; j++) {
                int shareCount = 0;
                int missingCount = 0;
                for (int k = 0; k < gt.numberOfSites(); k++) {
                    if (gt.genotype(taxaIndex[i], k)>-1 && gt.genotype(taxaIndex[j], k)>-1) shareCount++;
                    else if (gt.genotype(taxaIndex[i], k)<0 && gt.genotype(taxaIndex[j], k)<0) missingCount++;
                }
                shareList.add((double)shareCount/gt.numberOfSites());
                missingList.add((double)missingCount/gt.numberOfSites());
                
            }
        }
        double averShare = 0;
        for (int i = 0; i < shareList.size(); i++) {
            averShare += shareList.get(i)/shareList.size();
        }
        double averMissing = 0;
        for (int i = 0; i < missingList.size(); i++) {
            averMissing += missingList.get(i)/missingList.size();
        }
        System.out.println(averShare);
        System.out.println(averMissing);
    }
    
    public void skimApproachCoverage () {
        String inputfileS = "M:\\CML247_soap_5x.scafSeq_2.txt";
        Fasta f = new Fasta(inputfileS);
        f.sortRecordByLengthDescending();
        System.out.println(f.getTotalSeqLength());
        System.out.println(f.getL50());
        System.out.println(f.getSeqLength(f.getSeqNumber()-1));
    }
    
    public void mkGWASAccuracyTable () {
        //String predictionFileS = "E:\\Research\\geneticMapping\\accuracyTableGWAS\\box_M5Rules_Min50_Num_G.txt";
        //String tableFileS = "E:\\Research\\geneticMapping\\accuracyTableGWAS\\accuracyTableGWAS.txt";
        //String predictionFileS = "E:\\Research\\geneticMapping\\accuracyTableGWAS\\box_M5Rules_Min50_Num_GJ.txt";
        //String tableFileS = "E:\\Research\\geneticMapping\\accuracyTableGWAS\\accuracyTableGWASJoint.txt";
        String predictionFileS = "E:\\Research\\geneticMapping\\accuracyTableGWAS\\box_M5Rules_Min50_Num_J.txt";
        String tableFileS = "E:\\Research\\geneticMapping\\accuracyTableGWAS\\accuracyTableJoint.txt";
        int[] predictionCut = {10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000};
        int[] accuracyCut = {10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000, 20000000, 50000000};
        double[] logPreCut = new double[predictionCut.length];
        for (int i = 0; i < logPreCut.length; i++) logPreCut[i] = Math.log10(predictionCut[i]);
        double[] logAccCut = new double[accuracyCut.length];
        for (int i = 0; i < logAccCut.length; i++) logAccCut[i] = Math.log10(accuracyCut[i]);
        double[] actualValue = null;
        double[] predictionValue = null;
        try {
            BufferedReader br = new BufferedReader(new FileReader(predictionFileS), 65536);
            br.readLine();
            String temp;
            String[] tem;
            int cnt = 0;
            while ((temp = br.readLine()) != null) cnt++;
            actualValue = new double[cnt];
            predictionValue = new double[cnt];
            br = new BufferedReader(new FileReader(predictionFileS), 65536);
            br.readLine();
            for (int i = 0; i < cnt; i++) {
                tem = br.readLine().split("\\s+");
                actualValue[i] = Double.valueOf(tem[2]);
                predictionValue[i] = Double.valueOf(tem[3]);
            }
            br.close();
            BufferedWriter bw = new BufferedWriter (new FileWriter(tableFileS), 65536);
            bw.write("PredictionCutoff\tProportionRemain");
            for (int i = 0; i < accuracyCut.length; i++) {
                bw.write("\t<"+String.valueOf(accuracyCut[i]/1000)+ "Kb");
            }
            bw.newLine();
            for (int i = 0; i < logPreCut.length; i++) {
                double[] ratio = new double[accuracyCut.length];
                cnt = 0;
                for (int j = 0; j < actualValue.length; j++) {
                    if (predictionValue[j] > logPreCut[i]) continue;
                    for (int k = 0; k < logAccCut.length; k++) {
                        if (actualValue[j] < logAccCut[k]) ratio[k]++;
                    }
                    cnt++;
                }
                bw.write(String.valueOf(predictionCut[i]/1000)+"Kb\t"+String.valueOf((double)cnt/actualValue.length));
                for (int j = 0; j < logAccCut.length; j++) {
                    bw.write("\t"+String.valueOf((double)ratio[j]/cnt));
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("Accuracy table is at " + tableFileS);
    }
}
