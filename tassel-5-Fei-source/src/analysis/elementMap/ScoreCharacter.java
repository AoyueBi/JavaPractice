/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.elementMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.Arrays;
import java.util.List;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class ScoreCharacter {
    
    public ScoreCharacter() {
        this.getSubset();
        this.getScoreDistribution();
    }
    
    
    private void getScoreDistribution () {
        int size = 20000;
        String infileS = "M:\\production\\elementMap\\combinedScore\\chr010_combinedScore_16.txt";
        String outfileS = "M:\\production\\elementMap\\scoreCharacter\\chr10_scoreDistribution_16.txt";
        int binNumber = 1000;
        double interval = 0.05;
        double[] binStart = new double[binNumber];
        int[] binCounts = new int[binNumber];
        for (int i = 0; i < binNumber; i++) binStart[i] = i*interval;
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            String temp = br.readLine();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                List<String> l = FStringUtils.fastSplit(temp, "\t");
                if (l.get(2).startsWith("N")) continue;
                double value = Double.valueOf(l.get(2));
                int index = Arrays.binarySearch(binStart, value);
                if (index < 0) index = -index -2;
                if (index == binNumber) index = binNumber-1;
                binCounts[index]++;
                cnt++;
                if (cnt%10000000 == 0) System.out.println(cnt);
            }
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("CpScore\tCount\tProportion");
            bw.newLine();
            for (int i = 0; i < binNumber; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(binStart[i]).append("\t").append(binCounts[i]).append("\t").append((double)binCounts[i]/cnt);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void getSubset () {
        String inputFileS = "M:\\production\\elementMap\\combinedScore\\chr010_combinedScore_16.txt";
        String outfileS = "M:\\production\\elementMap\\scoreCharacter\\subset\\chr10_subset_16";
        try {
            
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
