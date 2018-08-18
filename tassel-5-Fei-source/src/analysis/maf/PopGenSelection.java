/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maf;

import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.List;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class PopGenSelection {
    
    public PopGenSelection () {
        this.assignCMToHmp();
    }
    
    void assignCMToHmp () {
        String geneticMapFileS = "M:\\production\\maf\\popgen\\source\\NAM_phasedImputed_1cM_CONVERTED_TO_V3.txt";
        String inputDirS = "M:\\production\\maf\\annotations\\siftScore\\001_hmp321Info\\gerpAncestral\\";
        String outFileDirS = "M:\\production\\maf\\popgen\\selection\\Hmp321_cm";
        int chrNum = 10;
        TIntArrayList[] cmList = new TIntArrayList[chrNum];
        TIntArrayList[] posList = new TIntArrayList[chrNum];
        int[][] cm = new int[chrNum][];
        int[][] poss = new int[chrNum][];
        for (int i = 0; i < chrNum; i++) {
            cmList[i] = new TIntArrayList();
            posList[i] = new TIntArrayList();
        }
        try {
            BufferedReader br = IoUtils.getTextReader(geneticMapFileS);
            String temp = br.readLine();
            while ((temp = br.readLine()) != null) {
                List<String> l = FStringUtils.fastSplit(temp, " ");
                if (l.get(6).startsWith("N")) continue;
                int chrIndex = Integer.valueOf(l.get(3)) - 1;
                cmList[chrIndex].add(Integer.valueOf(l.get(5)));
                posList[chrIndex].add(Integer.valueOf(l.get(6)));
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        for (int i = 0; i < chrNum; i++) {
            cm[i] = cmList[i].toArray();
            poss[i] = posList[i].toArray();
        }
        File[] fs = new File(inputDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            int chrIndex = Integer.valueOf(f.getName().replaceFirst("hmp32Info_chr", "").replaceFirst(".txt.gz", ""))-1;
            String outfileS = new File (outFileDirS, "hmp321_cm_"+FStringUtils.getNDigitNumber(3, chrIndex+1)+".txt.gz").getAbsolutePath();
            try {
                BufferedReader br = IoUtils.getTextGzipReader(f.getAbsolutePath());
                String temp = br.readLine();
                BufferedWriter bw = IoUtils.getTextGzipWriter(outfileS);
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    List<String> l = FStringUtils.fastSplit(temp);
                    String name = "s"+l.get(0)+"_"+l.get(1);
                    int chr = Integer.valueOf(l.get(0));
                    int pos = Integer.valueOf(l.get(1));
                    int index = Arrays.binarySearch(poss[chr-1], pos);
                    if (index < 0) index = -index-1;
                    if (index == poss[chr-1].length) index = poss[chr-1].length-1;
                    float cmValue = (float)((double)cm[chr-1][index]/100);
                    StringBuilder sb = new StringBuilder(name);
                    sb.append(" ").append(chr).append(" ").append(cmValue).append(" ").append(pos).append(" ").append(l.get(2)).append(" ").append(l.get(3));
                    bw.write(sb.toString());
                    bw.newLine();
                    cnt++;
                    if (cnt%1000000 ==0) System.out.println(String.valueOf(cnt)+"\t"+f.getName());
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
}
