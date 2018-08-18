/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.cassava.wgs;

import format.Table;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import utils.FArrayUtils;

/**
 *
 * @author fl262
 */
class BurdenCount {
    
    public BurdenCount () {
        ArrayList<String> list = new ArrayList();
        for (int i = 0; i < 10; i++) {
            int[] value = this.test();
            String s = "Pro: "+String.valueOf(value[0])+"\tLAC: "+String.valueOf(value[1])+"\tAfrican: "+String.valueOf(value[2]);
            list.add(s);
        }
        for (int i = 0; i < list.size(); i++) {
            System.out.println(list.get(i));
        }
    }
    
    public int[] test () {
        String africanVCF = "O:\\Cassava\\fixedLoad\\africaData_afterFiltering.recode.vcf";
        String lacVCF = "O:\\Cassava\\fixedLoad\\latinAmericaData_afterFiltering.recode.vcf";
        String proVCF = "O:\\Cassava\\fixedLoad\\progenitorsData_afterFiltering.recode.vcf";
        String dafFileS = "O:\\Cassava\\fixedLoad\\info_delMut_DAFrubber_GERPandTree2.txt";
        int chrNum = 18;
        int sampleSize = 16;
        HashMap<Integer, String>[] posDeMap = new HashMap[chrNum];
        for (int i = 0; i < chrNum; i++) {
            posDeMap[i] = new HashMap();
        }
        Table t = new Table (dafFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            int chrIndex = t.getIntValue(i, 0)-1;
            int pos = t.getIntValue(i, 1);
            posDeMap[chrIndex].put(pos, t.content[i][5]);
        }
        int afrCount = this.getFixed(africanVCF, sampleSize, posDeMap);
        int lacCount = this.getFixed(lacVCF, sampleSize, posDeMap);
        int proCount = this.getFixed(proVCF, sampleSize, posDeMap);
        int[] value = {proCount, lacCount, afrCount};
        return value;
    }
    
    private int getFixed (String infileS, int sampleSize, HashMap<Integer, String>[] posDeMap) {
        Table t = new Table (infileS);
        int[] indices = new int[t.getColumnNumber()-9];
        for (int i = 0; i < indices.length; i++) {
            indices[i] = i+9;
        }
        FArrayUtils.shuffleArray(indices);
        int[] index = new int[sampleSize];
        for (int i = 0; i < index.length; i++) {
            index[i] = indices[i];
        }
        Arrays.sort(index);
        int cnt = 0;
        for (int i = 0; i < t.getRowNumber(); i++) {
            int chrIndex = t.getIntValue(i, 0)-1;
            int pos = t.getIntValue(i, 1);
            String deA = posDeMap[chrIndex].get(pos);
            int deV = -1;
            if (deA.equals(t.content[i][3])) {
                deV = 0;
            }
            else if (deA.equals(t.content[i][4])) {
                deV = 1;
            }
            int count = 0;
            for (int j = 0; j < index.length; j++) {
                String[] temp = t.content[i][index[j]].split(":");
                temp = temp[0].split("\\|");
                for (int k = 0; k < temp.length; k++) {
                    if (Integer.valueOf(temp[k]) == deV) {
                        count++;
                    }
                }
            }
            if (count == 2*sampleSize) cnt++;   
        }
        return cnt;
    }
}
