/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.cassava.wgs;


import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class Reformat {
    
    public Reformat () {
        this.Reformat1();
    }
    
    public void Reformat1 () {
        String infileS = "E:\\Research\\cassava\\revision\\format\\cultiProCombinedSIFT.txt";
        String outfileS = "E:\\Research\\cassava\\revision\\format\\cultiProCombinedSIFT_order.txt";
        String outfileS2 = "E:\\Research\\cassava\\revision\\format\\cultiProCombinedSIFT_pos.txt";
        HashMap<Integer, String> posRecordMap = new HashMap ();
        String header = null;
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            header = br.readLine();
            String temp = null;
            TIntArrayList posList = new TIntArrayList();
            while ((temp = br.readLine()) != null) {
                List<String> l = FStringUtils.fastSplit(temp);
                int posValue = Integer.valueOf(l.get(0))*100_000_000+Integer.valueOf(l.get(1));
                posRecordMap.put(posValue, temp);
                posList.add(posValue);
            }
            br.close();
            int[] pos = posList.toArray();
            Arrays.sort(pos);
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write(header);
            bw.newLine();
            for (int i = 0; i < pos.length; i++) {
                bw.write(posRecordMap.get(pos[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            bw = IoUtils.getTextWriter(outfileS2);
            bw.write("Chr\tPos");
            bw.newLine();
            for (int i = 0; i < pos.length; i++) {
                List<String> l = FStringUtils.fastSplit(posRecordMap.get(pos[i]));
                bw.write(l.get(0)+"\t"+l.get(1));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
}
