/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.collaboration;

import format.Range;
import format.Table;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class Randy {
    
    public Randy () {
        this.selectRegion();
    }
    
    public void selectRegion () {
        String positionFileS = "M:\\collaboration\\randy\\position.txt";
        String outputDirS = "M:\\collaboration\\randy\\region\\";
        String scoreDirS = "M:\\production\\maf\\annotations\\siteUniqueness\\data\\combined_6Genomes\\";
        Table t = new Table (positionFileS);
        HashMap<String, Range> regionRangeMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            Range r = new Range(t.getIntValue(i, 1), Integer.valueOf(t.content[i][2].replaceAll(",", "")), Integer.valueOf(t.content[i][3].replaceAll(",", ""))+1);
            regionRangeMap.put(t.content[i][0], r);
        }
        Set<Map.Entry<String, Range>> set = regionRangeMap.entrySet();
        set.parallelStream().forEach(s -> {
            Range r = s.getValue();
            String regionName = s.getKey();
            int chr = r.getRangeChromosome();
            int start = r.getRangeStart();
            int end = r.getRangeEnd();
            String infileS = "chr"+FStringUtils.getNDigitNumber(3, chr)+"_uniqueness.txt.gz";
            infileS = new File(scoreDirS, infileS).getAbsolutePath();
            String outfileS = new File(outputDirS, regionName).getAbsolutePath();
            try {
                BufferedReader br = IoUtils.getTextGzipReader(infileS);
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                bw.write(br.readLine());
                bw.newLine();
                String temp = null;
                int cnt=0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%10000000 == 0) System.out.println(String.valueOf(cnt)+"\tchr\t"+ String.valueOf(chr));
                    List<String> l = FStringUtils.fastSplit(temp, "\t");
                    int pos = Integer.valueOf(l.get(1));
                    if (pos < start) continue;
                    else if (pos < end) {
                        bw.write(temp);
                        bw.newLine();
                    }
                    else break;
                }
                br.close();
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
    
}
