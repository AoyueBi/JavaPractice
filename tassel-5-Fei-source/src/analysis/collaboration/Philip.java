/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.collaboration;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class Philip {
    
    public Philip () {
        this.convertToBed();
    }
    
    public void convertToBed () {
        String infileS = "M:\\collaboration\\philip\\debug\\raw\\pos.txt";
        String outfileS = "M:\\collaboration\\philip\\debug\\raw\\pos_bed.txt";
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                tem = tem[0].split(":");
                StringBuilder sb = new StringBuilder(tem[0]);
                sb.append("\t").append(tem[1]).append("\t").append(Integer.valueOf(tem[1])+1);
                bw.write(sb.toString());
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
        
}
