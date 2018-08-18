/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.cassava.wgs;

import format.Table;
import java.io.BufferedWriter;
import java.util.Arrays;
import java.util.HashSet;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class Location {
    
    public Location () {
        this.summarizeLocation();
    }
    
    public void summarizeLocation () {
        String infileS = "E:\\Research\\cassava\\revision\\geographic\\location.txt";
        String outfileS = "E:\\Research\\cassava\\revision\\geographic\\location_summary.txt";
        Table t = new Table (infileS);
        HashSet<String> loSet = new HashSet();
        for (int i = 0; i < t.getRowNumber(); i++) {
            loSet.add(t.content[i][0]);
        }
        String[] locations = loSet.toArray(new String[loSet.size()]);
        int[] count = new int[locations.length];
        Arrays.sort(locations);
        for (int i = 0; i < t.getRowNumber(); i++ ) {
            int index = Arrays.binarySearch(locations, t.content[i][0]);
            count[index]++;
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Location\tSampleSize"); 
            bw.newLine();
            for (int i = 0 ; i < locations.length; i++) {
                bw.write(locations[i]+"\t"+String.valueOf(count[i]));
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
