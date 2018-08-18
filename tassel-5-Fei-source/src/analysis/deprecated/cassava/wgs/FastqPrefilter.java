/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.cassava.wgs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class FastqPrefilter {
    
    public FastqPrefilter() {
//        String r1Source = "O:\\Cassava\\QC\\samExtraction\\3422_7332_10687_HA63NADXX_I000070_ACAGTG_R1.fastq";
//        String r2Source = "O:\\Cassava\\QC\\samExtraction\\3422_7332_10687_HA63NADXX_I000070_ACAGTG_R2.fastq";
//        String samSource = "O:\\Cassava\\QC\\samExtraction\\I000070_HA63NADXX_unsortedt.sam";
//        String outputFileS1 = "O:\\Cassava\\QC\\samExtraction\\output\\temp_3422_7332_10687_HA63NADXX_I000070_ACAGTG_R1.fastq";
//        String outputFileS2 = "O:\\Cassava\\QC\\samExtraction\\output\\temp_3422_7332_10687_HA63NADXX_I000070_ACAGTG_R2.fastq";
//        String[] args = {r1Source, r2Source, samSource, outputFileS1, outputFileS2};
//        this.process(args);
        
        String r1Source = "O:\\Cassava\\QC\\samExtraction\\3422_7332_10687_HA63NADXX_I000070_ACAGTG_R1_special.fastq";
        String r2Source = "O:\\Cassava\\QC\\samExtraction\\3422_7332_10687_HA63NADXX_I000070_ACAGTG_R2_special.fastq";
        String samSource = "O:\\Cassava\\QC\\samExtraction\\I000070_HA63NADXX_unsortedt_special.sam";
        String outputFileS1 = "O:\\Cassava\\QC\\samExtraction\\output\\temp_3422_7332_10687_HA63NADXX_I000070_ACAGTG_R1.fastq";
        String outputFileS2 = "O:\\Cassava\\QC\\samExtraction\\output\\temp_3422_7332_10687_HA63NADXX_I000070_ACAGTG_R2.fastq";
        String[] args = {r1Source, r2Source, samSource, outputFileS1, outputFileS2};
        this.processSpecial(args);
    }
    
    public FastqPrefilter (String[] args) {
        //this.process(args);
        this.processSpecial(args);
    }
    
    //Rare format in Fastq, where R1/R2 reads are represented by /1 or /2
    public void processSpecial (String[] args) {
        ArrayList<String> nameList = new ArrayList();
        try {
            BufferedReader br = IoUtils.getTextReader(args[2]);
            String temp;
            ArrayList<String[]> sameList = new ArrayList();
            while (!(temp = br.readLine()).startsWith("@PG")) {}
            temp = br.readLine();
            String[] tem = temp.substring(0, 100).split("\t");
            String currentName = tem[0];
            if (Integer.valueOf(tem[1]) < 2000) sameList.add(tem);
            while ((temp = br.readLine()) != null) {
                tem = temp.substring(0, 100).split("\t");
                if (tem[0].equals(currentName)) {
                    if (Integer.valueOf(tem[1]) < 2000) sameList.add(tem);
                    //sameList.add(tem);
                }
                else {
                    boolean flag = true;
                    for (int i = 0; i < sameList.size(); i++) {
                        if (sameList.get(i)[2].equals("*")) {
                            nameList.add(sameList.get(i)[0]);
                            flag = false;
                            break;
                        }
                    }
                    if (flag) {
                        for (int i = 0; i < sameList.size(); i++) {
                            int chr = Integer.valueOf(sameList.get(i)[2]);
                            if (chr < 23) continue;
                            if (Integer.valueOf(sameList.get(i)[4]) < 20) {
                                nameList.add(sameList.get(i)[0]);
                                break;
                            }
                        }
                    }
                    currentName = tem[0];
                    sameList.clear();
                    sameList.add(tem);
                }
            }
            boolean flag = true;
            for (int i = 0; i < sameList.size(); i++) {
                if (sameList.get(i)[2].equals("*")) {
                    nameList.add(sameList.get(i)[0]);
                    flag = false;
                    break;
                }
            }
            if (flag) {
                for (int i = 0; i < sameList.size(); i++) {
                    int chr = Integer.valueOf(sameList.get(i)[1]);
                    if (chr < 23) continue;
                    if (Integer.valueOf(sameList.get(i)[1]) < 20) {
                        nameList.add(sameList.get(i)[0]);
                        break;
                    }
                }
            }
            br.close();
            BufferedReader br1 = IoUtils.getTextReader(args[0]);
            BufferedReader br2 = IoUtils.getTextReader(args[1]);
            BufferedWriter bw1 = IoUtils.getTextWriter(args[3]);
            BufferedWriter bw2 = IoUtils.getTextWriter(args[4]);
            String[] name = nameList.toArray(new String[nameList.size()]);
            Arrays.sort(name);
            String temp1 = null;
            String temp2 = null;
            while ((temp1 = br1.readLine()) != null) {
                temp2 = br2.readLine();
                String query = temp1.split("/")[0].replaceFirst("@", "");
                int index = Arrays.binarySearch(name, query);
                if (index < 0) {
                    for (int i = 0; i < 3; i++) {
                        br1.readLine();
                        br2.readLine();
                    }
                }
                else {
                    bw1.write(temp1);bw1.newLine();
                    bw2.write(temp2);bw2.newLine();
                    for (int i = 0; i < 3; i++) {
                        bw1.write(br1.readLine());bw1.newLine();
                        bw2.write(br2.readLine());bw2.newLine();
                    }
                }
            }
            bw1.flush();
            bw2.flush();
            bw1.close();
            bw2.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        } 
    }
    
    public void process (String[] args) {
        ArrayList<String> nameList = new ArrayList();
        try {
            BufferedReader br = IoUtils.getTextReader(args[2]);
            String temp;
            ArrayList<String[]> sameList = new ArrayList();
            while (!(temp = br.readLine()).startsWith("@PG")) {}
            temp = br.readLine();
            String[] tem = temp.substring(0, 100).split("\t");
            String currentName = tem[0];
            if (Integer.valueOf(tem[1]) < 2000) sameList.add(tem);
            while ((temp = br.readLine()) != null) {
                tem = temp.substring(0, 100).split("\t");
                if (tem[0].equals(currentName)) {
                    if (Integer.valueOf(tem[1]) < 2000) sameList.add(tem);
                    //sameList.add(tem);
                }
                else {
                    boolean flag = true;
                    for (int i = 0; i < sameList.size(); i++) {
                        if (sameList.get(i)[2].equals("*")) {
                            nameList.add(sameList.get(i)[0]);
                            flag = false;
                            break;
                        }
                    }
                    if (flag) {
                        for (int i = 0; i < sameList.size(); i++) {
                            int chr = Integer.valueOf(sameList.get(i)[2]);
                            if (chr < 23) continue;
                            if (Integer.valueOf(sameList.get(i)[4]) < 20) {
                                nameList.add(sameList.get(i)[0]);
                                break;
                            }
                        }
                    }
                    currentName = tem[0];
                    sameList.clear();
                    sameList.add(tem);
                }
            }
            boolean flag = true;
            for (int i = 0; i < sameList.size(); i++) {
                if (sameList.get(i)[2].equals("*")) {
                    nameList.add(sameList.get(i)[0]);
                    flag = false;
                    break;
                }
            }
            if (flag) {
                for (int i = 0; i < sameList.size(); i++) {
                    int chr = Integer.valueOf(sameList.get(i)[1]);
                    if (chr < 23) continue;
                    if (Integer.valueOf(sameList.get(i)[1]) < 20) {
                        nameList.add(sameList.get(i)[0]);
                        break;
                    }
                }
            }
            br.close();
            BufferedReader br1 = IoUtils.getTextReader(args[0]);
            BufferedReader br2 = IoUtils.getTextReader(args[1]);
            BufferedWriter bw1 = IoUtils.getTextWriter(args[3]);
            BufferedWriter bw2 = IoUtils.getTextWriter(args[4]);
            String[] name = nameList.toArray(new String[nameList.size()]);
            Arrays.sort(name);
            String temp1 = null;
            String temp2 = null;
            while ((temp1 = br1.readLine()) != null) {
                temp2 = br2.readLine();
                String query = temp1.split(" ")[0].replaceFirst("@", "");
                int index = Arrays.binarySearch(name, query);
                if (index < 0) {
                    for (int i = 0; i < 3; i++) {
                        br1.readLine();
                        br2.readLine();
                    }
                }
                else {
                    bw1.write(temp1);bw1.newLine();
                    bw2.write(temp2);bw2.newLine();
                    for (int i = 0; i < 3; i++) {
                        bw1.write(br1.readLine());bw1.newLine();
                        bw2.write(br2.readLine());bw2.newLine();
                    }
                }
            }
            bw1.flush();
            bw2.flush();
            bw1.close();
            bw2.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        } 
    }
    
    public static void main (String[] args) {
        new FastqPrefilter(args);
    }
}
