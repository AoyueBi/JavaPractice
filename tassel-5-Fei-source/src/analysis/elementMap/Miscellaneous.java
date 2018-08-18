/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.elementMap;

import format.Table;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import utils.FArrayUtils;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class Miscellaneous {
    
    Miscellaneous () {
        
    }
    
    void getHighCoverageTaxa () {
        String inputFileS = "M:\\production\\elementMap\\taxa\\taxaFqMap.txt";
        String depthFileS = "M:\\production\\elementMap\\taxaDepth\\taxaDepth.txt";
        String outputFileS = "M:\\production\\elementMap\\taxa\\taxaFqMap_5x.txt";
        Table t = new Table (depthFileS);
        ArrayList<String> highList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (Double.valueOf(t.content[i][1]) < 5) continue;
            highList.add(t.content[i][0]);
        }
        String[] highTaxa = highList.toArray(new String[highList.size()]);
        Arrays.sort(highTaxa);
        try {
            BufferedReader br = IoUtils.getTextReader(inputFileS);
            BufferedWriter bw = IoUtils.getTextWriter(outputFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                int index = Arrays.binarySearch(highTaxa, tem[0]);
                if (index < 0) continue;
                bw.write(temp);
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
    
    void getSingleCopySite () {
        String inputFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\chr010_uniqueness_6Genomes.txt";
        String outputFileS = "M:\\production\\elementMap\\taxaDepth\\singleCopySite.txt";
        TIntArrayList sList = new TIntArrayList();
        int siteNumber = 200000;
        try{
            BufferedReader br = IoUtils.getTextReader(inputFileS);
            String temp = br.readLine();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnt++;
                if (cnt%10000000==0) System.out.println(cnt);
                List<String> l = FStringUtils.fastSplit(temp);
                if (l.get(2).startsWith("N")) continue;
                float c = Float.valueOf(l.get(2));
                float v = Float.valueOf(l.get(3));
                if (c != 1) continue;
                if (v != 0) continue;
                sList.add(Integer.valueOf(l.get(1)));
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        int[] sites = sList.toArray();
        FArrayUtils.shuffleArray(sites);
        if (sites.length < siteNumber) siteNumber = sites.length;
        int[] subSites = new int[siteNumber];
        for (int i = 0; i < subSites.length; i++) subSites[i] = sites[i];
        Arrays.sort(subSites);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outputFileS);
            bw.write("Chr\tPos");
            bw.newLine();
            for (int i = 0; i < siteNumber; i++) {
                bw.write("10\t"+String.valueOf(subSites[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    void getTaxaNamList () {
        String folderOldFileS = "M:\\production\\elementMap\\taxa\\source\\282Folders_old.txt"; 
        String folderNovogeneFileS = "M:\\production\\elementMap\\taxa\\source\\282Folders_Nonogene.txt"; 
        String baseFolder = "R:\\WGS\\Zea\\";
        String taxaFqFileS = "M:\\production\\elementMap\\taxa\\taxaFqMap.txt";
        String taxaFqSizeFileS = "M:\\production\\elementMap\\taxa\\taxaFqSize.txt";
        Table t = new Table (folderOldFileS);
        ArrayList<File> oldDirList = new ArrayList();
        ArrayList<String> oldFileNameList = new ArrayList();
        HashMap<String, Double> nameSizeMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            File f = new File(baseFolder, t.content[i][0]);
            oldDirList.add(f);
            File[] fs = f.listFiles();
            fs = IoUtils.listFilesEndsWith(fs, "fastq.gz");
            for (int j = 0; j < fs.length; j++) {
                oldFileNameList.add(fs[j].getName());
                double size = (double)fs[j].length()/1024/1024/1024;
                nameSizeMap.put(fs[j].getName(), size);
            }
        }
        t = new Table (folderNovogeneFileS);
        ArrayList<File> novogeneDirList = new ArrayList();
        ArrayList<String> novogeneFileNameList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            File f = new File(baseFolder, t.content[i][0]);
            novogeneDirList.add(f);
            File[] fs = f.listFiles();
            fs = IoUtils.listFilesEndsWith(fs, "fq.gz");
            for (int j = 0; j < fs.length; j++) {
                novogeneFileNameList.add(fs[j].getName());
                double size = (double)fs[j].length()/1024/1024/1024;
                nameSizeMap.put(fs[j].getName(), size);
            }
        }
        HashSet<String> taxaSet = new HashSet();
        for (int i = 0; i < oldFileNameList.size(); i++) {
            String[] temp = oldFileNameList.get(i).split("_");
            taxaSet.add(temp[4]);
        }
        for (int i = 0; i < novogeneFileNameList.size(); i++) {
            String[] temp = novogeneFileNameList.get(i).split("_");
            taxaSet.add(temp[0]);
        }
        String[] taxa = taxaSet.toArray(new String[taxaSet.size()]);
        Arrays.sort(taxa);
        try {
            BufferedWriter bw1 = IoUtils.getTextWriter(taxaFqFileS);
            BufferedWriter bw2 = IoUtils.getTextWriter(taxaFqSizeFileS);
            bw2.write("Taxa\tNumberOfFile\tTotalSize(Gb)");
            bw2.newLine();
            for (int i = 0; i < taxa.length; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(taxa[i]);
                double fileSize = 0;
                int n = 0;
                for (int j = 0; j < oldFileNameList.size(); j++) {
                    String name = oldFileNameList.get(j);
                    String[] temp = name.split("_");
                    if (temp[4].equals(taxa[i])) {
                        sb.append("\t").append(name);
                        fileSize+=nameSizeMap.get(name);
                        n++;
                    }
                }
                for (int j = 0; j < novogeneFileNameList.size(); j++) {
                    String name = novogeneFileNameList.get(j);
                    String[] temp = name.split("_");
                    if (temp[0].equals(taxa[i])) {
                        sb.append("\t").append(name);
                        fileSize+=nameSizeMap.get(name);
                        n++;
                    }
                }
                bw1.write(sb.toString());
                bw1.newLine();
                bw2.write(taxa[i]+"\t"+String.valueOf(n)+"\t"+String.valueOf(fileSize));
                bw2.newLine();
            }
            bw1.flush();
            bw1.close();
            bw2.flush();
            bw2.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
}
