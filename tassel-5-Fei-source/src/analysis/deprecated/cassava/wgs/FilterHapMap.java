/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.cassava.wgs;

import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.List;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class FilterHapMap {
    
    public FilterHapMap () {
        this.filterByDepthIndel();
    }
    
    public void filterByDepthIndel () {
//        String inputHapMapDirS = "M:\\production\\cassava\\hapmap\\raw\\";
//        String validSiteDirS = "M:\\pipelineTest\\cassava\\wgs\\badRefSite\\validSite\\";
//        String outputHapMapDirS = "M:\\production\\cassava\\hapmap\\stage_01\\";
        String inputHapMapDirS = "/workdir/fl262/raw/";
        String validSiteDirS = "/workdir/fl262/validSite/";
        String outputHapMapDirS = "/workdir/fl262/stage_01/";
        File[] fs = new File (inputHapMapDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String infileS = new File(validSiteDirS, f.getName().replaceFirst(".VCF.txt", ".validSite.txt")).getAbsolutePath();
            String outfileS = new File(outputHapMapDirS, f.getName().replaceFirst(".VCF.txt", ".stage01.VCF.txt")).getAbsolutePath();
            TIntArrayList siteList = new TIntArrayList();
            try {
                BufferedReader br = IoUtils.getTextReader(infileS);
                String temp = br.readLine();
                while ((temp = br.readLine()) != null) {
                    siteList.add(Integer.valueOf(temp.split("\t")[1]));
                }
                br.close();
                int[] sites = siteList.toArray();
                br = IoUtils.getTextReader(f.getAbsolutePath());
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                bw.write(br.readLine());
                bw.newLine();
                while ((temp = br.readLine()) != null) {
                    int index = temp.indexOf("GT");
                    String[] tem = temp.substring(0, index).split("\t");
                    int depth = Integer.valueOf(tem[4].split(";")[0].replaceFirst("DP=", ""));
                    if (depth < 7000) continue;
                    if (depth > 15000) continue;
                    if (Arrays.binarySearch(sites, Integer.valueOf(tem[1])) < 0) continue;
                    tem = tem[4].split(";");
                    tem = tem[5].replaceFirst("DI=", "").split(",");
                    int cnt = 0;
                    for (int i = 0; i < tem.length; i++) cnt+=Integer.valueOf(tem[i]);
                    if (cnt>1) continue;
                    bw.write(temp);
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            int[] sites = siteList.toArray();
            
        });
    }
}
