/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.maf;

import format.Bins;
import format.Table;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author Fei Lu
 */
class InvariantSiteDistribution {
    
    public InvariantSiteDistribution () {
        //this.mkStatistics();
        //this.mkChr10Dis();
    }
    
    public void mkChr10Dis () {
        int binSize = 100000;
        String infoFileS = "E:\\Database\\InfoFile\\ChrLenCentPosi_agpV3.txt";
        String depthFileS = "M:\\production\\maf\\wgs\\invariantSite\\invariant_hmp31_ad_10.bin";
        String outfileS = "E:\\Research\\wgs_maf\\invariantSite_distribution\\chr10_dis.txt";
        Table t  = new Table(infoFileS);
        int chrLength = Integer.valueOf(t.content[9][1]);
        AlleleDepth ad  = new AlleleDepth(depthFileS, IOFileFormat.Binary);
        int[] pos = new int[ad.getSiteNumber()];
        double[] v = new double[ad.getSiteNumber()];
        for (int i = 0; i < ad.getSiteNumber(); i++) {
            pos[i] = ad.getPosition(i);
            v[i] = 1;
        }
        Bins b = new Bins (1, chrLength,binSize, pos, v);
        double[] vs = new double[b.getBinNum()];
        for (int i = 0; i < vs.length; i++) {
            vs[i] = b.getNumValues(i);
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("BinStart\tNumberOfInvariant");
            bw.newLine();
            for (int i = 0; i < b.getBinNum(); i++) {
                bw.write(String.valueOf(b.getBinStart(i))+"\t"+String.valueOf(b.getNumValues(i)));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkStatistics () {
        String infileDirS = "M:\\production\\maf\\wgs\\invariantSite\\";
        String outfileS = "E:\\Research\\wgs_maf\\invariantSite_distribution\\statistics.txt";
        File[] fs = new File(infileDirS).listFiles();
        AlleleDepth ad  = null;
        int[] size = new int[fs.length];
        int[] chr = new int[fs.length];
        long total = 0;
        for (int i = 0; i < fs.length; i++) {
            ad = new AlleleDepth(fs[i].getAbsolutePath(), IOFileFormat.Binary);
            size[i] = ad.getSiteNumber();
            chr[i] = ad.getChromosome(0);
            total += size[i];
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Total:\t" + String.valueOf(total));
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                bw.write("chromosome"+String.valueOf(chr[i])+":\t"+String.valueOf(size[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    public static AlleleDepth[] getInvariantSiteByChrome () {
        String infileDirS = "M:\\production\\maf\\wgs\\invariantSite";
        File[] fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        ArrayList<AlleleDepth> adList = new ArrayList();
        fsList.parallelStream().forEach(e -> {
            adList.add(new AlleleDepth(e.getAbsolutePath(), IOFileFormat.Binary));
        });
        AlleleDepth[] ads = adList.toArray(new AlleleDepth[adList.size()]);
        Arrays.sort(ads);
        return ads;
    }
}
