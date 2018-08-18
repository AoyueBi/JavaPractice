/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.gbs;

import java.io.BufferedWriter;
import java.io.File;
import net.maizegenetics.dna.map.TagGWASMap;
import net.maizegenetics.dna.map.TagGWASMapInfo;
import net.maizegenetics.dna.map.TagsOnGeneticMap;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import utils.IoUtils;

/**
 *
 * @author Fei Lu
 */
public class PanAnchorUtils {
    
    public PanAnchorUtils () {
        
    }
    
    public void resolutionEstimation () {
        String tagMapFileS = "M:\\production\\geneticMapping\\tagMap\\tagGWASMap.h5";
        String v1RigidFileS = "M:\\production\\panAnchorV2\\v1\\v1.panAnchor.rigid.togm.txt";
        String v1RelaxedFileS = "M:\\production\\panAnchorV2\\v1\\v1.panAnchor.relaxed.togm.txt";
        String v2RigidFileS = "M:\\production\\panAnchorV2\\rigid\\v2.panAnchor.rigid.togm.txt";
        String v2RelaxedFileS = "M:\\production\\panAnchorV2\\relaxed\\v2.panAnchor.relaxed.togm.txt";
        String resolutionFileS = "M:\\production\\panAnchorV2\\resolution\\resolution.txt";
        int[] accuracyCut = {10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000, 20000000, 50000000};
        double[] proportion = new double[accuracyCut.length];
        TagGWASMap tgm = new TagGWASMap (tagMapFileS);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(resolutionFileS);
            bw.write("FileName\tAnchorNumber");
            for (int i = 0; i < accuracyCut.length; i++) {
                bw.write("\t<"+String.valueOf(accuracyCut[i]/1000)+ " kb");
            }
            bw.newLine();
            TagsOnGeneticMap togm = new TagsOnGeneticMap (v1RigidFileS, FilePacking.Text);
            int[] v1Rigid = this.getActualDistance(tgm, togm);
            proportion = this.getProportion(accuracyCut, v1Rigid);
            bw.write(new File(v1RigidFileS).getName()+"\t"+String.valueOf(togm.getTagCount()));
            for (int i = 0; i < proportion.length; i++) {
                bw.write("\t"+String.valueOf(proportion[i]));
            }
            bw.newLine();
            togm = new TagsOnGeneticMap (v2RigidFileS, FilePacking.Text);
            int[] v2Rigid = this.getActualDistance(tgm, togm);
            proportion = this.getProportion(accuracyCut, v2Rigid);
            bw.write(new File(v2RigidFileS).getName()+"\t"+String.valueOf(togm.getTagCount()));
            for (int i = 0; i < proportion.length; i++) {
                bw.write("\t"+String.valueOf(proportion[i]));
            }
            bw.newLine();
            togm = new TagsOnGeneticMap (v1RelaxedFileS, FilePacking.Text);
            int[] v1Relaxed = this.getActualDistance(tgm, togm);
            proportion = this.getProportion(accuracyCut, v1Relaxed);
            bw.write(new File(v1RelaxedFileS).getName()+"\t"+String.valueOf(togm.getTagCount()));
            for (int i = 0; i < proportion.length; i++) {
                bw.write("\t"+String.valueOf(proportion[i]));
            }
            bw.newLine();
            togm = new TagsOnGeneticMap (v2RelaxedFileS, FilePacking.Text);
            int[] v2Relaxed = this.getActualDistance(tgm, togm);
            proportion = this.getProportion(accuracyCut, v2Relaxed);
            bw.write(new File(v2RelaxedFileS).getName()+"\t"+String.valueOf(togm.getTagCount()));
            for (int i = 0; i < proportion.length; i++) {
                bw.write("\t"+String.valueOf(proportion[i]));
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private double[] getProportion (int[] accuracyCut, int[] distance) {
        double[] proportion = new double[accuracyCut.length];
        for (int i = 0; i < distance.length; i++) {
            for (int j = 0; j < accuracyCut.length; j++) {
                if (distance[i] < accuracyCut[j]) proportion[j]++;
            }
        }
        for (int i = 0; i < proportion.length; i++) {
            proportion[i] = (double)proportion[i]/distance.length;
        }
        return proportion;
    }
    
    private int[] getActualDistance (TagGWASMap tgm, TagsOnGeneticMap togm) {
        int size = 30000;
        int[] distance = new int[size];
        int cnt = 0;
        for (int i = 0; i < tgm.getTagCount(); i++) {
            TagGWASMapInfo tgmi = tgm.getTagGWASMapInfo(i);
            if (!tgmi.isUniqueRef()) continue;
            long[] t = tgm.getTag(i);
            int index = togm.getTagIndex(t);
            if (index < 0) continue;
            int d = 1000000000;
            if (togm.getGChr(index) == tgmi.pChr) d = Math.abs(togm.getGPos(index) - tgmi.pPos);
            else d=d+Math.abs(togm.getGPos(index) - tgmi.pPos);
            distance[cnt] = d;
            cnt++;
            if (cnt == size) break;
        }
        return distance;
    }
}
