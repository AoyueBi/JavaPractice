/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.collaboration;

import format.ShortreadAlignment;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.Arrays;
import java.util.HashSet;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.map.TagsOnGeneticMap;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class Jianbing {
    
    public Jianbing () {
        //this.filterAlignment();
        this.getSubAnchorFileS();
    }
    
    public void getSubAnchorFileS () {
        String anchorFileS = "M:\\production\\panAnchorV2\\v1\\v1.panAnchor.rigid.togm.txt";
        String saFileS = "M:\\collaboration\\jianbing\\novelgene\\anchorOnNovel.filter.sam";
        String subAnchorFileS = "M:\\collaboration\\jianbing\\novelgene\\subAnchor.txt";
        TagsOnGeneticMap togm = new TagsOnGeneticMap(anchorFileS, FilePacking.Text);
        HashSet<Integer> indexSet = new HashSet();
        try {
            BufferedReader br = IoUtils.getTextReader(saFileS);
            String temp;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("@")) continue;
                String[] tem = temp.substring(0, 20).split("\t");
                indexSet.add(Integer.valueOf(tem[0]));
            }
            Integer[] index = indexSet.toArray(new Integer[indexSet.size()]);
            Arrays.sort(index);
            BufferedWriter bw = IoUtils.getTextWriter(subAnchorFileS);
            bw.write("AnchorIndex\tSeq\tGChr\tGPos");
            bw.newLine();
            for (int i = 0; i < index.length; i++) {
                bw.write(String.valueOf(index[i])+"\t"+BaseEncoder.getSequenceFromLong(togm.getTag(i))+"\t"+String.valueOf(togm.getGChr(index[i]))+"\t"+String.valueOf(togm.getGPos(index[i])));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void filterAlignment () {
        String samFileS = "M:\\collaboration\\jianbing\\novelgene\\anchorOnNovel.sam";
        String saFileS = "M:\\collaboration\\jianbing\\novelgene\\anchorOnNovel.filter.sam";
        try {
            BufferedReader br = IoUtils.getTextReader(samFileS);
            BufferedWriter bw = IoUtils.getTextWriter(saFileS);
            String temp;
            while ((temp = br.readLine()) != null) {
                if (temp.substring(0, 20).contains("*")) continue;
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
    }
    
    
    
}
