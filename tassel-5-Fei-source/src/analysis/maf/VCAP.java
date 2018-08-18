/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maf;

import format.Range;
import format.Ranges;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import utils.FStringUtils;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class VCAP {
    
    public VCAP () {
        this.mkAllHmpPositions();
    }
    
    public static void toPositionListInJson (int[] chr, int[] pos, String outfileS, IOFileFormat format) {
        try {
            BufferedWriter bw = null;
            if (format == IOFileFormat.Text) {
                bw = IoUtils.getTextWriter(outfileS);
            }
            else if (format == IOFileFormat.TextGzip) {
                bw = IoUtils.getTextGzipWriter(outfileS);
            }
            else {
                throw new UnsupportedOperationException("Not supported yet.");
            }        
            bw.newLine();
            bw.write("{\n");
            bw.write("    \"PositionList\":[\n");
            for (int i = 0; i < chr.length; i++) {
                bw.write("        {\n");
                bw.write("            \"chr\":{\n");
                StringBuilder sb = new StringBuilder();
                sb.append("                \"name\":\"").append(chr[i]).append("\"\n");
                bw.write(sb.toString());
                bw.write("            },\n");
                sb = new StringBuilder();
                sb.append("            \"position\":").append(pos[i]).append(",\n");
                bw.write(sb.toString());
                bw.write("            \"strand\":\"+\"\n");
                bw.write("        }");
                if (i == chr.length-1) {
                    bw.newLine();
                    bw.write("    ]\n");
                }
                else {
                    bw.write(",\n");
                }
            }
            bw.write("}");
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public static void toPositionListInJson (Ranges hmpRs, int[] chr, int[] pos, String outfileS) {
        TIntArrayList indexList = new TIntArrayList();
        hmpRs.sortByStartPosition();
        for (int i = 0; i < chr.length; i++) {
            if (!hmpRs.isInRanges(chr[i], pos[i])) continue;
            indexList.add(i);
        }
        int[] index = indexList.toArray();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.newLine();
            bw.write("{\n");
            bw.write("    \"PositionList\":[\n");
            for (int i = 0; i < index.length; i++) {
                bw.write("        {\n");
                bw.write("            \"chr\":{\n");
                StringBuilder sb = new StringBuilder();
                sb.append("                \"name\":\"").append(chr[index[i]]).append("\"\n");
                bw.write(sb.toString());
                bw.write("            },\n");
                sb = new StringBuilder();
                sb.append("            \"position\":").append(pos[index[i]]).append(",\n");
                bw.write(sb.toString());
                bw.write("            \"strand\":\"+\"\n");
                bw.write("        }");
                if (i == index.length-1) {
                    bw.newLine();
                    bw.write("    ]\n");
                }
                else {
                    bw.write(",\n");
                }
            }
            bw.write("}");
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void mkAllHmpPositions () {
        String infileS = "/workdir/fl262/vcap/genotypes/unzip/namrils_projected_hmp31_MAF02mnCnt2500.hmp.txt";
        String outfileS = "/workdir/fl262/vcap/genotypes/unzip/namrils_projected_hmp31_MAF02mnCnt2500.hmp.pos.range2.txt";
        ArrayList<Range> rList = new ArrayList();
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            String temp = br.readLine();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                temp = temp.substring(0, 100);
                List<String> l = FStringUtils.fastSplit(temp, "\t");
                int chr = Integer.valueOf(l.get(2));
                int pos = Integer.valueOf(l.get(3));
                Range r = new Range (chr, pos, pos+1);
                rList.add(r);
                cnt++;
                if (cnt%1000000 == 0) System.out.println(cnt);
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        Ranges rs = new Ranges(rList, new File(infileS).getName());
        rs.writeFile(outfileS, IOFileFormat.Text);
    }
    
}
