/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.hmp3sv;

import format.Fasta;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class Hmp3SVGo {
    int[] chrLength = {301476924, 237917468, 232245527, 242062272, 217959525, 169407836, 176826311, 175377492, 157038028, 149632204};
    
    public Hmp3SVGo () {
        this.depthPipe();
    }
    
    public void depthPipe () {
        this.pileup();
    }
    
    public void pileup () {
        int coreNumber = 23;
        int chromosome = 6;
        int binSize = 100;
        int regionSize = 10000000;
        String infileDirS = "/workdir/mingh/source/";
        String outDirS = "/workdir/mingh/depthInBin/";
        String pileupDirS = "/workdir/mingh/pileup/";
        String reference = "/workdir/mingh/ref/maize3.fa";
        File[] fs = IoUtils.listFilesEndsWith(new File(infileDirS).listFiles(), ".bam");
        new File(outDirS).mkdir();
        new File(pileupDirS).mkdir();
        int left = chrLength[chromosome-1]%regionSize;
        int regionNumber;
        if (left == 0) regionNumber = chrLength[chromosome-1]/regionSize;
        else regionNumber = chrLength[chromosome-1]/regionSize+1;
        int[] start = new int[regionNumber];
        int[] end = new int[regionNumber];
        if (left == 0) {
            for (int i = 0; i < regionNumber; i++) {
                start[i] = i*regionSize+1;
                end[i] = (i+1)*regionSize;
            }
        }
        else {
            for (int i = 0; i < regionNumber-1; i++) {
                start[i] = i*regionSize+1;
                end[i] = (i+1)*regionSize;
            }
            start[regionNumber-1] = chrLength[chromosome-1] - left + 1;
            end[regionNumber-1] = chrLength[chromosome-1];
        }
        
        left = chrLength[chromosome-1]%binSize;
        int binNumber = 0;
        if (left == 0) binNumber = chrLength[chromosome-1]/binSize;
        else binNumber = chrLength[chromosome-1]/binSize+1;
        int[] binStart = new int[binNumber];
        int[] binEnd = new int[binNumber];
        if (left == 0) {
            for (int i = 0; i < binNumber; i++) {
                binStart[i] = i*binSize+1;
                binEnd[i] = (i+1)*binSize;
            }
        }
        else {
            for (int i = 0; i < binNumber-1; i++) {
                binStart[i] = i*binSize+1;
                binEnd[i] = (i+1)*binSize;
            }
            binStart[binNumber-1] = chrLength[chromosome-1] - left + 1;
            binEnd[binNumber-1] = chrLength[chromosome-1];
        }
        
        
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        String[] taxaName = nameSet.toArray(new String[nameSet.size()]);
        
        int[][] depth = new int[binNumber][taxaName.length];
        Arrays.sort(taxaName);
        TIntArrayList[] indexList = new TIntArrayList[taxaName.length];
        int[][] taxaToFile = new int[taxaName.length][];
        for (int i = 0; i < indexList.length; i++) indexList[i] = new TIntArrayList();
        for (int i = 0; i < fs.length; i++) {
            String name = fs[i].getName().split("_")[0];
            int index = Arrays.binarySearch(taxaName, name);
            indexList[index].add(i);
        }
        for (int i = 0; i < indexList.length; i++) {
            taxaToFile[i] = indexList[i].toArray();
            Arrays.sort(taxaToFile[i]);
        }
        for (int i = 0; i < regionNumber; i++) {
            String[] cmd = new String[taxaName.length];
            for (int j = 0; j < taxaName.length; j++) {
                String outfileS = new File(pileupDirS, taxaName[j]+".pileup.txt").getAbsolutePath();
                StringBuilder sb = new StringBuilder();
                sb.append("samtools mpileup -A -B -q 30 -f ").append(reference);
                for (int k = 0; k < taxaToFile[j].length; k++) {
                    sb.append(" ").append(fs[taxaToFile[j][k]].getAbsolutePath());
                }
                sb.append(" -r ").append(chromosome).append(":").append(start[i]).append("-").append(end[i]).append(" -o ").append(outfileS);
                cmd[j] = sb.toString();
            }
            int cnt = 0;
            for (int j = 0; j < cmd.length; j+=coreNumber) {
                int endIndex = j+coreNumber;
                if (endIndex > cmd.length) endIndex = cmd.length;
                ArrayList<String> cmdList = new ArrayList();
                for (int k = j; k < endIndex; k++) {
                    cmdList.add(cmd[k]);
                }
                cmdList.parallelStream().forEach(element -> {
                    try {
                        Runtime run = Runtime.getRuntime();
                        Process p = run.exec(element);
                        p.waitFor();
                    }
                    catch (Exception e) {
                        e.printStackTrace();
                    }
                });
                System.out.println(endIndex + " " + String.valueOf((double)(endIndex/cmd.length)) + " taxa finished pileup in this regon.");
                cnt++;
            }
            System.out.println("Pileup is finished for region:\t" +String.valueOf(chromosome)+"\t"+String.valueOf(start[i])+"\t"+String.valueOf(end[i]));
            File[] piles = new File(pileupDirS).listFiles();           
            List<File> pileList = Arrays.asList(piles);
            pileList.parallelStream().forEach(element -> {
                int taxaIndex = Arrays.binarySearch(taxaName, element.getName().replaceFirst(".pileup.txt", ""));
                try {
                    BufferedReader br = IoUtils.getTextReader(element.getAbsolutePath());
                    String temp = br.readLine();
                    String[] tem = temp.split("\t", -1);
                    int n = tem.length/3-1;
                    br = IoUtils.getTextReader(element.getAbsolutePath());
                    while ((temp = br.readLine()) != null) {
                        tem = temp.split("\t");
                        int pos = Integer.valueOf(tem[1]);
                        int binIndex = Arrays.binarySearch(binStart, pos);
                        if (binIndex < 0) {
                            binIndex = -binIndex-2;
                        }
                        int d = 0;
                        for (int j = 0; j < n; j++) {
                            d+=Integer.valueOf(tem[j*3+3]);
                        }
                        depth[binIndex][taxaIndex]+=d;
                    }
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            });
        }
        
        String outfileS = FStringUtils.getNDigitNumber(2, chromosome)+".binDepth.txt";
        outfileS = new File(outDirS, outfileS).getAbsolutePath();
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Chromosome\tBinStart\tBinEnd");
            for (int i = 0; i < taxaName.length; i++) {
                bw.write("\t"+taxaName[i]);
            }
            bw.newLine();
            for (int i = 0; i < binStart.length; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(chromosome).append("\t").append(binStart[i]).append("\t").append(binEnd[i]);
                for (int j = 0; j < depth[i].length; j++) {
                    sb.append("\t").append(depth[i][j]);
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        } 
    }
    
    public static void main (String[] args) {
        new Hmp3SVGo();
    }
    
}
