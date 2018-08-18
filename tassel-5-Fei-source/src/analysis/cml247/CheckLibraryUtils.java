/*
* To change this license header, choose License Headers in Project Properties.
* To change this template file, choose Tools | Templates
* and open the template in the editor.
*/

package analysis.cml247;

import format.Fasta;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.TreeSet;
import net.maizegenetics.util.MultiMemberGZIPInputStream;

/**
 *
 * @author fl262
 */
public class CheckLibraryUtils {
    
    public CheckLibraryUtils() {
        
    }
    
    public void mkRatioTable (String checkSamFileS, String checkRatioTable) {
        try {
            TreeSet<String> laneSet = new TreeSet();
            BufferedReader br = new BufferedReader (new FileReader(checkSamFileS), 65536);
            String temp;
            String[] tem;
            while ((temp = br.readLine())!=null) {
                if (temp.startsWith("@PG"))break;
            }
            while ((temp = br.readLine()) != null) {
                tem = temp.split("\t");
                String[] te = tem[0].split("_");
                String laneName = "";
                for (int i = 0; i < te.length-1; i++) {
                    laneName = laneName+te[i]+"_";
                }
                laneName = laneName.substring(0, laneName.length()-1);
                laneSet.add(laneName);
            }
            String[] lanes = laneSet.toArray(new String[laneSet.size()]);
            Arrays.sort(lanes);
            double[][] ratio = new double[lanes.length][2];
            for (int i = 0; i < ratio.length; i++) {
                for (int j = 0; j < 2; j++) {
                    ratio[i][j] = 0;
                }
            }
            br = new BufferedReader (new FileReader(checkSamFileS), 65536);
            while ((temp = br.readLine())!=null) {
                if (temp.startsWith("@PG"))break;
            }
            while ((temp = br.readLine()) != null) {
                tem = temp.split("\t");
                String[] te = tem[0].split("_");
                String laneName = "";
                for (int i = 0; i < te.length-1; i++) {
                    laneName = laneName+te[i]+"_";
                }
                laneName = laneName.substring(0, laneName.length()-1);
                int index = Arrays.binarySearch(lanes, laneName);
                if (tem[2].startsWith("maize")) {
                    ratio[index][0]++;
                }
                else if (tem[2].startsWith("sorghum")) {
                    ratio[index][1]++;
                } 
            }
            BufferedWriter bw = new BufferedWriter (new FileWriter(checkRatioTable), 65536);
            bw.write("Lane\tMaizePercentage\tSorghumPercentage");
            bw.newLine();
            for (int i = 0; i < lanes.length; i++) {
                double sum = ratio[i][0]+ratio[i][1];
                bw.write(lanes[i]+"\t"+String.valueOf(ratio[i][0]/sum)+"\t"+String.valueOf(ratio[i][1]/sum));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mergeQuery (String queryDirS, String queryFileS) {
        File[] fs = new File(queryDirS).listFiles();
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(queryFileS));
            for (int i = 0; i < fs.length; i++) {
                BufferedReader br = new BufferedReader (new FileReader(fs[i]), 65536);
                String temp;
                while ((temp = br.readLine()) != null) {
                    bw.write(temp);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkPacBioQuery (String pacBioDirS, String queryDirS) {
        File[] fs = new File(pacBioDirS).listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith("fasta");
            }
        });
        int cnt = 0;
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(queryDirS+"/pacBio_cml247.txt"), 65536);
            for (int i = 0; i < fs.length; i++) {
                Fasta f = new Fasta(fs[i].getAbsolutePath());
                for (int j = 0; j < f.getSeqNumber(); j++) {
                    bw.write(">PacBio_"+String.valueOf(cnt));
                    bw.newLine();
                    int length = 150;
                    if (f.getSeq(j).length() < length) length = f.getSeq(j).length();
                    bw.write(f.getSeq(j).substring(0, length));
                    bw.newLine();
                    cnt++;
                }
                if (cnt>20000) break;
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    public void mkQueries (String fastqDir, String queryDirS) {
        File[] fs = new File(fastqDir).listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                if (name.endsWith("gz")) return true;
                else if (name.endsWith("txt")) return true;
                return false;
            }
        });
        for (int i = 0; i < fs.length; i++) {
            System.out.println(fs[i].getAbsoluteFile());
            try {
                String outfileS = queryDirS+"/"+fs[i].getName()+".txt";
                BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
                BufferedReader br = null;
                if (fs[i].getName().endsWith(".gz")) {
                    br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(fs[i]))));
                } else {
                    br = new BufferedReader(new FileReader(fs[i]), 65536);
                }
                for (int j = 0; j < 100000; j++) br.readLine();
                for (int j = 0; j < 2000; j++) {
                    br.readLine();
                    bw.write(">"+fs[i].getAbsolutePath()+"_"+String.valueOf(j));
                    bw.newLine();
                    bw.write(br.readLine());
                    bw.newLine();
                    br.readLine();br.readLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
    }
    
    public void checkSeqByBlast(String fastqDir, String queryFileS, String dataBase, String alignmentFileS, String resultFileS) {
        File[] fs = new File(fastqDir).listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith("gz");
            }
        });
        for (int i = 0; i < fs.length; i++) {
            System.out.println(fs[i].getAbsolutePath());
            try {
                BufferedWriter bw = new BufferedWriter (new FileWriter(queryFileS), 65536);
                BufferedReader br = null;
                if (fs[i].getName().endsWith(".gz")) {
                    br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(fs[i]))));
                } else {
                    br = new BufferedReader(new FileReader(fs[i]), 65536);
                }
                for (int j = 0; j < 10000; j++) br.readLine();
                for (int j = 0; j < 2000; j++) {
                    br.readLine();
                    bw.write(">"+String.valueOf(j));
                    bw.newLine();
                    bw.write(br.readLine());
                    bw.newLine();
                    br.readLine();br.readLine();
                }
                bw.flush();
                bw.close();
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            try {
                Runtime run = Runtime.getRuntime();
                String cmd = "cmd /c blastn -evalue 1e-60 -max_target_seqs 1 -num_alignments 1 -num_threads 3 -outfmt 6 -query " + queryFileS + " -db " + dataBase + " -out " + alignmentFileS;
                System.out.println(cmd);
                Process p = run.exec(cmd);
                p.waitFor();
                System.out.println("Alignment is made at " + alignmentFileS);
            }
            catch (Exception e) {
                System.out.println(e.toString());
                System.exit(1);
            }
        }
    }
    
    public void mergeReference(String maizeRef, String sorghumRef, String mergeRef) {
        Fasta maize = new Fasta(maizeRef);
        Fasta sorghum = new Fasta(sorghumRef);
        for (int i = 0; i < maize.getSeqNumber(); i++) {
            maize.setName("maize-"+maize.getName(i), i);
        }
        for (int i = 0; i < sorghum.getSeqNumber(); i++) {
            sorghum.setName("sorghum-"+sorghum.getName(i), i);
        }
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(mergeRef), 65536);
            for (int i = 0; i < maize.getSeqNumber(); i++) {
                bw.write(">"+maize.getName(i));
                bw.newLine();
                bw.write(maize.getSeq(i));
                bw.newLine();
            }
            for (int i = 0; i < sorghum.getSeqNumber(); i++) {
                bw.write(">"+sorghum.getName(i));
                bw.newLine();
                bw.write(sorghum.getSeq(i));
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
