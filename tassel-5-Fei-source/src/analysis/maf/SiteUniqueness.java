/*
* To change this license header, choose License Headers in Project Properties.
* To change this template file, choose Tools | Templates
* and open the template in the editor.
*/
package analysis.maf;

import format.Bins;
import format.Fasta;
import format.GeneFeature;
import format.Range;
import format.RangeAttribute;
import format.Ranges;
import format.Table;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;
import gnu.trove.map.hash.TIntFloatHashMap;
import gnu.trove.map.hash.TIntIntHashMap;
import graphcis.r.DensityPlot;
import graphcis.r.Histogram;
import graphcis.r.ScatterPlot;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import net.maizegenetics.dna.BaseEncoder;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import utils.Benchmark;
import utils.FArrayUtils;
import utils.FStringUtils;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class SiteUniqueness {
    HashMap<Byte, Byte> ascIIByteMap = null;

    public SiteUniqueness (String referenceGenome, String testGenome, String resultDirS) {

        this.buildUniquenessTable(referenceGenome, testGenome, resultDirS);
    }

    public SiteUniqueness () {
        this.buildUniquenessTableBatch();
        //this.combineUniqueness();
        //this.analysisPipe();
        //this.vcapPipe();
        //this.vcapPipe2();
    }
    
    public void vcapPipe2 () {
        //this.mkTwoGroupPositionList();
        //this.mkVCAPScript2();
        //this.mkSNPNumberTwoGroup();
        this.summerizeVCAPTwoGroup();
    }
    
    private void summerizeVCAPTwoGroup () {
        String inDirS = "M:\\production\\maf\\annotations\\siteUniqueness\\vcap\\positionList\\twoGroups\\results\\";
        String snpCountFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\vcap\\positionList\\twoGroups\\snpCount.txt";
        String totalVarianceFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\vcap\\positionList\\twoGroups\\totalVariance.txt";
        String heritabilityFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\vcap\\positionList\\twoGroups\\heritability.txt";
        String enrichmentFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\vcap\\positionList\\twoGroups\\enrichment.txt";
        Table t = new Table(snpCountFileS);
        int[][] snpCount = new int[t.getRowNumber()][t.getColumnNumber()];
        String[] rangeName = new String[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            rangeName[i] = t.content[i][0];
            snpCount[i][0] = t.getIntValue(i, 1);
            snpCount[i][1] = t.getIntValue(i, 2);
        }
        String infileDirS = new File(inDirS, rangeName[0]).getAbsolutePath();
        File[] fs = IoUtils.listFilesEndsWith(new File(infileDirS).listFiles(), "reml");
        String[] fileName = new String[fs.length];
        for (int i = 0; i < fileName.length; i++) {
            fileName[i] = fs[i].getName();
        }
        double[][] variance = new double[fileName.length][rangeName.length];
        double[][] h2 = new double[fileName.length][rangeName.length];
        double[][] enrich = new double[fileName.length][rangeName.length];
        for (int i = 0; i < fileName.length; i++) {
            for (int j = 0; j < rangeName.length; j++) {
                String infileS = new File(inDirS, rangeName[j]).getAbsolutePath();
                infileS = new File(infileS, fileName[i]).getAbsolutePath();
                try {
                    BufferedReader br = IoUtils.getTextReader(infileS);
                    String temp;
                    while (!(temp = br.readLine()).startsWith("Component Heritability")) {}
                    double v1 = Double.valueOf(br.readLine().split(" ")[1]);
                    double v2 = Double.valueOf(br.readLine().split(" ")[2]);
                    double sum = v1+v2;
                    variance[i][j] = sum;
                    h2[i][j] = v1/sum;
                    enrich[i][j] = (v1/snpCount[j][0])/(sum/(snpCount[j][0]+snpCount[j][1]));
                    br.close();
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(totalVarianceFileS);
            bw.write("Traint");
            for (int i = 0; i < rangeName.length; i++) {
                bw.write("\t"+rangeName[i]);
            }
            bw.newLine();
            for (int i = 0; i < fileName.length; i++) {
                bw.write(fileName[i].replaceFirst(".reml", ""));
                for (int j = 0; j < rangeName.length; j++) {
                    bw.write("\t"+String.valueOf(variance[i][j]));
                }
                bw.newLine();
            }
            bw.flush();bw.close();
            bw = IoUtils.getTextWriter(heritabilityFileS);
            bw.write("Traint");
            for (int i = 0; i < rangeName.length; i++) {
                bw.write("\t"+rangeName[i]);
            }
            bw.newLine();
            for (int i = 0; i < fileName.length; i++) {
                bw.write(fileName[i].replaceFirst(".reml", ""));
                for (int j = 0; j < rangeName.length; j++) {
                    bw.write("\t"+String.valueOf(h2[i][j]));
                }
                bw.newLine();
            }
            bw.flush();bw.close();
            bw = IoUtils.getTextWriter(enrichmentFileS);
            bw.write("Traint");
            for (int i = 0; i < rangeName.length; i++) {
                bw.write("\t"+rangeName[i]);
            }
            bw.newLine();
            for (int i = 0; i < fileName.length; i++) {
                bw.write(fileName[i].replaceFirst(".reml", ""));
                for (int j = 0; j < rangeName.length; j++) {
                    bw.write("\t"+String.valueOf(enrich[i][j]));
                }
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkSNPNumberTwoGroup () {
        double[] scoreThresh = {0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5};
        String inDirS = "M:\\production\\maf\\annotations\\siteUniqueness\\vcap\\positionList\\twoGroups\\positionList\\";
        String outfileS = "M:\\production\\maf\\annotations\\siteUniqueness\\vcap\\positionList\\twoGroups\\snpCount.txt";
        String[] dirName = new String[scoreThresh.length-1];
        int[] k1 = new int[scoreThresh.length-1];
        int[] k2 = new int[scoreThresh.length-1];
        ArrayList<Integer> indexList = new ArrayList();
        for (int i = 0; i < scoreThresh.length-1; i++) indexList.add(i);
        indexList.parallelStream().forEach(index -> {
            dirName[index] = "range_"+String.valueOf(scoreThresh[index])+"_"+String.valueOf(scoreThresh[index+1]);
            String dir = new File(inDirS, dirName[index]).getAbsolutePath();
            String k1FileS = new File(dir, "case.json.gz").getAbsolutePath();
            String k2FileS = new File(dir, "control.json.gz").getAbsolutePath();
            try {
                BufferedReader br = IoUtils.getTextGzipReader(k1FileS);
                int cnt = 0;
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    if (temp.contains("chr")) cnt++;
                }
                br.close();
                k1[index] = cnt;
                br = IoUtils.getTextGzipReader(k2FileS);
                System.out.println("finished " + k1FileS);
                cnt = 0;
                temp = null;
                while ((temp = br.readLine()) != null) {
                    if (temp.contains("chr")) cnt++;
                }
                br.close();
                k2[index] = cnt;
                System.out.println("finished " + k2FileS);
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Range\tK1Count\tK2Count");
            bw.newLine();
            for (int i = 0; i < dirName.length; i++) {
                bw.write(dirName[i]+"\t"+String.valueOf(k1[i])+"\t"+String.valueOf(k2[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void mkVCAPScript2 () {
        String runFileS = "/local/workdir/fl262/vcap/TASSEL5/run_pipeline.pl";
        String hmpFileS = "/local/workdir/fl262/vcap/genotypes/namrils_projected_hmp31_MAF02mnCnt2500.hmp.txt.gz.lix";
        String positionListDirS = "/local/workdir/fl262/vcap/positionList/";
        String kinshipDirS = "/local/workdir/fl262/vcap/kinship";
        String scriptFileS = "/local/workdir/fl262/script.pl";
        String resultDirS = "/local/workdir/fl262/vcap/result/";
        String phenotypeDirS = "/local/workdir/fl262/vcap/phenotypes/NAM/";
        String ldakPath = "/programs/LDAK/ldak.4.9.fast";
        new File (kinshipDirS).mkdir();
        ArrayList<String> commandList = new ArrayList();
        File[] fs = new File(positionListDirS).listFiles();
        double[] vs = new double[fs.length];
        for (int i = 0; i < vs.length; i++) {
            vs[i] = Double.valueOf(fs[i].getName().split("_")[1]);
        }
        for (int i = 0; i < fs.length-1; i++) {
            for (int j = i; j < fs.length; j++) {
                if (vs[i] > vs[j]) {
                    double t = vs[i];
                    vs[i] = vs[j];
                    vs[j] = t;
                    File tf = fs[i];
                    fs[i] = fs[j];
                    fs[j] = tf;
                }
            }
        }
        File[] phenofs = IoUtils.listFilesEndsWith(new File(phenotypeDirS).listFiles(), ".txt");
        Arrays.sort(phenofs);
        for (int i = 0; i < fs.length; i++) {
            String kinshipOutDirS = fs[i].getName();
            String resultOutDirS = fs[i].getName(); 
            kinshipOutDirS = new File(kinshipDirS, kinshipOutDirS).getAbsolutePath();
            resultOutDirS = new File(resultDirS, resultOutDirS).getAbsolutePath();
            new File(kinshipOutDirS).mkdir();
            new File(resultOutDirS).mkdir();
            StringBuilder sb = new StringBuilder();
            String k1positionFileS = new File (fs[i], "case.json.gz").getAbsolutePath();
            String k2positionFileS = new File (fs[i], "control.json.gz").getAbsolutePath();
            String k1KinshipFileS = new File(kinshipOutDirS, "k1").getAbsolutePath();
            String k2KinshipFileS = new File(kinshipOutDirS, "k2").getAbsolutePath();
            sb.append("system (\"").append("perl ").append(runFileS).append(" -Xmx120g -importGuess ").append(hmpFileS).append(" -FilterSiteBuilderPlugin -positionList ").append(k1positionFileS);
            sb.append(" -endPlugin -KinshipPlugin -method Normalized_IBS -endPlugin -export ").append(k1KinshipFileS);
            sb.append(" -exportType SqrMatrixBin \");");
            commandList.add(sb.toString());
            sb = new StringBuilder();
            sb.append("system (\"").append("perl ").append(runFileS).append(" -Xmx120g -importGuess ").append(hmpFileS).append(" -FilterSiteBuilderPlugin -positionList ").append(k2positionFileS);
            sb.append(" -endPlugin -KinshipPlugin -method Normalized_IBS -endPlugin -export ").append(k2KinshipFileS);
            sb.append(" -exportType SqrMatrixBin \");");
            commandList.add(sb.toString());
            String kinshipPathFileS = new File(kinshipOutDirS, "kinshipList.txt").getAbsolutePath();
            try {
                BufferedWriter bw = IoUtils.getTextWriter(kinshipPathFileS);
                bw.write(k1KinshipFileS); bw.newLine();
                bw.write(k2KinshipFileS); bw.newLine();
                bw.flush();bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            for (int j = 0; j < phenofs.length; j++) {
                sb = new StringBuilder();
                String outPrefixS = phenofs[j].getName().replaceFirst(".txt", "");
                 
                outPrefixS = new File (resultOutDirS, outPrefixS).getAbsolutePath();
                sb.append("system (\"");
                sb.append(ldakPath).append(" --mem-save YES --reml ").append(outPrefixS).append(" --mgrm ").append(kinshipPathFileS);
                sb.append(" --pheno ").append(phenofs[j].getAbsolutePath()).append(" --kinship-details NO").append("\");");
                commandList.add(sb.toString());
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(scriptFileS);
            for (int i = 0; i < commandList.size(); i++) {
                bw.write(commandList.get(i));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void mkTwoGroupPositionList () {
        double[] scoreThresh = {0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5};
        String infileS  = "M:\\production\\maf\\annotations\\siteUniqueness\\vcap\\hmp31WithScore.txt";
        String positionDirS = "M:\\production\\maf\\annotations\\siteUniqueness\\vcap\\positionList\\twoGroups\\positionList\\";
        new File(positionDirS).mkdir();
        try {
            String[][] outfileS = new String[scoreThresh.length-1][2];
            TIntArrayList[][] chrList = new TIntArrayList[scoreThresh.length-1][2];
            TIntArrayList[][] posList = new TIntArrayList[scoreThresh.length-1][2];
            for (int i = 0; i < scoreThresh.length-1; i++) {
                String outDirS = "range_"+String.valueOf(scoreThresh[i])+"_"+String.valueOf(scoreThresh[i+1]);
                outDirS = new File(positionDirS,outDirS).getAbsolutePath();
                new File(outDirS).mkdir();
                outfileS[i][0] = new File(outDirS, "case.json.gz").getAbsolutePath();
                outfileS[i][1] = new File(outDirS, "control.json.gz").getAbsolutePath();
                for (int j = 0; j < 2; j++) {
                    chrList[i][j] = new TIntArrayList();
                    posList[i][j] = new TIntArrayList();
                }
            }
            BufferedReader br = IoUtils.getTextReader(infileS);
            String temp = br.readLine();
            int cnt = 0;
            while ((temp = br.readLine())!=null) {
                List<String> l = FStringUtils.fastSplit(temp, "\t");
                if (l.get(2).startsWith("N")) continue;
                double value = Double.valueOf(l.get(2));
                int chr = Integer.valueOf(l.get(0));
                int pos = Integer.valueOf(l.get(1));
                for (int i = 0; i < scoreThresh.length-1; i++) {
                    if (value>=scoreThresh[i] && value<scoreThresh[i+1]) {
                        chrList[i][0].add(chr);
                        posList[i][0].add(pos);
                    }
                    else {
                        chrList[i][1].add(chr);
                        posList[i][1].add(pos);
                    }
                }
                cnt++;
                if (cnt%1000000 == 0) System.out.println(cnt);
            }
            br.close();
            for (int i = 0; i < scoreThresh.length-1; i++) {
                for (int j = 0; j < 2; j++) {
                    VCAP.toPositionListInJson(chrList[i][j].toArray(), posList[i][j].toArray(), outfileS[i][j], IOFileFormat.TextGzip);
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void vcapPipe () {
        //this.getSingleCopyAllChr();
        //this.mkPositionList();
        //this.mkHapMapWithUniqueScore2();
        //this.mkPositionListByScore();
        //this.mkVCAPScript();
        //this.summarizeVCAP();
    }
    
    public void summarizeVCAP () {
        String categaryFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\vcap\\positionList\\set03\\kinship\\mkKinship.pl";
        String positionListDirS = "M:\\production\\maf\\annotations\\siteUniqueness\\vcap\\positionList\\set03\\positionList";
        String h2DirS = "M:\\production\\maf\\annotations\\siteUniqueness\\vcap\\positionList\\set03\\results\\";
        String h2SummaryFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\vcap\\positionList\\set03\\h2Summary.txt";
        String enrichFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\vcap\\positionList\\set03\\h2Enrich.txt";
        String snpNumberFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\vcap\\positionList\\set03\\snpNumber.txt";
        HashMap<String, String> kinshipFileMap = new HashMap();
        try {
            BufferedReader br = IoUtils.getTextReader(categaryFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    String[] tem = temp.substring(1, temp.length()).split("\t");
                    kinshipFileMap.put(tem[1], tem[0]);
                }
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        int[] snpNumber = new int[kinshipFileMap.size()];
        Set<Map.Entry<String,String>> set = kinshipFileMap.entrySet();
        set.parallelStream().forEach(e -> {
            int index = Integer.valueOf(e.getKey().replaceFirst("k", ""))-1;
            try {
                String infileS = new File(positionListDirS, e.getValue()).getAbsolutePath();
                BufferedReader br = IoUtils.getTextGzipReader(infileS);
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.contains("chr")) cnt++;
                }
                snpNumber[index] = cnt;
                br.close();
            }
            catch (Exception ee) {
                ee.printStackTrace();
            }
        });
        int totalSNP = 0;
        try {
            BufferedWriter bw = IoUtils.getTextWriter(snpNumberFileS);
            bw.write(snpNumberFileS);
            bw.newLine();
            for (int i = 0; i < snpNumber.length; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(kinshipFileMap.get("k"+String.valueOf(i+1))).append("\t").append(snpNumber[i]);
                totalSNP+=snpNumber[i];
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.write("AllSNP"+"\t"+String.valueOf(totalSNP));
            bw.newLine();
            bw.flush();
            bw.close();
            System.out.println("Total SNP number:\t" + String.valueOf(totalSNP));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
        File[] fs = IoUtils.listFilesEndsWith(new File(h2DirS).listFiles(), "reml");
        Arrays.sort(fs);
        double[][] h2Sub = new double[fs.length][kinshipFileMap.size()];
        double[][] enrichment = new double[fs.length][kinshipFileMap.size()];
        double[] h2 = new double[fs.length];
        String[] traitName = new String[fs.length];
        for (int i = 0; i < fs.length; i++) {
            try {
                BufferedReader br = IoUtils.getTextReader(fs[i].getAbsolutePath());
                String temp = br.readLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("Her_")) {
                        String[] tem = temp.split(" ");
                        h2Sub[i][cnt] = Double.valueOf(tem[1]);
                        cnt++;
                    }
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            h2[i] = 0;
            for (int j = 0; j < h2Sub[i].length; j++) {
                h2[i]+=h2Sub[i][j];
            }
            String name = fs[i].getName().replaceFirst("NAM_", "").replaceFirst("_multiblup.txt", "");
            traitName[i] = name;
            double mean = h2[i]/totalSNP;
            for (int j = 0; j < h2Sub[i].length; j++) {
                enrichment[i][j] = h2Sub[i][j]/snpNumber[j]/mean;
                h2Sub[i][j] = h2Sub[i][j]/h2[i];
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(h2SummaryFileS);
            bw.write("Trait\th2");
            for (int i = 0; i < snpNumber.length; i++) {
                bw.write("\t"+kinshipFileMap.get("k"+String.valueOf(i+1)));
            }
            bw.newLine();
            for (int i = 0; i < traitName.length; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(traitName[i]).append("\t").append(h2[i]);
                for (int j = 0; j < snpNumber.length; j++) {
                    sb.append("\t").append(h2Sub[i][j]);
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            bw = IoUtils.getTextWriter(enrichFileS);
            bw.write("Trait");
            for (int i = 0; i < snpNumber.length; i++) {
                bw.write("\t"+kinshipFileMap.get("k"+String.valueOf(i+1)));
            }
            bw.newLine();
            for (int i = 0; i < traitName.length; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(traitName[i]);
                for (int j = 0; j < snpNumber.length; j++) {
                    sb.append("\t").append(enrichment[i][j]);
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
 
    
    public void mkVCAPScript () {
        String runFileS = "/local/workdir/fl262/vcap/TASSEL5/run_pipeline.pl";
        String hmpFileS = "/local/workdir/fl262/vcap/genotypes/namrils_projected_hmp31_MAF02mnCnt2500.hmp.txt.gz.lix";
        String positionListDirS = "/local/workdir/fl262/vcap/positionList/";
        String kinshipDirS = "/local/workdir/fl262/vcap/kinship";
        String mkKinshipFileS = "/local/workdir/fl262/mkKinship.pl";
        String kinshipPathFileS = "/local/workdir/fl262/kinshipList.txt";
        String resultDirS = "/local/workdir/fl262/vcap/result/";
        String pheonotypeDirS = "/local/workdir/fl262/vcap/phenotypes/NAM/";
        String ldakPath = "/programs/LDAK/ldak.4.9.fast";
        String estimateH2FileS = "/local/workdir/fl262/h2.pl";
        new File(kinshipDirS).mkdirs();
        new File(resultDirS).mkdirs();
        File[] fs = new File(positionListDirS).listFiles();
        double[] inter = new double[fs.length];
        for (int i = 0; i < fs.length; i++) {
            inter[i] = Double.valueOf(fs[i].getName().split("_")[1]);
        }
        for (int i = 0; i < inter.length-1; i++) {
            for (int j = i+1; j < inter.length; j++) {
                if (inter[i]>inter[j]) {
                    double tDouble = inter[i];
                    inter[i] = inter[j];
                    inter[j] = tDouble;
                    File f = fs[i];
                    fs[i] = fs[j];
                    fs[j] = f;
                }
            }
        }
        
        String[] command = new String[fs.length];
        String[] kinshipPath = new String[fs.length];
        for (int i = 0; i < fs.length; i++) {
            StringBuilder sb = new StringBuilder();
            String outfileS = new File(kinshipDirS, "k"+String.valueOf(i+1)).getAbsolutePath();
            sb.append("system (\"").append("perl ").append(runFileS).append(" -Xmx120g -importGuess ").append(hmpFileS).append(" -FilterSiteBuilderPlugin -positionList ").append(fs[i].getAbsolutePath());
            sb.append(" -endPlugin -KinshipPlugin -method Normalized_IBS -endPlugin -export ").append(outfileS);
            sb.append(" -exportType SqrMatrixBin \");");
            command[i] = sb.toString();
            kinshipPath[i] = outfileS;
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(mkKinshipFileS);
            for (int i = 0; i <  fs.length; i++) {
                bw.write("#"+fs[i].getName()+"\t"+"k"+String.valueOf(i+1));
                bw.newLine();
            }
            for (int i = 0; i < command.length; i++) {
                bw.write(command[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
            bw = IoUtils.getTextWriter(kinshipPathFileS);
            for (int i = 0; i < fs.length; i++) {
                bw.write(kinshipPath[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        fs = new File(pheonotypeDirS).listFiles();
        fs = IoUtils.listFilesEndsWith(fs, ".txt");
        Arrays.sort(fs);
        command = new String[fs.length];
        for (int i = 0; i < fs.length; i++) {
            StringBuilder sb = new StringBuilder();
            String outPrefixS = fs[i].getName().replaceFirst(".txt", "");
            outPrefixS = new File (resultDirS, outPrefixS).getAbsolutePath();
            sb.append("system (\"");
            sb.append(ldakPath).append(" --reml ").append(outPrefixS).append(" --mgrm ").append(kinshipPathFileS);
            sb.append(" --pheno ").append(fs[i].getAbsolutePath()).append(" --kinship-details NO").append("\");");
            command[i] = sb.toString();
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(estimateH2FileS);
            for (int i = 0; i < command.length; i++) {
                bw.write(command[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkPositionListByScore () {
        double[] scoreThresh = {0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5};
        //String infileS  = "/local/workdir/fl262/hapmapScore/hmp31WithScore.txt";
        //String outfileDirS = "/local/workdir/fl262/positionList";
        String infileS  = "M:\\production\\maf\\annotations\\siteUniqueness\\vcap\\hmp31WithScore.txt";
        String outfileDirS = "M:\\production\\maf\\annotations\\siteUniqueness\\vcap\\positionList\\";
        
        String[] outfiles = new String[scoreThresh.length];
        TIntArrayList[] chrList = new TIntArrayList[outfiles.length];
        TIntArrayList[] posList = new TIntArrayList[outfiles.length];
        for (int i = 0; i < outfiles.length; i++) {
            String outS = null;
            if (i == outfiles.length-1) {
                outS = "pos_"+String.valueOf(scoreThresh[i])+"_Plus"+".json.gz";
            }
            else {
                outS = "pos_"+String.valueOf(scoreThresh[i])+"_"+String.valueOf(scoreThresh[i+1])+".json.gz";
                
            }
            outS = new File(outfileDirS, outS).getAbsolutePath();
            outfiles[i] = outS;
            chrList[i] = new TIntArrayList();
            posList[i] = new TIntArrayList();
        }
        
        Table t = new Table(infileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.content[i][2].startsWith("N")) continue;
            double value = t.getDoubleValue(i, 2);
            int index = Arrays.binarySearch(scoreThresh, value);
            if (index < 0) index = -index-2;
            if (index >= outfiles.length-1) index = outfiles.length-1;
            chrList[index].add(t.getIntValue(i, 0));
            posList[index].add(t.getIntValue(i, 1));
        }
        for (int i = 0; i < outfiles.length; i++) {
            VCAP.toPositionListInJson(chrList[i].toArray(), posList[i].toArray(), outfiles[i], IOFileFormat.TextGzip);
        }
    }
    
    public void mkHapMapWithUniqueScore2 () {
        String infileDirS = "M:\\production\\maf\\annotations\\siteUniqueness\\data\\combined_6Genomes";
        String hmpRangeFileS = "M:\\Database\\maize\\vcap\\hmp31PosList\\namrils_projected_hmp31_MAF02mnCnt2500.hmp.pos.range.txt";
        String outfileS = "M:\\production\\maf\\annotations\\siteUniqueness\\vcap\\hmp31WithScore.txt";
        Ranges rs = new Ranges(hmpRangeFileS, IOFileFormat.Text);
        File[] fs = new File(infileDirS).listFiles();
        
        List<File> fList = Arrays.asList(fs);
        ArrayList<String>[] outList = new ArrayList[fs.length];
        for (int i = 0; i < fs.length; i++) {
            outList[i] = new ArrayList();
        }
        fList.parallelStream().forEach(f -> {
            try {
                BufferedReader br = IoUtils.getTextGzipReader(f.getAbsolutePath());
                String temp = br.readLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%10000000 == 0) System.out.println(String.valueOf(cnt)+" lines from " + f.getName());
                    List<String> l = FStringUtils.fastSplit(temp, "\t");
                    int chr = Integer.valueOf(l.get(0));
                    int pos = Integer.valueOf(l.get(1));
                    if (!rs.isInRanges(chr, pos)) continue;
                    outList[chr-1].add(temp);
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Chromosome\tPos\tSiteUniqueScore\tSiteUniqueScoreSD");
            bw.newLine();
            for (int i = 0; i < outList.length; i++) {
                for (int j = 0; j < outList[i].size(); j++) {
                    bw.write(outList[i].get(j));
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
    
    public void mkHapMapWithUniqueScore () {
        String infileDirS = "/local/workdir/fl262/combined_6Genomes/";
        String hmpRangeFileS = "/local/workdir/fl262/hapmapScore/namrils_projected_hmp31_MAF02mnCnt2500.hmp.pos.range.txt";
        String outfileS = "/local/workdir/fl262/hapmapScore/hmp31WithScore.txt";
        Ranges rs = new Ranges(hmpRangeFileS, IOFileFormat.Text);
        File[] fs = new File(infileDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        ArrayList<String>[] outList = new ArrayList[fs.length];
        for (int i = 0; i < fs.length; i++) {
            outList[i] = new ArrayList();
        }
        fList.parallelStream().forEach(f -> {
            try {
                BufferedReader br = IoUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%10000000 == 0) System.out.println(String.valueOf(cnt)+" lines from " + f.getName());
                    List<String> l = FStringUtils.fastSplit(temp, "\t");
                    int chr = Integer.valueOf(l.get(0));
                    int pos = Integer.valueOf(l.get(1));
                    if (!rs.isInRanges(chr, pos)) continue;
                    outList[chr-1].add(temp);
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Chromosome\tPos\tSiteUniqueScore\tSiteUniqueScoreSD");
            bw.newLine();
            for (int i = 0; i < outList.length; i++) {
                for (int j = 0; j < outList[i].size(); j++) {
                    bw.write(outList[i].get(j));
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
    
    /**
     * @deprecated
     */
    private void mkPositionList () {
        String rangeFileS = "M:\\Database\\maize\\vcap\\hmp31PosList\\namrils_projected_hmp31_MAF02mnCnt2500.hmp.pos.range.txt";
        String infileDirS = "M:\\production\\maf\\annotations\\siteUniqueness\\singleCopy\\allChrs\\set01\\";
        String outfileS = "M:\\production\\maf\\annotations\\siteUniqueness\\vcap\\positionLists\\singleCopy.json.txt";
        File[] fs = new File(infileDirS).listFiles();
        TIntArrayList chrList = new TIntArrayList();
        TIntArrayList posList = new TIntArrayList();
        for (int i = 0; i < fs.length; i++) {
            System.out.println(fs[i]);
            Ranges rs = new Ranges(fs[i].getAbsolutePath(), IOFileFormat.Text);
            rs.sortByStartPosition();
            for (int j = 0; j < rs.getRangeNumber(); j++) {
                int start = rs.getRangeStart(j);
                int end = rs.getRangeEnd(j);
                for (int k = start; k < end; k++) {
                    chrList.add(rs.getRangeChromosome(j));
                    posList.add(k);
                }
            }
        }
        Ranges rs = new Ranges(rangeFileS, IOFileFormat.Text);
        VCAP.toPositionListInJson(rs, chrList.toArray(), posList.toArray(), outfileS);
    }
    
    /**
     * @deprecated
     */
    private void getSingleCopyAllChr () {
        String infileDirS = "/local/workdir/fl262/combined_6Genomes/";
        String outfileDirS = "/local/workdir/fl262/singleCopy/";
        double lowThresh = 0.9;
        double highThresh = 1.1;
        double sdThresh = 0.2;
        File[] fs = new File(infileDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            int currentChr = Integer.valueOf(f.getName().split("_")[0].replaceFirst("chr", ""));
            String outfileS = f.getName().split("_")[0]+".singleCopy.range.txt";
            outfileS = new File(outfileDirS, outfileS).getAbsolutePath();
            ArrayList<Range> rList = new ArrayList();
            try {
                BufferedReader br = IoUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                int start = Integer.MIN_VALUE;
                int end = Integer.MAX_VALUE;
                int cnt = 0;
                
                while ((temp = br.readLine()) != null) {
                    List<String> l = FStringUtils.fastSplit(temp, "\t");
                    if (l.get(2).startsWith("N")) continue;
                    double value = Double.valueOf(l.get(2));
                    double sd = Double.valueOf(l.get(3));
                    if (value>lowThresh && value < highThresh && sd < sdThresh) {
                        if (start == Integer.MIN_VALUE) {
                            start = Integer.valueOf(l.get(1));
                            end = start+1;
                        }
                        else {
                            
                        }
                    }
                    else {
                        if (start == Integer.MIN_VALUE) {
                            
                        }
                        else {
                            end = Integer.valueOf(l.get(1));
                            Range r = new Range (currentChr, start, end);
                            rList.add(r);
                            start = Integer.MIN_VALUE;
                            end = Integer.MIN_VALUE;
                        }
                    }
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println(cnt);
                }
                if (end != Integer.MIN_VALUE) {
                    Range r = new Range (currentChr, start, cnt);
                    rList.add(r);
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            StringBuilder sb = new StringBuilder();
            sb.append("Chr10_singleCopy_LowThresh:").append(lowThresh).append("_HighThresh:").append(highThresh).append("_SDThresh:").append(sdThresh);
            Ranges rs = new Ranges(rList, sb.toString());
            rs.writeFile(outfileS, IOFileFormat.Text);
        });
    }
    
    public void analysisPipe () {
        //this.compareScore();
        //this.compareScorePlot();
        //this.getScoreDistribution();
        //this.getLogScoreDistribution();
        //this.getScoreAndSDScatter();
        //this.getScoreAndSDScatterPlot();
        //this.getScoreDistributionOnChr();
        //this.getSingleCopyDistributionOnChr();
        //this.getSingleCopyChr10();
        //this.getRangeLengthDistributionChr10();
        this.geneFootprintChr10();
        //this.geneFootprintAllChr();
        //this.mnaseFootprintChr10();
        //this.mnaseFootprintAllChr();
        //this.functionalOverlap();
        this.getAdh1();
    }
    
    
    private void getAdh1 () {
        String adh1Name = "GRMZM2G442658"; 
        int startPos = 274047276;
        int endPos = 274056225;
        String scoreFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\data\\combined_6Genomes\\chr001_uniqueness.txt.gz";
        String subScoreFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\adh1\\adh1_uniqueness.txt";
        String annotationFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\adh1\\adh1_uniqueness_annotation.txt";
        String geneFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
//        try {
//            BufferedReader br = IoUtils.getTextGzipReader(scoreFileS);
//            BufferedWriter bw = IoUtils.getTextWriter(subScoreFileS);
//            String temp = br.readLine();
//            bw.write(temp);
//            bw.newLine();
//            int cnt = 0;
//            while ((temp = br.readLine()) != null) {
//                cnt++;
//                if (cnt%1000000 == 0) System.out.println(cnt);
//                //List<String> l = FStringUtils.fastSplit(temp);
//                //int pos = Integer.valueOf(l.get(1));
//                if (cnt < startPos) continue;
//                if (cnt > endPos) break;
//                bw.write(temp);
//                bw.newLine();
//            }
//            bw.flush();
//            bw.close();
//        }
//        catch (Exception e) {
//            e.printStackTrace();
//        }
        GeneFeature gf = new GeneFeature(geneFileS);
        ArrayList<Range> cdsList = null;
        ArrayList<Range> utr5List = null;
        ArrayList<Range> utr3List = null;
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            if (gf.getGeneName(i).equals(adh1Name)) {
                cdsList = gf.getCDSList(i, 0);
                utr5List = gf.get5UTRList(i, 0);
                utr3List = gf.get3UTRList(i, 0);
            }
        }
        Ranges cdsRange = new Ranges(cdsList, "CDS");
        Ranges utr5Range = new Ranges(utr5List, "utr5");
        Ranges utr3Range = new Ranges(utr3List, "utr3");
        try {
            BufferedReader br = IoUtils.getTextReader(subScoreFileS);
            BufferedWriter bw = IoUtils.getTextWriter(annotationFileS);
            bw.write(br.readLine()+"\tAnnotation");
            bw.newLine();
            String temp = null;
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                int pos = Integer.valueOf(tem[1]);
                String annotation = "Intron";
                if (cdsRange.isInRanges(1, pos)) annotation ="CDS";
                if (utr5Range.isInRanges(1, pos)) annotation ="5'UTR";
                if (utr3Range.isInRanges(1, pos)) annotation ="3'UTR";
                if (pos >= 274053225) annotation ="Upstream";
                if (pos < 274050254) annotation ="Downstream";
                bw.write(temp+"\t"+annotation);
                bw.newLine();
            }
            bw.flush();
            bw.close();
            
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void functionalOverlap () {
        String mnaseFileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\MNase_shoot.ra.txt";
        String singleFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\singleCopy\\chr010_singleCopy.range.txt";
        String geneFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
        String overlapFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\functionOverlap\\overlap.txt";
        int chr = 10;
        GeneFeature gf = new GeneFeature(geneFileS);
        ArrayList<Range> rList =  new ArrayList();
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            if (gf.getGeneChromosome(i) != chr) continue;
            List<Range> cdsList = gf.getCDSList(i, 0);
            for (int j = 0; j < cdsList.size(); j++) {
                rList.add(cdsList.get(j));
            }
        }
        Ranges cds = new Ranges(rList, "func");
        Ranges srs = new Ranges(singleFileS, IOFileFormat.Text);
        Ranges mrs = new Ranges(mnaseFileS, IOFileFormat.Text);
        rList =  new ArrayList();
        for (int i = 0; i < srs.getRangeNumber(); i++) {
            Range r = srs.getRange(i);
            if (r.getRangeChromosome() != chr) continue;
            rList.add(r);
        }
        Ranges single = new Ranges(rList, "func");
        rList =  new ArrayList();
        for (int i = 0; i < mrs.getRangeNumber(); i++) {
            Range r = mrs.getRange(i);
            if (r.getRangeChromosome() != chr) continue;
            rList.add(r);
        }
        Ranges mnase = new Ranges(rList, "func");
        Ranges mnaseAndSingle = mnase.merge(single);
        mnaseAndSingle.collapse();
        Ranges mnaseAndCDS = mnase.merge(cds);
        mnaseAndCDS.collapse();
        Ranges cdsSingle = cds.merge(single);
        cdsSingle.collapse();
        Ranges singleMnaseCDS = single.merge(mnase).merge(cds);
        singleMnaseCDS.collapse();
        System.out.println("CDS length:\t"+String.valueOf(cds.getTotalRangeSize()));
        System.out.println("MNase length:\t"+String.valueOf(mnase.getTotalRangeSize()));
        System.out.println("Single copy length:\t"+String.valueOf(single.getTotalRangeSize()));
        System.out.println("MnaseAndSingle:\t"+String.valueOf(mnaseAndSingle.getTotalRangeSize()));
        System.out.println("mnaseAndCDS:\t"+String.valueOf(mnaseAndCDS.getTotalRangeSize()));
        System.out.println("cdsSingle:\t"+String.valueOf(cdsSingle.getTotalRangeSize()));
        System.out.println("singleMnaseCDS:\t"+String.valueOf(singleMnaseCDS.getTotalRangeSize()));
        int cnt = 0;
        for (int i = 0; i < single.getRangeNumber(); i++) {
            int start = single.getRangeStart(i);
            int end = single.getRangeEnd(i);
            int ch = single.getRangeChromosome(i);
            for (int j = start; j < end; j++) {
                if (!mnase.isInRanges(ch, j)) continue;
                if (!cds.isInRanges(ch, j)) continue;
                cnt++;
            }
        }
        System.out.println(cnt);
    }
    
    public void mnaseFootprintAllChr () {
        int streamSize = 2000;
        int binLength = 100;
        String rangeFileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\MNase_shoot.ra.txt";
        String infileDirS = "M:\\production\\maf\\annotations\\siteUniqueness\\data\\combined_6Genomes";
        String mnaseFootprintFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\footprint\\mnaseFootprintAllChr.txt";
        Ranges rs = new Ranges(rangeFileS, IOFileFormat.Text);
        int mnaseNumber = 0;
        for (int i = 0; i < rs.getRangeNumber(); i++) {
            mnaseNumber++;
        }
        float[][] upstream =new float[mnaseNumber][];
        float[][] mnase =new float[mnaseNumber][];
        float[][] downstream =new float[mnaseNumber][];
        System.out.println("MnaseNumber:\t " + String.valueOf(mnaseNumber));
        File[] fs = new File(infileDirS).listFiles();
        ArrayList<Integer> chrList = new ArrayList();
        for (int i = 0; i < fs.length; i++) {
            chrList.add(i+1);
        }
        AtomicInteger ac = new AtomicInteger(0);
        chrList.stream().forEach(chr -> {
            String infileS = new File (infileDirS, "chr"+FStringUtils.getNDigitNumber(3, chr)+"_uniqueness.txt.gz").getAbsolutePath();
            int chrIndex = chr  - 1;
            TIntFloatHashMap posScoreMap = new TIntFloatHashMap();
            try {
                BufferedReader br = IoUtils.getTextGzipReader(infileS);
                String temp = br.readLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    List<String> l = FStringUtils.fastSplit(temp, "\t");
                    if (l.get(2).startsWith("N")) continue;
                    posScoreMap.put(Integer.valueOf(l.get(1)), Float.valueOf(l.get(2)));
                    cnt++;
                    if (cnt%10000000==0) System.out.println(cnt+"\tChr"+String.valueOf(chr));
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            int geneCnt = ac.intValue();
            for (int i = 0; i < rs.getRangeNumber(); i++) {
                int c = rs.getRangeChromosome(i);
                if (c != chr) continue;
                byte strand = 1;
                int start = rs.getRangeStart(i)-streamSize;
                int end = rs.getRangeStart(i);
                upstream[geneCnt] = this.addScoreValue(posScoreMap, start, end, strand);
                start = rs.getRangeStart(i);
                end = rs.getRangeEnd(i);
                mnase[geneCnt] = this.addScoreValue(posScoreMap, start, end, strand);
                start = rs.getRangeEnd(i);
                end = rs.getRangeEnd(i)+streamSize;
                downstream[geneCnt] = this.addScoreValue(posScoreMap, start, end, strand);
                geneCnt++;
                ac.incrementAndGet();
            }
        });
        double[][][] meanSD = new double[3][binLength][2];
        meanSD[0] = this.getMeanSD(upstream, binLength);
        meanSD[1] = this.getMeanSD(mnase, binLength);
        meanSD[2] = this.getMeanSD(downstream, binLength);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(mnaseFootprintFileS);
            bw.write("Class\tPosition\tScoreMean\tScoreSD");
            bw.newLine();
            for (int i = 0; i < meanSD[0].length; i++) {
                bw.write("Upstream\t"+String.valueOf(i)+"\t"+String.valueOf(meanSD[0][i][0])+"\t"+String.valueOf(meanSD[0][i][1]));
                bw.newLine();
            }
            for (int i = 0; i < meanSD[1].length; i++) {
                bw.write("Mnase\t"+String.valueOf(i)+"\t"+String.valueOf(meanSD[1][i][0])+"\t"+String.valueOf(meanSD[1][i][1]));
                bw.newLine();
            }
            for (int i = 0; i < meanSD[2].length; i++) {
                bw.write("Downstream\t"+String.valueOf(i)+"\t"+String.valueOf(meanSD[2][i][0])+"\t"+String.valueOf(meanSD[2][i][1]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mnaseFootprintChr10 () {
        int streamSize = 2000;
        int binLength = 100;
        String rangeFileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\MNase_shoot.ra.txt";
        String infileS = "M:\\production\\maf\\annotations\\siteUniqueness\\chr010_uniqueness_6Genomes.txt";
        String mnaseFootprintFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\footprint\\mnaseFootprint.txt";
        Ranges rs = new Ranges(rangeFileS, IOFileFormat.Text);
        int mnaseNumber = 0;
        for (int i = 0; i < rs.getRangeNumber(); i++) {
            if (rs.getRangeChromosome(i)!= 10)continue;
            mnaseNumber++;
        }
        System.out.println("MnaseNumber:\t " + String.valueOf(mnaseNumber));
        float[][] upstream =new float[mnaseNumber][];
        float[][] mnase =new float[mnaseNumber][];
        float[][] downstream =new float[mnaseNumber][];
        TIntFloatHashMap posScoreMap = new TIntFloatHashMap();
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            String temp = br.readLine();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                List<String> l = FStringUtils.fastSplit(temp, "\t");
                if (l.get(2).startsWith("N")) continue;
                posScoreMap.put(Integer.valueOf(l.get(1)), Float.valueOf(l.get(2)));
                cnt++;
                if (cnt%10000000==0) System.out.println(cnt);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        int cnt = 0;
        for (int i = 0; i < rs.getRangeNumber(); i++) {
            int chr = rs.getRangeChromosome(i);
            if (chr != 10) continue;
            byte strand = 1;
            int start = rs.getRangeStart(i)-streamSize;
            int end = rs.getRangeStart(i);
            upstream[cnt] = this.addScoreValue(posScoreMap, start, end, strand);
            start = rs.getRangeStart(i);
            end = rs.getRangeEnd(i);
            mnase[cnt] = this.addScoreValue(posScoreMap, start, end, strand);
            start = rs.getRangeEnd(i);
            end = rs.getRangeEnd(i)+streamSize;
            downstream[cnt] = this.addScoreValue(posScoreMap, start, end, strand);
            cnt++;
        }
        double[][][] meanSD = new double[3][binLength][2];
        meanSD[0] = this.getMeanSD(upstream, binLength);
        meanSD[1] = this.getMeanSD(mnase, binLength);
        meanSD[2] = this.getMeanSD(downstream, binLength);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(mnaseFootprintFileS);
            bw.write("Class\tPosition\tScoreMean\tScoreSD");
            bw.newLine();
            for (int i = 0; i < meanSD[0].length; i++) {
                bw.write("Upstream\t"+String.valueOf(i)+"\t"+String.valueOf(meanSD[0][i][0])+"\t"+String.valueOf(meanSD[0][i][1]));
                bw.newLine();
            }
            for (int i = 0; i < meanSD[1].length; i++) {
                bw.write("Mnase\t"+String.valueOf(i)+"\t"+String.valueOf(meanSD[1][i][0])+"\t"+String.valueOf(meanSD[1][i][1]));
                bw.newLine();
            }
            for (int i = 0; i < meanSD[2].length; i++) {
                bw.write("Downstream\t"+String.valueOf(i)+"\t"+String.valueOf(meanSD[2][i][0])+"\t"+String.valueOf(meanSD[2][i][1]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void geneFootprintAllChr () {
        int streamSize = 2000;
        int binLength = 100;
        String gffFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
        String infileDirS = "M:\\production\\maf\\annotations\\siteUniqueness\\data\\combined_6Genomes\\";
        String GeneFootprintFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\footprint\\GeneFootprintAllChr.txt";
        GeneFeature gf = new GeneFeature(gffFileS);
        int geneNumber = 0;
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            if (!gf.isThere5UTR(i, 0)) continue;
            if (!gf.isThere3UTR(i, 0)) continue;
            geneNumber++;
        }
        System.out.println(String.valueOf(geneNumber)+" genes will be analyzed");
        float[][] upstream =new float[geneNumber][];
        float[][] utr5 =new float[geneNumber][];
        float[][] cds =new float[geneNumber][];
        float[][] intron =new float[geneNumber][];
        float[][] utr3 =new float[geneNumber][];
        float[][] downstream =new float[geneNumber][];
        File[] fs = new File(infileDirS).listFiles();
        ArrayList<Integer> chrList = new ArrayList();
        for (int i = 0; i < fs.length; i++) {
            chrList.add(i+1);
        }
        AtomicInteger ac = new AtomicInteger(0);
        chrList.stream().forEach(chr -> {
            String infileS = new File (infileDirS, "chr"+FStringUtils.getNDigitNumber(3, chr)+"_uniqueness.txt.gz").getAbsolutePath();
            int chrIndex = chr  - 1;
            TIntFloatHashMap posScoreMap = new TIntFloatHashMap();
            try {
                BufferedReader br = IoUtils.getTextGzipReader(infileS);
                String temp = br.readLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    List<String> l = FStringUtils.fastSplit(temp, "\t");
                    if (l.get(2).startsWith("N")) continue;
                    posScoreMap.put(Integer.valueOf(l.get(1)), Float.valueOf(l.get(2)));
                    cnt++;
                    if (cnt%10000000==0) System.out.println(cnt+"\tchr"+String.valueOf(chr));
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            int geneCnt = ac.intValue();
            for (int i = 0; i < gf.getGeneNumber(); i++) {
                if (gf.getGeneChromosome(i) != chr) continue;
                if (!gf.isThere5UTR(i, 0)) continue;
                if (!gf.isThere3UTR(i, 0)) continue;
                byte strand = gf.getTranscriptStrand(i, 0);
                int start = 0;
                int end = 0;
                if (strand == 1) {
                    end = gf.getAllGene5UTRRange().getRangeEnd(0);
                    start = end - streamSize;
                }
                else {
                    ArrayList<Range> rList = gf.get5UTRList(i, 0);
                    start = rList.get(rList.size()-1).getRangeEnd();
                    end = start + streamSize;
                }
                upstream[geneCnt] = this.addScoreValue(posScoreMap, start, end, strand);
                ArrayList<Range> rList = gf.get5UTRList(i, 0);
                utr5[geneCnt] = this.addScoreValue(posScoreMap, rList, strand);
                rList = gf.getCDSList(i, 0);
                cds[geneCnt] = this.addScoreValue(posScoreMap, rList, strand);
                rList = gf.getIntronList(i, 0);
                intron[geneCnt] = this.addScoreValue(posScoreMap, rList, strand);
                rList = gf.get3UTRList(i, 0);
                utr3[geneCnt] = this.addScoreValue(posScoreMap, rList, strand);
                if (strand == 1) {
                    rList = gf.get3UTRList(i, 0);
                    start = rList.get(rList.size()-1).getRangeEnd();
                    end = start + streamSize;
                }
                else {
                    rList = gf.get3UTRList(i, 0);
                    end = rList.get(0).getRangeStart();
                    start = end - streamSize;
                }
                downstream[geneCnt] = this.addScoreValue(posScoreMap, start, end, strand);
                geneCnt++;
                ac.incrementAndGet();
                if (geneCnt%100 == 0) System.out.println(String.valueOf(geneCnt)+" genes has been analyzed");
            }
        });
        
        double[][][] meanSD = new double[6][binLength][2];
        meanSD[0] = this.getMeanSD(upstream, binLength);
        meanSD[1] = this.getMeanSD(utr5, binLength);
        meanSD[2] = this.getMeanSD(cds, binLength);
        meanSD[3] = this.getMeanSD(intron, binLength);
        meanSD[4] = this.getMeanSD(utr3, binLength);
        meanSD[5] = this.getMeanSD(downstream, binLength);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(GeneFootprintFileS);
            bw.write("Class\tPosition\tScoreMean\tScoreSD");
            bw.newLine();
            for (int i = 0; i < meanSD[0].length; i++) {
                bw.write("Upstream\t"+String.valueOf(i)+"\t"+String.valueOf(meanSD[0][i][0])+"\t"+String.valueOf(meanSD[0][i][1]));
                bw.newLine();
            }
            for (int i = 0; i < meanSD[1].length; i++) {
                bw.write("UTR5\t"+String.valueOf(i)+"\t"+String.valueOf(meanSD[1][i][0])+"\t"+String.valueOf(meanSD[1][i][1]));
                bw.newLine();
            }
            for (int i = 0; i < meanSD[2].length; i++) {
                bw.write("CDS\t"+String.valueOf(i)+"\t"+String.valueOf(meanSD[2][i][0])+"\t"+String.valueOf(meanSD[2][i][1]));
                bw.newLine();
            }
            for (int i = 0; i < meanSD[3].length; i++) {
                bw.write("Intron\t"+String.valueOf(i)+"\t"+String.valueOf(meanSD[3][i][0])+"\t"+String.valueOf(meanSD[3][i][1]));
                bw.newLine();
            }
            for (int i = 0; i < meanSD[4].length; i++) {
                bw.write("UTR3\t"+String.valueOf(i)+"\t"+String.valueOf(meanSD[4][i][0])+"\t"+String.valueOf(meanSD[4][i][1]));
                bw.newLine();
            }
            for (int i = 0; i < meanSD[5].length; i++) {
                bw.write("Downstream\t"+String.valueOf(i)+"\t"+String.valueOf(meanSD[5][i][0])+"\t"+String.valueOf(meanSD[5][i][1]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void geneFootprintChr10 () {
        int streamSize = 2000;
        int binLength = 100;
        String gffFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
//        String infileS = "M:\\production\\maf\\annotations\\siteUniqueness\\chr010_uniqueness_6Genomes.txt";
//        String GeneFootprintFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\footprint\\GeneFootprint.txt";
        String infileS = "M:\\production\\elementMap\\combinedScore\\chr010_combinedScore_16.txt";
        String GeneFootprintFileS = "M:\\production\\elementMap\\combinedScore\\GeneFootprint_chr10.txt";
        GeneFeature gf = new GeneFeature(gffFileS);
        int geneNumber = 0;
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            if (gf.getGeneChromosome(i)!= 10)continue;
            if (!gf.isThere5UTR(i, 0)) continue;
            if (!gf.isThere3UTR(i, 0)) continue;
            geneNumber++;
        }
        System.out.println(String.valueOf(geneNumber)+" genes will be analyzed");
        float[][] upstream =new float[geneNumber][];
        float[][] utr5 =new float[geneNumber][];
        float[][] cds =new float[geneNumber][];
        float[][] intron =new float[geneNumber][];
        float[][] utr3 =new float[geneNumber][];
        float[][] downstream =new float[geneNumber][];
        TIntFloatHashMap posScoreMap = new TIntFloatHashMap();
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            String temp = br.readLine();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                List<String> l = FStringUtils.fastSplit(temp, "\t");
                if (l.get(2).startsWith("N")) continue;
                posScoreMap.put(Integer.valueOf(l.get(1)), Float.valueOf(l.get(2)));
                cnt++;
                if (cnt%10000000==0) System.out.println(cnt);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        int cnt = 0;
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            if (gf.getGeneChromosome(i)!= 10)continue;
            if (!gf.isThere5UTR(i, 0)) continue;
            if (!gf.isThere3UTR(i, 0)) continue;
            byte strand = gf.getTranscriptStrand(i, 0);
            int start = 0;
            int end = 0;
            if (strand == 1) {
                end = gf.getAllGene5UTRRange().getRangeEnd(0);
                start = end - streamSize;
            }
            else {
                ArrayList<Range> rList = gf.get5UTRList(i, 0);
                start = rList.get(rList.size()-1).getRangeEnd();
                end = start + streamSize;
            }
            upstream[cnt] = this.addScoreValue(posScoreMap, start, end, strand);
            ArrayList<Range> rList = gf.get5UTRList(i, 0);
            utr5[cnt] = this.addScoreValue(posScoreMap, rList, strand);
            rList = gf.getCDSList(i, 0);
            cds[cnt] = this.addScoreValue(posScoreMap, rList, strand);
            rList = gf.getIntronList(i, 0);
            intron[cnt] = this.addScoreValue(posScoreMap, rList, strand);
            rList = gf.get3UTRList(i, 0);
            utr3[cnt] = this.addScoreValue(posScoreMap, rList, strand);
            if (strand == 1) {
                rList = gf.get3UTRList(i, 0);
                start = rList.get(rList.size()-1).getRangeEnd();
                end = start + streamSize;
            }
            else {
                rList = gf.get3UTRList(i, 0);
                end = rList.get(0).getRangeStart();
                start = end - streamSize;
            }
            downstream[cnt] = this.addScoreValue(posScoreMap, start, end, strand);
            cnt++;
            if (cnt%100 == 0) System.out.println(String.valueOf(cnt)+" genes has been analyzed");
        }
        double[][][] meanSD = new double[6][binLength][2];
        meanSD[0] = this.getMeanSD(upstream, binLength);
        meanSD[1] = this.getMeanSD(utr5, binLength);
        meanSD[2] = this.getMeanSD(cds, binLength);
        meanSD[3] = this.getMeanSD(intron, binLength);
        meanSD[4] = this.getMeanSD(utr3, binLength);
        meanSD[5] = this.getMeanSD(downstream, binLength);
        try {
            BufferedWriter bw = IoUtils.getTextWriter(GeneFootprintFileS);
            bw.write("Class\tPosition\tScoreMean\tScoreSD");
            bw.newLine();
            for (int i = 0; i < meanSD[0].length; i++) {
                bw.write("Upstream\t"+String.valueOf(i)+"\t"+String.valueOf(meanSD[0][i][0])+"\t"+String.valueOf(meanSD[0][i][1]));
                bw.newLine();
            }
            for (int i = 0; i < meanSD[1].length; i++) {
                bw.write("UTR5\t"+String.valueOf(i)+"\t"+String.valueOf(meanSD[1][i][0])+"\t"+String.valueOf(meanSD[1][i][1]));
                bw.newLine();
            }
            for (int i = 0; i < meanSD[2].length; i++) {
                bw.write("CDS\t"+String.valueOf(i)+"\t"+String.valueOf(meanSD[2][i][0])+"\t"+String.valueOf(meanSD[2][i][1]));
                bw.newLine();
            }
            for (int i = 0; i < meanSD[3].length; i++) {
                bw.write("Intron\t"+String.valueOf(i)+"\t"+String.valueOf(meanSD[3][i][0])+"\t"+String.valueOf(meanSD[3][i][1]));
                bw.newLine();
            }
            for (int i = 0; i < meanSD[4].length; i++) {
                bw.write("UTR3\t"+String.valueOf(i)+"\t"+String.valueOf(meanSD[4][i][0])+"\t"+String.valueOf(meanSD[4][i][1]));
                bw.newLine();
            }
            for (int i = 0; i < meanSD[5].length; i++) {
                bw.write("Downstream\t"+String.valueOf(i)+"\t"+String.valueOf(meanSD[5][i][0])+"\t"+String.valueOf(meanSD[5][i][1]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private double[][] getMeanSD (float[][] feature, int binLength) {
        double[] bound = new double[binLength];
        double[][] meanSD = new double[binLength][2];
        TDoubleArrayList[] valueList = new TDoubleArrayList[binLength];
        for (int i = 0; i < binLength; i++) {
            bound[i] = (double)i/binLength;
            valueList[i] = new TDoubleArrayList();
        }
        for (int i = 0; i < feature.length; i++) {
            for (int j = 0; j < feature[i].length; j++) {
                double v = (double)j/feature[i].length;
                int index = Arrays.binarySearch(bound, v);
                if (index < 0) index = -index -2;
                if (index > bound.length-1) index = bound.length-1;
                valueList[index].add(feature[i][j]);
            }
        }
        for (int i = 0; i < valueList.length; i++) {
            double[] value = valueList[i].toArray();
            DescriptiveStatistics d = new DescriptiveStatistics(value);
            meanSD[i][0] = d.getMean();
            meanSD[i][1] = d.getStandardDeviation();
        }
        return meanSD;
    }
    
    private float[] addScoreValue (TIntFloatHashMap posScoreMap, List<Range> rList, byte strand) {
        TFloatArrayList fList = new TFloatArrayList();
        for (int i = 0; i < rList.size(); i++) {
            int start = rList.get(i).getRangeStart();
            int end = rList.get(i).getRangeEnd();
            for (int j = start; j < end; j++) {
                float value = posScoreMap.get(j); 
                if (value == 0) continue;
                fList.add(value);
            }
        }
        if (strand == -1) fList.reverse();
        float[] value = fList.toArray();
        return value;
    }
    
    private float[] addScoreValue (TIntFloatHashMap posScoreMap, int start, int end, byte strand) {
        TFloatArrayList fList = new TFloatArrayList();
        for (int j = start; j < end; j++) {
            float value = posScoreMap.get(j);
            if (value == 0) continue;
            fList.add(value);
        }
        if (strand == -1) fList.reverse();
        float[] value = fList.toArray();
        return value;
    }
    
    private void getRangeLengthDistributionChr10 () {
        String infileS = "M:\\production\\maf\\annotations\\siteUniqueness\\singleCopy\\chr010_singleCopy.range.txt";
        String outfileS = "M:\\production\\maf\\annotations\\siteUniqueness\\singleCopy\\chr010_singleCopy.distribution.txt";
        Ranges rs = new Ranges(infileS, IOFileFormat.Text);
        int[] len = new int[rs.getRangeNumber()];
        int max = Integer.MIN_VALUE;
        for (int i = 0; i < rs.getRangeNumber(); i++) {
            len[i] = rs.getRangeSize(i);
            if (len[i] > max) max = len[i];
        }
        int[] count = new int[max+1];
        for (int i = 0; i < len.length; i++) count[len[i]]++;
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Length\tCount\tBaseNumber");
            bw.newLine();
            for (int i = 0; i <  count.length; i++) {
                bw.write(String.valueOf(i)+"\t"+String.valueOf(count[i])+"\t"+String.valueOf(i*count[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void getSingleCopyChr10 () {
        String infileS = "M:\\production\\maf\\annotations\\siteUniqueness\\chr010_uniqueness_6Genomes.txt";
        String outfileS = "M:\\production\\maf\\annotations\\siteUniqueness\\singleCopy\\chr010_singleCopy.range.txt";
        double lowThresh = 0.8;
        double highThresh = 1.2;
        double sdThresh = Double.MAX_VALUE;
        ArrayList<Range> rList = new ArrayList();
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            String temp = br.readLine();
            int start = Integer.MIN_VALUE;
            int end = Integer.MAX_VALUE;
            int cnt = 0;
            int totalLength = 0;
            while ((temp = br.readLine()) != null) {
                List<String> l = FStringUtils.fastSplit(temp, "\t");
                if (l.get(2).startsWith("N")) {
                    cnt++;
                    continue;
                }
                double value = Double.valueOf(l.get(2));
                double sd = Double.valueOf(l.get(3));
                if (value>lowThresh && value < highThresh && sd < sdThresh) {
                    if (start == Integer.MIN_VALUE) {
                        start = Integer.valueOf(l.get(1));
                        end = start+1;
                    }
                    else {
                        
                    }
                }
                else {
                    if (start == Integer.MIN_VALUE) {
                        
                    }
                    else {
                        end = Integer.valueOf(l.get(1));
                        Range r = new Range (10, start, end);
                        rList.add(r);
                        start = Integer.MIN_VALUE;
                        end = Integer.MIN_VALUE;
                        totalLength +=r.getRangeSize();
                    }
                }
                cnt++;
                if (cnt%1000000 == 0) System.out.println(cnt);
            }
            if (end != Integer.MIN_VALUE) {
                Range r = new Range (10, start, cnt);
                rList.add(r);
                totalLength +=r.getRangeSize();
            }
            System.out.println("Total length:\t"+String.valueOf(totalLength));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        StringBuilder sb = new StringBuilder();
        sb.append("Chr10_singleCopy_LowThresh:").append(lowThresh).append("_HighThresh:").append(highThresh).append("_SDThresh:").append(sdThresh);
        Ranges rs = new Ranges(rList, sb.toString());
        rs.writeFile(outfileS, IOFileFormat.Text);
    }
    
    public void getSingleCopyDistributionOnChr () {
        int chrLength = 149627545;
        double lowThresh = 0.5;
        double highThresh = 1.5;
        double sdThresh = Double.MIN_VALUE;
        String infileS = "M:\\production\\maf\\annotations\\siteUniqueness\\chr010_uniqueness_6Genomes.txt";
        String binFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\scoreDistribution\\chr010_singleCopy_bin.txt";
        int binSize = 100000;
        TDoubleArrayList valueList = new TDoubleArrayList();
        TIntArrayList posList = new TIntArrayList();
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            String temp = br.readLine();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                List<String> list = FStringUtils.fastSplit(temp, "\t");
                if (list.get(2).startsWith("N")) continue;
                double value = Double.valueOf(list.get(2));
                double sd = Double.valueOf(list.get(3));
                if (value > highThresh || value < lowThresh || sd > sdThresh) {
                    continue;
                }
                valueList.add(Double.valueOf(list.get(2)));
                posList.add(Integer.valueOf(list.get(1)));
                cnt++;
                if (cnt%1000000 == 0) System.out.println(cnt);
            }
            double[] value = valueList.toArray();
            int[] pos = posList.toArray();
            Bins bin = new Bins(1, chrLength, binSize, pos, value);
            BufferedWriter bw = IoUtils.getTextWriter(binFileS);
            bw.write("BinStart\tNumberOfSingleCopySite");
            bw.newLine();
            for (int i = 0; i < bin.getBinNum(); i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(bin.getBinStart(i)).append("\t").append(bin.getNumValues(i));
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
    
    public void getScoreDistributionOnChr () {
        int chrLength = 149627545;
        String infileS = "M:\\production\\maf\\annotations\\siteUniqueness\\chr010_uniqueness_6Genomes.txt";
        String binFileS = "M:\\production\\maf\\annotations\\siteUniqueness\\scoreDistribution\\chr010_bin.txt";
        int binSize = 100000;
        TDoubleArrayList valueList = new TDoubleArrayList();
        TIntArrayList posList = new TIntArrayList();
        
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            String temp = br.readLine();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                List<String> list = FStringUtils.fastSplit(temp, "\t");
                if (list.get(2).startsWith("N")) continue;
                valueList.add(Double.valueOf(list.get(2)));
                posList.add(Integer.valueOf(list.get(1)));
                cnt++;
                if (cnt%1000000 == 0) System.out.println(cnt);
            }
            double[] value = valueList.toArray();
            int[] pos = posList.toArray();
            Bins bin = new Bins(1, chrLength, binSize, pos, value);
            BufferedWriter bw = IoUtils.getTextWriter(binFileS);
            bw.write("BinStart\tMean\tSD");
            bw.newLine();
            for (int i = 0; i < bin.getBinNum(); i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(bin.getBinStart(i)).append("\t").append(bin.getBinAverage(i)).append("\t").append(bin.getBinSD(i));
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
    
    public void getScoreAndSDScatterPlot () {
        double lowThresh = 0;
        double highThresh = 3;
        String infileS = "M:\\production\\maf\\annotations\\siteUniqueness\\scoreDistribution\\scoreVsSD.scatter.txt";
        String outfileS = "M:\\production\\maf\\annotations\\siteUniqueness\\scoreDistribution\\scoreVsSD.scatter.pdf";
        TDoubleArrayList scoreList = new TDoubleArrayList();
        TDoubleArrayList sdList = new TDoubleArrayList();
        Table t = new Table (infileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            double score = t.getDoubleValue(i, 0);
            double sd = t.getDoubleValue(i, 1);
            if (score < lowThresh || score > highThresh) continue;
            scoreList.add(score);
            sdList.add(sd);
        }
        double[] x = scoreList.toArray();
        double[] y = sdList.toArray();
        ScatterPlot d = new ScatterPlot(x, y);
        d.setColor(255, 0, 0, 5);
        d.setPlottingCharacter(16);
        d.setXLab("Uniqueness score");
        d.setYLab("SD of uniquess score");
        d.saveGraph(outfileS);
        //d.showGraph();
    }
    
    public void getScoreAndSDScatter () {
        int size = 10000;
        double lowThresh = 0.5;
        double highThresh = 3;
        String infileS = "M:\\production\\maf\\annotations\\siteUniqueness\\chr010_uniqueness_6Genomes.txt";
        String outfileS = "M:\\production\\maf\\annotations\\siteUniqueness\\scoreDistribution\\scoreVsSD.scatter.txt";
        TDoubleArrayList scoreList = new TDoubleArrayList();
        TDoubleArrayList sdList = new TDoubleArrayList();
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            String temp = br.readLine();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                List<String> l = FStringUtils.fastSplit(temp, "\t");
                if (l.get(2).startsWith("N")) continue;
                double score = Double.valueOf(l.get(2));
                if (score < lowThresh || score > highThresh) {
                    cnt++;
                    continue;
                }
                double sd = Double.valueOf(l.get(3));
                scoreList.add(score);
                sdList.add(sd);
                cnt++;
                if (cnt%1000000 == 0) System.out.println(cnt);
            }
            br.close();
            double[] score = scoreList.toArray();
            double[] sd = sdList.toArray();
            int[] index = FArrayUtils.getRandomIntArray(score.length, size);
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("UniqueScore\tSD");
            bw.newLine();
            for (int i = 0; i < index.length; i++) {
                bw.write(String.valueOf(score[index[i]])+"\t"+sd[index[i]]);
                bw.newLine();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void getLogScoreDistribution () {
        int chrLength = 149627545;
        int size = 10000;
        String infileS = "M:\\production\\maf\\annotations\\siteUniqueness\\chr010_uniqueness.txt";
        int[] index = FArrayUtils.getRandomIntArray(chrLength, size);
        Arrays.sort(index);
        TDoubleArrayList vList = new TDoubleArrayList();
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            String temp = br.readLine();
            int cnt = -1;
            while ((temp = br.readLine()) != null) {
                cnt++;
                if (cnt%1000000 == 0) System.out.println(cnt);
                if (Arrays.binarySearch(index, cnt) < 0) continue;
                List<String> l = FStringUtils.fastSplit(temp, "\t");
                if (l.get(2).startsWith("N")) continue;
                double value = Math.log10(Double.valueOf(l.get(2)));
                vList.add(value);
                
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        double[] value = vList.toArray();
        Histogram h = new Histogram(value);
        h.setBreakNumber(100);
        h.showGraph();
    }
    
    public void getScoreDistribution () {
        int size = 10000;
        String infileS = "M:\\production\\maf\\annotations\\siteUniqueness\\chr010_uniqueness_6Genomes.txt";
        String outfileS = "M:\\production\\maf\\annotations\\siteUniqueness\\scoreDistribution\\chr10_scoreDistribution.txt";
        int binNumber = 100;
        double interval = 0.5;
        double[] binStart = new double[binNumber];
        int[] binCounts = new int[binNumber];
        for (int i = 0; i < binNumber; i++) binStart[i] = i*interval;
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            String temp = br.readLine();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                List<String> l = FStringUtils.fastSplit(temp, "\t");
                if (l.get(2).startsWith("N")) continue;
                double value = Double.valueOf(l.get(2));
                int index = Arrays.binarySearch(binStart, value);
                if (index < 0) index = -index -2;
                if (index == binNumber) index = binNumber-1;
                binCounts[index]++;
                cnt++;
                if (cnt%1000000 == 0) System.out.println(cnt);
            }
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("SiteUniqueScore\tCount\tProportion");
            bw.newLine();
            for (int i = 0; i < binNumber; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(binStart[i]).append("\t").append(binCounts[i]).append("\t").append((double)binCounts[i]/cnt);
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
    private void compareScorePlot () {
        String infileS = "M:\\production\\maf\\annotations\\siteUniqueness\\qualityCheck\\genome4VsGenome6.txt";
        String outfileS = "M:\\production\\maf\\annotations\\siteUniqueness\\qualityCheck\\genome4VsGenome6.pdf";
        double lowThresh = 0.5;
        double highThresh = 3;
        Table t = new Table (infileS);
        TDoubleArrayList v1List = new TDoubleArrayList();
        TDoubleArrayList v2List = new TDoubleArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            double v1 = t.getDoubleValue(i, 0);
            double v2 = t.getDoubleValue(i, 1);
            if ((v1>lowThresh && v1 < highThresh) || (v2>lowThresh && v2 < highThresh)) {
                v1List.add(v1);
                v2List.add(v2);
            }
        }
        ScatterPlot s = new ScatterPlot(v1List.toArray(), v2List.toArray());
        s.setXLim(0.5, 5);
        s.setYLim(0.5, 5);
        s.setXLab("UScore from 4 genomes");
        s.setYLab("UScore from 6 genomes");
        s.setPlottingCharacter(16);
        s.setColor(255, 0, 0, 5);
        s.saveGraph(outfileS);
    }
    
    private void compareScore () {
        String infile1S = "M:\\production\\maf\\annotations\\siteUniqueness\\chr010_uniqueness_4Genomes.txt";
        String infile2S = "M:\\production\\maf\\annotations\\siteUniqueness\\chr010_uniqueness_6Genomes.txt";
        String outfileS = "M:\\production\\maf\\annotations\\siteUniqueness\\qualityCheck\\genome4VsGenome4.txt";
        int chrLength = 149627545;
        int size = 10000;
        int[] index = FArrayUtils.getRandomIntArray(chrLength, size);
        Arrays.sort(index);
        TDoubleArrayList v1List = new TDoubleArrayList();
        TDoubleArrayList v2List = new TDoubleArrayList();
        try {
            BufferedReader br1 = IoUtils.getTextReader(infile1S);
            BufferedReader br2 = IoUtils.getTextReader(infile2S);
            String temp1 = br1.readLine();
            String temp2 = br2.readLine();
            int cnt = -1;
            while ((temp1 = br1.readLine()) != null) {
                cnt++;
                if (cnt%1000000 == 0) System.out.println(cnt);
                temp2 = br2.readLine();
                List<String> list1 = FStringUtils.fastSplit(temp1, "\t");
                List<String> list2 = FStringUtils.fastSplit(temp2, "\t");
                if (Arrays.binarySearch(index, cnt) < 0) continue;
                if (list1.get(2).startsWith("N")) continue;
                v1List.add(Double.valueOf(list1.get(2)));
                v2List.add(Double.valueOf(list2.get(2)));
                
            }
            br1.close();
            br2.close();
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("UScore4Genome\tUScore6Genomes");
            bw.newLine();
            double[] v1 = v1List.toArray();
            double[] v2 = v2List.toArray();
            for (int i = 0; i < v1.length; i++) {
                bw.write(String.valueOf(v1[i])+"\t"+String.valueOf(v2[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void combineUniqueness () {
        String infileDirS = "/workdir/fl262/taxa/";
        String outfileDirS = "/workdir/fl262/combined/";
//        String infileDirS = "M:\\production\\maf\\annotations\\siteUniqueness\\taxa\\";
//        String outfileDirS = "M:\\production\\maf\\annotations\\siteUniqueness\\test\\";
        
        new File(outfileDirS).mkdir();
        File[] taxaDir = new File(infileDirS).listFiles();
        File[] tempFs = taxaDir[0].listFiles();
        ArrayList<Integer> chrList = new ArrayList();
        for (int i = 0; i < tempFs.length; i++) {
            chrList.add(Integer.valueOf(tempFs[i].getName().split("_")[0].replaceFirst("chr", "")));
        }
        chrList.parallelStream().forEach(chr -> {
            String fileName = "chr"+FStringUtils.getNDigitNumber(3, chr)+"_uniqueness.txt.gz";
            String outfileS = new File (outfileDirS, fileName).getAbsolutePath();
            BufferedReader[] brs = new BufferedReader[taxaDir.length];
            for (int i = 0; i < taxaDir.length; i++) {
                String infileS = new File(taxaDir[i], fileName).getAbsolutePath();
                brs[i] = IoUtils.getTextGzipReader(infileS);
            }
            try {
                BufferedWriter bw = IoUtils.getTextGzipWriter(outfileS);
                bw.write("Chromosome\tPos\tCopyNumberScore\tCopyNumberScoreSD");
                bw.newLine();
                for (int i = 0; i < brs.length; i++) brs[i].readLine();
                String[] temps = new String[brs.length];
                int cnt = 0;
                while ((temps[0] = brs[0].readLine()) != null) {
                    for (int i = 1; i <  brs.length; i++) temps[i] = brs[i].readLine();
                    List<String> l = FStringUtils.fastSplit(temps[0], "\t");
                    String chrS = l.get(0);
                    String posS = l.get(1);
                    if (l.get(2).startsWith("N")) {
                        bw.write(chrS+"\t"+posS+"\tNA\tNA");
                        bw.newLine();
                        continue;
                    }
                    double[] value = new double[brs.length];
                    value[0] = Double.valueOf(l.get(2));
                    for (int i = 1; i < temps.length; i++) {
                        l = FStringUtils.fastSplit(temps[i], "\t");
                        value[i] = Double.valueOf(l.get(2));
                    }
                    DescriptiveStatistics d = new DescriptiveStatistics(value);
                    StringBuilder sb = new StringBuilder();
                    sb.append(chrS).append("\t").append(posS).append("\t").append((float)d.getMean()).append("\t").append((float)d.getStandardDeviation());
                    bw.write(sb.toString());
                    bw.newLine();
                    cnt++;
                    if (cnt%10000000 == 0) System.out.println("Finished " + String.valueOf(cnt) + " sites on chromosome " + chrS);
                }
                bw.flush();
                bw.close();
            }
            catch(Exception e) {
                e.printStackTrace();
            }
        });
    }
    
    public void buildUniquenessTableBatch () {
        //String referenceGenome = "/workdir/fl262/reference/maize.agpV3.fa";
        //String referenceGenome = "/workdir/fl262/reference/cassavaV6_23Chr.fa";
        String referenceGenome = "/workdir/fl262/reference/Sorghum_v3.0.fa";
        String testGenomeDirS = "/workdir/fl262/genomes/";
        String outfileDirS = "/workdir/fl262/output/";
        new File(outfileDirS).mkdir();
        File[] fs = new File(testGenomeDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            String taxa = fs[i].getName().replaceFirst(".fa", "");
            String outDirS = new File(outfileDirS, "uniqueness_"+taxa).getAbsolutePath();
            this.buildUniquenessTable(referenceGenome, fs[i].getAbsolutePath(), outDirS);
        }
    }
    
    public void buildUniquenessTable (String referenceGenome, String testGenome, String outDirS) {
//        String referenceGenome = "M:\\production\\maf\\annotations\\siteUniqueness\\reference.fa";
//        String testGenome = "M:\\production\\maf\\annotations\\siteUniqueness\\test.fa";
//        String outfileDirS = "M:\\production\\maf\\annotations\\siteUniqueness\\output\\";
        new File(outDirS).mkdir();
        //fragmentSize should >= than kmerLength
        int kmerLength = 32;
        int fragmentSize = 100000;
        this.buildAscIIByteMap();
        ConcurrentHashMap<Long, Integer> kmerCountMap = new ConcurrentHashMap();
        Fasta f = new Fasta (testGenome);
        HashMap<Integer, String> chrSeqMap = new HashMap();
        for (int i = 0; i < f.getSeqNumber(); i++) {
            if (f.getSeqLength(i) < kmerLength) continue;
            chrSeqMap.put(i+1, f.getSeq(i));
        }
        Set<Map.Entry<Integer, String>> chrSeqset = chrSeqMap.entrySet();
        chrSeqset.parallelStream().forEach(e -> {
            int chr = e.getKey();
            String seq = e.getValue();
            byte[] bArray = this.convertByteArray(seq.getBytes());
            long[] seqL = new long[2];
            for (int i = 0; i < bArray.length-kmerLength+1; i++) {
                boolean flag = false;
                for (int j = i; j < i+kmerLength; j++) {
                    if (bArray[j] >3) {
                        i = j;
                        flag = true;
                        break;
                    }
                }
                if (flag) continue;
                byte[] subArray = Arrays.copyOfRange(bArray, i, i + kmerLength);
                seqL[0] = BaseEncoder.getLongSeqFromByteArray(subArray);
                seqL[1] = BaseEncoder.getLongSeqFromByteArray(this.getReverseComplementary(subArray));
                for (int j = 0; j < seqL.length; j++) {
                    if (kmerCountMap.containsKey(seqL[j])) {
                        kmerCountMap.computeIfPresent(seqL[j], (k, v) -> v + 1);
                    }
                    else {
                        kmerCountMap.put(seqL[j], 1);
                    }
                }
                int s = kmerCountMap.size();
                if (s%1000000 == 0) System.out.println("KmerCountMap size: "+String.valueOf(s));
            }
        });
        //this.printKmerMap(kmerCountMap);
        f = new Fasta(referenceGenome);
        System.out.println("Writing site uniqueness by chromosomes...");
        chrSeqMap = new HashMap();
        for (int i = 0; i < f.getSeqNumber(); i++) {
            chrSeqMap.put(Integer.valueOf(f.getName(i)), f.getSeq(i));
        }
        chrSeqset = chrSeqMap.entrySet();
        chrSeqset.parallelStream().forEach(entry -> {
            int chr = entry.getKey();
            String seq = entry.getValue();
            byte[] bArray = this.convertByteArray(seq.getBytes());
            StringBuilder sb  = new StringBuilder();
            sb.append("chr").append(FStringUtils.getNDigitNumber(3, chr)).append("_uniqueness.txt");
            String outfileS = new File(outDirS, sb.toString()).getAbsolutePath();
            int[][] bound = FArrayUtils.getSubsetsIndicesBySubsetSize(seq.length(), fragmentSize);
            try {
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                bw.write("Chromosome\tPos\tKmerUniqueness");
                bw.newLine();
                for (int i = 0; i < bound.length; i++) {
                    int intervalSize = bound[i][1] - bound[i][0];
                    TIntArrayList[] kmerCountList = new TIntArrayList[intervalSize];
                    for (int j = 0; j < kmerCountList.length; j++) kmerCountList[j] = new TIntArrayList();
                    int startIndex = bound[i][0]-kmerLength+1;
                    if (startIndex < 0) startIndex = 0;
                    int endIndex = bound[i][1];
                    if (endIndex-1+kmerLength > seq.length()) endIndex = seq.length()-kmerLength+1;
                    for (int j = startIndex; j < endIndex; j++) {
                        boolean flag = false;
                        for (int k = j; k < j+kmerLength; k++) {
                            if (bArray[k] >3) {
                                j = k;
                                flag = true;
                                break;
                            }
                        }
                        if (flag) continue;
                        
                        long query = BaseEncoder.getLongSeqFromByteArray(Arrays.copyOfRange(bArray, j, j + kmerLength));
                        Integer count = kmerCountMap.get(query);
                        if (count == null) count = 0;
                        int offSet = bound[i][0]-j;
                        
                        if (offSet > 0) {
                            for (int k = j; k < j+kmerLength-offSet; k++) {
                                kmerCountList[offSet+k-bound[i][0]].add(count);
                            }
                        }
                        else {
                            int end = j+kmerLength;
                            if (end > bound[i][1]) end = bound[i][1];
                            for (int k = j; k < end; k++) {
                                kmerCountList[k-bound[i][0]].add(count);
                            }
                        }
                    }
                    for (int j = 0; j < kmerCountList.length; j++) {
                        int[] kmerCount = kmerCountList[j].toArray();
                        sb = new StringBuilder();
                        sb.append(chr).append("\t").append(bound[i][0]+j+1).append("\t");
                        if (kmerCount.length == 0) {
                            sb.append("NA");
                        }
                        else {
                            float aver = 0;
                            for (int k = 0; k < kmerCount.length; k++) {
                                aver+=(float)kmerCount[k]/kmerCount.length;
                            }
                            sb.append(aver);
                        }
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.gc();
        });
        System.out.println("Results output to" + outDirS);
        System.out.println("Finished");
    }
    
    private byte[] getReverseComplementary (byte[] a) {
        byte[] b = new byte[a.length];
        for (int i = 0; i < a.length; i++) {
            byte t;
            if (a[i] == 0) t = 3;
            else if (a[i] == 1) t = 2;
            else if (a[i] == 2) t = 1;
            else if (a[i] == 3) t = 0;
            else t = a[i];
            b[a.length-i-1] = t;
        }
        return b;
    }
    
    private void buildAscIIByteMap () {
        this.ascIIByteMap = new HashMap();
        for (int i = 0; i < 128; i++) {
            ascIIByteMap.put((byte)i, (byte)4);
        }
        ascIIByteMap.put((byte)65, (byte)0);
        ascIIByteMap.put((byte)97, (byte)0);
        ascIIByteMap.put((byte)67, (byte)1);
        ascIIByteMap.put((byte)99, (byte)1);
        ascIIByteMap.put((byte)71, (byte)2);
        ascIIByteMap.put((byte)103, (byte)2);
        ascIIByteMap.put((byte)84, (byte)3);
        ascIIByteMap.put((byte)116, (byte)3);
        ascIIByteMap.put((byte)78, (byte)4);
        ascIIByteMap.put((byte)110, (byte)4);
    }
    
    private byte[] convertByteArray (byte[] bArray) {
        byte[] cArray = new byte[bArray.length];
        for (int i = 0; i < cArray.length; i++) {
            cArray[i] = ascIIByteMap.get(bArray[i]);
        }
        return cArray;
    }
    
    private void printKmerMap (ConcurrentHashMap<Long, Integer> kmerCountMap) {
        for (Map.Entry<Long, Integer> e : kmerCountMap.entrySet()) {
            System.out.println(BaseEncoder.getSequenceFromLong(e.getKey())+"\t"+String.valueOf(e.getValue()));
        }
    }
}
