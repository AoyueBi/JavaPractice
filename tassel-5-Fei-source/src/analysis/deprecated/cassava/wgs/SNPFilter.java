/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.cassava.wgs;

import format.Fasta;
import format.JPileUp;
import format.Table;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import graphcis.r.DensityPlot;
import graphcis.r.ScatterPlot;
import graphcis.r.ScatterPlotMultiClass;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import utils.FArrayUtils;
import utils.FStringUtils;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class SNPFilter {
    
    public SNPFilter () {
        //this.cloneDepthPipe();
        this.filterTestPipe();
        //this.badRefSitePipe();
    }
    
    public void badRefSitePipe () {
        //this.mkChrDepthFiles();
        //this.estimateTaxaDepth();
        //this.estimateRelativeDepthAndVariance();
        //this.getSampleDepthAndSD();
        //this.mkDepthAndSDFigure();
        //this.mkValidSites();
    }
    
    public void mkValidSites () {
        String infileDirS = "M:\\pipelineTest\\cassava\\wgs\\badRefSite\\relativeDepth\\";
        String validSiteDirS = "M:\\pipelineTest\\cassava\\wgs\\badRefSite\\validSite\\";
        File[] fs = new File(infileDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String outfileS = new File(validSiteDirS, f.getName().replaceFirst(".DAndV.txt", ".validSite.txt")).getAbsolutePath();
            try {
                BufferedReader br = IoUtils.getTextReader(f.getAbsolutePath());
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                bw.write("Chr\tPos");
                bw.newLine();
                String temp = br.readLine();
                while ((temp = br.readLine())!=null) {
                    String[] tem = temp.split("\t");
                    float mean = Float.valueOf(tem[2]);
                    float sd = Float.valueOf(tem[3]);
                    if (mean < 0.68) continue;
                    else if (mean < 0.764) {
                        double y = 5*mean-3.4;
                        if (sd > y) continue;
                    }
                    else if (mean < 1.416) {
                        if (sd > 0.42) continue;
                    }
                    else if (mean  < 1.5) {
                        double y = 7.5*mean-5;
                        if (sd > y) continue;
                    }
                    else continue;
                    bw.write(tem[0]+"\t"+tem[1]);
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        
    }
    
    public void mkDepthAndSDFigure () {
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\badRefSite\\depthAndSD.txt";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\badRefSite\\depthAndSD.pdf";
        Table t = new Table (infileS);
        double[] mValue = new double[t.getRowNumber()];
        double[] sValue = new double[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            mValue[i] = t.getDoubleValue(i, 1);
            sValue[i] = t.getDoubleValue(i, 2);
        }
        ScatterPlot sp = new ScatterPlot(mValue, sValue);
        sp.setXLab("Relative Depth");
        sp.setYLab("SD");
        sp.setColor(255, 0, 0, 50);
        sp.setXLim(0, 2);
        sp.setYLim(0, 1);
        sp.showGraph();
    }
    
    public void getSampleDepthAndSD () {
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\badRefSite\\relativeDepth\\chr001.DAndV.txt";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\badRefSite\\depthAndSD.txt";
        TFloatArrayList meanList = new TFloatArrayList();
        TFloatArrayList sdList = new TFloatArrayList();
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            String temp = br.readLine();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                meanList.add(Float.valueOf(tem[2]));
                sdList.add(Float.valueOf(tem[3]));
                cnt++;
                if (cnt%1000000 == 0) System.out.println(cnt);
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        float[] mean = meanList.toArray();
        float[] sd = sdList.toArray();
        int size = 10000;
        double[] mValue = new double[size];
        double[] sValue = new double[size];
        for (int i = 0; i < size; i++) {
            int index = (int)(Math.random()*mean.length);
            if (mean[index] == 0) {
                i--;
                continue;
            }
            mValue[i] = mean[index];
            sValue[i] = sd[index];
        }
//        ScatterPlot sp = new ScatterPlot(mValue, sValue);
//        sp.setXLab("Relative Depth");
//        sp.setYLab("SD");
//        sp.setXLim(0, 2);
//        sp.setYLim(0, 1);
//        sp.showGraph();
        try{
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("SiteID\tRelativeDepth\tSD");
            bw.newLine();
            for (int i = 0; i < mValue.length; i++) {
                bw.write(String.valueOf(i+1)+"\t"+String.valueOf(mValue[i])+"\t"+String.valueOf(sValue[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void estimateRelativeDepthAndVariance() {
        String chrDepthDirS = "/workdir/fl262/chrDepth/";
        String taxaDepthFileS = "/workdir/fl262/taxaDepth.txt";
        String outDirS = "/workdir/fl262/relativeDepth/";
        Table t = new Table(taxaDepthFileS);
        int[] taxaDepth = new int[t.getColumnNumber()];
        for (int i = 0; i < taxaDepth.length; i++) taxaDepth[i] = Integer.valueOf(t.content[0][i]);
        File[] fs = new File(chrDepthDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String outfileS = new File(outDirS, f.getName().replaceFirst("taxaDepth.txt", "DAndV.txt")).getAbsolutePath();
            try {
                BufferedReader br = IoUtils.getTextReader(f.getAbsolutePath());
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                String header = br.readLine();
                bw.write("Chr\tPos\tMeanDepth\tSD");
                bw.newLine();
                String temp = null;
                double[] depth = new double[taxaDepth.length];
                double mean;
                double sd;
                
                while ((temp = br.readLine()) != null) {
                    String[] tem = temp.split("\t");
                    for (int i = 0; i < depth.length; i++) {
                        depth[i] = Double.valueOf(tem[i+2])/taxaDepth[i];
                    }
                    DescriptiveStatistics st = new DescriptiveStatistics(depth);
                    mean = st.getMean();
                    sd = st.getStandardDeviation();
                    StringBuilder sb = new StringBuilder();
                    sb.append(tem[0]).append("\t").append(tem[1]).append("\t").append((float)mean).append("\t").append((float)sd);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
    
    public void estimateTaxaDepth () {
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\badRefSite\\chromDepth\\chr001.taxaDepth.txt";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\badRefSite\\taxaDepth\\taxaDepth.txt";
        try {
            int size = 1000000;
            BufferedReader br = IoUtils.getTextReader(infileS);
            String temp = br.readLine();
            String[] tem = temp.split("\t");
            String[] taxa = new String[tem.length-2];;
            for (int i = 0; i < taxa.length; i++) {
                taxa[i] = tem[i+2];
            }
            int[][] depth = new int[taxa.length][size];
            for (int i = 0; i < size; i++) {
                tem = br.readLine().split("\t");
                for (int j = 0; j < taxa.length; j++) {
                    depth[j][i] = Integer.valueOf(tem[j+2]);
                }
            }
            int[] peakDepth = new int[taxa.length];
            for (int i = 0; i < taxa.length; i++) {
                int[] d = new int[200];
                for (int j = 0; j < size; j++) {
                    if (depth[i][j] >= d.length) continue;
                    d[depth[i][j]]++;
                }
                int max = 0;
                int maxDepth = 0;
                for (int j = 1; j < d.length; j++) {
                    if (d[j] > max) {
                        max = d[j];
                        maxDepth = j;
                    }
                }
                peakDepth[i] = maxDepth;
            }
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            StringBuilder sb1 = new StringBuilder();
            StringBuilder sb2 = new StringBuilder();
            for (int i = 0; i < taxa.length; i++) {
                sb1.append(taxa[i]).append("\t");
                sb2.append(peakDepth[i]).append("\t");
            }
            sb1.deleteCharAt(sb1.length()-1);
            sb2.deleteCharAt(sb2.length()-1);
            bw.write(sb1.toString()); bw.newLine();
            bw.write(sb2.toString()); bw.newLine();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    } 
    
    private void mkChrDepthFiles () {
        String bamDirS = "/workdir/fl262/bam/";
        String pileupDirS = "/workdir/fl262/pileup/";
        String chrDepthChunkDirS = "/workdir/fl262/chrDepthChunk/";
        String chrDepthDirS = "/workdir/fl262/chrDepth/";
        String referenceFileS = "/workdir/fl262/reference/cassavaV6_23Chr.fa";
        int regionSize = 100000;
        new File (pileupDirS).mkdir();
        new File (chrDepthChunkDirS).mkdir();
        new File (chrDepthDirS).mkdir();
        HashMap<String, String[]> taxaBamMap = this.getTaxaBamMap(bamDirS);
        HashMap<String, String> bamPileupMap = this.getBamPileupMap(pileupDirS, bamDirS);
        Set<String> taxaSet = taxaBamMap.keySet();
        String[] taxaNames = taxaSet.toArray(new String[taxaSet.size()]);
        Arrays.sort(taxaNames);
        Fasta f = new Fasta(referenceFileS);
        int chrNum = 20;
        for (int i = 0; i < chrNum; i++) {
            int currentChr = i + 1;
            int chrLength = f.getSeqLength(i);
            int[][] regionBound = this.creatRegions(currentChr, regionSize, 1, chrLength);
            for (int j = 0; j < regionBound.length; j++) {
                String[] commands = this.creatPileupCommand(currentChr, regionBound[j][0], regionBound[j][1], taxaNames, taxaBamMap, bamPileupMap, referenceFileS, bamDirS, pileupDirS);
                this.performPileUp(commands);
                this.outputDepthMatrix(regionBound[j][0], regionBound[j][1], taxaNames, taxaBamMap, bamPileupMap, chrDepthChunkDirS, currentChr);
                System.out.println("Chromosome: "+String.valueOf(currentChr)+"\tPosition:"+String.valueOf(regionBound[j][0])+" is finished.");
            }
            this.mergeDepth(taxaNames, chrDepthChunkDirS, chrDepthDirS, currentChr);
            try {
                FileUtils.cleanDirectory(new File(pileupDirS));
                FileUtils.cleanDirectory(new File(chrDepthChunkDirS));
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
    }  
    
    private void mergeDepth (String[] taxaNames, String chrDepthChunkDirS, String chrDepthDirS, int currentChr) {
        System.out.println("\nMerging depth chunks");
        int cnt = 0;
        try {
            String outfileS = "chr"+FStringUtils.getNDigitNumber(3, currentChr)+".taxaDepth.txt";
            outfileS = new File(chrDepthDirS, outfileS).getAbsolutePath();
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Chr\tPos");
            for (int i = 0; i < taxaNames.length; i++) {
                bw.write("\t"+taxaNames[i]);
            }
            bw.newLine();
            File[] fs = new File (chrDepthChunkDirS).listFiles();
            HashMap<Integer, File> posFileMap = new HashMap();
            for (int i = 0; i < fs.length; i++) {
                int key = Integer.valueOf(fs[i].getName().split("_")[1]);
                posFileMap.put(key, fs[i]);
            }
            Set<Integer> posSet = posFileMap.keySet();
            Integer[] startPoses = posSet.toArray(new Integer[posSet.size()]);
            Arrays.sort(startPoses);
            for (int i = 0; i < startPoses.length; i++) {
                BufferedReader br = IoUtils.getTextReader(posFileMap.get(startPoses[i]).getAbsolutePath());
                String temp = br.readLine();
                while ((temp = br.readLine()) != null) {
                    bw.write(temp);
                    bw.newLine();
                    cnt++;
                }
                br.close();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println(String.valueOf(cnt)+" depth written");
    }
    
    private void outputDepthMatrix (int startPos, int endPos, String[] taxaNames, HashMap<String, String[]> taxaBamMap, HashMap<String, String> bamPileupMap, String chrDepthDirS, int currentChr) {
        int regionLength = endPos - startPos +1;
        ConcurrentHashMap<String, int[]> bamDepthMap = new ConcurrentHashMap();
        List<String> taxaList = Arrays.asList(taxaNames);
        taxaList.parallelStream().forEach(taxa -> {
            int taxaIndex = Arrays.binarySearch(taxaNames, taxa);
            String[] bams = taxaBamMap.get(taxa);
            List<String> bamList = Arrays.asList(bams);
            bamList.parallelStream().forEach(bam -> {
                String pileupFileS = bamPileupMap.get(bam);
                int[] depth = new int[regionLength];
                try {
                    BufferedReader br = IoUtils.getTextReader(pileupFileS);
                    String temp;
                    while ( (temp = br.readLine()) != null) {
                        String[] tem = temp.split("\t");
                        int siteIndex = Integer.valueOf(tem[1])-startPos;
                        int d = Integer.valueOf(tem[3]);
                        depth[siteIndex] = d;
                    }
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
                bamDepthMap.put(bam, depth);
            });
        });
        String outfileS = this.buildFileName(currentChr, startPos, endPos, chrDepthDirS, ".taxaDepth.txt");
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("chr\tPos");
            for (int i=  0; i < taxaNames.length; i++) {
                bw.write("\t"+taxaNames[i]);
            }
            bw.newLine();
            int[][] td = new int[taxaNames.length][regionLength];
            for (int i = 0; i < taxaNames.length; i++) {
                String[] bams = taxaBamMap.get(taxaNames[i]);
                for (int j = 0; j < bams.length; j++) {
                    int[] d = bamDepthMap.get(bams[j]);
                    for (int k = 0; k < d.length; k++) {
                        td[i][k]+=d[k];
                    }
                }
            }
            for (int i = 0; i < regionLength; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(currentChr).append("\t").append(startPos+i);
                for (int j = 0; j < td.length; j++) {
                    sb.append("\t").append(td[j][i]);
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
    
    private String buildFileName (int currentChr, int startPos, int endPos, String dirS, String suffix) {
        StringBuilder sb = new StringBuilder();
        sb.append("chr").append(FStringUtils.getNDigitNumber(3, currentChr)).append("_").append(startPos).append("_").append(endPos).append(suffix);
        return new File (dirS, sb.toString()).getAbsolutePath();
    }
    
    private int[][] creatRegions (int currentChr, int regionSize, int startPos, int endPos) {
        int[][] regionBound = FArrayUtils.getSubsetsIndicesBySubsetSize(endPos-startPos+1, regionSize);
        for (int i = 0; i < regionBound.length; i++) {
            regionBound[i][0]+=startPos;
            regionBound[i][1]+=startPos;
            regionBound[i][1]--;
        }
        return regionBound;
    }
    
    private void performPileUp (String[] commands) {
        int coreNumber = Runtime.getRuntime().availableProcessors();
        int[][] batchIndex = FArrayUtils.getSubsetsIndicesBySubsetSize(commands.length, coreNumber);
        List<String> commandList = Arrays.asList(commands);
        for (int i = 0; i < batchIndex.length; i++) {
            List<String> subList = commandList.subList(batchIndex[i][0], batchIndex[i][1]);
            subList.parallelStream().forEach(element -> {
                try {
                    Runtime run = Runtime.getRuntime();
                    Process p = run.exec(element);
                    p.waitFor();
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            });
        }
    }
    
    private String[] creatPileupCommand (int currentChr, int startPos, int endPos, String[] taxaNames, HashMap<String, String[]> taxaBamMap, HashMap<String, String> bamPileupMap, String referenceFileS, String bamDirS, String pileupDirS) {
        ArrayList<String> commandList = new ArrayList();
        for (int i = 0; i < taxaNames.length; i++) {
            String[] bams = taxaBamMap.get(taxaNames[i]);
            for (int j = 0; j < bams.length; j++) {
                StringBuilder sb = new StringBuilder();     
                String pileupFileS = bamPileupMap.get(bams[j]);
                sb.append("samtools mpileup -A -B -q 30 -Q 10 -f ").append(referenceFileS).append(" ").append(bams[j]).append(" -r ");
                sb.append(currentChr).append(":").append(startPos).append("-").append(endPos).append(" -o ").append(pileupFileS);
                commandList.add(sb.toString());
            }
        }
        return commandList.toArray(new String[commandList.size()]);
    }
    
    private HashMap<String, String> getBamPileupMap (String pileupDirS, String bamDirS) {
        File[] bams = IoUtils.listFilesEndsWith(new File(bamDirS).listFiles(), "bam");
        HashMap<String, String> bamPileupMap = new HashMap();
        for (int i = 0; i < bams.length; i++) {
            bamPileupMap.put(bams[i].getAbsolutePath(), new File (pileupDirS, bams[i].getName().replace(".bam", ".pileup.txt")).getAbsolutePath());
        }
        return bamPileupMap;
    }
    
    private HashMap<String, String[]> getTaxaBamMap (String bamDirS) {
        File[] bams = IoUtils.listFilesEndsWith(new File(bamDirS).listFiles(), "bam");
        HashSet<String> taxaSet = new HashSet();
        for (int i = 0; i < bams.length; i++) {
            taxaSet.add(bams[i].getName().split("_")[0]);
        }
        String[] taxa = taxaSet.toArray(new String[taxaSet.size()]);
        Arrays.sort(taxa);
        ArrayList<String>[] fileList = new ArrayList[taxa.length];
        for (int i = 0; i < fileList.length; i++) fileList[i] = new ArrayList();
        for (int i = 0; i < bams.length; i++) {
            String query = bams[i].getName().split("_")[0];
            int index = Arrays.binarySearch(taxa, query);
            fileList[index].add(bams[i].getAbsolutePath());
        }
        HashMap<String, String[]> taxaBamMap = new HashMap();
        for (int i = 0; i < taxa.length; i++) {
            taxaBamMap.put(taxa[i], fileList[i].toArray(new String[fileList[i].size()]));
        }
        return taxaBamMap;
    }
    
    public void cloneDepthPipe () {
        this.convertToJPileup();
        this.mkPrefilterDepth();
    }
    
    public void mkPrefilterDepth () {
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\snpFilter\\prefilterDepth\\I000070_chr001.jp.bin";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\snpFilter\\prefilterDepth\\I000070_prefilter.depth.pdf";
        int size = 20000;
        JPileUp jp = new JPileUp(infileS, IOFileFormat.Binary);
        double[] depth = new double[size];
        for (int i = 0; i < size; i++) {
            int index = (int)(Math.random()*jp.getSiteNumber());
            depth[i] = jp.getSiteDepth(index);
        }
        DensityPlot d = new DensityPlot(depth);
        d.setTitle("Read depth of I000070 aftet prefiltering");
        d.setXLab("Depth");
        d.setYLab("Density");
        d.setXLim(0, 60);
        d.setYLim(0, 0.1);
        d.setSmoothN(10000);
        d.saveGraph(outfileS);
    }
    
    public void convertToJPileup () {
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\snpFilter\\prefilterDepth\\I000070_chr001.pileup.txt";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\snpFilter\\prefilterDepth\\I000070_chr001.jp.bin";
        JPileUp jp = new JPileUp();
        jp.readFromPileup(infileS);
        jp.writeFile(outfileS, IOFileFormat.Binary);
    }
    
    public void filterTestPipe () {
        //this.makeRandomSite();
        this.depthDis(); 
        //this.makePCAPlot();
    }
    
    public void makePCAPlot () {
        String pcInfileS = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\pca\\cassava_pcs.txt";
        String annotationFileS = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\pca\\cassavaWGSlist_annotation.txt";
        String pdfFileS = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\pca\\cassava_pca.pdf";
        Table t = new Table (annotationFileS);
        HashMap<String,String> taxaClassMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaClassMap.put(t.content[i][0], t.content[i][3]);
        }
        t = new Table (pcInfileS);
        TDoubleArrayList xList = new TDoubleArrayList();
        TDoubleArrayList yList = new TDoubleArrayList();
        ArrayList<String> classList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            String key = t.content[i][0];
            String value = taxaClassMap.get(t.content[i][0]);
            System.out.println(key+"\t"+value);
            if (value == null) continue;
            xList.add(t.getDoubleValue(i, 1));
            yList.add(t.getDoubleValue(i, 2));
            classList.add(value);
        }
        ScatterPlotMultiClass sc = new ScatterPlotMultiClass(xList.toArray(), yList.toArray(), classList.toArray(new String[classList.size()]));
        sc.setXLab("PC1");
        sc.setYLab("PC2");
        sc.setTitle("PCA plot of cassava WGS clones");
        sc.saveGraph(pdfFileS);
    }
    
    public void depthDis () {
//        String infileS = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\vcfChr1\\chr001_random.raw.VCF.txt";
//        String pdfFileS = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\snpFilter\\depth\\chr001_random.raw.depth.pdf";
        String infileS = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\vcfChr1\\chr001_random.stage01.VCF.txt";
        String pdfFileS = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\snpFilter\\depth\\chr001_random.stage01.depth.pdf";
        Table t = new Table (infileS);
        int size = 10000;
        double[] v = new double[size];
        for (int i = 0; i < size; i++) {
            v[i] = Double.valueOf(t.content[i][4].split(";")[0].replaceFirst("DP=", ""));
        }
        DensityPlot d =new DensityPlot(v);
        d.setSmoothN(10000);
        d.setXLim(0, 20000);
        d.setXLab("Depth");
        d.setYLab("Density");
        d.saveGraph(pdfFileS);
    }
    
    public void makeRandomSite () {
//        String infileS = "M:\\production\\cassava\\hapmap\\raw\\chr001.VCF.txt";
//        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\vcfChr1\\chr001_random.raw.VCF.txt";
        String infileS = "M:\\production\\cassava\\hapmap\\stage_01\\chr001.stage01.VCF.txt";
        String outfileS = "M:\\pipelineTest\\cassava\\wgs\\snpDiscovery\\vcfChr1\\chr001_random.stage01.VCF.txt";
        int size = 20000;  
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            int cnt = 0;
            String temp = br.readLine();
            while ((temp = br.readLine()) != null) {
                cnt++;
                if (cnt%100000 == 0) System.out.println(cnt);
            }
            int snpNum = cnt;
            br.close();
            int[][] subsets = FArrayUtils.getSubsetsIndicesBySubsetNumber(snpNum, size);
            int[] index = new int[size];
            for (int i = 0; i < index.length; i++) index[i] = subsets[i][0];
            br = IoUtils.getTextReader(infileS);
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write(br.readLine());
            bw.newLine();
            cnt = -1;
            while ((temp = br.readLine()) != null) {
                cnt++;
                if (Arrays.binarySearch(index, cnt) < 0) {
                    continue;
                }
                bw.write(temp);
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
