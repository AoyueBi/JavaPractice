/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maf;

import format.Table;
import gnu.trove.list.array.TDoubleArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import stats.FSimpleRegression;
import stats.Transformation;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class PopGenAdaptation {
    
    public PopGenAdaptation () {
        //this.mergeFstAndSNPStatistics();
        //this.mergeTwoFile();
        //this.splitByGroup();
        //this.foldEnrichmentSIFT();
        //this.foldEnrichmentGERP();
        //this.correctFst4Maf();
        this.batchPdf();
    }
    
    public void batchPdf () {
        String fstDirS = "M:\\production\\maf\\popgen\\adaptation\\byGroup\\corrected";
        String rScirptDirS = "M:\\production\\maf\\popgen\\adaptation\\byGroup\\script";
        String tempRFileS = "M:\\production\\maf\\popgen\\adaptation\\byGroup\\tempR.R";
        File[] fstFs = new File(fstDirS).listFiles();
        File[] rFs = new File(rScirptDirS).listFiles();
        for (int i = 0; i < rFs.length; i++) {
            ArrayList<String> lList = new ArrayList();
            try {
                BufferedReader br = IoUtils.getTextReader(rFs[i].getAbsolutePath());
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    lList.add(temp);
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            for (int j = 0; j < fstFs.length; j++) {
                try {
                    String groupName = fstFs[j].getName().replaceFirst(".txt", "");
                    BufferedWriter bw = IoUtils.getTextWriter(tempRFileS);
                    for (int k = 0; k < lList.size(); k++) {
                        bw.write(lList.get(k).replaceFirst("FileName", groupName));
                        bw.newLine();
                    }
                    bw.flush();
                    bw.close();
                    StringBuilder sb = new StringBuilder("rscript");
                    sb.append(" ").append(tempRFileS);
                    String command = sb.toString();
                    System.out.println(command);
                    Runtime rt = Runtime.getRuntime();
                    Process p = rt.exec(command);
                    p.waitFor();
                    System.out.println(String.valueOf(i)+"\t"+String.valueOf(j));
                    System.out.println(rFs[i].getName()+"\t"+fstFs[i].getName());
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
            
        }
    }
    
    public void correctFst4Maf () {
        String infileDirS = "M:\\production\\maf\\popgen\\adaptation\\byGroup\\selected\\";
        String outDirS = "M:\\production\\maf\\popgen\\adaptation\\byGroup\\corrected\\";
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            Table t = new Table(fs[i].getAbsolutePath());
            double[] fst = t.getDoubleArrayByColumn(2);
            double[] maf = t.getDoubleArrayByColumn(5);
            FSimpleRegression fsr = new FSimpleRegression(maf, fst);

            double[] predictFst = fsr.getResidual();
            double[] zs = Transformation.getZScore(predictFst);
            try {
                String outfileS = new File (outDirS, fs[i].getName()).getAbsolutePath();
                BufferedReader br = IoUtils.getTextReader(fs[i].getAbsolutePath());
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                String temp = br.readLine();
                bw.write(temp+"\tCorrectFst\tZScore");
                bw.newLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    bw.write(temp+String.valueOf(predictFst[cnt])+"\t"+String.valueOf(zs[cnt]));
                    bw.newLine();
                    cnt++;
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
    }
    
    public void foldEnrichmentGERP () {
        String inputDirS = "M:\\production\\maf\\popgen\\adaptation\\byGroup\\selected\\";
        String outfileS = "M:\\production\\maf\\popgen\\adaptation\\foldEnrichment\\foldFstGerp.txt";
        File[] fs = new File(inputDirS).listFiles();
        double step = 0.5;
        int numOfGroup =(int) Math.ceil(4.81/step);
        double[] bound = new double[numOfGroup+1] ;
        for (int i = 0; i < bound.length; i++) {
            bound[i] = i*step;
        }
        double[][] folds = new double[fs.length][];
        StringBuilder sb = new StringBuilder("GERPGroup");
        for (int i = 0; i < fs.length; i++) {
            Table t = new Table(fs[i].getAbsolutePath());
            String population = fs[i].getName().replaceFirst(".txt", "");
            folds[i] = this.getFstFoldGERP(t, bound);
            sb.append("\t").append(population);
        }
        sb.append("\tMean\tError");
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < numOfGroup; i++) {
                sb = new StringBuilder();
                sb.append(bound[i]);
                TDoubleArrayList dList = new TDoubleArrayList();
                for (int j = 0; j < folds.length; j++) {
                    sb.append("\t").append(folds[j][i]);
                    if (Double.isNaN(folds[j][i])) continue;
                    if (Double.isInfinite(folds[j][i])) continue;
                    dList.add(folds[j][i]);
                }
                double[] d = dList.toArray();
                DescriptiveStatistics ds = new DescriptiveStatistics(d);
                double mean = ds.getMean();
                double sd = ds.getStandardDeviation();
                double se = sd/Math.sqrt(d.length);
                sb.append("\t").append(mean).append("\t").append(se);
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
    
    double[] getFstFoldGERP (Table t, double[] bound) {
        double[] fold = new double[bound.length-1];
        double[] gerpMean = new double[fold.length];
        double[] gerpError = new double[fold.length];
        int totalHighFst = 0;
        int[] groupHighFst = new int[fold.length];
        int[] groupSize = new int[fold.length];
        double[] values = new double[t.getRowNumber()];
        TDoubleArrayList[] gerpList = new TDoubleArrayList[fold.length];
        for (int i = 0; i < gerpList.length; i++) {
            gerpList[i] = new TDoubleArrayList();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            values[i] = t.getDoubleValue(i, 2);
        }
        Arrays.sort(values);
        double fstThresh = values[(int)((double)values.length*0.99)];
        System.out.println(fstThresh);
        int cnt = 0;
        for (int i = 0; i < t.getRowNumber(); i++) {
            double fst = t.getDoubleValue(i, 2);
            double gerp = t.getDoubleValue(i, 11);
            if (gerp < 0) continue;
            cnt++;
            int index = Arrays.binarySearch(bound, gerp);
            if (index < 0) index = -index-2;
            if (index >= fold.length) {
                index--;
            }
            gerpList[index].add(fst);
            groupSize[index]++;
            if (fst<fstThresh) continue;
            groupHighFst[index]++;
            totalHighFst++;
        }
        double base = (double)totalHighFst/cnt;
        for (int i = 0; i < fold.length; i++) {
            fold[i] = Math.log((double)groupHighFst[i]/groupSize[i]/base)/Math.log(2);
            values = gerpList[i].toArray();
            DescriptiveStatistics d = new DescriptiveStatistics(values);
            gerpMean[i] = d.getMean();
            gerpError[i] = d.getStandardDeviation()/Math.sqrt(values.length);
        }
//        return fstMean;
        return fold;
    }
    
    public void foldEnrichmentSIFT () {
        String inputDirS = "M:\\production\\maf\\popgen\\adaptation\\byGroup\\selected\\";
        String outfileS = "M:\\production\\maf\\popgen\\adaptation\\foldEnrichment\\foldFstSIFT.txt";
        File[] fs = new File(inputDirS).listFiles();
        double step = 0.05;
        int numOfGroup =(int) Math.ceil(1/step);
        double[] bound = new double[numOfGroup+1] ;
        for (int i = 0; i < bound.length; i++) {
            bound[i] = i*step;
        }
        double[][] folds = new double[fs.length][];
        StringBuilder sb = new StringBuilder("SIFTGroup");
        for (int i = 0; i < fs.length; i++) {
            Table t = new Table(fs[i].getAbsolutePath());
            String population = fs[i].getName().replaceFirst(".txt", "");
            folds[i] = this.getFstFoldSIFT(t, bound);
            sb.append("\t").append(population);
        }
        sb.append("\tMean\tError");
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < numOfGroup; i++) {
                sb = new StringBuilder();
                sb.append(bound[i]);
                TDoubleArrayList dList = new TDoubleArrayList();
                for (int j = 0; j < folds.length; j++) {
                    sb.append("\t").append(folds[j][i]);
                    if (Double.isNaN(folds[j][i])) continue;
                    if (Double.isInfinite(folds[j][i])) continue;
                    dList.add(folds[j][i]);
                }
                double[] d = dList.toArray();
                DescriptiveStatistics ds = new DescriptiveStatistics(d);
                double mean = ds.getMean();
                double sd = ds.getStandardDeviation();
                double se = sd/Math.sqrt(d.length);
                sb.append("\t").append(mean).append("\t").append(se);
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
    
    double[] getFstFoldSIFT (Table t, double[] bound) {
        double[] fold = new double[bound.length-1];
        double[] fstMean = new double[fold.length];
        double[] fstError = new double[fold.length];
        int totalHighFst = 0;
        int[] groupHighFst = new int[fold.length];
        int[] groupSize = new int[fold.length];
        double[] values = new double[t.getRowNumber()];
        TDoubleArrayList[] fstList = new TDoubleArrayList[fold.length];
        for (int i = 0; i < fstList.length; i++) {
            fstList[i] = new TDoubleArrayList();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            values[i] = t.getDoubleValue(i, 2);
        }
        Arrays.sort(values);
        double fstThresh = values[(int)((double)values.length*0.99)];
        System.out.println(fstThresh);
        for (int i = 0; i < t.getRowNumber(); i++) {
            double fst = t.getDoubleValue(i, 2);
            double sift = t.getDoubleValue(i, 9);
            int index = Arrays.binarySearch(bound, sift);
            if (index < 0) index = -index-2;
            if (index >= fold.length) {
                index--;
            }
            fstList[index].add(fst);
            groupSize[index]++;
            if (fst<fstThresh) continue;
            groupHighFst[index]++;
            totalHighFst++;
        }
        double base = (double)totalHighFst/t.getRowNumber();
        for (int i = 0; i < fold.length; i++) {
            fold[i] = Math.log((double)groupHighFst[i]/groupSize[i]/base)/Math.log(2);
            values = fstList[i].toArray();
            DescriptiveStatistics d = new DescriptiveStatistics(values);
            fstMean[i] = d.getMean();
            fstError[i] = d.getStandardDeviation()/Math.sqrt(values.length);
        }
//        return fstMean;
        return fold;
    }
    
    public void splitByGroup () {
        String infileS = "M:\\production\\maf\\popgen\\adaptation\\byClass\\all.txt";
        String outDirS = "M:\\production\\maf\\popgen\\adaptation\\byGroup\\";
        Table t = new Table (infileS);
        String[] info = new String[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            StringBuilder sb = new StringBuilder();
            for (int j = 17; j < t.getColumnNumber(); j++) {
                sb.append("\t").append(t.content[i][j]);
            }
            info[i] = sb.toString();
        }
        for (int i = 0; i < 15; i++) {
            String outfileS = new File (outDirS, t.header[i+2]+".txt").getAbsolutePath();
            try {
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                bw.write("Chr\tPos\tFst\tRatioFstVsMAF\tDAF	MinorAlleleFrequency	SiteDepth	HetCount	MutationClass	Sift	GerpTreeLength	Gerp	UScore	Transcripts");
                bw.newLine();
                for (int j = 0; j < t.getRowNumber(); j++) {
                    if (t.content[j][i+2].startsWith("N")) continue;
                    StringBuilder sb = new StringBuilder(t.content[j][0]);
                    double fst = Double.valueOf(t.content[j][i+2]);
                    List<String> l = FStringUtils.fastSplit(info[j]);
                    double ratio = fst/Double.valueOf(l.get(2));
                    sb.append("\t").append(t.content[j][1]).append("\t").append(t.content[j][i+2]).append("\t").append(ratio).append(info[j]);
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
    }
    
    public void mergeTwoFile () {
        String synFileS = "M:\\production\\maf\\popgen\\adaptation\\byClass\\synFst.txt";
        String nonSynFileS = "M:\\production\\maf\\popgen\\adaptation\\byClass\\nonSynTolerentFst.txt";
        String delSFileS = "M:\\production\\maf\\popgen\\adaptation\\byClass\\delSFst.txt";
        String mergedFileS = "M:\\production\\maf\\popgen\\adaptation\\byClass\\all.txt";
        Table t = new Table (synFileS);
        Table at = new Table (delSFileS);
        Table nt = t.getMergeTableByRow(at);
        nt = nt.getMergeTableByRow(new Table(nonSynFileS));
        nt.writeTable(mergedFileS);
    }
    
    public void mergeFstAndSNPStatistics () {
        String synFstFileS = "M:\\production\\maf\\popgen\\paremeters\\fst\\fstSynonymous\\siteFst.txt";
        String nonSynFstFileS = "M:\\production\\maf\\popgen\\paremeters\\fst\\fstNonSynonymousTolerent\\siteFst.txt";
        String delSFstFileS = "M:\\production\\maf\\popgen\\paremeters\\fst\\fstDeleteriousS\\siteFst.txt";
        String delSGFstFileS = "M:\\production\\maf\\popgen\\paremeters\\fst\\fstDeleteriousSG\\siteFst.txt";
        String annotationDirS = "M:\\production\\maf\\annotations\\siftScore\\003_hmp321SiftGerpUScore\\";
        String outSynFileS = "M:\\production\\maf\\popgen\\adaptation\\byClass\\synFst.txt";
        String outNonSynFileS = "M:\\production\\maf\\popgen\\adaptation\\byClass\\nonSynTolerentFst.txt";
        String outDelSFileS = "M:\\production\\maf\\popgen\\adaptation\\byClass\\delSFst.txt";
        String outDelSGFileS = "M:\\production\\maf\\popgen\\adaptation\\byClass\\delSGFst.txt";
        File[] fs = new File(annotationDirS).listFiles();
        HashMap<Integer, String>[] posInfoMap = new HashMap[fs.length];
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            int chrIndex = Integer.valueOf(f.getName().replaceFirst(".hmpSiftGerpUScore.txt", "").replaceFirst("chr", ""))-1;
            posInfoMap[chrIndex] = new HashMap();
            try {
                BufferedReader br = IoUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                while ((temp = br.readLine()) != null) {
                    List<String> l = FStringUtils.fastSplit(temp);
                    int pos = Integer.valueOf(l.get(1));
                    String daf = "NA";
                    if (l.get(4).equals(l.get(5))) daf = l.get(8);
                    else if (l.get(4).equals(l.get(6))) daf = l.get(7);
                    StringBuilder sb = new StringBuilder(daf);
                    for (int i = 8; i < l.size(); i++) {
                        sb.append("\t").append(l.get(i));
                    }
                    posInfoMap[chrIndex].put(pos, sb.toString());
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        this.outputFstAndSNPStatistics(posInfoMap, synFstFileS, outSynFileS);
        this.outputFstAndSNPStatistics(posInfoMap, nonSynFstFileS, outNonSynFileS);
        this.outputFstAndSNPStatistics(posInfoMap, delSFstFileS, outDelSFileS);
        this.outputFstAndSNPStatistics(posInfoMap, delSGFstFileS, outDelSGFileS);
    }
    
    private void outputFstAndSNPStatistics (HashMap<Integer, String>[] posInfoMap, String fstFileS, String ouputFileS) {
        try {
            BufferedReader br = IoUtils.getTextReader(fstFileS);
            BufferedWriter bw = IoUtils.getTextWriter(ouputFileS);
            String header = "Chr	Pos	China_specific_VS_Mixed	China_specific_VS_Non_stiff_stalk	China_specific_VS_Stiff_stalk	China_specific_VS_Teosinte	China_specific_VS_Tropical_subtropical	Mixed_VS_Non_stiff_stalk	Mixed_VS_Stiff_stalk	Mixed_VS_Teosinte	Mixed_VS_Tropical_subtropical	Non_stiff_stalk_VS_Stiff_stalk	Non_stiff_stalk_VS_Teosinte	Non_stiff_stalk_VS_Tropical_subtropical	Stiff_stalk_VS_Teosinte	Stiff_stalk_VS_Tropical_subtropical	Teosinte_VS_Tropical_subtropical";
            header = header+"\tDAF\t"+"MinorAlleleFrequency	SiteDepth	HetCount	MutationClass	Sift	GerpTreeLength	Gerp	UScore	Transcripts";
            bw.write(header);
            bw.newLine();
            String temp = br.readLine();
            while ((temp = br.readLine()) != null) {
                List<String> l = FStringUtils.fastSplit(temp);
                int chrIndex = Integer.valueOf(l.get(0))-1;
                StringBuilder sb = new StringBuilder(temp);
                int pos = Integer.valueOf(l.get(1));
                String info = posInfoMap[chrIndex].get(pos);
                if (info == null) continue;
                sb.append("\t").append(info);
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
}
