/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.cassava.wgs;

import format.Range;
import format.Ranges;
import format.Table;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import graphcis.r.LineChart;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
class Heterozygosity {
    
    public Heterozygosity () {
        //this.getSubTaxaHMP();
        //this.calculateHeterozygosity();
        //this.calWholegenomeHeterozygosity();
        
        //this.filterHighCoverageSite();
        this.calHetsInRegions();
        //this.hetsInRegionFigures();
    }
    
    public void hetsInRegionFigures () {
        String inDirS = "M:\\production\\cassava\\hmp\\highCoverageHmp\\homoTaxaInRegion\\result\\";
        String outDirS = "M:\\production\\cassava\\hmp\\highCoverageHmp\\homoTaxaInRegion\\pdf\\";
        File[] fs = new File(inDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            Table t = new Table(fs[i].getAbsolutePath());
            double[][] value = new double[t.getColumnNumber()-3][t.getRowNumber()];
            double[] x = new double[t.getRowNumber()];
            String[] taxaNames = new String[t.getColumnNumber()-3];
            for (int j = 0; j < taxaNames.length; j++) {
                taxaNames[j] = t.header[j+3];
            }
            for (int j = 0; j < t.getRowNumber(); j++) {
                x[j] = t.getDoubleValue(j, 1);
                for (int k = 0; k < value.length; k++) {
                    value[k][j] = t.getDoubleValue(j, 3+k);
                }
            }
            LineChart l  =new LineChart(x, value, taxaNames);
            l.showOnlyLine();
            l.setIfSmoothLine(true);
            String outfileS = new File(outDirS, fs[i].getName().replaceFirst(".txt", ".pdf")).getAbsolutePath();
            l.saveGraph(outfileS);
        }
    }
    
    public void calHetsInRegions () {
        String infileDirS = "Q:\\D\\Cassava\\hapmapFiles_V6\\hapmapCalls\\filtered_highCov_MAF5\\";
        String outDirS = "M:\\production\\cassava\\hmp\\highCoverageHmp\\homoTaxaInRegion\\result\\";
        String regionFileS = "M:\\production\\cassava\\hmp\\highCoverageHmp\\homoTaxaInRegion\\regionSelected.txt";
        int windowSize = 1000;
        int step = 500;
        int numberOfSelTaxa = 5;
        Table t = new Table(regionFileS);
        ArrayList<Range> periList = new ArrayList();
        ArrayList<Range> armList = new ArrayList();
        String[] anno = new String[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            Range r = new Range (Integer.valueOf(t.content[i][1]), t.getIntValue(i, 2), t.getIntValue(i, 3)+1);
            anno[i] = t.content[i][0];
            if (t.content[i][0].startsWith("Pe")) {
                periList.add(r);
            }
            else {
                armList.add(r);
            }
        }
        File[] fs = IoUtils.listFilesEndsWith(new File(infileDirS).listFiles(), "hmp.txt");
        Ranges someRange = new Ranges(periList, "pericentromere");
        calHetsInRegion(fs, someRange, numberOfSelTaxa,  windowSize, step, outDirS, "pericentromere");
        someRange = new Ranges(armList, "arm");
        calHetsInRegion(fs, someRange, numberOfSelTaxa,  windowSize, step, outDirS, "chromosome_arm");
    }
    
     public void calHetsInRegion(File[] fs, Ranges someRange, int numberOfSelTaxa, int windowSize, int step, String outDirS, String category) {
        int[] selChr = someRange.getChromosomes();
        for (int i = 0; i < fs.length; i++) {
            GenotypeTable gt = ImportUtils.readFromHapmap(fs[i].getAbsolutePath());
            int cuChr = gt.chromosome(0).getChromosomeNumber();
            int hit = Arrays.binarySearch(selChr, gt.chromosome(0).getChromosomeNumber());
            if (hit < 0) continue;
            for (int j = 0; j < someRange.getRangeNumber(); j++) {
                if (cuChr != someRange.getRangeChromosome(j)) continue;
                int startSiteIndex = 0;
                int endSiteIndex = 0;
                for (int k = 0; k < gt.numberOfSites(); k++) {
                    if (gt.physicalPositions()[k] > someRange.getRangeStart(j)) {
                        startSiteIndex = k;
                        break;
                    }
                }
                for (int k = startSiteIndex; k < gt.numberOfSites(); k++) {
                    if (gt.physicalPositions()[k] > someRange.getRangeEnd(j)) {
                        endSiteIndex = k;
                        break;
                    }
                }
                TaxaAndHet[] th = new TaxaAndHet[gt.numberOfTaxa()];
                for (int k = 0; k < gt.numberOfTaxa(); k++) {
                    int hetCnt = 0;
                    for (int u = startSiteIndex; u < endSiteIndex; u++) {
                        if (gt.isHeterozygous(k, u)) hetCnt++;
                    }
                    th[k] = new TaxaAndHet(gt.taxa().get(k), (double)hetCnt/(endSiteIndex-startSiteIndex));
                }
                Arrays.sort(th);
                TDoubleArrayList[] hetList = new TDoubleArrayList[numberOfSelTaxa];
                TIntArrayList[] startList = new TIntArrayList[numberOfSelTaxa];
                TIntArrayList[] endList = new TIntArrayList[numberOfSelTaxa];
                for (int k = 0; k < numberOfSelTaxa; k++) {
                    hetList[k] = new TDoubleArrayList();
                    startList[k] = new TIntArrayList();
                    endList[k] = new TIntArrayList();
                    int taxonIndex = gt.taxa().indexOf(th[k].t);
                    for (int u = startSiteIndex; u < endSiteIndex; u+=step) {
                        int end = u+windowSize;
                        if (end > endSiteIndex) end = endSiteIndex;
                        int hetCnt = 0;
                        for (int v = u; v < end; v++) {
                            if (gt.isHeterozygous(taxonIndex, v)) hetCnt++;
                        }
                        hetList[k].add((double)hetCnt/(end-u));
                        startList[k].add(gt.physicalPositions()[u]);
                        endList[k].add(gt.physicalPositions()[end]);
                    }
                }
                double[][] het = new double[numberOfSelTaxa][];
                int[] start = startList[0].toArray();
                int[] end = endList[0].toArray();
                for (int k = 0; k < het.length; k++) het[k] = hetList[k].toArray();
                String outfileS = category+"_"+String.valueOf(someRange.getRangeChromosome(j))+"_"+String.valueOf(someRange.getRangeStart(j)+"_"+String.valueOf(someRange.getRangeEnd(j)))+"_het.txt";
                try {
                    BufferedWriter bw = IoUtils.getTextWriter(new File(outDirS, outfileS).getAbsolutePath());
                    bw.write("Chr\tStart\tEnd");
                    for (int k = 0; k < numberOfSelTaxa; k++) {
                        bw.write("\t"+th[k].t.getName());
                    }
                    bw.newLine();
                    for (int k = 0; k < het[0].length; k++) {
                        bw.write(gt.chromosomeName(0)+"\t"+String.valueOf(start[k])+"\t"+String.valueOf(end[k]));
                        for (int u = 0; u < het.length; u++) {
                            bw.write("\t"+String.valueOf(het[u][k]));
                        }
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
        
    }
     
    class TaxaAndHet implements Comparable<TaxaAndHet>{
        Taxon t;
        double het;
        public TaxaAndHet (Taxon t, double het) {
            this.t = t;
            this.het = het;
        }

        @Override
        public int compareTo(TaxaAndHet o) {
            if (het < o.het) return -1;
            else if (het > o.het) return 1;
            return 0;
        }
    }
    
    public void filterHighCoverageSite () {
        boolean ifSelectChr = true;
        double mafCut = 0.05;
        double callRateCut = 1;
        String infileDirS = "Q:\\D\\Cassava\\hapmapFiles_V6\\hapmapCalls\\";
        String outfileDirS = "M:\\production\\cassava\\hmp\\highCoverageHmp\\";
        File[] fs = IoUtils.listFilesEndsWith(new File(infileDirS).listFiles(), "vcf");
        int[] chrSelected = {1,2,3,4,5};
        
        for (int i = 0; i < fs.length; i++) {
            String c = fs[i].getName().replace(".vcf", "");
            c = c.replaceFirst(".*_c", "");
            if (ifSelectChr) {
                int hit = Arrays.binarySearch(chrSelected, Integer.valueOf(c));
                if (hit < 0) continue;
            }
            GenotypeTable gt = ImportUtils.readFromVCF(fs[i].getAbsolutePath(), null, true);
            TIntArrayList siteList = new TIntArrayList();
            for (int j = 0; j < gt.numberOfSites(); j++) {
                if (gt.minorAlleleFrequency(j) < mafCut) continue;
                double callrate = (double)gt.totalNonMissingForSite(j)/gt.numberOfTaxa();
                if (callrate < callRateCut) continue;
                siteList.add(j);
            }
//            FilterGenotypeTable sgt = FilterGenotypeTable.getInstance(gt, siteList.toArray());
//            String outfileS = new File(outfileDirS, fs[i].getName().replace(".vcf", ".hmp.txt")).getAbsolutePath();
//            ExportUtils.writeToHapmap(sgt, outfileS);
//            System.out.println(outfileS+" exported");
        }
    }
    
     
    public void calWholegenomeHeterozygosity () {
        String infileDirS = "M:\\production\\cassava\\hmp\\lowHetsTaxa\\";
        String outfileS = "M:\\production\\cassava\\hmp\\lowHetsTaxa\\taxaList\\taxaHeterogosity.txt";
        File[] fs = IoUtils.listFilesEndsWith(new File(infileDirS).listFiles(), "hmp.txt");
        GenotypeTable gt = ImportUtils.readFromHapmap(fs[0].getAbsolutePath());
        String[] taxaNames = new String[gt.numberOfTaxa()];
        int[] hetCount = new int[taxaNames.length];
        for (int i = 0; i < taxaNames.length; i++) {
            taxaNames[i] = gt.taxaName(i);
        }
        int totalSite = 0;
        for (int i = 0; i < fs.length; i++) {
            gt = ImportUtils.readFromHapmap(fs[i].getAbsolutePath());
            for (int j = 0; j < gt.numberOfTaxa(); j++) {
                for (int k = 0; k < gt.numberOfSites(); k++) {
                    if (gt.isHeterozygous(j, k)) hetCount[j]++;
                }
            }
            totalSite+=gt.numberOfSites();
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Taxa\tHeterozygosity");
            bw.newLine();
            for (int i = 0; i < taxaNames.length; i++) {
                bw.write(taxaNames[i]+"\t"+String.valueOf((double)hetCount[i]/totalSite));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void getSubTaxaHMP () {
        double mafCut = 0.05;
        double callRateCut = 0.99;
        String infileDirS = "Q:\\D\\Cassava\\hapmapFiles_V6\\hapmapCalls\\";
        String taxaNameFileS = "M:\\production\\cassava\\hmp\\lowHetsTaxa\\taxaList\\lowHetsTaxa.txt";
        String outfileDirS = "M:\\production\\cassava\\hmp\\lowHetsTaxa\\";
        Table ta = new Table(taxaNameFileS);
        String[] taxaNames = new String[ta.getRowNumber()];
        for (int i = 0; i < taxaNames.length; i++) taxaNames[i] = ta.content[i][0];
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (int i = 0; i < taxaNames.length; i++) {
            Taxon t = new Taxon (taxaNames[i]);
            tlb.add(t);
        }
        TaxaList tList = tlb.build();
        File[] fs = IoUtils.listFilesEndsWith(new File(infileDirS).listFiles(), "vcf");
        for (int i = 0; i < fs.length; i++) {
            GenotypeTable gt = ImportUtils.readFromVCF(fs[i].getAbsolutePath(), null, true);
            TIntArrayList siteList = new TIntArrayList();
            for (int j = 0; j < gt.numberOfSites(); j++) {
                if (gt.minorAlleleFrequency(j) < mafCut) continue;
                double callrate = (double)gt.totalNonMissingForSite(j)/gt.numberOfTaxa();
                if (callrate < callRateCut) continue;
                siteList.add(j);
            }
//            FilterGenotypeTable sgt = FilterGenotypeTable.getInstance(gt, siteList.toArray());
//            GenotypeTable tgt = FilterGenotypeTable.getInstance(sgt, tList);
//            String outfileS = new File(outfileDirS, fs[i].getName().replace(".vcf", ".hmp.txt")).getAbsolutePath();
//            ExportUtils.writeToHapmap(tgt, outfileS);
        }
    }
    
    public void calculateHeterozygosity () {
        int windowSize = 1000;
        int step = 500;
        String infileDirS = "M:\\production\\cassava\\hmp\\lowHetsTaxa\\";
        String outfileDirS = "M:\\production\\cassava\\hmp\\lowHetsTaxa\\heterozygosity\\";
        File[] fs = IoUtils.listFilesEndsWith(new File(infileDirS).listFiles(), "txt");
        for (int i = 0; i < fs.length; i++) {
            String outfileS = new File(outfileDirS, fs[i].getName().replace(".hmp.", ".het.")).getAbsolutePath();
            GenotypeTable gt = ImportUtils.readFromHapmap(fs[i].getAbsolutePath());
            TDoubleArrayList[] hetList = new TDoubleArrayList[gt.numberOfTaxa()];
            ArrayList<String> chrList = new ArrayList();
            TIntArrayList startList = new TIntArrayList();
            TIntArrayList endList = new TIntArrayList();
            for (int j = 0; j < hetList.length; j++) {
                hetList[j] = new TDoubleArrayList();
            }
            int cnt = 0;
            for (int j = 0; j < gt.numberOfSites(); j += step) {
                int end = j+windowSize;
                if (end>gt.numberOfSites()) end = gt.numberOfSites();
                for (int k = 0; k < gt.numberOfTaxa(); k++) {
                    double v = 0;
                    for (int u = j; u < end; u++) {
                        if (gt.isHeterozygous(k, u)) v++;
                    }
                    v = v/(end-j);
                    hetList[k].add(v);
                }
                chrList.add(gt.chromosome(j).getName());
                startList.add(gt.physicalPositions()[j]);
                endList.add(gt.physicalPositions()[end-1]);
            }
            try {
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                bw.write("Region\tStart\tEnd");
                for (int j = 0 ; j < gt.numberOfTaxa(); j++) {
                    bw.write("\t"+gt.taxaName(j));
                }
                bw.newLine();
                String[] chr = chrList.toArray(new String[chrList.size()]);
                int[] start = startList.toArray();
                int[] end = endList.toArray();
                double[][] hets = new double[hetList.length][];
                for (int j = 0; j < hets.length; j++) hets[j] = hetList[j].toArray();
                for (int j = 0; j < chr.length; j++) {
                    bw.write(chr[j]+"\t"+String.valueOf(start[j])+"\t"+String.valueOf(end[j]));
                    for (int k = 0; k < hets.length; k++) {
                        bw.write("\t"+hets[k][j]);
                    }
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
    
    
}
