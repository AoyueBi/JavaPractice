/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.maf;

import format.Bins;
import format.Range;
import format.RangeAttribute;
import format.Ranges;
import format.Table;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import graphcis.r.ScatterPlot;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import utils.FArrayUtils;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author Fei Lu
 */
class InvariantSiteAnnotation {
    
    public InvariantSiteAnnotation () {
        //this.invariantInRange();
        //this.correlationWithAnnotation();
        //this.mkCorrelationImages();
    }
    
    public void mkCorrelationImages () {
        String infileS = "E:\\Research\\wgs_maf\\invariantSite_annotation\\invariantSiteWithAnnotation.txt";
        String imageDirS = "E:\\Research\\wgs_maf\\invariantSite_annotation\\scatterPlot\\";
        File imageDir = new File(imageDirS);
        imageDir.mkdir();
        Table t = new Table (infileS);
        double[][] value = new double[t.getColumnNumber()][t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            for (int j = 0; j < t.getColumnNumber(); j++) {
                value[j][i] = Double.valueOf(t.content[i][j]);
            }
        }
        for (int i = 1; i < value.length; i++) {
            double[][] vs = new double[2][];
            vs[0] = value[0];
            vs[1] = value[i];
            vs = FArrayUtils.removeNaN(value);
            String outfileS = new File (imageDir, t.header[i]+".pdf").getAbsolutePath();
            ScatterPlot sp = new ScatterPlot(vs[0], vs[i]);
            sp.setTitle("Annotation relationship in 100kb bin");
            sp.setXLab("Invariant site frequency");
            sp.setYLab(t.header[i]);
            sp.addTrendLine();
            sp.saveGraph(outfileS);
        }
    }
    
    public void correlationWithAnnotation () {
        int chr = 10;
        int binSize = 100000;
        int sampleInterval = 10;
        String infoFileS = "E:\\Database\\InfoFile\\ChrLenCentPosi_agpV3.txt";
        String infileDirS = "M:\\production\\maf\\wgs\\invariantSite\\";
        String outfileS = "E:\\Research\\wgs_maf\\invariantSite_annotation\\invariantSiteWithAnnotation.txt";
        String outfileS2 = "E:\\Research\\wgs_maf\\invariantSite_annotation\\invariantSiteCorrelation.txt";
        ArrayList<RangeAttribute> rList = new ArrayList();
        RangeAttribute raCHG = RangeAnnotation.getCHG();    rList.add(raCHG);
        RangeAttribute raCpG = RangeAnnotation.getCpG();    rList.add(raCpG);
        RangeAttribute raCHH = RangeAnnotation.getCHH();    rList.add(raCHH);
        RangeAttribute raCentro = RangeAnnotation.getCentromereDistance();  rList.add(raCentro);
        RangeAttribute raRecombination = RangeAnnotation.getRecombination();    rList.add(raRecombination);
        RangeAttribute[] ras = rList.toArray(new RangeAttribute[rList.size()]);
        
        Table t = new Table (infoFileS);
        int chrLength = Integer.valueOf(t.content[9][1]);
        File[] fs = new File(infileDirS).listFiles();
        AlleleDepth prime = new AlleleDepth();
        for (int i = 0; i <  fs.length; i++) {
            AlleleDepth another = new AlleleDepth(fs[i].getAbsolutePath(), IOFileFormat.Binary);
            prime = prime.merge(another);
        }
        TIntArrayList posList = new TIntArrayList();
        for (int i = 0; i < prime.getSiteNumber(); i++) {
            if (prime.getChromosome(i) != chr) continue;
            posList.add(prime.getPosition(i));
        }
        int[] invariantPos = posList.toArray();
        double[] v = new double[invariantPos.length];
        for (int i = 0; i < v.length; i++) v[i] = 1;
        Bins b = new Bins (1, chrLength,binSize, invariantPos, v);
        double[] invariantValue = new double[b.getBinNum()];
        for (int i = 0; i < b.getBinNum(); i++) {
            invariantValue[i] = (double)b.getNumValues(i)/(b.getBinEnd(i)-b.getBinStart(i));
        }
        int[] pos = new int[chrLength/sampleInterval];
        for (int i = 0; i < pos.length; i++) {
            pos[i] = i*sampleInterval+1;
        }
        double[][] annotationValues = new double[ras.length][b.getBinNum()];
        for (int i = 0; i < ras.length; i++) {
            TDoubleArrayList vList = new TDoubleArrayList();
            TIntArrayList pList = new TIntArrayList();
            for (int j = 0; j < pos.length; j++) {
                int index = ras[i].getRangeIndex(chr, pos[j]);
                if (index < 0) continue;
                pList.add(pos[j]);
                vList.add(ras[i].getValue(index));
            }
            double[] vs = vList.toArray();
            int[] ps = pList.toArray();
            b = new Bins (1, chrLength,binSize, ps, vs);
            for (int j = 0; j < b.getBinNum(); j++) {
                annotationValues[i][j] = b.getBinAverage(j);
            }
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("InvariantInBin");
            for (int i = 0; i < ras.length; i++) {
                bw.write("\t"+ras[i].getAnnotation());
            }
            bw.newLine();
            for (int i = 0; i < invariantValue.length; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(invariantValue[i]);
                for (int j = 0; j < ras.length; j++) {
                    sb.append("\t").append(annotationValues[j][i]);
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
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS2);
            bw.write("Class\tr\tr2");
            bw.newLine();
            for (int i = 0; i < ras.length; i++) {
                TDoubleArrayList v1List = new TDoubleArrayList();
                TDoubleArrayList v2List = new TDoubleArrayList();
                for (int j = 0; j < invariantValue.length; j++) {
                    if (Double.isNaN(invariantValue[j])) continue;
                    if (Double.isNaN(annotationValues[i][j])) continue;
                    v1List.add(invariantValue[j]);
                    v2List.add(annotationValues[i][j]);
                }
                double[] v1 = v1List.toArray();
                double[] v2 = v2List.toArray();
                double r = new PearsonsCorrelation().correlation(v1, v2);
                bw.write(ras[i].getAnnotation()+"\t"+String.valueOf(r)+"\t"+Math.pow(r, 2));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void invariantInRange () {
        String infileDirS = "M:\\production\\maf\\wgs\\invariantSite\\";
        String outfileS = "E:\\Research\\wgs_maf\\invariantSite_annotation\\invariantSiteByRange.txt";
        String uniqueGenomeFileS = "M:\\production\\maf\\wgs\\uniqueKmer\\position\\uniqueGenome.pos.r.txt";
        Ranges uniqueGenome = new Ranges(uniqueGenomeFileS, IOFileFormat.Text);
        long uniqueGenomeLength = uniqueGenome.getTotalRangeSize();
        ArrayList<RangeAttribute> rList = new ArrayList();
        RangeAttribute ra3UTR = RangeAnnotation.get3UTR();  rList.add(ra3UTR);
        RangeAttribute ra5UTR = RangeAnnotation.get5UTR();  rList.add(ra5UTR);
        RangeAttribute raAllRepeat = RangeAnnotation.getAllRepeat();    rList.add(raAllRepeat);
        RangeAttribute raCDS = RangeAnnotation.getCDS();    rList.add(raCDS);
        RangeAttribute raCopia = RangeAnnotation.getCopiaRepeat();  rList.add(raCopia);
        RangeAttribute raDNARepeat = RangeAnnotation.getDNARepeat(); rList.add(raDNARepeat);
        RangeAttribute raDownStream = RangeAnnotation.getGeneDownStream5000();  rList.add(raDownStream);
        RangeAttribute raUpStream = RangeAnnotation.getGeneUpStream5000();  rList.add(raUpStream);
        RangeAttribute raGypsyRepeat = RangeAnnotation.getGypsyRepeat();    rList.add(raGypsyRepeat);
        RangeAttribute raMutator = RangeAnnotation.getMutatorRepeat();    rList.add(raMutator);
        RangeAttribute raHelitron = RangeAnnotation.getHelitronRepeat();    rList.add(raHelitron);
        RangeAttribute raIntergenic = RangeAnnotation.getIntergenic();  rList.add(raIntergenic);
        RangeAttribute raIntron = RangeAnnotation.getIntron();  rList.add(raIntron);
        RangeAttribute raLTR = RangeAnnotation.getLTRRepeat();  rList.add(raLTR);
        RangeAttribute raLineRepeat = RangeAnnotation.getLineRepeat();  rList.add(raLineRepeat);
        RangeAttribute raMNaseRoot = RangeAnnotation.getMNaseRoot();    rList.add(raMNaseRoot);
        RangeAttribute raMNaseShoot = RangeAnnotation.getMNaseShoot();  rList.add(raMNaseShoot);
        RangeAttribute raSineRepeat = RangeAnnotation.getSineRepeat();  rList.add(raSineRepeat);
        RangeAttribute raTranscript = RangeAnnotation.getTranscript();  rList.add(raTranscript);
        RangeAttribute[] ras = rList.toArray(new RangeAttribute[rList.size()]);
        
        File[] fs = new File(infileDirS).listFiles();
        AlleleDepth prime = new AlleleDepth();
        for (int i = 0; i <  fs.length; i++) {
            AlleleDepth another = new AlleleDepth(fs[i].getAbsolutePath(), IOFileFormat.Binary);
            prime = prime.merge(another);
        }
        
        double[] valueOnRangeLength = new double[ras.length];
        double[] valueOnUniqueRangeLength = new double[ras.length];
        long[] ls = new long[ras.length];
        long[] uls = new long[ras.length];
        int[] cnt = new int[ras.length];
        for (int i = 0; i < ras.length; i++) {
            long totalLength = ras[i].getTotalRangeSize();
            uniqueGenome.annotation = ras[i].getAnnotation();
            Ranges nr = uniqueGenome.merge(ras[i]);
            nr.collapse();
            long nrLength = nr.getTotalRangeSize();
            long overlapLength = totalLength+uniqueGenomeLength-nrLength;
            
            int hit = 0;
            for (int j = 0; j < prime.getSiteNumber(); j++) {
                if (ras[i].isInRanges(prime.getChromosome(j), prime.getPosition(j))) {
                    hit++;
                }
                if (j%10000000 == 0) System.out.println(String.valueOf(j)+" sites scaned");
            }
            valueOnRangeLength[i] = (double)hit/totalLength;
            valueOnUniqueRangeLength[i] = (double)hit/overlapLength;
            ls[i] = totalLength;
            uls[i] = overlapLength;
            cnt[i] = hit;
            System.out.println(String.valueOf(i+1)+"th annotation\t"+ras[i].annotation+"\t:"+String.valueOf(valueOnRangeLength[i])+"\t"+String.valueOf(valueOnUniqueRangeLength[i]));
        }
        try {
            BufferedWriter bw = IoUtils.getTextWriter(outfileS);
            bw.write("Class\tInvariant site frequency\tInvariant site frequency in unique genome\tInvariant site number\tAnnotation length\tAnnotation length in unique genome");
            bw.newLine();
            for (int i = 0; i < ras.length; i++) {
                bw.write(ras[i].getAnnotation()+"\t"+String.valueOf(valueOnRangeLength[i])+"\t"+String.valueOf(valueOnUniqueRangeLength[i])+"\t"+String.valueOf(cnt[i])+"\t"+String.valueOf(ls[i])+"\t"+String.valueOf(uls[i]));
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
