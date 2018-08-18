/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.elementMap;

import utils.BaseCoder;
import format.Fasta;
import format.GeneFeature;
import format.Table;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TShortByteHashMap;
import gnu.trove.map.hash.TShortIntHashMap;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import net.openhft.koloboke.collect.map.hash.HashByteByteMap;
import net.openhft.koloboke.collect.map.hash.HashShortByteMap;
import net.openhft.koloboke.collect.map.hash.HashShortByteMaps;
import utils.FStringUtils;
import utils.IoUtils;

/**
 *
 * @author fl262
 */
public class UpstreamExpression {
    
    public UpstreamExpression () {
        this.mkGeneCpScoreTable();
        //this.mkGeneCpScoreKmerTable();
        //test();
    }
    
    private void test () {
        byte[] a = {0,0,0, 0}; //a, c, g
        byte[] b = BaseCoder.getReverseComplementary(a);
        short as = BaseCoder.getShortSeqFromByteArray(a);
        short bs = BaseCoder.getShortSeqFromByteArray(b);
        System.out.println(BaseCoder.getSequenceFromShort(as));
        System.out.println(BaseCoder.getSequenceFromShort(bs));
    }
    
    private void mkGeneCpScoreKmerTable () {
        String geneCpScoreDirS = "M:\\production\\elementMap\\upstreamExpression\\geneCpScore\\";
        String genomeFileS = "M:\\Database\\maize\\agpV3\\maize.agpV3.fa";
        String geneCpScoreKmerDirS = "M:\\production\\elementMap\\upstreamExpression\\geneKmer\\";
        int kmerLength = 8;
        int upstreamLength = 500;
        double scoreHighCut = 1;
        double scoreLowCut = 0;
        new File(geneCpScoreKmerDirS).mkdir();
        short[] kmers = this.getKmerSetFromKmerLength(kmerLength);
        String header = this.getHeader(kmers);
        byte[][] chrSeq = this.getGenomeSequenceByte(genomeFileS);
        File[] fs = new File(geneCpScoreDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            String outfileS = f.getName();
            outfileS = outfileS.replaceFirst("upstreamScore.txt.gz", "upstreamKmer.txt");
            outfileS = new File(geneCpScoreKmerDirS, outfileS).getAbsolutePath();
            try{
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                bw.write(header);
                bw.newLine();
                BufferedReader br = IoUtils.getTextGzipReader(f.getAbsolutePath());
                String temp = br.readLine();
                List<String> l = FStringUtils.fastSplit(temp);
                int oriLength = Integer.valueOf(l.get(l.size()-1));
                while ((temp = br.readLine()) != null) {
                    TShortIntHashMap[] kmerCountMaps = new TShortIntHashMap[2];
                    StringBuilder[] sbs = new StringBuilder[2];
                    for (int i = 0; i < kmerCountMaps.length; i++) {
                        int[] nullCount = new int[kmers.length];
                        kmerCountMaps[i] = new TShortIntHashMap(kmers, nullCount);
                        sbs[i] = new StringBuilder();
                    }
                    l = FStringUtils.fastSplit(temp);
                    for (int i = 0; i < sbs.length; i++) {
                        sbs[i].append(l.get(1)).append("\t").append(l.get(2)).append("\t").append(l.get(3)).append("\t").append(l.get(4)).append("\t").append(l.get(0)).append("\t").append(i+1);
                    }
                    float[] cpScore = new float[upstreamLength];
                    for (int i = 0; i < cpScore.length; i++) {
                        cpScore[i] = Float.valueOf(l.get(i+5+oriLength-upstreamLength));
                    }
                    int chr = Integer.valueOf(l.get(2));
                    int startPos = Integer.valueOf(l.get(3));
                    int strand = Integer.valueOf(l.get(4));
                    byte[] chrSeqByte = new byte[upstreamLength];
                    if (strand == 1) {
                        System.arraycopy(chrSeq[chr-1], startPos-upstreamLength, chrSeqByte, 0, upstreamLength);
                    }
                    else {
                        System.arraycopy(chrSeq[chr-1], startPos, chrSeqByte, 0, upstreamLength);
                        chrSeqByte = BaseCoder.getReverseComplementary(chrSeqByte);
                    }
                    int start = -1;
                    int end = -1;
                    for (int i = 0; i < cpScore.length; i++) {
                        if (cpScore[i] < scoreHighCut && cpScore[i] >= scoreLowCut) {
                            end = i+1;
                            if (end == cpScore.length) {
                                int length = end-start-1;
                                if (length>kmerLength) {
                                    byte[] bArray = new byte[length] ;
                                    System.arraycopy(chrSeqByte, start+1, bArray, 0, length);
                                    for (int j = 0; j < bArray.length-kmerLength; j++) {
                                        boolean flag = false;
                                        for (int k = j; k < j+kmerLength; k++) {
                                            if (bArray[k] > 3) {
                                                j = k+1;
                                                flag = true;
                                                break;
                                            }
                                        }
                                        if (flag) continue;
                                        short kmerS = BaseCoder.getShortSeqFromByteArray(Arrays.copyOfRange(bArray, j, j + kmerLength));
                                        short kmerRS = BaseCoder.getShortReverseComplement(kmerS, kmerLength);
                                        kmerCountMaps[0].adjustValue(kmerS, 1);
                                        kmerCountMaps[1].adjustValue(kmerS, 1);
                                        kmerCountMaps[1].adjustValue(kmerRS, 1);
                                    }
                                }
                            }
                        }
                        else {
                            if (start != end) {
                                int length = end-start;
                                if (length>kmerLength) {
                                    byte[] bArray = new byte[length] ;
                                    System.arraycopy(chrSeqByte, start+1, bArray, 0, length);
                                    for (int j = 0; j < bArray.length-kmerLength; j++) {
                                        boolean flag = false;
                                        for (int k = j; k < j+kmerLength; k++) {
                                            if (bArray[k] > 3) {
                                                j = k+1;
                                                flag = true;
                                                break;
                                            }
                                        }
                                        if (flag) continue;
                                        short kmerS = BaseCoder.getShortSeqFromByteArray(Arrays.copyOfRange(bArray, j, j + kmerLength));
                                        short kmerRS = BaseCoder.getShortReverseComplement(kmerS, kmerLength);
                                        kmerCountMaps[0].adjustValue(kmerS, 1);
                                        kmerCountMaps[1].adjustValue(kmerS, 1);
                                        kmerCountMaps[1].adjustValue(kmerRS, 1);
                                    }
                                }
                            }
                            start = i;
                            end = i;
                        }
                    }
                    for (int i = 0; i < kmers.length; i++) {
                        sbs[0].append("\t").append(kmerCountMaps[0].get(kmers[i]));
                        sbs[1].append("\t").append(kmerCountMaps[1].get(kmers[i]));
                    }
                    for (int i = 0; i < sbs.length; i++) {
                        bw.write(sbs[i].toString());
                        bw.newLine();
                    }
                     
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
    
    private String getHeader(short[] kmers) {
        StringBuilder sb = new StringBuilder("Gene");
        sb.append("\tChromosome\tStartPos\tStrand\tTaxa\tCountStrand");
        for (int i = 0; i < kmers.length; i++) {
            sb.append("\t").append(BaseCoder.getSequenceFromShort(kmers[i]));
        }
        return sb.toString();
    }
    
    private byte[][] getGenomeSequenceByte (String genomeFileS) {
        HashByteByteMap ascIIByteMap = BaseCoder.getAscIIByteMap();
        Fasta f = new Fasta(genomeFileS);
        f.sortRecordByNameValue();
        byte[][] chrSeq = new byte[f.getSeqNumber()][];
        for (int i = 0; i < f.getSeqNumber(); i++) {
            String seqS = f.getSeq(i);
            byte[] seqByte = f.getSeq(i).getBytes();
            for (int j = 0; j < seqByte.length; j++) {
                seqByte[j] = ascIIByteMap.get(seqByte[j]);
            }
            chrSeq[i] = seqByte;
        }
        return chrSeq;
    }
    
    private short[] getKmerSetFromKmerLength (int kmerLength) {
        byte[] baseByte = {0,1,2,3};
        int maxNumber = (int)Math.pow(4, kmerLength);
        short[] kmers = new short[maxNumber];
        int[] is = new int[kmerLength];
        int cnt = 0;
        for (is[0] = 0; is[0] < baseByte.length; is[0]++){
            for (is[1] = 0; is[1] < baseByte.length; is[1]++){
                for (is[2] = 0; is[2] < baseByte.length; is[2]++){
                    for (is[3] = 0; is[3] < baseByte.length; is[3]++){
                        for (is[4] = 0; is[4] < baseByte.length; is[4]++){
                            for (is[5] = 0; is[5] < baseByte.length; is[5]++){
                                for (is[6] = 0; is[6] < baseByte.length; is[6]++){
                                    for (is[7] = 0; is[7] < baseByte.length; is[7]++){
                                        byte[] bArray = new byte[kmerLength];
                                        for (int i = 0; i < bArray.length; i++) {
                                            bArray[i] = (byte)is[i];
                                        }
                                        short kmer = BaseCoder.getShortSeqFromByteArray(bArray);
                                        kmers[cnt] = kmer;
                                        cnt++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        Arrays.sort(kmers);
        return kmers;
    }
    
    private void mkGeneCpScoreTable () {
//        String geneinfoFileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
//        String taxaScoreDirS = "M:\\production\\elementMap\\upstreamExpression\\taxaScore\\";
//        String depthFileS = "M:\\production\\elementMap\\taxaDepth\\taxaDepth_16.txt";
//        String geneCpScoreDirS = "M:\\production\\elementMap\\upstreamExpression\\geneCpScore\\";
        
        String geneinfoFileS = "/workdir/fl262/Zea_mays.AGPv3.26.gf.txt";
        String taxaScoreDirS = "/workdir/fl262/taxaScore_16/";
        String depthFileS = "/workdir/fl262/taxaDepth/taxaDepth_16.txt";
        String geneCpScoreDirS = "/workdir/fl262/geneCpScore";
        
        int upstreamLength = 2000;
        double minDepth = 5;
        new File(geneCpScoreDirS).mkdir();
        ArrayList<String> taxaList = new ArrayList();
        Table t = new Table (depthFileS);
        HashMap<String, Double> taxaDepthMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (Double.valueOf(t.content[i][1]) < minDepth) continue;
            taxaList.add(t.content[i][0]);
            taxaDepthMap.put(t.content[i][0], Double.valueOf(t.content[i][1]));
        }
        String[] taxa = taxaList.toArray(new String[taxaList.size()]);
        Arrays.sort(taxa);
        taxaList = new ArrayList();
        File[] fs = new File (taxaScoreDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
           if (Arrays.binarySearch(taxa, fs[i].getName()) < 0) continue;
           taxaList.add(fs[i].getName());
        }
        taxa = taxaList.toArray(new String[taxaList.size()]);
        Arrays.sort(taxa);
        HashMap<String, Integer> taxaIndexMap = new HashMap();
        for (int i = 0; i < taxa.length; i++) {
            taxaIndexMap.put(taxa[i], i);
        }
        GeneFeature gf = new GeneFeature (geneinfoFileS);
        String[] genes = new String[gf.getGeneNumber()];
        HashSet<Integer> chrSet = new HashSet();
        for (int i = 0; i < genes.length; i++) {
            genes[i] = gf.getTranscriptName(i, 0);
            chrSet.add(gf.getGeneChromosome(i));
        }
        Arrays.sort(genes);
        int chrNum = chrSet.size();
        int[][] start = new int[chrNum][]; 
        byte[][] strand = new byte[chrNum][];
        String[][] geneName = new String[chrNum][];
        String[][] scoreStr = new String[gf.getGeneNumber()][taxa.length];
        TIntArrayList[] startList = new TIntArrayList[chrNum];
        TByteArrayList[] strandList = new TByteArrayList[chrNum];
        ArrayList<String>[] geneNameList = new ArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) {
            startList[i] = new TIntArrayList();
            strandList[i] = new TByteArrayList();
            geneNameList[i] = new ArrayList();
        }
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            int index = gf.getGeneChromosome(i)-1;
            startList[index].add(gf.getTranscriptStart(i, 0));//*************************************************************************problem
            strandList[index].add(gf.getTranscriptStrand(i, 0));
            geneNameList[index].add(gf.getTranscriptName(i, 0));
        }
        for (int i = 0; i < chrNum; i++) {
            start[i] = startList[i].toArray();
            strand[i] = strandList[i].toArray();
            geneName[i] = geneNameList[i].toArray(new String[geneNameList[i].size()]);
        }
        taxaList.parallelStream().forEach(taxon -> {
            String inDirS = new File(taxaScoreDirS, taxon).getAbsolutePath();
            double taxonDepth = taxaDepthMap.get(taxon);
            for (int i = 0; i < chrNum; i++) {
                int chr = i+1;
                StringBuilder sb  = new StringBuilder();
                sb.append(taxon).append("_chr").append(FStringUtils.getNDigitNumber(3, chr)).append("_CpScore.txt.gz");
                TFloatArrayList scoreList = new TFloatArrayList();
                try {
                    BufferedReader br = IoUtils.getTextGzipReader(new File(inDirS, sb.toString()).getAbsolutePath());
                    String temp = br.readLine();
                    int cnt = 0;
                    while ((temp = br.readLine()) != null) {
                        cnt++;
                        if (cnt%100000000 == 0) System.out.println(taxon+"\t"+String.valueOf(chr)+"\t"+String.valueOf(cnt));
                        float value;
                        if (temp.startsWith("N")) {
                            value = -1;
                        }
                        else {
                            value = (float)(Double.valueOf(temp)/taxonDepth);
                        }
                        scoreList.add(value);
                    }
                    br.close();
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
                float[] scores = scoreList.toArray();
                for (int j = 0; j < start[i].length; j++) {
                    sb = new StringBuilder(taxon);
                    sb.append("\t").append(geneName[i][j]).append("\t").append(chr).append("\t").append(start[i][j]).append("\t").append(strand[i][j]);
                    if (strand[i][j] == 1) {
                        for (int k = start[i][j]-upstreamLength; k < start[i][j]; k++) {
                            sb.append("\t").append(scores[k-1]);
                        }
                    }
                    else {
                        for (int k = start[i][j]+upstreamLength; k > start[i][j]; k--) {
                            sb.append("\t").append(scores[k-1]);
                        }
                    }
                    int geneIndex = Arrays.binarySearch(genes, geneName[i][j]);
                    int taxaIndex = taxaIndexMap.get(taxon);
                    scoreStr[geneIndex][taxaIndex] = sb.toString();
                }
                System.out.println(taxon+"\tchr"+String.valueOf(chr)+" is finished");
                scores = null;
                scoreList = null;
            }
        });
        StringBuilder sb = new StringBuilder();
        sb.append("Taxa\tTranscript\tChr\tStart\tStrand");
        for (int i = 0; i < upstreamLength; i++) {
            sb.append("\t").append(i+1);
        }
        String header = sb.toString();
        for (int i = 0; i < scoreStr.length; i++) {
            String outfileS = new File (geneCpScoreDirS, genes[i]+".upstreamScore.txt").getAbsolutePath();
            try {
                BufferedWriter bw = IoUtils.getTextWriter(outfileS);
                bw.write(header);
                bw.newLine();
                for (int j = 0; j < scoreStr[i].length; j++) {
                    if (scoreStr[i][j] == null) continue; 
                    bw.write(scoreStr[i][j]);
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
