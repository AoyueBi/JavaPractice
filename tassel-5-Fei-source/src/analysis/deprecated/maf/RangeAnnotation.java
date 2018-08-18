/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.deprecated.maf;

import format.Bins;
import format.GeneFeature;
import format.Range;
import format.RangeAttribute;
import format.Ranges;
import format.Table;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.THashSet;
import gnu.trove.set.hash.TIntHashSet;
import graphcis.r.Histogram;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import utils.IOFileFormat;
import utils.IoUtils;

/**
 *
 * @author Fei Lu
 */
class RangeAnnotation {
    
    public RangeAnnotation () {
        this.convertGFF();
        //this.mkGeneFeatureRange();
        //this.mkCentromereDistance();
        //this.mkRepeatRange();
        //this.mkMNaseRange();
        //this.mkRecombinationRange();
        //this.mkMethylationRange();
    }
    
    public void mkMethylationRange () {
        String infileS = "Q:\\Zea\\Genotypes\\Annotations\\methylation\\GSE39232_B73_sequence.rmdup.v3.mvavg.window1000.step10.dat.gz";
        String cpgFileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\CpG.methylation.ra.txt";
        String chgFileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\CHG.methylation.ra.txt";
        String chhFileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\CHH.methylation.ra.txt";
        ArrayList<Range> rList = new ArrayList();
        TFloatArrayList cpgVList = new TFloatArrayList();
        TFloatArrayList chgVList = new TFloatArrayList();
        TFloatArrayList chhVList = new TFloatArrayList();
        try {
            BufferedReader br = IoUtils.getTextGzipReader(infileS);
            String temp = br.readLine();
            String[] tem;
            while ((temp = br.readLine()) != null) {
                tem = temp.split("\t");
                int chr = Integer.valueOf(tem[0]);
                int start = Integer.valueOf(tem[1]);
                int end = Integer.valueOf(tem[2]);
                rList.add(new Range(chr, start, end));
                cpgVList.add(Float.valueOf(tem[3]));
                chgVList.add(Float.valueOf(tem[4]));
                chhVList.add(Float.valueOf(tem[5]));
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        Range[] rs = rList.toArray(new Range[rList.size()]);
        byte[] strand = new byte[rs.length];
        for (int i = 0; i < strand.length; i++) {
            strand[i] = Byte.MIN_VALUE;
        }
        RangeAttribute cpg = new RangeAttribute(rs, "CpG", strand, cpgVList.toArray());
        RangeAttribute chg = new RangeAttribute(rs, "CHG", strand, chgVList.toArray());
        RangeAttribute chh = new RangeAttribute(rs, "CHH", strand, chhVList.toArray());
        cpg.addNote("stepSize=10;windowSize=1000");
        chg.addNote("stepSize=10;windowSize=1000");
        chh.addNote("stepSize=10;windowSize=1000");
        cpg.writeFile(cpgFileS, IOFileFormat.Text);
        chg.writeFile(chgFileS, IOFileFormat.Text);
        chh.writeFile(chhFileS, IOFileFormat.Text);
    }
    
    public void mkRecombinationRange () {
        String infileS = "E:\\Database\\maize\\recombination\\crossoverInterval.nam.cnnam.agpv3.txt";
        String infoFileS = "E:\\Database\\InfoFile\\ChrLenCentPosi_agpV3.txt";
        String outfileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\recombination.ra.txt";
        Table t = new Table (infoFileS);
        int[] chrLength = new int[t.getRowNumber()];
        for (int i = 0; i < chrLength.length; i++) {
            chrLength[i] = Integer.valueOf(t.content[i][1]);
        }
        t = new Table (infileS);
        RangeAttribute[] ras = new RangeAttribute[chrLength.length];
        int step = 10000;
        int window = 500000;
        for (int i = 0; i < chrLength.length; i++) {
            TIntArrayList posList = new TIntArrayList();
            TDoubleArrayList valueList = new TDoubleArrayList();
            for (int j = 0; j < t.getRowNumber(); j++) {
                int chrIndex = Integer.valueOf(t.content[j][0]) - 1;
                if (chrIndex!=i) continue;
                int pos = (Integer.valueOf(t.content[j][4]) + Integer.valueOf(t.content[j][5]))/2;
                posList.add(pos);
                valueList.add(1);
            }
            int[] position = posList.toArray();
            Arrays.sort(position);
            double[] value = valueList.toArray();
            ras[i] = this.getStepAvgValueFromWindow(i+1, chrLength[i], step, window, position, value, "Recombination rate");
        }
        RangeAttribute ra = ras[0];
        for (int i = 1; i < ras.length; i++) ra = ra.merge(ras[i]);
        ra.addNote("stepSize="+String.valueOf(step)+";windowSize="+String.valueOf(window));
        ra.writeFile(outfileS, IOFileFormat.Text);
    }
    
    /**
     * Return the average value of a step (value/bp) from a window
     * @param chromLength
     * @param stepSize
     * @param windowSize
     * @param position need to be sorted
     * @param value
     * @return 
     */
    private RangeAttribute getStepAvgValueFromWindow (int chromosome, int chromLength, int stepSize, int windowSize, int[] position, double[] value, String annotation) {
        int half = windowSize/2;
        Bins b = new Bins(1, chromLength, stepSize);
        TFloatArrayList vList = new TFloatArrayList();
        ArrayList<Range> rList = new ArrayList();
        for (int i = 0; i < b.getBinNum(); i++) {
            int start = b.getBinStart(i);
            int end = b.getBinEnd(i);
            int mid = (start+end)/2;
            start = mid-half;
            end = mid + half;
            if (start < 1) {
                end = end-start;
                start = 1;
            }
            if (end > chromLength+1) {
                start = start - (end-chromLength);
                end = chromLength+1;
            }
            int startIndex = Arrays.binarySearch(position, start);
            if (startIndex < 0) startIndex = -startIndex-1;
            int endIndex = Arrays.binarySearch(position, end);
            if (endIndex < 0) endIndex = -endIndex-1;
            double sum = 0;
            for (int j = startIndex; j < endIndex; j++) {
                sum+=value[j];
            }
            vList.add((float)(sum/windowSize));
            rList.add(new Range(chromosome, b.getBinStart(i), b.getBinEnd(i)));
        }
        Range[] r = rList.toArray(new Range[rList.size()]);
        float[] values = vList.toArray();
        byte[] strand = new byte[r.length];
        for (int i = 0; i < strand.length; i++) strand[i] = Byte.MIN_VALUE;
        RangeAttribute ra = new RangeAttribute(r, annotation, strand, values);
        ra.sortByStartPosition();
        return ra;
    }
    
    public void mkMNaseRange () {
        String infileS = "Q:\\Zea\\Genotypes\\Annotations\\HankBass_MNase\\AGPv3_mapping\\bowtie_output\\RP.bfthresh1.1.MNaseHS.Ranges.dat";
        String identifier = "MNase_root";
        String outfileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\MNase_root.ra.txt";
        this.writeMNaseRange(infileS, outfileS, identifier);
        infileS = "Q:\\Zea\\Genotypes\\Annotations\\HankBass_MNase\\AGPv3_mapping\\bowtie_output\\AP.bfthresh1.1.MNaseHS.Ranges.dat";
        identifier = "MNase_shoot";
        outfileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\MNase_shoot.ra.txt";
        this.writeMNaseRange(infileS, outfileS, identifier);
    }
    
    private void writeMNaseRange (String infileS, String outfileS, String identifier) {
        ArrayList<Range> rList = new ArrayList();
        try {
            BufferedReader br = IoUtils.getTextReader(infileS);
            String temp;
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                rList.add(new Range(Integer.valueOf(tem[0]), Integer.valueOf(tem[1]), Integer.valueOf(tem[2])));
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        Ranges r = new Ranges(rList.toArray(new Range[rList.size()]), identifier);
        r.sortByStartPosition();
        byte[] strand = new byte[r.getRangeNumber()];
        float[] value = new float[r.getRangeNumber()];
        for (int i = 0; i < strand.length; i++) {
            strand[i] = Byte.MIN_VALUE;
            value[i] = Float.NaN;
        }
        RangeAttribute ra = new RangeAttribute(r, strand, value);
        ra.writeFile(outfileS, IOFileFormat.Text);
    }
    
    public void mkRepeatRange () {
        String[] allRepeatType = {"Dust", "LTRs", "Low complexity regions", "Other repeats", "RNA repeats", "Satellite repeats", 
            "Tandem repeats", "Type I Transposons/LINE", "Type I Transposons/SINE", "Type II Transposons", "Unknown"};
        String[] allRepeatClass = {"DNA", "DNA/En-Spm", "DNA/Harbinger", "DNA/MuDR", "DNA/TcMar", "DNA/TcMar-Pogo", "DNA/TcMar-Tc1", 
            "DNA/hAT", "DNA/hAT-Ac", "DNA/hAT-Tip100", "LINE", "LTR", "LTR/", "LTR/Copia", "LTR/Gypsy", "LTR/abiri", "LTR/afad", "LTR/afeke", 
            "LTR/afuv", "LTR/agep", "LTR/ahoru", "LTR/ahov", "LTR/ajajog", "LTR/ajipe", "LTR/ajit", "LTR/alaw", "LTR/amiin", "LTR/anar", 
            "LTR/aneas", "LTR/anim", "LTR/ansuya", "LTR/anysaf", "LTR/apil", "LTR/arar", "LTR/atej", "LTR/atop", "LTR/avahi", "LTR/awuhe", 
            "LTR/baha", "LTR/baso", "LTR/bavav", "LTR/bawigu", "LTR/beboso", "LTR/beby", "LTR/bene", "LTR/beva", "LTR/bihar", "LTR/bobeg", 
            "LTR/bobobo", "LTR/bogu", "LTR/boha", "LTR/boja", "LTR/bomevy", "LTR/bori", "LTR/bosohe", "LTR/bovo", "LTR/bowuow", "LTR/bs1", "LTR/buire", 
            "LTR/bula", "LTR/bumy", "LTR/bygum", "LTR/cinful-zeon", "LTR/crm", "LTR/dabe", "LTR/dadeir", "LTR/dady", "LTR/dagaf", "LTR/daju", "LTR/dala", "LTR/debeh", 
            "LTR/defub", "LTR/demo", "LTR/depuw", "LTR/dijap", "LTR/doba", "LTR/doke", "LTR/dolovu", "LTR/donuil", "LTR/dugiab", "LTR/ebel", "LTR/ehahu", "LTR/ekoj", "LTR/elalal", 
            "LTR/emuh", "LTR/eninu", "LTR/epiil", "LTR/epohi", "LTR/epom", "LTR/etug", "LTR/eugene", "LTR/ewib", "LTR/ewigyw", "LTR/ewiut", "LTR/ewog", "LTR/ewot", "LTR/fajy", 
            "LTR/fanuab", "LTR/fara", "LTR/fate", "LTR/fege", "LTR/fehod", "LTR/finaij", "LTR/fipi", "LTR/flip", "LTR/fosu", "LTR/fourf", "LTR/fuved", "LTR/fuvej", "LTR/gate", "LTR/gati", 
            "LTR/gekog", "LTR/giepum", "LTR/gilovu", "LTR/giream", "LTR/gofi", "LTR/grande", "LTR/guafa", "LTR/gudyeg", "LTR/gufa", "LTR/guhis", "LTR/guvi", "LTR/guwiot", 
            "LTR/guwo", "LTR/gylu", "LTR/gyma", "LTR/gyte", "LTR/habu", "LTR/hago", "LTR/halo", "LTR/hani", "LTR/hera", "LTR/hesa", "LTR/hiimam", "LTR/hiri", "LTR/hoda", 
            "LTR/homy", "LTR/hooni", "LTR/hopscotch", "LTR/huck", "LTR/hute", "LTR/huti", "LTR/hutu", "LTR/ibulaf", "LTR/ifab", "LTR/ijaat", "LTR/ijiret", "LTR/ikal", 
            "LTR/ilofaw", "LTR/ilyl", "LTR/ipiki", "LTR/iseb", "LTR/iwik", "LTR/iwim", "LTR/jakek", "LTR/janoov", "LTR/japov", "LTR/jaws", "LTR/jelat", "LTR/jeli", "LTR/ji", 
            "LTR/joemon", "LTR/juta", "LTR/kahoba", "LTR/kaise", "LTR/kake", "LTR/kase", "LTR/kawivo", "LTR/kinosi", "LTR/kise", "LTR/kubi", "LTR/kupu", "LTR/kuvi", "LTR/labe", 
            "LTR/labu", "LTR/lafa", "LTR/laiwa", "LTR/lamyab", "LTR/lata", "LTR/lenu", "LTR/leso", "LTR/lise", "LTR/loba", "LTR/loukuv", "LTR/lowy", "LTR/lusi", "LTR/lute", 
            "LTR/luteja", "LTR/lyna", "LTR/lyruom", "LTR/lywy", "LTR/machiavelli", "LTR/mafigi", "LTR/mafogo", "LTR/magellan", "LTR/mako", "LTR/maono", "LTR/maro", "LTR/mauky", 
            "LTR/mewu", "LTR/mibaab", "LTR/mijuw", "LTR/milt", "LTR/miva", "LTR/moorud", "LTR/mopin", "LTR/muekeh", "LTR/mufeub", "LTR/mulaf", "LTR/muusi", "LTR/mywur", 
            "LTR/naadira", "LTR/naasuj", "LTR/nabu", "LTR/naiba", "LTR/naijaj", "LTR/nakovu", "LTR/nakuuv", "LTR/name", "LTR/nana", "LTR/naseup", "LTR/nasi", "LTR/neafu", 
            "LTR/neha", "LTR/nene", "LTR/neteut", "LTR/niki", "LTR/nisow", "LTR/nitat", "LTR/niypo", "LTR/nobe", "LTR/nopip", "LTR/notu", "LTR/nowuv", "LTR/nuhan", "LTR/nyjuvy", 
            "LTR/odip", "LTR/ogiv", "LTR/oguod", "LTR/ojah", "LTR/ojam", "LTR/ojav", "LTR/ojokat", "LTR/okoj", "LTR/okopam", "LTR/okur", "LTR/olepo", "LTR/omoha", "LTR/omud", 
            "LTR/onal", "LTR/onub", "LTR/opie", "LTR/osed", "LTR/ovamef", "LTR/oveah", "LTR/ovev", "LTR/oviis", "LTR/ovikoh", "LTR/oweiw", "LTR/owiit", "LTR/owume", "LTR/pagof", 
            "LTR/panen", "LTR/pebi", "LTR/petopi", "LTR/pibo", "LTR/pifo", "LTR/piube", "LTR/poarow", "LTR/pope", "LTR/prem1", "LTR/puck", "LTR/pute", "LTR/raga", "LTR/raider", 
            "LTR/reina", "LTR/rely", "LTR/riiryl", "LTR/rijuep", "LTR/rimaar", "LTR/rowi", "LTR/ruda", "LTR/rufefu", "LTR/ruhi", "LTR/rulo", "LTR/ruugu", "LTR/ruwi", 
            "LTR/saahol", "LTR/sari", "LTR/satulo", "LTR/sawujo", "LTR/sehoad", "LTR/seko", "LTR/sela", "LTR/seufyt", "LTR/seuwe", "LTR/sido", "LTR/small", "LTR/soefes", 
            "LTR/sofi", "LTR/soger", "LTR/sokiit", "LTR/sowu", "LTR/stonor", "LTR/suda", "LTR/sywu", "LTR/taname", "LTR/taro", "LTR/tata", "LTR/tatu", "LTR/tekay", "LTR/teki",
            "LTR/teuta", "LTR/tisy", "LTR/tituer", "LTR/tiwe", "LTR/tojena", "LTR/toro", "LTR/totu", "LTR/tufe", "LTR/tuku", "LTR/tuteh", "LTR/tywo", "LTR/ubat", "LTR/ubel", 
            "LTR/ubep", "LTR/ubid", "LTR/ubow", "LTR/udav", "LTR/udokup", "LTR/ufonah", "LTR/ugano", "LTR/ugog", "LTR/ugymos", "LTR/uhun", "LTR/ujinas", "LTR/ukov", "LTR/ulik", 
            "LTR/uloh", "LTR/ulon", "LTR/uluil", "LTR/ulyg", "LTR/umojev", "LTR/uper", "LTR/upus", "LTR/ures", "LTR/urogor", "LTR/urum", "LTR/usif", "LTR/usuf", "LTR/utar", 
            "LTR/uvet", "LTR/uvis", "LTR/uwaf", "LTR/uwew", "LTR/uwub", "LTR/uwum", "LTR/uwuw", "LTR/vafim", "LTR/vaofen", "LTR/vedi", "LTR/vegu", "LTR/victim", "LTR/vodida", 
            "LTR/volo", "LTR/vora", "LTR/votaed", "LTR/vufe", "LTR/vufi", "LTR/vuijon", "LTR/vuna", "LTR/vusu", "LTR/waepo", "LTR/wamenu", "LTR/waneer", "LTR/wawo", "LTR/wawu",
            "LTR/weaniv", "LTR/weki", "LTR/wemu", "LTR/wihov", "LTR/wiolus", "LTR/wiru", "LTR/witi", "LTR/wiwa", "LTR/wugaab", "LTR/wuwe", "LTR/wuywu", "LTR/wyly", "LTR/xilon-diguus", 
            "LTR/ydut", "LTR/yemi", "LTR/yfages", "LTR/ypel", "LTR/yraj", "LTR/yrer", "LTR/yreud", "LTR/ytar", "LTR/ytub", "LTR/yvoj", "LTR/ywely", "LTR/ywuv", "LTR/ywyt", 
            "Low_complexity", "MobileElement", "Other", "Other/Centromeric", "Other/Simple", "RC/Helitron", "Retroelement", "SINE", "Satellite", "Unknown", "dust", "nonLTR", "rRNA", "trf"};
        String[] selectedRepeatType = {"LTRs", "Satellite repeats", 
            "Tandem repeats", "Type I Transposons/LINE", "Type I Transposons/SINE", "Type II Transposons"};
        String[] selectedRepeatClass = {"DNA", "DNA/En-Spm", "DNA/Harbinger", "DNA/MuDR", "DNA/TcMar", "DNA/TcMar-Pogo", "DNA/TcMar-Tc1", 
            "DNA/hAT", "DNA/hAT-Ac", "DNA/hAT-Tip100", "LINE", "LTR/Copia", "LTR/Gypsy", "RC/Helitron", "SINE"};
        
        String infileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gff3.gz";
        HashSet<String> infoList = new HashSet();
        try {
            BufferedReader br = IoUtils.getTextGzipReader(infileS);
            String temp = br.readLine();
            while ((temp = br.readLine()) != null) {
                infoList.add(temp);
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        String[] info = infoList.toArray(new String[infoList.size()]);
        String identifier = "type";
        String outDirS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\repeat\\type\\";
        this.writeRepeatRange(info, identifier, selectedRepeatType, outDirS);
        
        identifier = "class";
        outDirS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\repeat\\class\\";
        this.writeRepeatRange(info, identifier, selectedRepeatClass, outDirS);
        
        String outfileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\allRepeat.ra.txt";
        this.writeAllRepeatRange(info, outfileS);
    }
    
    private void writeAllRepeatRange(String[] info, String outfileS) {
        ArrayList<Range> rList = new ArrayList();
        String[] temp;
        for (int i = 0; i < info.length; i++) {
            char s = info[i].charAt(0);
            if ((int)s < 48 || (int)s > 57) continue;
            temp = info[i].split("\t");
            if (!temp[2].startsWith("repeat")) continue;
            Range r = new Range(Integer.valueOf(temp[0]),Integer.valueOf(temp[3]), Integer.valueOf(temp[4]));
            rList.add(r);
        }
        Ranges r = new Ranges(rList.toArray(new Range[rList.size()]), "allRepeat");
        r.collapse();
        byte[] strand = new byte[r.getRangeNumber()];
        float[] values = new float[r.getRangeNumber()]; 
        for (int j = 0; j < strand.length; j++) {
            strand[j] = Byte.MIN_VALUE;
            values[j] = Float.NaN ;
        }    
        RangeAttribute ra = new RangeAttribute(r, strand, values);
        ra.sortByStartPosition();
        ra.writeFile(outfileS, IOFileFormat.Text);
    }
    
    private void writeRepeatRange (String[] info, String identifier, String[] type, String outDirS) {
        String[] renamedType = new String[type.length];
        System.arraycopy(type, 0, renamedType, 0, type.length);
        for (int i = 0; i < renamedType.length; i++) {
            renamedType[i] = renamedType[i].replaceAll("\\/", "_");
        }
        new File(outDirS).mkdir();
        ArrayList<Range>[] rs = new ArrayList[type.length];
        for (int i = 0; i < rs.length; i++) rs[i] = new ArrayList();
        String[] temp;
        for (int i = 0; i < info.length; i++) {
            char s = info[i].charAt(0);
            if ((int)s < 48 || (int)s > 57) continue;
            temp = info[i].split("\t");
            if (!temp[2].startsWith("repeat")) continue;
            String[] tem = temp[8].split(";");
            for (int j = tem.length-1; j > -1; j--) {
                if (tem[j].startsWith(identifier)) {
                    String query = tem[j].split("=")[1];
                    int index = Arrays.binarySearch(type, query);
                    if (index < 0) break;
                    Range r = new Range(Integer.valueOf(temp[0]),Integer.valueOf(temp[3]), Integer.valueOf(temp[4]));
                    rs[index].add(r);
                    break;
                }
            }
        }
        for (int i = 0; i < rs.length; i++) {
            String outfileS = new File(outDirS, renamedType[i]+".ra.txt").getAbsolutePath();
            Ranges r = new Ranges(rs[i].toArray(new Range[rs[i].size()]), type[i]);
            r.collapse();
            byte[] strand = new byte[r.getRangeNumber()];
            float[] values = new float[r.getRangeNumber()]; 
            for (int j = 0; j < strand.length; j++) {
                strand[j] = Byte.MIN_VALUE;
                values[j] = Float.NaN;
            }
            RangeAttribute ra = new RangeAttribute(r, strand, values);
            ra.sortByStartPosition();
            ra.writeFile(outfileS, IOFileFormat.Text);
        }
    }
    
    public void mkCentromereDistance () {
        String infileS = "E:\\Database\\InfoFile\\ChrLenCentPosi_agpV3.txt";
        String outfileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\centromereDistance.ra.txt";
        int intervalSize = 10000;
        Table t = new Table (infileS);
        ArrayList<Range> rList = new ArrayList();
        TFloatArrayList vList = new TFloatArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            int chr = i+1;
            int cent = (Integer.valueOf(t.content[i][2])+Integer.valueOf(t.content[i][3]))/2;
            int current = 1;
            int chrLength = Integer.valueOf(t.content[i][1]);
            int firstArmLength = cent-1;
            int secondArmLength = chrLength-cent+1;
            float v;
            while (current < chrLength) {
                int next = current+intervalSize;
                if (next > chrLength) next = chrLength+1;
                int mid = (current+next)/2;
                if (mid < cent) v = (float)(cent-mid)/firstArmLength;
                else v = (float)(mid-cent+1)/secondArmLength;
                rList.add(new Range(chr, current, next));
                vList.add(v);
                current+=intervalSize;
            }
        }
        Range[] r = rList.toArray(new Range[rList.size()]);
        float[] value = vList.toArray(new float[vList.size()]);
        byte[] strand = new byte[r.length];
        for (int i = 0; i < strand.length; i++) strand[i] = Byte.MIN_VALUE;
        RangeAttribute ra = new RangeAttribute(r, "Distance to centromere", strand, value);
        ra.writeFile(outfileS, IOFileFormat.Text);
    }
    
    public void mkGeneFeatureRange () {
        String infileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
        String outfileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\transcript.ra.txt";
        GeneFeature gf = new GeneFeature(infileS);
        gf.getAllGeneTranscriptRange().writeFile(outfileS, IOFileFormat.Text);
        outfileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\geneUpstream5000.ra.txt";
        gf.getAllGeneUpstreamRange(5000).writeFile(outfileS, IOFileFormat.Text);
        outfileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\5UTR.ra.txt";
        gf.getAllGene5UTRRange().writeFile(outfileS, IOFileFormat.Text);
        outfileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\3UTR.ra.txt";
        gf.getAllGene3UTRRange().writeFile(outfileS, IOFileFormat.Text);
        outfileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\CDS.ra.txt";
        gf.getAllGeneCDSRange().writeFile(outfileS, IOFileFormat.Text);
        outfileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\intron.ra.txt";
        gf.getAllIntronRange().writeFile(outfileS, IOFileFormat.Text);
        outfileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\geneDownstream5000.ra.txt";
        gf.getAllGeneDownstreamRange(5000).writeFile(outfileS, IOFileFormat.Text);
        outfileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\intergenic.ra.txt";
        gf.getIntergeneicRange().writeFile(outfileS, IOFileFormat.Text);
    }
    
    public void convertGFF () {
        String infileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gff3.gz";
        String outfileS = "E:\\Database\\maize\\agpv3\\gene\\Zea_mays.AGPv3.26.gf.txt";
        GeneFeature gf = new GeneFeature();
        gf.readFromMaizeGFF(infileS);
        gf.writeFile(outfileS);
    }
    
    public static RangeAttribute get3UTR () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\3UTR.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute get5UTR () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\5UTR.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getCDS() {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\CDS.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getCentromereDistance () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\centromereDistance.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getCpG () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\CpG.methylation.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getCHG () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\CHG.methylation.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getCHH () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\CHH.methylation.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getGeneDownStream5000 () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\geneDownStream5000.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getGeneUpStream5000 () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\geneUpStream5000.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getIntron () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\intron.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getMNaseRoot () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\MNase_root.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getMNaseShoot () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\MNase_shoot.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getRecombination () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\recombination.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getTranscript () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\transcript.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getIntergenic () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\intergenic.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getAllRepeat () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\allRepeat.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getLTRRepeat () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\repeat\\type\\LTRs.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getLineRepeat () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\repeat\\type\\Type I Transposons_LINE.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getSineRepeat () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\repeat\\type\\Type I Transposons_SINE.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getDNARepeat () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\repeat\\type\\Type II Transposons.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getMutatorRepeat () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\repeat\\class\\DNA_MuDR.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getHelitronRepeat () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\repeat\\class\\RC_Helitron.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getCopiaRepeat () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\repeat\\class\\LTR_Copia.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
    
    public static RangeAttribute getGypsyRepeat () {
        String infileS = "M:\\Database\\maize\\annotation\\rangeAnnotation\\repeat\\class\\LTR_Gypsy.ra.txt";
        return new RangeAttribute(infileS, IOFileFormat.Text);
    }
}

