package analysis.maf;

public class UniquenessScoreGo {

    public UniquenessScoreGo (String referenceGenome, String testGenome, String outputDirS) {
//        referenceGenome = "/Users/feilu/Documents/analysisL/pipelineTest/uniqueScore/maize_chr12.fa";
//        testGenome = "/Users/feilu/Documents/analysisL/pipelineTest/uniqueScore/maize_chr12.fa";
//        outputDirS = "/Users/feilu/Documents/analysisL/pipelineTest/uniqueScore/result";
        new SiteUniqueness(referenceGenome, testGenome, outputDirS);
    }

    public static void main (String[] args) {
        new UniquenessScoreGo(args[0], args[1], args[2]);
    }
}
