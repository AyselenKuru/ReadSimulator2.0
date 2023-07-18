package Assignment2;

import Assignment1.Gene;
import Assignment1.GtfRow;
import Assignment1.Transcript;
import Utils.FileUtils;
import org.apache.commons.math3.distribution.NormalDistribution;




import java.io.File;

import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class ReadSimulator {

    public static void main(String[] args) {

        //parameter treatment
        if (args.length == 0) {
            System.out.println("execute the program with the following parameters:");
            System.out.println("-length <int> : read length");
            System.out.println("-frlength <int>  : fragment length distribution");
            System.out.println("- SD <int> : fragment length distribution");
            System.out.println("-readcounts <file> : table of gene_id, transcript_id, count tuples");
            System.out.println("mutationsrate <int> : mutation rate in percent (between 0.0 and 100.0");
            System.out.println("-fasta <FASTA file> : genome fasta");
            System.out.println("-fidx <FASTA file index> : index for genome fasta");
            System.out.println("-gtf <GTF file> : genomic annotation");
            System.out.println("-od <output file path>: path to the output directory");
            System.exit(1);
        }

        int readLength = 0;
        int mean = 0;
        int standardDeviation = 0;
        String readcountsFilename = null;
        double mutationsrate = 0.0;
        String fastaFilename = null;
        String fastaidxFilename = null;
        String gtfFilename = null;
        String outputpath = null;

        for (int i = 0; i < args.length - 1; i++) {

            if (args[i].equals("-length")) {
                if (!args[i + 1].startsWith("-")) {
                    readLength = Integer.parseInt(args[i + 1]);
                }
            } else if (args[i].equals("-frlength")) {
                if (!args[i + 1].startsWith("-")) {
                    mean = Integer.parseInt(args[i + 1]);
                }
            } else if (args[i].equals("-SD")) {
                if (!args[i + 1].startsWith("-")) {
                    standardDeviation = Integer.parseInt(args[i + 1]);
                }
            } else if (args[i].equals("-readcounts")) {
                if (!args[i + 1].startsWith("-")) {
                    readcountsFilename = args[i + 1];
                }
            } else if (args[i].equals("-mutationrate")) {
                if (!args[i + 1].startsWith("-")) {
                    mutationsrate = Double.parseDouble(args[i + 1]);
                }
            } else if (args[i].equals("-fasta")) {
                if (!args[i + 1].startsWith("-")) {
                    fastaFilename = args[i + 1];
                }
            } else if (args[i].equals("-fidx")) {
                if (!args[i + 1].startsWith("-")) {
                    fastaidxFilename = args[i + 1];
                }
            } else if (args[i].equals("-gtf")) {
                if (!args[i + 1].startsWith("-")) {
                    gtfFilename = args[i + 1];
                }
            } else if (args[i].equals("-od")) {
                if (!args[i + 1].startsWith("-")) {
                    outputpath = args[i + 1];
                }
            }
        }
        //notifications
        if (readcountsFilename == null) {
            System.out.println("readcounts file is missing");
            System.exit(1);
        } else if (fastaidxFilename == null) {
            System.out.println("fasta index file path is missing");
            System.exit(1);
        } else if (gtfFilename == null) {
            System.out.println("gtf file path is missing");
            System.exit(1);
        } else if (fastaFilename == null) {
            System.out.println("fasta  file path is missing");
            System.exit(1);
        } else if (readLength == 0) {
            System.out.println("readLength is missing");
            System.exit(1);
        } else if (mean == 0) {
            System.out.println("mean is missing");
            System.exit(1);
        } else if (standardDeviation == 0) {
            System.out.println("standard deviation is missing");
            System.exit(1);
        } else if (mutationsrate == 0) {
            System.out.println("mutationrate is missing");
            System.exit(1);
        } else if (outputpath == null) {
            System.out.println("outputhpath is missing");
            System.exit(1);
        }

        //Dateien einlesen

        //readcounts file for transcripts who should be simulated
        ArrayList<String[]> readcounts;
        readcounts = FileUtils.readCounts(readcountsFilename);

        HashMap<String, Gene> gtfMap;
        gtfMap = FileUtils.readGTFA3(gtfFilename, readcounts); //reads rows with gene and exon as feature


        try {
            GenomeSequenceExtractor extractor = new GenomeSequenceExtractor(new File(fastaFilename), new File(fastaidxFilename));

            String geneSequence = "";
            int geneStart = 0;

            int rmember = 0;
            HashMap<String, String> geneMap = new HashMap<>();

            FileWriter writefileFW = new FileWriter(outputpath + "/fw.fastq");
            FileWriter writefileRW = new FileWriter(outputpath + "/rw.fastq");
            FileWriter writefileRM = new FileWriter(outputpath + "/read.mappinginfo");
            //FileWriter writefilePlotFL = new FileWriter(outputpath + "/fragmentlengthDistribution");
            //FileWriter writefileMD = new FileWriter(outputpath + "/mutationDistribution");


            writefileRM.write("readid\tchr\tgene\ttranscript\tt_fw_regvec\tt_rw_regvec\tfw_regvec\trw_regvec\tfw_mut\trw_mut" + "\n");


            for (int i = 0; i < readcounts.size(); i++) {

                String geneID = readcounts.get(i)[0];
                String transcriptID = readcounts.get(i)[1];
                StringBuilder transcirptSequenceBuilder = new StringBuilder();
                String chr = "";

                ArrayList<Integer[]> exonsInTrans = new ArrayList<>();

                Gene g = gtfMap.get(geneID);
                //if (g.data.strand.equals("-")) {
                // continue;
                //}
                if (!geneMap.containsKey(geneID)) {
                    geneSequence = extractor.getSequence(g.data.seqname, g.data.start, g.data.end); //get gene Sequence
                    geneMap.put(geneID, geneSequence);
                } else {
                    geneSequence = geneMap.get(geneID);
                }
                chr = g.data.seqname;
                geneStart = g.data.start;
                Transcript t = g.transcriptMap.get(transcriptID);
                String transcriptSeq = "";

                for (GtfRow exon : t.exons) {
                    String exonSeq = geneSequence.substring(exon.start - geneStart, exon.end - geneStart + 1);
                    Integer[] eChrPos = new Integer[2];
                    eChrPos[0] = exon.start;
                    eChrPos[1] = exon.end;
                    exonsInTrans.add(eChrPos);
                    chr = exon.seqname;
                    transcirptSequenceBuilder.append(exonSeq);
                }
                transcriptSeq = transcirptSequenceBuilder.toString();
                /*if (g.data.strand.equals("-")) {
                    transcriptSeq = reverseCompliment(transcriptSeq);
                }*/
                //transcript bearbeitung
                if(transcriptID.equals("ENST00000581687")){
                    transcriptSeq = transcirptSequenceBuilder.toString();
                }

                int readcount = Integer.parseInt(readcounts.get(i)[2]);
                //int rmember =0;

                for (int r = 0; r < readcount; r++) {
                    int transcriptLen = transcirptSequenceBuilder.length();
                    int fragmentLen = getSampleFragLen(readLength, mean, standardDeviation, transcriptLen);
                    int pos = selectRandomPos(transcriptLen, fragmentLen);
                    String fragmentSeq = getFragmentSeq(transcriptSeq, pos, (pos + fragmentLen));
                    if(g.data.strand.equals("-")){
                        fragmentSeq = reverseCompliment(fragmentSeq);
                    }
                    String[] seq = getReadSeq(fragmentSeq, readLength);
                    //String fw = reverseCompliment(seq[0]);
                    //String rw = reverseCompliment(seq[1]);

                    //mutationen
                    SequenceError fw = simulateMutation(seq[0], mutationsrate);
                    SequenceError rw = simulateMutation(seq[1], mutationsrate);

                    FileUtils.writeFastq(rmember, fw.seq, writefileFW);
                    FileUtils.writeFastq(rmember, rw.seq, writefileRW);

                    String fwVectors = getVector(exonsInTrans, pos, readLength);    //mit 810 getestet

                    int rwPos = (pos + fragmentLen) - readLength;
                    String rwVectors = getVector(exonsInTrans, rwPos, readLength); //mit 858 getestet

                    if (g.data.strand.equals("-")){
                        int rwStart = transcriptLen-pos-readLength;

                        FileUtils.writeReadMapping(getReadMappingInformation(g.data.strand,rmember, chr, geneID, transcriptID, rwStart, readLength,
                                fragmentLen, rwVectors,fwVectors, fw.positions, rw.positions), writefileRM);
                    }else {
                        FileUtils.writeReadMapping(getReadMappingInformation(g.data.strand,rmember, chr, geneID, transcriptID, pos, readLength,
                                fragmentLen, fwVectors, rwVectors, fw.positions, rw.positions), writefileRM);
                    }
                    //FileUtils.writeFragmentlength(rmember,fragmentLen,writefilePlotFL);
                    //FileUtils.writeMutationDis(rmember, fw.positions, rw.positions, writefileMD);
                    rmember += 1;

                }
            }
            writefileFW.close();
            writefileRW.close();
            writefileRM.close();
            //writefilePlotFL.close();
            //writefileMD.close();

        } catch (IOException e) {
            System.out.println("extractions failed: " + e.getMessage());
        }
    }

    private static int getSampleFragLen(int readLen, int mean, int sd, int transcriptLen) {

        NormalDistribution nd = new NormalDistribution(mean, sd);
        int value = 0;
        while (value < readLen || value > transcriptLen) {
            value = (int) Math.round(nd.sample());
        }
        return value;
    }

    private static int selectRandomPos(int transcriptLen, int fragmentLen) {
        double value = Math.random() * (transcriptLen - fragmentLen);
        return (int) Math.floor(value);
    }

    private static String getFragmentSeq(String transcriptseq, int start, int end) {
        String fragmentSeq = transcriptseq.substring(start, end);
        return fragmentSeq;
    }

    private static String[] getReadSeq(String fragmentSeq, int readLen) {

        String forwardSeq = fragmentSeq.substring(0, readLen); //fw
        StringBuilder reversSeq = new StringBuilder();
        for (int i = fragmentSeq.length() - 1; i > fragmentSeq.length() - readLen - 1; i--) {
            char nuk = fragmentSeq.charAt(i);
            if (nuk == 'A') {
                reversSeq.append('T');
            }
            if (nuk == 'T') {
                reversSeq.append('A');
            }
            if (nuk == 'G') {
                reversSeq.append('C');
            }
            if (nuk == 'C') {
                reversSeq.append('G');
            }
        }
        return new String[]{forwardSeq, reversSeq.toString()};
    }

    private static SequenceError simulateMutation(String transcriptSeq, double mutationrate) {

        int transciptLen = transcriptSeq.length();//fragmente len eigentlich

        StringBuilder mutatedSeq = new StringBuilder();
        List<Integer> mutationPositions = new ArrayList<>();

        List<Character> aList = Arrays.asList('T', 'C', 'G');
        List<Character> tList = Arrays.asList('A', 'C', 'G');
        List<Character> cList = Arrays.asList('A', 'T', 'G');
        List<Character> gList = Arrays.asList('A', 'C', 'T');

        Random random = new Random();
        int randomChar = 0;

        for (int i = 0; i < transciptLen; i++) {

            if (random.nextDouble() * 100 < mutationrate) {
                mutationPositions.add(i);
                if (transcriptSeq.charAt(i) == 'A') {
                    randomChar = random.nextInt(aList.size());
                    mutatedSeq.append(aList.get(randomChar));
                } else if (transcriptSeq.charAt(i) == 'T') {
                    randomChar = random.nextInt(tList.size());
                    mutatedSeq.append(tList.get(randomChar));
                } else if (transcriptSeq.charAt(i) == 'C') {
                    randomChar = random.nextInt(cList.size());
                    mutatedSeq.append(cList.get(randomChar));
                } else if (transcriptSeq.charAt(i) == 'G') {
                    randomChar = random.nextInt(gList.size());
                    mutatedSeq.append(gList.get(randomChar));
                } else {
                    mutatedSeq.append(transcriptSeq.charAt(i)); //if there is any unknown letter
                }
            } else {
                mutatedSeq.append(transcriptSeq.charAt(i));
            }
        }

        SequenceError mutation = new SequenceError();
        mutation.seq = mutatedSeq.toString();
        mutation.positions = mutationPositions;

        return mutation;

    }

    public static String getReadMappingInformation(String strand,int readnum, String chr, String geneID, String transcriptID, int start,
                                                   int readlength, int fragmentlength, String fwVector, String rwVector,
                                                   List<Integer> fwMpos, List<Integer> rwMpos) { //start == pos, start from fragment/read
        int end = start + fragmentlength;
        String fw = start + "-" + (start + readlength);
        String rw = (end - readlength) + "-" + end;

        if (strand.equals("-")){
            rw = start + "-" + (start + readlength);

            int rwEnd = (start + readlength);
            int fwStart = rwEnd -fragmentlength;
            int fwEnd = fwStart +readlength;

            fw = fwStart+ "-" +fwEnd;
        }


        String mappingLine = readnum + "\t" + chr + "\t" + geneID + "\t" + transcriptID + "\t" + fw + "\t" + rw + "\t" + fwVector + "\t"
                + rwVector + "\t" + mutPos(fwMpos) + "\t" + mutPos(rwMpos);

        return mappingLine;
    }

    public static String mutPos(List<Integer> pos) {
        StringBuilder builder = new StringBuilder();
        for (Integer e : pos) {
            if (!builder.isEmpty()) {
                builder.append(",");
            }
            builder.append(e);
        }
        return builder.toString();
    }

    public static String getVector(ArrayList<Integer[]> exonsInTrans, int fragmentStart, int readlength) {

        int startFRest = fragmentStart;
        ArrayList<Integer[]> result = new ArrayList<>();
        boolean foundstart = false;
        int startPos = 0;

        for (Integer[] startEnd : exonsInTrans) {
            Integer[] startEndresult = new Integer[2];

            int startE = startEnd[0];
            int endE = startEnd[1];
            int lenE = endE - startE + 1;

            if (startFRest >= lenE && !foundstart) { //restlen als erstews so groß wie fragmentstart
                startFRest = startFRest - lenE;
                continue;
            } else {
                if (!foundstart) {  //gets startposition
                    foundstart = true;
                    startPos = (startE + startFRest);

                    startEndresult[0] = startPos; // start vom fragment in einem exon gefunden

                    if (startPos + readlength > endE) { //if region bigger than exon

                        startEndresult[1] = (endE + 1);
                        readlength = readlength - ((endE + 1) - startPos);
                    } else {    //if region is in one exon
                        startEndresult[1] = startPos + readlength;
                        readlength = 0;
                    }
                } else if (readlength >= lenE) {    //full exon is part of fragment
                    startEndresult[0] = startE;
                    startEndresult[1] = endE + 1;
                    readlength = readlength - lenE;

                } else {    //fragment ends in this exon

                    startEndresult[0] = startE;
                    startEndresult[1] = startE + readlength;
                    readlength = 0;

                }
            }
            //only if exon is part of fragment
            if (startEndresult[0] != null && startEndresult[1] != null) {
                result.add(startEndresult);
            }

            if (readlength == 0) {
                break;
            }
        }

        StringBuilder resultS = new StringBuilder();
        for (Integer[] e : result) {
            if (!resultS.isEmpty()) {
                resultS.append("|");
            }
            resultS.append(e[0]).append("-").append(e[1]);

        }
        return resultS.toString();
    }

    public static String reverseCompliment(String seq) {
        StringBuilder compliment = new StringBuilder();
        for (int i = 0; i < seq.length(); i++) {
            char nuk = seq.charAt(i);
            if (nuk == 'A') {
                compliment.append('T');
            } else if (nuk == 'T') {
                compliment.append('A');
            } else if (nuk == 'G') {
                compliment.append('C');
            } else if (nuk == 'C') {
                compliment.append('G');
            } else {
                compliment.append(nuk);
                System.out.println("blöd gelaufen " + nuk);
            }
        }
        String ret = compliment.reverse().toString();
        return ret;
    }

    public static HashMap<String, Integer> simulateDifferentialExpression(HashMap<String, Integer> transcriptMap){

        int numGenes = 10000; // total number of genes
        double mean = 50; // mean expression level
        double foldChange = 2; // fold change between conditions
        double dispersion = 0.5; // dispersion parameter


        //NegativeBinomial nbDistCond1 = new NegativeBinomialDistribution(mean, dispersion);
        //NegativeBinomialDistribution nbDistCond2 = new NegativeBinomialDistribution(mean * foldChange, dispersion);


        //int[] geneCountsCond1 = nbDistCond1.sample(numGenes);
       return transcriptMap;

    }




}
