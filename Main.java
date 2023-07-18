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

public class ReadSimulator2 {

    public static void main(String[] args) {

        //parameter treatment
        if (args.length == 0) {
            System.out.println("execute the program with the following parameters:");
            System.out.println("-gene_id <String> : gene id");
            System.out.println("-length <int> : read length");
            System.out.println("-frlength <int>  : fragment length distribution");
            System.out.println("- SD <int> : fragment length distribution");
            System.out.println("-readcounts <file> : table of gene_id, transcript_id, count tuples");
            System.out.println("-non_unique_frq <int> : frequency how many non unique reads should be simulate, default: 20 %");
            System.out.println("-sequence_error_rate <int> : sequence error rate in percent");
            System.out.println("-dot_mutation_rate <int[4]> : dot mutation rate in percent for all bases A,C,G,T , default: 1,1,1,1");
            System.out.println("-insertion_mutation_rate <int[4]> : insertion mutation rate in percent for all bases A,C,G,T , default: 1,1,1,1");
            System.out.println("-deletion_mutation_rate <int[4]> : deletion mutation rate in percent for all bases A,C,G,T , default: 1,1,1,1");
            System.out.println("-junction_mutation_rate <int> : junction mutation rate in percent , default: 5");
            System.out.println("-junction_size <int> : affected bases at junction boundaries, default: 3");
            System.out.println("-overlap_on_more_than_1_T <boolean> : to adjust difficulty, default: false");
            System.out.println("-minoverlapsize <int> : to adjust difficulty, default: 50");
            System.out.println("-fasta <FASTA file> : genome fasta");
            System.out.println("-fidx <FASTA file index> : index for genome fasta");
            System.out.println("-gtf <GTF file> : genomic annotation");
            System.out.println("-od <output file path>: path to the output directory");
            System.exit(1);
        }

        String gene_id = null;
        int readLength = 0;
        int mean = 0;
        int standardDeviation = 0;
        String readcountsFilename = null;
        double errorate = 0.0;
        String fastaFilename = null;
        String fastaidxFilename = null;
        String gtfFilename = null;
        String outputpath = null;
        double dotMutationRates[] = new double[]{1, 1, 1, 1};
        double inMutationRates[] = new double[]{1, 1, 1, 1};
        double delMutationRates[] = new double[]{1, 1, 1, 1};
        double juncMutation = 5;
        int junctionSize = 3;
        boolean overlapOnMoreThan1T = false;
        int minoverlapsize = 0;
        int nonuniqueFrequency = 20;

        int uniqueCount = 0;
        int nonUniqueCount = 0;
        int unKnownCount = 0;

        int transcriptCountforGene = 0;
        int totalReadCount = 0;



        for (int i = 0; i < args.length - 1; i++) {

            if (args[i].equals("-length")) {
                if (!args[i + 1].startsWith("-")) {
                    readLength = Integer.parseInt(args[i + 1]);
                }
            } else if (args[i].equals("-gene_id")) {
                if (!args[i + 1].startsWith("-")) {
                    gene_id = (args[i + 1]);
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
            } else if (args[i].equals("-non_unique_frq")) {
                if (!args[i + 1].startsWith("-")) {
                    nonuniqueFrequency = Integer.parseInt(args[i + 1]);
                }
            } else if (args[i].equals("-dot_mutation_rate")) {
                if (!args[i + 1].startsWith("-")) {
                    dotMutationRates = readArrayInput(args[i + 1]);
                }
            } else if (args[i].equals("-insertion_mutation_rate")) {
                if (!args[i + 1].startsWith("-")) {
                    inMutationRates = readArrayInput(args[i + 1]);
                }
            } else if (args[i].equals("-deletion_mutation_rate")) {
                if (!args[i + 1].startsWith("-")) {
                    delMutationRates = readArrayInput(args[i + 1]);
                }
            } else if (args[i].equals("-junction_mutation_rate")) {
                if (!args[i + 1].startsWith("-")) {
                    juncMutation = Double.parseDouble(args[i + 1]);
                }
            } else if (args[i].equals("-junction_size")) {
                if (!args[i + 1].startsWith("-")) {
                    junctionSize = Integer.parseInt(args[i + 1]);
                }
            } else if (args[i].equals("-sequence_error_rate")) {
                if (!args[i + 1].startsWith("-")) {
                    errorate = Double.parseDouble(args[i + 1]);
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
            } else if (args[i].equals("-overlap_on_more_than_1_T")) {
                if (!args[i + 1].startsWith("-")) {
                    if (args[i + 1].equals("true")) {
                        overlapOnMoreThan1T = true;
                    }
                }
            } else if (args[i].equals("-minoverlapsize")) {
                if (!args[i + 1].startsWith("-")) {
                    minoverlapsize = Integer.parseInt(args[i + 1]);
                }
            } else if (args[i].equals("-transcript_count_for_gene")) {
                if (!args[i + 1].startsWith("-")) {
                    transcriptCountforGene = Integer.parseInt(args[i + 1]);
                }
            } else if (args[i].equals("-total_read_count")) {
                if (!args[i + 1].startsWith("-")) {
                    totalReadCount = Integer.parseInt(args[i + 1]);
                }
            }
        }
        //notifications
        if (fastaidxFilename == null) {
            System.out.println("fasta index file path is missing");
            System.exit(1);
        } else if (gtfFilename == null) {
            System.out.println("gtf file path is missing");
            System.exit(1);
        } else if (gene_id == null) {
            System.out.println("gene id is missing");
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
        } else if (outputpath == null) {
            System.out.println("outputhpath is missing");
            System.exit(1);
        }/* else if (totalReadCount == 0) {
            System.out.println("total read count is missing");
            System.exit(1);
        }*/

        //Dateien einlesen

        //readcounts file for transcripts who should be simulated
        ArrayList<String[]> readcounts = new ArrayList<>();
        if (readcountsFilename != null) {
            readcounts = FileUtils.readCounts(readcountsFilename);
        }
        boolean needOriginal = false;
        if(minoverlapsize!=0 || overlapOnMoreThan1T== true){
            needOriginal = true;
        }

        HashMap<String, Gene> gtfMap;
        gtfMap = FileUtils.readGTFGeneBased(gtfFilename, gene_id, readcounts); //reads rows with gene and exon as feature

        readcounts = completeReadCount(gtfMap, gene_id, transcriptCountforGene, readcounts, totalReadCount);
        addOverlapsforGene(gtfMap, gene_id, overlapOnMoreThan1T, minoverlapsize);
        addNonOverlapsforGene(gtfMap, gene_id);
        if(needOriginal){
            addOriginalOverlapsforGene(gtfMap,gene_id);
            addOriginalNonOverlapsforGene(gtfMap,gene_id);
        }

        try {
            GenomeSequenceExtractor extractor = new GenomeSequenceExtractor(new File(fastaFilename), new File(fastaidxFilename));

            String geneSequence = "";
            int geneStart = 0;

            int rmember = 0;
            HashMap<String, String> geneMap = new HashMap<>();
            HashMap<String, String> original_geneMap = new HashMap<>();

            FileWriter writefileFW = new FileWriter(outputpath + "/fw.fastq");
            FileWriter writefileRW = new FileWriter(outputpath + "/rw.fastq");
            FileWriter writefileRM = new FileWriter(outputpath + "/read.mappinginfo");
            FileWriter writefileAdditonal = new FileWriter(outputpath + "/read.additionalinfo");

            writefileRM.write("readid\tchr\tgene\ttranscript\tt_fw_regvec\tt_rw_regvec\tfw_regvec\trw_regvec\tfw_mut\trw_mut" + "\n");
            String[] inoutAdd = new String[]{"gene id: " + gene_id + "\tread length: " + readLength + "\tfragment length average: " + mean};
            FileUtils.writeTsvLine(writefileAdditonal, inoutAdd);
            writefileAdditonal.write("Transkript IDS from reads\toverlap region in percent\n");
            ArrayList<String[]> otherGenes = new ArrayList<>();


            for (int i = 0; i < readcounts.size(); i++) {
                //System.out.println(readcounts.get(i)[0]);
                String geneID = readcounts.get(i)[0];
                String transcriptID = readcounts.get(i)[1];
                StringBuilder transcirptSequenceBuilder = new StringBuilder();
                String chr = "";

                ArrayList<Integer[]> exonsInTrans = new ArrayList<>();

                Gene g = gtfMap.get(geneID);
                if (!geneMap.containsKey(geneID)) {
                    geneSequence = extractor.getSequence(g.data.seqname, g.data.start, g.data.end); //get gene Sequence
                    original_geneMap.put(geneID, geneSequence);
                    //apply mutation
                    Mutation mutation = simulateJunctionsMutation(g, geneSequence, juncMutation, junctionSize);
                    geneSequence = mutation.seq;
                    mutation = simulateDotMutation(g, geneSequence, dotMutationRates);
                    geneSequence = mutation.seq;
                    mutation = simulateIndelMutation(g, geneSequence, inMutationRates, delMutationRates);
                    geneSequence = mutation.seq;
                    geneMap.put(geneID, geneSequence);
                } else {
                    geneSequence = geneMap.get(geneID);
                }
                //System.out.println(gene_id);
                //System.out.println(geneSequence);
                chr = g.data.seqname;
                geneStart = g.data.start;
                Transcript t = g.transcriptMap.get(transcriptID);
                //System.out.println(transcriptID);
                String transcriptSeq = "";
                String originalTranscriptSeq = "";

                transcriptSeq = getTranscriptSequence(t, geneSequence, geneStart, exonsInTrans);
                originalTranscriptSeq = getTranscriptSequence(t, original_geneMap.get(geneID), geneStart, new ArrayList<>());

                //System.out.println(transcriptID);
                //System.out.println(transcriptSeq);

                int readcount = 0;

                if (readcounts.get(i).length > 2) {
                    readcount = Integer.parseInt(readcounts.get(i)[2]);
                }
                //System.out.println(readcount);

                //int rmember =0;

                Random random = new Random();
                String type = "";

                for (int r = 0; r < readcount; r++) {

                    int transcriptLen = transcriptSeq.length();
                    int fragmentLen = getSampleFragLen(readLength, mean, standardDeviation, transcriptLen);
                    //int pos =selectRandomPos(transcriptLen, fragmentLen);

                    int pos = -1;
                    if(geneID.equals(gene_id)){
                        if (random.nextDouble() * 100 < nonuniqueFrequency) {
                            pos = selectPosFromOverlap(transcriptLen, fragmentLen, t);
                            //type = "NonUnique area";
                            //nonUniqueCount += 1;
                        } else {
                            pos = selectPosFromNonOverlap(transcriptLen, fragmentLen, t);
                            //type = "Unique areas included";
                            //uniqueCount += 1;
                        }
                        if (pos == -2) {
                            pos = selectPosFromOverlap(transcriptLen, fragmentLen, t);
                            //type = "NonUnique because including unique was not possible";
                            //uniqueCount -= 1;
                            //nonUniqueCount += 1;
                        }
                        if (pos == -1) {
                            pos = selectRandomPos(transcriptLen, fragmentLen);

                        }
                        type = checkPos(t,pos, needOriginal);
                        if(type.equals("Unique area")){
                            uniqueCount += 1;
                        }else if(type.equals("NonUnique area")){
                            nonUniqueCount += 1;
                        }

                    }else{
                        pos = selectRandomPos(transcriptLen, fragmentLen);
                        type = "Random";
                    }

                    String fragmentSeq = getFragmentSeq(transcriptSeq, pos, (pos + fragmentLen));
                    String orginalFragmentSeq = getFragmentSeq(originalTranscriptSeq, pos, (pos + fragmentLen));
                    if (g.data.strand.equals("-")) {
                        fragmentSeq = reverseCompliment(fragmentSeq);
                        orginalFragmentSeq = reverseCompliment(orginalFragmentSeq);
                    }
                    String[] seq = getReadSeq(fragmentSeq, readLength);
                    String[] originalSeq = getReadSeq(orginalFragmentSeq, readLength);

                    int diffFW = compareReads(seq[0], originalSeq[0]);
                    int diffRW = compareReads(seq[1], originalSeq[1]);

                    //error rate
                    SequenceError fw = simulateSequenceErrors(seq[0], errorate);
                    SequenceError rw = simulateSequenceErrors(seq[1], errorate);


                    FileUtils.writeFastq(rmember, fw.seq, writefileFW);
                    FileUtils.writeFastq(rmember, rw.seq, writefileRW);

                    /*String fwVectors = getVector(exonsInTrans, pos, readLength);    //mit 810 getestet
                    int rwPos = (pos + fragmentLen) - readLength;
                    String rwVectors = getVector(exonsInTrans, rwPos, readLength); //mit 858 getestet

                     */

                    /*String fwVectors = vector_pos(g.data.strand, exonsInTrans,pos,readLength);
                    int rwPos = (pos + fragmentLen) - readLength;
                    String rwVectors = vector_pos(g.data.strand,exonsInTrans, rwPos, readLength);

                     */



                    if (g.data.strand.equals("-")) {
                        int r_rev_s= transcriptSeq.length()-pos;
                        int fw_start= r_rev_s-readLength;
                        int rw_start=r_rev_s-fragmentLen;
                        //int rwPos = (pos + fragmentLen) - readLength;
                        String fwVectors = getVector(exonsInTrans, fw_start, readLength);    //mit 810 getestet
                        String rwVectors = getVector(exonsInTrans, rw_start, readLength); //mit 858 getestet

                        int rwStart = transcriptLen - pos - readLength;

                        FileUtils.writeReadMapping(getReadMappingInformation(g.data.strand, rmember, chr, geneID, transcriptID, rwStart, readLength,
                                fragmentLen, rwVectors, fwVectors, fw.positions, rw.positions, type, diffFW, diffRW), writefileRM);
                    } else {
                        String fwVectors = getVector(exonsInTrans, pos, readLength);    //mit 810 getestet
                        int rwPos = (pos + fragmentLen) - readLength;
                        String rwVectors = getVector(exonsInTrans, rwPos, readLength); //mit 858 getestet

                        FileUtils.writeReadMapping(getReadMappingInformation(g.data.strand, rmember, chr, geneID, transcriptID, pos, readLength,
                                fragmentLen, fwVectors, rwVectors, fw.positions, rw.positions, type, diffFW, diffRW), writefileRM);
                    }
                    //FileUtils.writeFragmentlength(rmember,fragmentLen,writefilePlotFL);
                    //FileUtils.writeMutationDis(rmember, fw.positions, rw.positions, writefileMD);
                    rmember += 1;

                }

                String[] input = new String[2];
                if(needOriginal){
                    input = new String[]{transcriptID, String.valueOf(t.originaloverlapinPercentage)};
                }else{
                     input = new String[]{transcriptID, String.valueOf(t.overlapinPercentage)};
                }

                if(geneID.equals(gene_id)){
                    FileUtils.writeTsvLine(writefileAdditonal, input);
                }else{
                    String[] otherGene = new String[]{geneID, transcriptID, readcounts.get(i)[2] };
                    otherGenes.add(otherGene);
                }

            }
            String[] unique_RC = new String[]{"reads in unique areas: ", String.valueOf(uniqueCount)};
            FileUtils.writeTsvLine(writefileAdditonal, unique_RC);
            String[] nonunique_RC = new String[]{"reads in non unique areas: ", String.valueOf(nonUniqueCount)};
            FileUtils.writeTsvLine(writefileAdditonal, nonunique_RC);
            String[] unknown_RC = new String[]{"random area: ", String.valueOf(unKnownCount)};
            FileUtils.writeTsvLine(writefileAdditonal, unknown_RC);
            String[] seqErrorR = new String[]{"sequence error rate:", String.valueOf(errorate)};
            FileUtils.writeTsvLine(writefileAdditonal, seqErrorR);
            String[] dotMutationRate = new String[]{"dotmutation rate:", "A: ", String.valueOf(dotMutationRates[0]), "\tC: ", String.valueOf(dotMutationRates[1]), "\tG: ", String.valueOf(dotMutationRates[2]), "\tT: ", String.valueOf(dotMutationRates[3])};
            FileUtils.writeTsvLine(writefileAdditonal, dotMutationRate);
            String[] inMutationRate = new String[]{"insertion rate:", "A: ", String.valueOf(inMutationRates[0]), "\tC: ", String.valueOf(inMutationRates[1]), "\tG: ", String.valueOf(inMutationRates[2]), "\tT: ", String.valueOf(inMutationRates[3])};
            FileUtils.writeTsvLine(writefileAdditonal, inMutationRate);
            String[] delMutationRate = new String[]{"deletion rate:", "A: ", String.valueOf(delMutationRates[0]), "\tC: ", String.valueOf(delMutationRates[1]), "\tG: ", String.valueOf(delMutationRates[2]), "\tT: ", String.valueOf(delMutationRates[3])};
            FileUtils.writeTsvLine(writefileAdditonal, delMutationRate);
            writefileAdditonal.write("junction area mutation: " + String.valueOf(juncMutation) + "\t\tjunction area size: " + String.valueOf(junctionSize) + "\n");
            writefileAdditonal.write("overlap reads fit to more than two transcripts: " + String.valueOf(overlapOnMoreThan1T) + "\n");
            writefileAdditonal.write("reads from other genes: \n" );
            for(String[] ele : otherGenes){
                writefileAdditonal.write(ele[0]+"\t "+ ele[1]+"\t "+ ele[2] );
            }



            writefileFW.close();
            writefileRW.close();
            writefileRM.close();
            writefileAdditonal.close();
            //writefilePlotFL.close();
            //writefileMD.close();

        } catch (IOException e) {
            System.out.println("extractions failed: " + e.getMessage());
        }
    }

    private static String getTranscriptSequence(Transcript t, String geneSequence, int geneStart, ArrayList<Integer[]> exonsInTrans) {

        StringBuilder transcriptSequenceBuilder = new StringBuilder();
        String transcriptSeq = "";
        for (GtfRow exon : t.exons) {
            String exonSeq = "";
            if (exon.start - geneStart < geneSequence.length()) {
                int end = exon.end - geneStart + 1;
                if (end > geneSequence.length()) {
                    end = geneSequence.length();
                }
                exonSeq = geneSequence.substring(exon.start - geneStart, end);
            }

            Integer[] eChrPos = new Integer[2];
            eChrPos[0] = exon.start;
            eChrPos[1] = exon.end;
            exonsInTrans.add(eChrPos);

            transcriptSequenceBuilder.append(exonSeq);
        }
        transcriptSeq = transcriptSequenceBuilder.toString();

        if (t.exons.get(1).strand.equals("-")) {
            StringBuilder rev = new StringBuilder();
            for (int k = transcriptSeq.length() - 1; k >= 0; k--) {
                char c = transcriptSeq.charAt(k);
                if (c == 'T') {
                    rev.append('A');
                }
                if (c == 'A') {
                    rev.append('T');
                }
                if (c == 'C') {
                    rev.append('G');
                }
                if (c == 'G') {
                    rev.append('C');
                }
            }
            //rev.reverse();
            transcriptSeq = rev.toString();

        }
        return transcriptSeq;

    }

    private static int compareReads(String r1, String r2) {
        int diff = 0;
        for (int i = 0; i < r1.length(); i++) {
            if (r1.charAt(i) != r2.charAt(i)) {
                diff += 1;
            }
        }
        return diff;
    }

    private static double[] readArrayInput(String input) {
        String[] elements = input.split(",");
        if (elements.length == 1) {
            double p = Double.parseDouble(elements[0]);
            return new double[]{p, p, p, p};
        }
        if (elements.length == 4) {
            double[] result = new double[4];
            for (int i = 0; i < elements.length; i++) {
                result[i] = Double.parseDouble(elements[i]);
            }
            return result;
        }
        System.out.println("number of input values are wrong: " + input);
        System.exit(1);
        return null;
    }

    private static ArrayList<String[]> completeReadCount(HashMap<String, Gene> gtfMap, String gene_id, int transcriptCountForGene,
                                                         ArrayList<String[]> readcounts, int totalreadCount) {


        HashSet<String> transcriptsOfGivenGene = new HashSet<>();
        int sumReadCounts = 0;
        Gene givenGene = gtfMap.get(gene_id);
        boolean tCisZero = false;
        if (transcriptCountForGene == 0) {
            tCisZero = true;
        }
        for (String[] row : readcounts) {
            if (row[0].equals(gene_id)) {
                transcriptsOfGivenGene.add(row[1]);
                if (row.length > 2) {
                    sumReadCounts += Integer.parseInt(row[2]);
                    //transcriptCountForGene-=1;
                }
                if (tCisZero) {       //so that transcript count is optional
                    transcriptCountForGene += 1;
                }
            }
        }


        totalreadCount = totalreadCount - sumReadCounts;
        int meanCounts = totalreadCount / transcriptCountForGene;
        //NormalDistribution nd = new NormalDistribution(meanCounts, (meanCounts / 10));


        int entries = readcounts.size();
        for (int i = 0; i < readcounts.size(); i++) {
            String[] row = readcounts.get(i);
            if (row[0].equals(gene_id)) {
                transcriptsOfGivenGene.add(row[1]);
                if (row.length == 2 || Integer.parseInt(row[2]) == 0) {
                    if (i == entries - 1) {
                        String[] newElement = {row[0], row[1], String.valueOf(totalreadCount)};
                        readcounts.set(i, newElement);
                    } else {
                        NormalDistribution nd = new NormalDistribution(meanCounts, (meanCounts / 10));
                        int count = (int) nd.sample();
                        String[] newElement = {row[0], row[1], String.valueOf(count)};
                        readcounts.set(i, newElement);
                        totalreadCount -= count;
                    }
                }
            }
        }


        transcriptCountForGene = transcriptCountForGene - transcriptsOfGivenGene.size();
        ArrayList<String> transcriptList = new ArrayList<>(givenGene.transcriptMap.keySet());
        Collections.shuffle(transcriptList);

        for (String t : transcriptList) {
            if (transcriptCountForGene <= 0) {
                break;
            }
            if (!transcriptsOfGivenGene.contains(t)) {
                NormalDistribution nd = new NormalDistribution(meanCounts, (meanCounts / 10));
                int count = (int) nd.sample();
                if (transcriptCountForGene == 1) {
                    count = totalreadCount;
                }
                String[] newElement = {gene_id, t, String.valueOf(count)}; // meanCount add standard deviation
                readcounts.add(newElement);
                transcriptCountForGene--;
                totalreadCount -= count;
            }
        }
        return readcounts;
    }

    public static String vector_pos(String strand, ArrayList<Integer[]> exonsInTrans,  int fr_start,  int len){

        ArrayList<String> exons = new ArrayList<>();
        int ex_lens=0;
        int read_len=0;

        if(strand.equals("-")){

            for(Integer[] e: exonsInTrans){
                int startE = e[0];
                int endE = e[1];
                int lenE = endE - startE + 1;
                ex_lens+=lenE;
                if(ex_lens>=fr_start){
                    int ex_start= endE-(ex_lens-fr_start)+read_len+1;
                    int increment=(ex_lens-fr_start);
                    int ex_end=ex_start+increment;
                    if(ex_end>endE+1){
                        ex_end=endE+1;
                    }
                    read_len+=(ex_end-ex_start);
                    if(read_len>=len){
                        ex_end=ex_end-(read_len-len);
                        if(ex_end-ex_start>0){
                            String exs= ex_start+"-"+ex_end;
                            exons.add(exs);
                        }
                        break;
                    }
                    else {
                        if(ex_end-ex_start>0){
                            String exs= ex_start+"-"+ex_end;
                            exons.add(exs);
                        }
                    }
                }
            }
        }
        else {
            for(Integer[] e: exonsInTrans){
                int startE = e[0];
                int endE = e[1];
                int lenE = endE - startE + 1;
                ex_lens+=lenE;
                if(ex_lens>=fr_start){
                    int ex_start= startE+lenE-(ex_lens-fr_start)+read_len;
                    int increment=(ex_lens-fr_start);
                    int ex_end=ex_start+increment;
                    if(ex_end>endE+1){
                        ex_end=endE+1;
                    }
                    read_len+=(ex_end-ex_start);
                    if(read_len>=len){
                        ex_end=ex_end-(read_len-len);
                        if(ex_end-ex_start>0){
                            String exs= ex_start+"-"+ex_end;
                            exons.add(exs);
                        }
                        break;
                    }
                    else {
                        if(ex_end-ex_start>0){
                            String exs= ex_start+"-"+ex_end;
                            exons.add(exs);
                        }
                    }
                }
            }
        }

        StringBuilder resultS = new StringBuilder();
        for (String e : exons) {
            if (!resultS.isEmpty()) {
                resultS.append("|");
            }
            resultS.append(e);

        }
        return resultS.toString();

    }

    public static void addOverlapsforGene(HashMap<String, Gene> gtfMap, String gene_id, boolean overlapOnMoreThan1T, int minoverlapsize) {
        Gene g = gtfMap.get(gene_id);
        for (Transcript t1 : g.transcriptMap.values()) {
            for (Transcript t2 : g.transcriptMap.values()) {
                if (!t1.equals(t2)) {
                    t1.addOverlapsFor(t2, minoverlapsize);
                }
            }
            Collections.sort(t1.overlaps, new Comparator<Overlap>() {
                @Override
                public int compare(Overlap o1, Overlap o2) {
                    return o1.startPos - o2.startPos;
                }
            });
            if (overlapOnMoreThan1T) {
                int i = 1;
                while (i < t1.overlaps.size()) {
                    Overlap olOld = t1.overlaps.get(i - 1);
                    Overlap olNew = t1.overlaps.get(i);

                    if (olOld.startPos <= olNew.startPos && olOld.endPos >= olNew.endPos) {
                        //------------
                        //   -----
                        t1.overlaps.remove(i - 1);
                    } else if (olOld.startPos >= olNew.startPos && olOld.endPos <= olNew.endPos) {
                        //   -----
                        //----------
                        t1.overlaps.remove(i);
                    } else {
                        i++;
                    }
                }

            }
        }


    }
    public static void addOriginalOverlapsforGene(HashMap<String, Gene> gtfMap, String gene_id ) {
        Gene g = gtfMap.get(gene_id);
        for (Transcript t1 : g.transcriptMap.values()) {
            for (Transcript t2 : g.transcriptMap.values()) {
                if (!t1.equals(t2)) {
                    t1.addOriginalOverlapsFor(t2);
                }
            }
            Collections.sort(t1.orignialOverlaps, new Comparator<Overlap>() {
                @Override
                public int compare(Overlap o1, Overlap o2) {
                    return o1.startPos - o2.startPos;
                }
            });

        }


    }

    public static void addNonOverlapsforGene(HashMap<String, Gene> gtfMap, String gene_id) {
        Gene g = gtfMap.get(gene_id);
        for (Transcript t : g.transcriptMap.values()) {
            int start = 0;
            int transcriptLen = 0;
            for (GtfRow exon : t.exons) {
                transcriptLen += exon.end - exon.start + 1;
            }
            //System.out.println(t.exons.get(1).attribute.get("transcript_id"));
            //System.out.println("overlaps region: ");
            int nonOverlapLen = 0;
            for (Overlap o : t.overlaps) {
                //System.out.println(o.startPos + " - " + o.endPos);
                if (o.startPos <= start) {
                    if (start < o.endPos + 1) {
                        start = o.endPos + 1;
                    }

                } else {
                    Nonoverlap nonoverlaps = new Nonoverlap();
                    nonoverlaps.startPos = start;
                    nonoverlaps.endPos = o.startPos - 1;
                    t.nonoverlaps.add(nonoverlaps);
                    start = o.endPos + 1;
                    nonOverlapLen += (nonoverlaps.endPos - nonoverlaps.startPos);
                }
            }

            if (start < transcriptLen) {
                Nonoverlap nonoverlaps = new Nonoverlap();
                nonoverlaps.startPos = start;
                nonoverlaps.endPos = transcriptLen - 1;
                t.nonoverlaps.add(nonoverlaps);
                nonOverlapLen += (nonoverlaps.endPos - nonoverlaps.startPos);
            }
            /*System.out.println("nonoverlap: ");
            for (Nonoverlap no : t.nonoverlaps) {

                System.out.println(no.startPos + " - " + no.endPos);
            }
            System.out.println(transcriptLen + " : " + nonOverlapLen);

             */

            t.overlapinPercentage = Math.round((100 - ((double) nonOverlapLen * 100) / transcriptLen) * 1000) / 1000;

        }

    }

    public static void addOriginalNonOverlapsforGene(HashMap<String, Gene> gtfMap, String gene_id) {
        Gene g = gtfMap.get(gene_id);
        for (Transcript t : g.transcriptMap.values()) {
            int start = 0;
            int transcriptLen = 0;
            for (GtfRow exon : t.exons) {
                transcriptLen += exon.end - exon.start + 1;
            }
            //System.out.println(t.exons.get(1).attribute.get("transcript_id"));
            //System.out.println("overlaps region: ");
            int nonOverlapLen = 0;
            for (Overlap o : t.orignialOverlaps) {
                //System.out.println(o.startPos + " - " + o.endPos);
                if (o.startPos <= start) {
                    if (start < o.endPos + 1) {
                        start = o.endPos + 1;
                    }

                } else {
                    Nonoverlap nonoverlaps = new Nonoverlap();
                    nonoverlaps.startPos = start;
                    nonoverlaps.endPos = o.startPos - 1;
                    t.orignialNonoverlaps.add(nonoverlaps);
                    start = o.endPos + 1;
                    nonOverlapLen += (nonoverlaps.endPos - nonoverlaps.startPos);
                }
            }

            if (start < transcriptLen) {
                Nonoverlap nonoverlaps = new Nonoverlap();
                nonoverlaps.startPos = start;
                nonoverlaps.endPos = transcriptLen - 1;
                t.orignialNonoverlaps.add(nonoverlaps);
                nonOverlapLen += (nonoverlaps.endPos - nonoverlaps.startPos);
            }

            t.originaloverlapinPercentage = Math.round((100 - ((double) nonOverlapLen * 100) / transcriptLen) * 1000) / 1000;

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

    private static int selectPosFromOverlap(int transcriptLen, int fragmentLen, Transcript t) {
        Collections.shuffle(t.overlaps);

        double random = Math.random();
        for (int i = 2; i < 5; i++) {
            int shift = (int) (random * (fragmentLen / i));

            for (Overlap o : t.overlaps) {
                if ((o.startPos + shift + fragmentLen) <= transcriptLen) {
                    return o.startPos + shift;
                } else if ((o.endPos - fragmentLen - shift >= 0)) {
                    return o.endPos - fragmentLen - shift;
                }
            }
        }

        //System.out.println("xxx  - " + t.overlaps.size());

        return -1;

    }

    private static int selectPosFromNonOverlap(int transcriptLen, int fragmentLen, Transcript t) {
        Collections.shuffle(t.nonoverlaps);
        double random = Math.random();
        for (int i = 10; i > 8; i--) {
            int shift = (int) (random * (fragmentLen / i));
            for (Nonoverlap no : t.nonoverlaps) {
                if (((no.startPos + shift) < no.endPos) && (no.startPos + shift + fragmentLen) <= transcriptLen) {
                    return no.startPos + shift;
                } else if ((no.endPos <= fragmentLen) && (no.endPos - fragmentLen - shift >= 0)) {
                    return no.endPos - fragmentLen - shift;
                }
            }
        }
        //System.out.println("xxx  - " + t.overlaps.size());

        return -2;

    }

    private static String getFragmentSeq(String transcriptseq, int start, int end) {
        int trans_len = transcriptseq.length();
        String fragmentSeq = "";
        if(trans_len< end){
            fragmentSeq = transcriptseq.substring(start, trans_len);
        }else{
            fragmentSeq = transcriptseq.substring(start, end);
        }

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

    private static SequenceError simulateSequenceErrors(String transcriptSeq, double errorRate) {

        int transciptLen = transcriptSeq.length();//fragmente len eigentlich

        StringBuilder seq = new StringBuilder();
        List<Integer> positions = new ArrayList<>();

        List<Character> aList = Arrays.asList('T', 'C', 'G');
        List<Character> tList = Arrays.asList('A', 'C', 'G');
        List<Character> cList = Arrays.asList('A', 'T', 'G');
        List<Character> gList = Arrays.asList('A', 'C', 'T');

        Random random = new Random();
        int randomChar = 0;

        for (int i = 0; i < transciptLen; i++) {

            if (random.nextDouble() * 100 < errorRate) {
                positions.add(i);
                if (transcriptSeq.charAt(i) == 'A') {
                    randomChar = random.nextInt(aList.size());
                    seq.append(aList.get(randomChar));
                } else if (transcriptSeq.charAt(i) == 'T') {
                    randomChar = random.nextInt(tList.size());
                    seq.append(tList.get(randomChar));
                } else if (transcriptSeq.charAt(i) == 'C') {
                    randomChar = random.nextInt(cList.size());
                    seq.append(cList.get(randomChar));
                } else if (transcriptSeq.charAt(i) == 'G') {
                    randomChar = random.nextInt(gList.size());
                    seq.append(gList.get(randomChar));
                } else {
                    seq.append(transcriptSeq.charAt(i)); //if there is any unknown letter
                }
            } else {
                seq.append(transcriptSeq.charAt(i));
            }
        }

        SequenceError sequenceError = new SequenceError();
        sequenceError.seq = seq.toString();
        sequenceError.positions = positions;

        return sequenceError;

    }

    private static Mutation simulateDotMutation(Gene g, String geneSeq, double[] mutationRate) {

        HashMap<Character, Double> mutationRateMap = new HashMap<>();
        mutationRateMap.put('A', mutationRate[0]);
        mutationRateMap.put('C', mutationRate[1]);
        mutationRateMap.put('G', mutationRate[2]);
        mutationRateMap.put('T', mutationRate[3]);

        double transitionsRate = 65;

        int genLen = geneSeq.length();

        StringBuilder seq = new StringBuilder();
        List<Integer> positions = new ArrayList<>();

        List<Character> aList = Arrays.asList('T', 'C');
        List<Character> tList = Arrays.asList('A', 'G');
        List<Character> cList = Arrays.asList('A', 'G');
        List<Character> gList = Arrays.asList('C', 'T');

        Random random = new Random();
        int randomChar = 0;

        for (int i = 0; i < genLen; i++) {
            Double rate = mutationRateMap.get(geneSeq.charAt(i));
            if (rate == null) {
                rate = Double.valueOf(0.0);
            }
            if (random.nextDouble() * 100 < rate) {
                positions.add(i);
                boolean isTransition = random.nextDouble() * 100 < transitionsRate;

                if (geneSeq.charAt(i) == 'A') {
                    if (isTransition) {
                        seq.append('G');
                    } else {
                        randomChar = random.nextInt(aList.size());
                        seq.append(aList.get(randomChar));
                    }
                } else if (geneSeq.charAt(i) == 'T') {
                    if (isTransition) {
                        seq.append('C');
                    } else {
                        randomChar = random.nextInt(tList.size());
                        seq.append(tList.get(randomChar));
                    }
                } else if (geneSeq.charAt(i) == 'C') {
                    if (isTransition) {
                        seq.append('T');
                    } else {
                        randomChar = random.nextInt(cList.size());
                        seq.append(cList.get(randomChar));
                    }
                } else if (geneSeq.charAt(i) == 'G') {
                    if (isTransition) {
                        seq.append('A');
                    } else {
                        randomChar = random.nextInt(gList.size());
                        seq.append(gList.get(randomChar));
                    }
                } else {
                    seq.append(geneSeq.charAt(i)); //if there is any unknown letter
                }
            } else {    //no mutation
                seq.append(geneSeq.charAt(i));
            }
        }

        Mutation mutationSequence = new Mutation();
        mutationSequence.seq = seq.toString();
        mutationSequence.positions = positions;

        return mutationSequence;

    }

    private static Mutation simulateJunctionsMutation(Gene g, String geneSeq, double mutationRate, int junctionSize) {

        HashSet<Integer> junctions = new HashSet<>();
        for (Transcript t : g.transcriptMap.values()) {
            for (GtfRow exon : t.exons) {

                for (int i = 0; i < junctionSize; i++) {
                    junctions.add(exon.start + i - g.getStart());
                    junctions.add(exon.end - i - g.getStart());
                }
            }
        }

        double transitionsRate = 65;
        int genLen = geneSeq.length();

        StringBuilder seq = new StringBuilder();
        List<Integer> positions = new ArrayList<>();

        List<Character> aList = Arrays.asList('T', 'C');
        List<Character> tList = Arrays.asList('A', 'G');
        List<Character> cList = Arrays.asList('A', 'G');
        List<Character> gList = Arrays.asList('C', 'T');

        Random random = new Random();
        int randomChar = 0;

        for (int i = 0; i < genLen; i++) {
            if (junctions.contains(i) && random.nextDouble() * 100 < mutationRate) {
                positions.add(i);
                boolean isTransition = random.nextDouble() * 100 < transitionsRate;

                if (geneSeq.charAt(i) == 'A') {
                    if (isTransition) {
                        seq.append('G');
                    } else {
                        randomChar = random.nextInt(aList.size());
                        seq.append(aList.get(randomChar));
                    }
                } else if (geneSeq.charAt(i) == 'T') {
                    if (isTransition) {
                        seq.append('C');
                    } else {
                        randomChar = random.nextInt(tList.size());
                        seq.append(tList.get(randomChar));
                    }
                } else if (geneSeq.charAt(i) == 'C') {
                    if (isTransition) {
                        seq.append('T');
                    } else {
                        randomChar = random.nextInt(cList.size());
                        seq.append(cList.get(randomChar));
                    }
                } else if (geneSeq.charAt(i) == 'G') {
                    if (isTransition) {
                        seq.append('A');
                    } else {
                        randomChar = random.nextInt(gList.size());
                        seq.append(gList.get(randomChar));
                    }
                } else {
                    seq.append(geneSeq.charAt(i)); //if there is any unknown letter
                }
            } else {    //no mutation
                seq.append(geneSeq.charAt(i));
            }
        }

        Mutation mutationSequence = new Mutation();
        mutationSequence.seq = seq.toString();
        mutationSequence.positions = positions;

        return mutationSequence;

    }

    private static Mutation simulateIndelMutation(Gene g, String geneSeq, double[] inRates, double[] delRates) {

        HashMap<Character, Double> insertRateMap = new HashMap<>();
        insertRateMap.put('A', inRates[0]);
        insertRateMap.put('C', inRates[1]);
        insertRateMap.put('G', inRates[2]);
        insertRateMap.put('T', inRates[3]);

        HashMap<Character, Double> deletionRateMap = new HashMap<>();
        deletionRateMap.put('A', delRates[0]);
        deletionRateMap.put('C', delRates[1]);
        deletionRateMap.put('G', delRates[2]);
        deletionRateMap.put('T', delRates[3]);

        int genLen = geneSeq.length();

        StringBuilder seq = new StringBuilder();
        List<Integer> positions = new ArrayList<>();

        List<Character> list = Arrays.asList('T', 'C', 'A', 'G');

        Random random = new Random();
        int randomChar = 0;

        for (int i = 0; i < genLen; i++) {
            Double insertionRate = insertRateMap.get(geneSeq.charAt(i));
            Double deletionRate = deletionRateMap.get(geneSeq.charAt(i));
            if (insertionRate == null) {
                insertionRate = Double.valueOf(0.0);
            }
            if (deletionRate == null) {
                deletionRate = Double.valueOf(0.0);
            }
            if (random.nextDouble() * 100 < insertionRate) {
                positions.add(i);
                seq.append(geneSeq.charAt(i));
                randomChar = random.nextInt(list.size());
                seq.append(list.get(randomChar));

            } else if (random.nextDouble() * 100 < deletionRate) {
                positions.add(i);

            } else {    //no mutation
                seq.append(geneSeq.charAt(i));
            }
        }

        Mutation mutationSequence = new Mutation();
        mutationSequence.seq = seq.toString();
        mutationSequence.positions = positions;

        return mutationSequence;

    }

    public static String getReadMappingInformation(String strand, int readnum, String chr, String geneID, String transcriptID, int start,
                                                   int readlength, int fragmentlength, String fwVector, String rwVector,
                                                   List<Integer> fwMpos, List<Integer> rwMpos, String type, int diffFW, int diffRW) { //start == pos, start from fragment/read
        int end = start + fragmentlength;
        String fw = start + "-" + (start + readlength);
        String rw = (end - readlength) + "-" + end;

        if (strand.equals("-")) {
            rw = start + "-" + (start + readlength);

            int rwEnd = (start + readlength);
            int fwStart = rwEnd - fragmentlength;
            int fwEnd = fwStart + readlength;

            fw = fwStart + "-" + fwEnd;
        }


        String mappingLine = readnum + "\t" + chr + "\t" + geneID + "\t" + transcriptID + "\t" + fw + "\t" + rw + "\t" + fwVector + "\t"
                + rwVector + "\t" + mutPos(fwMpos) + "\t" + mutPos(rwMpos) + "\t" + type + "\t" + diffFW + "\t" + diffRW;

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

    public static String checkPos(Transcript t, int pos, boolean needOriginal ){

        if(needOriginal){
            for(Overlap o : t.orignialOverlaps){
                if(o.startPos<= pos && pos <= o.endPos){
                    return "NonUnique area";
                }
            }
            for(Nonoverlap no : t.orignialNonoverlaps){
                if(no.startPos<= pos && pos <= no.endPos){
                    return "Unique area";
                }
            }
        }else{
            for(Overlap o : t.overlaps){
                if(o.startPos<= pos && pos <= o.endPos){
                    return "NonUnique area";
                }
            }
            for(Nonoverlap no : t.nonoverlaps){
                if(no.startPos<= pos && pos <= no.endPos){
                    return "Unique area";
                }
            }
        }



        return "random";
    }


}
