package Assignment2;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import org.apache.commons.math3.distribution.NormalDistribution;




public class TranscriptFragmentDistribution {
    HashMap<String,ArrayList<Integer>> transcriptFragments = new HashMap<>();
    HashMap<String,Double> transcript_mm = new HashMap<>();
    HashMap<String,Double> transcript_tot_len = new HashMap<>();
    HashMap<String, Double> transcript_readcons = new HashMap<String, Double>();
    public static void main(String[] args) {

        // To Tanja, I will send you a ordner that has the samples, you need to call them I called it in main
        //after reading the files call the calc_1T that is implemented so that you only need to call it with the list
        // number of mutation are given in the file, if you want to implement a rate divide it by the num
        // of reads . Love you <3
        // the return of the calc_1T is a normal distrubiton so you can sample n amount of fragments for n amount of reads
        // call the calc1t only if the fragment starts has more than 1 start if not use that start
        String directoryPath = "/Users/test/Documents/gene_Samples";

        // Get all files in directory
   /*
        File[] files = new File(directoryPath).listFiles();
        for (File file : files) {
            if (file.isFile() && file.getName().endsWith(".txt")) {
                readFile(file);
            }
        }
        calc();

    */
    }
    public TranscriptFragmentDistribution(String directoryPath){
        File[] files = new File(directoryPath).listFiles();
        for (File file : files) {
            if (file.isFile() && file.getName().endsWith(".txt")) {
                readFile(file);
            }
        }
    }

    public  NormalDistribution getFragStart(ArrayList<Integer> fragmentStarts){
        double mean = getMean(fragmentStarts);
        if(fragmentStarts.size()>1){
            double sd = getStandardDeviation(fragmentStarts,mean);
            NormalDistribution normDist = new NormalDistribution(mean, sd);
            for (int i = 0; i < 2; i++) {
                double sample = normDist.sample();
                double probability = normDist.density(sample);
                // System.out.println("Sampled start position: " + sample + ", Probability: " + probability);
            }
            return normDist;
        }else {
            double sd =50;
            NormalDistribution r = new NormalDistribution(mean,sd);
            return r;

        }


    }
    public  NormalDistribution getFragLen(ArrayList<Integer> fragmentStarts){
        double mean = getMean(fragmentStarts);
        if(fragmentStarts.size()>1){
            double sd = getStandardDeviation(fragmentStarts,mean);
            NormalDistribution normDist = new NormalDistribution(mean, sd);
            for (int i = 0; i < 2; i++) {
                double sample = normDist.sample();
                double probability = normDist.density(sample);
                System.out.println("Sampled start position: " + sample + ", Probability: " + probability);
            }
            return normDist;
        }else {
            System.out.println("only one sample "+ fragmentStarts.get(0));
            System.err.print("I told you to call this with more than one frag start");
            return null;

        }


    }
    public  void calc(){
        for(String transcriptId : transcriptFragments.keySet()){
            ArrayList<Integer> fragmentStarts = transcriptFragments.get(transcriptId);
            double totalLength = transcript_tot_len.get(transcriptId);
            double numReads = transcript_readcons.get(transcriptId)/2.0;
            double mean = getMean(fragmentStarts);
            if(fragmentStarts.size()>1){
                double sd = getStandardDeviation(fragmentStarts,mean);
                NormalDistribution normDist = new NormalDistribution(mean, sd);
                System.out.println("For Transcript: "+transcriptId);
                for (int i = 0; i < 2; i++) {
                    double sample = normDist.sample();
                    double probability = normDist.density(sample);
                    System.out.println("Sampled start position: " + sample + ", Probability: " + probability);
                }
            }else {
                System.out.println("only one sample "+ fragmentStarts.get(0));

            }



        }

    }


    public  void readFile(File filename){
        try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
            String line;
            while ((line = br.readLine()) != null) {
                if(!line.startsWith("#")){
                    String[] data = line.split("\t");
                    String transcriptId = data[1];
                    Double tot_read = Double.parseDouble(data[2]);
                    Double tot_len= Double.parseDouble(data[4]);
                    Double mm = Double.parseDouble(data[3]);
                    String[] fragmentStarts = data[5].split("-");
                    ArrayList<Integer> starts = new ArrayList<Integer>();
                    for (String start : fragmentStarts) {
                        starts.add(Integer.parseInt(start));
                    }
                    if(transcript_readcons.containsKey(transcriptId)){
                        double tmp = transcript_readcons.get(transcriptId);
                        tmp+=tot_read;
                        transcript_readcons.put(transcriptId,tmp);
                    }else {
                        transcript_readcons.put(transcriptId,tot_read);
                    }
                    if(transcript_mm.containsKey(transcriptId)){
                        double tmp = transcript_mm.get(transcriptId);
                        tmp+=mm;
                        transcript_mm.put(transcriptId,tmp);
                    }else {
                        transcript_mm.put(transcriptId,mm);
                    }
                    if(transcriptFragments.containsKey(transcriptId)){
                        transcriptFragments.get(transcriptId).addAll(starts);

                    }else {
                        transcriptFragments.put(transcriptId, starts);
                    }
                    if(transcript_tot_len.containsKey(transcriptId)){
                        double tmp = transcript_tot_len.get(transcriptId);
                        tmp+=tot_len;
                        transcript_tot_len.put(transcriptId,tmp);
                    }else {
                        transcript_tot_len.put(transcriptId,tot_len);
                    }
                }


            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    private static double getMean(ArrayList<Integer> numbers) {
        double sum = 0.0;
        for (double num : numbers) {
            sum += num;
        }
        int n = numbers.size();
        return sum / n;
    }

    // helper method to calculate the standard deviation of a list of numbers
    private static double getStandardDeviation(List<Integer> numbers,double mean) {
        double sumSquaredDiff = 0.0;
        for (double num : numbers) {
            sumSquaredDiff += Math.pow(num - mean, 2);
        }
        int n = numbers.size();
        return Math.sqrt(sumSquaredDiff / (n - 1));
    }
}
/*
     ArrayList<Integer> fragmentStarts = transcriptFragments.get(transcriptId);
              int fragmentCount = fragmentStarts.size();
              Collections.sort(fragmentStarts);
              if(fragmentCount>1){
                  double totalLength = fragmentStarts.get(fragmentCount - 1) - fragmentStarts.get(0);
                  double lambda = totalLength / (double) fragmentCount;
                  PoissonDistribution distribution = new PoissonDistribution(lambda);
                  transcriptDistributions.put(transcriptId, distribution);
              }else {

                  double lambda = 250;
                  PoissonDistribution distribution = new PoissonDistribution(lambda);
                  transcriptDistributions.put(transcriptId, distribution);
              }

 */