package Assignment2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

public class ExtractorTest {

    public static void main(String[] args) throws IOException {


        GenomeSequenceExtractor extractor = new GenomeSequenceExtractor(new File("/Users/test/IdeaProjects/GOBI-TanjaPock/data_A2/test.fa"), new File("/Users/test/IdeaProjects/GOBI-TanjaPock/data_A2/Homo_sapiens.GRCh37.75.dna.toplevel.fa.fai"));
        String result = extractor.getSequence("1", 0, 5000);
        String str = result.substring(611, 776);

        System.out.println(reverseCompliment(str));
    }

    public static String reverseCompliment(String seq){
        StringBuilder compliment = new StringBuilder();
        for(int i = 0; i < seq.length(); i++){
            char nuk = seq.charAt(i);
            if (nuk == 'A'){
                compliment.append('T');
            }else if  (nuk == 'T'){
                compliment.append('A');
            }else if (nuk == 'U'){
                compliment.append('C');
            }else if (nuk == 'C'){
                compliment.append('U');
            }
        }
        return seq.toString();
    }
}
