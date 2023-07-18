package Assignment2;

public class CheckPos {

    public static void main(String[] args) {
        String seq = "GGCGGAGGTTGCAGTTGCAGTTGCAGTGAGCTGAGATCACGCCACTGCACTCCACTTGGGCGACAGAGCGAGACTCTGTCTCCAAAAAATAAATTAGTAAATAAAAGATATGAGTAAAGATTGCCAAGAAGTTCATTGGCGGCCTCTGTTTTGTTTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTTGAACAGTTGCCAGTAGTCATCAGAGTACCAAAACTGATCTATTTCTCAGTAATGAGTGTGTGCCTCATGCCTGTTTCAATATTGGGTTTTGGAGACATTATTGTACCAGGCCTGTTGATTGCATACTGTAGAAGATTT";
        System.out.println(reverseCompliment(seq));
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
}
