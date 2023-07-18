package Assignment2;

import Utils.FileUtils;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

public class GenomeSequenceExtractor {

    ArrayList<String[]> index = new ArrayList<>();
    RandomAccessFile raf;

    public GenomeSequenceExtractor(File fasta, File indexFile) throws FileNotFoundException {
        index = readIndex(indexFile);
        raf = new RandomAccessFile(fasta, "r");
    }

    static ArrayList<String[]> readIndex(File indexFile) {
        ArrayList<String[]> readIndex = new ArrayList<>();
        try {
            BufferedReader reader = new BufferedReader(new FileReader(indexFile));
            String line;
            while ((line = reader.readLine()) != null) {
                String[] cols = line.split("\t");
                readIndex.add(cols);
            }
            reader.close();
        } catch (IOException e) {
            System.out.println("reading readcounts file failed: " + e.getMessage());
        }
        return readIndex;
    }

    public String getSequence(String chr, int start, int end) throws IOException {
        start--;
        end--;
        for (int i = 0; i < index.size(); i++) {
            if (index.get(i)[0].equals(chr)) {
                int rowLen = Integer.parseInt(index.get(i)[3]);
                long chrStart = Long.parseLong(index.get(i)[2]);
                long numNewLines = (start / rowLen);
                //System.out.println("rowlen: " + rowLen + "chrStart: " + chrStart + "numNewLines: " + numNewLines);

                long pos = chrStart + numNewLines + start;
                raf.seek(pos);
                int numNewLinesLen = end/rowLen - start/rowLen;
                //int numNewLinesLen = (int)Math.ceil((end - start) / rowLen);
                int len = (end - start) +1 + numNewLinesLen;
                byte[] bytes = new byte[len];
                raf.read(bytes);
                StringBuilder s = new StringBuilder();
                int begin = 0;
                for (int k = 0; k < bytes.length; k++) {
                    if (bytes[k] == '\n' && k > begin ) {
                        s.append(new String(bytes, begin , k - begin, StandardCharsets.UTF_8));
                        begin = k+1;
                    }
                }
                if (bytes.length > begin ) {
                    s.append(new String(bytes, begin, bytes.length - begin, StandardCharsets.UTF_8));
                }
                //String s = new String(bytes, StandardCharsets.UTF_8);
                return s.toString();
            }

        }
        return null;
    }



}
