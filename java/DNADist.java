//package java;
import java.io.*;
import java.util.*;
import java.text.DecimalFormat;

class Sequence {
    String name;
    String data;

    Sequence(String name, String data) {
        this.name = name;
        this.data = data;
    }
}

public class DNADist {

    public static double computeDistance(String seq1, String seq2) {
        int differences = 0;
        int total = 0;

        for (int i = 0; i < seq1.length(); ++i) {
            char base1 = Character.toUpperCase(seq1.charAt(i));
            char base2 = Character.toUpperCase(seq2.charAt(i));

            // Only compare valid nucleotide pairs (A, C, G, T)
            if ((base1 == 'A' || base1 == 'C' || base1 == 'G' || base1 == 'T') &&
                (base2 == 'A' || base2 == 'C' || base2 == 'G' || base2 == 'T')) {
                total++;

                if (base1 != base2) {
                    differences++;
                }
            }
        }

        if (total == 0) {
            return -1.0;  // No valid overlap
        }

        double p = (double) differences / total;

        // Check for the maximum p value for which the model is defined
        if (p >= 0.75) {
            return -1.0;  // Infinite distance
        }

        double distance = -0.75 * Math.log(1.0 - (4.0 / 3.0) * p);
        return distance;
    }

    public static void main(String[] args) {
        String input = "seq";
        BufferedReader inputFile = null;
        try {
            inputFile = new BufferedReader(new FileReader(input));
            // Read the number of sequences and their length
            String line = inputFile.readLine();
            if (line == null) {
                System.err.println("Error: Input file is empty.");
                return;
            }
            String[] tokens = line.trim().split("\\s+");
            if (tokens.length < 2) {
                System.err.println("Error: Invalid header line.");
                return;
            }
            int numSequences = Integer.parseInt(tokens[0]);
            int sequenceLength = Integer.parseInt(tokens[1]);

            List<Sequence> sequences = new ArrayList<>();

            // Read the sequences
            for (int i = 0; i < numSequences; ++i) {
                line = inputFile.readLine();
                if (line == null) {
                    System.err.println("Error: Not enough sequences in the input file.");
                    return;
                }
                tokens = line.trim().split("\\s+");
                if (tokens.length < 2) {
                    System.err.println("Error: Invalid sequence line.");
                    return;
                }
                String name = tokens[0];
                String data = tokens[1];
                sequences.add(new Sequence(name, data));
            }
            inputFile.close();

            // Initialize the distance matrix
            double[][] distanceMatrix = new double[numSequences][numSequences];

            // Compute distances between each pair of sequences
            for (int i = 0; i < numSequences; ++i) {
                for (int j = i + 1; j < numSequences; ++j) {
                    double distance = computeDistance(sequences.get(i).data, sequences.get(j).data);
                    distanceMatrix[i][j] = distance;
                    distanceMatrix[j][i] = distance;
                }
            }

            // Output the distance matrix
            DecimalFormat df = new DecimalFormat("0.0000");
            for (int i = 0; i < numSequences; ++i) {
                System.out.printf("%10s ", sequences.get(i).name);
                for (int j = 0; j < numSequences; ++j) {
                    if (distanceMatrix[i][j] == -1.0) {
                        System.out.printf("%8s", "N/A");  // No valid distance
                    } else {
                        System.out.printf("%8s", df.format(distanceMatrix[i][j]));
                    }
                }
                System.out.println();
            }

        } catch (FileNotFoundException e) {
            System.err.println("Error: Could not open the input file 'seq'.");
            e.printStackTrace();
        } catch (IOException e) {
            System.err.println("Error reading the input file.");
            e.printStackTrace();
        } catch (NumberFormatException e) {
            System.err.println("Error: Invalid number format in input file.");
            e.printStackTrace();
        } finally {
            if (inputFile != null) {
                try {
                    inputFile.close();
                } catch (IOException e) {
                    // Ignore
                }
            }
        }
    }
}
