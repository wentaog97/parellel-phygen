import mpi.*;
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

    public static void initializeFullMatrix(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                matrix[i][j] = -1.0;
            }
        }
    }

    public static void main(String[] args) throws MPIException {
        // initialize MPI environment
        MPI.Init(args);

        // get mpi data
        int rank = MPI.COMM_WORLD.Rank();
        int size = MPI.COMM_WORLD.Size();

        // timer setup
        long startTime = 0; 
        long endTime = 0;
        if (rank == 0) startTime = System.currentTimeMillis();

        // argument handling
        if (args.length == 0) {
            System.out.println("Proper Usage is: java program filename");
            System.exit(0);
        }

        // input file handling
        String input = args[0];
        BufferedReader inputFile = null;
        try {
            inputFile = new BufferedReader(new FileReader(input));
            // read number of sequences and length
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

            // read sequences
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

            // divide work among processes based on rank
            int rowsPerProcess = numSequences / size;
            int remainder = numSequences % size;
            int startRow = rank * rowsPerProcess + Math.min(rank, remainder);
            int endRow = (rank + 1) * rowsPerProcess + Math.min(rank + 1, remainder);
            if (rank == size - 1) endRow = numSequences;

            // initialize local matrix to store computed distances for assigned rows of each rank
            double[][] localMatrix = new double[numSequences][numSequences];
            for (int i = startRow; i < endRow; ++i) {
                for (int j = i + 1; j < numSequences; ++j) {
                    double distance = computeDistance(sequences.get(i).data, sequences.get(j).data);
                    localMatrix[i][j] = distance;
                    localMatrix[j][i] = distance;
                }
            }

            // allocate full matrix on root process for collecting results
            double[][] fullMatrix = null;
            if (rank == 0) {
                fullMatrix = new double[numSequences][numSequences];
                initializeFullMatrix(fullMatrix);
            }

            // non-blocking communication
            // send computed rows from non-root processes to root process
            if (rank != 0) {
                for (int i = 0; i < (endRow - startRow); i++) {
                    MPI.COMM_WORLD.Isend(localMatrix[i], 0, numSequences, MPI.DOUBLE, 0, startRow + i);
                }
            } else {
                // root process copies its own data to full matrix
                for (int i = startRow; i < endRow; i++) {
                    fullMatrix[i] = Arrays.copyOf(localMatrix[i - startRow], numSequences);
                }
                // root process receives data from other processes
                Request[] requests = new Request[numSequences];
                for (int src = 1; src < size; src++) {
                    int startRowRecv = src * rowsPerProcess + Math.min(src, remainder);
                    int endRowRecv = (src + 1) * rowsPerProcess + Math.min(src + 1, remainder);
                    for (int i = startRowRecv; i < endRowRecv; i++) {
                        double[] rowBuffer = new double[numSequences];
                        requests[i] = MPI.COMM_WORLD.Irecv(rowBuffer, 0, numSequences, MPI.DOUBLE, src, i);
                        requests[i].Wait();
                        fullMatrix[i] = rowBuffer;
                    }
                }
            }

            // root process prints final distance matrix
            if (rank == 0) {
                DecimalFormat df = new DecimalFormat("0.0000");
                for (int i = 0; i < numSequences; ++i) {
                    System.out.printf("%10s ", sequences.get(i).name);
                    for (int j = 0; j < numSequences; ++j) {
                        if (fullMatrix[i][j] == -1.0) {
                            System.out.printf("%8s", "N/A");
                        } else {
                            System.out.printf("%8s", df.format(fullMatrix[i][j]));
                        }
                    }
                    System.out.println();
                }
            }
        } catch (Exception e) {
            System.err.println("Error handling file or data: " + e.getMessage());
        } finally {
            if (inputFile != null) try { inputFile.close(); } catch (IOException e) { }
        }

        // measure and print elapsed time
        if (rank == 0) {
            endTime = System.currentTimeMillis();
            long elapsedTime = endTime - startTime;
            System.err.printf("Elapsed time: %.3f seconds %n", elapsedTime / 1000.0);
        }

        // end all processes
        MPI.Finalize();
    }
}