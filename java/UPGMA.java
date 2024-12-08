import mpi.*;
import java.io.*;
import java.util.*;
import java.text.DecimalFormat;

public class UPGMA {

    static final double LARGE_DISTANCE = 1e6;

    // Class to represent a node in the phylogenetic tree
    static class Node {
        String name;
        Node left, right;
        double branchLengthLeft, branchLengthRight;
        double height;
        int size;

        Node(String name) {
            this.name = name;
            this.left = null;
            this.right = null;
            this.branchLengthLeft = 0.0;
            this.branchLengthRight = 0.0;
            this.height = 0.0;
            this.size = 1;
        }
    }

    // Read the distance matrix from file (same as before)
    public static void readDistanceMatrix(String filename, List<String> taxaNames, List<List<Double>> distances) {
        try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                String[] tokens = line.trim().split("\\s+");
                String name = tokens[0];
                taxaNames.add(name);
                List<Double> row = new ArrayList<>();
                for (int i = 1; i < tokens.length; i++) {
                    String token = tokens[i];
                    double dist = token.equalsIgnoreCase("N/A") ? LARGE_DISTANCE : Double.parseDouble(token);
                    row.add(dist);
                }
                distances.add(row);
            }
        } catch (IOException e) {
            System.err.println("Error: Could not open the distance matrix file '" + filename + "'.");
            System.exit(1);
        }
    }

    // Build Newick format tree (same as before)
    public static String buildNewick(Node node) {
        if (node == null) return "";
        if (node.left == null && node.right == null) return node.name;
        else {
            String leftSubtree = buildNewick(node.left);
            String rightSubtree = buildNewick(node.right);
            DecimalFormat df = new DecimalFormat("0.######");
            return "(" + leftSubtree + ":" + df.format(node.branchLengthLeft) + "," +
                   rightSubtree + ":" + df.format(node.branchLengthRight) + ")";
        }
    }

    // UPGMA algorithm with MPI parallelization
    public static void UPGMA(List<String> taxaNames, List<List<Double>> distances) throws MPIException {
        int n = taxaNames.size();
        List<Node> nodes = new ArrayList<>();

        // Initialize clusters with individual taxa
        for (String name : taxaNames) {
            nodes.add(new Node(name));
        }

        // Create distance matrix to share between processes
        int matrixSize = distances.size();
        double[][] D = new double[matrixSize][matrixSize];
        int rank = MPI.COMM_WORLD.Rank();
        int size = MPI.COMM_WORLD.Size();
        
        if (rank == 0) {
            for (int i = 0; i < matrixSize; i++) {
                for (int j = 0; j < matrixSize; j++) {
                    D[i][j] = distances.get(i).get(j);
                }
            }
        }

        // Distribute the distance matrix among all processes
        MPI.COMM_WORLD.Bcast(D, 0, matrixSize * matrixSize, MPI.DOUBLE, 0);

        List<Integer> activeIndices = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            activeIndices.add(i);
        }

        while (activeIndices.size() > 1) {
            int m = activeIndices.size();

            // Each process computes its local minimum distance
            double minDist = Double.MAX_VALUE;
            int minI = -1, minJ = -1;
            for (int ii = rank; ii < m - 1; ii += size) {
                int idxI = activeIndices.get(ii);
                for (int jj = ii + 1; jj < m; jj++) {
                    int idxJ = activeIndices.get(jj);
                    if (D[idxI][idxJ] < minDist) {
                        minDist = D[idxI][idxJ];
                        minI = idxI;
                        minJ = idxJ;
                    }
                }
            }

            // Perform reduction to get the global minimum pair
            double[] globalMinDist = new double[1];
            int[] globalMinI = new int[1];
            int[] globalMinJ = new int[1];
            MPI.COMM_WORLD.Allreduce(new double[] {minDist}, 0, globalMinDist, 0, 1, MPI.DOUBLE, MPI.MIN);
            MPI.COMM_WORLD.Allreduce(new int[] {minI}, 0, globalMinI, 0, 1, MPI.INT, MPI.MIN);
            MPI.COMM_WORLD.Allreduce(new int[] {minJ}, 0, globalMinJ, 0, 1, MPI.INT, MPI.MIN);

            // Root process merges clusters and updates distance matrix
            if (rank == 0) {
                Node clusterI = nodes.get(globalMinI[0]);
                Node clusterJ = nodes.get(globalMinJ[0]);
                Node newCluster = new Node("");
                newCluster.left = clusterI;
                newCluster.right = clusterJ;

                double newHeight = globalMinDist[0] / 2.0;
                newCluster.height = newHeight;
                newCluster.branchLengthLeft = newHeight - clusterI.height;
                newCluster.branchLengthRight = newHeight - clusterJ.height;
                newCluster.size = clusterI.size + clusterJ.size;

                nodes.add(newCluster);
                int newIdx = nodes.size() - 1;

                // Broadcast updated tree structure to all processes
                MPI.COMM_WORLD.Bcast(newCluster, 0, newCluster.size, MPI.OBJECT, 0);
                // Update active indices and distance matrix (not shown fully)

                activeIndices.remove(Integer.valueOf(globalMinI[0]));
                activeIndices.remove(Integer.valueOf(globalMinJ[0]));
                activeIndices.add(newIdx);
            }
        }

        // Final result at root
        if (rank == 0) {
            Node root = nodes.get(0);
            String newickTree = buildNewick(root) + ";";
            System.out.println(newickTree);
        }
    }

    public static void main(String[] args) throws MPIException {
        // Initialize MPI
        MPI.Init(args);
        int rank = MPI.COMM_WORLD.Rank();
        int size = MPI.COMM_WORLD.Size();

        if (args.length == 0) {
            if (rank == 0) {
                System.out.println("Proper Usage is: java program filename");
            }
            MPI.Finalize();
            System.exit(0);
        }

        List<String> taxaNames = new ArrayList<>();
        List<List<Double>> distances = new ArrayList<>();

        String filePath = args[0];
        if (rank == 0) {
            readDistanceMatrix(filePath, taxaNames, distances);
        }

        // Perform UPGMA with MPI
        UPGMA(taxaNames, distances);

        // Finalize MPI
        MPI.Finalize();
    }
}
