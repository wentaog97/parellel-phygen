import java.io.*;
import mpi.*;
import java.util.*;
import java.text.DecimalFormat;

/**
 * Parallelized version of UPGMA using MPI Java to compute the phylogenic tree.
 * 
 * By Davis Zhong
 * 12/09/2024
 */

public class UPGMA {

    static final double LARGE_DISTANCE = 1e6;  // A large value to replace "N/A" distances

    // Class to represent a node in the phylogenetic tree
    static class Node {
        String name;               // Name of the taxon or cluster
        Node left;                 // Left child
        Node right;                // Right child
        double branchLengthLeft;   // Branch length to left child
        double branchLengthRight;  // Branch length to right child
        double height;             // Height of the node in the tree
        int size;                  // Number of taxa in the cluster

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

    // Function to read the distance matrix from the file
    public static void readDistanceMatrix(String filename, List<String> taxaNames, List<List<Double>> distances) {
        try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.trim().isEmpty())
                    continue;

                String[] tokens = line.trim().split("\\s+");
                if (tokens.length == 0)
                    continue;

                String name = tokens[0];
                taxaNames.add(name);

                List<Double> row = new ArrayList<>();
                for (int i = 1; i < tokens.length; i++) {
                    String token = tokens[i];
                    double dist;
                    if (token.equalsIgnoreCase("N/A")) {
                        dist = LARGE_DISTANCE;  // Replace "N/A" with a large distance
                    } else {
                        try {
                            dist = Double.parseDouble(token);
                        } catch (NumberFormatException e) {
                            System.err.println("Error: Invalid distance value '" + token + "' in the distance matrix.");
                            System.exit(1);
                            return;
                        }
                    }
                    row.add(dist);
                }
                distances.add(row);
            }
        } catch (IOException e) {
            System.err.println("Error: Could not open the distance matrix file '" + filename + "'.");
            System.exit(1);
        }
    }

    // Function to output the tree in Newick format
    public static String buildNewick(Node node) {
        if (node == null)
            return "";

        if (node.left == null && node.right == null) {
            // Leaf node
            return node.name;
        } else {
            // Internal node
            String leftSubtree = buildNewick(node.left);
            String rightSubtree = buildNewick(node.right);
            StringBuilder sb = new StringBuilder();
            DecimalFormat df = new DecimalFormat("0.######");  // Format branch lengths
            sb.append("(")
              .append(leftSubtree)
              .append(":")
              .append(df.format(node.branchLengthLeft))
              .append(",")
              .append(rightSubtree)
              .append(":")
              .append(df.format(node.branchLengthRight))
              .append(")");
            return sb.toString();
        }
    }

    // UPGMA algorithm implementation, parallelized
    public static void UPGMA(List<String> taxaNames, List<List<Double>> distances) throws MPIException {
        int n = taxaNames.size();
        List<Node> nodes = new ArrayList<>();

        // initialize clusters with individual taxa
        for (String name : taxaNames) {
            nodes.add(new Node(name));
        }
        // copy of distance matrix
        int size = distances.size();
        double[][] D = new double[size][size];
        // Populate the distance matrix D
        for (int i = 0; i < size; i++) {
            List<Double> row = distances.get(i);
            for (int j = 0; j < row.size(); j++) {
                D[i][j] = row.get(j);
            }
            // Ensure D is square
            for (int j = row.size(); j < size; j++) {
                D[i][j] = 0.0;
            }
        }
        // Ensure D is symmetric
        for (int i = 0; i < size; i++) {
            for (int j = i + 1; j < size; j++) {
                D[j][i] = D[i][j];
            }
        }
        // initialize MPI
        MPI.Init(new String[0]);
        int worldRank = MPI.COMM_WORLD.Rank();
        int worldSize = MPI.COMM_WORLD.Size();

        // timer setup
        long startTime = 0; 
        long endTime = 0;
        if (worldRank == 0) startTime = System.currentTimeMillis();
        // divide work among MPI processes
        int chunkSize = (n + worldSize - 1) / worldSize;  // round up to ensure each process gets some work
        // create list to hold active indices of taxa for each process
        List<Integer> activeIndices = new ArrayList<>();
        for (int i = worldRank * chunkSize; i < Math.min((worldRank + 1) * chunkSize, n); i++) {
            activeIndices.add(i);
        }
        
        while (activeIndices.size() > 1) {
            int m = activeIndices.size();
            
            // find pair of clusters with smallest distance in the local data
            double minDist = Double.MAX_VALUE;
            int minI = -1, minJ = -1;
            for (int ii = 0; ii < m - 1; ii++) {
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
           
            double[] minDistArr = {minDist};
            int[] minIArr = {minI};
            int[] minJArr = {minJ};
            // perform Allreduce to get minimum across all processes
            MPI.COMM_WORLD.Allreduce(minDistArr, 0, minDistArr, 0, 1, MPI.DOUBLE, MPI.MIN);
            MPI.COMM_WORLD.Allreduce(minIArr, 0, minIArr, 0, 1, MPI.INT, MPI.MIN);
            MPI.COMM_WORLD.Allreduce(minJArr, 0, minJArr, 0, 1, MPI.INT, MPI.MIN);
            
            // continue merging process with parallelized communication
            minDist = minDistArr[0];
            minI = minIArr[0];
            minJ = minJArr[0];
            // merge clusters[minI] and clusters[minJ]
            Node clusterI = nodes.get(minI);
            Node clusterJ = nodes.get(minJ);
            Node newCluster = new Node("");
            newCluster.left = clusterI;
            newCluster.right = clusterJ;
            double newHeight = minDist / 2.0;
            newCluster.height = newHeight;
            newCluster.branchLengthLeft = newHeight - clusterI.height;
            newCluster.branchLengthRight = newHeight - clusterJ.height;
            newCluster.size = clusterI.size + clusterJ.size;
            nodes.add(newCluster);
            int newIdx = nodes.size() - 1;
            // resize distance matrix D
            double[][] newD = new double[D.length + 1][D.length + 1];
            for (int i = 0; i < D.length; i++) {
                for (int j = 0; j < D.length; j++) {
                    newD[i][j] = D[i][j];
                }
            }
            // initialize new row and column for new cluster
            for (int i = 0; i < newD.length; i++) {
                newD[i][newD.length - 1] = 0.0;
                newD[newD.length - 1][i] = 0.0;
            }
            // update new distance matrix for all existing clusters
            for (int idxK : activeIndices) {
                if (idxK == minI || idxK == minJ) continue;
                Node clusterK = nodes.get(idxK);
                double dist = (D[minI][idxK] * clusterI.size + D[minJ][idxK] * clusterJ.size) / (clusterI.size + clusterJ.size);
                // update distance between new cluster and existing clusters
                newD[newIdx][idxK] = dist;
                newD[idxK][newIdx] = dist;
            }
            // set updated distance matrix to D
            D = newD;
            // remove merged clusters and add new one
            activeIndices.remove(Integer.valueOf(minI));
            activeIndices.remove(Integer.valueOf(minJ));
            activeIndices.add(newIdx);
            
            // synchronize after each merge
            MPI.COMM_WORLD.Barrier();
        }
        // The last remaining cluster is the root
        int rootIdx = activeIndices.get(0);
        Node root = nodes.get(rootIdx);
        // Output the tree in Newick format
        if (worldRank == 0) {
            String newickTree = buildNewick(root) + ";";
            System.out.println(newickTree);
        }
        // measure and print elapsed time
        if (worldRank == 0) {
            endTime = System.currentTimeMillis();
            long elapsedTime = endTime - startTime;
            System.err.printf("Elapsed time: %.3f seconds %n", elapsedTime / 1000.0);
        }
        MPI.Finalize();
    }


    public static void main(String[] args) throws MPIException {
        if (args.length == 0) {
            System.out.println("Proper Usage is: java program filename");
            System.exit(0);
        }
        
        List<String> taxaNames = new ArrayList<>();
        List<List<Double>> distances = new ArrayList<>();

        String filePath = args[0];

        // Read the distance matrix from the file "DistanceMatrix"
        readDistanceMatrix(filePath, taxaNames, distances);

        // Perform UPGMA algorithm
        UPGMA(taxaNames, distances);
    }
}
