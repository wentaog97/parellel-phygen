import java.io.*;
import java.util.*;
import java.text.DecimalFormat;

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

    // UPGMA algorithm implementation
    public static void UPGMA(List<String> taxaNames, List<List<Double>> distances) {
        int n = taxaNames.size();
        List<Node> nodes = new ArrayList<>();

        // Initialize clusters with individual taxa
        for (String name : taxaNames) {
            nodes.add(new Node(name));
        }

        // Copy of the distance matrix
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

        List<Integer> activeIndices = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            activeIndices.add(i);
        }

        while (activeIndices.size() > 1) {
            int m = activeIndices.size();

            // Find the pair of clusters with the smallest distance
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

            // Merge clusters[minI] and clusters[minJ]
            Node clusterI = nodes.get(minI);
            Node clusterJ = nodes.get(minJ);

            // Create new cluster
            Node newCluster = new Node("");
            newCluster.left = clusterI;
            newCluster.right = clusterJ;

            // The height of the new cluster is half the distance between merged clusters
            double newHeight = minDist / 2.0;
            newCluster.height = newHeight;

            // Set branch lengths
            newCluster.branchLengthLeft = newHeight - clusterI.height;
            newCluster.branchLengthRight = newHeight - clusterJ.height;

            // Update size
            newCluster.size = clusterI.size + clusterJ.size;

            // Add new cluster to nodes
            nodes.add(newCluster);
            int newIdx = nodes.size() - 1;

            // Resize D to add new row and column
            int DSize = D.length;
            // Expand D to (DSize+1) x (DSize+1)
            double[][] newD = new double[DSize + 1][DSize + 1];
            for (int i = 0; i < DSize; i++) {
                for (int j = 0; j < DSize; j++) {
                    newD[i][j] = D[i][j];
                }
            }
            D = newD;

            // Initialize new row and column
            for (int i = 0; i <= DSize; i++) {
                D[i][DSize] = 0.0;
                D[DSize][i] = 0.0;
            }

            // Compute distances between the new cluster and other active clusters
            for (int idxK : activeIndices) {
                if (idxK == minI || idxK == minJ)
                    continue;
                Node clusterK = nodes.get(idxK);
                double dist = (D[minI][idxK] * clusterI.size + D[minJ][idxK] * clusterJ.size)
                               / (clusterI.size + clusterJ.size);

                // Update D[newIdx][idxK] and D[idxK][newIdx]
                D[newIdx][idxK] = dist;
                D[idxK][newIdx] = dist;
            }

            // Mark distances of merged clusters as inactive
            for (int i = 0; i < D.length; i++) {
                D[minI][i] = D[i][minI] = Double.MAX_VALUE;
                D[minJ][i] = D[i][minJ] = Double.MAX_VALUE;
            }

            // Update active indices
            activeIndices.remove(Integer.valueOf(minI));
            activeIndices.remove(Integer.valueOf(minJ));
            activeIndices.add(newIdx);
        }

        // The last remaining cluster is the root
        int rootIdx = activeIndices.get(0);
        Node root = nodes.get(rootIdx);

        // Output the tree in Newick format
        String newickTree = buildNewick(root) + ";";
        System.out.println(newickTree);

        // No need to explicitly delete the tree in Java (garbage collected)
    }

    public static void main(String[] args) {
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
