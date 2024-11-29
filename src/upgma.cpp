#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <limits>
#include <map>
#include <algorithm>

using namespace std;

// Structure to represent a node in the phylogenetic tree
struct Node {
    string name;         // Name of the taxon or cluster
    Node* left;          // Left child
    Node* right;         // Right child
    double branch_length_left;   // Branch length to left child
    double branch_length_right;  // Branch length to right child
    double height;       // Height of the node in the tree
    int size;            // Number of taxa in the cluster

    Node(string n) : name(n), left(nullptr), right(nullptr),
                     branch_length_left(0.0), branch_length_right(0.0),
                     height(0.0), size(1) {}
};

// Function to read the distance matrix from the file
void readDistanceMatrix(const string& filename,
                        vector<string>& taxa_names,
                        vector<vector<double>>& distances) {
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Error: Could not open the distance matrix file '" << filename << "'." << endl;
        exit(1);
    }

    string line;
    while (getline(infile, line)) {
        if (line.empty())
            continue;

        istringstream iss(line);
        string name;
        iss >> name;
        taxa_names.push_back(name);

        vector<double> row;
        double dist;
        while (iss >> dist) {
            row.push_back(dist);
        }
        distances.push_back(row);
    }
    infile.close();
}

// Function to output the tree in Newick format
string buildNewick(Node* node) {
    if (!node)
        return "";

    if (!node->left && !node->right) {
        // Leaf node
        return node->name;
    } else {
        // Internal node
        string left_subtree = buildNewick(node->left);
        string right_subtree = buildNewick(node->right);
        ostringstream oss;
        oss << "(" << left_subtree << ":" << fixed << setprecision(6) << node->branch_length_left
            << "," << right_subtree << ":" << fixed << setprecision(6) << node->branch_length_right << ")";
        return oss.str();
    }
}

// Function to delete the tree recursively
void deleteTree(Node* node) {
    if (node == nullptr)
        return;
    deleteTree(node->left);
    deleteTree(node->right);
    delete node;
}

// UPGMA algorithm implementation
void UPGMA(vector<string>& taxa_names,
           vector<vector<double>>& distances) {
    int n = taxa_names.size();
    vector<Node*> clusters;

    // Initialize clusters with individual taxa
    for (int i = 0; i < n; ++i) {
        clusters.push_back(new Node(taxa_names[i]));
    }

    vector<vector<double>> D = distances;  // Copy of the distance matrix

    while (clusters.size() > 1) {
        int m = clusters.size();

        // Find the pair of clusters with the smallest distance
        double min_dist = numeric_limits<double>::max();
        int min_i = -1, min_j = -1;
        for (int i = 0; i < m - 1; ++i) {
            for (int j = i + 1; j < m; ++j) {
                if (D[i][j] < min_dist) {
                    min_dist = D[i][j];
                    min_i = i;
                    min_j = j;
                }
            }
        }

        // Merge clusters[min_i] and clusters[min_j]
        Node* cluster_i = clusters[min_i];
        Node* cluster_j = clusters[min_j];

        // Create new cluster
        Node* new_cluster = new Node("");
        new_cluster->left = cluster_i;
        new_cluster->right = cluster_j;

        // The height of the new cluster is half the distance between merged clusters
        double new_height = min_dist / 2.0;
        new_cluster->height = new_height;

        // Set branch lengths
        new_cluster->branch_length_left = new_height - cluster_i->height;
        new_cluster->branch_length_right = new_height - cluster_j->height;

        // Update size
        new_cluster->size = cluster_i->size + cluster_j->size;

        // Remove the merged clusters and add the new cluster
        clusters.erase(clusters.begin() + max(min_i, min_j));
        clusters.erase(clusters.begin() + min(min_i, min_j));
        clusters.push_back(new_cluster);

        // Update the distance matrix
        // Remove rows and columns corresponding to clusters[min_i] and clusters[min_j]
        D.erase(D.begin() + max(min_i, min_j));
        D.erase(D.begin() + min(min_i, min_j));
        for (auto& row : D) {
            row.erase(row.begin() + max(min_i, min_j));
            row.erase(row.begin() + min(min_i, min_j));
        }

        // Calculate distances between the new cluster and the remaining clusters
        vector<double> new_dist_row;
        for (int k = 0; k < clusters.size() - 1; ++k) {
            double dist = (D[k][k] * clusters[k]->size + D[k][k] * clusters[k]->size) / (clusters[k]->size + clusters[k]->size);
            dist = (D[k][clusters.size() - 1] * clusters[k]->size + D[k][clusters.size() - 1] * clusters[k]->size) / (clusters[k]->size + clusters[k]->size);
            dist = (D[k][k] * clusters[k]->size + min_dist * clusters[k]->size) / (clusters[k]->size + new_cluster->size);

            // Calculate average distance
            double dist_ik = (D[k][min_i] * cluster_i->size + D[k][min_j] * cluster_j->size) / (cluster_i->size + cluster_j->size);
            new_dist_row.push_back(dist_ik);
        }
        new_dist_row.push_back(0.0);  // Distance to itself is zero

        // Add the new distances to the matrix
        for (int i = 0; i < D.size(); ++i) {
            D[i].push_back(new_dist_row[i]);
        }
        D.push_back(new_dist_row);
    }

    // The last remaining cluster is the root
    Node* root = clusters[0];

    // Output the tree in Newick format
    string newick_tree = buildNewick(root) + ";";
    cout << newick_tree << endl;

    // Clean up memory
    deleteTree(root);
}

int main() {
    vector<string> taxa_names;
    vector<vector<double>> distances;

    // Read the distance matrix from the file "DistanceMatrix"
    readDistanceMatrix("DistanceMatrix", taxa_names, distances);

    // Perform UPGMA
    UPGMA(taxa_names, distances);

    return 0;
}
