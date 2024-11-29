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

const double LARGE_DISTANCE = 1e6;

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

    const double LARGE_DISTANCE = 1e6; // Define a large distance value

    string line;
    while (getline(infile, line)) {
        if (line.empty())
            continue;

        istringstream iss(line);
        string name;
        iss >> name;
        taxa_names.push_back(name);

        vector<double> row;
        string token;
        while (iss >> token) {
            double dist;
            if (token == "N/A" || token == "n/a") {
                dist = LARGE_DISTANCE; // Replace "N/A" with a large distance
            } else {
                // Convert the token to a double
                istringstream token_stream(token);
                if (!(token_stream >> dist)) {
                    cerr << "Error: Invalid distance value '" << token << "' in the distance matrix." << endl;
                    exit(1);
                }
            }
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

void UPGMA(vector<string>& taxa_names,
           vector<vector<double>>& distances) {
    int n = taxa_names.size();
    vector<Node*> nodes;

    // Initialize clusters with individual taxa
    for (int i = 0; i < n; ++i) {
        nodes.push_back(new Node(taxa_names[i]));
    }

    vector<vector<double>> D = distances;  // Copy of the distance matrix

    // Ensure D is square and symmetric
    for (int i = 0; i < n; ++i) {
        D[i].resize(n, 0.0);
    }

    vector<int> active_indices;
    for (int i = 0; i < n; ++i) {
        active_indices.push_back(i);
    }

    while (active_indices.size() > 1) {
        int m = active_indices.size();

        // Find the pair of clusters with the smallest distance
        double min_dist = numeric_limits<double>::max();
        int min_i = -1, min_j = -1;
        for (int ii = 0; ii < m - 1; ++ii) {
            int idx_i = active_indices[ii];
            for (int jj = ii + 1; jj < m; ++jj) {
                int idx_j = active_indices[jj];
                if (D[idx_i][idx_j] < min_dist) {
                    min_dist = D[idx_i][idx_j];
                    min_i = idx_i;
                    min_j = idx_j;
                }
            }
        }

        // Merge clusters[min_i] and clusters[min_j]
        Node* cluster_i = nodes[min_i];
        Node* cluster_j = nodes[min_j];

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

        // Add new cluster to nodes
        nodes.push_back(new_cluster);
        int new_idx = nodes.size() - 1;

        // Resize D to add new row and column
        int D_size = D.size();
        for (int i = 0; i < D_size; ++i) {
            D[i].resize(D_size + 1, 0.0);
        }
        D.push_back(vector<double>(D_size + 1, 0.0));

        // Compute distances between the new cluster and other active clusters
        for (int idx_k : active_indices) {
            if (idx_k == min_i || idx_k == min_j)
                continue;
            Node* cluster_k = nodes[idx_k];
            double dist = (D[min_i][idx_k] * cluster_i->size + D[min_j][idx_k] * cluster_j->size) / (cluster_i->size + cluster_j->size);

            // Update D[new_idx][idx_k] and D[idx_k][new_idx]
            D[new_idx][idx_k] = dist;
            D[idx_k][new_idx] = dist;
        }

        // Mark distances of merged clusters as inactive
        for (int i = 0; i < D.size(); ++i) {
            D[min_i][i] = D[i][min_i] = numeric_limits<double>::max();
            D[min_j][i] = D[i][min_j] = numeric_limits<double>::max();
        }

        // Update active indices
        active_indices.erase(remove(active_indices.begin(), active_indices.end(), min_i), active_indices.end());
        active_indices.erase(remove(active_indices.begin(), active_indices.end(), min_j), active_indices.end());
        active_indices.push_back(new_idx);
    }

    // The last remaining cluster is the root
    int root_idx = active_indices[0];
    Node* root = nodes[root_idx];

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
