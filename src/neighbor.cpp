#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <limits>
#include <map>

using namespace std;

// Structure to represent a node in the phylogenetic tree
struct Node {
    string name;         // Name of the taxon (for leaf nodes)
    Node* left;          // Left child
    Node* right;         // Right child
    double branch_length_left;   // Branch length to left child
    double branch_length_right;  // Branch length to right child
    double height;       // Height of the node in the tree
    int index;           // Index of the node in the distance matrix

    Node(string n, int idx) : name(n), left(nullptr), right(nullptr),
                              branch_length_left(0.0), branch_length_right(0.0),
                              height(0.0), index(idx) {}
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

// Neighbor-Joining algorithm implementation
void neighborJoining(vector<string>& taxa_names,
                     vector<vector<double>>& distances) {
    int n = taxa_names.size();
    vector<Node*> nodes;

    // Initialize nodes
    for (int i = 0; i < n; ++i) {
        nodes.push_back(new Node(taxa_names[i], i));
    }

    vector<vector<double>> D = distances;  // Copy of the distance matrix
    vector<int> active_indices(n);
    for (int i = 0; i < n; ++i) {
        active_indices[i] = i;
    }

    int next_node_index = n;  // Index for new internal nodes

    while (active_indices.size() > 2) {
        int m = active_indices.size();
        vector<double> R(m, 0.0);

        // Compute total distances R
        for (int i = 0; i < m; ++i) {
            int idx_i = active_indices[i];
            for (int j = 0; j < m; ++j) {
                int idx_j = active_indices[j];
                R[i] += D[idx_i][idx_j];
            }
        }

        // Compute the Q-matrix and find the pair with minimum Q-value
        double min_Q = numeric_limits<double>::max();
        int min_i = -1, min_j = -1;
        for (int i = 0; i < m - 1; ++i) {
            int idx_i = active_indices[i];
            for (int j = i + 1; j < m; ++j) {
                int idx_j = active_indices[j];
                double Q = (m - 2) * D[idx_i][idx_j] - R[i] - R[j];
                if (Q < min_Q) {
                    min_Q = Q;
                    min_i = i;
                    min_j = j;
                }
            }
        }

        int idx_i = active_indices[min_i];
        int idx_j = active_indices[min_j];

        // Compute branch lengths
        double delta = (R[min_i] - R[min_j]) / (m - 2);
        double limb_length_i = 0.5 * D[idx_i][idx_j] + 0.5 * delta;
        double limb_length_j = 0.5 * D[idx_i][idx_j] - 0.5 * delta;

        // Create new node
        Node* new_node = new Node("", next_node_index++);
        new_node->left = nodes[idx_i];
        new_node->right = nodes[idx_j];
        new_node->branch_length_left = limb_length_i;
        new_node->branch_length_right = limb_length_j;

        // Update distances
        nodes.push_back(new_node);
        int new_idx = nodes.size() - 1;
        active_indices[min_i] = new_idx;  // Replace idx_i with new_idx
        active_indices.erase(active_indices.begin() + min_j);  // Remove idx_j

        // Update the distance matrix
        for (int k = 0; k < D.size(); ++k) {
            if (k == idx_i || k == idx_j)
                continue;
            D[k].push_back(0.0);
        }
        D.push_back(vector<double>(D.size() + 1, 0.0));

        for (int k = 0; k < D.size(); ++k) {
            if (k == new_idx)
                continue;
            int idx_k = k;
            double dik = D[idx_i][idx_k];
            double djk = D[idx_j][idx_k];
            double dikj = (dik + djk - D[idx_i][idx_j]) / 2.0;
            D[new_idx][idx_k] = dikj;
            D[idx_k][new_idx] = dikj;
        }

        // Mark old distances as invalid
        for (int k = 0; k < D.size(); ++k) {
            D[idx_i][k] = D[k][idx_i] = 0.0;
            D[idx_j][k] = D[k][idx_j] = 0.0;
        }

        nodes[idx_i] = new_node;
        nodes[idx_j] = nullptr;
    }

    // Join the last two nodes
    int idx_a = active_indices[0];
    int idx_b = active_indices[1];
    double limb_length_a = D[idx_a][idx_b] / 2.0;
    double limb_length_b = D[idx_a][idx_b] / 2.0;

    Node* root = new Node("", next_node_index++);
    root->left = nodes[idx_a];
    root->right = nodes[idx_b];
    root->branch_length_left = limb_length_a;
    root->branch_length_right = limb_length_b;

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

    // Perform Neighbor-Joining
    neighborJoining(taxa_names, distances);

    return 0;
}
