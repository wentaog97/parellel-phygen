#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

struct Sequence {
    string name;
    string data;
};

double computeDistance(const string& seq1, const string& seq2) {
    int differences = 0;
    int total = 0;

    for (size_t i = 0; i < seq1.size(); ++i) {
        char base1 = toupper(seq1[i]);
        char base2 = toupper(seq2[i]);

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

    double p = static_cast<double>(differences) / total;

    // Check for the maximum p value for which the model is defined
    if (p >= 0.75) {
        return -1.0;  // Infinite distance
    }

    double distance = -0.75 * log(1.0 - (4.0 / 3.0) * p);
    return distance;
}

int main() {
    // Open the input file
    string input = "seq";
    ifstream inputFile(input);
    if (!inputFile.is_open()) {
        cerr << "Error: Could not open the input file 'seqs'." << endl;
        return 1;
    }

    // Read the number of sequences and their length
    int numSequences, sequenceLength;
    inputFile >> numSequences >> sequenceLength;

    vector<Sequence> sequences(numSequences);

    // Read the sequences
    for (int i = 0; i < numSequences; ++i) {
        inputFile >> sequences[i].name >> sequences[i].data;
    }
    inputFile.close();

    // Initialize the distance matrix
    vector< vector<double> > distanceMatrix(numSequences, vector<double>(numSequences, 0.0));

    // Compute distances between each pair of sequences
    for (int i = 0; i < numSequences; ++i) {
        for (int j = i + 1; j < numSequences; ++j) {
            double distance = computeDistance(sequences[i].data, sequences[j].data);
            distanceMatrix[i][j] = distance;
            distanceMatrix[j][i] = distance;
        }
    }

    // Output the distance matrix
    cout << fixed << setprecision(4);
    for (int i = 0; i < numSequences; ++i) {
        cout << setw(10) << sequences[i].name << " ";
        for (int j = 0; j < numSequences; ++j) {
            if (distanceMatrix[i][j] == -1.0) {
                cout << setw(8) << "N/A";  // No valid distance
            } else {
                cout << setw(8) << distanceMatrix[i][j];
            }
        }
        cout << endl;
    }

    return 0;
}
