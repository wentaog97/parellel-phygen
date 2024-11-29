#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Function to generate a random nucleotide, gap, or ambiguous base
char randomCharacter(double gapProb, double ambiguousProb) {
    double randValue = (double)rand() / RAND_MAX; // Random value between 0 and 1

    if (randValue < gapProb) {
        return '-'; // Gap character
    } else if (randValue < gapProb + ambiguousProb) {
        return 'N'; // Ambiguous base
    } else {
        char nucleotides[] = {'A', 'T', 'C', 'G'};
        return nucleotides[rand() % 4]; // Random nucleotide
    }
}

// Function to generate sequences
void generateSequences(int n, int m, double gapProb, double ambiguousProb) {
    printf("%d %d\n", n, m); // Print header in PHYLIP format
    for (int i = 0; i < n; i++) {
        // Print organism ID (e.g., Org1, Org2, ...)
        printf("Org%-6d", i + 1); 
        // Generate and print a sequence of length m
        for (int j = 0; j < m; j++) {
            printf("%c", randomCharacter(gapProb, ambiguousProb));
        }
        printf("\n"); // Newline after each sequence
    }
}

int main(int argc, char *argv[]) {
    if (argc != 5 && argc!= 3) {
        printf("Usage: %s n m gapProb ambiguousProb\n", argv[0]);
        printf("n = number of organisms\n");
        printf("m = length of gene sequence\n");
        printf("Optional:\n");
        printf("gapProb = probability of a gap (0-1)\n");
        printf("ambiguousProb = probability of an ambiguous base (0-1)\n");
        return 1;
    }

    // Parse input arguments
    int n = atoi(argv[1]); // Number of organisms
    int m = atoi(argv[2]); // Length of gene sequence

    double gapProb = 0;
    double ambiguousProb = 0;

    if(argc == 5){
        gapProb = atof(argv[3]); // Probability of a gap
        ambiguousProb = atof(argv[4]); // Probability of an ambiguous base
    }

    if (n <= 0 || m <= 0 || gapProb < 0 || gapProb > 1 || ambiguousProb < 0 || ambiguousProb > 1) {
        printf("Error: Invalid inputs. Ensure n, m > 0 and 0 <= gapProb, ambiguousProb <= 1.\n");
        return 1;
    }

    // Seed the random number generator
    srand(time(NULL));

    // Generate sequences
    generateSequences(n, m, gapProb, ambiguousProb);

    return 0;
}
