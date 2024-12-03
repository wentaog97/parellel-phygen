# parellel-phygen
Construct Phylogenetic Tree with parallel implementations.

The purpose of this program is calculating genetic distances between input organisms and then clustering them using UPGMA method.

Seqgen is used to generate testing data, hypothetical organisms with generated DNA sequences.

Use the output of Seqgen and put it into dnadist to calculate the distance matrix.

Use the distance matrix in UPGMA and get the newick tree. 

---

To Run the java program
1. Have input ready, the sample input file is "genes" 
2. Compile DNADist.java with "Javac DNADist.java"
3. Run DNADist with "Java DNADist genes > distMat"
4. Compile UPGMA.java with "Javac UPGMA.java"
5. Run UPGMA with "Java UPGMA distMat > result"
