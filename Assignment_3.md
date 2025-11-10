# Assignment 3



### <center><u>FASTA</u></center>
### <u>GOAL </u>
Implement a simplified FASTA like algorithm to find the the best matching sequence,

### <u>Data</u>
**1)** A fasta file containing the query sequence.
**2)** A multifasta file containing multiple possible matches of the query sequence.

### <u>Approach</u>
**1)** index the query and the database sequences to build a table containing a list of kmers and their corresponding index in the sequence.
**2)** Find overlaps between the two that correspond to diagonals of the matrix query x sequence. The highest scoring diagonals are the diagonals with the most overlaps(matches).
**3)** Diagonals with gap smaller than a predefined distance g can be merged to form a bigger sequence.

### <u>Implementantion</u>


<u>Load the libraries</u>
```python
from tqdm import tqdm as tqdm
import os
```


<u>Define a function to load the multifasta file as a dictionary</u>
```python
def read_data(path): #reads a textfile with multiple genes in fasta format and returns a dictionary withentries genes and their sequences
    file = open(path)
    gene_dic = {}
    genes = file.read()
    genes = genes.split( '\n')
    for gene_name in range(0,len(genes) -1,2):
        gene_dic[genes[gene_name]] = genes[gene_name + 1]

    return(gene_dic)
```

        
<u>Define a function that builds a kmmer index out of a given sequence.</u>
```python
def indexKmers(sequence, k):
    positions = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        if kmer not in positions:
            positions[kmer] = []
        positions[kmer].append(i)
    return positions
```

        
<u>Load the data</u>
```python
working_directory = os.getcwd()
Q = read_data(f'{working_directory}/data_for_assignments/query.fa')
S = read_data(f'{working_directory}/data_for_assignments/all_yeast_genes_minplus1k.fa')
```

        
<u>Define the parameters of the algorithm and run the algorithm.</u>
```python
threshold = 20
g = 3
k = 5

query_kmers = indexKmers(list(Q.values())[0], k)
max_score = 0


for seq_name, seq in tqdm(S.items()):
    sequence_kmers = indexKmers(seq,k)

    #Score each diagonal based on the number of the matching kmers
    diagonal_dic = {}
    for kmer, indices in sequence_kmers.items():
        if kmer in query_kmers:
            for index in indices:
                for j in query_kmers[kmer]:
                    diagonal_index = (index - j)
                    # increment the score for that diagonal
                    diagonal_dic[diagonal_index] = diagonal_dic.get(diagonal_index, 0) + 1


    #Filter out diagonals with score less than the threshold
    filtered_diagonal_dic = {diagonal : score for diagonal, score in diagonal_dic.items() if score >= threshold}



    #Sort the diagonal based on their index to make the merging more efficient
    sorted_diagonals = sorted(filtered_diagonal_dic.items())  # list of (diagonal_index, score)

    #Merge the diagonals based on the gap(g) parameter
    if sorted_diagonals:
        merged_diagonals = []
        current_diag, current_score = sorted_diagonals[0]

        for diag, score in sorted_diagonals[1:]:
            if diag - current_diag <= g:
                # If the two diagonals are closer than the gap parameter then merge them
                current_score += score
            else:
                merged_diagonals.append((current_diag, current_score))
                current_diag, current_score = diag, score

        merged_diagonals.append((current_diag, current_score))  # append the last merged diagonal
        best_diagonal = max(merged_diagonals, key=lambda x: x[1])
        if best_diagonal[1] > max_score:
            max_score = best_diagonal[1]
            best_sequence = seq_name

print(best_sequence)
```

<br>
<br>

### Plot the query sequence against the best matching sequence.


<br><br>
<div style="text-align: center;">
    <img src="figures/fasta_array.png" alt="Nonmers plot" width="700">

</div>
<br><br>
