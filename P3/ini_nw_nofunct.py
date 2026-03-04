seq_i = "FAST"
seq_j = "FAT"
match = 1
mismatch = -1
gap = -2

# Store the lengths of the sequences
n_i, n_j = len(seq_i), len(seq_j)

# Initialize scores and traceback matrices
scores = [[0] * (n_j+1) for _ in range(n_i +1)]
traceback  = [[0] * (n_j+1) for _ in range(n_i +1)]

# Initialize edges with gaps
for i in range(1, n_i +1):
    scores[i][0] = i * gap
    traceback[i][0] = 1

for j in range(1,n_j +1):
    scores[0][j] = j * gap
    traceback[0][j] = -1

# Fill the matrices
for i in range(1, n_i +1):
    for j in range(,1 n_j +1):
    
        # Calculate scores for all possibilities
        subst_score = match if seq_i[i-1] == seq_j[j-1] else mismatch
        diag = scores[i-1][j-1] + subst_score
        left = scores[i][j-1] + gap
        up = scores[i-1][j] + gap
        # Choose the best scores
        scores[i][j] = max(diag, left, up)

        # Set traceback pointer
        if diag > left and diag > up:
            traceback[i][j] = 0  # diagonal (match/mismatch)
        elif left > up:
            traceback[i][j] = -1  # left (insertion)
        else:
            traceback[i][j] = 1  # up (deletion)


# Print scores matrix (for debugging)
header = 8 * " " + " ".join(f"{c:>3}" for c in seq_j)
print(header)
for aa, row in zip(f" {seq_i}", scores):
    print(f"{aa:>3}", " ".join(f"{v:>3}" for v in row))

# Print traceback matrix (for debugging))
header = 8 * " " + " ".join(f"{c:>3}" for c in seq_j)
print("\n" + header)
for aa, row in zip(f" {seq_i}", traceback):
    print(f"{aa:>3}", " ".join(f"{v:>3}" for v in row))

# Print optimal score
print("Optimal score:  __")  # fill in

# Prepare for traceback
aln_i, aln_j = [], []
i, j = n_i, n_j  # fill in; this is the starting point for traceback

# Traceback
while i > 0 or j > 0:  # fill in
    if traceback[i][j] == 0:  # diagonal
        i -= 1
        j -= 1
        aln_i.append(seq_i[i])
        aln_j.append(seq_j[j])
    elif traceback[i][j] == -1:  # left
        j -= 1
        aln_i.append('-')
        aln_j.append(seq_j[j])       
    else:  # up (traceback[i][j] == 1)
        i -= 1
        aln_i.append(seq_i[i])
        aln_j.append('-')

# Print the alignment
# a few lines here