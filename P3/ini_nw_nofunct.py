seq_i = "FAST"
seq_j = "FAT"
match = 1
mismatch = -1
gap = -2

# Store the lengths of the sequences
n_i, n_j = len(seq_i), len(seq_j)

# Initialize scores and traceback matrices
scores = __  # fill in
traceback = __  # fill in

# Initialize edges with gaps
# a few lines here

# Fill the matrices
# a few lines here

        # Calculate scores for all possibilities
        diag = __  # fill in
        left = __  # fill in  
        up = __   # fill in
        # Choose the best scores
        scores[i][j] = max(diag, left, up)

        # Set traceback pointer
        if diag > left and diag > up:
            traceback[i][j] = 0  # diagonal (match/mismatch)
            # a few lines here
        elif left > up:
            traceback[i][j] = -1  # left (insertion)
            # a few lines here
        else:
            traceback[i][j] = 1  # up (deletion)
            # a few lines here

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
i, j = __, __  # fill in; this is the starting point for traceback

# Traceback
while i __ or j __:  # fill in
    if traceback[i][j] == 0:  # diagonal
            i -= 1
            j -= 1
            aln_i.append(seq_i[i])
            aln_j.append(seq_j[j])
    elif traceback[i][j] == -1:  # left
        # a few lines here
    else:  # up (traceback[i][j] == 1)
        # a few lines here

# Print the alignment
# a few lines here