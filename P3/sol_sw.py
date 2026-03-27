import argparse
from Bio.Align import substitution_matrices

subst_mat = substitution_matrices.load("BLOSUM62")


def sw(seq_i, seq_j, gap):
    """
    Local alignment of 2 sequences using Smith-Waterman algorithm
    with BLOSUM62 substitution matrix and a specified gap penalty.

    >>> sw('THEFASTCAT', 'THERAT', -4)
    Optimal score: 20.0
    ('THEFAS', 'THERAT')
    >>> sw('FAT', 'THEFASTCAT', -4)
    Optimal score: 11.0
    ('FAT', 'FAS')
    >>> sw('THECATISFAT', 'AFASTCAT', -4)
    Optimal score: 18.0
    ('CAT', 'CAT')
    """

    n_i, n_j = len(seq_i), len(seq_j)

    # Initialize scores and traceback matrices
    scores = [[0] * (n_j + 1) for _ in range(n_i + 1)]
    traceback = [[0] * (n_j + 1) for _ in range(n_i + 1)]

    # Initialize edges with gaps
    for i in range(1, n_i + 1):
        traceback[i][0] = 2  # deletion (gap in seq_j)

    for j in range(1, n_j + 1):
        traceback[0][j] = 2  # insertion (gap in seq_i)

    # Fill the matrices
    best_score = 0
    for i in range(1, n_i + 1):
        for j in range(1, n_j + 1):
            # Calculate scores for all possibilities
            match_score = subst_mat[seq_i[i - 1], seq_j[j - 1]]
            diag = scores[i - 1][j - 1] + match_score  # match/mismatch
            left = scores[i][j - 1] + gap  # insertion (gap in seq_i)
            up = scores[i - 1][j] + gap  # deletion (gap in seq_j)

            # Choose the best score
            scores[i][j] = max(0, diag, left, up)

            # Set traceback pointer
            if diag > left and diag > up and diag > 0:
                traceback[i][j] = 0  # diagonal (match/mismatch)
            elif left > up and left > 0:
                traceback[i][j] = -1  # left (insertion)
            elif up > 0:
                traceback[i][j] = 1  # up (deletion)
            else:
                traceback[i][j] = 2

            if scores[i][j] > best_score:
                best_score = scores[i][j]
                best_score_cell = (i, j)
    # Print optimal score
    print("Optimal score:", best_score)

    # Traceback
    aln_i, aln_j = [], []
    i, j = best_score_cell

    while scores[i][j] != 0:
        if traceback[i][j] == 0:  # diagonal
            i -= 1
            j -= 1
            aln_i.append(seq_i[i])
            aln_j.append(seq_j[j])
        elif traceback[i][j] == -1:  # left
            j -= 1
            aln_i.append("-")
            aln_j.append(seq_j[j])
        elif traceback[i][j] == 1:  # up
            i -= 1
            aln_i.append(seq_i[i])
            aln_j.append("-")

    return "".join(reversed(aln_i)), "".join(reversed(aln_j))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Smith-Waterman local sequence alignment"
    )
    parser.add_argument("--seq1", required=True, help="First sequence")
    parser.add_argument("--seq2", required=True, help="Second sequence")
    parser.add_argument("--gap", type=int, required=True, help="Gap penalty")
    args = parser.parse_args()

    # a1, a2 = nw("FAT", "FAST", 2, -1, -1)
    # a1, a2 = nw("THEBIGCAT", "THERAT", 1, -1, -2)
    a1, a2 = sw(args.seq1, args.seq2, args.gap)
    print(a1)
    print(a2)
    print()
