import argparse


def nw(seq_i, seq_j, match, mismatch, gap):
    """
    Global alignment of 2 sequences using Needleman-Wunsch.

    >>> nw('FAT', 'FAST', 2, -1, -1)
    Optimal score: 5
    ('FA-T', 'FAST')
    >>> nw('THEFASTCAT', 'THERAT', 8, -8, -4)
    Optimal score: 16
    ('THE-FASTCAT', 'THER-A-T---')
    >>> nw('THERAT', 'THEFASTCAT', 8, -8, -4)
    Optimal score: 16
    ('THE-RA-T---', 'THEF-ASTCAT')
    """
    n_i, n_j = len(seq_i), len(seq_j)

    # Initialize scores and traceback matrices
    scores = [[0] * (n_j + 1) for _ in range(n_i + 1)]
    traceback = [[0] * (n_j + 1) for _ in range(n_i + 1)]

    # Initialize edges with gaps
    for i in range(1, n_i + 1):
        scores[i][0] = i * gap
        traceback[i][0] = 1  # deletion (gap in seq_j)

    for j in range(1, n_j + 1):
        scores[0][j] = j * gap
        traceback[0][j] = -1  # insertion (gap in seq_i)

    # Fill the matrices
    for i in range(1, n_i + 1):
        for j in range(1, n_j + 1):
            # Calculate scores for all possibilities
            match_score = match if seq_i[i - 1] == seq_j[j - 1] else mismatch
            diag = scores[i - 1][j - 1] + match_score  # match/mismatch
            left = scores[i][j - 1] + gap  # insertion (gap in seq_i)
            up = scores[i - 1][j] + gap  # deletion (gap in seq_j)

            # Choose the best score
            scores[i][j] = max(diag, left, up)

            # Set traceback pointer
            if diag > left and diag > up:
                traceback[i][j] = 0  # diagonal (match/mismatch)
            elif left > up:
                traceback[i][j] = -1  # left (insertion)
            else:
                traceback[i][j] = 1  # up (deletion)

    # Print optimal score
    print("Optimal score:", scores[-1][-1])

    # Traceback
    aln_i, aln_j = [], []
    i, j = n_i, n_j

    while i > 0 or j > 0:
        if traceback[i][j] == 0:  # diagonal
            i -= 1
            j -= 1
            aln_i.append(seq_i[i])
            aln_j.append(seq_j[j])
        elif traceback[i][j] == -1:  # left
            j -= 1
            aln_i.append("-")
            aln_j.append(seq_j[j])
        else:  # up (traceback[i][j] == 1)
            i -= 1
            aln_i.append(seq_i[i])
            aln_j.append("-")

    return "".join(reversed(aln_i)), "".join(reversed(aln_j))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Needleman-Wunsch global sequence alignment"
    )
    parser.add_argument("--seq1", required=True, help="First sequence")
    parser.add_argument("--seq2", required=True, help="Second sequence")
    parser.add_argument("--match", type=int, required=True, help="Match score")
    parser.add_argument("--mismatch", type=int, required=True, help="Mismatch penalty")
    parser.add_argument("--gap", type=int, required=True, help="Gap penalty")
    args = parser.parse_args()

    # a1, a2 = nw("FAT", "FAST", 2, -1, -1)
    # a1, a2 = nw("THEBIGCAT", "THERAT", 1, -1, -2)
    a1, a2 = nw(args.seq1, args.seq2, args.match, args.mismatch, args.gap)
    print(a1)
    print(a2)
    print()
