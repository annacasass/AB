import argparse

out_of_band = -1000  # very negative value for out-of-band cells


def nw_banded(seq_i, seq_j, match, mismatch, gap, k=None):
    """
    Global alignment of 2 sequences using a banded Needleman-Wunsch.

    Only cells dp[i][j] where |i - j| <= k are filled.  When k is None it
    defaults to max(len(seq_i), len(seq_j)) // 4, which works well in
    practice.  If the optimal alignment requires more than k gaps, k is
    widened to abs(len(seq_i) - len(seq_j)) to ensure the endpoint is
    reachable.

    >>> nw_banded('FAT', 'FAST', 2, -1, -1, k=1)
    Optimal score: 5
    ('FA-T', 'FAST')
    >>> nw_banded('THEBIGCAT', 'THERAT', 8, -8, -4, k=3)
    Optimal score: 20
    ('THEBIGCAT', 'THE---RAT')
    >>> nw_banded('THERAT', 'THEBIGCAT', 8, -8, -4, k=3)
    Optimal score: 20
    ('THE---RAT', 'THEBIGCAT')
    """

    n_i, n_j = len(seq_i), len(seq_j)

    # When k is None it defaults to max(len(seq_i), len(seq_j)) // 4,
    if k is None:
        k = max(1, max(n_i, n_j) // 4)

    # k is widened to abs(len(seq_i) - len(seq_j)) if necessary to ensure the endpoint is reachable.
    min_k = abs(n_i - n_j)
    if min_k > k:
        k = min_k

    scores = [[out_of_band] * (n_j + 1) for _ in range(n_i + 1)]
    traceback = [[0] * (n_j + 1) for _ in range(n_i + 1)]

    scores[0][0] = 0

    for i in range(1, min(n_i, k) + 1):
        scores[i][0] = i * gap
        traceback[i][0] = 1

    for j in range(1, min(n_j, k) + 1):
        scores[0][j] = j * gap
        traceback[0][j] = -1

    for i in range(1, n_i + 1):
        j_lo = max(1, i - k)
        j_hi = min(n_j, i + k)

        for j in range(j_lo, j_hi + 1):
            match_score = match if seq_i[i - 1] == seq_j[j - 1] else mismatch

            diag_val = scores[i - 1][j - 1] + match_score
            up_val = scores[i - 1][j] + gap
            left_val = scores[i][j - 1] + gap

            best = max(diag_val, left_val, up_val)
            scores[i][j] = best

            if diag_val >= left_val and diag_val >= up_val:
                traceback[i][j] = 0
            elif left_val >= up_val:
                traceback[i][j] = -1
            else:
                traceback[i][j] = 1

    print("Optimal score:", scores[n_i][n_j])

    # Uncomment to display all scores
    # for score in scores:
    #    print(score)

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
        else:  # up
            i -= 1
            aln_i.append(seq_i[i])
            aln_j.append("-")

    return "".join(reversed(aln_i)), "".join(reversed(aln_j))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Banded Needleman-Wunsch global sequence alignment"
    )
    parser.add_argument("--seq1", required=True, help="First sequence")
    parser.add_argument("--seq2", required=True, help="Second sequence")
    parser.add_argument("--match", type=int, required=True, help="Match score")
    parser.add_argument("--mismatch", type=int, required=True, help="Mismatch penalty")
    parser.add_argument("--gap", type=int, required=True, help="Gap penalty")
    parser.add_argument(
        "--band",
        type=int,
        default=None,
        help="Half-width k of the diagonal band (default: len/4)",
    )
    args = parser.parse_args()

    a1, a2 = nw_banded(
        args.seq1,
        args.seq2,
        args.match,
        args.mismatch,
        args.gap,
        k=args.band,
    )
    print(a1)
    print(a2)
    print()
