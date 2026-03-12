from Bio.Align import substitution_matrices
from fasta2dict import fasta2dict

subst_mat = substitution_matrices.load("BLOSUM62")


def align_profiles(seqs, names_i, names_j, gap):
    """Makes a profile alignment of sequences with names names_i
    and names_j, whose sequences are in seqs, a dictionary of sequences
    with names (as keys)
    >>> align_profiles(fasta2dict('cats_and_rats.fasta'), ['FAST_CAT', 'FAT_CAT'], ['A_RAT', 'THE_RATS'], -2)
    (3.42, {'FAST_CAT': 'thefastcat', 'FAT_CAT': 'a--fa-tcat', 'A_RAT': 'a--ra-t---', 'THE_RATS': 'thera-t--s'})
    >>> align_profiles(fasta2dict('hmgb.fasta'), ['hmgt_mouse'], ['hmgl_trybr', 'hmgl_wheat', 'hmgb_chite'], -2)[0]
    9.33
    """
    n_i, n_j = len(seqs[names_i[0]]), len(seqs[names_j[0]])

    # Initialize scores and traceback matrices
    scores = [[0] * (n_j + 1) for _ in range(n_i + 1)]
    traceback = [[0] * (n_j + 1) for _ in range(n_i + 1)]

    # Initialize edges with gaps
    for i in range(1, n_i + 1):
        scores[i][0] = i * gap
        traceback[i][0] = 1  # up (gap in profile_j)

    for j in range(1, n_j + 1):
        scores[0][j] = j * gap
        traceback[0][j] = -1  # left (gap in profile_i)

    # Fill the matrices
    for i in range(1, n_i + 1):
        for j in range(1, n_j + 1):
            # Average substitution score across all sequence pairs, including gap penalties
            score, n_pairs = 0, 0
            for name_i in names_i:
                for name_j in names_j:
                    aa_i, aa_j = seqs[name_i][i - 1], seqs[name_j][j - 1]

                    # FILL IN THE SCORE CALCULATION

            score = score / n_pairs if n_pairs > 0 else 0

            diag = scores[i - 1][j - 1] + score  # match/mismatch
            left = scores[i][j - 1] + gap  # insertion (gap in profile_i)
            up = scores[i - 1][j] + gap  # deletion (gap in profile_j)

            scores[i][j] = max(diag, left, up)

            if diag >= left and diag >= up:
                traceback[i][j] = 0  # diagonal
            elif left >= up:
                traceback[i][j] = -1  # left
            else:
                traceback[i][j] = 1  # up

    # Traceback
    i, j = n_i, n_j
    aln = {name: [] for name in names_i + names_j}

    while i > 0 or j > 0:
        if traceback[i][j] == 0:  # diagonal
            i -= 1
            j -= 1
            for name_i in names_i:
                aln[name_i].append(seqs[name_i][i])
            for name_j in names_j:
                aln[name_j].append(seqs[name_j][j])
        elif traceback[i][j] == -1:  # left (gap in profile_i)
            j -= 1
            for name_i in names_i:
                aln[name_i].append("-")
            for name_j in names_j:
                aln[name_j].append(seqs[name_j][j])
        else:  # up (gap in profile_j)
            i -= 1
            for name_i in names_i:
                aln[name_i].append(seqs[name_i][i])
            for name_j in names_j:
                aln[name_j].append("-")

    # Write aligned sequences back into seqs
    for name, chars in aln.items():
        seqs[name] = "".join(reversed(chars))

    return round(scores[-1][-1], 2), seqs


if __name__ == "__main__":
    print("== cats ==")
    print(
        align_profiles(
            fasta2dict("cats_and_rats.fasta"),
            ["FAST_CAT", "FAT_CAT"],
            ["A_RAT", "THE_RATS"],
            -4,
        )
    )
    print("== hmgb ==")
    print(
        align_profiles(
            fasta2dict("hmgb.fasta"),
            ["hmgt_mouse"],
            ["hmgl_trybr", "hmgl_wheat", "hmgb_chite"],
            -4,
        )
    )
