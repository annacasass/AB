from Bio.Align import substitution_matrices


def nw(seq_i, seq_j, gep):
    """Returns the global sequence alignment of seq_i and seq_j
    by using the match, mismatch and gep parameters
    nw('FAT', 'FAST', -1)
    (14.0, 'FA-T', 'FAST')
    >>> nw('THEFASTCAT', 'THERAT', -4)
    (10.0, 'THEFASTCAT', 'THE---R-AT')
    >>> nw('THERAT', 'THEFASTCAT', -4)
    (10.0, 'THE---R-AT', 'THEFASTCAT')
    >>> nw('T', 'FAST', -1)
    (2.0, '---T', 'FAST')
    """
    subst_mat = substitution_matrices.load("BLOSUM62")

    # We initialize the matrices
    score_mat = [[0 for _ in range(len(seq_j) + 1)] for _ in range(len(seq_i) + 1)]
    traceback = [[0 for _ in range(len(seq_j) + 1)] for _ in range(len(seq_i) + 1)]

    for i in range(1, len(seq_i) + 1):
        score_mat[i][0] = i * gep
        traceback[i][0] = 1

    for j in range(1, len(seq_j) + 1):
        score_mat[0][j] = j * gep
        traceback[0][j] = -1

    # We fill the matrices
    for i in range(1, len(seq_i) + 1):
        for j in range(1, len(seq_j) + 1):
            subst = score_mat[i - 1][j - 1] + subst_mat[seq_i[i - 1]][seq_j[j - 1]]
            inser = score_mat[i][j - 1] + gep
            delet = score_mat[i - 1][j] + gep

            if subst >= inser and subst >= delet:
                score_mat[i][j] = subst
                traceback[i][j] = 0
            elif inser >= delet:
                score_mat[i][j] = inser
                traceback[i][j] = -1
            else:
                score_mat[i][j] = delet
                traceback[i][j] = 1

    # We do the traceback
    aln_i = []
    aln_j = []

    i = len(seq_i)
    j = len(seq_j)
    while i != 0 or j != 0:
        if traceback[i][j] == 0:
            i -= 1
            j -= 1
            aln_i.append(seq_i[i])
            aln_j.append(seq_j[j])
        elif traceback[i][j] == -1:
            j -= 1
            aln_i.append("-")
            aln_j.append(seq_j[j])
        else:
            i -= 1
            aln_i.append(seq_i[i])
            aln_j.append("-")

    return score_mat[-1][-1], "".join(reversed(aln_i)), "".join(reversed(aln_j))
