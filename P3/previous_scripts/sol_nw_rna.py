def nw(seq_i, seq_j, match, mismatch, gep):
    """Returns the global sequence alignment of seq_i and seq_j
    by using the match, mismatch and gep parameters
    >>> nw('T', 'U', 1, -3, -1)
    ('T', 'U')
    >>> nw('T', 'A', 1, -3, -1)
    ('T-', '-A')
    >>> nw('ATCGATGCTATGCTAAATACGAT', 'UCGAUGCUAUCUAAAUAACGAU', 1, -3, -1)
    ('ATCGATGCTATGCTAAAT-ACGAT', '-UCGAUGCUAU-CUAAAUAACGAU')
    >>> nw('UCGAUGCUAUCUAAAUAACGAU', 'ATCGATGCTATGCTAAATACGAT', 1, -3, -1)
    ('-UCGAUGCUAU-CUAAAUAACGAU', 'ATCGATGCTATGCTAAAT-ACGAT')
    """
    # change to DNA alphabet
    dna_i = seq_i.replace("U", "T")
    dna_j = seq_j.replace("U", "T")

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
            score = match if dna_i[i - 1] == dna_j[j - 1] else mismatch
            subst = score_mat[i - 1][j - 1] + score
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

    return "".join(reversed(aln_i)), "".join(reversed(aln_j))
