def nw3d(seq1, seq2, seq3, match=1, mismatch=-1, gap=-1):
    """
    Global alignment of 3 sequences using Needleman-Wunsch (3D).

    >>> nw3d("ABC", "ABC", "ABC")
    ('ABC', 'ABC', 'ABC')
    >>> nw3d("AB", "AC", "BC", 1, -1, -1)
    ('AB-', 'A-C', '-BC')
    """
    n1, n2, n3 = len(seq1), len(seq2), len(seq3)

    # Initialize score and traceback matrices
    score = [[[0] * (n3 + 1) for _ in range(n2 + 1)] for _ in range(n1 + 1)]
    traceback = [[[0] * (n3 + 1) for _ in range(n2 + 1)] for _ in range(n1 + 1)]

    # Initialize edges with gaps
    for i in range(1, n1 + 1):
        score[i][0][0] = i * gap
        traceback[i][0][0] = (i - 1, 0, 0)
    for j in range(1, n2 + 1):
        score[0][j][0] = j * gap
        traceback[0][j][0] = (0, j - 1, 0)
    for k in range(1, n3 + 1):
        score[0][0][k] = k * gap
        traceback[0][0][k] = (0, 0, k - 1)

    # Edges in planes
    for i in range(1, n1 + 1):
        for j in range(1, n2 + 1):
            score[i][j][0] = max(
                score[i - 1][j - 1][0]
                + (match if seq1[i - 1] == seq2[j - 1] else mismatch),
                score[i - 1][j][0] + gap,
                score[i][j - 1][0] + gap,
            )
            # Traceback: choose which gave max
            if score[i][j][0] == score[i - 1][j - 1][0] + (
                match if seq1[i - 1] == seq2[j - 1] else mismatch
            ):
                traceback[i][j][0] = (i - 1, j - 1, 0)
            elif score[i][j][0] == score[i - 1][j][0] + gap:
                traceback[i][j][0] = (i - 1, j, 0)
            else:
                traceback[i][j][0] = (i, j - 1, 0)
    for i in range(1, n1 + 1):
        for k in range(1, n3 + 1):
            score[i][0][k] = max(
                score[i - 1][0][k - 1]
                + (match if seq1[i - 1] == seq3[k - 1] else mismatch),
                score[i - 1][0][k] + gap,
                score[i][0][k - 1] + gap,
            )
            if score[i][0][k] == score[i - 1][0][k - 1] + (
                match if seq1[i - 1] == seq3[k - 1] else mismatch
            ):
                traceback[i][0][k] = (i - 1, 0, k - 1)
            elif score[i][0][k] == score[i - 1][0][k] + gap:
                traceback[i][0][k] = (i - 1, 0, k)
            else:
                traceback[i][0][k] = (i, 0, k - 1)
    for j in range(1, n2 + 1):
        for k in range(1, n3 + 1):
            score[0][j][k] = max(
                score[0][j - 1][k - 1]
                + (match if seq2[j - 1] == seq3[k - 1] else mismatch),
                score[0][j - 1][k] + gap,
                score[0][j][k - 1] + gap,
            )
            if score[0][j][k] == score[0][j - 1][k - 1] + (
                match if seq2[j - 1] == seq3[k - 1] else mismatch
            ):
                traceback[0][j][k] = (0, j - 1, k - 1)
            elif score[0][j][k] == score[0][j - 1][k] + gap:
                traceback[0][j][k] = (0, j - 1, k)
            else:
                traceback[0][j][k] = (0, j, k - 1)

    # Fill 3D matrix
    for i in range(1, n1 + 1):
        for j in range(1, n2 + 1):
            for k in range(1, n3 + 1):
                # 7 possibilities: match/mismatch in 3, gaps in 1/2/3, 2-gap combos
                scores = []
                pointers = []

                # match/mismatch all three
                s = score[i - 1][j - 1][k - 1] + (
                    match if seq1[i - 1] == seq2[j - 1] == seq3[k - 1] else mismatch
                )
                scores.append(s)
                pointers.append((i - 1, j - 1, k - 1))

                # gaps in one sequence
                scores.append(score[i - 1][j][k] + gap)
                pointers.append((i - 1, j, k))
                scores.append(score[i][j - 1][k] + gap)
                pointers.append((i, j - 1, k))
                scores.append(score[i][j][k - 1] + gap)
                pointers.append((i, j, k - 1))

                # gaps in two sequences
                scores.append(score[i - 1][j - 1][k] + gap)
                pointers.append((i - 1, j - 1, k))
                scores.append(score[i - 1][j][k - 1] + gap)
                pointers.append((i - 1, j, k - 1))
                scores.append(score[i][j - 1][k - 1] + gap)
                pointers.append((i, j - 1, k - 1))

                best = max(scores)
                score[i][j][k] = best
                traceback[i][j][k] = pointers[scores.index(best)]

    # Traceback
    aln1, aln2, aln3 = [], [], []
    i, j, k = n1, n2, n3
    while i > 0 or j > 0 or k > 0:
        prev = traceback[i][j][k]
        pi, pj, pk = prev
        if i > pi:
            aln1.append(seq1[i - 1])
        else:
            aln1.append("-")
        if j > pj:
            aln2.append(seq2[j - 1])
        else:
            aln2.append("-")
        if k > pk:
            aln3.append(seq3[k - 1])
        else:
            aln3.append("-")
        i, j, k = pi, pj, pk

    return "".join(reversed(aln1)), "".join(reversed(aln2)), "".join(reversed(aln3))


# Example
if __name__ == "__main__":
    a1, a2, a3 = nw3d("FASTCAT", "FATCAT", "THECAT", match=1, mismatch=-1, gap=-2)
    print(a1)
    print(a2)
    print(a3)
