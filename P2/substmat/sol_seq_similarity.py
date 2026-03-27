from Bio.Align import substitution_matrices

# Load the matrix once at module level instead of every function call
BLOSUM62 = substitution_matrices.load("BLOSUM62")


def seq_similarity(seq1, seq2):
    """Computes the similarity between two sequences as a percentage based on the BLOSUM62 matrix.
    >>> seq_similarity("FASTCAT", "FATCAT")

    >>> seq_similarity("FASTCAT", "FASTCAT")
    100.0
    >>> seq_similarity("FASTCAT", "FASTRAT")
    85.7
    >>> seq_similarity("-FASTCAT", "-FASTRAT")
    85.7
    >>> seq_similarity("FASTCAT", "FA-TCAT")
    100.0
    >>> seq_similarity("FASTCAT", "FAT-CAT")
    100.0
    >>> seq_similarity("AFASTCAT", "-FASTRAT")
    85.7
    >>> seq_similarity("FASTCAT", "AAAAAAA")
    42.9
    >>> seq_similarity("FASTCAT", "AFAAAFA")
    14.3
    """
    if len(seq1) != len(seq2):
        return None

    # Use sum with generator expression for more efficient counting
    matches = sum(
        1
        for char1, char2 in zip(seq1, seq2)
        if char1 != "-" and char2 != "-" and BLOSUM62[char1][char2] > 0
    )

    # Count non-gap positions
    tot_len = sum(
        1 for char1, char2 in zip(seq1, seq2) if char1 != "-" and char2 != "-"
    )

    if tot_len == 0:
        return 0.0

    return round((matches / tot_len) * 100, 1)


if __name__ == "__main__":
    seq1 = "-FASTCAT"
    seq2 = "-FASTRAT"
    print(seq1)
    print(seq2)
    print(seq_similarity(seq1, seq2))
    85.7
