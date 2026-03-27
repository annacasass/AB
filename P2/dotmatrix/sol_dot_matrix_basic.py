def dot_matrix(seq1, seq2):
    """
    >>> dot_matrix("FASTCAT", "FASTCAT")
    [['o', ' ', ' ', ' ', ' ', ' ', ' '], [' ', 'o', ' ', ' ', ' ', 'o', ' '], [' ', ' ', 'o', ' ', ' ', ' ', ' '], [' ', ' ', ' ', 'o', ' ', ' ', 'o'], [' ', ' ', ' ', ' ', 'o', ' ', ' '], [' ', 'o', ' ', ' ', ' ', 'o', ' '], [' ', ' ', ' ', 'o', ' ', ' ', 'o']]
    >>> dot_matrix("FASTCAT", "TACTSAF")
    [[' ', ' ', ' ', 'o', ' ', ' ', 'o'], [' ', 'o', ' ', ' ', ' ', 'o', ' '], [' ', ' ', ' ', ' ', 'o', ' ', ' '], [' ', ' ', ' ', 'o', ' ', ' ', 'o'], [' ', ' ', 'o', ' ', ' ', ' ', ' '], [' ', 'o', ' ', ' ', ' ', 'o', ' '], ['o', ' ', ' ', ' ', ' ', ' ', ' ']]
    """
    dot_mat = [[" " if char1 != char2 else "o" for char1 in seq1] for char2 in seq2]

    # dot_mat = [[" " for _ in range(len(seq1))] for _ in range(len(seq2))]
    # for i, char2 in enumerate(seq2):
    #    for j, char1 in enumerate(seq1):
    #        if char1 == char2:
    #            dot_mat[i][j] = "o"
    return dot_mat


if __name__ == "__main__":
    seq1 = "FASTCAT"
    seq2 = "FATRAT"
    dot_mat = dot_matrix(seq1, seq2)
    print(f" {seq1}")
    for char2, row in zip(seq2, dot_mat):
        print(char2 + "".join(row))
