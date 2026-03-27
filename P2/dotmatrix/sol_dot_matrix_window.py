def dot_matrix(seq1, seq2, window=1):
    """
    >>> dot_matrix("FASTCAT", "FATRAT", 2)
    [['o', ' ', ' ', ' ', ' ', ' ', ' '], [' ', ' ', ' ', ' ', ' ', 'o', ' '], [' ', ' ', ' ', ' ', ' ', ' ', ' '], [' ', ' ', ' ', ' ', ' ', ' ', ' '], [' ', ' ', ' ', ' ', ' ', 'o', ' '], [' ', ' ', ' ', ' ', ' ', ' ', ' ']]
    >>> dot_matrix("FASTCAT", "FATRAT")
    [['o', ' ', ' ', ' ', ' ', ' ', ' '], [' ', 'o', ' ', ' ', ' ', 'o', ' '], [' ', ' ', ' ', 'o', ' ', ' ', 'o'], [' ', ' ', ' ', ' ', ' ', ' ', ' '], [' ', 'o', ' ', ' ', ' ', 'o', ' '], [' ', ' ', ' ', 'o', ' ', ' ', 'o']]
    >>> dot_matrix("FASTCAT", "FATRAT", 0)
    Window should be 1 or larger
    """
    if window < 1:
        print("Window should be 1 or larger")
        return None
    dot_mat = [[" " for _ in seq1] for _ in seq2]
    # create the matrix
    for i in range(len(seq2) - window + 1):
        for j in range(len(seq1) - window + 1):
            if all(seq1[j + k] == seq2[i + k] for k in range(window)):
                dot_mat[i][j] = "o"
    return dot_mat


if __name__ == "__main__":
    # print the matrix
    seq1 = "FASTCAT"
    seq2 = "FATRAT"
    dot_mat = dot_matrix(seq1, seq2, 2)
    print(f" {seq1}")
    for char2, row in zip(seq2, dot_mat):
        print(char2 + "".join(row))
