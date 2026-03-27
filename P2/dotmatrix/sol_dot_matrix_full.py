def dot_matrix(seq1, seq2, window=1, threshold=1):
    """
    >>> dot_matrix("FASTCAT", "FATRAT", 4, 2)
    [['o', ' ', ' ', ' ', ' ', ' ', ' '], [' ', 'o', 'o', ' ', ' ', ' ', ' '], [' ', ' ', 'o', 'o', ' ', ' ', ' '], [' ', ' ', ' ', 'o', 'o', ' ', ' '], [' ', ' ', ' ', ' ', ' ', 'o', ' '], [' ', ' ', ' ', ' ', ' ', ' ', 'o']]
    >>> dot_matrix("FASTCAT", "FATRAT", 2, 4)
    A threshold larger than the window does not make sense
    >>> dot_matrix("FASTCAT", "FATRAT")
    [['o', ' ', ' ', ' ', ' ', ' ', ' '], [' ', 'o', ' ', ' ', ' ', 'o', ' '], [' ', ' ', ' ', 'o', ' ', ' ', 'o'], [' ', ' ', ' ', ' ', ' ', ' ', ' '], [' ', 'o', ' ', ' ', ' ', 'o', ' '], [' ', ' ', ' ', 'o', ' ', ' ', 'o']]
    """
    if threshold > window:
        print("A threshold larger than the window does not make sense")
        return None

    dot_mat = [[" " for _ in seq1] for _ in seq2]

    for i in range(len(seq2) - window + 1):
        for j in range(len(seq1) - window + 1):
            matches = sum(seq1[j + k] == seq2[i + k] for k in range(window))
            if matches >= threshold:
                for k in range(window):
                    dot_mat[i + k][j + k] = "o"
    return dot_mat


if __name__ == "__main__":
    seq1 = "FASTCAT"
    seq2 = "FATRAT"
    dot_mat = dot_matrix(seq1, seq2, 4, 2)
    # print the matrix
    print(f" {seq1}")
    for char2, row in zip(seq2, dot_mat):
        print(char2 + "".join(row))
