def read_substmat_half(filename):
    """Reads a substitution matrix from a file."""

    with open(filename) as f:
        mat = {}
        lines = f.readlines()
        aa_names = lines[-1].strip().split()
        for line in lines[:-1]:
            aa_name1, *scores = line.split()
            scores = map(int, scores)
            mat[aa_name1] = {}
            for aa_name2, score in zip(aa_names, scores):
                mat[aa_name1][aa_name2] = score
                mat[aa_name2][aa_name1] = score
    return mat
