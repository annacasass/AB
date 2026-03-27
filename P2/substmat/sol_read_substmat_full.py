def read_substmat_full(filename):
    """Reads a substitution matrix from a file."""

    with open(filename) as f:
        mat = {}
        aa_names = f.readline().strip().split()
        for line in f:
            aa_name, *scores = line.split()
            mat[aa_name] = {aa_names[i]: int(score) for i, score in enumerate(scores)}
    return mat


if __name__ == "__main__":
    mat = read_substmat_full("blosum62_full.mat")
    print(mat)
