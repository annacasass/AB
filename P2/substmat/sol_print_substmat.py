from read_substmat_full import read_substmat_full
from read_substmat_half import read_substmat_half


def print_substmat(filename, format_type="full"):
    """
    Prints a substitution matrix.

    Args:
        filename: Path to the substitution matrix file
        format_type: "full" or "half" to specify the matrix format
    """
    # Select the appropriate reading function based on the format type
    if format_type == "full":
        smat = read_substmat_full(filename)
    elif format_type == "half":
        smat = read_substmat_half(filename)
    else:
        raise ValueError("format_type must be 'full' or 'half'")

    # Generate the formatted string for the substitution matrix
    aa_names = list(smat.keys())
    result = "     " + "    ".join(aa_names) + "\n"
    for num_scores, aa1 in enumerate(aa_names, 1):
        line = [aa1]
        for i in range(num_scores):
            aa2 = aa_names[i]
            line.append(str(smat[aa1][aa2]).rjust(2))
        result += "   ".join(line) + "\n"
    print(result, end="")


if __name__ == "__main__":
    print_substmat("blosum62_full.mat", "full")
    print("\n")
    print_substmat("blosum45_half.mat", "half")
