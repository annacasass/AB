from print_substmat import print_substmat
from read_substmat_half import read_substmat_half
import sys


def diff_substmat(subsmat_filename1, subsmat_filename2):
    """Compares two substitution matrices and prints the differences.

    Args:
        subsmat_filename1: First substitution matrix file
        subsmat_filename2: Second substitution matrix file
        format_type: "full" or "half" to specify the matrix format
        print_result: If True, prints the difference matrix; if False, only returns it

    Returns:
        Dictionary containing the differences between the two matrices
    """

    mat1 = read_substmat_half(subsmat_filename1)
    mat2 = read_substmat_half(subsmat_filename2)
    aa_names = list(mat1.keys())

    mat_difs = {}
    for aa1 in aa_names:
        mat_difs[aa1] = {}
        for aa2 in aa_names:
            if aa2 in mat1[aa1]:
                mat_difs[aa1][aa2] = mat2[aa1][aa2] - mat1[aa1][aa2]

    # Generate the formatted string for the substitution matrix
    result = "     " + "    ".join(aa_names) + "\n"
    for num_scores, aa1 in enumerate(aa_names, 1):
        line = [aa1]
        for i in range(num_scores):
            aa2 = aa_names[i]
            line.append(str(mat_difs[aa1][aa2]).rjust(2))
        result += "   ".join(line) + "\n"
    print(result, end="")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python diff_substmat.py <file1> <file2> [format_type]")
        print(
            "Example: python diff_substmat.py blosum62_half.mat blosum45_half.mat half"
        )
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    format_type = sys.argv[3] if len(sys.argv) > 3 else "half"

    # Directly prints the result
    diff_substmat(file1, file2, format_type)
