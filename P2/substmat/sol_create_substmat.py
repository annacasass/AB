from Bio import SeqIO
import math
import sys


def create_substmat(fasta_file, aa_alphabet="ACDEFGHIKLMNPQRSTVWY"):
    """Create a substitution matrix from a multiple sequence alignment in FASTA format."""

    # Read the sequences from the FASTA file
    records = list(SeqIO.parse(fasta_file, "fasta"))

    # Initialize amino acid frequencies and observed counts
    aa_freqs = {}
    obs = {}
    for aa in aa_alphabet:
        obs[aa] = {}
        aa_freqs[aa] = 0
        for aa2 in aa_alphabet:
            obs[aa][aa2] = 0

    # Count amino acid frequencies
    for record in records:
        for aa in record.seq:
            if aa == "-":
                continue
            if aa not in aa_alphabet:
                print(
                    "The alignment contains amino acids that are not in the provided alphabet"
                )
                return
            aa_freqs[aa] += 1

    total_count = sum(aa_freqs.values())
    for aa in aa_alphabet:
        aa_freqs[aa] /= total_count

    # Count observed substitutions
    total_obs = 0
    nseq = len(records)
    for i in range(nseq):  # i= index of seq1
        for j in range(i + 1, nseq):  # j= index of seq2
            for k in range(len(records[i].seq)):  # column index
                aa1 = records[i].seq[k]
                aa2 = records[j].seq[k]
                if aa1 != "-" and aa2 != "-":
                    obs[aa1][aa2] += 1
                    total_obs += 1

    subst_mat = {}
    for aa1 in aa_alphabet:
        subst_mat[aa1] = {}
        for aa2 in aa_alphabet:
            if aa1 == aa2:
                obs_freq = obs[aa1][aa2] / total_obs
                expected = aa_freqs[aa1] * aa_freqs[aa2]
            else:
                obs_freq = (obs[aa1][aa2] + obs[aa2][aa1]) / total_obs
                expected = 2 * aa_freqs[aa1] * aa_freqs[aa2]

            if obs_freq == 0 or aa_freqs[aa1] == 0 or aa_freqs[aa2] == 0:
                subst_mat[aa1][aa2] = -999
            else:
                subst_mat[aa1][aa2] = int(
                    round(math.log10(obs_freq / expected) * 10, 0)
                )

    aa_names = list(subst_mat.keys())
    result = "     " + "    ".join(aa_names) + "\n"
    for num_scores, aa1 in enumerate(aa_names, 1):
        line = aa1
        for i in range(num_scores):
            aa2 = aa_names[i]
            line += str(subst_mat[aa1][aa2]).rjust(5)
        result += line + "\n"
    print(result, end="")


if __name__ == "__main__":
    create_substmat("baba.fasta", "ABC")
