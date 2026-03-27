#!/usr/bin/env python3
import sys
import os
from collections import Counter
from itertools import combinations

# BLOSUM62 simplified dictionary (scores for standard amino acids)
BLOSUM62 = {
    "A": {
        "A": 4,
        "R": -1,
        "N": -2,
        "D": -2,
        "C": 0,
        "Q": -1,
        "E": -1,
        "G": 0,
        "H": -2,
        "I": -1,
        "L": -1,
        "K": -1,
        "M": -1,
        "F": -2,
        "P": -1,
        "S": 1,
        "T": 0,
        "W": -3,
        "Y": -2,
        "V": 0,
    },
    "R": {
        "A": -1,
        "R": 5,
        "N": 0,
        "D": -2,
        "C": -3,
        "Q": 1,
        "E": 0,
        "G": -2,
        "H": 0,
        "I": -3,
        "L": -2,
        "K": 2,
        "M": -1,
        "F": -3,
        "P": -2,
        "S": -1,
        "T": -1,
        "W": -3,
        "Y": -2,
        "V": -3,
    },
    "N": {
        "A": -2,
        "R": 0,
        "N": 6,
        "D": 1,
        "C": -3,
        "Q": 0,
        "E": 0,
        "G": 0,
        "H": 1,
        "I": -3,
        "L": -3,
        "K": 0,
        "M": -2,
        "F": -3,
        "P": -2,
        "S": 1,
        "T": 0,
        "W": -4,
        "Y": -2,
        "V": -3,
    },
    "D": {
        "A": -2,
        "R": -2,
        "N": 1,
        "D": 6,
        "C": -3,
        "Q": 0,
        "E": 2,
        "G": -1,
        "H": -1,
        "I": -3,
        "L": -4,
        "K": -1,
        "M": -3,
        "F": -3,
        "P": -1,
        "S": 0,
        "T": -1,
        "W": -4,
        "Y": -3,
        "V": -3,
    },
    "C": {
        "A": 0,
        "R": -3,
        "N": -3,
        "D": -3,
        "C": 9,
        "Q": -3,
        "E": -4,
        "G": -3,
        "H": -3,
        "I": -1,
        "L": -1,
        "K": -3,
        "M": -1,
        "F": -2,
        "P": -3,
        "S": -1,
        "T": -1,
        "W": -2,
        "Y": -2,
        "V": -1,
    },
    "Q": {
        "A": -1,
        "R": 1,
        "N": 0,
        "D": 0,
        "C": -3,
        "Q": 5,
        "E": 2,
        "G": -2,
        "H": 0,
        "I": -3,
        "L": -2,
        "K": 1,
        "M": 0,
        "F": -3,
        "P": -1,
        "S": 0,
        "T": -1,
        "W": -2,
        "Y": -1,
        "V": -2,
    },
    "E": {
        "A": -1,
        "R": 0,
        "N": 0,
        "D": 2,
        "C": -4,
        "Q": 2,
        "E": 5,
        "G": -2,
        "H": 0,
        "I": -3,
        "L": -3,
        "K": 1,
        "M": -2,
        "F": -3,
        "P": -1,
        "S": 0,
        "T": -1,
        "W": -3,
        "Y": -2,
        "V": -2,
    },
    "G": {
        "A": 0,
        "R": -2,
        "N": 0,
        "D": -1,
        "C": -3,
        "Q": -2,
        "E": -2,
        "G": 6,
        "H": -2,
        "I": -4,
        "L": -4,
        "K": -2,
        "M": -3,
        "F": -3,
        "P": -2,
        "S": 0,
        "T": -2,
        "W": -2,
        "Y": -3,
        "V": -3,
    },
    "H": {
        "A": -2,
        "R": 0,
        "N": 1,
        "D": -1,
        "C": -3,
        "Q": 0,
        "E": 0,
        "G": -2,
        "H": 8,
        "I": -3,
        "L": -3,
        "K": -1,
        "M": -2,
        "F": -1,
        "P": -2,
        "S": -1,
        "T": -2,
        "W": -2,
        "Y": 2,
        "V": -3,
    },
    "I": {
        "A": -1,
        "R": -3,
        "N": -3,
        "D": -3,
        "C": -1,
        "Q": -3,
        "E": -3,
        "G": -4,
        "H": -3,
        "I": 4,
        "L": 2,
        "K": -3,
        "M": 1,
        "F": 0,
        "P": -3,
        "S": -2,
        "T": -1,
        "W": -3,
        "Y": -1,
        "V": 3,
    },
    "L": {
        "A": -1,
        "R": -2,
        "N": -3,
        "D": -4,
        "C": -1,
        "Q": -2,
        "E": -3,
        "G": -4,
        "H": -3,
        "I": 2,
        "L": 4,
        "K": -2,
        "M": 2,
        "F": 0,
        "P": -3,
        "S": -2,
        "T": -1,
        "W": -2,
        "Y": -1,
        "V": 1,
    },
    "K": {
        "A": -1,
        "R": 2,
        "N": 0,
        "D": -1,
        "C": -3,
        "Q": 1,
        "E": 1,
        "G": -2,
        "H": -1,
        "I": -3,
        "L": -2,
        "K": 5,
        "M": -1,
        "F": -3,
        "P": -1,
        "S": 0,
        "T": -1,
        "W": -3,
        "Y": -2,
        "V": -2,
    },
    "M": {
        "A": -1,
        "R": -1,
        "N": -2,
        "D": -3,
        "C": -1,
        "Q": 0,
        "E": -2,
        "G": -3,
        "H": -2,
        "I": 1,
        "L": 2,
        "K": -1,
        "M": 5,
        "F": 0,
        "P": -2,
        "S": -1,
        "T": -1,
        "W": -1,
        "Y": -1,
        "V": 1,
    },
    "F": {
        "A": -2,
        "R": -3,
        "N": -3,
        "D": -3,
        "C": -2,
        "Q": -3,
        "E": -3,
        "G": -3,
        "H": -1,
        "I": 0,
        "L": 0,
        "K": -3,
        "M": 0,
        "F": 6,
        "P": -4,
        "S": -2,
        "T": -2,
        "W": 1,
        "Y": 3,
        "V": -1,
    },
    "P": {
        "A": -1,
        "R": -2,
        "N": -2,
        "D": -1,
        "C": -3,
        "Q": -1,
        "E": -1,
        "G": -2,
        "H": -2,
        "I": -3,
        "L": -3,
        "K": -1,
        "M": -2,
        "F": -4,
        "P": 7,
        "S": -1,
        "T": -1,
        "W": -4,
        "Y": -3,
        "V": -2,
    },
    "S": {
        "A": 1,
        "R": -1,
        "N": 1,
        "D": 0,
        "C": -1,
        "Q": 0,
        "E": 0,
        "G": 0,
        "H": -1,
        "I": -2,
        "L": -2,
        "K": 0,
        "M": -1,
        "F": -2,
        "P": -1,
        "S": 4,
        "T": 1,
        "W": -3,
        "Y": -2,
        "V": -2,
    },
    "T": {
        "A": 0,
        "R": -1,
        "N": 0,
        "D": -1,
        "C": -1,
        "Q": -1,
        "E": -1,
        "G": -2,
        "H": -2,
        "I": -1,
        "L": -1,
        "K": -1,
        "M": -1,
        "F": -2,
        "P": -1,
        "S": 1,
        "T": 5,
        "W": -2,
        "Y": -2,
        "V": 0,
    },
    "W": {
        "A": -3,
        "R": -3,
        "N": -4,
        "D": -4,
        "C": -2,
        "Q": -2,
        "E": -3,
        "G": -2,
        "H": -2,
        "I": -3,
        "L": -2,
        "K": -3,
        "M": -1,
        "F": 1,
        "P": -4,
        "S": -3,
        "T": -2,
        "W": 11,
        "Y": 2,
        "V": -3,
    },
    "Y": {
        "A": -2,
        "R": -2,
        "N": -2,
        "D": -3,
        "C": -2,
        "Q": -1,
        "E": -2,
        "G": -3,
        "H": 2,
        "I": -1,
        "L": -1,
        "K": -2,
        "M": -1,
        "F": 3,
        "P": -3,
        "S": -2,
        "T": -2,
        "W": 2,
        "Y": 7,
        "V": -1,
    },
    "V": {
        "A": 0,
        "R": -3,
        "N": -3,
        "D": -3,
        "C": -1,
        "Q": -2,
        "E": -2,
        "G": -3,
        "H": -3,
        "I": 3,
        "L": 1,
        "K": -2,
        "M": 1,
        "F": -1,
        "P": -2,
        "S": -2,
        "T": 0,
        "W": -3,
        "Y": -1,
        "V": 4,
    },
}


def detect_format(file_path):
    """Detects alignment format based on file extension."""
    ext = os.path.splitext(file_path)[1].lower()
    if ext in [".fa", ".fasta"]:
        return "fasta"
    elif ext == ".sto":
        return "sto"
    elif ext in [".aln", ".clustal"]:
        return "aln"
    else:
        raise ValueError(f"Cannot detect format from extension: {ext}")


def read_alignment(file_path):
    """
    Reads an alignment file and returns a dictionary {seq_id: sequence}.
    Format is detected automatically.
    """
    fmt = detect_format(file_path)
    sequences = {}

    if fmt == "fasta":
        with open(file_path) as f:
            seq_id = None
            seq = []
            for line in f:
                line = line.strip()
                if not line:
                    continue  # skip empty lines
                if line.startswith(">"):
                    if seq_id:  # save previous sequence
                        sequences[seq_id] = "".join(seq)
                    seq_id = line[1:].split()[0]
                    seq = []
                else:  # append sequence lines
                    seq.append(line)
            if seq_id:
                sequences[seq_id] = "".join(seq)

    elif fmt == "sto":
        with open(file_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or line == "//":
                    continue  # skip comments and end marker
                parts = line.split()
                if len(parts) < 2:
                    continue  # skip malformed lines
                seq_id, seq_part = parts[0], parts[1]
                sequences[seq_id] = (
                    f"{sequences.get(seq_id, '')}{seq_part.replace('.', '-')}"
                )

    elif fmt == "aln":
        with open(file_path) as f:
            for line in f:
                line = line.rstrip()
                if not line or line.startswith("CLUSTAL") or line.startswith(" "):
                    continue  # skip headers and empty lines
                parts = line.split()
                if len(parts) < 2:
                    continue  # skip malformed lines
                seq_id, seq_part = parts[0], parts[1]
                sequences[seq_id] = sequences.get(seq_id, "") + seq_part

    return sequences


def make_consensus_line(seq_block, thresh_high=3, thresh_low=1):
    """
    Generates the consensus line for a Clustal block using BLOSUM62 scores.
    seq_block: list of equal-length sequence fragments
    thresh_high: threshold for ':'
    thresh_low: threshold for '.'
    """
    consensus = []
    aln_len = len(seq_block[0])
    for i in range(aln_len):
        column = [seq[i] for seq in seq_block]
        if "-" in column:
            consensus.append(" ")
            continue
        if all(res == column[0] for res in column):
            consensus.append("*")
            continue
        # calculate average BLOSUM62 score
        scores = []
        for a, b in combinations(column, 2):
            try:
                scores.append(BLOSUM62[a][b])
            except KeyError:
                scores.append(0)  # unknown residue
        avg_score = sum(scores) / len(scores)
        if avg_score >= thresh_high:
            consensus.append(":")
        elif avg_score >= thresh_low:
            consensus.append(".")
        else:
            consensus.append(" ")
    return " " * 31 + "".join(consensus)


def write_alignment(sequences, file_path):
    """Writes sequences to file in fasta/sto/aln with 60 residues per line."""
    fmt = detect_format(file_path)
    BLOCK_SIZE = 60

    with open(file_path, "w") as out:
        if fmt == "fasta":
            for seq_id, seq in sequences.items():
                out.write(f">{seq_id}\n")
                for i in range(0, len(seq), BLOCK_SIZE):
                    out.write(f"{seq[i:i+BLOCK_SIZE]}\n")

        elif fmt == "sto":
            out.write("# STOCKHOLM 1.0\n")
            for seq_id, seq in sequences.items():
                out.write(f"{seq_id.ljust(30)} {seq}\n")
            out.write("//\n")

        elif fmt == "aln":
            out.write("CLUSTAL W alignment\n\n")  # only once at the top
            aln_len = max(len(seq) for seq in sequences.values())
            ids = list(sequences.keys())
            for start in range(0, aln_len, BLOCK_SIZE):
                # get the current block of sequences
                seq_block = [
                    sequences[seq_id][start : start + BLOCK_SIZE] for seq_id in ids
                ]
                # write each sequence line
                for seq_id, fragment in zip(ids, seq_block):
                    out.write(f"{seq_id.ljust(30)} {fragment}\n")
                # write consensus line
                out.write(
                    make_consensus_line(seq_block) + "\n\n"
                )  # blank line after block


# Example usage
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input_file output_file")
        sys.exit(1)

    in_file = sys.argv[1]
    out_file = sys.argv[2]

    seqs = read_alignment(in_file)
    write_alignment(seqs, out_file)
    print(f"Converted {in_file} -> {out_file}")
