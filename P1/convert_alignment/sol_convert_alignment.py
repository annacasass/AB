import sys
import os


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
                out.write("\n\n")  # blank line after block


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
