#!/usr/bin/env python3
import sys
import os

# Constants
FASTA_LINE_WIDTH = 60  # Standard line width for FASTA and CLUSTALformat
ALN_ID_WIDTH = 30  # Width for sequence IDs in CLUSTAL and STOCKHOLM formats


def detect_format(file_path):
    """Detects alignment format based on file extension.

    >>> detect_format("test.fasta")
    'fasta'
    >>> detect_format("test.fa")
    'fasta'
    >>> detect_format("test.sto")
    'sto'
    >>> detect_format("test.aln")
    'aln'
    >>> detect_format("test.clustal")
    'aln'
    >>> detect_format("TEST.FASTA")
    'fasta'
    """
    ext = os.path.splitext(file_path)[1].lower()
    format_map = {
        ".fa": "fasta",
        ".fasta": "fasta",
        ".sto": "sto",
        ".aln": "aln",
        ".clustal": "aln",
    }

    if ext not in format_map:
        raise ValueError("Unsupported format: {}".format(ext))
    return format_map[ext]


def read_alignment(file_path):
    """Reads an alignment file and returns a dictionary {seq_id: sequence}."""
    fmt = detect_format(file_path)

    try:
        with open(file_path) as f:
            if fmt == "fasta":
                return _read_fasta(f)
            elif fmt == "sto":
                return _read_stockholm(f)
            elif fmt == "aln":
                return _read_clustal(f)
    except IOError as e:
        raise IOError("Error reading file {}: {}".format(file_path, e))


def _read_fasta(f):
    """Read FASTA format efficiently using lists."""
    sequences = {}
    seq_id = None
    seq_parts = []

    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seq_id:
                sequences[seq_id] = "".join(seq_parts)
            seq_id = line[1:].split()[0]
            seq_parts = []
        else:
            seq_parts.append(line)

    if seq_id:
        sequences[seq_id] = "".join(seq_parts)

    return sequences


def _read_stockholm(f):
    """Read Stockholm format efficiently using lists."""
    seq_parts = {}  # dict of lists for efficient building

    for line in f:
        line = line.strip()
        if not line or line.startswith("#") or line == "//":
            continue

        parts = line.split(None, 1)  # Split on first whitespace only
        if len(parts) < 2:
            continue

        seq_id, seq_fragment = parts
        if seq_id not in seq_parts:
            seq_parts[seq_id] = []
        seq_parts[seq_id].append(seq_fragment.replace(".", "-"))

    # Join all parts into final sequences
    sequences = {}
    for seq_id, parts in seq_parts.items():
        sequences[seq_id] = "".join(parts)

    return sequences


def _read_clustal(f):
    """Read Clustal format efficiently using lists."""
    seq_parts = {}

    for line in f:
        line = line.rstrip()
        if (
            not line
            or line.startswith("CLUSTAL")
            or line.startswith(" ")
            or line.startswith("*")
        ):
            continue

        parts = line.split(None, 1)
        if len(parts) < 2:
            continue

        seq_id, seq_fragment = parts
        if seq_id not in seq_parts:
            seq_parts[seq_id] = []
        seq_parts[seq_id].append(seq_fragment)

    # Join all parts into final sequences
    sequences = {}
    for seq_id, parts in seq_parts.items():
        sequences[seq_id] = "".join(parts)

    return sequences


def validate_alignment(sequences):
    """Validates that all sequences have the same length.

    >>> validate_alignment({'seq1': 'ACGT', 'seq2': 'TGCA'})
    >>> validate_alignment({'seq1': 'ACGT', 'seq2': 'TGC'})
    Traceback (most recent call last):
        ...
    ValueError: Sequences have different lengths: {3, 4}
    >>> validate_alignment({})
    Traceback (most recent call last):
        ...
    ValueError: No sequences found
    """
    if not sequences:
        raise ValueError("No sequences found")

    lengths = set(len(seq) for seq in sequences.values())
    if len(lengths) > 1:
        raise ValueError("Sequences have different lengths: {}".format(lengths))


def write_alignment(sequences, file_path):
    """Writes sequences to file in the detected format."""
    fmt = detect_format(file_path)
    validate_alignment(sequences)

    try:
        with open(file_path, "w") as out:
            if fmt == "fasta":
                _write_fasta(sequences, out)
            elif fmt == "sto":
                _write_stockholm(sequences, out)
            elif fmt == "aln":
                _write_clustal(sequences, out)
    except IOError as e:
        raise IOError("Error writing file {}: {}".format(file_path, e))


def _write_fasta(sequences, out):
    """Write FASTA format."""
    for seq_id, seq in sequences.items():
        out.write(">{}\n".format(seq_id))
        for i in range(0, len(seq), FASTA_LINE_WIDTH):
            out.write("{}\n".format(seq[i : i + FASTA_LINE_WIDTH]))


def _write_stockholm(sequences, out):
    """Write Stockholm format."""
    out.write("# STOCKHOLM 1.0\n")
    for seq_id, seq in sequences.items():
        out.write("{} {}\n".format(seq_id.ljust(ALN_ID_WIDTH), seq))
    out.write("//\n")


def _write_clustal(sequences, out):
    """Write Clustal format."""
    out.write("CLUSTAL W alignment\n\n")

    if not sequences:
        return

    # Get alignment length from first sequence
    aln_len = len(list(sequences.values())[0])
    ids = list(sequences.keys())

    for start in range(0, aln_len, FASTA_LINE_WIDTH):
        for seq_id in ids:
            fragment = sequences[seq_id][start : start + FASTA_LINE_WIDTH]
            out.write("{} {}\n".format(seq_id.ljust(ALN_ID_WIDTH), fragment))
        out.write("\n")  # Single blank line between blocks


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: {} input_file output_file".format(sys.argv[0]))
        sys.exit(1)

    # Run doctests if --test flag is provided
    if "--test" in sys.argv:
        import doctest

        doctest.testmod()
        print("All doctests passed!")
        sys.exit(0)

    try:
        seqs = read_alignment(sys.argv[1])
        write_alignment(seqs, sys.argv[2])
        print("Converted {} -> {}".format(sys.argv[1], sys.argv[2]))
    except Exception as e:
        print("Error: {}".format(e), file=sys.stderr)
        sys.exit(1)
