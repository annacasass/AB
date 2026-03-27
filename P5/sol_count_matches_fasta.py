from Bio import SeqIO


def count_matches_fasta(fasta_file):
    """
    >>> count_matches_fasta("fastcats.fasta")
    5
    """
    records = SeqIO.parse(fasta_file, "fasta")
    seqs = [record.seq for record in records]
    # matches = 0
    #    for i in range(len(seqs[0])):  # all columns
    #        for j in range(len(seqs) - 1):  # all sequences
    #            match_pos = True  # the column contains the same residue in all seqs
    #            if seqs[j][i] != seqs[j + 1][i]:  # compare with next seq
    #                match_pos = False
    #                break  # if not a match we can go to the next column
    #        if match_pos:
    #            matches += 1  # one more column with matches in each sequence
    return sum(
        all(seqs[j][i] == seqs[j + 1][i] for j in range(len(seqs) - 1))
        for i in range(len(seqs[0]))
    )


if __name__ == "__main__":
    print(count_matches_fasta("fastcats.fasta"))
