from Bio import SeqIO


def count_matches_stars(fasta_file):
    """
    >>> count_matches_stars("fastcats.fasta")
    '       ** * **'
    """
    records = SeqIO.parse(fasta_file, "fasta")
    seqs = [record.seq for record in records]
    # for seq in seqs:
    #     print(seq)
    """
    matches = 0
    stars = []  # we will add "*" or "-" for each column
    for i in range(len(seqs[0])):
        for j in range(len(seqs) - 1):
            match_pos = True
            if seqs[j][i] != seqs[j + 1][i]:
                match_pos = False
                break
        if match_pos:
            matches += 1
            stars.append("*")
        else:
            stars.append(" ")
    return "".join(stars)  # convert list to str of stars
    """

    return "".join(
        [
            (
                "*"
                if all(seqs[j][i] == seqs[j + 1][i] for j in range(len(seqs) - 1))
                else " "
            )
            for i in range(len(seqs[0]))
        ]
    )


if __name__ == "__main__":
    print(count_matches_stars("fastcats.fasta"))
