from Bio import SeqIO


def fasta2dict(filename):
    """
    Convert a FASTA file to a dictionary mapping sequence IDs to sequences.

    >>> fasta2dict('cats.fasta')
    {'FAST_CAT': Seq('thefastcat'), 'FAT_CAT': Seq('a--fa-tcat')}
    >>> fasta2dict('rats.fasta')
    {'A_RAT': Seq('a--rat-'), 'THE_RATS': Seq('therats')}
    """
    records = SeqIO.parse(filename, "fasta")

    # FILL IN
    # return dict


if __name__ == "__main__":
    print(fasta2dict("cats.fasta"))
    print(fasta2dict("rats.fasta"))
