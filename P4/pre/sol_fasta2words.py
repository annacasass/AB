from Bio import SeqIO


def fasta2words(filename, k):
    """Creates a dictionary of sequence.id as keys and sorted lists of words of length k from a SeqRecord object as values
    >>> fasta2words("thefastcat.fasta", 3)
    {'fast_cat': ['AST', 'CAT', 'EFA', 'FAS', 'HEF', 'STC', 'TCA', 'THE']}
    >>> len(fasta2words("protein_database.fasta", 3))
    10
    >>> prots = fasta2words("protein_database.fasta", 3)
    >>> len(prots['K7D4N0_PANTR'])
    178
    """
    records = SeqIO.parse(filename, "fasta")
    seqs_as_words = {}
    for record in records:
        seq = str(record.seq)
        seq_id = record.id.split("|")[2] if "|" in record.id else record.id
        seqs_as_words[seq_id] = list(
            sorted({seq[i : i + k] for i in range(len(seq) - k + 1)})
        )
    return seqs_as_words
