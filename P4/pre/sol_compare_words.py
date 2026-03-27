from Bio.Align import substitution_matrices
from fasta2words import fasta2words


def compare_words(seq_fasta, db_fasta, t=8, k=3):
    """Scores one sequnce in the seq_fasta file with all sequences in the db_fasta file.
    The score should be the sum of the individuall scores of all word of length k with
    a threshold t using a substitution matrix.
    >>> compare_words("SGMR2.fasta", "protein_database.fasta", t=8, k=3)[:3]
    [(6568.0, 'A0A2K5EYW3_AOTNA'), (6448.0, 'K7D4N0_PANTR'), (6447.0, 'A0A2I2ZIR7_GORGO')]
    """
    subst_mat = substitution_matrices.load("BLOSUM62")
    protein_as_words = fasta2words(seq_fasta, k)
    database_as_words = fasta2words(db_fasta, k)

    # this dictionary contains only one value (one sequence)
    prot_words = list(protein_as_words.values())[0]

    results = []
    # all proteins in the database (many sequences)
    for seq_id, words in database_as_words.items():
        # all words in the database
        seq_score = 0
        for word in words:
            # against all words in the query
            for prot_word in prot_words:
                # compute the score of each character and sum
                word_score = sum(subst_mat[prot_word[i]][word[i]] for i in range(k))
                if word_score >= t:
                    seq_score += word_score
        results.append((seq_score, seq_id))

    # sort, reverse and display
    return list(reversed(sorted(results)))
