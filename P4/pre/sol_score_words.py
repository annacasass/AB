from kmers import *

from Bio.Align import substitution_matrices


def score_words(words):
    """
    Write a code that scores each of the words in the query against all the words of the query. For this you will use a BLOSUM62 substitution matrix. The output should be something similar to the cell below: one line per word and a list of tuples with the score and the aligned word sorted from larger to smaller score.
    >>> score_words(['AAA', 'AAS', 'ASA', 'SSA', 'SSS'])
    {'AAA': [(12.0, 'AAA'), (9.0, 'ASA'), (9.0, 'AAS'), (6.0, 'SSA'), (3.0, 'SSS')], 'AAS': [(12.0, 'AAS'), (9.0, 'AAA'), (6.0, 'SSS'), (6.0, 'ASA'), (3.0, 'SSA')], 'ASA': [(12.0, 'ASA'), (9.0, 'SSA'), (9.0, 'AAA'), (6.0, 'SSS'), (6.0, 'AAS')], 'SSA': [(12.0, 'SSA'), (9.0, 'SSS'), (9.0, 'ASA'), (6.0, 'AAA'), (3.0, 'AAS')], 'SSS': [(12.0, 'SSS'), (9.0, 'SSA'), (6.0, 'ASA'), (6.0, 'AAS'), (3.0, 'AAA')]}
    >>> score_words(kmers("QLNFQLMSAGQLQ", 10))
    {'QLNFQLMSAG': [(49.0, 'QLNFQLMSAG'), (-8.0, 'NFQLMSAGQL'), (-13.0, 'LNFQLMSAGQ'), (-15.0, 'FQLMSAGQLQ')], 'LNFQLMSAGQ': [(49.0, 'LNFQLMSAGQ'), (-3.0, 'FQLMSAGQLQ'), (-13.0, 'QLNFQLMSAG'), (-13.0, 'NFQLMSAGQL')], 'NFQLMSAGQL': [(49.0, 'NFQLMSAGQL'), (-8.0, 'QLNFQLMSAG'), (-12.0, 'FQLMSAGQLQ'), (-13.0, 'LNFQLMSAGQ')], 'FQLMSAGQLQ': [(48.0, 'FQLMSAGQLQ'), (-3.0, 'LNFQLMSAGQ'), (-12.0, 'NFQLMSAGQL'), (-15.0, 'QLNFQLMSAG')]}
    """
    subst_mat = substitution_matrices.load("BLOSUM62")
    scores = {}
    # words in the query
    for prot_word in words:
        scores_and_words = []

        # against words in the database
        for db_word in words:
            # compute the score of each character and sum
            score = sum(
                subst_mat[prot_word[i]][db_word[i]] for i in range(len(words[0]))
            )
            scores_and_words.append((score, db_word))
        sorted_scores = list(reversed(sorted(scores_and_words)))
        scores[prot_word] = sorted_scores
    return scores
