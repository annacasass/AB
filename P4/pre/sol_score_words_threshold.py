from kmers import *

from Bio.Align import substitution_matrices


def score_words(words, t):
    """
    Write a code that scores each of the words in the query against all the words of the query. For this you will use a BLOSUM62 substitution matrix. The output should be something similar to the cell below: one line per word and a list of tuples with the score and the aligned word sorted from larger to smaller score.
    >>> score_words(['AAA', 'AAS', 'ASA', 'SSA', 'SSS'], 8)
    {'AAA': [(12.0, 'AAA'), (9.0, 'ASA'), (9.0, 'AAS')], 'AAS': [(12.0, 'AAS'), (9.0, 'AAA')], 'ASA': [(12.0, 'ASA'), (9.0, 'SSA'), (9.0, 'AAA')], 'SSA': [(12.0, 'SSA'), (9.0, 'SSS'), (9.0, 'ASA')], 'SSS': [(12.0, 'SSS'), (9.0, 'SSA')]}
    >>> score_words(kmers("QLNFQLMSAGQLQ", 3), 8)
    {'QLN': [(15.0, 'QLN'), (9.0, 'QLQ')], 'LNF': [(16.0, 'LNF')], 'NFQ': [(17.0, 'NFQ')], 'FQL': [(15.0, 'FQL')], 'QLM': [(14.0, 'QLM'), (9.0, 'QLQ')], 'LMS': [(13.0, 'LMS')], 'MSA': [(13.0, 'MSA')], 'SAG': [(14.0, 'SAG')], 'AGQ': [(15.0, 'AGQ')], 'GQL': [(15.0, 'GQL')], 'QLQ': [(14.0, 'QLQ'), (9.0, 'QLN'), (9.0, 'QLM')]}
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
            if score >= t:
                scores_and_words.append((score, db_word))
        sorted_scores = list(reversed(sorted(scores_and_words)))
        scores[prot_word] = sorted_scores
    return scores
