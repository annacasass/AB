from Bio.Align import substitution_matrices

# from build_index import build_index

subst_mat = substitution_matrices.load("BLOSUM62")


def find_seeds(query, db, k, t):
    """
    Return all seeds whose BLOSUM62 score against a query k-mer meets
    threshold t, using a prebuilt k-mer index for lookups.

    >>> db = {'KTM': [(1, 0)], 'MKT': [(1, 2)], 'MRT': [(0, 0)], 'RTA': [(0, 1)], 'TAY': [(0, 2)], 'TMK': [(1, 1)]}
    >>> find_seeds("MKT", db, 3, 11)
    [{'db_iseq': 1, 'db_pos': 2, 'q_pos': 0, 'q_kmer': 'MKT', 'db_kmer': 'MKT', 'score': 15.0}, {'db_iseq': 0, 'db_pos': 0, 'q_pos': 0, 'q_kmer': 'MKT', 'db_kmer': 'MRT', 'score': 12.0}]
    """
    query_kmers = [(query[i : i + k], i) for i in range(len(query) - k + 1)]
    seeds = []
    for q_word, q_pos in query_kmers:
        for db_word, positions in db.items():
            score = sum(subst_mat[q_word[i]][db_word[i]] for i in range(k))
            if score >= t:
                for db_iseq, db_pos in positions:
                    seeds.append(
                        {
                            "db_iseq": db_iseq,
                            "db_pos": db_pos,
                            "q_pos": q_pos,
                            "q_kmer": q_word,  # original query k-mer
                            "db_kmer": db_word,  # matched db k-mer (may differ)
                            "score": score,
                        }
                    )
    return seeds


if __name__ == "__main__":
    k = 3
    t = 11
    db = {
        "KTM": [(1, 0)],
        "MKT": [(1, 2)],
        "MRT": [(0, 0)],
        "RTA": [(0, 1)],
        "TAY": [(0, 2)],
        "TMK": [(1, 1)],
    }
    print(find_seeds("MKT", db, k, t))
