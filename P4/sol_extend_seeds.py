from Bio.Align import substitution_matrices

subst_mat = substitution_matrices.load("BLOSUM62")


def extend_seeds(seeds, query, db_seqs, k, X=15):
    """Extend seeds in both directions using the BLOSUM62 matrix. Return a list of hits.

    >>> seeds = [{'db_iseq': 1, 'db_pos': 2, 'q_pos': 0, 'q_kmer': 'MKT', 'db_kmer': 'MKT', 'score': 15.0}, {'db_iseq': 0, 'db_pos': 0, 'q_pos': 0, 'q_kmer': 'MKT', 'db_kmer': 'MRT', 'score': 12.0}]
    >>> db_seqs = ["MRTAY", "KTLKT"]
    >>> extend_seeds(seeds, "MKTAY", db_seqs, 3, X=15)
    [{'db_iseq': 1, 'query_start': 0, 'query_end': 2, 'db_start': 2, 'db_end': 4, 'query_seq': 'MKT', 'db_seq': 'LKT', 'seed': 'MKT'}, {'db_iseq': 0, 'query_start': 0, 'query_end': 4, 'db_start': 0, 'db_end': 4, 'query_seq': 'MKTAY', 'db_seq': 'MRTAY', 'seed': 'MKT'}]
    """
    hits = []  # will contain actual hits
    seen = (
        set()
    )  # to avoid duplicates, stores tuples (db_iseq, q_start, q_end, db_start, db_end)

    for s in seeds:
        q_pos = s["q_pos"]
        db_iseq = s["db_iseq"]
        db_pos = s["db_pos"]
        db_seq = db_seqs[db_iseq]

        # --- Extend right ---
        score, best_score, best_right = 0, 0, 0
        right = 0
        while (
            q_pos + k + right < len(query)
            and db_pos + k + right < len(db_seq)
            and score >= best_score - X
        ):
            score += subst_mat[query[q_pos + k + right]][db_seq[db_pos + k + right]]
            if score > best_score:
                best_score = score
                best_right = right + 1
            right += 1

        # --- Extend left ---
        score, best_score, best_left = 0, 0, 0
        left = 0
        while (
            q_pos - left - 1 >= 0 and db_pos - left - 1 >= 0 and score >= best_score - X
        ):
            score += subst_mat[query[q_pos - left - 1]][db_seq[db_pos - left - 1]]
            if score > best_score:
                best_score = score
                best_left = left + 1
            left += 1

        q_start = q_pos - best_left
        q_end = q_pos + k + best_right - 1
        db_start = db_pos - best_left
        db_end = db_pos + k + best_right - 1

        key = (db_iseq, q_start, q_end, db_start, db_end)
        if key not in seen:
            seen.add(key)
            hits.append(
                {
                    "db_iseq": db_iseq,
                    "query_start": q_start,
                    "query_end": q_end,
                    "db_start": db_start,
                    "db_end": db_end,
                    "query_seq": query[q_start : q_end + 1],
                    "db_seq": db_seq[db_start : db_end + 1],
                    "seed": s["q_kmer"],
                }
            )

    return hits


if __name__ == "__main__":
    db_seqs = ["MRTAY", "KTLKT"]
    seeds = [
        {
            "db_iseq": 1,
            "db_pos": 2,
            "q_pos": 0,
            "q_kmer": "MKT",
            "db_kmer": "MKT",
            "score": 15.0,
        },
        {
            "db_iseq": 0,
            "db_pos": 0,
            "q_pos": 0,
            "q_kmer": "MKT",
            "db_kmer": "MRT",
            "score": 12.0,
        },
    ]
    hits = extend_seeds(seeds, "MKTAY", db_seqs, 3, X=15)
    for hit in hits:
        print(hit)
