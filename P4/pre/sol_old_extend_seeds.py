from Bio.Align import substitution_matrices
from build_index import build_index
from find_seeds import find_seeds

subst_mat = substitution_matrices.load("BLOSUM62")


def extend_seeds(seeds, query, db, k, X=15):
    """
    For each seed, extend left and right using score-based extension (BLOSUM62).
    Extension continues as long as the cumulative score does not drop more
    than X below the running maximum score seen so far.

    >>> index = build_db_iseq(["MKTAYIAK"], 3)
    >>> seeds = find_seeds("MKTAYIAK", index, 3, 11)
    >>> hits = extend_seeds(seeds, "MKTAYIAK", ["MKTAYIAK"], 3)
    >>> hits[0]["query_seq"]
    'MKTAYIAK'
    >>> hits[0]["db_seq"]
    'MKTAYIAK'
    """
    hits = []
    seen = set()

    for s in seeds:
        q_pos = s["q_pos"]
        db_iseq = s["db_iseq"]
        db_pos = s["db_pos"]
        db_seq = db[db_iseq]  # look up by integer index into sequence list

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
        if key in seen:
            continue
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
                "seed": s["kmer"],
            }
        )

    return hits
