from Bio import SeqIO
from Bio.Align import substitution_matrices


def build_db_index(database, k):
    """
    Preprocess the database into a k-mer index.

    Returns a dict:
        { db_kmer: [(db_index, db_pos), ...], ... }

    """
    index = {}
    for db_index, db_kmer in enumerate(database):
        for db_pos in range(len(db_kmer) - k + 1):
            db_word = db_kmer[db_pos : db_pos + k]
            if db_word not in index:
                index[db_word] = []
            index[db_word].append((db_index, db_pos))
    return dict(sorted(index.items()))


def find_seeds(query_kmers, db_index, k, t):
    """
    Return all seeds whose BLOSUM62 score against a query k-mer meets
    threshold t, using a prebuilt k-mer index.

    Instead of scanning every db position, we iterate over the (much
    smaller) set of unique db k-mers and look up their positions directly.
    """
    subst_mat = substitution_matrices.load("BLOSUM62")
    seeds = []
    for q_word, q_pos in query_kmers:
        for db_word, positions in db_index.items():
            # Score the query k-mer against this db k-mer using BLOSUM62
            score = sum(subst_mat[q_word[i]][db_word[i]] for i in range(k))
            if score >= t:
                for db_kmer_index, db_pos in positions:
                    seeds.append(
                        {
                            "db_index": db_kmer_index,
                            "db_pos": db_pos,
                            "q_pos": q_pos,
                            "kmer": q_word,  # original query k-mer
                            "db_kmer": db_word,  # matched db k-mer (may differ)
                            "score": score,
                        }
                    )
    return seeds


def extend_seeds(seeds, query, db, k, X=15):
    """
    For each seed, extend left and right using score-based extension (BLOSUM62).

    Extension continues as long as the cumulative score does not drop more
    than X below the running maximum score seen so far. This mirrors real
    BLAST behaviour: conservative substitutions are allowed mid-extension,
    but the alignment is trimmed back to the highest-scoring endpoint.

    X: dropoff threshold (default 15, same as BLAST ungapped default)
    """
    subst_mat = substitution_matrices.load("BLOSUM62")
    hits = []
    seen = set()

    for s in seeds:
        q_pos = s["q_pos"]
        db_index = s["db_index"]
        db_pos = s["db_pos"]
        db_kmer = db[db_index]

        # --- Extend right ---
        score, best_score, best_right = 0, 0, 0
        right = 0
        while q_pos + k + right < len(query) and db_pos + k + right < len(db_kmer):
            score += subst_mat[query[q_pos + k + right]][db_kmer[db_pos + k + right]]
            if score > best_score:
                best_score = score
                best_right = right + 1
            if score < best_score - X:
                break
            right += 1

        # --- Extend left ---
        score, best_score, best_left = 0, 0, 0
        left = 0
        while q_pos - left - 1 >= 0 and db_pos - left - 1 >= 0:
            score += subst_mat[query[q_pos - left - 1]][db_kmer[db_pos - left - 1]]
            if score > best_score:
                best_score = score
                best_left = left + 1
            if score < best_score - X:
                break
            left += 1

        q_start = q_pos - best_left
        q_end = q_pos + k + best_right - 1
        db_start = db_pos - best_left
        db_end = db_pos + k + best_right - 1

        key = (db_index, q_start, q_end, db_start, db_end)
        if key in seen:
            continue
        seen.add(key)

        hits.append(
            {
                "db_index": db_index,
                "query_start": q_start,
                "query_end": q_end,
                "db_start": db_start,
                "db_end": db_end,
                "query_seq": query[q_start : q_end + 1],
                "db_kmer": db_kmer[db_start : db_end + 1],
                "seed": s["kmer"],
            }
        )

    return hits


def merge_overlapping_hits(hits, query, db):
    """
    Merge overlapping hits per (db sequence, diagonal).

    FIX 2: group by diagonal (db_start - query_start) so that only hits
    on the same alignment path are merged. Hits on different diagonals
    represent different alignment positions and must NOT be merged.
    """
    merged_hits = []

    # Group by (db_index, diagonal)
    groups = {}
    for h in hits:
        diagonal = h["db_start"] - h["query_start"]
        key = (h["db_index"], diagonal)
        groups.setdefault(key, []).append(h)

    for (db_index, _diagonal), group in groups.items():
        group.sort(key=lambda x: x["query_start"])
        current = dict(group[0])  # copy so we don't mutate original

        for h in group[1:]:
            # Overlapping or adjacent on the same diagonal
            if h["query_start"] <= current["query_end"] + 1:
                current["query_end"] = max(current["query_end"], h["query_end"])
                current["db_end"] = max(current["db_end"], h["db_end"])
                current["query_seq"] = query[
                    current["query_start"] : current["query_end"] + 1
                ]
                # FIX 3: db slice is safe because we're on the same diagonal
                current["db_kmer"] = db[db_index][
                    current["db_start"] : current["db_end"] + 1
                ]
            else:
                merged_hits.append(current)
                current = dict(h)

        merged_hits.append(current)

    return merged_hits


if __name__ == "__main__":
    import sys

    k = 3  # k-mer length
    t = 11  # BLOSUM62 score threshold

    query_file = sys.argv[1] if len(sys.argv) > 1 else "query.fasta"
    db_file = sys.argv[2] if len(sys.argv) > 2 else "database.fasta"

    query = str(next(SeqIO.parse(query_file, "fasta")).seq).upper()
    db_file_lst = [str(r.seq).upper() for r in SeqIO.parse(db_file, "fasta")]

    print(f"Query ({len(query_file)} aa): {query[:60]}...")
    print(f"Database: {len(db_file_lst)} sequence(s)")

    # Preprocess database into k-mer indexs
    print(f"Preprocessing database (k={k})...")
    db = build_db_index(db_file_lst, k)
    print(f"  Index contains {len(db)} unique {k}-mers")

    # Find seeds in the queryand extend to hits
    query_kmers = [(query[i : i + k], i) for i in range(len(query) - k + 1)]
    seeds = find_seeds(query_kmers, db, k, t)
    print(f"Found seeds (threshold t={t}):")
    for s in seeds:
        print(
            f"  Query[{s['q_pos']}:{s['q_pos']+k}]='{s['kmer']}'"
            f" ~ DB[{s['db_index']}][{s['db_pos']}:{s['db_pos']+k}]='{s['db_kmer']}'"
            f"  score={s['score']}"
        )

    # Extend seeds to hits
    hits = extend_seeds(seeds, query, db, k)
    print("\nExtended hits (deduplicated):")
    for h in hits:
        print(
            f"  DB[{h['db_index']}] Q[{h['query_start']}:{h['query_end']}]"
            f" -> DB[{h['db_start']}:{h['db_end']}]"
        )
        print(f"    Seed: {h['seed']}")
        print(f"    Query: {h['query_seq']}")
        print(f"    DB:    {h['db_kmer']}")
"""
    # Merge overlapping hits on the same diagonal
    merged = merge_overlapping_hits(hits, query, db)
    print("\nMerged hits (longest HSPs per diagonal):")
    for h in merged:
        diagonal = h["db_start"] - h["query_start"]
        print(
            f"  DB[{h['db_index']}] diag={diagonal}"
            f" Q[{h['query_start']}:{h['query_end']}]"
            f" -> DB[{h['db_start']}:{h['db_end']}]"
        )
        print(f"    Query: {h['query_seq']}")
        print(f"    DB:    {h['db_kmer']}")
"""
