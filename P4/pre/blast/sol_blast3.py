from Bio.Align import substitution_matrices


def read_fasta_sequence(filepath):
    """
    Read the first sequence from a FASTA file and return it as a plain string.
    Lines starting with '>' are headers and are skipped.
    """
    seq_lines = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_lines.append(line.upper())
    return "".join(seq_lines)


def read_fasta_database(filepath):
    """
    Read all sequences from a FASTA file and return a list of plain strings,
    one per record.
    """
    sequences = []
    current = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current).upper())
                    current = []
            else:
                current.append(line)
    if current:
        sequences.append("".join(current).upper())
    return sequences


def get_kmers(seq, k):
    """Return a list of tuples (k-mer, position) from the sequence."""
    return [(seq[i : i + k], i) for i in range(len(seq) - k + 1)]


def score_kmer_pair(word1, word2, subst_mat):
    """Return the BLOSUM62 score between two same-length k-mers."""
    return sum(subst_mat[word1[i]][word2[i]] for i in range(len(word1)))


def find_seeds(query_kmers, database, k, t):
    """
    Return all seeds whose BLOSUM62 score against the query k-mer meets
    threshold t. Each seed records both the query word and the matching
    db word (which may differ).
    """
    subst_mat = substitution_matrices.load("BLOSUM62")
    seeds = []
    for db_index, db_seq in enumerate(database):
        for db_pos in range(len(db_seq) - k + 1):
            db_word = db_seq[db_pos : db_pos + k]
            for q_word, q_pos in query_kmers:
                score = score_kmer_pair(q_word, db_word, subst_mat)
                if score >= t:
                    seeds.append(
                        {
                            "db_index": db_index,
                            "db_pos": db_pos,
                            "q_pos": q_pos,
                            "kmer": q_word,  # original query k-mer
                            "db_kmer": db_word,  # matched db k-mer (may differ)
                            "score": score,
                        }
                    )
    return seeds


def extend_seeds(seeds, query, database, k):
    """For each seed, extend left and right to create HSPs."""
    hits = []
    seen = set()  # FIX 1: deduplicate identical hits before merging

    for s in seeds:
        q_pos = s["q_pos"]
        db_index = s["db_index"]
        db_pos = s["db_pos"]
        db_seq = database[db_index]

        left = 0
        right = k

        while (
            q_pos + right < len(query)
            and db_pos + right < len(db_seq)
            and query[q_pos + right] == db_seq[db_pos + right]
        ):
            right += 1

        while (
            q_pos - left - 1 >= 0
            and db_pos - left - 1 >= 0
            and query[q_pos - left - 1] == db_seq[db_pos - left - 1]
        ):
            left += 1

        q_start = q_pos - left
        q_end = q_pos + right - 1
        db_start = db_pos - left
        db_end = db_pos + right - 1

        # FIX 1: skip exact duplicate extended hits
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
                "db_seq": db_seq[db_start : db_end + 1],
                "seed": s["kmer"],
            }
        )

    return hits


def merge_overlapping_hits(hits, query, database):
    """
    Merge overlapping hits per (database sequence, diagonal).

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
                current["db_seq"] = database[db_index][
                    current["db_start"] : current["db_end"] + 1
                ]
            else:
                merged_hits.append(current)
                current = dict(h)

        merged_hits.append(current)

    return merged_hits


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print("Usage: python blast_fixed.py <query.fasta> <database.fasta>")
        sys.exit(1)

    query_file, db_file = sys.argv[1], sys.argv[2]

    query = read_fasta_sequence(query_file)
    database = read_fasta_database(db_file)

    print(f"Query ({len(query)} aa): {query}")
    print(f"Database: {len(database)} sequence(s)")

    k = 3
    t = 11  # BLOSUM62 score threshold

    query_kmers = get_kmers(query, k)

    seeds = find_seeds(query_kmers, database, k, t)
    print(f"Found seeds (threshold t={t}):")
    for s in seeds:
        print(
            f"  Query[{s['q_pos']}:{s['q_pos']+k}]='{s['kmer']}'"
            f" ~ DB[{s['db_index']}][{s['db_pos']}:{s['db_pos']+k}]='{s['db_kmer']}'"
            f"  score={s['score']}"
        )

    hits = extend_seeds(seeds, query, database, k)
    print("\nExtended hits (deduplicated):")
    for h in hits:
        print(
            f"  DB[{h['db_index']}] Q[{h['query_start']}:{h['query_end']}]"
            f" -> DB[{h['db_start']}:{h['db_end']}]"
        )
        print(f"    Seed: {h['seed']}")
        print(f"    Query: {h['query_seq']}")
        print(f"    DB:    {h['db_seq']}")

    merged = merge_overlapping_hits(hits, query, database)
    print("\nMerged hits (longest HSPs per diagonal):")
    for h in merged:
        diagonal = h["db_start"] - h["query_start"]
        print(
            f"  DB[{h['db_index']}] diag={diagonal}"
            f" Q[{h['query_start']}:{h['query_end']}]"
            f" -> DB[{h['db_start']}:{h['db_end']}]"
        )
        print(f"    Query: {h['query_seq']}")
        print(f"    DB:    {h['db_seq']}")
