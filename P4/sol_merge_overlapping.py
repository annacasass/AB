def merge_overlapping(hits, query, db_seqs):
    """
    Merge overlapping hits per (database sequence, diagonal).

    Group by diagonal (db_start - query_start) so that only hits on the
    same alignment path are merged. Hits on different diagonals represent
    different alignment positions and must NOT be merged.

    >>> h1 = {"db_iseq": 0, "query_start": 0, "query_end": 4, "db_start": 0, "db_end": 4, "query_seq": "MKTAY", "db_seq": "MKTAY", "seed": "MKT"}
    >>> h2 = {"db_iseq": 0, "query_start": 2, "query_end": 6, "db_start": 2, "db_end": 6, "query_seq": "TAYIA", "db_seq": "TAYIA", "seed": "TAY"}
    >>> merged = merge_overlapping([h1, h2], "MKTAYIAK", ["MKTAYIAK"])
    >>> merged[0]['query_seq']
    'MKTAYIA'
    >>> h1 = {"db_iseq": 1, "query_start": 0, "query_end": 2, "db_start": 2, "db_end": 4, "query_seq": "MKT", "db_seq": "LKT", "seed": "MKT"}
    >>> h2 = {"db_iseq": 0, "query_start": 0, "query_end": 4, "db_start": 0, "db_end": 4, "query_seq": "MKTAY", "db_seq": "MRTAY", "seed": "MKT"}
    >>> merged = merge_overlapping([h1, h2], "MKTAY", ["MKTAY", "KTMKT"])
    >>> len(merged)
    2
    """
    merged_hits = []

    groups = {}
    for h in hits:
        diagonal = h["db_start"] - h["query_start"]
        key = (h["db_iseq"], diagonal)
        if key not in groups:
            groups[key] = []
        groups[key].append(h)

    for (db_iseq, _diagonal), group in groups.items():
        group.sort(key=lambda x: x["query_start"])
        current = dict(group[0])

        for h in group[1:]:
            if h["query_start"] <= current["query_end"] + 1:
                current["query_end"] = max(current["query_end"], h["query_end"])
                current["db_end"] = max(current["db_end"], h["db_end"])
                current["query_seq"] = query[
                    current["query_start"] : current["query_end"] + 1
                ]
                current["db_seq"] = db_seqs[db_iseq][
                    current["db_start"] : current["db_end"] + 1
                ]
            else:
                merged_hits.append(current)
                current = dict(h)

        merged_hits.append(current)

    return merged_hits


if __name__ == "__main__":
    h1 = {
        "db_iseq": 0,
        "query_start": 0,
        "query_end": 4,
        "db_start": 0,
        "db_end": 4,
        "query_seq": "MKTAY",
        "db_seq": "MKTAY",
        "seed": "MKT",
    }
    h2 = {
        "db_iseq": 0,
        "query_start": 2,
        "query_end": 6,
        "db_start": 2,
        "db_end": 6,
        "query_seq": "TAYIA",
        "db_seq": "TAYIA",
        "seed": "TAY",
    }
    print(merge_overlapping([h1, h2], "MKTAYIAK", ["MKTAYIAK"]))

    h1 = {
        "db_iseq": 1,
        "query_start": 0,
        "query_end": 2,
        "db_start": 2,
        "db_end": 4,
        "query_seq": "MKT",
        "db_seq": "LKT",
        "seed": "MKT",
    }
    h2 = {
        "db_iseq": 0,
        "query_start": 0,
        "query_end": 4,
        "db_start": 0,
        "db_end": 4,
        "query_seq": "MKTAY",
        "db_seq": "MRTAY",
        "seed": "MKT",
    }
    print(merge_overlapping([h1, h2], "MKTAYIAK", ["MKTAYIAK"]))
