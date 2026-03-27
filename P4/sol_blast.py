from Bio import SeqIO
from Bio.Align import substitution_matrices
import sys

from build_index import build_index
from find_seeds import find_seeds
from extend_seeds import extend_seeds
from merge_overlapping import merge_overlapping


if __name__ == "__main__":

    k = 3  # k-mer length
    t = 11  # BLOSUM62 score threshold

    query_file = sys.argv[1] if len(sys.argv) > 1 else "query.fasta"
    db_file = sys.argv[2] if len(sys.argv) > 2 else "database.fasta"

    query = str(next(SeqIO.parse(query_file, "fasta")).seq).upper()
    sequences = [str(r.seq).upper() for r in SeqIO.parse(db_file, "fasta")]

    print(f"\nPreprocessing database (k={k})...")
    index = build_index(sequences, k)

    print(f"\nFinding seeds (threshold t={t})...")
    seeds = find_seeds(query, index, k, t)

    print(f"\nExtending seeds...")
    hits = extend_seeds(seeds, query, sequences, k)

    print(f"\nMerging overlapping hits...")
    merged = merge_overlapping(hits, query, sequences)
    merged.sort(key=lambda h: h["query_end"] - h["query_start"], reverse=True)

    print(f"\nMerged hits sorted by length (longest first): {len(merged)}")
    for h in merged[:10]:  # print top 10 longest hits
        diagonal = h["db_start"] - h["query_start"]
        print(
            f"  DB[{h['db_iseq']}] diag={diagonal}"
            f" Q[{h['query_start']}:{h['query_end']}]"
            f" -> DB[{h['db_start']}:{h['db_end']}]"
        )
        print(f"    Query: {h['query_seq']}")
        print(f"    DB:    {h['db_seq']}")
