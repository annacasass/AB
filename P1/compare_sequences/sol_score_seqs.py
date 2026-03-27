def score_seqs(seq1, seq2, match, mismatch, gap):
    """
    Calculate the score of two sequences based on a matching score,
    a mismatching score, and a gap penalty.

    >>> score_seqs("THEFASTCAT", "THEFATCAT-", 1, -1, -2)
    -1
    >>> score_seqs("THEFASTCAT", "THEFATCA-T", 1, -1, -2)
    1
    >>> score_seqs("THEFA-TCAT", "THEFASTCAT", 1, -1, -2)
    7
    >>> score_seqs("THEFASTCAT", "THE", 1, -1, -2)
    >>> score_seqs("THE-FASTCAT", "THE-FASTCAT", 1, -1, -2)
    10
    """
    if len(seq1) == len(seq2):
        score = 0
        for char1, char2 in zip(seq1, seq2):
            if char1 == "-" and char2 == "-":
                pass  # gap-gap alignment scores 0, so skip
            elif char1 == "-" or char2 == "-":
                score += gap  # gap penalty
            elif char1 == char2:
                score += match
            else:
                score += mismatch
        return score


if __name__ == "__main__":
    print(score_seqs("THEFASTCAT", "THEFATCAT", 1, -1, -2))
