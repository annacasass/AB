def move_gaps(seq1, seq2):
    """
    Generate alignments by moving gaps in one sequence
    in another sequence.
    >>> move_gaps("THEFASTCAT", "THEFATCAT")
    ['THEFASTCAT', '-THEFATCAT', 'T-HEFATCAT', 'TH-EFATCAT', 'THE-FATCAT', 'THEF-ATCAT', 'THEFA-TCAT', 'THEFAT-CAT', 'THEFATC-AT', 'THEFATCA-T', 'THEFATCAT-']
    >>> move_gaps("THEFASTCAT", "AFASTCAT")
    ['THEFASTCAT', '--AFASTCAT', 'A--FASTCAT', 'AF--ASTCAT', 'AFA--STCAT', 'AFAS--TCAT', 'AFAST--CAT', 'AFASTC--AT', 'AFASTCA--T', 'AFASTCAT--']
    >>> move_gaps("THEFATCAT", "THEFASTCAT")
    ['THEFASTCAT', '-THEFATCAT', 'T-HEFATCAT', 'TH-EFATCAT', 'THE-FATCAT', 'THEF-ATCAT', 'THEFA-TCAT', 'THEFAT-CAT', 'THEFATC-AT', 'THEFATCA-T', 'THEFATCAT-']
    >>> move_gaps("AFASTCAT", "THEFASTCAT")
    ['THEFASTCAT', '--AFASTCAT', 'A--FASTCAT', 'AF--ASTCAT', 'AFA--STCAT', 'AFAS--TCAT', 'AFAST--CAT', 'AFASTC--AT', 'AFASTCA--T', 'AFASTCAT--']
    """
    num_gaps = len(seq1) - len(seq2)
    if num_gaps >= 0:
        long_seq = seq1
        short_seq = seq2
    else:
        long_seq = seq2
        short_seq = seq1

    num_gaps = abs(num_gaps)
    aligns = [long_seq]
    aligns.extend(
        (short_seq[:i] + "-" * num_gaps + short_seq[i:] for i in range(len(short_seq)))
    )
    aligns.append(short_seq + "-" * num_gaps)
    return aligns


if __name__ == "__main__":
    for align in move_gaps("THEFASTCAT", "AFASTCAT"):
        print(align)
