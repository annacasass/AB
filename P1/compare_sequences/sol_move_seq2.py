from score_seqs import score_seqs


def move_seq2(seq1, seq2):
    """
    >>> move_seq2("THEFASTCAT", "THEFATCAT")
    ['THEFASTCAT---------', 'THEFATCAT----------', '-THEFATCAT---------', '--THEFATCAT--------', '---THEFATCAT-------', '----THEFATCAT------', '-----THEFATCAT-----', '------THEFATCAT----', '-------THEFATCAT---', '--------THEFATCAT--', '---------THEFATCAT-', '----------THEFATCAT']
    >>> move_seq2("THEFASTCAT", "AFASTCAT")
    ['THEFASTCAT--------', 'AFASTCAT----------', '-AFASTCAT---------', '--AFASTCAT--------', '---AFASTCAT-------', '----AFASTCAT------', '-----AFASTCAT-----', '------AFASTCAT----', '-------AFASTCAT---', '--------AFASTCAT--', '---------AFASTCAT-', '----------AFASTCAT']
    >>> move_seq2("THEFASTCAT", "THECAT")
    ['THEFASTCAT------', 'THECAT----------', '-THECAT---------', '--THECAT--------', '---THECAT-------', '----THECAT------', '-----THECAT-----', '------THECAT----', '-------THECAT---', '--------THECAT--', '---------THECAT-', '----------THECAT']
    """
    # Total length both sequences need to be
    total_length = len(seq1) + len(seq2)

    # First alignment: seq1 with gaps at the end
    alignments = [seq1 + "-" * len(seq2)]

    # Generate all positions where seq2 can start (0 to len(seq1))
    for start_position in range(len(seq1) + 1):
        leading_gaps = "-" * start_position
        trailing_gaps = "-" * (total_length - start_position - len(seq2))
        aligned_seq2 = leading_gaps + seq2 + trailing_gaps
        alignments.append(aligned_seq2)

    return alignments


def print_scores(seq1, seq2, match, mismatch, gap):
    """
    Score and print the alignments of two sequences with match, mismatch and gap.
    >>> print_scores('THEFASTCAT', 'AFASTCAT', 1, -1, -2)
    THEFASTCAT--------
    AFASTCAT---------- -12
    -AFASTCAT--------- -12
    --AFASTCAT-------- 2
    ---AFASTCAT------- -15
    ----AFASTCAT------ -16
    -----AFASTCAT----- -19
    ------AFASTCAT---- -22
    -------AFASTCAT--- -27
    --------AFASTCAT-- -28
    ---------AFASTCAT- -33
    ----------AFASTCAT -36
    """
    alignments = move_seq2(seq1, seq2)
    seq1_padded = alignments[0]  # First element is padded seq1
    seq2_alignments = alignments[1:]  # Rest are seq2 in different positions

    print(seq1_padded)
    for seq2_aligned in seq2_alignments:
        score = score_seqs(seq1_padded, seq2_aligned, match, mismatch, gap)
        print(seq2_aligned, score)


if __name__ == "__main__":
    print_scores("THEFASTCAT", "THEFATCAT", 1, 0, -1)
