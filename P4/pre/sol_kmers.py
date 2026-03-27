def kmers(query, k):
    """Returns a list of all the words of lenght k contained in the query.
    >>> kmers("QLNFQLMSAGQLQ", 3)
    ['QLN', 'LNF', 'NFQ', 'FQL', 'QLM', 'LMS', 'MSA', 'SAG', 'AGQ', 'GQL', 'QLQ']
    >>> kmers("QLNFQLMSAGQLQ", 4)
    ['QLNF', 'LNFQ', 'NFQL', 'FQLM', 'QLMS', 'LMSA', 'MSAG', 'SAGQ', 'AGQL', 'GQLQ']
    """
    return [query[i : i + k] for i in range(len(query) - k + 1)]
