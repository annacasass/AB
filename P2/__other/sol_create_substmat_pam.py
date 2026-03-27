from Bio import SeqIO
import math
import numpy as np
from collections import defaultdict

def create_pam_matrix(fasta_file, pam_distance=1, min_identity=0.85,
                      aa_alphabet="ACDEFGHIKLMNPQRSTVWY"):
    """
    Create a PAM (Point Accepted Mutation) substitution matrix using
    the original Dayhoff methodology (without phylogenetic trees).
    
    Original PAM methodology (Dayhoff et al., 1978):
    1. Use closely related sequences (≥85% identity)
    2. Count mutations directly (parsimony assumption)
    3. Calculate mutation probability matrix (PAM1)
    4. Extrapolate to higher PAM distances by matrix multiplication
    
    Args:
        fasta_file (str): Path to FASTA alignment of closely related sequences
        pam_distance (int): PAM distance (1 = 1% divergence, 250 = 250%)
        min_identity (float): Minimum sequence identity to include (default 0.85)
        aa_alphabet (str): Valid amino acids
        
    Returns:
        dict: PAM substitution matrix with log-odds scores
    """
    records = list(SeqIO.parse(fasta_file, "fasta"))
    
    print(f"Loaded {len(records)} sequences")
    
    # Step 1: Count mutations from closely related sequences
    mutation_counts, aa_counts = count_mutations_dayhoff(
        records, aa_alphabet, min_identity
    )
    
    # Step 2: Calculate amino acid frequencies
    aa_freq = {aa: count / sum(aa_counts.values()) 
               for aa, count in aa_counts.items()}
    
    # Step 3: Create PAM1 mutation probability matrix
    pam1_matrix = create_pam1_matrix(mutation_counts, aa_counts, aa_alphabet)
    
    # Step 4: Extrapolate to desired PAM distance
    if pam_distance == 1:
        pam_matrix = pam1_matrix
    else:
        pam_matrix = extrapolate_pam(pam1_matrix, pam_distance, aa_alphabet)
    
    # Step 5: Convert to log-odds scores
    log_odds_matrix = convert_to_log_odds(pam_matrix, aa_freq, aa_alphabet)
    
    return log_odds_matrix


def count_mutations_dayhoff(records, aa_alphabet, min_identity):
    """
    Count mutations using Dayhoff's original approach.
    
    - Only use closely related sequence pairs (≥85% identity)
    - Assume parsimony: each difference is one mutation
    - Count both directions (symmetric)
    
    Returns:
        mutation_counts: dict[aa1][aa2] = number of observed aa1→aa2 mutations
        aa_counts: dict[aa] = total occurrences of each amino acid
    """
    mutation_counts = defaultdict(lambda: defaultdict(float))
    aa_counts = defaultdict(float)
    
    # Initialize
    for aa in aa_alphabet:
        aa_counts[aa] = 0
        for aa2 in aa_alphabet:
            mutation_counts[aa][aa2] = 0
    
    n = len(records)
    pairs_used = 0
    
    # Compare all sequence pairs
    for i in range(n):
        for j in range(i + 1, n):
            seq1 = records[i].seq
            seq2 = records[j].seq
            
            if len(seq1) != len(seq2):
                continue
            
            # Calculate identity
            identity = calculate_identity(seq1, seq2)
            
            # Only use closely related sequences (original PAM used ≥85%)
            if identity < min_identity:
                continue
            
            pairs_used += 1
            
            # Count amino acids and mutations
            for k in range(len(seq1)):
                aa1 = str(seq1[k]).upper()
                aa2 = str(seq2[k]).upper()
                
                # Skip gaps
                if aa1 == "-" or aa2 == "-":
                    continue
                
                # Check valid amino acids
                if aa1 not in aa_alphabet or aa2 not in aa_alphabet:
                    continue
                
                # Count amino acid occurrences (for both sequences)
                aa_counts[aa1] += 0.5
                aa_counts[aa2] += 0.5
                
                # Count mutations (differences)
                if aa1 != aa2:
                    # Count mutation in both directions (symmetric)
                    # Each mutation counted once, but we don't know direction
                    # So we split it: 0.5 for aa1→aa2 and 0.5 for aa2→aa1
                    mutation_counts[aa1][aa2] += 0.5
                    mutation_counts[aa2][aa1] += 0.5
    
    print(f"Used {pairs_used} closely related sequence pairs (≥{min_identity:.0%} identity)")
    
    return mutation_counts, aa_counts


def create_pam1_matrix(mutation_counts, aa_counts, aa_alphabet):
    """
    Create the PAM1 mutation probability matrix.
    
    PAM1 represents 1 Point Accepted Mutation per 100 amino acids.
    
    For each amino acid i:
        M[i][j] = probability that i mutates to j in one PAM unit
    
    The matrix is row-stochastic (each row sums to 1).
    """
    M = {}
    
    for aa in aa_alphabet:
        M[aa] = {}
        
        # Calculate total mutations FROM this amino acid
        total_mutations_from_aa = sum(mutation_counts[aa].values())
        total_occurrences = aa_counts[aa]
        
        if total_occurrences == 0:
            # No data for this amino acid
            for aa2 in aa_alphabet:
                M[aa][aa2] = 1.0 / len(aa_alphabet)  # Uniform distribution
            continue
        
        # Calculate mutation probability for each target amino acid
        for aa2 in aa_alphabet:
            if aa == aa2:
                # Diagonal: probability of NO mutation
                # This will be set after we calculate off-diagonal elements
                M[aa][aa2] = 0
            else:
                # Off-diagonal: probability of mutation to aa2
                # Proportional to observed mutations normalized by frequency
                if mutation_counts[aa][aa2] > 0:
                    M[aa][aa2] = mutation_counts[aa][aa2] / total_occurrences
                else:
                    M[aa][aa2] = 0
        
        # Calculate mutation rate (lambda) for this amino acid
        mutation_rate = total_mutations_from_aa / total_occurrences
        
        # Normalize so that the total mutation probability matches PAM1
        # PAM1 = 1% accepted point mutations per 100 residues = 0.01
        # Scale the off-diagonal elements
        off_diagonal_sum = sum(M[aa][aa2] for aa2 in aa_alphabet if aa2 != aa)
        
        if off_diagonal_sum > 0:
            # Scale to achieve 1% mutation rate for PAM1
            scale_factor = 0.01 / off_diagonal_sum
            for aa2 in aa_alphabet:
                if aa2 != aa:
                    M[aa][aa2] *= scale_factor
        
        # Set diagonal (probability of no change)
        M[aa][aa] = 1.0 - sum(M[aa][aa2] for aa2 in aa_alphabet if aa2 != aa)
        
        # Ensure probability is valid
        if M[aa][aa] < 0:
            # Renormalize
            total = sum(M[aa][aa2] for aa2 in aa_alphabet if aa2 != aa)
            if total > 0:
                for aa2 in aa_alphabet:
                    if aa2 != aa:
                        M[aa][aa2] /= (total / 0.99)
            M[aa][aa] = 0.01
    
    return M


def extrapolate_pam(pam1_matrix, n, aa_alphabet):
    """
    Extrapolate from PAM1 to PAMn using matrix multiplication.
    
    PAM_n = PAM1^n
    
    This models n evolutionary time units.
    """
    print(f"Extrapolating from PAM1 to PAM{n}...")
    
    # Convert dict to numpy array
    aa_list = sorted(aa_alphabet)
    size = len(aa_list)
    aa_to_idx = {aa: i for i, aa in enumerate(aa_list)}
    
    # Create numpy matrix
    M = np.zeros((size, size))
    for i, aa1 in enumerate(aa_list):
        for j, aa2 in enumerate(aa_list):
            M[i][j] = pam1_matrix[aa1][aa2]
    
    # Compute matrix power: M^n
    M_n = np.linalg.matrix_power(M, n)
    
    # Convert back to dict
    result = {}
    for i, aa1 in enumerate(aa_list):
        result[aa1] = {}
        for j, aa2 in enumerate(aa_list):
            result[aa1][aa2] = M_n[i][j]
    
    return result


def convert_to_log_odds(prob_matrix, aa_freq, aa_alphabet):
    """
    Convert mutation probability matrix to log-odds scores.
    
    Score(i,j) = 10 * log10(M[i][j] / f[j])
    
    where:
        M[i][j] = probability of amino acid i mutating to j
        f[j] = background frequency of amino acid j
    """
    log_odds = {}
    
    for aa1 in aa_alphabet:
        log_odds[aa1] = {}
        for aa2 in aa_alphabet:
            observed = prob_matrix[aa1][aa2]
            expected = aa_freq[aa2]
            
            if observed > 1e-10 and expected > 1e-10:
                # PAM uses log10 scaled by 10
                score = 10 * math.log10(observed / expected)
            else:
                score = -10  # Large negative score
            
            log_odds[aa1][aa2] = int(round(score, 0))
    
    return log_odds


def calculate_identity(seq1, seq2):
    """Calculate pairwise sequence identity (ignoring gaps)."""
    if len(seq1) != len(seq2):
        return 0.0
    
    matches = 0
    total = 0
    
    for aa1, aa2 in zip(seq1, seq2):
        if aa1 != "-" and aa2 != "-":
            total += 1
            if aa1 == aa2:
                matches += 1
    
    return matches / total if total > 0 else 0.0


def print_matrix_sample(matrix, aa_list="ARNDCQEGHILKMFPSTWYV"):
    """Pretty print a sample of the matrix."""
    sample = aa_list[:8]  # Show first 8 amino acids
    print("\n" + " " * 4 + "  ".join(sample))
    for aa1 in sample:
        scores = [f"{matrix[aa1][aa2]:3d}" for aa2 in sample]
        print(f"{aa1}: " + "  ".join(scores))


# Example usage
if __name__ == "__main__":
    # Create PAM1 (very conservative, for closely related sequences)
    print("Creating PAM1 matrix...")
    pam1 = create_pam_matrix("closely_related.fasta", pam_distance=1)
    print_matrix_sample(pam1)
    
    # Create PAM250 (permissive, for distantly related sequences)
    print("\nCreating PAM250 matrix...")
    pam250 = create_pam_matrix("closely_related.fasta", pam_distance=250)
    print_matrix_sample(pam250)