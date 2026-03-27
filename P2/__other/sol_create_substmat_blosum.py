from Bio import SeqIO
import math
from collections import defaultdict

def create_blosum_matrix(fasta_file, identity_threshold=0.62, 
                         aa_alphabet="ACDEFGHIKLMNPQRSTVWY"):
    """
    Create a BLOSUM-style substitution matrix.
    
    The methodology follows Henikoff & Henikoff (1992):
    1. Cluster sequences at identity threshold
    2. Weight each cluster equally (not each sequence)
    3. Count substitutions between clusters
    4. Compute log-odds scores
    
    Args:
        fasta_file (str): Path to FASTA alignment (should be ungapped blocks ideally)
        identity_threshold (float): Clustering threshold (e.g., 0.62 for BLOSUM62)
        aa_alphabet (str): Valid amino acids
        
    Returns:
        dict: Substitution matrix
    """
    records = list(SeqIO.parse(fasta_file, "fasta"))
    
    # Step 1: Cluster sequences
    clusters = cluster_sequences_blosum(records, identity_threshold)
    
    # Step 2: Initialize counters
    pair_counts = defaultdict(float)  # (aa1, aa2) -> weighted count
    total_pairs = 0
    
    # Step 3: Count pairs between and within clusters
    num_clusters = len(clusters)
    
    for i in range(num_clusters):
        for j in range(i, num_clusters):  # Include i==j for within-cluster pairs
            cluster1 = clusters[i]
            cluster2 = clusters[j]
            
            # Weight for this cluster pair
            # Each cluster pair comparison gets weight 1
            # Divided by number of sequence pairs between/within clusters
            n1 = len(cluster1)
            n2 = len(cluster2)
            
            if i == j:  # Within same cluster
                num_pairs = n1 * (n1 - 1) / 2 if n1 > 1 else 0
                if num_pairs == 0:
                    continue
                weight = 1.0 / num_pairs
                
                # Compare all pairs within cluster
                for idx1 in range(n1):
                    for idx2 in range(idx1 + 1, n1):
                        seq1 = records[cluster1[idx1]].seq
                        seq2 = records[cluster1[idx2]].seq
                        count_pairs(seq1, seq2, pair_counts, weight, aa_alphabet)
                        total_pairs += weight
            else:  # Between different clusters
                num_pairs = n1 * n2
                weight = 1.0 / num_pairs
                
                # Compare all pairs between clusters
                for idx1 in cluster1:
                    for idx2 in cluster2:
                        seq1 = records[idx1].seq
                        seq2 = records[idx2].seq
                        count_pairs(seq1, seq2, pair_counts, weight, aa_alphabet)
                        total_pairs += weight
    
    # Step 4: Compute amino acid frequencies from pair counts
    aa_counts = defaultdict(float)
    for (aa1, aa2), count in pair_counts.items():
        aa_counts[aa1] += count
        aa_counts[aa2] += count
    
    total_aa = sum(aa_counts.values())
    aa_freq = {aa: count / total_aa for aa, count in aa_counts.items()}
    
    # Normalize pair counts
    pair_freq = {pair: count / total_pairs for pair, count in pair_counts.items()}
    
    # Step 5: Compute log-odds matrix
    subst_mat = {}
    for aa1 in aa_alphabet:
        subst_mat[aa1] = {}
        for aa2 in aa_alphabet:
            # Get observed frequency
            pair = tuple(sorted([aa1, aa2]))
            observed = pair_freq.get(pair, 0)
            
            # Compute expected frequency
            f1 = aa_freq.get(aa1, 1e-10)
            f2 = aa_freq.get(aa2, 1e-10)
            
            if aa1 == aa2:
                expected = f1 * f2
            else:
                expected = 2 * f1 * f2
            
            # Compute log-odds score (BLOSUM uses log base 2, scaled by 2)
            if observed == 0 or expected == 0:
                log_odds = -10
            else:
                # Using log2 and multiplying by 2 is standard for BLOSUM
                # This gives scores in "half-bits"
                log_odds = math.log2(observed / expected) * 2
            
            subst_mat[aa1][aa2] = int(round(log_odds, 0))
    
    return subst_mat


def count_pairs(seq1, seq2, pair_counts, weight, aa_alphabet):
    """Count amino acid pairs between two sequences."""
    for k in range(len(seq1)):
        aa1 = str(seq1[k])
        aa2 = str(seq2[k])
        
        # Skip gaps
        if aa1 == "-" or aa2 == "-":
            continue
        
        # Check valid amino acids
        if aa1 not in aa_alphabet or aa2 not in aa_alphabet:
            continue
        
        # Store in canonical order for symmetry
        pair = tuple(sorted([aa1, aa2]))
        pair_counts[pair] += weight


def cluster_sequences_blosum(records, identity_threshold):
    """
    Cluster sequences using BLOSUM-style clustering.
    
    Uses a greedy approach where sequences are clustered if they share
    >= identity_threshold with any sequence already in the cluster.
    
    This is a simplified version; true BLOSUM uses more sophisticated clustering.
    """
    n = len(records)
    clusters = []
    assigned = [False] * n
    
    for i in range(n):
        if assigned[i]:
            continue
        
        # Start new cluster
        cluster = [i]
        assigned[i] = True
        
        # Try to add remaining sequences
        for j in range(i + 1, n):
            if assigned[j]:
                continue
            
            # Check if j is similar to ANY sequence in the cluster
            seq_j = records[j].seq
            for idx in cluster:
                seq_idx = records[idx].seq
                identity = calculate_identity(seq_idx, seq_j)
                
                if identity >= identity_threshold:
                    cluster.append(j)
                    assigned[j] = True
                    break  # Only need to match one sequence in cluster
        
        clusters.append(cluster)
    
    return clusters


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