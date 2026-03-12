import sys
from tree_nodes import Node, newick2nodes
from fasta2dict import fasta2dict
from align_profiles_names import align_profiles


def node2splits_align(i_node, nodes, seqs, gep):
    """
        Aligns the profiles of the two subtrees of the node i_node, and returns
      the list of names of the sequences in the subtree rooted at i_node.
      
    >>> node2splits_align(0, newick2nodes("hmgb.dnd"), fasta2dict("hmgb.fasta"), -4)
    ['hmgb_chite', 'hmgl_wheat', 'hmgl_trybr', 'hmgt_mouse']
    >>> node2splits_align(1, newick2nodes("hmgb.dnd"), fasta2dict("hmgb.fasta"), -4)
    ['hmgb_chite', 'hmgl_wheat', 'hmgl_trybr']
    >>> node2splits_align(2, newick2nodes("hmgb.dnd"), fasta2dict("hmgb.fasta"), -4)
    ['hmgb_chite', 'hmgl_wheat']
    """

    # a few lines missing

        _, seqs = align_profiles(seqs, left_list, right_list, gep)
    
    # a few lines missing



if __name__ == "__main__":
    nodes = newick2nodes(sys.argv[1])  # "hmgb.dnd"
    id2seq = fasta2dict(sys.argv[2])  # "hmgb.fasta"
    node2splits_align(0, nodes, id2seq, -4)
    for name, seq in id2seq.items():
        print(f">{name}\n{seq}")
