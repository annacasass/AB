from Bio import Phylo
from io import StringIO
import os


class Node:
    def __init__(self, name="", left=0, right=0, parent=0, distance=0.0):
        self.name = name
        self.left = left
        self.right = right
        self.parent = parent
        self.distance = distance

    def __repr__(self):
        return f"Node(name={self.name!r}, left={self.left}, right={self.right}, parent={self.parent}), distance={self.distance})"


def newick2nodes(newick):
    # newick can be a file path or a newick string
    if os.path.isfile(newick):
        # read from file
        with open(newick) as f:
            tree = Phylo.read(f, "newick")
    else:
        # read from string
        tree = Phylo.read(StringIO(newick), "newick")
    return _tree2nodes(tree)


def _tree2nodes(tree):
    nodes = {}
    counter = [0]

    def new_id():
        nid = counter[0]
        counter[0] += 1
        return nid

    def build(clade, parent_id):
        nid = new_id()
        nodes[nid] = Node(parent=parent_id, distance=clade.branch_length or 0.0)
        if clade.is_terminal():
            nodes[nid].name = clade.name or ""
        else:
            children = clade.clades
            left_id = build(children[0], nid)
            right_id = build(children[1], nid)
            nodes[nid].left = left_id
            nodes[nid].right = right_id

        return nid

    root_id = build(tree.root, -1)
    nodes[root_id].parent = -1
    return nodes
