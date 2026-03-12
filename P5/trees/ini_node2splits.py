# Figure out how to to print *leaves*. A *leaf* is a node with no descendants.

from tree_nodes import Node, newick2nodes


def node2splits(i_node, nodes):
    """This function should call itself in a recursive way. The idea is to start
    exploring left then rigth. Takes as input the node id to start
    the recursion and a dictionary of node objects.
    >>> node2splits(0, newick2nodes("tree_abcdef.dnd"))
    ['A'] ['B']
    ['D'] ['E']
    ['C'] ['D', 'E']
    ['C', 'D', 'E'] ['F']
    ['A', 'B'] ['C', 'D', 'E', 'F']
    ['A', 'B', 'C', 'D', 'E', 'F']
    >>> node2splits(0, newick2nodes("hmgb.dnd"))
    ['hmgb_chite'] ['hmgl_wheat']
    ['hmgb_chite', 'hmgl_wheat'] ['hmgl_trybr']
    ['hmgb_chite', 'hmgl_wheat', 'hmgl_trybr'] ['hmgt_mouse']
    ['hmgb_chite', 'hmgl_wheat', 'hmgl_trybr', 'hmgt_mouse']
    """

    # FILL IN THE CODE

    print(left_list, right_list)
    return left_list + right_list


if __name__ == "__main__":
    nodes = newick2nodes("((A,B),((C,(D,E)),F));")
    node2splits(0, nodes)
