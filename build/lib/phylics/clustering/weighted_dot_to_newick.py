#!/usr/bin/envs python
import re
import sys
import json
import argparse
from Bio import Phylo
from gen_tree import gen_tree 
from io import StringIO
import matplotlib.pyplot as plt
from pprint import pprint

def from_dot_to_dict(dotfile):
    """
        dot file example:
        
        NODES SECTION

        "0: [0,1]";
        "1:[0,0.52],0.0016";
        "5:[0.0000,0.4201],0.0017";
        ...

        EDGES SECTION
        "0: [0,1]" -> "1:[0,0.52],0.0016";
        "0: [0,1]" -> "2:[0.52,1],0.0008";
        "1:[0,0.52],0.0016" -> "5:[0.0000,0.4201],0.0017";
        "1:[0,0.52],0.0016" -> "6:[0.4201,0.5235],0.0021";
        "5:[0.0000,0.4201],0.0017" -> "19:[0.0000,0.2712],0.0000";
        "5:[0.0000,0.4201],0.0017" -> "20:[0.2712,0.4201],0.0007";
        ...
    """
    """
        json example:

        node_to_children = {
            'root': {'b': 3, 'c': 5},
            'a': {'leaf1': 12, 'leaf2': 32},
            'b': {'a': 2, 'leaf3': 21, 'leaf4': 3},
            'c': {'leaf5': 5, 'leaf6': 7}
        }
    """

    tree = {}
    root_node = ""
    with open(dotfile, 'r') as f:
        lines = f.read().splitlines()
        for l in lines:
            if l.startswith("digraph") or l.startswith("}"):
                continue
            nodes = l.split("->")
            if len(nodes) == 1:
                if root_node == "": #nodes section 
                    root_node = re.sub('[\\t\"\\s\;]+', '', nodes[0]).split(":")[0]
            else:               #edges section
                parent_id = re.sub('[\\t\"\\s\;]+', '', nodes[0]).split(":")[0] # strip special chars
                child_id = re.sub('[\\t\"\\s\;]+', '', nodes[1]).split(":")[0]
                edge_weight = re.sub('[\\t\"\\s\;]+', '', nodes[1]).split(",")[2]
                if parent_id not in tree:
                    tree[parent_id] = {}
                tree[parent_id][child_id] = edge_weight 
                """
                    {
                        ...
                        "parent_id": {"child_id": edge_weight}
                    }
                """
    return tree, root_node

def from_json_to_newick(node_to_children, root_node) -> str:
    visited_nodes = set()

    def newick_render_node(name, distance: float) -> str:
        assert name not in visited_nodes, "Error: The tree may not be circular!"

        if name not in node_to_children:
            # Leafs
            return F'{name}:{distance}'
        else:
            # Nodes
            visited_nodes.add(name)
            children = node_to_children[name]
            children_strings = [newick_render_node(child, children[child]) for child in children.keys()]
            children_strings = ",".join(children_strings)
            return F'({children_strings}){name}:{distance}'

    newick_string = newick_render_node(root_node, 0) + ';'

    # Ensure no entries in the dictionary are left unused.
    assert visited_nodes == set(node_to_children.keys()), "Error: some nodes aren't in the tree"

    return newick_string


def main():
    parser = argparse.ArgumentParser(description="Weighted dot to json converter.")

    parser.add_argument("dotfile", metavar="tree.dot", action="store", type=str, help="Dot file.")
    parser.add_argument("outdir", metavar="out/dir/path", action="store", type=str, help="Output directory path.")

    args = parser.parse_args()

    dotfile = args.dotfile
    outdir = args.outdir

    tree, root_node = from_dot_to_dict(dotfile)
    with open(outdir+'/tree.json', 'w') as f:
        json.dump(tree, f)
    #pprint(tree)
    #print("root node: {}".format(root_node))
    newick = from_json_to_newick(tree, root_node=root_node)
    
    # draw newick tree
    handle = StringIO(newick)
    tree = Phylo.read(handle, "newick")
    tree.ladderize()
    Phylo.draw(tree)
    plt.gcf().set_size_inches(21, 40)
    plt.savefig(outdir+"/tree.newick.png")
    plt.clf()

    with open(outdir+"/tree.newick", 'w+') as f:
        f.write(newick)


if __name__ == "__main__":
    sys.exit(main())








