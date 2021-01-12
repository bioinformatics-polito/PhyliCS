import sys, os
import argparse
import pandas as pd
from ete3 import Tree, ProfileFace

def remove_prefix(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix):]
    return text  

def plot_tree(tree, data, out_path):
    for lf in tree.iter_leaves():
        lf.add_features(profile = data[int(lf.name)])
        lf.add_features(deviation = [0 for x in range(len(data))])
        lf.add_face(ProfileFace(max_v=12., min_v=0.0, center_v=2.0, width=len(data), height=20, style='heatmap', colorscheme=2), column=1, position="aligned")
    tree.render(os.path.join(out_path, "nwk_heatmap.png"))


def main():
    parser = argparse.ArgumentParser(description="Plot newick tree corresponding heatmap.")

    parser.add_argument("-t", "--tree", required=True, help="Newick tree file.")
    parser.add_argument("-d", "--data", required=True, help="CNV profile matrix.")

    parser.add_argument("-o", "--out_path", required=True, help="Output folder path.")

    args = parser.parse_args()

    t = Tree(args.tree)
    data = pd.read_csv(args.data, sep="\t", usecols = lambda column : column not in ['CHR', 'START', 'END'])

    names = []
    for c in data.columns:
        names.append(int(remove_prefix(c, "cell")))
    data.columns = names
    plot_tree(t, data, args.out_path)

if __name__ == "__main__":
    sys.exit(main())