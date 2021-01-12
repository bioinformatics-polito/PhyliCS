#!/usr/bin/envs python
import sys
import numpy
import argparse
from gen_tree import gen_tree
import numpy as np
from anytree.exporter import DotExporter


def main():
    parser = argparse.ArgumentParser(description="From tree to dot.")

    parser.add_argument("npy", metavar="tree.npy", action="store", type=str, help="Tree file.")
    parser.add_argument("outdir", metavar="out/dir/path", action="store", type=str, help="Output directory path.")

    args = parser.parse_args()

    npy = args.npy
    outdir = args.outdir


    tree = numpy.load(npy, allow_pickle=True)
    DotExporter(tree[0]).to_dotfile(outdir+"/tree.dot")


if __name__ == "__main__":
    sys.exit(main())