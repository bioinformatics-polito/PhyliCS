
#!/usr/bin/env python

import argparse
import sys
import pandas as pd
import itertools
from Bio import Phylo

def main():
    parser = argparse.ArgumentParser(description="From newick tree to distance matrix.")

    parser.add_argument("newick", metavar="tree.newick", help="Newick tree file", action="store", type=str)
    parser.add_argument("output", metavar="output_prefix", help="Output prefix", action="store", type=str)

    args = parser.parse_args()

    newick = args.newick
    output_prefix = args.output

    
    t = Phylo.read(newick, 'newick')

    d = {}   
    for x, y in itertools.combinations(t.get_terminals(), 2):
        v = t.distance(x, y)
        d[x.name] = d.get(x.name, {})
        d[x.name][y.name] = v
        d[y.name] = d.get(y.name, {})
        d[y.name][x.name] = v
    for x in t.get_terminals():
        d[x.name][x.name] = 0

    m = pd.DataFrame(d)
    m.to_csv(output_prefix+"/tree.mat", sep='\t')

if __name__ == "__main__":
    sys.exit(main())
    
    