import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="noisy cells filtering")


parser.add_argument("input1", metavar='SegCopy', action='store',
        help='cnvs file.', type=str)
parser.add_argument("input2", metavar='per_cell_summary_metrics.csv', action='store',
        help='10x csv file.', type=str)

parser.add_argument("outdir", metavar='outdir', action='store',
        help='Output directory path.', type=str)

args=parser.parse_args()

input1 = args.input1
input2 = args.input2
outdir = args.outdir

df = pd.read_csv(input1, sep="\t")
df2 = pd.read_csv(input2, usecols=['barcode', 'is_noisy'])

# select the barcodes corresponding to noisy cells
barcodes = df2['barcode'][df2['is_noisy'] == 1].str[:-2]

# drop columns corresponding to noisy cells
df = df.drop(columns=barcodes, axis=1)

df.to_csv(outdir+'/SegCopy_filtered', sep='\t', index=False)
