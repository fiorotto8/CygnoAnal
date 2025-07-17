import pandas as pd
import ROOT
from array import array
import argparse

parser = argparse.ArgumentParser(description="Convert CSV to ROOT TTree.")
parser.add_argument("--basedir", default="ParamsScan", help="Base directory containing the CSV file")
parser.add_argument("--csv", default="global_output.csv", help="Path to the global output CSV file")
args = parser.parse_args()

BASEDIR = args.basedir
GLOBAL_OUTPUT = BASEDIR + "/" + args.csv
df = pd.read_csv(GLOBAL_OUTPUT, header=0)
df = df.sort_values(by="run")
#print(df)
# Create a ROOT file
root_file = ROOT.TFile("output.root", "RECREATE")

# Create a TTree
tree = ROOT.TTree("data", "Data from CSV")

# Prepare storage for each column
branches = {}
data_types = {
    'int64': 'L',
    'float64': 'D',
    'object': 'C'
}

for col in df.columns:
    dtype = str(df[col].dtype)
    if dtype.startswith('int'):
        arr = array('l', [0])
        tree.Branch(col, arr, f"{col}/L")
    elif dtype.startswith('float'):
        arr = array('d', [0.])
        tree.Branch(col, arr, f"{col}/D")
    else:
        arr = ROOT.std.string()
        tree.Branch(col, arr)
    branches[col] = arr

# Fill the tree
for _, row in df.iterrows():
    for col in df.columns:
        dtype = str(df[col].dtype)
        if dtype.startswith('int'):
            branches[col][0] = int(row[col])
        elif dtype.startswith('float'):
            branches[col][0] = float(row[col])
        else:
            branches[col].replace(0, len(branches[col]), str(row[col]))
    tree.Fill()

# Write and close the file
tree.Write()
root_file.Close()