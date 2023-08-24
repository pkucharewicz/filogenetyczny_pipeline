from Bio import Phylo
from Bio.Phylo import Consensus
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("trees", help="file with newick trees")
parser.add_argument("outdir", help="dir to write files to")           
args = parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

trees = list(Phylo.parse(args.trees, "newick"))

for cutoff in [0.4,0.3,0.2,0.1,0]:
    majority_tree = Consensus.majority_consensus(trees, cutoff)
    #majority_tree.root_with_outgroup({"name":'Gloeobacter_violaceus_PCC_7421'})
    #Phylo.draw(majority_tree)
    Phylo.write(majority_tree, args.outdir+"/consensus_"+str(float(cutoff)),'newick')