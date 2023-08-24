import re
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("dir", help="dir with trees")
parser.add_argument("outdir", help="dir to write files to")           
args = parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

name0="Synechococcus_sp._JA-2-3B'a2-13"
name1="'Synechococcus_sp._JA-2-3B'a(2-13)'"
name_alt="Synechococcus_sp._JA-2-3Ba2-13"

out_file=args.outdir+"/"+'all_trees.nwk'
out_unrooted=args.outdir+"/"+'all_trees_unrooted.nwk'
all_trees=open(out_file,"w")
unrooted_trees=open(out_unrooted,"w")
for filename in os.listdir(args.dir):
   f = os.path.join(args.dir, filename)
   tree0=open(f,"r").read()
   if re.search("\)-|:-",tree0):
      continue
   else:
      tree1=re.sub(name0,name_alt,tree0)
      tree2=re.sub('_id_\d+','',tree1)
      tree3=re.sub(r"'Synechococcus_sp._JA-2-3B\\'a\(2-13\)'",name_alt,tree2)
      tree3=re.sub('\.','',tree3)
      tree3=re.sub('\=|\-','_',tree3)
      # tree3=re.sub('Inner\d+','',tree2)
      all_trees.write(tree3)
      unrooted_trees.write("[&U]"+tree3)

all_trees.close()
unrooted_trees.close()


