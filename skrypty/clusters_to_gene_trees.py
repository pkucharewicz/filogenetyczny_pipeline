import argparse
import os
import subprocess
from Bio import AlignIO
from Bio.Phylo.Applications import PhymlCommandline
from Bio import pairwise2
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

def mmseqs_to_files(file,n_genomes,outdir):
    clusters=[]
    cluster=set()
    prots=[]
    id=None
    first_cluster=True
    with open(file) as f:
        for l in f.readlines():
            if first_cluster:
                cluster_file=outdir+"/"+l.strip().split(" ")[0][1:]+".faa"
                first_cluster=False
                continue
            if l.startswith(">"):
                if id==l.strip().split(" ")[0]:
                    clusters.append(cluster)
                    if len(cluster)>=n_genomes:
                        clust_file=open(cluster_file,"w")
                        for genome,seq in prots:
                            clust_file.write(f'>{genome}\n{seq}\n')
                    cluster=set()
                    prots=[]
                    cluster_file=outdir+"/"+id[1:]+".faa"
                id=l.strip().split(" ")[0]
                genome=l.strip().split('[')[-1][:-1]
            else:
                cluster.add(genome)
                prots.append((genome.replace(" ","_"),l.strip()))

    #the last cluster
    clusters.append(cluster)
    if len(cluster)>=n_genomes:
        clust_file=open(cluster_file,"w")
        for genome,seq in prots:
            clust_file.write(f'>{genome}\n{seq}\n')
   
    return clusters


def get_MSA(seq_file, msa_file):
    mafft_out = open(msa_file, 'w')
    subprocess.call(["mafft","--auto",seq_file], stdout=mafft_out)
    mafft_out.close()

def cluster_to_bijective(file,out_file):
    repr_seq=False
    prots={}
    with open(file) as f:
        for l in f.readlines():
            l=l.strip()
            if l.startswith(">"):
                genome=l[1:]
                if genome not in prots:
                    prots[genome]=[]
            elif l:
                if not repr_seq:
                    repr_seq=l
                prot_seq=l.strip()
                prots[genome].append(prot_seq)
    new_file=open(out_file,"w")
    for genome, seqs in prots.items():
        if len(seqs)==1:
            new_file.write(f">{genome}\n{seqs[0]}\n")
        else:
            max_score=pairwise2.align.globalxx(repr_seq,seqs[0],score_only=True)
            best_seq=seqs[0]
            for s in seqs[1:]:
                aln_score=pairwise2.align.globalxx(repr_seq,s,score_only=True)
                if aln_score>max_score:
                    max_score=aln_score
                    best_seq=s

            new_file.write(f">{genome}\n{best_seq}\n")
    new_file.close()
   
def process_cluster_non_bijective(cluster_file,out_file):
    prots={}
    lines_n=0
    with open(cluster_file) as f:
        for l in f.readlines():
            lines_n+=1
            l=l.strip()
            if l.startswith(">"):
                genome=l[1:]
                if genome not in prots:
                    prots[genome]=[]
            elif l:
                prot_seq=l.strip()
                prots[genome].append(prot_seq)

    #check if fasta file has more than one sequence
    if lines_n>2:
        new_file=open(out_file,"w")
        for genome in prots:
            for i in range(len(prots[genome])):
                #genomes can repeat, so there is number added at the end of them
                prot_seq = prots[genome][i]
                new_file.write(f">{genome}_id_{i}\n{prot_seq}\n")
        new_file.close()

def construct_NJ(msa_file,out_file):
    aln=AlignIO.read(msa_file,'fasta')
    constructor = DistanceTreeConstructor()
    calculator = DistanceCalculator('identity')
    dist_M = calculator.get_distance(aln)
    NJ_tree = constructor.nj(dist_M)
    Phylo.write(NJ_tree, out_file,'newick')


parser = argparse.ArgumentParser()
parser.add_argument("clusters", help="fasta with all clusters (fasta mmseqs output)")
parser.add_argument("outdir", help="dir to write output to")
parser.add_argument("n", help="minimal number of genomes in cluster") 
parser.add_argument("--bijective", action='store_true', help="if clusters chould be unique (1 genome - 1 protein)")           
args = parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

clusters_dir=args.outdir+"/clusters"
if not os.path.exists(clusters_dir):
   os.makedirs(clusters_dir)

clusters_to_align_dir=args.outdir+"/clusters_to_align"
if not os.path.exists(clusters_to_align_dir):
   os.makedirs(clusters_to_align_dir)

MSA_dir=args.outdir+"/alignments"
if not os.path.exists(MSA_dir):
   os.makedirs(MSA_dir)

trees_dir=args.outdir+"/gene_trees"
if not os.path.exists(trees_dir):
   os.makedirs(trees_dir)

mmseqs_to_files(args.clusters,int(args.n), clusters_dir)

for filename in os.listdir(clusters_dir):
   f = os.path.join(clusters_dir, filename)
   out_file=clusters_to_align_dir+"/"+filename
   if args.bijective:
      cluster_to_bijective(f,out_file)
   else:
      process_cluster_non_bijective(f,out_file)

#MSA for every filtered cluster
n_files=len(os.listdir(clusters_to_align_dir))
MSA_done=0
for filename in os.listdir(clusters_to_align_dir):
    f = os.path.join(clusters_to_align_dir, filename)
    out_file=MSA_dir+"/"+filename
    get_MSA(f,out_file)
    MSA_done+=1
    print(f'{MSA_done}/{n_files} MSA done')

print('All MSA done')

#tree for every MSA

trees_done=0
for filename in os.listdir(MSA_dir):
    f = os.path.join(MSA_dir, filename)
    out_file=trees_dir+"/"+filename+'_tree.nwk'
    try:
        construct_NJ(f,out_file)
        trees_done+=1
        print(f'{trees_done}/{n_files} trees constructed') 
    except:
        print("Tree construction failed for msa file: "+f)
        print(open(f).read())

print("Trees done")