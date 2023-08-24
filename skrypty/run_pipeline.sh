#bash download_prot_acc.sh
#mmseqs easy-cluster all_proteoms.fasta clusters_mmseqs temp
#python clusters_to_gene_trees.py clusters_mmseqs_all_seqs.fasta bijective 32 --bijective
#python clusters_to_gene_trees.py clusterRes_all_seqs.fasta all 1

#python trees_to_one_file.py bijective/gene_trees final_bijective
#python trees_to_one_file.py all/gene_trees final_all

Rscript constr_consensus.r
python consensus_trees.py final_bijective/all_trees.nwk python_consensus
#./duptree -i final_all/all_trees_unrooted.nwk -o duptree_out_all.nwk
#./duptree -i final_bijective/all_trees_unrooted.nwk -o duptree_out_bijective.nwk

