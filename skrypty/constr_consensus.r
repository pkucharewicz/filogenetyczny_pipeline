library(ape)
library(TreeTools)

trees_file <- '/home/pampuch/studia/genomika_porównawcza/GP_projekt/final_pipeline/final_bijective/all_trees.nwk'
trees <- read.tree(file=trees_file)
consensus_tree <- consensus(trees, p = 0.5, check.labels = TRUE, rooted = FALSE)
plot(consensus_tree)
write.tree(consensus_tree,"/home/pampuch/studia/genomika_porównawcza/GP_projekt/final_pipeline/consensus_tree_05.nwk")
