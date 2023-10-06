library(igraph)

args = commandArgs(trailingOnly=TRUE)
outfile <- args[1]

prob_symm <- read.csv('pref_matrix.csv', header=F, sep=",")
verteces_clusters <- read.csv('block_sizes.csv', header=F, sep=",")
vert <- sum(as.integer(verteces_clusters))

g_sbm <- sample_sbm(n=vert, pref.matrix=as.matrix(prob_symm), block.sizes=as.integer(verteces_clusters), directed=F, loops=F)
df <- as.data.frame(as_edgelist(g_sbm))

write.table(df, file=outfile, row.names=FALSE, col.names=FALSE, sep=" ", quote=FALSE)