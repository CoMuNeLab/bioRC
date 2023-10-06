rm(list=ls())

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=6) {
  stop("Usage: script.R <N> <edgeslist> <signal> <SNR> <MC> <outfile.rds>", call.=FALSE)
} 

source("./RCLib.R")
library(igraph)

N <- as.numeric(args[1])
edges <- read.table(args[2], header=F)
data <- as.matrix(read.table(args[3], header=F))
SNR <- as.numeric(args[4])
MC <- as.numeric(args[5])
outfile <- args[6]

N <- max(c(edges[,1],edges[,2])+1)

print(N)
  
edges[,1] <- edges[,1]+1
edges[,2] <- edges[,2]+1
g <- igraph::graph_from_data_frame(edges, directed=F, vertices=1:N)

data <- addWhiteNoise(data, SNRdb=SNR)
data <- normSignal(data)

######################
# Reservoir setup
######################

#PARAMETERS

inSize = outSize = 1
trainLen = 2500
testLen = 5000
initLen = 100

#COMPLEX NETWORK
cat(paste(" <k>:", mean(degree(g)), "\n"))

W.net <- as.matrix(get.adjacency(g))

#add self-loops
diag(W.net) <- 1

res <- list()
for(m in 1:MC){
  W.wei <- matrix(runif(N*N,-1,1),N)

  W <- W.net * W.wei

  RC <- runReservoir(Signal=data, trainLen=trainLen, testLen=testLen, initLen=initLen, 
                     W, inSize=inSize, outSize=outSize, 
                     RC_leakRate = 0.3, RC_mode = "PREDICTIVE", RC_reg = 1e-8, RC_bias_intensity = 1)

  MSE <- getMSE(RC, trainLen, errorLen=500, norm=T)
  print( paste( 'MSE = ', MSE ) )

  res[[m]] <- RC
}

saveRDS(res, outfile)
