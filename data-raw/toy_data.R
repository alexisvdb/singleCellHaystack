
set.seed(1)


# load data
raw.dat <- read.csv(file = "data-raw/icon.csv", header=F)
raw.dat[is.na(raw.dat)] <- 0
heatmap(as.matrix(raw.dat), Rowv=NA, Colv=NA, scale="none")

raw.coords <- which(raw.dat==1,arr.ind = T)

tsne.x <- rep(NA,nrow(raw.coords))
tsne.y <- rep(NA,nrow(raw.coords))
for(r in 1:nrow(raw.coords)){
  x <- raw.coords[r,1]
  y <- raw.coords[r,2]
  tsne.y[r] <- x + rnorm(1)
  tsne.x[r] <- y + rnorm(1)
}

tsne.x <- tsne.x - mean(tsne.x)
tsne.y <- -(tsne.y - mean(tsne.y))

plot(tsne.x,tsne.y)

dat.tsne <- data.frame(tSNE1=tsne.x,tSNE2=tsne.y, row.names = paste0("cell_",1:601))


###########################
### generate the artificial gene expression data

set.seed(1)
# make artificial expression data
dat.expression <- matrix(exp(rnorm(n=601*500, mean=-3)),ncol=601)*1.5
dat.expression <- dat.expression*runif(min=0,max=10,n=601)
colnames(dat.expression) <- paste0("cell_",1:601)
rownames(dat.expression) <- paste0("gene_",1:500)
hist(apply(dat.expression>1,1,sum))

# add 50 "interesting" genes
selection.rows <- sample(1:nrow(dat.expression), 50)
for(g in 1:50){
  index <- selection.rows[g]

  center.x <- runif(min(tsne.x),max(tsne.x), n=1)
  center.y <- runif(min(tsne.y),max(tsne.y), n=1)

  dists <- (tsne.x-center.x)^2 + (tsne.y-center.y)^2
  o <- order(dists)
  dat.expression[index,] <- runif(min=min(dat.expression),max=.99,n=601)
  expressing.cell.count <- runif(min=10,max=100,n=1)
  dat.expression[index,o[1:expressing.cell.count]] <- runif(min=1.01,max=max(dat.expression),expressing.cell.count)
}

dat.expression <- round(dat.expression-.5,digits = 0)
dat.detection <- dat.expression >= 1
plot(dat.tsne$tSNE1,dat.tsne$tSNE2, col=dat.detection[selection.rows[1],]+1)

save(dat.tsne, dat.expression, file="data/toy.rda")











