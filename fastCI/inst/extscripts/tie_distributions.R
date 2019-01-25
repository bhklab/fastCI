
plotTriangle <- function(x, y, sort.order = c("y", "x")){
  sort.order <- match.arg(sort.order)
  if(sort.order == "x"){
    y.temp <- y
    y <- x
    x <- y.temp
  }
  toRemove <- is.na(x) | is.na(y)
  x <- x[!toRemove]
  y <- y[!toRemove]
  
  y.order <- order(y)
  x <- x[y.order]
  y <- y[y.order]
  incomparable.y <- list()
  for(dup in unique(y[duplicated(y)])){
    idx <- which(y==dup)
    incomparable.y <- c(incomparable.y,list(apply(combinations(length(idx),2), 1, function(x){
      idx[x]
    })))
  }
  incomparable.x <- list()
  for(dup in unique(x[duplicated(x)])){
    idx <- which(x==dup)
    incomparable.x <- c(incomparable.x,list(apply(combinations(length(idx),2), 1, function(x){
      idx[x]
    })))
  }
  
  rect.y <- matrix(0,length(y),length(y))
  for(ii in 1:length(y)){
    rect.y[ii,ii] <- 1
  }
  for(kk in seq_along(incomparable.y)){
    curMatrix <- incomparable.y[[kk]]
    for(colNum in seq_len(NCOL(curMatrix))){
      rect.y[curMatrix[1,colNum],curMatrix[2,colNum]] <- 1
    }
  }
  
  for(kk in seq_along(incomparable.x)){
    curMatrix <- incomparable.x[[kk]]
    for(colNum in seq_len(NCOL(curMatrix))){
      rect.y[curMatrix[1,colNum],curMatrix[2,colNum]] <- rect.y[curMatrix[1,colNum],curMatrix[2,colNum]] + 2
    }
  }
  
  image(rect.y, col=c("white","blue", "red", "purple"), zlim=c(0,3), xaxt='n', yaxt='n', bty="n")
  if(sort.order == "y") {
    legend("bottomright", legend = c("Tied Y", "Tied X", "Tied Both"), fill=c("blue", "red", "purple"))
  } else {
    legend("bottomright", legend = c("Tied X", "Tied Y", "Tied Both"), fill=c("blue", "red", "purple"))
  }
  
}



library(PharmacoGx)




load("~/Code/pSetUploadR/PSets/CCLE.CTRPv2.RData")

drug.sens <- summarizeSensitivityProfiles(CCLE.CTRPv2)

apply(drug.sens, 1, function(x) return(sum(duplicated(na.omit(x)))))

hist(apply(drug.sens, 1, function(x) return(sum(duplicated(na.omit(x))))),breaks=50)


my.features <- rownames(featureInfo(CCLE.CTRPv2, "rnaseq"))[featureInfo(CCLE.CTRPv2, "rnaseq")$GeneBioType=="protein_coding"]

rnaseq <- summarizeMolecularProfiles(CCLE.CTRPv2, "rnaseq", features = my.features)

### Standard variance prefilter: remove genes with variance less than .1 on log scale

myx <- apply(rnaseq, 1, var, na.rm=TRUE) > 0.1

rnaseq <- rnaseq[myx,]

hist(apply(rnaseq, 1, function(x) return(sum(duplicated(na.omit(x))))),breaks=50)

test <- apply(rnaseq, 1, function(x) {
  tt <- table(x)
  tt <- tt[tt>1]
})


names(test) <- NULL

test <- unique(test)


## For arguments sake, taking largest tissue
my.cells <- na.omit(cellNames(CCLE.CTRPv2)[cellInfo(CCLE.CTRPv2)$tissueid == "lung"])

## We will look at the "triangles" for representitive genes from each quadrant of variance, and also for a broad action and targeted therapy

x.var <- apply(rnaseq, 1, var, na.rm=TRUE)

summary(x.var)

my.order <- order(x.var)

rnaseq <- rnaseq[my.order,]

cur.x <- Biobase::exprs(rnaseq[sample(floor(length(my.order)/4),1),])

cur.y <- drug.sens["paclitaxel",]


cur.x <- cur.x[,my.cells]
cur.y <- cur.y[my.cells]

pdf("lung_first_quartile_gene_paclitaxel_tie_dist.pdf", onefile=TRUE)
plotTriangle(cur.x, cur.y)
plotTriangle(cur.x, cur.y, "x")
dev.off()



cur.x <- Biobase::exprs(rnaseq[sample(floor(length(my.order)/4),1)*2,])

cur.y <- drug.sens["paclitaxel",]


cur.x <- cur.x[,my.cells]
cur.y <- cur.y[my.cells]

pdf("lung_second_quartile_gene_paclitaxel_tie_dist.pdf")
plotTriangle(cur.x, cur.y)
plotTriangle(cur.x, cur.y, "x")
dev.off()


cur.x <- Biobase::exprs(rnaseq[sample(floor(length(my.order)/4),1)*3,])

cur.y <- drug.sens["paclitaxel",]


cur.x <- cur.x[,my.cells]
cur.y <- cur.y[my.cells]

pdf("lung_third_quartile_gene_paclitaxel_tie_dist.pdf")
plotTriangle(cur.x, cur.y)
plotTriangle(cur.x, cur.y, "x")
dev.off()



cur.x <- Biobase::exprs(rnaseq[sample(floor(length(my.order)/4),1)*4,])

cur.y <- drug.sens["paclitaxel",]


cur.x <- cur.x[,my.cells]
cur.y <- cur.y[my.cells]

pdf("lung_forth_quartile_gene_paclitaxel_tie_dist.pdf")
plotTriangle(cur.x, cur.y)
plotTriangle(cur.x, cur.y, "x")
dev.off()

## Now doing worst case scenerio

cur.x <- Biobase::exprs(rnaseq[sample(floor(length(my.order)/4),1),])

cur.y <- drug.sens["BRD8958",]


cur.x <- cur.x[,my.cells]
cur.y <- cur.y[my.cells]

pdf("lung_first_quartile_gene_BRD8958_tie_dist.pdf", onefile=TRUE)
plotTriangle(cur.x, cur.y)
plotTriangle(cur.x, cur.y, "x")
dev.off()



cur.x <- Biobase::exprs(rnaseq[sample(floor(length(my.order)/4),1)*2,])

cur.y <- drug.sens["BRD8958",]


cur.x <- cur.x[,my.cells]
cur.y <- cur.y[my.cells]

pdf("lung_second_quartile_gene_BRD8958_tie_dist.pdf")
plotTriangle(cur.x, cur.y)
plotTriangle(cur.x, cur.y, "x")
dev.off()


cur.x <- Biobase::exprs(rnaseq[sample(floor(length(my.order)/4),1)*3,])

cur.y <- drug.sens["BRD8958",]


cur.x <- cur.x[,my.cells]
cur.y <- cur.y[my.cells]

pdf("lung_third_quartile_gene_BRD8958_tie_dist.pdf")
plotTriangle(cur.x, cur.y)
plotTriangle(cur.x, cur.y, "x")
dev.off()



cur.x <- Biobase::exprs(rnaseq[sample(floor(length(my.order)/4),1)*3,])

cur.y <- drug.sens["BRD8958",]


cur.x <- cur.x[,my.cells]
cur.y <- cur.y[my.cells]

pdf("lung_forth_quartile_gene_BRD8958_tie_dist.pdf")
plotTriangle(cur.x, cur.y)
plotTriangle(cur.x, cur.y, "x")
dev.off()



