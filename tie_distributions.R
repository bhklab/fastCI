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
