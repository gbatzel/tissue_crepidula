library("edgeR")

# Load in count file and assign groups for each tissue sample
x <- read.delim("/Users/gbatzel/diff_expr/RSEM/gene_counts/matrix/counts_matrix.txt")
group <- c(1,1,1,2,2,2,3,3,3,4,4,4)
y <- DGEList(counts=x[,2:13], genes = x[,1] ,group=group)

# Matrix details
y$samples
dim(y)
str(y)

# Counts per million read
head(cpm(y))

# total gene counts per sample, 2=column of matrix y
apply(y$counts, 2, sum) 

# Filter out low reads
keep <- rowSums(cpm(y)>100) >= 2
y <- y[keep,]
dim(y)

# Reset the library size
y$samples$lib.size <- colSums(y$counts)
y$samples

## Normalizing data
y <- calcNormFactors(y)
y

# Data exploration
plotMDS(y, method="bcv", col=as.numeric(y$samples$group))
#legend("bottomleft", as.character(unique(y$samples$group)), col=1:3, pch=20) #optional legend

# Estimating the dispersion
y1 <- estimateCommonDisp(y, verbose=T)
## Disp = 0.12239 , BCV = 0.3498 
names(y1)

y1 <- estimateTagwiseDisp(y1)
names(y1)

plotBCV(y1)

## GLM estimates of dispersion

design.mat <- model.matrix(~ 0 + y$samples$group)
colnames(design.mat) <- levels(y$samples$group)
y2 <- estimateGLMCommonDisp(y,design.mat)
y2 <- estimateGLMTrendedDisp(y2,design.mat, method="power")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
y2 <- estimateGLMTagwiseDisp(y2,design.mat)
plotBCV(y2)

# differential expression
et21 <- exactTest(y1, pair=c(2,1)) # compare groups 2 and 1: gill vs foot
et31 <- exactTest(y1, pair=c(3,1)) # compare groups 3 and 1: head vs foot
et32 <- exactTest(y1, pair=c(3,2)) # compare groups 3 and 2: head vs gill
et41 <- exactTest(y1, pair=c(4,1)) # compare groups 4 and 1: mantle vs foot
et42 <- exactTest(y1, pair=c(4,2)) # compare groups 4 and 2: mantle vs gill
et43 <- exactTest(y1, pair=c(4,3)) # compare groups 4 and 3: mantle vs head


# Top hits for each of the differential expression comparisons
topTags(et41, n=10)
topTags(et42, n=10)
topTags(et43, n=10)
topTags(et31, n=10)
topTags(et32, n=10)
topTags(et21, n=10)

# Number of upregulated, neutral, and downregulated genes
de1 <- decideTestsDGE(et41, adjust.method="BH", p.value=0.001)
summary(de1)

de2 <- decideTestsDGE(et42, adjust.method="BH", p.value=0.001)
summary(de2)

de3 <- decideTestsDGE(et43, adjust.method="BH", p.value=0.001)
summary(de3)

de4 <- decideTestsDGE(et31, adjust.method="BH", p.value=0.001)
summary(de4)

de5 <- decideTestsDGE(et32, adjust.method="BH", p.value=0.001)
summary(de5)

de6 <- decideTestsDGE(et21, adjust.method="BH", p.value=0.001)
summary(de6)

# differentially expressed tags from the naive method in y1
de1tags41 <- rownames(y1)[as.logical(de1)] 
plotSmear(et41, de.tags=de1tags41)
abline(h = c(-1, 1), col = "blue")

de2tags42 <- rownames(y1)[as.logical(de2)] 
plotSmear(et42, de.tags=de2tags42)
abline(h = c(-1, 1), col = "blue")

de3tags43 <- rownames(y1)[as.logical(de3)] 
plotSmear(et43, de.tags=de3tags43)
abline(h = c(-1, 1), col = "blue")

de4tags31 <- rownames(y1)[as.logical(de4)] 
plotSmear(et31, de.tags=de4tags31)
abline(h = c(-1, 1), col = "blue")

de5tags32 <- rownames(y1)[as.logical(de5)] 
plotSmear(et32, de.tags=de5tags32)
abline(h = c(-1, 1), col = "blue")

de6tags21 <- rownames(y1)[as.logical(de6)] 
plotSmear(et21, de.tags=de6tags21)
abline(h = c(-1, 1), col = "blue")

#Export gene names of diff. expressed list
out41 <- topTags(et41, n=Inf, adjust.method="BH")
keep41 <- out41$table$FDR <= 0.001
gene_list_41 <- out41[keep41,]

out42 <- topTags(et42, n=Inf, adjust.method="BH")
keep42 <- out42$table$FDR <= 0.001
gene_list_42 <- out42[keep42,]

out43 <- topTags(et43, n=Inf, adjust.method="BH")
keep43 <- out43$table$FDR <= 0.001
gene_list_43 <- out43[keep43,]

out21 <- topTags(et21, n=Inf, adjust.method="BH")
keep21 <- out21$table$FDR <= 0.001
gene_list_21 <- out21[keep21,]

out31 <- topTags(et31, n=Inf, adjust.method="BH")
keep31 <- out31$table$FDR <= 0.001
gene_list_31 <- out31[keep31,]

out32 <- topTags(et32, n=Inf, adjust.method="BH")
keep32 <- out32$table$FDR <= 0.001
gene_list_32 <- out32[keep32,]

#Save to file
#write.table(gene_list_32, file="gene_list_32_test.txt")
