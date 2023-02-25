

library(data.table)
library(dplyr)
library(purrr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(viridis)
library(ggsci)
library(ggcorrplot)
library(pheatmap)
library(ggsignif)
library(plyr)
library(vsn)
library(tidyr)
library(DESeq2)
library(RColorBrewer)
library(PCAtools)
library(RUVSeq)
library(glmGamPoi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)

directory <- getwd()
directory <- getwd
directory <- getwd()
directory <- getwd()


# Prepare sample data
files <- list.files(directory, pattern ="*_fcounts.txt")
df <- map(files, .f=fread,header=TRUE, select = c(6:8))

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

id = c(1:length(df))
for (i in id){
  df[[id[i]]] <- df[[id[i]]] %>% dplyr::mutate(tpm(df[[id[i]]][[paste(colnames(df[[id[i]]][,3]))]], df[[id[i]]]$Length), .before=3)
  df[[id[i]]] <- df[[id[i]]] %>%  dplyr::select(2,3)
  df[[id[i]]] <- df[[id[i]]][!duplicated(df[[id[i]]][["gene_name"]]),]
  df[[id[i]]] <- drop_na(df[[id[i]]])
}

# Get count data
cts <- join_all(df, by = "gene_name", match = "all")
cts = cts %>%  remove_rownames %>% column_to_rownames(var="gene_name")
colnames(cts) = gsub("\\_.*","",files)


# Prepare metadata
sampleTable <- read.table("SraRunTable.txt", header=F, sep = ",", skip = 1) %>%  dplyr::select(V1, V23)
colnames(sampleTable) = c("sampleName", "X")
sampleFiles <- as.data.frame(list.files(directory, pattern = "*_fcounts.txt"))
colnames(sampleFiles) = c("sampleFiles")
sampleFiles$sampleName = gsub("\\_.*","", sampleFiles$sampleFiles)
sampleTable <- merge(sampleTable, sampleFiles, by.y = "sampleName")
coldata <- sampleTable %>% mutate(condition = case_when(startsWith(as.character(X), "N") ~ "HC",
                                                            startsWith(as.character(X), "g") ~ "DIS")) %>% dplyr::select("sampleName", "condition")

# Get counts for DIS 
t.cts <- as.data.frame(t(cts))
dis.cts <- tibble::rownames_to_column(t.cts, "sampleName")
cts.dis <- merge(coldata, dis.cts, by.y="sampleName") %>% filter(condition == "DIS") 
dis.cts <- as.data.frame(t(cts.dis))
colnames(dis.cts) <- dis.cts[1,]
dis.cts <- dis.cts[-c(1,2),]
write.csv(dis.cts, "04_Brain.cts.csv")



# Continue for degs
rownames(coldata) = coldata$sampleName
coldata$sampleName = NULL
coldata$condition <- factor(coldata$condition)
condition <- factor(colnames(cts))

# cross-check data
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))

# Pre-filter genes with low counts
dds <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Variance stabilized
dds1 <- dds
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd)
rv <- rowVars(assay(vsd))

#distance clustering
distance <- dist(t(vsd@assays@data@listData[[1]]), method = "maximum") 
clusters <- hclust(distance, method = "average") 
plot(clusters, labels=t(coldata))

#PCA
#pc <- prcomp(t(assay(vsd)[head(order(-rv),1000),]))
#idx <- pc$x[,2] > 5
#sum(idx)
#plot(pc$x[,1:2], col=idx+1, pch=20, asp=1)
#condition <- condition[!idx]
#dds <- dds[,!idx]
#vsd <- vst(dds, blind=FALSE)

#remove outliers from metadata
#c <- as.data.frame(condition)
#colnames(c) = c("c")
#coldata$c = rownames(coldata)
#coldata <- merge(c, coldata, by.y="c")
#rownames(coldata) = coldata$c
#coldata$c = NULL

#replot PCA and cluster
#plotPCA(vsd)
#distance <- dist(t(vsd@assays@data@listData[[1]]), method = "maximum") 
#clusters <- hclust(distance, method = "average") 
#plot(clusters, labels=t(coldata))

# determine degs
deg <- DESeq(dds)
res <- results(deg, format = c("DataFrame"), alpha = 0.05, 
               pAdjustMethod = "BH", lfcThreshold = 1,
               contrast=c("condition","DIS","HC"))
resOrdered <- res[order(res$padj),]
summary(res)


# Removing unwanted variation 
set <- newSeqExpressionSet(counts(deg))
set <- betweenLaneNormalization(set, which="upper")
not_sig <- rownames(res)[which(res$pvalue > .2)]
empirical <- rownames(set)[ rownames(set) %in% not_sig ]
set <- RUVg(set, empirical, k=1)
pdat <- pData(set)
pdat$condition <- condition
vsd$W1 <- pdat$W_1
plotPCA(vsd, intgroup="W1")
# add factors to design and retest for degs
colData(deg) <- cbind(colData(deg), pdat[,1:2])
design(deg) <- ~W_1 + condition
deg <- DESeq(deg, test="LRT", reduced=~W_1, fitType="glmGamPoi")
res <- results(deg)
resOrdered <- res[order(res$padj),]
summary(res)


# assessment plots
plotCounts(deg, gene=which.min(res$padj), intgroup="condition")
meanSdPlot(assay(vsd))

# Heatmap vsd
select <- order(rowMeans(counts(deg,normalized=TRUE)), decreasing=TRUE)[1:100]
pheatmap(assay(vsd)[select,], show_rownames=FALSE, annotation_col=coldata,
         angle_col = "45")

# Sample correlation
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Get degs
# subset degs
data <- as.data.frame(resOrdered) %>% dplyr::select(log2FoldChange, padj)
data <- tibble::rownames_to_column(data, "Gene")
rownames(data) <- NULL
df <- data %>% drop_na()

# volcano
df.v <- df
df.v$Change <- "No change"
df.v$Change[df.v$log2FoldChange > 1.5 & df.v$padj < 0.05] <- "Increased"
df.v$Change[df.v$log2FoldChange < -1.5 &  df.v$padj  < 0.05] <- "Decreased"
df.v$delabel <- NA
df.v$delabel[df.v$Change != "No change"] <- df.v$Gene[df.v$Change != "No change"]
#plot
#x = paste(name, "_volcano")
ggplot(df.v, aes(x=log2FoldChange, y=-log10(padj), col=Change, label=delabel)) +  geom_point() +  
  theme_classic() + geom_text_repel() + labs(x = "Fold change (log2FoldChange)", y = "Significance (-log10(adj.P.Val)) ") + 
  scale_color_manual(values=c("blue", "red", "black")) 


### Functional enrichment
# Get degs for functional enrichment
df.p <- dplyr::filter(df, padj < 0.05) %>% dplyr::filter(log2FoldChange < -1 | log2FoldChange > 1) %>%  dplyr::select(c(1,2))
df.p <- df.p[order(-df.p$log2FoldChange),]
colnames(df.p) = c("Genes", "logFC")

#GO analysis
original_gene_list <- df.p$logFC
names(original_gene_list) <- df.p$Gene
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

#GO analysis
GO <- gseGO( geneList = gene_list ,  ont = "BP",  OrgDb = org.Hs.eg.db,
             keyType = "SYMBOL",   exponent = 1,  minGSSize = 10,  maxGSSize = 200,
             eps = 1e-10,   pvalueCutoff = 0.5,   pAdjustMethod = "BH",   verbose = TRUE,
             seed = FALSE,  by = "fgsea" )

# visualizations
dotplot(GO, orderBy = "x", split=".sign", showCategory=10) + facet_grid(~.sign)
edox <- setReadable(GO, 'org.Hs.eg.db', 'SYMBOL')
heatplot(edox, foldChange=gene_list, showCategory=30)
p1 <- cnetplot(edox, foldChange=gene_list)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=gene_list)
p3 <- cnetplot(edox, foldChange=gene_list, circular = TRUE, colorEdge = TRUE) 
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
edox2 <- pairwise_termsim(edox)
treeplot(edox2, hclust_method = "average")
edo <- pairwise_termsim(GO)
emapplot(edo, cex_category=1.5,layout="kk") 
upsetplot(edo)

#Save files
#DEGs <- df %>% filter(padj < 0.05) %>% dplyr::select(-c(3))
write.csv(df.p, "04.Brain.deg.csv")
GOs <- as.data.frame(GO@result %>% dplyr::select("Description", "NES", "p.adjust", "qvalue", "core_enrichment") %>% filter(p.adjust < 0.05))
write.csv(GOs, "04.Brain.GO.csv")
