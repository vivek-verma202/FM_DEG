##%######################################################%##
#                                                          #
####                     01. Setup                      ####
#                                                          #
##%######################################################%##

load.lib<-c("AnnotationHub",
            "ensembldb",
            "DESeq2",
            "tidyverse",
            "RColorBrewer",
            "pheatmap",
            "DEGreport",
            "tximport",
            "ggplot2",
            "ggrepel",
            "readxl"
)
install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=T)
install.lib<-load.lib[!load.lib %in% installed.packages()]
if (!requireNamespace("BiocManager", quietly = T))
  install.packages("BiocManager")
BiocManager::install(install.lib)
sapply(load.lib,require,character=T)
theme_set(theme_bw())

##%######################################################%##
#                                                          #
####                02. Annotations                     ####
#                                                          #
##%######################################################%##

# Connect to AnnotationHub
ah <- AnnotationHub()
human_ens <- query(ah, c("Homo sapiens", "EnsDb"))
# latest on 31-Jul-2020: Ensembl 100 (AH79689)
# Extract annotations of interest
human_ens <- human_ens[["AH79689"]]
# Extract gene-level information
genes(human_ens, return.type = "data.frame") %>% View()
# Extract transcript-level information
# transcripts(human_ens, return.type = "data.frame") %>% View()
# Extract exon-level information
# exons(human_ens, return.type = "data.frame") %>% View()
# Create a gene-level dataframe 
annotations_ahb <- genes(human_ens, return.type = "data.frame")  %>%
  dplyr::select(gene_id, entrezid, gene_biotype, symbol)
length(which(map(annotations_ahb$entrezid, length) > 1)) #197
# Wait a second, we don't have one-to-one mappings!
# keep the first identifier multiple mapping cases
annotations_ahb$entrezid <- map(annotations_ahb$entrezid,1) %>%  unlist()
which(is.na(annotations_ahb$symbol)) %>% length()
which(duplicated(annotations_ahb$symbol)) %>% length()
# Return only the non-duplicated genes using index
annotations_ahb <- annotations_ahb[which(duplicated(annotations_ahb$symbol)
                                         == F),]
which(is.na(annotations_ahb$entrezid)) %>%  length()
# Create a transcript dataframe
txdb <- transcripts(human_ens, return.type = "data.frame") %>%
  dplyr::select(tx_id, gene_id)
txdb <- txdb[grep("ENST", txdb$tx_id),]
# Create a gene-level dataframe
genedb <- genes(human_ens, return.type = "data.frame")  %>%
  dplyr::select(gene_id, symbol)
# Merge the two dataframes together
annotations <- inner_join(txdb, genedb)
write.table(annotations,"tx2gene_grch38_ens94.txt",
            quote = F, sep = "\t",row.names = F)

##%######################################################%##
#                                                          #
####                 03. Get quant.sf                   ####
#                                                          #
##%######################################################%##

tx2gene <- read.delim("./tx2gene_grch38_ens100.txt")

samples <- list.files(path = "./salmon_output/",
                      full.names = T, pattern="salmon$")
files <- file.path(samples, "quant.sf")
names(files) <- str_replace(samples, "./salmon_output/", "") %>% 
  str_replace(".Fastq.gz.salmon", "")
all(file.exists(files))
# Run tximport
txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("tx_id", "symbol")], countsFromAbundance="lengthScaledTPM",
                ignoreTxVersion = T)
saveRDS(txi,"txi.RDS")
# txi <- readRDS("txi.RDS")
# Write the counts to an object
data <- txi$counts %>% 
  round() %>% 
  data.frame()
##%######################################################%##
#                                                          #
####                    04. Metadata                    ####
#                                                          #
##%######################################################%##

meta <- read_excel("./FM_pheno.xlsx")
meta <- meta[meta$ID %in% names(files),c(1:4)]
rna_batch <- read.delim("./rna_batch.txt")
meta <- dplyr::inner_join(meta, rna_batch)
meta <- as.data.frame(unclass(meta),stringsAsFactors = T)
rownames(meta) <- meta$ID
sapply(meta, function(x) sum(is.na(x)))
meta["F09",3] <- round(mean(meta$age, na.rm = T))
meta["F09",4] <- "Female"
meta$age <- scale(meta$age)
names(meta)[4] <- "sex"
meta$FM <- relevel(meta$FM, ref = "Control")
## QC
all(colnames(txi$counts) %in% rownames(meta))
all(colnames(txi$counts) == rownames(meta))
ggplot(data) +
  geom_histogram(aes(x = E05), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

mean_counts <- apply(data[,6:8], 1, mean)
variance_counts <- apply(data[,6:8], 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e8)) +
  scale_x_log10(limits = c(1,1e8)) +
  geom_abline(intercept = 0, slope = 1, color="red")

saveRDS(meta,"FM_metadata.RDS")

##%######################################################%##
#                                                          #
####                 05. Normalization                  ####
#                                                          #
##%######################################################%##

## Create DESeq2Dataset object
dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = meta,
                                design = ~ age + sex + batch + FM)
rm(list = c("meta","txi")); gc()
View(counts(dds))
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
write.table(counts(dds, normalized=TRUE),
            file="./normalized_counts.txt",
            sep="\t", quote=F, col.names=NA)

##%######################################################%##
#                                                          #
####                       06. QC                       ####
#                                                          #
##%######################################################%##

### pre-filter unexpressed genes:
genes <- row.names(counts(dds))
exp <- rowSums(counts(dds))
df <- as.data.frame(cbind(genes,exp))
df$exp <- as.numeric(df$exp)
length(which(df$exp==0)) # 12846 out of 35972
keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep,]
rm(list = setdiff(ls(),c("dds","meta"))); gc()
### Transform counts for data visualization
rld <- varianceStabilizingTransformation(dds, blind=T)
plotPCA(rld, intgroup="FM")
plotPCA(rld, intgroup="batch")
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
pheatmap(rld_cor, annotation = meta)
dotchart(rowMeans(rld_cor))
# looks like F62 is a weird sample
dds1 <- dds[,-84]
rld <- varianceStabilizingTransformation(dds, blind=T)
rld1 <- varianceStabilizingTransformation(dds1, blind=T)
plotPCA(rld, intgroup="FM")
plotPCA(rld1, intgroup="FM")
dds <- dds[,-84]
rm(list = setdiff(ls(),"dds")); gc()

##%######################################################%##
#                                                          #
####                     07. DeSeq                      ####
#                                                          #
##%######################################################%##

dds <- DESeq(dds)
sizeFactors(dds)
plotDispEsts(dds)
results(dds, name="FM_Case_vs_Control") %>% data.frame() %>% View()
plotMA(results(dds, name="FM_Case_vs_Control"), ylim=c(-1,1))
# shrunken LFC are better for downstream analyses:
#
resLFC <- lfcShrink(dds, coef="FM_Case_vs_Control", type="apeglm")
plotMA(resLFC, ylim=c(-0.01,0.01))
resLFC %>% data.frame() %>% View()
resLFC <- resLFC %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
write.table(resLFC,
            file="./deseq_res_all_shrunken_LFC.txt",
            sep="\t", quote=F, col.names=NA)
saveRDS(dds,"dds_all.RDS")

##%######################################################%##
#                                                          #
####                  08. Females only                  ####
#                                                          #
##%######################################################%##

rm(list = ls()); gc()
tx2gene <- read.delim("./tx2gene_grch38_ens100.txt")
samples <- list.files(path = "./salmon_output/",
                      full.names = T, pattern="salmon$")
files <- file.path(samples, "quant.sf")
names(files) <- str_replace(samples, "./salmon_output/", "") %>% 
  str_replace(".Fastq.gz.salmon", "")
all(file.exists(files))
meta <- readRDS("./FM_metadata.RDS")
meta <- meta %>% filter(sex == "Female") %>% 
  select(-sex) %>% filter(ID != "F62")
files <- files[names(files) %in% meta$ID]
txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("tx_id", "symbol")], countsFromAbundance="lengthScaledTPM",
                ignoreTxVersion = T)
saveRDS(txi,"txi_F_only.RDS")
dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = meta,
                                        design = ~ age + batch + FM)
dds <- estimateSizeFactors(dds)
write.table(counts(dds, normalized=TRUE),
            file="./normalized_counts_F_only.txt",
            sep="\t", quote=F, col.names=NA)
keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep,] # 12517 out of 35972
rm(list = setdiff(ls(),c("dds","meta"))); gc()
dds <- DESeq(dds)
plotDispEsts(dds)
results(dds, name="FM_Case_vs_Control") %>% data.frame() %>% View()
plotMA(results(dds, name="FM_Case_vs_Control"), ylim=c(-1,1))
resLFC <- lfcShrink(dds, coef="FM_Case_vs_Control", type="apeglm")
plotMA(resLFC, ylim=c(-0.00001,0.00001))
resLFC <- resLFC %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
write.table(resLFC,
            file="./deseq_res_F_only_shrunken_LFC.txt",
            sep="\t", quote=F, col.names=T,row.names = F)
saveRDS(dds,"dds_F_only.RDS")

#decided to proceed with F only
##%######################################################%##
#                                                          #
####                     09. Plots                      ####
#                                                          #
##%######################################################%##

plotCounts(dds, gene="CABIN1", intgroup="FM") 
plotCounts(dds, gene="TMEM8B", intgroup="FM")
plotCounts(dds, gene="JCHAIN", intgroup="FM")
heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap using the metadata data frame for the annotation
norm_count <- read.delim("./normalized_counts_F_only.txt")
norm_count <- norm_count %>% 
  filter(X %in% row.names(rowData(dds))) %>%
  column_to_rownames(var="X") %>% as.matrix()
pheatmap(norm_count, 
         color = c("grey40","white","orange"), 
         cluster_rows = T, 
         show_rownames = F,
         annotation = meta, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

##%######################################################%##
#                                                          #
####              10. Functional analysis               ####
#                                                          #
##%######################################################%##

# BiocManager::install(c("clusterProfiler","DOSE","pathview"))
library(DOSE)
library(pathview)
library(clusterProfiler)
library(tidyverse)
res <- read.delim("./deseq_res_F_only_shrunken_LFC.txt")
annotations_ahb <- read.csv("./annotations_ahb.csv")
res_ids <- inner_join(res, annotations_ahb, by=c("gene"="gene_name"))
sum(is.na(res_ids$entrezid))
hist(res_ids$pvalue,breaks = 10000, ylim = c(0,100))
allOE_genes <- as.character(res_ids$gene)
sigOE <- dplyr::filter(res_ids, pvalue < 0.05)
sigup <- dplyr::filter(sigOE, log2FoldChange > 0)
sigdn <- dplyr::filter(sigOE, log2FoldChange < 0)
sigup_genes <- as.character(sigup$gene)
sigdn_genes <- as.character(sigdn$gene)

egoup <- enrichGO(gene = sigup_genes, 
                universe = allOE_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)









##%######################################################%##
#                                                          #
####                     00. Extras                     ####
#                                                          #
##%######################################################%##

resultsNames(dds)


05. batch effect

vsd <- varianceStabilizingTransformation(dds)
plotPCA(vsd, "batch")
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
plotPCA(vsd, intgroup=c("FM", "batch", "ID"))

06. colData




write.table(names(files),"fn.txt", row.names = F, quote = F)

plot(-log(res$pvalue),res$log2FoldChange, 
     xlim = c(0,5),ylim = c(-0.00001,0.00001))




