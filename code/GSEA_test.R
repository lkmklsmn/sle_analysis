GPL10558<- read.ftable("./GPL10558-50081.txt")
scan('./GPL10558-50081.txt',nlines = 4)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

geoqry <- "GPL10558-50081"
geoqry <- "GPL10558"

if(length(grep(geoqry,list.files()))==0){
  #Download GDS file, put it in the current directory, and load it:
gds <- GEOquery::getGEO(geoqry, destdir=".")
}}
#Or, open an existing GDS file (even if its compressed):
gds <- GEOquery::getGEO(filename=paste0(geoqry,'.soft.gz'))

subGene <- gds@dataTable@table$Symbol
subGene <- grep('LOC|phage|permuted|negative|positive',unique(subGene),invert = T,value = T)

# library("enrichR")
# BiocManager::install("clusterProfiler")
# library("clusterProfiler")
# 
# # https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C7
# genelist_c7 <- read.csv("./data/c7.all.v2024.1.Hs.symbols.gmt",sep = '\t',header = F)
genelist_c7 <- read.gmt("./data/c7.all.v2024.1.Hs.symbols.gmt")
genelist_ha <- read.gmt("./data/h.all.v2024.1.Hs.symbols.gmt")
# 
# library(tidyverse);library("clusterProfiler")
# df <- top_table %>% mutate(diffexpressed= case_when(
#   logFC   > 0 & adj.P.Val < 0.05 ~ 'UP',
#   logFC   < 0 & adj.P.Val < 0.05 ~ 'DOWN',
#   adj.P.Val > 0.05 ~ 'NO'
# ))
# df <-df[df$diffexpressed!="NO",]
# 
# genelist_c7sub <- genelist_c7[genelist_c7$gene%in%df$gene,]
# deg_results_list <- split(df, df$diffexpressed)
# 
# # Run clusterProfiler on each sub-dataframe
# res <- lapply(names(deg_results_list),
#               function(x) enricher(gene = deg_results_list[[x]]$gene_symbol,
#                                    TERM2GENE = bg_genes))
# fgsea
library("fgsea")
# can use GMT pathways

# fgseaRes <- fgsea(pathways = examplePathways, 
#                   stats    = exampleRanks,
#                   minSize  = 15,
#                   maxSize  = 500)
genelist <- list(genelist_c7,genelist_ha)[[2]]
ranks <- top_table$logFC;names(ranks) <- top_table$gene

qry_term<- grep('LUP',unique(genelist$term),value=T); names(qry_term) <- unlist(lapply(qry_term,function(x)paste(strsplit(x[[1]],"_")[[1]][-1],collapse = ' ')))
lapply(names(qry_term),function(qry){
  plotEnrichment(genelist[genelist$term==qry_term[qry],'gene'],
                 ranks) + labs(title=qry)
})
#
egPathways <- split(genelist$gene,genelist$term)

fgseaRes <- fgsea(pathways = egPathways, 
                  stats    = ranks,
                  minSize  = 15,
                  maxSize  = 500)
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(egPathways[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)
