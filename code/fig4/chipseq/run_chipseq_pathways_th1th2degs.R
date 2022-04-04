#require(ontologyIndex)
require(clusterProfiler)


#load output from linear model and perform pathway analysis
#needs to process df and derive clean gene sets, then match to stat and TF targets

require(msigdbr)
require(readr)
require(dplyr)
require(stringr)
require(tidyr)
require(ggplot2)

source("code/paper_theme.R")

background <- read.csv("genesets_output/background_genes_cleaned.txt", header = F, col.names = "gene")
universe <- as.character(background$gene)

# or use as background all genes between Th1 and Th2 DEGs ??
df_deg <- read.csv("genesets_output/quantDEGS/degs_Th1_Th2_cor0.csv")

# read data
rep <- "rep2"

if(rep == "rep1"){
  df <- read_csv("data/peine_logfc_tidy.csv")
} else if(rep == "rep2"){
  df <- read_csv("data/replicate_2/peine_tidy_rep2_log2fc.csv")
}


df <- df[df$gene %in% df_deg$gene,]

# only keep th1 and th2 sample
df <- df %>% filter(celltype != "Th0") %>% filter(celltype != "Th1/2")

# clean isoforms
df <- df %>% separate(gene, into = c("gene", "iso"), sep = "\\.")

# get maximum fold change per gene
out <- df %>% group_by(gene) %>% summarize(maxFC = max(value))

# only focus on genes with high fold changes
df2 <- left_join(df, out) %>% filter(maxFC > 1) %>% filter(maxFC == value)

# genes where max fc is in th1 or th2 samplerespectively
th1_genes <- df2[df2$celltype=="Th1",] 
th2_genes <- df2[df2$celltype=="Th2",]

th1_genes <- unique(th1_genes$gene)
th2_genes <- unique(th2_genes$gene)


chipseq_genes <- read.csv("genesets_literature/references/chipseq_genes_all.csv")
chipseq_genes <- as.data.frame(chipseq_genes)


# remove tbet pos and neg

chipseq_genes <- chipseq_genes %>% filter(gs_name != "Tbet_pos") %>% filter(gs_name != "Tbet_neg")
chipseq_genes$gs_name[chipseq_genes$gs_name == "Tbet_large"] <- "Tbet"

res_s1 <- enricher(th1_genes, 
                   TERM2GENE = chipseq_genes, 
                   universe = universe, 
                   pAdjustMethod = "fdr", 
                   maxGSSize = 5000, 
                   pvalueCutoff = 1, 
                   qvalueCutoff = 1)
res1 <- res_s1@result
res1$model <- "Th1"

res_s2 <- enricher(th2_genes, 
                   TERM2GENE = chipseq_genes, 
                   universe = universe, 
                   pAdjustMethod = "fdr", 
                   maxGSSize = 5000, 
                   pvalueCutoff = 1, 
                   qvalueCutoff = 1)
res2 <- res_s2@result
res2$model <- "Th2"

res <- bind_rows(res1, res2)

res <- res %>% mutate(pval2 = -log10(pvalue))

g <- ggplot(data = res, aes(x = ID, y = pval2))


g2 <-  g + geom_bar(aes(fill=model), stat = "identity", position = position_dodge()) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "", y = "-log10(pvalue)") + 
  geom_hline(yintercept = -log10(0.05))+
  scale_y_continuous(limits = c(0,NA), expand = expansion(mult = c(0, .1))) +
  scale_fill_manual(values=c(myblue, myred))
print(g2)

sname <- c("figures/fig4/chipseq/chipseq_pathways_th1th2_")
ggsave(paste0(sname, rep, ".png"), width = 2.5, height = 2, bg = "white")
ggsave(paste0(sname, rep, ".svg"), width = 2.5, height = 2, bg = "white")


