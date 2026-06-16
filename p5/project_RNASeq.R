# Load the packges we need
# (These should already be installed on the RStudio server for R version 4.0.0)
library(edgeR)
library(baySeq)
library(data.table)
library(clusterProfiler)
library(tidyverse)
library(org.Hs.eg.db)
library(enrichR)
library(ggplot2)

# Read in the data
# You need to set the correct path!
# Skip the first 4 lines of each (they contain summary information, not expression data)
hbr1 <- as.data.frame(fread("~/p5/HBR_Rep1ReadsPerGene.out.tab", skip=4))
hbr2 <- as.data.frame(fread("~/p5/HBR_Rep2ReadsPerGene.out.tab", skip=4))
uhr1 <- as.data.frame(fread("~/p5/UHR_Rep1ReadsPerGene.out.tab", skip=4))
uhr2 <- as.data.frame(fread("~/p5/UHR_Rep2ReadsPerGene.out.tab", skip=4))

# Peek at the data
# There are 6 HBR samples (stored in hbr1 and hbr2)
# There are 6 UHR samples (stored in uhr1 and uhr2)
# Each has the same number of rows (60,710 transcripts)
# The first column has the gene/transcript name
# Columns 2-4 of each data frame contain the expression levels
dim(hbr1)
dim(hbr2)
dim(uhr1)
dim(uhr2)
head(hbr1)
head(hbr2)
head(uhr1)
head(uhr2)

# We need to put the data all into one matrix before applying the functions we need to apply
#   check that the row names are exactly the same first (they are)
# Column-bind the data, leaving out the first column of each (gene names)
# Call the combined data frame expr
# PROVIDE YOUR CODE TO DO THIS!!!
# You should have a gene/transcript per row (60710 rows) and a sample per column (12 samples)
all.equal(hbr1$V1, hbr2$V1) # Returns TRUE
all.equal(hbr1$V1, uhr1$V1) # Returns TRUE
all.equal(hbr1$V1, uhr2$V1) # Returns TRUE
expr <- cbind(hbr1[2:4], hbr2[2:4], uhr1[2:4], uhr2[2:4])
  
  # Convert from data frame to matrix
  # Replace the row names with the transcript names (which are in the right order)
  # Get rid of the column names (they're arbitrary and won't be used)
  expr <- as.matrix(expr)
rownames(expr) <- hbr1$V1
colnames(expr) <- NULL
dim(expr)
head(expr)

# Recall that the samples come from two groups: HBR and UHR
# The first 6 were from HBR, the second 6 from UHR
# Record the samples groups'
# (when we calculate differential expression between the two groups,
#   we need to tell the functions which group each sample belongs to)
data_groups <- c(rep("hbr",6), rep("uhr",6))
data_groups

# Create a DGEList object (Differential Gene Expression) from the data
# Supply the matrix of transcript counts per gene/sample & the group labels for the samples
# The DGEList will be used as input for the differential expression calculation functions
# See ?DGEList for details
d <- DGEList(counts=expr, group=factor(data_groups))
d
dim(d)

# Look at and filter the data first
# Use the cpm (counts per million) function to normalize the transcript counts
head(d$counts)
head(cpm(d))
apply(d$counts, 2, sum) # This calculate the total number of transcript counts per sample
keep <- rowSums(cpm(d)>100) >= 2 # This picks out row numbers of genes that are expressed (cpm > 100) in at least 2 samples
d <- d[keep,] # Keep only the genes expressed (cpm > 100) in at least 2 samples
dim(d) # 642 genes are left
d$samples$lib.size <- colSums(d$counts) # Update the "library size" (the total number of transcripts) for each sample
d$samples

# Normalize the reads and plot them
# We'll use MDS (multidimensional scaling) to plot the samples based on the similarity of their gene expression profiles
# (Think of the MDS plot as analogous to a principal components analysis plot)
d <- calcNormFactors(d) # This stores the scaling factors for normalization in d$samples$norm.factors
d
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:2, pch=20)

# Estimate and plot dispersion
# First plot dispersion (a measure similar to sd/mean) vs log(cpm)
# Since there is a trend in dispersion vs. log(cpm),
# a single dispersion value (red line) doesn't describe all the transcripts very well
# Below we try to find a better fit
d1 <- estimateCommonDisp(d, verbose=T)
names(d1)
d1 <- estimateTagwiseDisp(d1)
names(d1)
d1 <- estimateCommonDisp(d, verbose=T)
d1 <- estimateTagwiseDisp(d1)
plotBCV(d1)

# Estimate and plot again, using a generalize linear model to get better fit
design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)

# Now test for differential expression between the HBR and UHR groups
# de1 is a one-column matrix with gene names as rownames,
#   "-1"=down-regulated; "1"=up-regulated
# Extract the names of the DE genes from de1
et12 <- exactTest(d1, pair=c(1,2)) # compare groups 1 and 2 (HBR, UHR)
topTags(et12, n=10) # Look at the top 10 DE genes
de1 <- decideTests(et12, adjust.method="BH", p.value=0.05) # Find the DE genes; 
summary(de1) # How many are up- or down-regulated?

# Plot log(fold change) vs. log(cpm)
# Significantly differently-expressed genes are shown in red
de1tags12 <- rownames(d1)[as.logical(de1)] 
plotSmear(et12, de.tags=de1tags12)
abline(h = c(-2, 2), col = "blue")

# Clean differential expression data for enrichment analysis
df <- de1 %>% as.data.frame
df <- rownames_to_column(df, var = "Gene")
logfc <- et12$table %>% as.data.table
logfc$Gene <- rownames(et12$table)
df <- df %>% left_join(logfc, by = "Gene")
  
# run enrichment following this reference:
# https://ycl6.github.io/GO-Enrichment-Analysis-Demo/4_enrichR.html

all_genes <- unlist(logfc$Gene)
gene_list <- unlist(df$Gene)
# also have to prune the version from the ends of these
all_genes <- str_replace(all_genes, pattern = "\\.[0-9]+$", replacement = "")
gene_list <- str_replace(gene_list, pattern = "\\.[0-9]+$", replacement = "")
# lists need to be sorted for gseGO
all_genes <- sort(all_genes, decreasing = TRUE)
gene_list <- sort(gene_list, decreasing = TRUE)

db <- "org.Hs.eg.db"

# not all of our genes map back to a symbol
# we will filter these down
mapped <- clusterProfiler::bitr(
  gene_list,
  fromType = "ENSEMBL",
  toType = "SYMBOL", 
  OrgDb = db
)

# view the available databases
dbs <- listEnrichrDbs()
dbs <- dbs[order(dbs$libraryName),]

dbs[grep("GO",dbs$libraryName),]$libraryName

dbs_go <- c(
  "GO_Biological_Process_2025",
  "GO_Cellular_Component_2025",
  "GO_Molecular_Function_2025"
)

enriched_go <- enrichr(
  genes = mapped$SYMBOL,
  databases = dbs_go
)

plotEnrich(enriched_go[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")
plotEnrich(enriched_go[[2]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")
plotEnrich(enriched_go[[3]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")

# look at clinical/phenotypes
phe_dbs <- c(
  "PhenGenI_Association_2021",
  "ClinVar_2025"
)

enriched_phe <- enrichr(
  genes = mapped$SYMBOL,
  databases = phe_dbs
)

plotEnrich(enriched_phe[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")
plotEnrich(enriched_phe[[2]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")

# look at the brain atlas for up and down separately
brain_dbs <- c(
  "Allen_Brain_Atlas_up",
  "Allen_Brain_Atlas_down"
)

enriched_brain <- enrichr(
  genes = mapped$SYMBOL,
  databases = brain_dbs
)

plotEnrich(enriched_brain[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")
plotEnrich(enriched_brain[[2]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")

# generate a volcano plot
df$Gene <- str_replace(df$Gene, pattern = "\\.[0-9]+$", replacement = "")
df <- df %>% 
  mutate(
    NegLogP = -log10(PValue),
    Sig = if_else(
      PValue > 0.05,
      "Not",
      if_else(
        logFC > 0,
        "Up",
        "Down"
      )
    )
  ) %>% 
  left_join(
    mapped %>% rename(Gene = ENSEMBL),
    by = "Gene"
  ) %>% 
  filter(!is.na(SYMBOL)) %>% 
  mutate(
    Label = if_else(
      Sig != "Not",
      SYMBOL,
      NA
    )
  )

df %>% 
  ggplot() +
  aes(
    x = logFC,
    y = NegLogP,
    color = Sig,
    label = Label
  ) +
  geom_point(alpha = 0.5) +
  geom_label(hjust=0, vjust=0, size=3, show.legend = FALSE) +
  theme_minimal() +
  scale_color_manual(
    values = c(
      Up = "red",
      Down = "blue",
      Not = "black"
    )
  )
