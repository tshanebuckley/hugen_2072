#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#install.packages("biomaRt")
#BiocManager::install("GWENA")

library(GWENA)
library(magrittr) # Not mandatory, we use the pipe `%>%` to ease readability.

threads_to_use <- 2

# Q1. explain how the next line deals with the format difference from the example file. 
# Hint 1: what is the purpose of function t that is added to the example file. 
# Hint 2: you can try the next line without function t to see the difference.
UHRBHR_expr = t(read.table("./UHRBHRExpressions.csv", sep=',', header=TRUE, row.names=1))
ncol(UHRBHR_expr)
nrow(UHRBHR_expr)
UHRBHR_expr[1:4,1:4]
is_data_expr(UHRBHR_expr)

UHRBHR_traits = read.table("UHRBHR_traits.txt", header=TRUE, row.names=1)
UHRBHR_traits
unique(UHRBHR_traits$Condition)

# Q2. explain how the following module selects expressed genes to further anlayze. 
UHRBHR_expr_filtered <- UHRBHR_expr[,colSums(UHRBHR_expr)>1]

#Q3. explain what the next line does for the given gene names and why this is needed. 
# Hint: run the subsequent lines without line 29 and 30 to see the difference, especially in calling bio_enrich
genes<-colnames(UHRBHR_expr_filtered)
genes.adj <- substr(genes, 1, 15)
colnames(UHRBHR_expr_filtered) <- genes.adj

# Remaining number of genes
ncol(UHRBHR_expr_filtered)
UHRBHR_expr_filtered<-UHRBHR_expr_filtered[,1:200]
net <- build_net(UHRBHR_expr_filtered, cor_func = "spearman", n_threads = threads_to_use)

net$metadata$power
fit_power_table <- net$metadata$fit_power_table
fit_power_table[fit_power_table$Power == net$metadata$power, "SFT.R.sq"]

modules <- detect_modules(UHRBHR_expr_filtered, 
                          net$network, 
                          detailled_result = TRUE,
                          merge_threshold = 0.25)
length(unique(modules$modules_premerge))
length(unique(modules$modules))
layout_mod_merge <- plot_modules_merge(
  modules_premerge = modules$modules_premerge, 
  modules_merged = modules$modules)

ggplot2::ggplot(data.frame(modules$modules %>% stack), 
                ggplot2::aes(x = ind)) + ggplot2::stat_count() +
  ggplot2::ylab("Number of genes") +
  ggplot2::xlab("Module")

#Q4. Can you explain the enrichment result shown on the resulting plot 
enrichment <- bio_enrich(modules$modules, organism="hsapiens")
plot_enrichment(enrichment)
