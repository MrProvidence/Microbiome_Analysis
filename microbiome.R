###############################################################################
# RE-ANALYSIS OF THE 16S MICROBIOME DATA FROM A RESEARCH PAPER
# PAPER TITLE: Undermining the cry for help: the phytopathogenic fungus
# Verticillium dahliae secretes an antimicrobial effector protein to
# undermine host recruitment of antagonistic Pseudomonas bacteria
# DADA2 WORK FLOW (CLEANED & PUBLISHED-QUALITY OUTPUTS)
###############################################################################

########## PACKAGES USED IN THIS WORK FLOW
library(dada2)
library(curl)
library(Biostrings)
library(phyloseq)
library(ggplot2)
library(metagenomeSeq)
library(vegan)
library(DESeq2)
library(dplyr)
library(ggrepel)
library(pheatmap)
library(cowplot)     # for combined panels
library(tidyr)
library(stringr)

# Publication-style ggplot base
theme_set(theme_bw(base_size = 14))

########## INSTALLATION 
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("dada2","phyloseq","DESeq2","metagenomeSeq")) # uncomment if fresh install

###############################################################################
# ---- CREATE DIRECTORIES FOR ORGANIZATION ---- #
###############################################################################
dir.create("fastq_files", showWarnings = FALSE)
dir.create("Plots", showWarnings = FALSE)
dir.create("Tables", showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)  # extra location for publication fig

###############################################################################
# --------------The original experiement utilized 11 replicates ----------------
##### five will be downloaded here
###############################################################################



###############################
### 1. Create table of Run IDs and Submitter IDs
###############################
samples <- data.frame(
  Run = c(
    "ERR15119940","ERR15119941","ERR15119942","ERR15119943","ERR15119944",
    "ERR15119950","ERR15119951","ERR15119952","ERR15119953","ERR15119954",
    "ERR15119961","ERR15119962","ERR15119963","ERR15119964","ERR15119965"
  ),
  Submitter = c(
    "Mock_1","Mock_2","Mock_3","Mock_4","Mock_5",
    "TO22_1","TO22_2","TO22_3","TO22_4","TO22_5",
    "deltaAv2_1","deltaAv2_2","deltaAv2_3","deltaAv2_4","deltaAv2_5"
  )
)

###############################
# -----------------------------
# 2. Download function for ENA fastq
# -----------------------------
download_fastq <- function(run, submitter) {
  base <- paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=", run,
                 "&result=read_run&fields=fastq_ftp&format=tsv&download=true")
  
  message("Fetching download URLs for: ", run)
  
  tmp <- tempfile()
  download.file(base, tmp, quiet = TRUE)
  
  df <- read.table(tmp, sep="\t", header=TRUE, stringsAsFactors = FALSE)
  urls <- unlist(strsplit(df$fastq_ftp, ";"))
  
  # Loop over R1 and R2 URLs
  for (url in urls) {
    file_name <- basename(url)
    
    # detect F or R
    if (grepl("_1.fastq", file_name)) {
      new_name <- paste0(submitter, "_R1.fastq.gz")
    } else if (grepl("_2.fastq", file_name)) {
      new_name <- paste0(submitter, "_R2.fastq.gz")
    } else {
      next
    }
    
    full_url <- paste0("https://", url)
    
    message("Downloading: ", full_url)
    download.file(full_url, new_name, mode = "wb")
  }
}

# -----------------------------
# 3. Loop to download all files
# -----------------------------
for (i in 1:nrow(samples)) {
  download_fastq(samples$Run[i], samples$Submitter[i])
}


####### Adapter Trimming was carried out with Cutadapt in WSL Ubuntu
#path <- "C:/Users/sunda/Downloads/microbiom"
path <- "C:/Users/sunda/Downloads/microbiom/microbiom_trimmed" #new dir


list.files(path)

##########DADA2 wokflow 

fnFs <- sort(list.files(path, pattern = "_R1\\.trimmed\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2\\.trimmed\\.fastq\\.gz$", full.names = TRUE))
sample.names <- sub("_R1\\.trimmed\\.fastq\\.gz$", "", basename(fnFs))

###############################
### 4. Filtering and trimming
###############################
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen = c(240,150),
  maxN = 0, maxEE = c(2,2),
  truncQ = 2, rm.phix = TRUE,
  compress = TRUE, multithread = FALSE
)

###############################
### 5. Learn error rates
###############################
errF <- learnErrors(filtFs, nbases=1e7, multithread=FALSE, verbose=TRUE)
errR <- learnErrors(filtRs, nbases=1e7, multithread=FALSE, verbose=TRUE)

###############################
### 6. Denoising reads with DADA2
###############################
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)

###############################
### 7. Merge paired-end reads
###############################
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs)

###############################
### 8. Construct sequence table
###############################
seqtab <- makeSequenceTable(mergers)

###############################
### 9. Remove chimeras
###############################
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

###############################
### 10. Track read counts through pipeline
###############################
getN <- function(x) sum(getUniques(x))
track <- cbind(
  out,
  sapply(dadaFs, getN),
  sapply(dadaRs, getN),
  sapply(mergers, getN),
  rowSums(seqtab.nochim)
)
colnames(track) <- c("input","filtered","denoisedF","denoisedR","merged","nonchim")
rownames(track) <- sample.names

###############################
### 11. Assign taxonomy
###############################
taxa <- assignTaxonomy(seqtab.nochim, "fastq_files/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "fastq_files/silva_species_assignment_v132.fa.gz")

###############################
### 12. Build phyloseq object
###############################
samples.out <- rownames(seqtab.nochim)
parts <- strsplit(samples.out, "_")

treatment <- sapply(parts, `[`, 1)
replicate <- sapply(parts, `[`, 2)

samdf <- data.frame(
  Treatment = treatment,
  Replicate = replicate,
  stringsAsFactors = FALSE
)
rownames(samdf) <- samples.out

ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  sample_data(samdf),
  tax_table(taxa)
)

dna <- DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

###############################
### 13. Clean phyloseq (remove mitochondria & chloroplast)
###############################
ps.clean <- subset_taxa(ps, Kingdom=="Bacteria" & Family!="Mitochondria" & Order!="Chloroplast")
ps.clean <- prune_taxa(taxa_sums(ps.clean) > 0, ps.clean)

###############################
### 14 Alpha diversity
###############################
alpha <- estimate_richness(ps.clean, measures = c("Shannon", "Simpson", "Chao1", "Observed"))
alpha$Treatment <- sample_data(ps.clean)$Treatment
alpha$Replicate <- sample_data(ps.clean)$Replicate
write.csv(alpha, "Tables/AlphaDiversity_Table.csv", row.names = TRUE)

p_alpha <- plot_richness(ps.clean, x="Treatment", measures=c("Shannon","Simpson"), color="Treatment")
ggsave("Plots/AlphaDiversity_Shannon_Simpson.png", p_alpha, width=7, height=5)

p_shannon <- ggplot(alpha, aes(x = Treatment, y = Shannon, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  theme_bw() +
  labs(title = "Shannon Diversity Index", y = "Shannon Index", x = "")
ggsave("Plots/Shannon_Boxplot.png", p_shannon, width=7, height=5)

###############################
### 15. Beta diversity: NMDS + Bray-Curtis
###############################
ps.prop <- transform_sample_counts(ps.clean, function(x) x / sum(x))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
p_nmds <- plot_ordination(ps.prop, ord.nmds.bray, color="Treatment", title="Bray-Curtis NMDS") +
  theme_bw()
ggsave("Plots/Bray_NMDS.png", p_nmds, width=7, height=5)

###############################
### 16. Relative abundance barplot (PHYLA)
###############################
ps_rel <- transform_sample_counts(ps.clean, function(x) x / sum(x))
ps_phylum <- tax_glom(ps_rel, taxrank="Phylum")
df_phylum <- psmelt(ps_phylum)
relabun <- ggplot(df_phylum, aes(x=Treatment, y=Abundance, fill=Phylum)) +
  geom_bar(stat="identity") +
  theme_bw() +
  labs(y="Relative Abundance", x="") +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("Plots/RelativeAbundance_Phylum.png", relabun, width=8, height=5)

##############################################################################
# ---- 17. Save important tables ---- #
##############################################################################
write.csv(seqtab.nochim, "Tables/ASV_table_nochim.csv")
write.csv(as.data.frame(taxa), "Tables/taxonomy_table.csv")
write.csv(samdf, "Tables/sample_metadata.csv")
write.csv(track, "Tables/read_tracking.csv")

###############################################################################
# ---- 18. Beta statistical tests & PERMANOVA ----
###############################################################################
library(vegan)
pss <- ps.clean
sample_data(pss)$Treatment <- factor(sample_data(pss)$Treatment)

ps.prop <- transform_sample_counts(pss, function(x) x / sum(x))
dist_bray <- phyloseq::distance(ps.prop, method="bray")
samdf <- data.frame(sample_data(pss))

set.seed(42)
adonis_bray <- adonis2(dist_bray ~ Treatment, data=samdf, permutations=999)
write.csv(as.data.frame(adonis_bray), "Tables/PERMANOVA_bray.csv")

bd <- betadisper(dist_bray, samdf$Treatment)
anova_bd <- anova(bd)
permutest_bd <- permutest(bd)
sink("Tables/PERMDISP_bray.txt"); print(anova_bd); print(permutest_bd); sink()

###############################################################################
# ---- 19. DESeq2 differential abundance 
###############################################################################

library(DESeq2)

# extract count table
countdata <- as.data.frame(otu_table(ps.clean))
if (taxa_are_rows(ps.clean)) {
  countdata <- t(countdata)
}

# metadata
meta <- as.data.frame(sample_data(ps.clean))

# ensure Treatment is a factor AND make TO22 the reference
meta$Treatment <- factor(meta$Treatment, levels = c("TO22", "deltaAv2"))

### create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = meta,
  design = ~ Treatment
)

dds <- DESeq(dds)

res <- results(dds, contrast = c("Treatment", "TO22", "deltaAv2"))

levels(meta$Treatment)


res_df <- as.data.frame(res)
res_df$ASV <- rownames(res_df)

##-------------------------------
## ADD TAXONOMY
##-------------------------------

tax <- as.data.frame(tax_table(ps.clean))
tax$ASV <- rownames(tax)

merged <- merge(res_df, tax, by = "ASV", all.x = TRUE)


###############################################################################
# ---- ORDER-LEVEL DRIVERS (color by Phylum) ----
###############################################################################

order_fc <- merged %>%
  group_by(Order, Phylum) %>%
  summarize(
    log2FC = median(log2FoldChange, na.rm = TRUE),
    padj   = min(padj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(Order))

# Apply threshold
drivers <- order_fc %>%
  filter(abs(log2FC) >= 1.5) %>%
  arrange(log2FC)

drivers$Order <- factor(drivers$Order, levels = drivers$Order)

### ORDER PLOT
p_order <- ggplot(drivers, aes(x = log2FC, y = Order, color = Phylum)) +
  geom_vline(xintercept = 0, size = 1) +
  geom_point(size = 4) +
  scale_x_continuous(limits = c(-2.5, 2.5), breaks = seq(-2, 2, 1)) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_blank()
  ) +
  xlab("log2 fold change") +
  ggtitle("ΔAv2 vs WT(TO22)— Order")


###############################################################################
# ---- GENUS-LEVEL DRIVERS (color by Phylum) ----
###############################################################################

genus_fc <- merged %>%
  filter(Order %in% drivers$Order) %>%      # only genera in driver Orders
  group_by(Genus, Order, Phylum) %>%
  summarize(log2FC = median(log2FoldChange, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(Genus))

genus_fc$Genus <- factor(genus_fc$Genus, levels = genus_fc$Genus)

### GENUS PLOT
p_genus <- ggplot(genus_fc, aes(x = log2FC, y = Genus, color = Phylum)) +
  geom_vline(xintercept = 0, size = 1) +
  geom_point(size = 4) +
  scale_x_continuous(limits = c(-2.5, 2.5), breaks = seq(-2, 2, 1)) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right"
  ) +
  xlab("log2 fold change") +
  ggtitle("ΔAv2 vs WT(TO22) — Genus")


###############################################################################
# ---- SAVE TABLES ----
###############################################################################

write.csv(res_df,     "Tables/DESeq2_full_results.csv", row.names = FALSE)
write.csv(order_fc,   "Tables/Order_level_results.csv",  row.names = FALSE)
write.csv(genus_fc,   "Tables/Genus_level_results.csv",  row.names = FALSE)

###############################################################################
# ---- PRINT PLOTS ----
###############################################################################

print(p_order)
print(p_genus)






###############################################################################
# ---- 22. Combine all plots into a single multipanel figure  ----
###############################################################################
# Collect plot objects
plot_list <- list(
  "Alpha_Shannon" = p_shannon,
  "NMDS_Bray" = p_nmds,
  "RelAbun_Phylum" = relabun,
  "Order_Drivers" = p_order,
  "Genus_Drivers" = p_genus,
  "DESeq2_Volcano" = p_vol
)

library(cowplot)

##############################
### PAGE 1 — COMMUNITY STRUCTURE
##############################

# Adjust margin to prevent label collision
plot_list$RelAbun_Phylum <- plot_list$RelAbun_Phylum +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 40))

page1 <- plot_grid(
  plot_list$Alpha_Shannon,
  plot_list$NMDS_Bray,
  plot_list$RelAbun_Phylum,
  ncol = 3,
  labels = "AUTO",
  label_size = 14,
  rel_widths = c(1, 1, 1.3)
)

### PDF ###
pdf("Plots/Page1_CommunityStructure.pdf", width = 18, height = 6)
print(page1)
dev.off()

### PNG (FIXED) ###
png("Plots/Page1_CommunityStructure.png", width = 3600, height = 1200, res = 300)
print(page1)
dev.off()


##############################
### PAGE 2 — DIFFERENTIAL ABUNDANCE
##############################

page2 <- plot_grid(
  plot_list$Order_Drivers,
  plot_list$Genus_Drivers,
  ncol = 2,
  labels = "AUTO",
  label_size = 14,
  rel_widths = c(1, 1)
)

### PDF ###
pdf("Plots/Page2_DA.pdf", width = 18, height = 6)
print(page2)
dev.off()

### PNG (FIXED — print was missing) ###
png("Plots/Page2_DA.png", width = 3600, height = 1200, res = 300)
print(page2)   # ← THIS LINE WAS MISSING
dev.off()




###############################################################################
# ---- 23. Save final key tables ----
###############################################################################
write.csv(order_fc, "Tables/Order_level_results.csv", row.names = FALSE)
write.csv(genus_fc, "Tables/Genus_level_results.csv", row.names = FALSE)
write.csv(res_df, "Tables/DESeq2_full_results.csv", row.names = TRUE)

###############################################################################
# ---- 24. End of pipeline ----
###############################################################################





















