##########Re-analysis of publicly available 16s microbiome profiling dataset
###paper: 


library(dada2)



if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("dada2")


install.packages("curl")

library(curl)
library(dada2)
# Folder to save
dir.create("fastq_files", showWarnings = FALSE)


# -----------------------------
# 1. Create a table of Run IDs and Submitter IDs
# -----------------------------
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



#######Downloaded files were trimmed with cutadapter to remove adapters
######WSL Ubunt
#path <- "C:/Users/sunda/Downloads/microbiom"
path <- "C:/Users/sunda/Downloads/microbiom/microbiom_trimmed" #new dir


list.files(path)



fnFs <- sort(list.files(path, pattern = "_R1\\.trimmed\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2\\.trimmed\\.fastq\\.gz$", full.names = TRUE))

sample.names <- sub("_R1\\.trimmed\\.fastq\\.gz$", "", basename(fnFs))
###############removng primer
library(dada2)





filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names



#######Filter and Trim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(240,150),   # adjust based on quality plots
                     maxN = 0, maxEE = c(2,2),
                     truncQ = 2, rm.phix = TRUE,
                     compress=TRUE, multithread=FALSE)

####Learn error rate
errF <- learnErrors(filtFs, nbases=1e7, multithread=FALSE, verbose=TRUE)
errR <- learnErrors(filtRs, nbases=1e7, multithread=FALSE, verbose=TRUE)



dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)


dadaFs[[1]]
########### save created files as of now
save.image("dada2_workspace.RData")
save.image(".RData")

###Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

##Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

######no sequence out of the desired in read

##REMOVE CHIMERAS
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN),
               sapply(dadaRs, getN), 
               sapply(mergers, getN),
               rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

##Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "fastq_files/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "fastq_files/silva_species_assignment_v132.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


install.packages("phyloseq")
## [1] '1.30.0'
library(Biostrings);
library(phyloseq)
library(ggplot2)
theme_set(theme_bw())




samples.out <- rownames(seqtab.nochim)

# Split sample names by underscore
parts <- strsplit(samples.out, "_")

# Extract treatment and replicate info
treatment <- sapply(parts, `[`, 1)
replicate <- sapply(parts, `[`, 2)

# Build a metadata dataframe
samdf <- data.frame(Treatment = treatment,
                    Replicate = replicate,
                    stringsAsFactors = FALSE)

# Set rownames
rownames(samdf) <- samples.out
# Ensure sample_data ordering matches
samdf <- samdf[match(rownames(seqtab.nochim), rownames(samdf)), ]

# Build phyloseq (correct orientation)
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  sample_data(samdf),
  tax_table(taxa)
)


ps

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

sample_variables(ps)

plot_richness(ps, x="Replicate", measures=c("Shannon", "Simpson"), color="Treatment")


# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="Treatment", title="Bray NMDS")



ps.clean <- subset_taxa(ps, Kingdom == "Bacteria" & Family != "Mitochondria" & Order != "Chloroplast")
ps.clean <- prune_taxa(taxa_sums(ps.clean) > 0, ps.clean)  # drop taxa with zero after subsetting

top20 <- names(sort(taxa_sums(ps.clean), decreasing=TRUE))[1:30]
ps.top20 <- transform_sample_counts(ps.clean, function(x) x / sum(x))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Treatment", fill="Family") + facet_wrap(~Treatment, scales="free_x")
save.image(".RData")


##########Normalization, Permova and Differential abundance analysis

library(metagenomeSeq)
library(vegan)
library(DESeq2)


ordu_wu <- ordinate(ps, method="PCoA", distance="wunifrac")

plot_ordination(ps, ordu_wu, color="Treatment") +
  ggtitle("PCoA – Weighted UniFrac") +
  theme_bw()
ps.f <- prune_samples(sample_sums(ps) > 0, ps)
ps.ra <- transform_sample_counts(ps.f, function(x) x / sum(x))
bray.dist <- distance(ps.ra, method = "bray")
ordu <- ordinate(ps.ra, method = "NMDS", distance = bray.dist)



library(vegan)

metadata <- data.frame(sample_data(ps.ra))

adonis2(bray.dist ~ Treatment, data = metadata)

alpha <- estimate_richness(ps.f, measures = c("Shannon", "Simpson", "Chao1", "Observed"))
alpha$Treatment <- sample_data(ps.f)$Treatment  # Add metadata column

boxplot(Shannon ~ Treatment, data = alpha)

###########D


library(phyloseq)
library(ggplot2)
library(dplyr)

# -------------------------
# 1. Agglomerate to Phylum
library(phyloseq)
library(ggplot2)

# Transform to relative abundance
ps_rel <- transform_sample_counts(ps.clean, function(x) x / sum(x))

# Agglomerate to Phylum level
ps_phylum <- tax_glom(ps_rel, taxrank = "Phylum")

# Melt to long format
df_phylum <- psmelt(ps_phylum)

# Simple clean barplot (same as earlier)
relabun <- ggplot(df_phylum, aes(x = Treatment, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 12) +
  labs(y = "Relative Abundance", x = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )
relabun
