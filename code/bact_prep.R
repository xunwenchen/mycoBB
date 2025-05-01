# Construct phyloseq object and plot relative abundance bar chart ----

# Analyzing bacterial data using zotu (ASV) table by unoise3----
# Note: To calculate alpha diversity, a non-normalized zotu table should be better. 
# For PCoA, a normalized zotu table is preferred. 
# Read: https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html

# This analysis requires the packages in my_packages. Install them and proceed. 
# Install and load packages, set plot theme, and source own functions ----
# install 'Rtools', search website, download, and install
# install 'phyloseq' and 'DESeq2'

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install("DESeq2")

# install 'pairwiseAdonis'
install.packages("htmlwidgets") # needed for 'devtools'
install.packages("devtools") # use the package panel of Rstudio to install
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis") # require 'devtools'


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("apeglm")

remotes::install_github("mlverse/chattr")
devtools::install_github("zdk123/SpiecEasi")



my_packages <- c('readxl', 'phyloseq', 'ape', 'Biostrings', 'tidyverse', 'dplyr', 'pairwiseAdonis', 'reshape2', 'ggplot2', 'ggpubr', 'DESeq2', 'GGally', 'igraph', 'vegan', 'ecodist', 'agricolae', 'stats', 'writexl', 'dendextend', 'apeglm', 'pheatmap', 'GOplot','tibble', 'boot', 'chattr', 'SpiecEasi')

install.packages(my_packages) # 'phyloseq', 'DESeq2', and 'pairwiseAdonis' may not be installed using this code. But it is fine as they have been installed above. my_packages is also used for library attaching as follows.

 
# START -------------
# If all required packages are installed, start from here

# Clear the environment
rm(list = ls())

# Clear the console
cat("\014")  # This is equivalent to pressing Ctrl+L

# Clear all plots
if (!is.null(dev.list())) dev.off()

# save or not
save <- FALSE

# Load packages needed
my_packages <- c('readxl', 'phyloseq', 'ape', 'Biostrings', 'tidyverse', 'dplyr', 'pairwiseAdonis', 'reshape2', 'ggplot2', 'ggpubr', 'DESeq2', 'GGally', 'igraph', 'vegan', 'ecodist', 'agricolae', 'stats', 'writexl', 'dendextend', 'apeglm', 'pheatmap', 'GOplot','tibble', 'boot', 'chattr', 'SpiecEasi')

lapply(my_packages, library, character.only = TRUE) # load my_packages


theme_set(theme_bw()) # Set plot theme of ggplot2. It can be changed to other themes.
source('code/fun.R') # load own function(s)

# check versions of R and other essential packages
R.Version()
packageVersion('phyloseq') # version 1.44.0
packageVersion('DESeq2') # version 1.40.2
packageVersion('vegan') # version 2.6.4

# if use chattr for integrated ChatGPT, run the following
chattr_app()

# ~ Load raw data ----
# prepare 5 files according to required format before proceeding:
# 1. zotu table
# 2. tax info
# 3. tree file
# 4. metadata
# 5. seq data

# The following zotu raw table is from Zhengsheng Yu (co-author), i.e., zotutab_raw.txt. This file contains sample ID used during sequencing, zotu IDs, and reads. It is not scaled/normalized. This file has been edited in excel to change the sample ID to recognized ones, and the reads were scaled/normalized using the smallest read count. The zotutab_raw.txt was obtained using codes in the Supplementary Material "Suppl_Mater_unoise3".


# after obtaining the above 5 files in hands, proceed as follows
# ~~ 1. load zotu table ----
zotu_table_zhengsheng <- read_xlsx("data/zotu_table_raw.xlsx") # Zhengsheng's raw zotu table - un-rarefied. "data/zotu_table_raw.xlsx" was obtained using unoise3.

# python codes are shown as follows:
####################################
####################################


# python codes here


####################################
####################################
zotu_table_mengying <- read_xlsx("data/zotu_table.xlsx") # Mengying's zotu table - rarefied. Rarefied zotu table can be obtained using un-rarefied one by rarefy_even_depth() but note counts can be different every time and some data can be discarded due to the nature of rarefaction (random resampling).

# The un-rarefied table will be used to calculate alpha-diversity, while the rarefied one will be used to analyze beta-diversity

# check if the otu table has been normalized
see_if_norm <- zotu_table_zhengsheng %>% 
  summarise(across(hps_ck_1:hp_f_cd_8, ~ sum(.x, na.rm = TRUE)))
see_if_norm # not normalized 
min(see_if_norm); max(see_if_norm) # check the smallest read count
rm(see_if_norm) # remove the temp file 'see_if_norm'

see_if_norm <- zotu_table_mengying %>% 
  summarise(across(hps_ck_1:hp_f_cd_8, ~ sum(.x, na.rm = TRUE)))
see_if_norm # not normalized 
min(see_if_norm); max(see_if_norm) # check the smallest read count
rm(see_if_norm) # remove the temp file 'see_if_norm'


# ~~ 2. load tax info ----
tax <- read_xlsx('data/tax.xlsx') 

# ~~ 3. load tree file ----
tree <- read.tree('data/zotu_tree.nwk') # ape package is needed

# ~~ 4. load metadata ----
metadata <- read_xlsx('data/metadata.xlsx')

# ~~ 5. load seqs data ----
# 'dna-seq.fasta' is extracted from 'zotu_rep-seqs.qza') # obtained using qiime2
# The format is like >Zotu1
# ATTGGACAATGGGCGCAAGCCTGATCCAGCCATGCCG....
# There are 20861 zotus in total
seq <- readDNAStringSet('data/dna-seq.fasta', format="fasta",
                        nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)


# ~~ Prepare zotu table matrix ----
# should be matrix, rownames are zotu IDs, colnames are sample IDs
# it is required if need to use phyloseq package.
# convert to matrix
# get the column zotu_ID, then assign as row names
zotu_table_zhengsheng2 <- column_to_rownames(zotu_table_zhengsheng, var = "Zotu_ID")
zotu_table_mengying2 <- column_to_rownames(zotu_table_mengying, var = "Zotu_ID")

# convert data frame to matrix
zotu_zs_mat <- data.matrix(zotu_table_zhengsheng2)
zotu_my_mat <- data.matrix(zotu_table_mengying2)

rm(zotu_table_zhengsheng, zotu_table_mengying,
   zotu_table_zhengsheng2, zotu_table_mengying2) # remove temp files. 


# ~~ Prepare tax file matrix ----
# dealing with tax info
# get the column zotu_ID, then assign as row names
tax2 <- column_to_rownames(tax, var = "Zotu_ID")
# convert data frame to matrix
# use as.matrix rather than data.matrix
# otherwise, text will be deleted
tax_mat <- as.matrix(tax2) 
rm(tax2) # remove temp file 'tax2'
rm(tax) # remove temp file 'tax'

# ~~ Prepare metadata dataframe ----
metadata <- as.data.frame(metadata)

metadata_df <- column_to_rownames(metadata, var = "Sample_ID")

# convert things in metadata_df to factors
metadata_df$Cd <- as.factor(metadata_df$Cd)
metadata_df$Compartment <- as.factor(metadata_df$Compartment)
metadata_df$Mycorrhizal <- as.factor(metadata_df$Mycorrhizal)
metadata_df$Spatial_order <- as.factor(metadata_df$Spatial_order)

# set factor order. The later plotting shall follow the set order
# Assuming metadata_df is your dataframe
metadata_df$Group <- factor(metadata_df$Group, levels = c("Ctrl", "Cd", "M", "M+Cd"))

metadata_df$Compartment <- factor(metadata_df$Compartment, levels = c("Endosphere", "Rhizoplane", "Rhizosphere", "Hyphosphere", "Hyphae", "Bulk soil"))
metadata_df$Mycorrhizal <- factor(metadata_df$Mycorrhizal, levels = c("NM", "M"))
# Check the levels of the different factors
levels(metadata_df$Group)
levels(metadata_df$Compartment)
levels(metadata_df$Mycorrhizal)

rm(metadata) # remove temp file


# ~~ Prepare tree and seqs data ----
# no need to prepare(?)


# CHECK POINT 1 ---------------------------------------------------------------
# 1. zotu_mat - has otu id and sample id -- OK ******** it should be a matrix
head(zotu_zs_mat)
head(zotu_my_mat)
class(zotu_zs_mat)
class(zotu_my_mat)

# 2. tax_mat - has otu id and taxon info   ******** it should be a matrix
head(tax_mat)
class(tax_mat)

# 3. metadata - has sample id and treatments -- OK ******** it should be a data.frame
head(metadata_df)
class(metadata_df)


# 4. tree - has otu id -- OK(?) how to check(?)
taxa_names(tree)
class(tree)

# 5. seq - has DNA sequences and otu id -- OK(?) how to check(?)
taxa_names(seq)
class(seq)


### the above 5 items can be combined as a phylo object for further analysis


# IMPORTANT: telling how to combine the 5 items #### 
ZOTU_zhengsheng <- otu_table(zotu_zs_mat, taxa_are_rows = TRUE) # IMPORTANT step ***
rm(zotu_zs_mat)

ZOTU_mengying <- otu_table(zotu_my_mat, taxa_are_rows = TRUE) # IMPORTANT step ***
rm(zotu_my_mat)

TAX <- tax_table(tax_mat) # IMPORTANT step ***
rm(tax_mat)

metadata <- sample_data(metadata_df) # IMPORTANT step ***
rm(metadata_df) # remove temp file

# TREE <- phy_tree(tree) # actually no need, seems 'tree' =? 'TREE'

# using the following to check type

class(ZOTU_zhengsheng)
class(ZOTU_mengying)
class(TAX)
class(metadata)
class(tree)
class(seq)

# Now in the 'Global Environment' panel, there should be 5 objects: ZOTU, TAX, metadata, tree, seq

# Combine to construct phyloseq object ----------------------------------------
pl_un_rarefied <- phyloseq(ZOTU_zhengsheng, TAX, metadata, tree, seq)
pl <- phyloseq(ZOTU_mengying, TAX, metadata, tree, seq)

# if return error saying taxa/OTU names are not match. Use the followings to check
# a <- taxa_names(ZOTU)
# b <- taxa_names(TAX)
# c <- taxa_names(seq)
# d <- taxa_names(tree)
# setdiff(a, b); setdiff(a, c); setdiff(a, d) etc.

# seems probably the taxa_names of tree are different from other phyloseq objects
# zotu vs Zotu; but "Zotu1" in ZOTU vs "'Zotu1'" in tree seems OK. For example, check:
# setdiff(taxa_names(ZOTU),taxa_names(tree))

# see this post: https://github.com/joey711/phyloseq/issues/1044


# CHECK POINT 2 ---------------------------------------------------------------
# inspect different names and variables

pl_un_rarefied
sample_names(pl_un_rarefied)
rank_names(pl_un_rarefied)
sample_variables(pl_un_rarefied)
taxa_names(pl_un_rarefied)
refseq(pl_un_rarefied)

pl
sample_names(pl)
rank_names(pl)
sample_variables(pl)
taxa_names(pl)
refseq(pl)

# check again the min and max read counts
# If the data is not rarefied, min != max.
print(paste('min =', min(sample_sums(pl_un_rarefied)))); print(paste('max =', max(sample_sums(pl_un_rarefied)))) # min != max -- un-rarefied

print(paste('min =', min(sample_sums(pl)))); print(paste('max =', max(sample_sums(pl)))) # min = max -- rarefied

# remove unnecessary files
rm(metadata, tree, ZOTU_zhengsheng, ZOTU_mengying, seq)

# now there should be only two phyloseq objects left for further analyses. 
# check
pl_un_rarefied
pl
# OK

# PLOT ------------------------------------------------------------------------

# ~ Overall bacterial relative abundance vs all samples ----
# We can rarefy the 'pl_un_rarefied', but the resulting otu table will be different from 'pl' due to the nature of rarefaction. We use 'pl' instead of rarefying the 'pl_un_rarefied.' This is not a critical issue. 


# Plot the bar plot for all phyla
# ***Note it can take time to run depending on the computer spec.***
plot_bar(pl, "sample_name", fill="Phylum")+
  geom_bar(stat="identity", color = NA, width = 1)


# the bar plot was exported and edited in AI/CorelDraw to add group names and labels for axes. Relative zotu count was converted to % in AI/CorelDraw. 

# ~ calculate actual relative abundance at the phylum level ----

# Below code snippet demonstrate how to achieve this.
# The arrange(), rename() and select() are from the dplyr package and spread() is from the tidyr package, both packages are part of the tidyverse.
# The select() and spread() are used to convet the output from long to wide format.
# ref: https://github.com/joey711/phyloseq/issues/1521


# pl_relabun_phylum <- pl %>% 
#   tax_glom(taxrank = "Phylum") %>% 
#   transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>%
#   select(Phylum, Sample, Abundance) %>% 
#   spread(Sample, Abundance)

# subset different compartments
# endosphere
# endo <- subset_samples(pl, Compartment == 'Endosphere')
# sample_names(endo)
# plot_bar(endo, "Phylum", fill="Phylum", facet_grid = Mycorrhizal~Cd)+
#   geom_bar(aes(color="Phylum", fill="Phylum"), stat="identity", position="stack")

# Plot top phylum ----
# ***Note this step can be time consuming ***
# Transform to relative abundances
ps.rel <- transform_sample_counts(pl, function(x) x/sum(x))

# Melt to a data frame. ***a bit time-consuming***
df.m <- psmelt(ps.rel)

# Calculate the total abundance for each Phylum
df.sum <- df.m %>% group_by(Phylum) %>% summarise(Abundance = sum(Abundance))

# Get the top 20 Phyla
top10 <- df.sum %>% arrange(desc(Abundance)) %>% head(10) %>% pull(Phylum)

# Create a new column in the original data frame. If the Phylum is in the top 10, use the Phylum name. Otherwise, label as 'Others'.
df.m$phylum.top10 <- ifelse(df.m$Phylum %in% top10, df.m$Phylum, 'Others')

# Check if 'other' assigned
unique(df.m$phylum.top10)

# Now you can plot with ggplot2

colors <- c("Proteobacteria" = '#8dd3c7',
            "Firmicutes" = '#ffffb3',
            "Bacteroidetes"  = '#bebada',
            "Actinobacteria" = '#fb8072',
            "Planctomycetes" = '#80b1d3',
            "Patescibacteria" = '#fdb462',
            "Gemmatimonadetes"= '#b3de69',
            "Acidobacteria" = '#fccde5',
            "Verrucomicrobia" = '#d9d9d9',
            "Chloroflexi"='#bc80bd',
            "Others" = 'darkgray')

df.m$phylum.top10 <- factor(df.m$phylum.top10, levels = c("Proteobacteria",
                                                      "Firmicutes",
                                                      "Bacteroidetes",
                                                      "Actinobacteria",
                                                      "Planctomycetes",
                                                      "Patescibacteria",
                                                      "Gemmatimonadetes",
                                                      "Acidobacteria",
                                                      "Verrucomicrobia",
                                                      "Chloroflexi",
                                                      "Others"))

# ~~~~~ Fig. S. Rel abundance of bacterial comm (every rep) ----
rel_ab <- ggplot(df.m, aes(x = sample_name, y = Abundance, fill = phylum.top10)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = colors, breaks = c("Proteobacteria",
                                                "Firmicutes",
                                                "Bacteroidetes",
                                                "Actinobacteria",
                                                "Planctomycetes",
                                                "Patescibacteria",
                                                "Gemmatimonadetes",
                                                "Acidobacteria",
                                                "Verrucomicrobia",
                                                "Chloroflexi",
                                                "Others")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
rel_ab

if(save){
  ggsave(filename = "out/Fig. 1. Rel_ab_by_compt.jpg", plot = rel_ab, 
       width = 7680, height = 4320 , units = "px", dpi = 600)
}

# create rel ab plot using average value of 8 replicates
# create a new column for correct grouping, 6 compartments, 4 treatment groups, and 1 compartment contain only 2 groups, so 6*4-2 = 22 groups for the plot in total
df.m$Compartment_and_Group <- paste(df.m$Compartment, df.m$Group, sep = "_")
# check if 22 groups
length(unique(df.m$Compartment_and_Group)) == 22
# calculate average value for each group
df.m_avg <- df.m %>% 
  group_by(Compartment_and_Group, Phylum, phylum.top10, Compartment, Group) %>% 
  summarise(avg_rel_ab = sum(Abundance)/8)
length(unique(df.m_avg$Compartment_and_Group)) == 22

# set factor levels in the order of:
# for compartment like: endosphere, rhizoplane, rhizosphere, hyphosphere, hyphae, bulk soil
# for group like: Ctrl, Cd, M, M+Cd
df.m_avg$Compartment_and_Group <- factor(df.m_avg$Compartment_and_Group, 
  levels = c(
    "Endosphere_Ctrl", "Endosphere_Cd", "Endosphere_M", "Endosphere_M+Cd",
    "Rhizoplane_Ctrl", "Rhizoplane_Cd", "Rhizoplane_M", "Rhizoplane_M+Cd",
    "Rhizosphere_Ctrl", "Rhizosphere_Cd", "Rhizosphere_M", "Rhizosphere_M+Cd",
    "Hyphosphere_Ctrl", "Hyphosphere_Cd", "Hyphosphere_M", "Hyphosphere_M+Cd",
    "Hyphae_M", "Hyphae_M+Cd",
    "Bulk soil_Ctrl", "Bulk soil_Cd", "Bulk soil_M", "Bulk soil_M+Cd"))


# ~~~~~ Fig. 1a Rel abundance of bacterial comm (average of 8 reps) -----------
rel_ab_avg <- ggplot(df.m_avg, aes(x = Group, y = avg_rel_ab, fill = phylum.top10)) +
  geom_bar(stat = 'identity') +
  # facet wrap by Compartment
  facet_wrap(~Compartment, ncol = 6)+
  scale_fill_manual(values = colors, breaks = c("Proteobacteria",
                                                "Firmicutes",
                                                "Bacteroidetes",
                                                "Actinobacteria",
                                                "Planctomycetes",
                                                "Patescibacteria",
                                                "Gemmatimonadetes",
                                                "Acidobacteria",
                                                "Verrucomicrobia",
                                                "Chloroflexi",
                                                "Others")) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))+
  xlab("Treatment")+
  ylab("Relative abundance")+
  # make y-axis to be 100 not 1, # remove % in tick numbers
  scale_y_continuous(labels = scales::percent_format(scale = 100))+
  # remove grid lines
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  # change legend title to Top 10 Phyla
  labs(fill = "Top 10 Phyla")
  
rel_ab_avg

# END ----