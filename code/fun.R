plot_box <- function(df, x, y, z){
  ggplot2::theme_set(theme_bw())
  
  plot <- ggplot2::ggplot(df, aes(x = {{x}}, y = {{y}}), fill = {{z}})+
    geom_boxplot()+
    # geom_jitter(alpha = 0.5, width = 0.2)+
    stat_summary(fun=mean, geom="point", 
                 shape=20, size=3, color="indianred")
  
  return(plot)
}


plot_w_unifrac <- function(phylo_ob, by){
  library(phyloseq)
  ordu <-  ordinate(phylo_ob, "PCoA", "unifrac", weighted = TRUE)
  plot_ordination(phylo_ob, ordu, color=by, shape=by)+
    geom_point(size=2, alpha=0.5)+
    
    scale_shape_manual(values = c(16, 15, 17, 18)) +
    scale_color_manual(values = c("#bababa","#404040", "#f4a582", "#ca0020"))+
  
    #scale_colour_brewer(type="qual", palette="Set1")+
    ggtitle("PCoA on weighted UniFrac distance")
  
}


plot_w_unifrac_cpmt <- function(phylo_ob, by){
  library(phyloseq)
  ordu <-  ordinate(phylo_ob, "PCoA", "unifrac", weighted = TRUE)
  plot_ordination(phylo_ob, ordu, color=by, shape=by)+
    geom_point(size=2, alpha=0.5)+
    
    scale_shape_manual(values = c(1, 0, 16, 15, 18)) +
    scale_color_manual(values = c('sandybrown','dodgerblue','indianred4','hotpink','darkgrey'))+
    
    #scale_colour_brewer(type="qual", palette="Set1")+
    ggtitle("PCoA on weighted UniFrac distance")
  
}

plot_w_unifrac_cpmt_hy <- function(phylo_ob, by){
  library(phyloseq)
  ordu <-  ordinate(phylo_ob, "PCoA", "unifrac", weighted = TRUE)
  plot_ordination(phylo_ob, ordu, color=by, shape=by)+
    geom_point(size=2, alpha=0.5)+
    
    scale_shape_manual(values = c(1, 0, 16, 15, 17, 18)) +
    scale_color_manual(values = c('sandybrown','dodgerblue','indianred4','hotpink','purple', 'darkgrey'))+
    
    #scale_colour_brewer(type="qual", palette="Set1")+
    ggtitle("PCoA on weighted UniFrac distance")
  
}

plot_PCoA_bray <- function(phylo_ob, by){
  library(phyloseq)
  ordu <-  ordinate(phylo_ob, method = "PCoA", distance = "bray")
  plot_ordination(phylo_ob, ordu, color=by, shape=by)+
    geom_point(size=2, alpha=0.5)+
    scale_shape_manual(values = c(16, 15, 17, 18)) +
    scale_color_manual(values = c("#dfc27d","#a6611a", "#80cdc1", "#018571"))+
    #scale_colour_brewer(type="qual", palette="Set1")+
    ggtitle("PCoA on Bray-Curtis distance")
  
}


plot_uw_unifrac <- function(phylo_ob, by){
  library(phyloseq)
  ordu_uw = ordinate(phylo_ob, "PCoA", "unifrac", weighted = FALSE)
  plot_ordination(phylo_ob, ordu_uw, color=by, shape=by)+
  geom_point(size=2, alpha=0.5)+
    scale_shape_manual(values = c(16, 15, 17, 18)) +
    scale_color_manual(values = c("#dfc27d","#a6611a", "#80cdc1", "#018571"))+
  # scale_colour_brewer(type="qual", palette="Set1")+
  ggtitle("PCoA on unweighted UniFrac distance")
}


# # The following not working
# adonis <- function(phylo_ob, by){
#   dist <- phyloseq::distance(phylo_ob, method ='wunifrac')
#   library(pairwiseAdonis)
#   result <- pairwise.adonis(dist, sample_data(phylo_ob)$by, perm = 999)
#   
#   
# }


# # The following not working
# plot_uw_unifrac2 <- function(phylo_ob, var1, x, var2, y, by){
#   library(phyloseq)
#   subset1 <- subset_samples(phylo_ob, {{var1}} == x)
#   subset2 <- subset_samples(subset1, {{var2}} == y)
#   
#   ordu_uw = ordinate(subset2, "PCoA", "unifrac", weighted = FALSE)
#   plot_ordination(subset2, ordu_uw, color=by, shape=by)+
#     geom_point(size=5, alpha=0.5)+
#     scale_colour_brewer(type="qual", palette="Set1")+
#     ggtitle("PCoA on unweighted UniFrac distance")
# }


## Create function
## function to plot within-group beta diversity distance
beta_boxplot <- function(physeq, method = "bray", group) {
  
  # physeq: phyloseq-class object
  # method: beta-diversity metric. Default "bray", i.e., Bray-Curtis dissimilarity 
  # group: factorial variable to group
  
  ## Packages
  require("phyloseq") # v.1.30.0
  require("ggplot2") # v.3.3.2
  
  ## Identify the correspondence: group and samples
  group2samp <- list() # list to save the correspondence between group <--> samples
  group_list <- as.factor(get_variable(sample_data(physeq), group)) # list of group elements
  for (groups in levels(group_list)) { # loop over the no. of group levels
    target_group <- which(group_list == groups) # vct pos of the curr group variable 
    group2samp[[ groups ]] <- sample_names(physeq)[target_group] # matching samples: based on vct pos
  }  
  
  ## Calculate beta-diversity
  beta_div_dist <- phyloseq::distance(physeq = physeq, method = method)
  beta_div_dist <- as(beta_div_dist, "matrix")
  
  
  
  
  
  ## Coerce distance mtx into a tidy data frame by group
  dist_df <- data.frame() # save results in df 
  counter <- 1 
  for (groups in names(group2samp)) { # loop over group fct levels 
    sub_dist <- beta_div_dist[ group2samp[[groups]], group2samp[[groups]] ] # subset dist mtx to curr samples
    #print(sub_dist)
    no_samp_col <- ncol(sub_dist) # n cols: curr sub dist
    no_samp_row <- nrow(sub_dist) # n rows: curr sub dist
    for ( cols in seq(no_samp_col) ) { # loop over cols: curr sub_dist
      if ( cols > 1 ) {
        for ( rows in seq((cols-1)) ) { # loop over rows: curr sub_dist 
          ## Save results
          dist_df[ counter, "sample_pair" ] <- paste0( colnames(sub_dist)[cols], "-",  
                                                       rownames(sub_dist)[rows] ) # sample pair
          dist_df[ counter, "group" ] <- groups # group  
          dist_df[ counter, "beta_div_method" ] <- method # method
          dist_df[ counter, "beta_div_value" ] <- sub_dist[rows, cols] # beta-diversity for the sample pair     
          counter = counter + 1
        }
      }
    }
  }
  
  ## Create a ggplot2 boxplot
  ## set levels of grpup as 'Ctrl', 'Cd', 'M', 'M+Cd'
  dist_df$group <- factor(dist_df$group, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))
  plot_boxplot <- ggplot(data = dist_df, aes(x = group, y = beta_div_value)) + 
    geom_boxplot() + #geom_boxplot(outlier.shape=NA)
    geom_jitter(alpha = 0.5, width = 0.2) + 
    theme_bw() + 
    xlab(group) + ylab(method) 
    # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  ## Save df and boxplot into a list 
  # list_Out <- list("data" = dist_df, "plot" = plot_boxplot) 
  # 
  # return(list_Out)
}


# function to plot between-group beta diversity comparison

beta_boxplot_btw <- function(phylo_ob, dist_method, pair1, pair2){
  
  # calculate distance matrix
  dist <- phyloseq::distance(phylo_ob, method = dist_method)
  dist <- as.matrix(dist)
  
  # assign '0' to upper triangle of the dist matrix
  dist[upper.tri(dist)] <- 0
  
  # using melt() in 'reshape2' package to convert dist matrix to dataframe in common form and delete previously assigned '0'.
  library(reshape2)
  dist_melt <- subset(melt(dist), value!=0)
  
  # assign Sample_ID from metadata file to dist dataframe - dist_melt
  # Var1 and Var2 are column names of dist_melt, Sample_ID is a column name in metadata (load metadata first if have not)
  metadata.in.fun <- read_xlsx('data/metadata.xlsx')
  library(dplyr)
  temp <- dplyr::left_join(dist_melt, metadata.in.fun, by = c ('Var1'='Sample_ID'))
  temp2 <- dplyr::left_join(temp, metadata.in.fun, by = c ('Var2'='Sample_ID'))
  
  # generate new columns, and exclude unwanted columns
  temp3 <- temp2 %>% dplyr::transmute(distance = value, Compartment = Compartment.x, Pair = paste0(Group.x, sep = "_VS_", Group.y))
  
  # see what pairs can be used for comparison
  groups_to_select <- temp3 %>% group_by(Pair) %>% summarise()
  print('select two pairs for comparison:')
  print(groups_to_select)
  # plot 
  plot <- temp3 %>% 
    filter(Pair == pair1 | Pair == pair2) %>% 
    ggplot(aes(x=Pair, y = distance))+
    geom_boxplot()+
    geom_jitter(alpha = 0.5, width = 0.2)
  return(plot)
  
}


log2fc <- function(phylo_ob, Group, trt_group, ref_group){
  library("DESeq2")
  
  # since an error shows "In DESeqDataSet(se, design = design, ignoreRank): some variables in design formula are characters, converting to factors". The variables in design formula ("~ Group" in our case) are converted to factors using the following as.factor() function
  sample_data(phylo_ob)$Group <- as.factor(sample_data(phylo_ob)$Group)

  
  # The following two lines actually do all the complicated DESeq2 work. The function phyloseq_to_deseq2 converts the phyloseq-format microbiome data (pl_hps in our case) into a DESeqDataSet with dispersions estimated, using the experimental design formula, also shown (the ~Group term in our case). The DESeq function does the rest of the testing, in this case with default testing framework, but we can actually use alternatives. 
  deseq2_temp <- phyloseq_to_deseq2(phylo_ob, ~ Group)
  
  deseq2_temp2  <-  DESeq(deseq2_temp, test="Wald", fitType="parametric")
  
  
  # The following results function call creates a table of the results of the tests. Very fast. The hard work was already stored with the rest of the DESeq2-related data in our latest version of the deseq2_temp2 object (see above). I then order by the adjusted p-value, removing the entries with an NA value (seems we did not remove them). The rest of this example is just formatting the results table with taxonomic information for nice(ish) display in the HTML output.
  
  res <- results(deseq2_temp2, 
                 contrast=c(Group,trt_group,ref_group), 
                 cooksCutoff = FALSE)
  res # use this to check whether it is M+Cd vs M, or M vs M+Cd:
  # Group M.Cd vs M 
  # Wald test p-value: Group M.Cd vs M 
  # DataFrame with 20861 rows and 6 columns
  
  
  
  
  alpha = 0.01
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phylo_ob)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  dim(sigtab)
  
  # Let's look at the OTUs that were significantly different between the two groups (M vs M+Cd in our case? I think so) The following makes a nice ggplot2 summary of the results.
  
  library("ggplot2")
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  # Phylum level
  # x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  # x = sort(x, TRUE)
  # sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  
  # class level
  
  
  # Order level
  # x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
  # x = sort(x, TRUE)
  # sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))
  
  # Family level
  x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))
  
  
  # Genus level
  # x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
  # x = sort(x, TRUE)
  # sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
  
  # the following x also need to be changed 
  
  # point plot 
  ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) +
    geom_point(size=3, alpha = 0.5)
  
  
  # heatmap plot
  # ggplot(sigtab, aes(x = Family, y = Phylum, fill = log2FoldChange)) +
  #   geom_tile() +
  #   scale_fill_gradient(low = "#4870B1", high = "#D73735")  # Adjust the color scale as needed
    
  # + 
  #   theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
}



log2fc_genus <- function(phylo_ob, Group, trt_group, ref_group){
  library("DESeq2")
  
  # since an error shows "In DESeqDataSet(se, design = design, ignoreRank): some variables in design formula are characters, converting to factors". The variables in design formula ("~ Group" in our case) are converted to factors using the following as.factor() function
  sample_data(phylo_ob)$Group <- as.factor(sample_data(phylo_ob)$Group)
  
  
  # The following two lines actually do all the complicated DESeq2 work. The function phyloseq_to_deseq2 converts the phyloseq-format microbiome data (pl_hps in our case) into a DESeqDataSet with dispersions estimated, using the experimental design formula, also shown (the ~Group term in our case). The DESeq function does the rest of the testing, in this case with default testing framework, but we can actually use alternatives. 
  deseq2_temp <- phyloseq_to_deseq2(phylo_ob, ~ Group)
  
  deseq2_temp2  <-  DESeq(deseq2_temp, test="Wald", fitType="parametric")
  
  
  # The following results function call creates a table of the results of the tests. Very fast. The hard work was already stored with the rest of the DESeq2-related data in our latest version of the deseq2_temp2 object (see above). I then order by the adjusted p-value, removing the entries with an NA value (seems we did not remove them). The rest of this example is just formatting the results table with taxonomic information for nice(ish) display in the HTML output.
  
  res <- results(deseq2_temp2, 
                 contrast=c(Group,trt_group,ref_group), 
                 cooksCutoff = FALSE)
  res # use this to check whether it is M+Cd vs M, or M vs M+Cd:
  # Group M.Cd vs M 
  # Wald test p-value: Group M.Cd vs M 
  # DataFrame with 20861 rows and 6 columns
  
  
  
  
  alpha = 0.01
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phylo_ob)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  dim(sigtab)
  
  # Let's look at the OTUs that were significantly different between the two groups (M vs M+Cd in our case? I think so) The following makes a nice ggplot2 summary of the results.
  
  library("ggplot2")
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  # Phylum level
  # x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  # x = sort(x, TRUE)
  # sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  
  # class level
  
  
  # Order level
  # x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
  # x = sort(x, TRUE)
  # sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))
  
  # Family level
  # x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
  # x = sort(x, TRUE)
  # sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))
  
  
  # Genus level
  x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
  
  # the following x also need to be changed 
  
  # point plot 
  ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
    geom_point(size=3, alpha = 0.5)
  
  
  # heatmap plot
  # ggplot(sigtab, aes(x = Family, y = Phylum, fill = log2FoldChange)) +
  #   geom_tile() +
  #   scale_fill_gradient(low = "#4870B1", high = "#D73735")  # Adjust the color scale as needed
  
  # + 
  #   theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
}


# Function to perform hierarchical clustering based on the bacterial tree
dendro_counts <- function(phylo_ob, Group, trt_group, ref_group) {
  
  # load packages
  library(DESeq2)
  library(phyloseq)
  library(ggplot2)
  library(dendextend)
  library(ape)
  
  # Prepare the data for DESeq2
  sample_data(phylo_ob)$Group <- as.factor(sample_data(phylo_ob)$Group)
  deseq2_temp <- phyloseq_to_deseq2(phylo_ob, ~ Group)
  
  # Perform DESeq2 analysis
  deseq2_temp2 <- DESeq(deseq2_temp, test = "Wald", fitType = "parametric")
  res <- results(deseq2_temp2, contrast = c(Group, trt_group, ref_group), cooksCutoff = FALSE) # use Cooks distance to remove outliers or not
  
  alpha <- 0.01
  sigtab <- res[which(res$padj < alpha), ]
  sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(phylo_ob)[rownames(sigtab), ], "matrix"))
  
  # Apply DESeq2's variance stabilization transformation
  variance_stab_counts <- assay(varianceStabilizingTransformation(deseq2_temp2, blind = TRUE))
  rownames(variance_stab_counts) <- rownames(sigtab)
  
  # Extract the bacterial tree from the phyloseq object
  phy_tree <- phy_tree(phylo_ob)
  
  # Create a hierarchical clustering dendrogram based on the phylogenetic tree
   tree_dend <- as.dendrogram(as.hclust(phy_tree$edge[, 2:1])) 
  # By reversing the order of the columns with [, 2:1], you’re essentially flipping the direction of each edge in the tree.
  
  
  # Reorder the variance-stabilized counts based on the dendrogram
  variance_stab_counts_ordered <- variance_stab_counts[, order.dendrogram(tree_dend)]
  
  # Create a heatmap using ggplot2 with hierarchical clustering
  ggplot(data = melt(variance_stab_counts_ordered), aes(x = Family, y = variable, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "#4870B1", high = "#D73735") +
    labs(title = "Heatmap of Variance Stabilized Counts by Family and Sample",
         x = "Family",
         y = "Sample") +
    theme_minimal() +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}






flat_mt <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    Source = rownames(cormat)[row(cormat)[ut]],
    Target = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}



metabo_otu_cor_p_tab <- function(df_otu, df_metabo, which_sample, to_which_sample, cor_p_value, cor_R_value){
  # For each treatment in a compartment, only those OTUs detected in all 8 samples were retained. These more-stringent OTU tables were used to correlate with the metabolomic profiles. 
  # omit rows contain 0 otu number
  # df2 <- df_otu[!rowSums(df_otu==0),]
  df2 <- df_otu
  # omit rows with group-max value smaller than 5
  df2$rowMax <- apply(df2, 1, max, na.rm=TRUE) # find row max value assign to column 'rowMax'. apply(df,1,...) --- '1' means manipulation is performed on rows. '2' is on columns
  
  df3 <- df2[df2$rowMax >= 15, ]     # select only rows with rowMax more than x zotu counts (the larger the more stringent)
  df4 <- select(df3, -rowMax) # drop the rowMax column
  
  # remove rows with >=5 occurrence of 0 zotu count
  df4$count.0 <- apply(df4, 1, function(x) length(which(x==0))) # count the no. of occurrence of 0
  df5 <- df4[df4$count.0 <= 0, ]  # retain only those equal or less than x (the smaller the more stringent in this case)
  df6 <- select(df5, -count.0) # drop the count.0 column
  
  # transpose and convert to dataframe
  df7 <- as.data.frame(t(as.data.frame(df6)))
  
  
  # delete compound name which is not needed later
  metabo_df2 <- df_metabo[,-2]
  # transpose
  metabo_df3 <- t(metabo_df2)
  # take the first row and assign it as the colname
  colnames(metabo_df3) <- metabo_df3[1, ]
  # remove first row which is not needed (repeated)
  metabo_df4 <- metabo_df3[-1, ]
  # select a treatment group in the same compartment (e.g., Rhizosphere and Ctrl). Each compartment and each treatment has a network
  metabo_df5 <- as.data.frame(metabo_df4[c(which_sample:to_which_sample),])
  # check if the df is numeric, if not, need to conver to numeric as follow
  metabo_df6 <- as.data.frame(lapply(metabo_df5, as.numeric))
  
  # construct correlation matrix. This matrix is large. Maybe not needed
  # temp_cor_mt <- cor(df3, metabo_df6, 'pearson', use = "complete.obs")
  
  # construct R and P value tables
  library(Hmisc)
  # cor_p_mt is a very large file
  cor_p_mt <- rcorr(as.matrix(df7), as.matrix(metabo_df6), type = 'spearman')
  
  # flattern the tables using the function flat_mt() (own function)
  
  # using the function
  # flat_cor_p_tab is a large file
  flat_cor_p_tab <- flat_mt(cor_p_mt$r, cor_p_mt$P)
  
  
  library(stringr)
  # select only zotu paired with compound (not zotu vs zotu, or com vs com)
  zotu_com_tab <- flat_cor_p_tab %>% 
    filter(str_detect(Source, 'Zotu') & str_detect(Target, 'Com') | str_detect(Source, 'Com') & str_detect(Target, 'Zotu'))
  
  zotu_com_tab_sig <- zotu_com_tab[abs(zotu_com_tab$cor) >= cor_R_value &
                                     zotu_com_tab$p <= cor_p_value, ]
  
  return(zotu_com_tab_sig)
}


construct_node_file <- function(edge_raw, TAX, metabo_df){
  # edge_raw is the melted/flattern zotu-com correlation with p and r value; TAX is the phylo object of tax info. metabo_df is the pos and neg combined metabo file with all treatments and samples. 
  node_raw <- 
    as.data.frame(
      unique(
        c(edge_raw$Source, edge_raw$Target)))
  
  # assign "Name" as column name for later merging.
  colnames(node_raw)[1]  <- "Name" 
  
  # prepare for zotu_tax_tab file
  zotu_tax_tab <- as.data.frame(TAX)
  zotu_tax_tab <- rownames_to_column(zotu_tax_tab, "Name") # here used 'Name' is for later merging with node file 
  
  
  # prepare for metabo file
  metabo_Name_com_name <- metabo_df %>% select(com_id, com_name) %>% rename(Name = com_id)
  
  # first merge
  
  node_temp1 <- merge(node_raw, zotu_tax_tab, by = 'Name', all.x = TRUE)
  
  # last merge
  
  node_temp2 <- merge(node_temp1, metabo_Name_com_name, by = 'Name', all.x = TRUE)
  
  
  # assign a column 'ID' and order before 'Name'. It is needed for later VLookup in excel to create Source and Target numbers, to be used in gephi. 
  
  node_final <- node_temp2 %>% mutate(ID = 1:n(), .before = 'Name')
  
  
  return(node_final)
  
}





node_merge_otu <- function(node_file_address){
  
  node <- read_csv(node_file_address)
  zotu_tax_tab <- as.data.frame(TAX)
  zotu_tax_tab2 <- rownames_to_column(zotu_tax_tab, "Name") # here used 'Name' is for later merging with node file
  node2 <- merge(node, zotu_tax_tab2, by = 'Name', all.x = TRUE)
  node3 <- node2 %>% mutate(ID = 1:n(), .before = 'Name')
  return(node3)
  
}



plot_points_with_line <- function(df, x, y){
  
  ggplot2::theme_set(theme_bw())
  
  plot <- ggplot2::ggplot(df, aes(x = {{x}}, y = {{y}}), group=1)+
    # add group = 1 here is needed to create line connecting the points
    geom_line(color = 'grey', size = 1.5)+
    geom_point(shape=21, color="black", fill="#69b3a2", size=5)
  
  return(plot)
  
}




plot_pi_zi <- function(df){
  
  ggplot(df, aes(Pi, Zi))+
    geom_point(size = 2, alpha = 0.4)+
    geom_hline(yintercept=2.5, linetype="dashed", color = "darkgrey")+
    geom_vline(xintercept=0.62, linetype="dashed", color = "darkgrey")+
    xlab('Among-module connectivity (Pi)')+
    ylab('Within-module connectivity (Zi)')+
    xlim(0, 0.8)+
    ylim(0, 4)
  # set pi = 0.62, zi = 2.5
}




plot_node_degree <- function(df1, df2, df3, df4, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'Ctrl') 
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'Cd') 
  c <- df3 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M')
  d <- df4 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M+Cd')
  
  df_list <- list(a, b, c, d)
  df_merge <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
  df_merge$trt <- factor(df_merge$trt, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))
  ggplot(df_merge, aes({{tax_level}}, Accum_degree, color = trt, group = trt))+
    # geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)
    # + coord_flip()
  
}



plot_node_degree_2_trt <- function(df1, df2, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M')
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M+Cd')
  
  df_list <- list(a, b)
  
  
  df_merge <- Reduce(function(x, y) merge(x, y,
                                           all = T), df_list)
  
  ggplot(df_merge, aes({{tax_level}}, Accum_degree, color = trt, group = trt))+
    # geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)
    # + coord_flip()
  
}



plot_node_degree0 <- function(df1, df2, df3, df4, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'Ctrl') 
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'Cd') 
  c <- df3 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M')
  d <- df4 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M+Cd')
  
  df_list <- list(a, b, c, d)
  df_merge <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
  df_merge$trt <- factor(df_merge$trt, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))
  ggplot(df_merge, aes({{tax_level}}, Accum_degree, color = trt, group = trt))+
    geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)+
    coord_flip()
  
}



plot_node_degree_2_trt0 <- function(df1, df2, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M')
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M+Cd')
  
  df_list <- list(a, b)
  
  
  df_merge <- Reduce(function(x, y) merge(x, y,
                                          all = T), df_list)
  
  ggplot(df_merge, aes({{tax_level}}, Accum_degree, color = trt, group = trt))+
    geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)+ 
    coord_flip()
  
}



plot_btw0 <- function(df1, df2, df3, df4, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'Ctrl') 
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'Cd') 
  c <- df3 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'M')
  d <- df4 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'M+Cd')
  
  df_list <- list(a, b, c, d)
  df_merge <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
  df_merge$trt <- factor(df_merge$trt, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))
  ggplot(df_merge, aes({{tax_level}}, Betweenness, color = trt, group = trt))+
    geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)+
    coord_flip()
  
}

plot_btw_2_trt0 <- function(df1, df2, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'M')
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'M+Cd')
  
  df_list <- list(a, b)
  
  
  df_merge <- Reduce(function(x, y) merge(x, y,
                                          all = T), df_list)
  
  ggplot(df_merge, aes({{tax_level}}, Betweenness, color = trt, group = trt))+
    geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)+ 
    coord_flip()
  
}



plot_btw <- function(df1, df2, df3, df4, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'Ctrl') 
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'Cd') 
  c <- df3 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'M')
  d <- df4 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'M+Cd')
  
  df_list <- list(a, b, c, d)
  df_merge <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
  df_merge$trt <- factor(df_merge$trt, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))
  ggplot(df_merge, aes({{tax_level}}, Betweenness, color = trt, group = trt))+
    # geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)
  # + coord_flip()
  
}



plot_btw_2_trt <- function(df1, df2, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'M')
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'M+Cd')
  
  df_list <- list(a, b)
  
  
  df_merge <- Reduce(function(x, y) merge(x, y,
                                          all = T), df_list)
  
  ggplot(df_merge, aes({{tax_level}}, Betweenness, color = trt, group = trt))+
    # geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)
  # + coord_flip()
  
}





plot_module <- function(df1, df2, df3, df4, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'Ctrl') 
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'Cd') 
  c <- df3 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'M')
  d <- df4 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'M+Cd')
  
  df_list <- list(a, b, c, d)
  df_merge <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
  df_merge$trt <- factor(df_merge$trt, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))
  ggplot(df_merge, aes({{tax_level}}, Number_of_module, color = trt, group = trt))+
    # geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)
  # + coord_flip()
  
}


plot_module_2_trt <- function(df1, df2, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'M')
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'M+Cd')
  
  df_list <- list(a, b)
  
  
  df_merge <- Reduce(function(x, y) merge(x, y,
                                          all = T), df_list)
  
  ggplot(df_merge, aes({{tax_level}}, Number_of_module, color = trt, group = trt))+
    # geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)
  # + coord_flip()
  
}


plot_module0 <- function(df1, df2, df3, df4, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'Ctrl') 
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'Cd') 
  c <- df3 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'M')
  d <- df4 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'M+Cd')
  
  df_list <- list(a, b, c, d)
  df_merge <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
  df_merge$trt <- factor(df_merge$trt, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))
  ggplot(df_merge, aes({{tax_level}}, Number_of_module, color = trt, group = trt))+
    geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)+
    coord_flip()
  
}


plot_module_2_trt0 <- function(df1, df2, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'M')
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'M+Cd')
  
  df_list <- list(a, b)
  
  
  df_merge <- Reduce(function(x, y) merge(x, y,
                                          all = T), df_list)
  
  ggplot(df_merge, aes({{tax_level}}, Number_of_module, color = trt, group = trt))+
    geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)+
    coord_flip()
  
}


heatmap_vst <- function(phylo_obj, taxrank, trt_group, ref_group){
  
  sample_data(phylo_obj)$Group <- factor(sample_data(phylo_obj)$Group)
  phylo_obj_glom <- tax_glom(phylo_obj, taxrank = taxrank, NArm = TRUE)

# check if tax_glom() actually summing up the counts of different OTUs affiliated with the same taxonomic group. Note if NArm = TRUE (default), the following checking may fail. But should be okay. Total OTU counts of a sample should be the same before and after the tax_glom(), if no 'NA' removed. 

# sum(as.data.frame(otu_table(phylo_obj_glom))$rp_f_1)-
#   sum(as.data.frame(otu_table(phylo_obj))$rp_f_1)

# the number shows the number of counts removed.


deseq2_data <- phyloseq_to_deseq2(phylo_obj_glom, ~ Group)

# Run DESeq2 analysis
deseq2_results <- DESeq(deseq2_data)


# Extract log2FoldChange, p-values, and adjusted p-values
results_table <- results(deseq2_results, contrast=c('Group', trt_group, ref_group), cooksCutoff = FALSE)

# Extract ASV (zotu) counts using variance stabilization transformation
# can use vst() here instead of varianceStabilizingTransformation(), but note:
# In summary, while both functions perform a variance stabilizing transformation, vst() is optimized for speed by only considering a subset of variables to estimate the dispersion trend. https://www.biostars.org/p/459013/
variance_stab_counts <- assay(varianceStabilizingTransformation(deseq2_results, blind = TRUE))


vst_counts_t <- t(variance_stab_counts) # transpose first as scale works column-wise


vst_counts_t_scale <- scale(vst_counts_t) # scale column-wise (accross diff. samples)

variance_stab_counts_scale <- t(vst_counts_t_scale) # transpose back


# Combine the following tables

# class(results_table)
# class(variance_stab_counts)
# class(tax_table(phylo_obj_glom))

# Combine the following tables accoring to Zotu IDs which are rownames
# tax_table(pl_subset_glom) # [Zotu IDs, tax]
# variance_stab_counts # VST counts #[Zotu IDs, sample IDs]
# results_table # [Zotu IDs, log2fc related parameters like 'pvalue' and 'padj']

# Convert the tax_table to a data frame for merging
tax_df <- as.data.frame(tax_table(phylo_obj_glom))
# Create a new column named 'Zotu_ID' with the Zotu IDs
tax_df$Zotu_ID <- rownames(tax_df)

# Convert the variance_stab_counts matrix to a data frame
variance_stab_df <- as.data.frame(variance_stab_counts_scale)

# Create a new column named 'Zotu_ID' with the Zotu IDs
variance_stab_df$Zotu_ID <- rownames(variance_stab_df)

# Convert the results_table to a data frame
results_df <- as.data.frame(results_table)
# Create a new column named 'Zotu_ID' with the Zotu IDs
results_df$Zotu_ID <- rownames(results_df)


# Merge the tables based on the 'Zotu_ID' column
combined_table <- tax_df %>%
  left_join(variance_stab_df, by = 'Zotu_ID') %>%
  left_join(results_df, by = 'Zotu_ID')

# select only those padj < 0.05 and remove padj = NA
subset_combined_table <- combined_table[combined_table$padj < 0.05 & !is.na(combined_table$padj), ]

melted_data <- melt(subset_combined_table, 
                    id.vars = c(taxrank), 
                    measure.vars = c(8:23)) # value is VST count

# revert the sample order
melted_data$variable <- factor(melted_data$variable, levels = rev(levels(melted_data$variable)))

# plot heatmap
ggplot(melted_data, aes(x = melted_data[,1], y = variable, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colours=rev(c('#b2182b','#ef8a62','#fddbc7','#d1e5f0','#67a9cf','#2166ac'))) + # gradientn -- not gradient
  labs(x = taxrank, y = "Sample") +
  # theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

}

pheatmap_vst_scaled <- function(phylo_obj, trt_group, ref_group, low_count_threshold, top_x){

  # Remove low-count Zotu. The threshold low_count_threshold e.g. 10 can be adjusted. Those count sums more than 10 will be kept. 
  phylo_obj_filtered <- prune_taxa(taxa_sums(phylo_obj) > low_count_threshold, phylo_obj)
  
  
  # taxa_filtered <- tax_table(phylo_obj_filtered)
  
  
  deseq2_data <- phyloseq_to_deseq2(phylo_obj_filtered, ~ Group) # use filtered and non-glom data
  
  # Run DESeq2 analysis
  deseq2_results <- DESeq(deseq2_data)
  
  # save normalized read counts
  # nor_counts <- counts(deseq2_results, normalized = TRUE)
  # head(nor_counts)
  # write.csv(nor_counts, 'out/pl_unra_rzs_filt_DESeq_nor_counts.csv')
  
  
  # Extract log2FoldChange, p-values, and adjusted p-values
  results_table <- results(deseq2_results, contrast=c("Group", trt_group, ref_group), cooksCutoff = FALSE, alpha = 0.05) # default alpha = 0.1 
  # results_table
  # summary(results_table)
  
  
  res_tab_ordered <- results_table[order(results_table$padj), ]
  
  
  
  # res_tab_ordered
  
  top_zotu <- res_tab_ordered[1:top_x, ]
  top_zotu2 <- row.names(top_zotu)
  # top_zotu2
  
  # plotMA(results_table, cex = 0.7, ylim = c(-10, 10))
  # abline(h=c(-1, 1), col = 'indianred', lwd = 3)
  
  # plotMA with lfcShrink
  # resultsNames(deseq2_results)
  # resLFC <- lfcShrink(deseq2_results, coef = 'Group_M.Cd_vs_M', type = 'apeglm')
  # plotMA(resLFC, cex = 0.7, ylim = c(-10, 10))
  # abline(h=c(-1, 1), col = 'indianred', lwd = 3)
  
    # plotDispEsts(dds, main = 'Dispersion plot')

  # Extract ASV (zotu) counts using variance stabilization transformation
  # can use vst() here instead of varianceStabilizingTransformation(), but note:
  # In summary, while both functions perform a variance stabilizing transformation, vst() is optimized for speed by only considering a subset of variables to estimate the dispersion trend. https://www.biostars.org/p/459013/
  variance_stab_counts <- assay(varianceStabilizingTransformation(deseq2_results, blind = TRUE))
  
  vst_counts_t <- t(variance_stab_counts) # transpose first as scale works column-wise
  
  
  vst_counts_t_scale <- scale(vst_counts_t) # scale column-wise (accross diff. samples)
  
  variance_stab_counts_scale <- t(vst_counts_t_scale) # transpose back
  
  
  # pheatmap(t(variance_stab_counts[top_zotu2, ]),
  #          cluster_rows = FALSE,
  #          cluster_cols = TRUE)
  # 
  pheatmap(t(variance_stab_counts_scale[top_zotu2, ]),
           cluster_rows = FALSE,
           cluster_cols = TRUE)
  

}
  
  
  
pheatmap_vst <- function(phylo_obj, trt_group, ref_group, low_count_threshold, top_x){
  
  # Remove low-count Zotu. The threshold low_count_threshold e.g. 10 can be adjusted. Those count sums more than 10 will be kept. 
  phylo_obj_filtered <- prune_taxa(taxa_sums(phylo_obj) > low_count_threshold, phylo_obj)
  
  
  # taxa_filtered <- tax_table(phylo_obj_filtered)
  
  
  deseq2_data <- phyloseq_to_deseq2(phylo_obj_filtered, ~ Group) # use filtered and non-glom data
  
  # Run DESeq2 analysis
  deseq2_results <- DESeq(deseq2_data)
  
  # save normalized read counts
  # nor_counts <- counts(deseq2_results, normalized = TRUE)
  # head(nor_counts)
  # write.csv(nor_counts, 'out/pl_unra_rzs_filt_DESeq_nor_counts.csv')
  
  
  # Extract log2FoldChange, p-values, and adjusted p-values
  results_table <- results(deseq2_results, contrast=c("Group", trt_group, ref_group), cooksCutoff = FALSE, alpha = 0.05) # default alpha = 0.1 
  # results_table
  # summary(results_table)
  
  
  res_tab_ordered <- results_table[order(results_table$padj), ]
  
  
  
  # res_tab_ordered
  
  top_zotu <- res_tab_ordered[1:top_x, ]
  top_zotu2 <- row.names(top_zotu)
  # top_zotu2
  
  # plotMA(results_table, cex = 0.7, ylim = c(-10, 10))
  # abline(h=c(-1, 1), col = 'indianred', lwd = 3)
  
  # plotMA with lfcShrink
  # resultsNames(deseq2_results)
  # resLFC <- lfcShrink(deseq2_results, coef = 'Group_M.Cd_vs_M', type = 'apeglm')
  # plotMA(resLFC, cex = 0.7, ylim = c(-10, 10))
  # abline(h=c(-1, 1), col = 'indianred', lwd = 3)
  
  # plotDispEsts(dds, main = 'Dispersion plot')
  
  # Extract ASV (zotu) counts using variance stabilization transformation
  # can use vst() here instead of varianceStabilizingTransformation(), but note:
  # In summary, while both functions perform a variance stabilizing transformation, vst() is optimized for speed by only considering a subset of variables to estimate the dispersion trend. https://www.biostars.org/p/459013/
  variance_stab_counts <- assay(varianceStabilizingTransformation(deseq2_results, blind = TRUE))
  
  vst_counts_t <- t(variance_stab_counts) # transpose first as scale works column-wise
  
  
  vst_counts_t_scale <- scale(vst_counts_t) # scale column-wise (accross diff. samples)
  
  variance_stab_counts_scale <- t(vst_counts_t_scale) # transpose back
  
  
  pheatmap(t(variance_stab_counts[top_zotu2, ]),
           cluster_rows = FALSE,
           cluster_cols = TRUE)

  # pheatmap(t(variance_stab_counts_scale[top_zotu2, ]),
  #          cluster_rows = FALSE,
  #          cluster_cols = TRUE)
  
  
}


str_wk <- function(pl_ob, gl_tax, target_group, c.cutoff){
  
  # Agglomerate target taxa
  pl.gl <- tax_glom(pl_ob, taxrank = gl_tax, )
  
  # Extract the OTU table and convert to mx/df
  otu.tab.gl <- otu_table(pl.gl)
  ot.mx <- as.data.frame(otu.tab.gl)
  
  # transpose matrix
  ot.mx.t <- t(ot.mx)
  
  # subset target group
  tg <- ot.mx.t[grepl(target_group, rownames(ot.mx.t)), ]
  
  # filter OTUs with cutoff
  tg.f <- tg[, apply(tg, 2, function(col) any(col >= c.cutoff))]
  
  # # OR
  # Subset columns based on column sum (sum count of an OTU) threshold
  # c.cutoff = 10
  # hp.cd.f <- hp.cd[, colSums(hp.cd) >= c.cutoff]
  
  # Initialize an empty dataframe to store the results
  rs <- data.frame()
  
  # Linear regression fitting for each pair, using loop.
  # This step can be time-consuming.
  # use n*(n-1)/2 to calculate number of pairs
   for(i in 1:(ncol(tg.f)-1)){
    for(j in (i+1):ncol(tg.f)){
      fit <- lm(tg.f[,i] ~ tg.f[,j]) # y = a*x+b
      a <- coef(fit)[2]  # get the slope a. 
      b <- coef(fit)[1]  # get the intercept b. 
      r_squared <- summary(fit)$r.squared
      p_value <- summary(fit)$coefficients[2,4]
      cor_coefficient <- cor(tg.f[,i], tg.f[,j]) # default method is Pearson
      # for spearman: cor(tg.f[,i], tg.f[,j], method = "spearman")
      
      # Calculate the angle in degrees based on slope a
      angle <- atan(a) * (180 / pi)
      
      # Add the results to the dataframe
      rs <- rbind(rs, data.frame(avs_n = colnames(tg.f)[i], 
                                 asv_m = colnames(tg.f)[j], 
                                 slope_a = a, 
                                 intercept_b = b, 
                                 R_squared = r_squared, 
                                 p_value = p_value, 
                                 cor_coeff = cor_coefficient, 
                                 slope_angle = angle))
      
      # Print progress update
      if (i %% 10 == 0) {
        cat("Processing pair", i, "-", j, "\n")
      }
    }
  }
  

  return(rs)

}



plot_sw <- function(rs){

 
    g <- deparse(substitute(rs))
  
  
par(mfrow=c(1,3),mar=c(4,4,2,1))

hist(rs$slope_angle, breaks = seq(-100, 100, by = 10), ylim = c(0, 18000), 
     main = g,
     xlab = "Slope angle of linear regression line")

hist(rs$cor_coeff, breaks = seq(-1, 1, by = 0.05), ylim = c(0, 13000),
     main = g,
     xlab = "Pearson's rho")

hist(rs$p_value, breaks = seq(0, 1, by = 0.01), ylim = c(0, 2800),
     main = g,
     xlab = "p-value of Pearson's corr.")


}

plot_sw_sig <- function(rs){
  
  # Subset the dataframe
  rs.sig <- rs[rs$p_value <= 0.05, ]
  
  g <- deparse(substitute(rs))
  
  
  par(mfrow=c(1,3),mar=c(4,4,2,1))
  
  hist(rs.sig$slope_angle, breaks = seq(-100, 100, by = 10), ylim = c(0, 2000), 
       main = g,
       xlab = "Slope angle of linear regression line")
  
  hist(rs.sig$cor_coeff, breaks = seq(-1, 1, by = 0.05), ylim = c(0, 2000),
       main = g,
       xlab = "Pearson's rho")
  
  hist(rs.sig$p_value, breaks = seq(0, 0.05, by = 0.01), ylim = c(0, 2500),
       main = g,
       xlab = "p-value of Pearson's corr.")
  
  
}


sig.p.perc <- function(rs){
  
  sig.p <- sum(rs$p_value<=0.05)
  all.p <- length(rs$p_value)
  
  perc <- round(sig.p*100/all.p, 2)
  
  print(paste('Perc of p < 0.05 = ', perc, '%'))
  
}





angle.perc <- function(rs){
  
  a30 <- sum(rs$slope_angle > 0 & rs$slope_angle <= 30)
  a60 <- sum(rs$slope_angle > 30 & rs$slope_angle <= 60)
  a90 <- sum(rs$slope_angle > 60 & rs$slope_angle <= 90)
  
  ao30 <- sum(rs$slope_angle >= -30 & rs$slope_angle <= 0)
  ao60 <- sum(rs$slope_angle >= -60 & rs$slope_angle < -30)
  ao90 <- sum(rs$slope_angle >= -90 & rs$slope_angle < -60)
  
  check <- sum(a30, a60, a90,
               ao30, ao60, ao90) == length(rs$slope_angle)
  
  if (check) {
    print("Count sum = total number of slope angles. Result in %:")
  }
  
  if (!check) {
    stop("Count sum and number of slope angles not the same. Pls. check data.")
  }
  
  angle_group <- c(a30, a60, a90, ao30, ao60, ao90)
  angle_name <- c('a30', 'a60', 'a90', 'ao30', 'ao60', 'ao90')
  perc.rs <- numeric(length(angle_group)) 
  
  for(i in seq_along(angle_group)) {
    
    
    perc.rs[i] <- angle_group[i]*100/length(rs$slope_angle)
  }
  
  
  a.perc.tab <- as.data.frame(rbind(angle_name, perc.rs))
  a.perc.tab.t <- t(a.perc.tab)
  a.perc.tab.t.df <- as.data.frame(a.perc.tab.t)
  # colnames(a.perc.tab) <- a.perc.tab[1, ]
  # a.perc.tab <- a.perc.tab[-1, ]
  
  return(a.perc.tab.t.df)
  
  
}


# Create the pie chart
plot_pie <- function(df){
  df$angle_name <- factor(df$angle_name, levels = c("ao30", "ao60", "ao90", 'a90', 'a60', 'a30'))
  # df$angle_name <- factor(df$angle_name, levels = c("ao30", "ao60", "ao90", 'a30', 'a60', 'a90'))
  df$perc.rs <- as.numeric(df$perc.rs)
  
  custom_colors <- c("a30" = "#DFF2E9", "a60" = "#A7DBC3", "a90" = "#53B98B",
                     "ao30" = "#FBDEE4", "ao60" = "#F5B1BF", "ao90" = "#EF839B")
  
  
  ggplot(df, aes(x = '', y = perc.rs, fill = angle_name)) +
    # geom_bar(stat = 'identity') +
    geom_bar(stat = 'identity', width = 1, colour = "white") +
    coord_polar("y", start = 0) +  # Start at -pi/2 radians (bottom of the pie)
    # labs(title = NULL) +
    theme_void() +  # Optional: to remove axis lines and labels
    # geom_text(aes(label = angle_name)) +  # Add labels
    
    scale_fill_manual(values = custom_colors)+
    theme(legend.position = "none")
}



BC_tukey <- function(phylo_obj) {
  ## Identify the correspondence: group and samples
  group2samp <- list() # list to save the correspondence between group <--> samples
  group_list <- as.factor(get_variable(sample_data(phylo_obj), 'Group')) # list of group elements
  for (groups in levels(group_list)) { # loop over the no. of group levels
    target_group <- which(group_list == groups) # vct pos of the curr group variable 
    group2samp[[ groups ]] <- sample_names(phylo_obj)[target_group] # matching samples: based on vct pos
  }  
  
  ## Calculate beta-diversity
  beta_div_dist <- phyloseq::distance(physeq = phylo_obj, method = 'bray')
  beta_div_dist <- as(beta_div_dist, "matrix")
  
  ## Coerce distance mtx into a tidy data frame by group
  dist_df <- data.frame() # save results in df 
  counter <- 1 
  for (groups in names(group2samp)) { # loop over group fct levels 
    sub_dist <- beta_div_dist[ group2samp[[groups]], group2samp[[groups]] ] # subset dist mtx to curr samples
    #print(sub_dist)
    no_samp_col <- ncol(sub_dist) # n cols: curr sub dist
    no_samp_row <- nrow(sub_dist) # n rows: curr sub dist
    for ( cols in seq(no_samp_col) ) { # loop over cols: curr sub_dist
      if ( cols > 1 ) {
        for ( rows in seq((cols-1)) ) { # loop over rows: curr sub_dist 
          ## Save results
          dist_df[ counter, "sample_pair" ] <- paste0( colnames(sub_dist)[cols], "-",  
                                                       rownames(sub_dist)[rows] ) # sample pair
          dist_df[ counter, "group" ] <- groups # group  
          dist_df[ counter, "beta_div_method" ] <- 'bray' # method
          dist_df[ counter, "beta_div_value" ] <- sub_dist[rows, cols] # beta-diversity for the sample pair     
          counter = counter + 1
        }
      }
    }
  }
  dist_df$group <- factor(dist_df$group, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))
  plot <- ggplot(dist_df, aes(x = group, y = beta_div_value, color = group))+
    geom_boxplot()+
    geom_jitter(color = 'gray', width = 0.1)
  
  lm <- lm(beta_div_value ~ group, dist_df)
  lm.av <- aov(lm)
  tukey.result <- HSD.test(lm.av, trt = 'group')
  return(list(tukey.result, plot))
}


BC_t.test <- function(pl_ob_subset, group1, group2){
  dist <- phyloseq::distance(pl_ob_subset, method = 'bray')
  dist <- as.matrix(dist)
  
  # assign '0' to upper triangle of the dist matrix
  dist[upper.tri(dist)] <- 0
  
  # using melt() in 'reshape2' package to convert dist matrix to dataframe in common form and delete previously assigned '0'.
  library(reshape2)
  dist_melt <- subset(melt(dist), value!=0)
  
  # assign Sample_ID from metadata file to dist dataframe - dist_melt
  # Var1 and Var2 are column names of dist_melt, Sample_ID is a column name in metadata (load metadata first if have not)
  metadata.in.fun <- read_xlsx('data/metadata.xlsx')
  library(dplyr)
  temp <- dplyr::left_join(dist_melt, metadata.in.fun, by = c ('Var1'='Sample_ID'))
  temp2 <- dplyr::left_join(temp, metadata.in.fun, by = c ('Var2'='Sample_ID'))
  
  # generate new columns, and exclude unwanted columns
  temp3 <- temp2 %>% dplyr::transmute(distance = value, Compartment = Compartment.x, Pair = paste0(Mycorrhizal.x, sep = "_VS_", Mycorrhizal.y, sep='_', Cd.x, sep='_', Cd.y))
  
  cat('group to select:', unique(temp3$Pair))
  
  
  temp4 <- temp3 %>% filter(Pair == group1 | Pair == group2)
  
  
  plot <- temp4 %>% ggplot(aes(x=Pair, y = distance, color = Pair))+
    geom_boxplot()+
    geom_jitter(alpha = 0.5, width = 0.2)
  
  t_test_result <- t.test(distance ~ Pair, data = temp4)
  
  return(list(plot, t_test_result))
  
}


cor_coeff_perc <- function(rs){
  
  c1 <- sum(rs$cor_coeff >= -1 & rs$cor_coeff <= -0.7)
  c2 <- sum(rs$cor_coeff > -0.7 & rs$cor_coeff <= -0.3)
  c3 <- sum(rs$cor_coeff > -0.3 & rs$cor_coeff <= 0)
  
  c4 <- sum(rs$cor_coeff > 0 & rs$cor_coeff <=0.3)
  c5 <- sum(rs$cor_coeff > 0.3 & rs$cor_coeff <=0.7)
  c6 <- sum(rs$cor_coeff > 0.7 & rs$cor_coeff <=1)
  
  temp_check <- sum(unlist(mget(paste0("c", 1:6)))) == length(rs$cor_coeff)
  
  if (temp_check) {
    print("Count sum = total number of slope angles. Result in %:")
  }
  
  if (!temp_check) {
    stop("Count sum and number of slope angles not the same. Pls. check data.")
  }
  
  rho_group <- c(c1, c2, c3, c4, c5, c6)
  cor.coeff.group <- c('-1.0 ~ -0.7', '-0.7 ~ -0.3', '-0.3 ~ 0.0', 
                       '0.0 ~ 0.3', '0.3 ~ 0.7', '0.7 ~ 1.0')
  perc.rs <- numeric(length(rho_group)) 
  
  for(i in seq_along(rho_group)) {
    
    
    perc.rs[i] <- rho_group[i]*100/length(rs$cor_coeff)
  }
  
  
  a.perc.tab <- as.data.frame(rbind(cor.coeff.group, perc.rs))
  a.perc.tab.t <- t(a.perc.tab)
  a.perc.tab.t.df <- as.data.frame(a.perc.tab.t)
  
  new.colname <- deparse(substitute(rs))
  names(a.perc.tab.t.df)[names(a.perc.tab.t.df) == "perc.rs"] <- new.colname
  
  # colnames(a.perc.tab) <- a.perc.tab[1, ]
  # a.perc.tab <- a.perc.tab[-1, ]
  
  return(a.perc.tab.t.df)
  
}

cc.compare <- function(df, subset){
  df %>% filter(grepl(subset, group)) %>% 
    ggplot(aes(x = cor.coeff.group, y = as.numeric(perc.), fill = group)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = NULL,
         x = "Correlation Coefficient Group",
         y = "Percentage (%)",
         fill = "Group")
}


cal.even.inner <- function(data.vector){
  H <- vegan::diversity(data.vector, index = "shannon")
  S <- length(data.vector)
  J <- H / log(S)
  return(J)
  
}


group.angle <- function(angle.vector, break_by = 30) {
  breaks <- seq(-90, 90, by = break_by)
  grouped <- cut(angle.vector, breaks = breaks, include.lowest = TRUE)
  counts <- table(grouped)
  return(counts)
}






calculate_evenness <- function(data) {
  return(cal.even.inner(group.slope(data$slope_a)))
}

cal.chao1 <- function(data.vector){
  counts <- group.slope(data.vector$slope_a)
  # Convert the vector to a matrix with one row
  mt <- matrix(counts, nrow = 1)
  
  # Calculate Chao1 index
  rs <- estimateR(mt)
  
  # Extract Chao1 value
  chao1 <- rs["S.chao1", ]
  return(chao1)
}


group.slope <- function(vector, num_groups = 999){
  min <- min(vector)
  max <- max(vector)
  breaks <- seq(min, max, length.out = num_groups+1)
  grouped <- cut(vector, breaks = breaks, include.lowest = TRUE)
  counts <- table(grouped)
  return(counts)
}


cal.ace <- function(data.vector){
  counts <- group.slope(data.vector$slope_a)
  # Convert the vector to a matrix with one row
  mt <- matrix(counts, nrow = 1)
  
  # Calculate Chao1 index
  rs <- estimateR(mt)
  
  # Extract Chao1 value
  ACE <- rs["S.ACE", ]
  return(ACE)
}

boot_chao1 <- function(data, indices) {
  d <- data[indices]
  return(estimateR(d)["S.chao1"])
}

boot_even <- function(data, indices) {
  d <- data[indices]
  return(cal.even.inner(d))
}

boot_ace <- function(data, indices) {
  d <- data[indices]
  return(estimateR(d)["S.ACE"])
}


# cal_rob <- function(graph, type = "vertex", measure = "degree", N = 1000) {
#   robustness <- numeric(N)
#   for (i in 1:N) {
#     if (type == "vertex") {
#       if (measure == "degree") {
#         node_to_remove <- which.max(degree(graph))
#       } else if (measure == "btwn.cent") {
#         node_to_remove <- which.max(betweenness(graph))
#       } else {
#         node_to_remove <- sample(V(graph), 1)
#       }
#       graph <- delete_vertices(graph, node_to_remove)
#     } else {
#       if (measure == "random") {
#         edge_to_remove <- sample(E(graph), 1)
#       } else {
#         edge_to_remove <- which.max(edge_betweenness(graph))
#       }
#       graph <- delete_edges(graph, edge_to_remove)
#     }
#     robustness[i] <- mean(degree(graph))
#   }
#   return(mean(robustness))
# }


# Example usage:
# result <- node_rm_simu(pl_un_rarefied, 'Rhizosphere', 'Cd', c_cutoff = 5, rho_cutoff = 0.8)
# Removing node simulation ----
# Robustness simulation functions
rand.remov.once <- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) {
  id.rm <- sample(1:nrow(netRaw), round(nrow(netRaw) * rm.percent))
  net.Raw <- netRaw
  net.Raw[id.rm, ] <- 0
  net.Raw[, id.rm] <- 0
  if (abundance.weighted) {
    net.stength <- net.Raw * sp.ra
  } else {
    net.stength <- net.Raw
  }
  sp.meanInteration <- colMeans(net.stength)
  id.rm2 <- which(sp.meanInteration <= 0)
  remain.percent <- (nrow(netRaw) - length(id.rm2)) / nrow(netRaw)
  remain.percent
  }
  
rmsimu <- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm = 100) {
  t(sapply(rm.p.list, function(x) {
  remains <- sapply(1:nperm, function(i) {
  rand.remov.once(netRaw = netRaw, rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted)
  })
  remain.mean <- mean(remains)
  remain.sd <- sd(remains)
  remain.se <- sd(remains) / sqrt(nperm)
  result <- c(remain.mean, remain.sd, remain.se)
  names(result) <- c("remain.mean", "remain.sd", "remain.se")
  result
  }))
  }
  

  
ask_to_save <- function() {
  response <- readline(prompt = "Save the files such as plots and tables? (type 'yes' or 'no'): ")
  if (tolower(response) == "yes") {
    save <- TRUE
  } else {
    save <- FALSE
  }
  return(save)
}



