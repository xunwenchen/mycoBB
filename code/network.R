# prepare for network analysis
# prepare OTU table to be used in MENA
# http://ieg4.rccc.ou.edu/MENA/

# Endosphere compartment OTU tables ----

# pl_endo_cd <- otu_table(pl %>% subset_samples(Compartment =='Endosphere' & Group == 'Cd'))
# pl_endo_cd <- as.data.frame(pl_endo_cd)
# library(tibble)
# pl_endo_cd <- rownames_to_column(pl_endo_cd, "ID")
# pl_endo_cd[pl_endo_cd == 0] <- ''  
# write.table(pl_endo_cd, file = "out/pl_endo_cd.txt", sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)
# rm(pl_endo_cd)


# pl_endo_ctrl <- otu_table(pl %>% subset_samples(Compartment =='Endosphere' & Group == 'Ctrl'))
# pl_endo_ctrl[pl_endo_ctrl == 0] <- ''   
# write.table(pl_endo_ctrl, file = "out/pl_endo_ctrl.txt", sep = "\t",
#             row.names = TRUE, col.names = TRUE, quote = FALSE)
# rm(pl_endo_ctrl)
# 
# 
# pl_endo_m <- otu_table(pl %>% subset_samples(Compartment =='Endosphere' & Group == 'M'))
# pl_endo_m[pl_endo_m == 0] <- ''   
# write.table(pl_endo_m, file = "out/pl_endo_m.txt", sep = "\t",
#             row.names = TRUE, col.names = TRUE, quote = FALSE)
# rm(pl_endo_m)
# 
# 
# pl_endo_m_cd <- otu_table(pl %>% subset_samples(Compartment =='Endosphere' & Group == 'M+Cd'))
# pl_endo_m_cd[pl_endo_m_cd == 0] <- ''   
# write.table(pl_endo_m_cd, file = "out/pl_endo_m_cd.txt", sep = "\t",
#             row.names = TRUE, col.names = TRUE, quote = FALSE)
# rm(pl_endo_m_cd)


# Rhizoplane compartment OTU tables ----
# pl_rzp_cd <- otu_table(pl %>% subset_samples(Compartment =='Rhizoplane' & Group == 'Cd'))
# pl_rzp_cd <- as.data.frame(pl_rzp_cd)
# library(tibble)
# pl_rzp_cd <- rownames_to_column(pl_rzp_cd, "ID")
# pl_rzp_cd[pl_rzp_cd == 0] <- ''   
# write.table(pl_rzp_cd, file = "out/pl_rzp_cd.txt", sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)
# rm(pl_rzp_cd)
# 
# 
# 
# 
# pl_rzp_ctrl <- otu_table(pl %>% subset_samples(Compartment =='Rhizoplane' & Group == 'Ctrl'))
# pl_rzp_ctrl <- as.data.frame(pl_rzp_ctrl)
# library(tibble)
# pl_rzp_ctrl <- rownames_to_column(pl_rzp_ctrl, "ID")
# pl_rzp_ctrl[pl_rzp_ctrl == 0] <- ''   
# write.table(pl_rzp_ctrl, file = "out/pl_rzp_ctrl.txt", sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)
# rm(pl_rzp_ctrl)
# 
# 
# 
# pl_rzp_m <- otu_table(pl %>% subset_samples(Compartment =='Rhizoplane' & Group == 'M'))
# pl_rzp_m <- as.data.frame(pl_rzp_m)
# library(tibble)
# pl_rzp_m <- rownames_to_column(pl_rzp_m, "ID")
# pl_rzp_m[pl_rzp_m == 0] <- ''   
# write.table(pl_rzp_m, file = "out/pl_rzp_m.txt", sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)
# rm(pl_rzp_m)
# 
# 
# 
# pl_rzp_m_cd <- otu_table(pl %>% subset_samples(Compartment =='Rhizoplane' & Group == 'M+Cd'))
# pl_rzp_m_cd <- as.data.frame(pl_rzp_m_cd)
# library(tibble)
# pl_rzp_m_cd <- rownames_to_column(pl_rzp_m_cd, "ID")
# pl_rzp_m_cd[pl_rzp_m_cd == 0] <- ''   
# write.table(pl_rzp_m_cd, file = "out/pl_rzp_m_cd.txt", sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)
# rm(pl_rzp_m_cd)


# Rhizosphere compartment OTU tables ----
pl_rzs_cd <- otu_table(pl %>% subset_samples(Compartment =='Rhizosphere' & Group == 'Cd'))
pl_rzs_cd <- as.data.frame(pl_rzs_cd)
library(tibble)
pl_rzs_cd <- rownames_to_column(pl_rzs_cd, "ID")
pl_rzs_cd[pl_rzs_cd == 0] <- ''
write.table(pl_rzs_cd, file = "data/data_for_network_construction/pl_rzs_cd.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
rm(pl_rzs_cd)


pl_rzs_ctrl <- otu_table(pl %>% subset_samples(Compartment =='Rhizosphere' & Group == 'Ctrl'))
pl_rzs_ctrl <- as.data.frame(pl_rzs_ctrl)
pl_rzs_ctrl <- rownames_to_column(pl_rzs_ctrl, "ID")
pl_rzs_ctrl[pl_rzs_ctrl == 0] <- ''
write.table(pl_rzs_ctrl, file = "data/data_for_network_construction/pl_rzs_ctrl.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
rm(pl_rzs_ctrl)


pl_rzs_m <- otu_table(pl %>% subset_samples(Compartment =='Rhizosphere' & Group == 'M'))
pl_rzs_m <- as.data.frame(pl_rzs_m)
library(tibble)
pl_rzs_m <- rownames_to_column(pl_rzs_m, "ID")
pl_rzs_m[pl_rzs_m == 0] <- ''
write.table(pl_rzs_m, file = "data/data_for_network_construction/pl_rzs_m.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
rm(pl_rzs_m)


pl_rzs_m_cd <- otu_table(pl %>% subset_samples(Compartment =='Rhizosphere' & Group == 'M+Cd'))
pl_rzs_m_cd <- as.data.frame(pl_rzs_m_cd)
library(tibble)
pl_rzs_m_cd <- rownames_to_column(pl_rzs_m_cd, "ID")
pl_rzs_m_cd[pl_rzs_m_cd == 0] <- ''
write.table(pl_rzs_m_cd, file = "data/data_for_network_construction/pl_rzs_m_cd.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
rm(pl_rzs_m_cd)



# 
# 
# 
# 
# 
# 
# Hyphosphere compartment OTU tables ----
pl_hps_cd <- otu_table(pl %>% subset_samples(Compartment =='Hyphosphere' & Group == 'Cd'))
pl_hps_cd <- as.data.frame(pl_hps_cd)
pl_hps_cd <- rownames_to_column(pl_hps_cd, "ID")
pl_hps_cd[pl_hps_cd == 0] <- ''
write.table(pl_hps_cd, file = "data/data_for_network_construction/pl_hps_cd.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
rm(pl_hps_cd)


pl_hps_ctrl <- otu_table(pl %>% subset_samples(Compartment =='Hyphosphere' & Group == 'Ctrl'))
pl_hps_ctrl <- as.data.frame(pl_hps_ctrl)
pl_hps_ctrl <- rownames_to_column(pl_hps_ctrl, "ID")
pl_hps_ctrl[pl_hps_ctrl == 0] <- ''
write.table(pl_hps_ctrl, file = "data/data_for_network_construction/pl_hps_ctrl.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
rm(pl_hps_ctrl)


# ~~ the following two may have problem ----
pl_hps_m <- otu_table(pl %>% subset_samples(Compartment =='Hyphosphere' & Group == 'M'))
pl_hps_m <- as.data.frame(pl_hps_m)
pl_hps_m <- rownames_to_column(pl_hps_m, "ID")
pl_hps_m[pl_hps_m == 0] <- ""
write.table(pl_hps_m, file = "data/data_for_network_construction/pl_hps_m.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
rm(pl_hps_m)


pl_hps_m_cd <- otu_table(pl %>% subset_samples(Compartment =='Hyphosphere' & Group == 'M+Cd'))
pl_hps_m_cd <- as.data.frame(pl_hps_m_cd)
pl_hps_m_cd <- rownames_to_column(pl_hps_m_cd, "ID")
pl_hps_m_cd[pl_hps_m_cd == 0] <- ""
write.table(pl_hps_m_cd, file = "data/data_for_network_construction/pl_hps_m_cd.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
rm(pl_hps_m_cd)

# 
# 
# 
# 
# Hyphae compartment OTU tables ----
# no pl_hph_cd and pl_hph_ctrl since Cd and Ctrl treatments did not have mycorrhizal inoculation 

pl_hph_m <- otu_table(pl %>% subset_samples(Compartment =='Hyphae' & Group == 'M'))
pl_hph_m <- as.data.frame(pl_hph_m)
pl_hph_m <- rownames_to_column(pl_hph_m, "ID")
pl_hph_m[pl_hph_m == 0] <- ''
write.table(pl_hph_m, file = "data/data_for_network_construction/pl_hph_m.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
rm(pl_hph_m)

pl_hph_m_cd <- otu_table(pl %>% subset_samples(Compartment =='Hyphae' & Group == 'M+Cd'))
pl_hph_m_cd <- as.data.frame(pl_hph_m_cd)
pl_hph_m_cd <- rownames_to_column(pl_hph_m_cd, "ID")
pl_hph_m_cd[pl_hph_m_cd == 0] <- ''
write.table(pl_hph_m_cd, file = "data/data_for_network_construction/pl_hph_m_cd.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
rm(pl_hph_m_cd)


# 
# 
# Bulk soil compartment OTU tables ----
pl_bs_cd <- otu_table(pl %>% subset_samples(Compartment =='Bulk soil' & Group == 'Cd'))
pl_bs_cd <- as.data.frame(pl_bs_cd)
pl_bs_cd <- rownames_to_column(pl_bs_cd, "ID")
pl_bs_cd[pl_bs_cd == 0] <- ''   
write.table(pl_bs_cd, file = "data/data_for_network_construction/pl_bs_cd.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
rm(pl_bs_cd)


pl_bs_ctrl <- otu_table(pl %>% subset_samples(Compartment =='Bulk soil' & Group == 'Ctrl'))
pl_bs_ctrl <- as.data.frame(pl_bs_ctrl)
pl_bs_ctrl <- rownames_to_column(pl_bs_ctrl, "ID")
pl_bs_ctrl[pl_bs_ctrl == 0] <- ''   
write.table(pl_bs_ctrl, file = "data/data_for_network_construction/pl_bs_ctrl.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
rm(pl_bs_ctrl)


pl_bs_m <- otu_table(pl %>% subset_samples(Compartment =='Bulk soil' & Group == 'M'))
pl_bs_m <- as.data.frame(pl_bs_m)
pl_bs_m <- rownames_to_column(pl_bs_m, "ID")
pl_bs_m[pl_bs_m == 0] <- ''   
write.table(pl_bs_m, file = "data/data_for_network_construction/pl_bs_m.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
rm(pl_bs_m)


pl_bs_m_cd <- otu_table(pl %>% subset_samples(Compartment =='Bulk soil' & Group == 'M+Cd'))
pl_bs_m_cd <- as.data.frame(pl_bs_m_cd)
pl_bs_m_cd <- rownames_to_column(pl_bs_m_cd, "ID")
pl_bs_m_cd[pl_bs_m_cd == 0] <- ''   
write.table(pl_bs_m_cd, file = "data/data_for_network_construction/pl_bs_m_cd.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
rm(pl_bs_m_cd)


# the above txt files were used to obtain node files using MENA (http://ieg4.rccc.ou.edu/MENA/). The obtained files are stored in subfolders of 'data/network_data', depending on groups. Then, the node files were assigned with taxonomic info using the following codes for each group. # sif files were used in Cytoscape. We used other files.
# 
# construct zotu-tax table ----
zotu_tax_tab <- as.data.frame(tax_table(pl))
zotu_tax_tab <- rownames_to_column(zotu_tax_tab, "ID")
write.table(zotu_tax_tab, file = "data/data_for_network_construction/zotu_tax_tab.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
# change the column name "ID" to "Name"
# Change the column name "ID" to "Name"
colnames(zotu_tax_tab)[colnames(zotu_tax_tab) == "ID"] <- "Name"



# bulk soil compartment
# Ctrl group
write.csv(node_merge_otu('data/network_data/bs_majority_7/bs_ctrl_node.csv'),
          'data/network_data/bs_majority_7/bs_ctrl_node_final.csv',
          row.names = FALSE)

# Cd group
write.csv(node_merge_otu('data/network_data/bs_majority_7/bs_cd_node.csv'),
          'data/network_data/bs_majority_7/bs_cd_node_final.csv',
          row.names = FALSE)

# M group
write.csv(node_merge_otu('data/network_data/bs_majority_7/bs_m_node.csv'),
          'data/network_data/bs_majority_7/bs_m_node_final.csv',
          row.names = FALSE)

# M+Cd group
write.csv(node_merge_otu('data/network_data/bs_majority_7/bs_m_cd_node.csv'),
          'data/network_data/bs_majority_7/bs_m_cd_node_final.csv',
          row.names = FALSE)


# Hyphosphere -- ctrl treatment (majority = 6, tweak in MENA pipeline) ----
write.csv(node_merge_otu('data/network_data/hps_majority_6/hps_ctrl_node.csv'),
          'data/network_data/hps_majority_6/hps_ctrl_node_final.csv',
          row.names = FALSE)

# Hyphosphere -- Cd treatment (majority = 6) ----
write.csv(node_merge_otu('data/network_data/hps_majority_6/hps_cd_node.csv'),
          'data/network_data/hps_majority_6/hps_cd_node_final.csv',
          row.names = FALSE)

# Hyphosphere -- M treatment (majority = 6) ----
write.csv(node_merge_otu('data/network_data/hps_majority_6/hps_m_node.csv'),
          'data/network_data/hps_majority_6/hps_m_node_final.csv',
          row.names = FALSE)

# Hyphosphere -- M + Cd treatment (majority = 6) ----
write.csv(node_merge_otu('data/network_data/hps_majority_6/hps_m_cd_node.csv'),
          'data/network_data/hps_majority_6/hps_m_cd_node_final.csv',
          row.names = FALSE)

# Hyphosphere -- Ctrl treatment (majority = 7, more stringent than 6)----
write.csv(node_merge_otu('data/network_data/hps_majority_7/hps_ctrl_node.csv'),
          'data/network_data/hps_majority_7/hps_ctrl_node_final.csv',
          row.names = FALSE)

# Hyphosphere -- Cd treatment (majority = 7, more stringent than 6)----
write.csv(node_merge_otu('data/network_data/hps_majority_7/hps_cd_node.csv'),
          'data/network_data/hps_majority_7/hps_cd_node_final.csv',
          row.names = FALSE)

# Hyphosphere -- Ctrl treatment (majority = 8, more stringent than 6 and 7)----
write.csv(node_merge_otu('data/network_data/hps_majority_8/hps_ctrl_node.csv'),
          'data/network_data/hps_majority_8/hps_ctrl_node_final.csv',
          row.names = FALSE)

# Hyphosphere -- Cd treatment (majority = 8, more stringent than 6 and 7)----
write.csv(node_merge_otu('data/network_data/hps_majority_8/hps_cd_node.csv'),
          'data/network_data/hps_majority_8/hps_cd_node_final.csv',
          row.names = FALSE)

# Hyphosphere -- M treatment (majority = 5, less stringent than 6)----
write.csv(node_merge_otu('data/network_data/hps_majority_5/hps_m_node.csv'),
          'data/network_data/hps_majority_5/hps_m_node_final.csv',
          row.names = FALSE)
# Hyphosphere -- M+Cd treatment (majority = 5, less stringent than 6)----
write.csv(node_merge_otu('data/network_data/hps_majority_5/hps_m_cd_node.csv'),
          'data/network_data/hps_majority_5/hps_m_cd_node_final.csv',
          row.names = FALSE)


# Rhizosphere compartment

write.csv(node_merge_otu('data/network_data/rzs_majority_7/rzs_ctrl_node.csv'),
          'data/network_data/rzs_majority_7/rzs_ctrl_node_final.csv',
          row.names = FALSE)

write.csv(node_merge_otu('data/network_data/rzs_majority_7/rzs_cd_node.csv'),
          'data/network_data/rzs_majority_7/rzs_cd_node_final.csv',
          row.names = FALSE)

write.csv(node_merge_otu('data/network_data/rzs_majority_7/rzs_m_node.csv'),
          'data/network_data/rzs_majority_7/rzs_m_node_final.csv',
          row.names = FALSE)

write.csv(node_merge_otu('data/network_data/rzs_majority_7/rzs_m_cd_node.csv'),
          'data/network_data/rzs_majority_7/rzs_m_cd_node_final.csv',
          row.names = FALSE)



# hyphae-associated -- M treatment ----
# 
# construct_node_file <- function(address_of_raw_file){
#   df <- read.csv(address_of_raw_file)
#   
#   zotu_tax_tab <- as.data.frame(TAX)
#   zotu_tax_tab <- rownames_to_column(zotu_tax_tab, "Name") # here used 'Name' is for later merging with node file hps_m_node
#   df2 <- merge(df, zotu_tax_tab, by = 'Name', all.x = TRUE)
#   
# }

hph_m_node <- read.csv('data/network_data/hph_m_node_raw.csv')
hph_m_node_final <- merge(hph_m_node, zotu_tax_tab, by = 'Name', all.x = TRUE)


write.csv(hph_m_node_final,"data/network_data/hph_m_node_final.csv", row.names = FALSE)



# hyphae-associated -- M+Cd treatment ----
hph_m_cd_node <- read.csv('data/network_data/hph_m_cd_node_raw.csv')
hph_m_cd_node_final <- merge(hph_m_cd_node, zotu_tax_tab, by = 'Name', all.x = TRUE)


write.csv(hph_m_cd_node_final,"data/network_data/hph_m_cd_node_final.csv", row.names = FALSE)


# Rhizosphere - ctrl treatment ----
rzs_ctrl_node <- read.csv('data/network_data/rzs_ctrl_node_raw.csv')
rzs_ctrl_node_final <- merge(rzs_ctrl_node, zotu_tax_tab, by = 'Name', all.x = TRUE)

write.csv(rzs_ctrl_node_final,"data/network_data/rzs_ctrl_node_final.csv", row.names = FALSE)

# Rhizosphere - ctrl treatment 2 ----
rzs_ctrl_node2 <- read.csv('data/network_data/rzs_ctrl_node_raw2.csv')
rzs_ctrl_node_final2 <- merge(rzs_ctrl_node2, zotu_tax_tab, by = 'Name', all.x = TRUE)


write.csv(rzs_ctrl_node_final2,"data/network_data/rzs_ctrl_node_final2.csv", row.names = FALSE)


# SUMMARizing network centrality measures ----
# ~~ Zi-Pi plots -- 10 plots ----
# load data
df <- read.csv('data/network_data/rzs_majority_7/rzs_ctrl_node_final.csv')
p1 <- plot_pi_zi(df)+ggtitle('Rhizosphere - Ctrl')

df <- read.csv('data/network_data/rzs_majority_7/rzs_cd_node_final.csv')
p2 <- plot_pi_zi(df)+ggtitle('Rhizosphere - Cd')

df <- read.csv('data/network_data/rzs_majority_7/rzs_m_node_final.csv')
p3 <- plot_pi_zi(df)+ggtitle('Rhizosphere - M')

df <- read.csv('data/network_data/rzs_majority_7/rzs_m_cd_node_final.csv')
p4 <- plot_pi_zi(df)+ggtitle('Rhizosphere - M+Cd')

df <- read.csv('data/network_data/hps_majority_6/hps_ctrl_node_final.csv')
p5 <- plot_pi_zi(df)+ggtitle('Hyphosphere - Ctrl')

df <- read.csv('data/network_data/hps_majority_6/hps_cd_node_final.csv')
p6 <- plot_pi_zi(df)+ggtitle('Hyphosphere - Cd')

df <- read.csv('data/network_data/hps_majority_6/hps_m_node_final.csv')
p7 <- plot_pi_zi(df)+ggtitle('Hyphosphere - M')

df <- read.csv('data/network_data/hps_majority_6/hps_m_cd_node_final.csv')
p8 <- plot_pi_zi(df)+ggtitle('Hyphosphere - M+Cd')


df <- read.csv('data/network_data/hph_m_node_final.csv')
p9 <- plot_pi_zi(df)+ggtitle('Hyphae - M')

df <- read.csv('data/network_data/hph_m_cd_node_final.csv')
p10 <- plot_pi_zi(df)+ggtitle('Hyphae - M+Cd')

# ~~~~ Additional Figure. zi-pi plot ----
pi_zi_plot <- ggarrange(p1, p2, p3, p4,
          p5, p6, p7, p8,
          p9, p10,
          labels = c('a', 'b', 'c', 'd', 
                     'e', 'f', 'g', 'h', 
                     'i', 'j'),
          nrow = 3, ncol = 4)

ggsave('out/pi_zi.pdf', plot = pi_zi_plot,
       w = 88*3, h = 66*3.5, units = 'mm')
ggsave('out/pi_zi.jpg', plot = pi_zi_plot,
       w = 88*3, h = 66*3.5, units = 'mm', dpi = 300)
rm(p1, p2, p3, p4,
   p5, p6, p7, p8,
   p9, p10)


# node degree and betweenness corrected plot ----

# ~ load all data ----
# rhizosphere
df1 <- read.csv('data/network_data/rzs_majority_7/rzs_ctrl_node_final.csv')
df2 <- read.csv('data/network_data/rzs_majority_7/rzs_cd_node_final.csv')
df3<- read.csv('data/network_data/rzs_majority_7/rzs_m_node_final.csv')
df4 <- read.csv('data/network_data/rzs_majority_7/rzs_m_cd_node_final.csv')

# hyphosphere
df5 <- read.csv('data/network_data/hps_majority_6/hps_ctrl_node_final.csv')
df6 <- read.csv('data/network_data/hps_majority_6/hps_cd_node_final.csv')
df7 <- read.csv('data/network_data/hps_majority_6/hps_m_node_final.csv')
df8 <- read.csv('data/network_data/hps_majority_6/hps_m_cd_node_final.csv')

# hyphae-associated
df9 <- read.csv('data/network_data/hph_m_node_final.csv')
df10 <- read.csv('data/network_data/hph_m_cd_node_final.csv')



# ~~ Node degree ----
# ~~~ construct functions to combine tables ----
combine_node_deg_tables <- function(df1, df2, df3, df4,tax_level, tax_for_merge){
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'Ctrl') 
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'Cd') 
  c <- df3 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M')
  d <- df4 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M+Cd')
  
  df_list <- list(a, b, c, d)
  df_merge <- Reduce(function(x, y) merge(x, y, by = tax_for_merge, all=TRUE), df_list)
  return(df_merge)
}
combine_node_deg_tables_hph <- function(df1, df2, tax_level, tax_for_merge){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M')
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M+Cd')
  
  df_list <- list(a, b)
  df_merge <- Reduce(function(x, y) merge(x, y, by = tax_for_merge, all=TRUE), df_list)
  return(df_merge)
}

# combine tables. Family level has been finished so convert to comment using #
# rhizosphere
rzs_node_deg_phylum <- combine_node_deg_tables(df1, df2, df3, df4, Phylum, 'Phylum')
write.csv(rzs_node_deg_phylum, 'out/rzs_node_deg_phylum.csv')
# rzs_node_deg_Family <- combine_node_deg_tables(df1, df2, df3, df4, Family, 'Family')

# hyphosphere
hps_node_deg_phylum <- combine_node_deg_tables(df5, df6, df7, df8, Phylum, 'Phylum')
write.csv(hps_node_deg_phylum, 'out/hps_node_deg_phylum.csv')
# hps_node_deg_Family <- combine_node_deg_tables(df5, df6, df7, df8, Family, 'Family')

# hyphae-associated
hph_node_deg_phylum <- combine_node_deg_tables_hph(df9, df10, Phylum, 'Phylum')
write.csv(hph_node_deg_phylum, 'out/hph_node_deg_phylum.csv')
# hph_node_deg_Family <- combine_node_deg_tables_hph(df9, df10, Family, 'Family')



# Node degree corrected plots ----
# the data was obtained using R then modified in excel, then reloaded here
# load data
# rzs
rzs_family_node_degree_mod <- read.csv('data/network_data/Table S2.1. rzs_family_node_degree_mod.csv')
rzs_family_node_degree_mod$Treatment <- factor(rzs_family_node_degree_mod$Treatment, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))


rzs_phylum_node_degree_mod <- read.csv('data/network_data/rzs_node_deg_phylum_mod.csv')
rzs_phylum_node_degree_mod$Treatment <- factor(rzs_phylum_node_degree_mod$Treatment, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))


# hps
hps_family_node_degree_mod <- read.csv('data/network_data/Table S2.2. hps_family_node_degree_mod.csv')
hps_family_node_degree_mod$Treatment <- factor(hps_family_node_degree_mod$Treatment, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))


hps_phylum_node_degree_mod <- read.csv('data/network_data/hps_node_deg_phylum_mod.csv')
hps_phylum_node_degree_mod$Treatment <- factor(hps_phylum_node_degree_mod$Treatment, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))



# hph
hph_family_node_degree_mod <- read.csv('data/network_data/Table S2.3. hph_family_node_degree_mod.csv')
hph_family_node_degree_mod$Treatment <- factor(hph_family_node_degree_mod$Treatment, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))


hph_phylum_node_degree_mod <- read.csv('data/network_data/hph_node_deg_phylum_mod.csv')
hph_phylum_node_degree_mod$Treatment <- factor(hph_phylum_node_degree_mod$Treatment, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))




p1 <- ggplot(rzs_family_node_degree_mod, aes(Family, Accum_node_deg, color = Treatment, group = Treatment))+
  #geom_point(size = 4, alpha = 0.4)+
  geom_line()+
  ggtitle('Rhizosphere')+
  theme(legend.position = c(0.1, 0.6))+
  ylab('Accumulative node degree')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        legend.title = element_blank())


p2 <- ggplot(hps_family_node_degree_mod, aes(Family, Accum_node_deg, color = Treatment, group = Treatment))+
  #geom_point(size = 4, alpha = 0.4)+
  geom_line()+
  ggtitle('Hyphosphere')+
  theme(legend.position = c(0.1, 0.6))+
  ylab('Accumulative node degree')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        legend.title = element_blank())

p3 <- ggplot(hph_family_node_degree_mod, aes(Family, Accum_node_deg, color = Treatment, group = Treatment))+
  #geom_point(size = 4, alpha = 0.4)+
  geom_line()+
  ggtitle('Hyphae')+
  theme(legend.position = c(0.1, 0.8))+
  ylab('Accumulative node degree')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        legend.title = element_blank())

# ~~~~ Fig. S13. Accumulative node degree (family level) ----
ggarrange(p1, p2, p3,
          labels = c('a', 'b', 'c'),
          nrow = 3, ncol = 1)
if(save){
ggsave('out/Fig. S13. node_plot_family2.pdf',
       w = 88*2, h = 66*3, units = 'mm')

ggsave('out/Fig. S13. node_plot_family2.jpg',
       w = 88*2, h = 66*3, units = 'mm')
}
rm(p1, p2, p3)


t1 <- ggplot(rzs_phylum_node_degree_mod, aes(Phylum, Accum_degree, color = Treatment, group = Treatment))+
  geom_point(size = 4, alpha = 0.4)+
  geom_line()+
  ggtitle('Rhizosphere')+
  theme(legend.position = c(0.8, 0.5))+
  ylab('Accumulative node degree')+
  scale_x_discrete(limits=rev)+
  coord_flip()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank())


t2 <- ggplot(hps_phylum_node_degree_mod, aes(Phylum, Accum_degree, color = Treatment, group = Treatment))+
  geom_point(size = 4, alpha = 0.4)+
  geom_line()+
  ggtitle('Hyphosphere')+
  theme(legend.position = c(0.8, 0.5))+
  ylab('Accumulative node degree')+
  scale_x_discrete(limits=rev)+
  coord_flip()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank())

t3 <- ggplot(hph_phylum_node_degree_mod, aes(Phylum, Accum_degree, color = Treatment, group = Treatment))+
  geom_point(size = 4, alpha = 0.4)+
  geom_line()+
  ggtitle('Hyphae')+
  theme(legend.position = c(0.8, 0.5))+
  ylab('Accumulative node degree')+
  scale_x_discrete(limits=rev)+
  coord_flip()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank())



# ~~~~ Fig. S12. Accumulative node degree (phylum level) ----
ggarrange(t1, t2, t3,
          labels = c('a', 'b', 'c'),
          nrow = 1, ncol = 3)

ggsave('out/Fig. S12. node_plot_phylum2.pdf',
       w = 88*4, h = 88*1.5, units = 'mm')

ggsave('out/Fig. S12. node_plot_phylum2.jpg',
       w = 88*4, h = 88*1.5, units = 'mm')

rm(t1, t2, t3)

# Betweenness corrected plots ----
# ~~~ construct functions to combine tables 
combine_node_btw_tables <- function(df1, df2, df3, df4, tax_level, tax_for_merge){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_betw = sum(node.betw), trt = 'Ctrl') 
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_betw = sum(node.betw), trt = 'Cd') 
  c <- df3 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_betw = sum(node.betw), trt = 'M')
  d <- df4 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_betw = sum(node.betw), trt = 'M+Cd')
  
  df_list <- list(a, b, c, d)
  df_merge <- Reduce(function(x, y) merge(x, y, by = tax_for_merge, all=TRUE), df_list)
  return(df_merge)
}
combine_node_btw_tables_hph <- function(df1, df2, tax_level, tax_for_merge){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_betw = sum(node.betw), trt = 'M')
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_betw = sum(node.betw), trt = 'M+Cd')
  
  df_list <- list(a, b)
  df_merge <- Reduce(function(x, y) merge(x, y, by = tax_for_merge, all=TRUE), df_list)
  return(df_merge)
}

# combine tables. Family level has been finished so convert to comment using #
# rhizosphere
rzs_node_btw_phylum <- combine_node_btw_tables(df1, df2, df3, df4, Phylum, 'Phylum')
write.csv(rzs_node_btw_phylum, 'out/rzs_node_btw_phylum.csv')
# rzs_node_btw_Family <- combine_node_btw_tables(df1, df2, df3, df4, Family, 'Family')

# hyphosphere
hps_node_btw_phylum <- combine_node_btw_tables(df5, df6, df7, df8, Phylum, 'Phylum')
write.csv(hps_node_btw_phylum, 'out/hps_node_btw_phylum.csv')
# hps_node_btw_Family <- combine_node_btw_tables(df5, df6, df7, df8, Family, 'Family')

# hyphae-associated
hph_node_btw_phylum <- combine_node_btw_tables_hph(df9, df10, Phylum, 'Phylum')
write.csv(hph_node_btw_phylum, 'out/hph_node_btw_phylum.csv')
# hph_node_btw_Family <- combine_node_btw_tables_hph(df9, df10, Family, 'Family')


# the data was obtained using R then modified in excel, then reloaded here
# load data
# rzs
rzs_family_node_btw_mod <- read.csv('data/network_data/Table S3.1. rzs_family_node.betw_mod.csv')
rzs_family_node_btw_mod$Treatment <- factor(rzs_family_node_btw_mod$Treatment, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))


rzs_phylum_node_btw_mod <- read.csv('data/network_data/rzs_node_btw_phylum_mod.csv')
rzs_phylum_node_btw_mod$Treatment <- factor(rzs_phylum_node_btw_mod$Treatment, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))


# hps
hps_family_node_btw_mod <- read.csv('data/network_data/Table S3.2. hps_family_node.betw_mod.csv')
hps_family_node_btw_mod$Treatment <- factor(hps_family_node_btw_mod$Treatment, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))


hps_phylum_node_btw_mod <- read.csv('data/network_data/hps_node_btw_phylum_mod.csv')
hps_phylum_node_btw_mod$Treatment <- factor(hps_phylum_node_btw_mod$Treatment, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))


# hph
hph_family_node_btw_mod <- read.csv('data/network_data/Table S3.3. hph_family_node.betw_mod.csv')
hph_family_node_btw_mod$Treatment <- factor(hph_family_node_btw_mod$Treatment, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))

hph_phylum_node_btw_mod <- read.csv('data/network_data/hph_node_btw_phylum_mod.csv')
hph_phylum_node_btw_mod$Treatment <- factor(hph_phylum_node_btw_mod$Treatment, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))



p4 <- ggplot(rzs_family_node_btw_mod, aes(Family, Accum_btw, color = Treatment, group = Treatment))+
  #geom_point(size = 4, alpha = 0.4)+
  geom_line()+
  ggtitle('Rhizosphere')+
  theme(legend.position = c(0.1, 0.6))+
  ylab('Accumulative betweenness')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        legend.title = element_blank())


p5 <- ggplot(hps_family_node_btw_mod, aes(Family, Accum_btw, color = Treatment, group = Treatment))+
  #geom_point(size = 4, alpha = 0.4)+
  geom_line()+
  ggtitle('Hyphosphere')+
  theme(legend.position = c(0.1, 0.6))+
  ylab('Accumulative betweenness')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        legend.title = element_blank())

p6 <- ggplot(hph_family_node_btw_mod, aes(Family, Accum_btw, color = Treatment, group = Treatment))+
  #geom_point(size = 4, alpha = 0.4)+
  geom_line()+
  ggtitle('Hyphae')+
  theme(legend.position = c(0.1, 0.8))+
  ylab('Accumulative betweenness')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        legend.title = element_blank())

# ~~~~ Fig. S15. Accumulative betweenness (family level) ----
ggarrange(p4, p5, p6,
          labels = c('a', 'b', 'c'),
          nrow = 3, ncol = 1)
if(save){
ggsave('out/Fig. S15. betweenness_plot_family2.pdf',
       w = 88*2, h = 66*3, units = 'mm')

ggsave('out/Fig. S15. betweenness_plot_family2.jpg',
       w = 88*2, h = 66*3, units = 'mm')
}
rm(p4, p5, p6)



t4 <- ggplot(rzs_phylum_node_btw_mod, aes(Phylum, Accum_betw, color = Treatment, group = Treatment))+
  geom_point(size = 4, alpha = 0.4)+
  geom_line()+
  ggtitle('Rhizosphere')+
  theme(legend.position = c(0.8, 0.5))+
  ylab('Accumulative betweenness')+
  scale_x_discrete(limits=rev)+
  coord_flip()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank())


t5 <- ggplot(hps_phylum_node_btw_mod, aes(Phylum, Accum_betw, color = Treatment, group = Treatment))+
  geom_point(size = 4, alpha = 0.4)+
  geom_line()+
  ggtitle('Hyphosphere')+
  theme(legend.position = c(0.8, 0.5))+
  ylab('Accumulative betweenness')+
  scale_x_discrete(limits=rev)+
  coord_flip()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank())

t6 <- ggplot(hph_phylum_node_btw_mod, aes(Phylum, Accum_betw, color = Treatment, group = Treatment))+
  geom_point(size = 4, alpha = 0.4)+
  geom_line()+
  ggtitle('Hyphae')+
  theme(legend.position = c(0.8, 0.7))+
  ylab('Accumulative betweenness')+
  scale_x_discrete(limits=rev)+
  coord_flip()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank())


# ~~~~ Fig. S14. Accumulative betweenness (phylum level) ----
ggarrange(t4, t5, t6,
          labels = c('a', 'b', 'c'),
          nrow = 1, ncol = 3)
if(save){
ggsave('out/Fig. S14. btw_plot_phylum2.pdf',
       w = 88*4, h = 88*1.5, units = 'mm')

ggsave('out/Fig. S14. btw_plot_phylum2.jpg',
       w = 88*4, h = 88*1.5, units = 'mm')
}
rm(t4, t5, t6)



# ~~ NO. MODULES ----
# ~~~~~~ phylum level ----
p1 <- plot_module0(df1, df2, df3, df4, Phylum)+
  ggtitle('Rhizosphere')+
  theme(legend.position = c(0.8, 0.5))+
  ylab('Number of modules')+
  theme(legend.title = element_blank())

p2 <- plot_module0(df5, df6, df7, df8, Phylum)+
  ggtitle('Hyphosphere')+
  theme(legend.position = c(0.8, 0.5))+
  ylab('Number of modules')+
  theme(legend.title = element_blank())

p3 <- plot_module_2_trt0(df9, df10, Phylum)+
  ggtitle('Hyphae')+
  theme(legend.position = c(0.8, 0.7))+
  ylab('Number of modules')+
  theme(legend.title = element_blank())

# ~~~~ Additional Figure. Number of module (phylum level) ----
ggarrange(p1, p2, p3,
          labels = c('a', 'b', 'c'),
          nrow = 1, ncol = 3)

ggsave('out/no.module_plot_phylum.pdf',
       w = 88*4, h = 88*1.5, units = 'mm')

ggsave('out/no.module_plot_phylum.jpg',
       w = 88*4, h = 88*1.5, units = 'mm')

rm(p1, p2, p3)

# ~~~~~~ family level ----

p1 <- plot_module(df1, df2, df3, df4, Family)+
  ggtitle('Rhizosphere')+
  theme(legend.position = c(0.15, 0.65))+
  ylab('Number of modules')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        legend.title = element_blank())

p2 <- plot_module(df5, df6, df7, df8, Family)+
  ggtitle('Hyphosphere')+
  theme(legend.position = c(0.15, 0.65))+
  ylab('Number of modules')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        legend.title = element_blank())

p3 <- plot_module_2_trt(df9, df10, Family)+
  ggtitle('Hyphae')+
  theme(legend.position = c(0.15, 0.65))+
  ylab('Number of modules')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        legend.title = element_blank())

# ~~~~ Additional Figure. Number of module (family level) ----
ggarrange(p1, p2, p3,
          labels = c('a', 'b', 'c'),
          nrow = 3, ncol = 1)

ggsave('out/no.module_plot_family.pdf',
       w = 88*2, h = 66*3, units = 'mm')

ggsave('out/no.module_plot_family.jpg',
       w = 88*2, h = 66*3, units = 'mm')

rm(p1, p2, p3)




# # igraph to calculate network properties ----
# # ~~~~ Additional Fig. Modularity of networks ---- 
# # ~~~ rhizosphere ----
# png(filename = "out/add_Fig.1-4. Modularity.png", width = 2160, height = 2160, units = "px", res = 300)
# par(mfrow = c(2, 2), mar = c(2, 2, 2, 2), oma = c(1, 1, 1, 1)) # Adjust the fourth value (right margin) as needed
# 
# 
# rzs_ctrl_edge <- read.csv('data/network_data/rzs_majority_7/rzs_ctrl_edge.csv', header = TRUE)
# # prefix "Zotu" for each cell of column 'Source' and 'Target'
# rzs_ctrl_edge$Source <- paste0('Zotu', rzs_ctrl_edge$Source)
# rzs_ctrl_edge$Target <- paste0('Zotu', rzs_ctrl_edge$Target)
# head(rzs_ctrl_edge)
# 
# rzs_ctrl_node <- read.csv('data/network_data/rzs_majority_7/rzs_ctrl_node_final.csv', header = TRUE)
# # remove first column 'ID'
# rzs_ctrl_node <- rzs_ctrl_node[, -1]
# head(rzs_ctrl_node)
# 
# 
# rzs_ctrl_nw <- as_data_frame(d = rzs_ctrl_edge, vertices = rzs_ctrl_node, directed = FALSE)
# 
# rzs_ctrl_nw <- graph_from_data_frame(d = rzs_ctrl_edge, vertices = rzs_ctrl_node, directed = FALSE)
# 
# plot(rzs_ctrl_nw, vertex.size=5, vertex.label=NA)
# # title(main = "Rhizosphere - Ctrl")
# # modularity(cluster_fast_greedy(rzs_ctrl_nw))
# mtext(paste("Rhizosphere (Ctrl)\nmodularity =", round(modularity(cluster_fast_greedy(rzs_ctrl_nw)), 4)), side = 3, line = 0.5)
# 
# # Modularity is a measure of the structure of networks or graphs which measures the strength of division of a network into modules. Networks with high modularity have dense connections between the nodes within modules but sparse connections between nodes in different modules. Wikipedia
# 
# 
# rzs_cd_edge <- read.csv('data/network_data/rzs_majority_7/rzs_cd_edge.csv', 
#                           header = TRUE)
# rzs_cd_nw <- graph.data.frame(rzs_cd_edge, directed = FALSE)
# plot(rzs_cd_nw, vertex.size=5, vertex.label=NA)
# # title(main = "Rhizosphere - Cd")
# # modularity(cluster_fast_greedy(rzs_cd_nw))
# mtext(paste("Rhizosphere (Cd)\nmodularity =", round(modularity(cluster_fast_greedy(rzs_cd_nw)), 4)), side = 3, line = 0.5)
# 
# 
# rzs_m_edge <- read.csv('data/network_data/rzs_majority_7/rzs_m_edge.csv', 
#                         header = TRUE)
# rzs_m_nw <- graph.data.frame(rzs_m_edge, directed = FALSE)
# plot(rzs_m_nw, vertex.size=5, vertex.label=NA)
# # title(main = "Rhizosphere - M")
# # modularity(cluster_fast_greedy(rzs_m_nw))
# mtext(paste("Rhizosphere (M)\nmodularity =", round(modularity(cluster_fast_greedy(rzs_m_nw)), 4)), side = 3, line = 0.5)
# 
# 
# rzs_m_cd_edge <- read.csv('data/network_data/rzs_majority_7/rzs_m_cd_edge.csv', 
#                        header = TRUE)
# rzs_m_cd_nw <- graph.data.frame(rzs_m_cd_edge, directed = FALSE)
# plot(rzs_m_cd_nw, vertex.size=5, vertex.label=NA)
# # title(main = "Rhizosphere - M+Cd")
# # modularity(cluster_fast_greedy(rzs_m_cd_nw))
# mtext(paste("Rhizosphere (M+Cd)\nmodularity =", round(modularity(cluster_fast_greedy(rzs_m_cd_nw)), 4)), side = 3, line = 0.5)
# 
# 
# par(mfrow = c(1, 1))
# 
# dev.off()  
# 
# # ~~~ hyphosphere ----
# png(filename = "out/add_Fig.5-8. Modularity.png", width = 2160, height = 2160, units = "px", res = 300)
# par(mfrow = c(2, 2), mar = c(2, 2, 2, 2), oma = c(1, 1, 1, 1)) # Adjust the fourth value (right margin) as needed
# 
# 
# hps_ctrl_edge <- read.csv('data/network_data/hps_majority_6/hps_ctrl_edge.csv', 
#                           header = TRUE)
# hps_ctrl_nw <- graph.data.frame(hps_ctrl_edge, directed = FALSE)
# plot(hps_ctrl_nw, vertex.size=5, vertex.label=NA)
# # modularity(cluster_fast_greedy(hps_ctrl_nw))
# mtext(paste("Hyphosphere (Ctrl)\nmodularity =", round(modularity(cluster_fast_greedy(hps_ctrl_nw)), 4)), side = 3, line = 0.5)
# 
# 
# hps_cd_edge <- read.csv('data/network_data/hps_majority_6/hps_cd_edge.csv', 
#                           header = TRUE)
# hps_cd_nw <- graph.data.frame(hps_cd_edge, directed = FALSE)
# plot(hps_cd_nw, vertex.size=5, vertex.label=NA)
# # modularity(cluster_fast_greedy(hps_cd_nw))
# mtext(paste("Hyphosphere (Cd)\nmodularity =", round(modularity(cluster_fast_greedy(hps_cd_nw)), 4)), side = 3, line = 0.5)
# 
# 
# 
# hps_m_edge <- read.csv('data/network_data/hps_majority_6/hps_m_edge.csv', 
#                        header = TRUE)
# hps_m_nw <- graph.data.frame(hps_m_edge, directed = FALSE)
# plot(hps_m_nw, vertex.size=5, vertex.label=NA)
# # modularity(cluster_fast_greedy(hps_m_nw))
# mtext(paste("Hyphosphere (M)\nmodularity =", round(modularity(cluster_fast_greedy(hps_m_nw)), 4)), side = 3, line = 0.5)
# 
# 
# 
# hps_m_cd_edge <- read.csv('data/network_data/hps_majority_6/hps_m_cd_edge.csv', 
#                        header = TRUE)
# hps_m_cd_nw <- graph.data.frame(hps_m_cd_edge, directed = FALSE)
# plot(hps_m_cd_nw, vertex.size=5, vertex.label=NA)
# # modularity(cluster_fast_greedy(hps_m_cd_nw))
# mtext(paste("Hyphosphere (M+Cd)\nmodularity =", round(modularity(cluster_fast_greedy(hps_m_cd_nw)), 4)), side = 3, line = 0.5)
# 
# 
# par(mfrow = c(1, 1))
# dev.off()
# 
# # ~~~ hyphae-associated ----
# png(filename = "out/add_Fig.9-10. Modularity.png", width = 2160, height = 1080, units = "px", res = 300)
# par(mfrow = c(1, 2), mar = c(2, 2, 2, 2), oma = c(1, 1, 1, 1)) # Adjust the fourth value (right margin) as needed
# 
# hph_m_edge <- read.csv('data/network_data/hph_m_edge_final.csv', 
#                        header = TRUE)
# hph_m_nw <- graph.data.frame(hph_m_edge, directed = FALSE)
# plot(hph_m_nw, vertex.size=5, vertex.label=NA)
# # modularity(cluster_fast_greedy(hph_m_nw))
# mtext(paste("Hyphae (M)\nmodularity =", round(modularity(cluster_fast_greedy(hph_m_nw)), 4)), side = 3, line = 0.5)
# 
# 
# 
# hph_m_cd_edge <- read.csv('data/network_data/hph_m_cd_edge_final.csv', 
#                           header = TRUE)
# hph_m_cd_nw <- graph.data.frame(hph_m_cd_edge, 
#                                 # correlation=hph_m_cd_edge$Label,
#                                 directed = FALSE)
# plot(hph_m_cd_nw, vertex.size=5, vertex.label=NA)
# # modularity(cluster_fast_greedy(hph_m_cd_nw))
# mtext(paste("Hyphae (M+Cd)\nmodularity =", round(modularity(cluster_fast_greedy(hph_m_cd_nw)), 4)), side = 3, line = 0.5)
# 
# par(mfrow = c(1, 1))
# dev.off()




# End 



