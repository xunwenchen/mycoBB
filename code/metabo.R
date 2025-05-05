# metabolomic analysis -------
#  samples were collected from the rhizosphere (A, B, C, D) and hyphoshere (HA, HB, HC, HD)
# A is ctrl
# B is +AMF
# c is +AMF+Cd
# D is Cd
# H suffixed means Hyphosphere


# load metabo data ----
metabo_data_neg <- read_xlsx('data/all.exp.NEG_com_name.xlsx') 
metabo_data_pos <- read_xlsx('data/all.exp.POS_com_name.xlsx')
# the above tables were modified from all.exp.NEG.xlsx and all.exp.POS.xlsx by manually melting the MS1_name and MS2_name, with two new columns, com_name and name_type. The com_name shall be used to merge NEG and POS tables. 


# There were 566 compounds detected using the negative ion mode, and 1078 compounds detected using the positive ion mode. These two lists of compounds were combined for subsequent analysis. For the compounds presented in both modes, the ones with higher concentrations were retained while the other ones were omitted. After combining, the total number of compounds is 1,644 indicating that no overlapping compounds from the two different detection modes. 

# merge two tables of compound concentration derived from positive and negative ion modes ----
metabo_combined <- merge(metabo_data_neg, metabo_data_pos, 
                         by = 'com_name', all = TRUE)

metabo_df0 <- metabo_combined %>% 
  summarise(com_id = coalesce(`index.x`, `index.y`), 
            # Here should be summariSe, not summariZe!!!! see here: https://github.com/tidyverse/dplyr/issues/505
            com_name = `com_name`,
        # Rhizosphere compartment -> control
        rs_ck_1 = pmax(`A-1.x`, `A-1.y`, na.rm = T),
         rs_ck_2 = pmax(`A-2.x`, `A-2.y`, na.rm = T),
         rs_ck_3 = pmax(`A-3.x`, `A-3.y`, na.rm = T),
         rs_ck_4 = pmax(`A-4.x`, `A-4.y`, na.rm = T),
         rs_ck_5 = pmax(`A-5.x`, `A-5.y`, na.rm = T),
         rs_ck_6 = pmax(`A-6.x`, `A-6.y`, na.rm = T),
         rs_ck_7 = pmax(`A-7.x`, `A-7.y`, na.rm = T),
         rs_ck_8 = pmax(`A-8.x`, `A-8.y`, na.rm = T),
         # Rhizosphere compartment -> AMF
         rs_f_1 = pmax(`B-1.x`, `B-1.y`, na.rm = T),
         rs_f_2 = pmax(`B-2.x`, `B-2.y`, na.rm = T),
         rs_f_3 = pmax(`B-3.x`, `B-3.y`, na.rm = T),
         rs_f_4 = pmax(`B-4.x`, `B-4.y`, na.rm = T),
         rs_f_5 = pmax(`B-5.x`, `B-5.y`, na.rm = T),
         rs_f_6 = pmax(`B-6.x`, `B-6.y`, na.rm = T),
         rs_f_7 = pmax(`B-7.x`, `B-7.y`, na.rm = T),
         rs_f_8 = pmax(`B-8.x`, `B-8.y`, na.rm = T),
         # Rhizosphere compartment -> AMF+Cd
         rs_f_cd_1 = pmax(`C-1.x`, `C-1.y`, na.rm = T),
         rs_f_cd_2 = pmax(`C-2.x`, `C-2.y`, na.rm = T),
         rs_f_cd_3 = pmax(`C-3.x`, `C-3.y`, na.rm = T),
         rs_f_cd_4 = pmax(`C-4.x`, `C-4.y`, na.rm = T),
         rs_f_cd_5 = pmax(`C-5.x`, `C-5.y`, na.rm = T),
         rs_f_cd_6 = pmax(`C-6.x`, `C-6.y`, na.rm = T),
         rs_f_cd_7 = pmax(`C-7.x`, `C-7.y`, na.rm = T),
         rs_f_cd_8 = pmax(`C-8.x`, `C-8.y`, na.rm = T),
         # Rhizosphere compartment -> Cd
         rs_cd_1 = pmax(`D-1.x`, `D-1.y`, na.rm = T),
         rs_cd_2 = pmax(`D-2.x`, `D-2.y`, na.rm = T),
         rs_cd_3 = pmax(`D-3.x`, `D-3.y`, na.rm = T),
         rs_cd_4 = pmax(`D-4.x`, `D-4.y`, na.rm = T),
         rs_cd_5 = pmax(`D-5.x`, `D-5.y`, na.rm = T),
         rs_cd_6 = pmax(`D-6.x`, `D-6.y`, na.rm = T),
         rs_cd_7 = pmax(`D-7.x`, `D-7.y`, na.rm = T),
         rs_cd_8 = pmax(`D-8.x`, `D-8.y`, na.rm = T),
         
         
         # Hyphosphere compartment -> control
         hps_ck_1 = pmax(`HA-1.x`, `HA-1.y`, na.rm = T),
         hps_ck_2 = pmax(`HA-2.x`, `HA-2.y`, na.rm = T),
         hps_ck_3 = pmax(`HA-3.x`, `HA-3.y`, na.rm = T),
         hps_ck_4 = pmax(`HA-4.x`, `HA-4.y`, na.rm = T),
         hps_ck_5 = pmax(`HA-5.x`, `HA-5.y`, na.rm = T),
         hps_ck_6 = pmax(`HA-6.x`, `HA-6.y`, na.rm = T),
         hps_ck_7 = pmax(`HA-7.x`, `HA-7.y`, na.rm = T),
         hps_ck_8 = pmax(`HA-8.x`, `HA-8.y`, na.rm = T),
         # Hyphosphere compartment -> AMF
         hps_f_1 = pmax(`HB-1.x`, `HB-1.y`, na.rm = T),
         hps_f_2 = pmax(`HB-2.x`, `HB-2.y`, na.rm = T),
         hps_f_3 = pmax(`HB-3.x`, `HB-3.y`, na.rm = T),
         hps_f_4 = pmax(`HB-4.x`, `HB-4.y`, na.rm = T),
         hps_f_5 = pmax(`HB-5.x`, `HB-5.y`, na.rm = T),
         hps_f_6 = pmax(`HB-6.x`, `HB-6.y`, na.rm = T),
         hps_f_7 = pmax(`HB-7.x`, `HB-7.y`, na.rm = T),
         hps_f_8 = pmax(`HB-8.x`, `HB-8.y`, na.rm = T),
         # Hyphosphere compartment -> AMF+Cd
         hps_f_cd_1 = pmax(`HC-1.x`, `HC-1.y`, na.rm = T),
         hps_f_cd_2 = pmax(`HC-2.x`, `HC-2.y`, na.rm = T),
         hps_f_cd_3 = pmax(`HC-3.x`, `HC-3.y`, na.rm = T),
         hps_f_cd_4 = pmax(`HC-4.x`, `HC-4.y`, na.rm = T),
         hps_f_cd_5 = pmax(`HC-5.x`, `HC-5.y`, na.rm = T),
         hps_f_cd_6 = pmax(`HC-6.x`, `HC-6.y`, na.rm = T),
         hps_f_cd_7 = pmax(`HC-7.x`, `HC-7.y`, na.rm = T),
         hps_f_cd_8 = pmax(`HC-8.x`, `HC-8.y`, na.rm = T),
         # Hyphosphere compartment -> Cd
         hps_cd_1 = pmax(`HD-1.x`, `HD-1.y`, na.rm = T),
         hps_cd_2 = pmax(`HD-2.x`, `HD-2.y`, na.rm = T),
         hps_cd_3 = pmax(`HD-3.x`, `HD-3.y`, na.rm = T),
         hps_cd_4 = pmax(`HD-4.x`, `HD-4.y`, na.rm = T),
         hps_cd_5 = pmax(`HD-5.x`, `HD-5.y`, na.rm = T),
         hps_cd_6 = pmax(`HD-6.x`, `HD-6.y`, na.rm = T),
         hps_cd_7 = pmax(`HD-7.x`, `HD-7.y`, na.rm = T),
         hps_cd_8 = pmax(`HD-8.x`, `HD-8.y`, na.rm = T))
   if(save){write.csv(metabo_df0, 'data/metabo_data/metabo_df0.csv')}    

# ANOVA for each compound across treatments ----
View(metabo_df0)

# Find the most different compounds ----

# Find the most similar compounds ----

# compound correlation and cluster ----


# pattern analysis, but our treatment is not a ladder such as dosage 1 ppm, 5 ppm, 10 ppm. So cannot do correlation between each compound with dosage? may be we can take ctrl, AMF, Cd, AMF+Cd like a dosage? should be not. 

# PCA/PCoA ----

# ~~~~ metabo PCoA - Rhizosphere ----
metabo_df_for_PCoA <- subset(metabo_df0, select = c(3:34))

     
# construct bray curtis matrix
    metabo_bc_dist <- vegdist(t(metabo_df_for_PCoA), method = "bray")
     metabo_bc_pcoa <- pco(metabo_bc_dist)
     metabo_bc_pcoa_df <- data.frame(pcoa1 = metabo_bc_pcoa$vectors[,1], 
                                       pcoa2 = metabo_bc_pcoa$vectors[,2])
     
     
     metabo_bc_pcoa_df_trt <- cbind(metabo_bc_pcoa_df,
                                    trt = c(
                                      'Ctrl', 'Ctrl','Ctrl','Ctrl',
                                      'Ctrl','Ctrl','Ctrl','Ctrl',
                                      'M','M','M','M',
                                      'M','M','M','M', 
                                      'M+Cd','M+Cd','M+Cd','M+Cd',
                                      'M+Cd','M+Cd','M+Cd','M+Cd', 
                                      'Cd', 'Cd', 'Cd', 'Cd', 
                                      'Cd', 'Cd', 'Cd', 'Cd'))
     metabo_bc_pcoa_df_trt$trt <- factor(metabo_bc_pcoa_df_trt$trt,
                                         levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))

# Extract the loading scores of the first two axes
loadings_rhz <- eigenvals(metabo_bc_pcoa)[1:2]*100
loadings_rhz

p1 <- ggplot(metabo_bc_pcoa_df_trt, aes(x=pcoa1, y=pcoa2, 
                                       color = trt, shape = trt)) +
       geom_point(size = 2) +
  labs(x = paste0("PC1 (", round(loadings_rhz[1], 1), "%)"),
       y = paste0("PC2 (", round(loadings_rhz[2], 1), "%)"), 
            title = "Rhizosphere") +
  theme(legend.position = 'none')+
  scale_color_manual(values = c("#bababa","#404040", "#f4a582", "#ca0020"))
 
# ~~~~ metabo PCoA - hyphosphere ----
     
metabo_df_for_PCoA_hps <- subset(metabo_df0, select = c(35:66))
metabo_bc_dist_hps <- vegdist(t(metabo_df_for_PCoA_hps), method = "bray")
metabo_bc_pcoa_hps <- pco(metabo_bc_dist_hps)
metabo_bc_pcoa_df_hps <- data.frame(pcoa1 = metabo_bc_pcoa_hps$vectors[,1], 
                                pcoa2 = metabo_bc_pcoa_hps$vectors[,2])

metabo_bc_pcoa_df_hps_trt <- cbind(metabo_bc_pcoa_df_hps,
                               trt = c(
                                 'Ctrl', 'Ctrl','Ctrl','Ctrl',
                                 'Ctrl','Ctrl','Ctrl','Ctrl',
                                 'M','M','M','M',
                                 'M','M','M','M', 
                                 'M+Cd','M+Cd','M+Cd','M+Cd',
                                 'M+Cd','M+Cd','M+Cd','M+Cd', 
                                 'Cd', 'Cd', 'Cd', 'Cd', 
                                 'Cd', 'Cd', 'Cd', 'Cd'))
metabo_bc_pcoa_df_hps_trt$trt <- factor(metabo_bc_pcoa_df_hps_trt$trt,
                                    levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))

# Extract the loading scores of the first two axes
loadings_hps <- eigenvals(metabo_bc_pcoa_hps)[1:2]
loadings_hps

p2 <- ggplot(metabo_bc_pcoa_df_hps_trt, aes(x=pcoa1, y=pcoa2, 
                                  color = trt, shape = trt)) +
  geom_point(size = 2) +
  labs(x = paste0("PC1 (", round(loadings_hps[1], 1), "%)"),
       y = paste0("PC2 (", round(loadings_hps[2], 1), "%)"),
       title = "Hyphosphere") +
  theme(legend.position = c(0.84, 0.65),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  scale_color_manual(values = c("#bababa","#404040", "#f4a582", "#ca0020"))

# ~~~~ Extended Data Fig. 3. The PCoA of metabolome B-C dist ----
p_metab_pcoa <- ggarrange(p1, p2,
          labels = c('a', 'b'),
          ncol = 2, nrow = 1)

if(save){
  ggsave('out/ExtDataFig. 3. metabo_PCoA.pdf', plot = p_metab_pcoa,
       w = 88*2, h = 66, units = 'mm')
  ggsave('out/ExtDataFig. 3. metabo_PCoA.jpg', plot = p_metab_pcoa,
       w = 88*2, h = 66, units = 'mm')
}

rm(p1, p2)

# ~~~~ Extended Data Fig. 3. Pairwise adonis - rhizosphere ----
pairwiseAdonis::pairwise.adonis(metabo_bc_dist, metabo_bc_pcoa_df_trt$trt, perm = 999)

# ~~~~ Extended Data Fig. 3. Pairwise adonis - hyphosphere ----
pairwiseAdonis::pairwise.adonis(metabo_bc_dist_hps, metabo_bc_pcoa_df_hps_trt$trt, perm = 999)



# PCoA for metabo hyphosphere vs rhizosphere ------
# build a function first
metabo_pco_plot <- function(data, title) {
  
  # Calculate distance matrix
  dist_matrix <- vegdist(t(data), method = "bray")
  
  # Perform PCoA
  pcoa_result <- pco(dist_matrix)
  
  # Extract the loading scores of the first two axes
  loadings <- eigenvals(pcoa_result)[1:2]
  
  # Create data frame with PCoA coordinates
  pcoa_df <- data.frame(pcoa1 = pcoa_result$vectors[,1], pcoa2 = pcoa_result$vectors[,2])
  
  # Add compartment column
  compartments <- c(rep('Rhizosphere', 8), rep('Hyphosphere', 8))
  pcoa_df <- cbind(pcoa_df, compartment = compartments)
  
  # Create plot
  plot <- ggplot(pcoa_df, aes(x = pcoa1, y = pcoa2, color = compartment)) + 
    geom_point(size = 2) + 
    labs(x = paste0("PC1 (", round(loadings[1]/sum(loadings)*100,2), "%)"), 
         y = paste0("PC2 (", round(loadings[2]/sum(loadings)*100,2), "%)"), 
         title = title)
  
  return(plot)
}

# extract data then plot 
metabo_df_temp_ctrl <- subset(metabo_df0, select = c(3:10, 35:42))
p1 <- metabo_pco_plot(metabo_df_temp_ctrl, 'Ctrl')

metabo_df_temp_cd <- subset(metabo_df0, select = c(27:(27+7), 59:(59+7)))
p2 <- metabo_pco_plot(metabo_df_temp_cd, 'Cd')

metabo_df_temp_m <- subset(metabo_df0, select = c(11:(11+7), 43:(43+7)))
p3 <- metabo_pco_plot(metabo_df_temp_m, 'M')

metabo_df_temp_m_cd <- subset(metabo_df0, select = c(19:(19+7), 51:(51+7)))
p4 <- metabo_pco_plot(metabo_df_temp_m_cd, 'M+Cd')

# ~~~~ Additional Figure. PCoA for metabo hyphosphere vs rhizosphere ----
p_pcoa_metabo_h_vs_r <- ggarrange(p1, p2, p3, p4,
          labels = c('a', 'b', 'c', 'd'),
          ncol = 2, nrow = 2, legend = "none")

if(save){
ggsave('out/metabo_pcoa_comp.jpg', plot = p_pcoa_metabo_h_vs_r,
       w =88*2, h=66*2, unit= 'mm', dpi=300)
ggsave('out/metabo_pcoa_comp.pdf', plot = p_pcoa_metabo_h_vs_r,
       w =88*2, h=66*2, unit= 'mm', dpi=300)
}
rm(p1, p2, p3, p4)

# TEST #####################
############################
dummy_df <- read_xlsx('data/dummy_data.xlsx')
dummy_df_dist <- vegdist(t(dummy_df), method = "bray")
dummy_bc_pcoa <- pco(dummy_df_dist)
dummy_bc_pcoa_df <- data.frame(pcoa1 = dummy_bc_pcoa$vectors[,1], 
                                    pcoa2 = dummy_bc_pcoa$vectors[,2])

dummy_bc_pcoa_df_trt <- cbind(dummy_bc_pcoa_df,
                                   trt = c(
                                     'Ctrl', 'Ctrl','Ctrl','Ctrl', 'Ctrl','Ctrl',
                                     'trt_a','trt_a','trt_a','trt_a', 'trt_a','trt_a',
                                     'trt_b','trt_b','trt_b','trt_b','trt_b','trt_b', 
                                     'trt_c', 'trt_c', 'trt_c', 'trt_c', 'trt_c', 'trt_c'))
dummy_bc_pcoa_df_trt$trt <- factor(dummy_bc_pcoa_df_trt$trt,
                                        levels = c('Ctrl', 'trt_a', 'trt_b', 'trt_c'))

ggplot(dummy_bc_pcoa_df_trt, aes(x=pcoa1, y=pcoa2, color = trt, shape = trt)) +
  geom_point(size = 2) +
  labs(x = "PC1", y = "PC2", title = "TEST")
# TEST END #####################
############################
     
# ~~~~ statistical test for PCoA ----

pairwiseAdonis::pairwise.adonis(dummy_df_dist, dummy_bc_pcoa_df_trt$trt, perm = 999)

########################################
########################################


    
# partial least squares-discriminant analysis (PLS-DA) ----
# use PLS-DA option to view the separation of the groups. PLS-DA 'rotates' the PCA axes to maximize separation. Look at the 2D PLS scores plot. Look at the Q^2 and R^2 values (cross validation). Use the VIP plot to ID (identify) important metabolites (VIP>1.2). 
# The above info is from Bioinformatics.ca the youtube video https://www.youtube.com/watch?v=HNzAe_SDYv0. 
# For later enrichment analysis and pathway analysis, seems MetaboAnylst can only do human stuffs. 
# how about soil? What method should be used then?



# chemodiversity correlate with microbial diversity? ----

# may use two compartments and four treatments (hyphal samples only two treatments)
# to do correlation analysis. 
# array 1: chemdiversity of samples: A1-8, B1-8, C1-8, D1-8, HA1-8, HB1-8, HC1-8, HD1-8 = 8*8 = 64 Samples
# array 2: bacterial diversity of samples: rhz

# metabo_df add H/C and O/C ratios data ----
metabo_hc_oc <- read_xlsx('data/metabo_data/HC_OC_ratio.xlsx')
# combine with metabo_df
     
metabo_df <- merge(metabo_hc_oc, metabo_df0, by = 'com_id', all = T)
write.csv(metabo_df, 'data/metabo_data/metabo_df.csv')

ggplot(metabo_df, aes(OC_ratio, HC_ratio))+
  geom_point(alpha = 0.3)+
  xlim(0, 1.2)+
  ylim(0, 3)

# compare with Zhang Peng's runzhang
ft_df1 <- read_xlsx('data/metabo_data/FT-ICR_MS_ALL_Formula_DATA.xlsx')

library(purrr)
ft_df1_z <- ft_df1 %>%
  mutate(across(c(starts_with("rs"), starts_with("hps")), list(z = ~ scale(.))))

ft_df1_z_avg <- ft_df1_z %>%
  mutate(rs_ck_avg_z = rowMeans(select(.,93:100), na.rm = TRUE),
         rs_f_avg_z = rowMeans(select(., 101:108), na.rm = TRUE),
         rs_f_cd_avg_z = rowMeans(select(., 109:116), na.rm = TRUE),
         rs_cd_avg_z = rowMeans(select(., 117:124), na.rm = TRUE),
         hps_ck_avg_z = rowMeans(select(., 125:132), na.rm = TRUE),
         hps_f_avg_z = rowMeans(select(., 133:140), na.rm = TRUE),
         hps_f_cd_avg_z = rowMeans(select(., 141:148), na.rm = TRUE),
         hps_cd_avg_z = rowMeans(select(., 149:156), na.rm = TRUE))


# ft_df1_avg <- ft_df1 %>% mutate(
#   rs_ck_avg = rowMeans(select(., 10:17), na.rm = TRUE),
#   rs_f_avg = rowMeans(select(., 18:25), na.rm = TRUE),
#   rs_f_cd_avg = rowMeans(select(., 26:33), na.rm = TRUE),
#   rs_cd_avg = rowMeans(select(., 34:41), na.rm = TRUE),
#   
#   hps_ck_avg = rowMeans(select(., 42:49), na.rm = TRUE),
#   hps_f_avg = rowMeans(select(., 50:57), na.rm = TRUE),
#   hps_f_cd_avg = rowMeans(select(., 58:65), na.rm = TRUE),
#   hps_cd_avg = rowMeans(select(., 66:73), na.rm = TRUE))




ggplot(ft_df1_z_avg, aes(`O/C`, `H/C`))+
  geom_point(alpha = 0.3)+
  xlim(0, 1.2)+
  ylim(0, 3)


# so those omitted points should be those that are not fulfill the FT-ICR MS data processing requirements. Zhang Peng had them omitted. - GOOD ~~ 


# draw frames to indicate types of organic compounds such as 

# ~~~~ Addtional Figure. V-K plot ----
ggplot(ft_df1_z_avg, aes(`O/C`, `H/C`))+
  geom_point(alpha = 0.3, size = 0.5)+
  xlim(-0.1, 1.2)+
  ylim(0, 2.5)+
  geom_rect(aes(xmin = 0, xmax = 0.2, ymin = 1.5, ymax = 2.3), alpha = 0, color = "black", linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.2, xmax = 0.52, ymin = 1.5, ymax = 2.2), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.52, xmax = 0.7, ymin = 1.5, ymax = 2.3), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.7, xmax = 1.1, ymin = 1.5, ymax = 2.4), alpha = 0, color = "black",  linetype='dotted',  linewidth = 0.2)+
  geom_rect(aes(xmin = 0, xmax = 0.25, ymin = 0.5, ymax = 1.25), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.25, xmax = 0.67, ymin = 0.75, ymax = 1.5), alpha = 0, color = "black", linetype='dotted',  linewidth = 0.2)+
  geom_rect(aes(xmin = 0.67, xmax = 0.97, ymin = 0.53, ymax = 1.5), alpha = 0, color = "black", linetype='dotted',  linewidth = 0.2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
if(save){
  ggsave('out/vk_plot.pdf',
       w = 88, h = 66, units = 'mm')
ggsave('out/vk_plot.jpg',
       w = 88, h = 66, units = 'mm')
}

vk_p1 <- ggplot(ft_df1_z_avg, aes(x=`O/C`, y=`H/C`, color=`rs_ck_avg_z`)) +
  geom_point(alpha=0.6) +
  scale_color_gradient(low = "blue", high = "red") +
  #scale_size(range = c(1, 20)) +# Adjust the range based on your preference
  xlim(-0.1, 1.2)+
  ylim(0, 2.5)+
  geom_rect(aes(xmin = 0, xmax = 0.2, ymin = 1.5, ymax = 2.3), alpha = 0, color = "black", linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.2, xmax = 0.52, ymin = 1.5, ymax = 2.2), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.52, xmax = 0.7, ymin = 1.5, ymax = 2.3), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.7, xmax = 1.1, ymin = 1.5, ymax = 2.4), alpha = 0, color = "black",  linetype='dotted',  linewidth = 0.2)+
  geom_rect(aes(xmin = 0, xmax = 0.25, ymin = 0.5, ymax = 1.25), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.25, xmax = 0.67, ymin = 0.75, ymax = 1.5), alpha = 0, color = "black", linetype='dotted',  linewidth = 0.2)+
  geom_rect(aes(xmin = 0.67, xmax = 0.97, ymin = 0.53, ymax = 1.5), alpha = 0, color = "black", linetype='dotted',  linewidth = 0.2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle('Rhizosphere (Ctrl)') 

vk_p2 <- ggplot(ft_df1_z_avg, aes(x=`O/C`, y=`H/C`, color=`rs_cd_avg_z`)) +
  geom_point(alpha=0.6) +
  scale_color_gradient(low = "blue", high = "red") +
  #scale_size(range = c(1, 20)) +# Adjust the range based on your preference
  xlim(-0.1, 1.2)+
  ylim(0, 2.5)+
  geom_rect(aes(xmin = 0, xmax = 0.2, ymin = 1.5, ymax = 2.3), alpha = 0, color = "black", linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.2, xmax = 0.52, ymin = 1.5, ymax = 2.2), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.52, xmax = 0.7, ymin = 1.5, ymax = 2.3), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.7, xmax = 1.1, ymin = 1.5, ymax = 2.4), alpha = 0, color = "black",  linetype='dotted',  linewidth = 0.2)+
  geom_rect(aes(xmin = 0, xmax = 0.25, ymin = 0.5, ymax = 1.25), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.25, xmax = 0.67, ymin = 0.75, ymax = 1.5), alpha = 0, color = "black", linetype='dotted',  linewidth = 0.2)+
  geom_rect(aes(xmin = 0.67, xmax = 0.97, ymin = 0.53, ymax = 1.5), alpha = 0, color = "black", linetype='dotted',  linewidth = 0.2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle('Rhizosphere (Cd)') 

vk_p3 <- ggplot(ft_df1_z_avg, aes(x=`O/C`, y=`H/C`, color=`rs_f_avg_z`)) +
  geom_point(alpha=0.6) +
  scale_color_gradient(low = "blue", high = "red") +
  # scale_size(range = c(1, 20)) +# Adjust the range based on your preference
  xlim(-0.1, 1.2)+
  ylim(0, 2.5)+
  geom_rect(aes(xmin = 0, xmax = 0.2, ymin = 1.5, ymax = 2.3), alpha = 0, color = "black", linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.2, xmax = 0.52, ymin = 1.5, ymax = 2.2), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.52, xmax = 0.7, ymin = 1.5, ymax = 2.3), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.7, xmax = 1.1, ymin = 1.5, ymax = 2.4), alpha = 0, color = "black",  linetype='dotted',  linewidth = 0.2)+
  geom_rect(aes(xmin = 0, xmax = 0.25, ymin = 0.5, ymax = 1.25), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.25, xmax = 0.67, ymin = 0.75, ymax = 1.5), alpha = 0, color = "black", linetype='dotted',  linewidth = 0.2)+
  geom_rect(aes(xmin = 0.67, xmax = 0.97, ymin = 0.53, ymax = 1.5), alpha = 0, color = "black", linetype='dotted',  linewidth = 0.2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle('Rhizosphere (M)') 

vk_p4 <- ggplot(ft_df1_z_avg, aes(x=`O/C`, y=`H/C`, color=`rs_f_cd_avg_z`)) +
  geom_point(alpha=0.6) +
  scale_color_gradient(low = "blue", high = "red") +
  # scale_size(range = c(1, 20)) +# Adjust the range based on your preference
  xlim(-0.1, 1.2)+
  ylim(0, 2.5)+
  geom_rect(aes(xmin = 0, xmax = 0.2, ymin = 1.5, ymax = 2.3), alpha = 0, color = "black", linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.2, xmax = 0.52, ymin = 1.5, ymax = 2.2), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.52, xmax = 0.7, ymin = 1.5, ymax = 2.3), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.7, xmax = 1.1, ymin = 1.5, ymax = 2.4), alpha = 0, color = "black",  linetype='dotted',  linewidth = 0.2)+
  geom_rect(aes(xmin = 0, xmax = 0.25, ymin = 0.5, ymax = 1.25), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.25, xmax = 0.67, ymin = 0.75, ymax = 1.5), alpha = 0, color = "black", linetype='dotted',  linewidth = 0.2)+
  geom_rect(aes(xmin = 0.67, xmax = 0.97, ymin = 0.53, ymax = 1.5), alpha = 0, color = "black", linetype='dotted',  linewidth = 0.2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle('Rhizosphere (M+Cd)') 





vk_p5 <- ggplot(ft_df1_z_avg, aes(x=`O/C`, y=`H/C`, color=`hps_ck_avg_z`)) +
  geom_point(alpha=0.6) +
  scale_color_gradient(low = "blue", high = "red") +
  # scale_size(range = c(1, 20)) +# Adjust the range based on your preference
  xlim(-0.1, 1.2)+
  ylim(0, 2.5)+
  geom_rect(aes(xmin = 0, xmax = 0.2, ymin = 1.5, ymax = 2.3), alpha = 0, color = "black", linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.2, xmax = 0.52, ymin = 1.5, ymax = 2.2), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.52, xmax = 0.7, ymin = 1.5, ymax = 2.3), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.7, xmax = 1.1, ymin = 1.5, ymax = 2.4), alpha = 0, color = "black",  linetype='dotted',  linewidth = 0.2)+
  geom_rect(aes(xmin = 0, xmax = 0.25, ymin = 0.5, ymax = 1.25), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.25, xmax = 0.67, ymin = 0.75, ymax = 1.5), alpha = 0, color = "black", linetype='dotted',  linewidth = 0.2)+
  geom_rect(aes(xmin = 0.67, xmax = 0.97, ymin = 0.53, ymax = 1.5), alpha = 0, color = "black", linetype='dotted',  linewidth = 0.2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle('Hyphosphere (Ctrl)') 

vk_p6 <- ggplot(ft_df1_z_avg, aes(x=`O/C`, y=`H/C`, color=`hps_cd_avg_z`)) +
  geom_point(alpha=0.6) +
  scale_color_gradient(low = "blue", high = "red") +
  # scale_size(range = c(1, 20)) +# Adjust the range based on your preference
  xlim(-0.1, 1.2)+
  ylim(0, 2.5)+
  geom_rect(aes(xmin = 0, xmax = 0.2, ymin = 1.5, ymax = 2.3), alpha = 0, color = "black", linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.2, xmax = 0.52, ymin = 1.5, ymax = 2.2), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.52, xmax = 0.7, ymin = 1.5, ymax = 2.3), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.7, xmax = 1.1, ymin = 1.5, ymax = 2.4), alpha = 0, color = "black",  linetype='dotted',  linewidth = 0.2)+
  geom_rect(aes(xmin = 0, xmax = 0.25, ymin = 0.5, ymax = 1.25), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.25, xmax = 0.67, ymin = 0.75, ymax = 1.5), alpha = 0, color = "black", linetype='dotted',  linewidth = 0.2)+
  geom_rect(aes(xmin = 0.67, xmax = 0.97, ymin = 0.53, ymax = 1.5), alpha = 0, color = "black", linetype='dotted',  linewidth = 0.2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle('Hyphosphere (Cd)') 

vk_p7 <- ggplot(ft_df1_z_avg, aes(x=`O/C`, y=`H/C`, color=`hps_f_avg_z`)) +
  geom_point(alpha=0.6) +
  scale_color_gradient(low = "blue", high = "red") +
  # scale_size(range = c(1, 20)) +# Adjust the range based on your preference
  xlim(-0.1, 1.2)+
  ylim(0, 2.5)+
  geom_rect(aes(xmin = 0, xmax = 0.2, ymin = 1.5, ymax = 2.3), alpha = 0, color = "black", linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.2, xmax = 0.52, ymin = 1.5, ymax = 2.2), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.52, xmax = 0.7, ymin = 1.5, ymax = 2.3), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.7, xmax = 1.1, ymin = 1.5, ymax = 2.4), alpha = 0, color = "black",  linetype='dotted',  linewidth = 0.2)+
  geom_rect(aes(xmin = 0, xmax = 0.25, ymin = 0.5, ymax = 1.25), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.25, xmax = 0.67, ymin = 0.75, ymax = 1.5), alpha = 0, color = "black", linetype='dotted',  linewidth = 0.2)+
  geom_rect(aes(xmin = 0.67, xmax = 0.97, ymin = 0.53, ymax = 1.5), alpha = 0, color = "black", linetype='dotted',  linewidth = 0.2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle('Hyphosphere (M)') 

vk_p8 <- ggplot(ft_df1_z_avg, aes(x=`O/C`, y=`H/C`, color=`hps_f_cd_avg_z`)) +
  geom_point(alpha=0.6) +
  scale_color_gradient(low = "blue", high = "red") +
  # scale_size(range = c(1, 20)) +# Adjust the range based on your preference
  xlim(-0.1, 1.2)+
  ylim(0, 2.5)+
  geom_rect(aes(xmin = 0, xmax = 0.2, ymin = 1.5, ymax = 2.3), alpha = 0, color = "black", linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.2, xmax = 0.52, ymin = 1.5, ymax = 2.2), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.52, xmax = 0.7, ymin = 1.5, ymax = 2.3), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.7, xmax = 1.1, ymin = 1.5, ymax = 2.4), alpha = 0, color = "black",  linetype='dotted',  linewidth = 0.2)+
  geom_rect(aes(xmin = 0, xmax = 0.25, ymin = 0.5, ymax = 1.25), alpha = 0, color = "black",  linetype='dotted', linewidth = 0.2)+
  geom_rect(aes(xmin = 0.25, xmax = 0.67, ymin = 0.75, ymax = 1.5), alpha = 0, color = "black", linetype='dotted',  linewidth = 0.2)+
  geom_rect(aes(xmin = 0.67, xmax = 0.97, ymin = 0.53, ymax = 1.5), alpha = 0, color = "black", linetype='dotted',  linewidth = 0.2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle('Hyphosphere (M+Cd)') 

# ~~~~ Additional Figure. V-K figure by rhiz and hypho, by group ----
p_vk_g <- ggarrange(vk_p1, vk_p2, vk_p3, vk_p4,
          vk_p5, vk_p6, vk_p7, vk_p8,
          labels = c('a', 'b', 'c', 'd',
                     'e', 'f', 'g', 'h'),
          ncol = 4, nrow = 2, legend = "none")

if(save){
  ggsave('out/vk_plot_by_cmpt_trt.pdf', plot = p_vk_g,
       w = 88*5, h = 66*3, units = 'mm')
ggsave('out/vk_plot_by_cmpt_trt.jpg', plot = p_vk_g,
       w = 88*5, h = 66*3, units = 'mm')
}
rm(vk_p1, vk_p2, vk_p3, vk_p4,
   vk_p5, vk_p6, vk_p7, vk_p8)


# O/C lines: 0, 0.2, 0.25, 0.52, 0.67, 0.7, 0.97, 1.1
# H/C lines: 0.5, 0.53, 0.75, 1.25, 1.5, 2.2, 2.3, 2.4, 

# Lipid-like :
# ---- 0 ≤ O/C ≤ 0.2 & 1.5 ≤ H/C ≤ 2.3 & N/C ≤ 0.04 & P/C ≤ 0.03;
# Protein-like: 
# ---- 0.2 ≤ O/C ≤ 0.52 & 1.5 ≤ H/C ≤ 2.2 & 0.178 ≤ N/C ≤ 0.44 & P/C ≤ 0.06;
# Amino sugar-like: 
# ---- 0.52 ≤ O/C ≤ 0.7 & 1.5 ≤ H/C ≤ 2.2 & 0.07 ＜ N/C ≤ 0.182 & P/C ≤ 0.167;
# Carbohydrate-like:
# ---- 0.7 ≤ O/C ≤ 1.1 & 1.5 ≤ H/C ≤ 2.4 & N = 0 & P = 0 ;
# Condensed aromatics-like:
# ---- 0 ≤ O/C ≤ 0.25 & 0.5 ≤ H/C ≤ 1.25;
# Lignin-like:
# ---- 0.25 ≤ O/C ≤ 0.67 & 0.75 ≤ H/C ≤1.5;
# Tannin-like:
# ---- 0.67 ≤ O/C ≤ 0.97 & 0.53 ≤ H/C ≤ 1.5;
# CRAM / 富羧酸脂环：
# ---- 0.30 ≤ DBE/C ≤ 0.68 & 0.20 ≤ DBE/H ≤ 0.95 & 0.77 ≤ DBE/O ≤ 1.75

  
  

# add formula info to metabo_df and send to Zhang Peng
com_id_formula <- read_xlsx('data/metabo_data/com_id_and_formula.xlsx')
metabo_df_Zhang_Peng <- merge(metabo_df, com_id_formula, by = 'com_id', all = T)
write_xlsx(metabo_df_Zhang_Peng, 'data/metabo_data/metabo_df_Zhang_Peng.xlsx')

# compound and zotu network construction -----------------

# Range of samples in the metabo_df4 within the function metabo_otu_cor_p_tab(). Not the range in the metabo_df file. 
# rs_ck  1-8
# rs_f     9-16
# rs_f_cd  17-24
# rs_cd    25-32
# 
# hps_ck  33-40
# hps_f     41-48
# hps_f_cd  49-56
# hps_cd    57-64



# find differently expressed compounds in POS and NEG modes

# load DEC of POS mode (comparing hps M and M+Cd)
# DEC_hps_M_vs_M_cd_POS <- read_xlsx('data/HB-vs-HC.filter.name.annot_POS.xlsx')
# DEC_hps_M_vs_M_cd_NEG <- read_xlsx('data/HB-vs-HC.filter.name.annot_NEG.xlsx')
# 
# # selection criteria: OPLS-DA VIP≥1 and T-test P<0.05
# dec_list_POS <- DEC_hps_M_vs_M_cd_POS$index
# dec_list_NEG <- DEC_hps_M_vs_M_cd_NEG$index
# dec_list <- c(dec_list_NEG, dec_list_POS)

###################
# hps_ck  33-40
# hps_f     41-48
# hps_f_cd  49-56
# hps_cd    57-64
###################

# Hyphosphere - Ctrl ----
# subsetting data
hps_Ctrl <- as.data.frame(otu_table(subset_samples(pl, Compartment == 'Hyphosphere' & Group == 'Ctrl')))

# construct melted spearman cor p and r table -- also contains the edge file. Numbers are sample range, p value, and r value selected. p and r values can be changed according to requirement/preference
zotu_COR_com_hps_Ctrl_edge <- metabo_otu_cor_p_tab(hps_Ctrl, metabo_df, 33, 40, 0.01, 0.5)

# construct the node file by mering three matrices: 1. edge file, TAX phyloseq object, and metrabo_df.
zotu_COR_com_hps_Ctrl_node <- construct_node_file(zotu_COR_com_hps_Ctrl_edge, TAX, metabo_df)

# save as .csv file for later processing in excel to be used in gephi. i.e., create Source and Target numbers for the edge file. And the edge file needs to be add columns like inter (pp or np - positive or negative) and Label (1 or -1) like those generated by MENA pipeline. 
write.csv(zotu_COR_com_hps_Ctrl_node,"data/metabo_data/zotu_COR_com_hps_Ctrl_node.csv", row.names = FALSE)

write.csv(zotu_COR_com_hps_Ctrl_edge,"data/metabo_data/zotu_COR_com_hps_Ctrl_edge.csv", row.names = FALSE)




# Hyphosphere - M ----
hps_M <- as.data.frame(otu_table(subset_samples(pl, Compartment == 'Hyphosphere' & Group == 'M')))
zotu_COR_com_hps_M_edge <- metabo_otu_cor_p_tab(hps_M, metabo_df, 41, 48, 0.01, 0.5)

zotu_COR_com_hps_M_node <- construct_node_file(zotu_COR_com_hps_M_edge, TAX, metabo_df)


write.csv(zotu_COR_com_hps_M_node,"data/metabo_data/zotu_COR_com_hps_M_node.csv", row.names = FALSE)

write.csv(zotu_COR_com_hps_M_edge,"data/metabo_data/zotu_COR_com_hps_M_edge.csv", row.names = FALSE)



# Hyphosphere - M+Cd ----
hps_M_Cd <- as.data.frame(otu_table(subset_samples(pl, Compartment == 'Hyphosphere' & Group == 'M+Cd'))) # subset target compartment and treatment group
zotu_COR_com_hps_M_Cd_edge <- metabo_otu_cor_p_tab(hps_M_Cd, metabo_df, 59, 66, 0.01, 0.5) # use the metabo_otu_cor_p_tab() function to construct edge file, 
# ***** CHECK if number of sample columns are correct, as the some columns were added to metabo_df so the index of the sample columns were changed. Use colnames(metabo_df) to check sample index. 
# *** NOTE the correct index number is important as incorrect index leads to error. 
# e.g. if type "59, 67", which include 9 samples instead of 8 --> error ***

zotu_COR_com_hps_M_Cd_node <- construct_node_file(zotu_COR_com_hps_M_Cd_edge, TAX, metabo_df)


write.csv(zotu_COR_com_hps_M_Cd_node,"data/metabo_data/zotu_COR_com_hps_M_Cd_node.csv", row.names = FALSE)

write.csv(zotu_COR_com_hps_M_Cd_edge,"data/metabo_data/zotu_COR_com_hps_M_Cd_edge.csv", row.names = FALSE)




# Hyphosphere - Cd ----
hps_Cd <- as.data.frame(otu_table(subset_samples(pl, Compartment == 'Hyphosphere' & Group == 'Cd')))
zotu_COR_com_hps_Cd_edge <- metabo_otu_cor_p_tab(hps_Cd, metabo_df, 57, 64, 0.01, 0.5)


zotu_COR_com_hps_Cd_node <- construct_node_file(zotu_COR_com_hps_Cd_edge, TAX, metabo_df)


write.csv(zotu_COR_com_hps_Cd_node,"data/metabo_data/zotu_COR_com_hps_Cd_node.csv", row.names = FALSE)

write.csv(zotu_COR_com_hps_Cd_edge,"data/metabo_data/zotu_COR_com_hps_Cd_edge.csv", row.names = FALSE)







# rzs_Ctrl  1-8
# rzs_M     9-16
# rzs_M_Cd  17-24
# rzs_Cd    25-32

###############################################################################

###############################################################################
# Rhizosphere - Ctrl ---------------
# subset target subgroups
rzs_Ctrl <- as.data.frame(otu_table(subset_samples(pl, Compartment == 'Rhizosphere' & Group == 'Ctrl')))
# the function: 1 is the first position of the group, 8 is the last position. 0.01 is a selected p value, 0.8 is a selected R value (spearman correlation). Other values can be selected, such as p = 0.05 and R = 0.6
zotu_COR_com_rzs_Ctrl_edge <- metabo_otu_cor_p_tab(rzs_Ctrl, metabo_df, 1, 8, 0.01, 0.5)

zotu_COR_com_rzs_Ctrl_node <- construct_node_file(zotu_COR_com_rzs_Ctrl_edge, TAX, metabo_df)

write.csv(zotu_COR_com_rzs_Ctrl_node,"data/metabo_data/zotu_COR_com_rzs_Ctrl_node.csv", row.names = FALSE)

write.csv(zotu_COR_com_rzs_Ctrl_edge,"data/metabo_data/zotu_COR_com_rzs_Ctrl_edge.csv", row.names = FALSE)



# Rhizosphere - M ----
rzs_M <- as.data.frame(otu_table(subset_samples(pl, Compartment == 'Rhizosphere' & Group == 'M')))
zotu_COR_com_rzs_M_edge <- metabo_otu_cor_p_tab(rzs_M, metabo_df, 9, 16, 0.01, 0.5)

zotu_COR_com_rzs_M_node <- construct_node_file(zotu_COR_com_rzs_M_edge, TAX, metabo_df)

write.csv(zotu_COR_com_rzs_M_node,"data/metabo_data/zotu_COR_com_rzs_M_node.csv", row.names = FALSE)

write.csv(zotu_COR_com_rzs_M_edge,"data/metabo_data/zotu_COR_com_rzs_M_edge.csv", row.names = FALSE)

# Rhizosphere - M+Cd ----
rzs_M_Cd <- as.data.frame(otu_table(subset_samples(pl, Compartment == 'Rhizosphere' & Group == 'M+Cd')))
zotu_COR_com_rzs_M_Cd_edge <- metabo_otu_cor_p_tab(rzs_M_Cd, metabo_df, 17, 24, 0.01, 0.5)

zotu_COR_com_rzs_M_Cd_node <- construct_node_file(zotu_COR_com_rzs_M_Cd_edge, TAX, metabo_df)

write.csv(zotu_COR_com_rzs_M_Cd_node,"data/metabo_data/zotu_COR_com_rzs_M_Cd_node.csv", row.names = FALSE)

write.csv(zotu_COR_com_rzs_M_Cd_edge,"data/metabo_data/zotu_COR_com_rzs_M_Cd_edge.csv", row.names = FALSE)


# Rhizosphere - Cd ----
rzs_Cd <- as.data.frame(otu_table(subset_samples(pl, Compartment == 'Rhizosphere' & Group == 'Cd')))
zotu_COR_com_rzs_Cd_edge <- metabo_otu_cor_p_tab(rzs_Cd, metabo_df, 25, 32, 0.01, 0.5)

zotu_COR_com_rzs_Cd_node <- construct_node_file(zotu_COR_com_rzs_Cd_edge, TAX, metabo_df)

write.csv(zotu_COR_com_rzs_Cd_node,"data/metabo_data/zotu_COR_com_rzs_Cd_node.csv", row.names = FALSE)

write.csv(zotu_COR_com_rzs_Cd_edge,"data/metabo_data/zotu_COR_com_rzs_Cd_edge.csv", row.names = FALSE)






# chemical diversity and ASV diversity calculation ----
# ~~ using the metabo_df and Shannon index ----
# H' = - sum of 1-s pi*lnPi


library(vegan) # using diversity() function in vegan package to calculate Shannon index (H') of metabo and ASV
tdf <- metabo_df %>% select(11:74)# select numeric columns
tdf_tp <- t(tdf) #transpose
chem_div <- vegan::diversity(tdf_tp, 'shannon') # cal Shannon index

# for H/C > 1.5 compounds: 
metabo_df_HC1.5 <- filter (metabo_df, HC_ratio < 1.5)
chem_div_HC1.5 <- vegan::diversity(t(metabo_df_HC1.5 %>% select(11:74)), 'shannon')

# cal ASV Shannon index
ZOTU_rzs_hps <- as.data.frame(otu_table(pl_un_rarefied)) %>% select(1:32, 129:(129+31))
ASV_table_tp <- t(otu_table(pl_un_rarefied))
ASV_div <- vegan::diversity(ASV_table_tp, 'shannon')

# combine calculated results as a df
# convert named vectors (chem_div and ASV_div) to dataframe, and then transpose

chem_div_df  <-  as.data.frame(t(data.frame(as.list(chem_div)))) %>%  rownames_to_column(var="sample_ID")


chem_div_df_HC1.5  <-  as.data.frame(t(data.frame(as.list(chem_div_HC1.5)))) %>%  rownames_to_column(var="sample_ID")

ASV_div_df <-  as.data.frame(t(data.frame(as.list(ASV_div)))) %>%  rownames_to_column(var="sample_ID")

chem_ASV_div <- merge(chem_div_df, ASV_div_df, by= 'sample_ID')

chem_ASV_div_HC1.5 <- merge(chem_div_df_HC1.5, ASV_div_df, by= 'sample_ID')

# correlation
# prepare data using excel
library("writexl")
if(save){
  write_xlsx(chem_ASV_div,"out/chem_ASV_div.xlsx")
  write_xlsx(chem_ASV_div_HC1.5,"out/chem_ASV_div_HC1.5.xlsx")
}
# the following chem_ASV_div_add_trt_compartment_info.xlsx file was modified from chem_ASV_div.xlsx by adding columns: trt and compartment, and change other column names

chem_ASV_div2 <- read_xlsx('out/chem_ASV_div_add_trt_compartment_info.xlsx')

chem_ASV_div2$trt <- factor(chem_ASV_div2$trt, 
                           levels=c('Ctrl', 'Cd', 'M', 'M+Cd'))


chem_ASV_div2_HC1.5 <- read_xlsx('out/chem_ASV_div_HC1.5_add_trt_compartment_info.xlsx')

chem_ASV_div2_HC1.5$trt <- factor(chem_ASV_div2_HC1.5$trt, 
                            levels=c('Ctrl', 'Cd', 'M', 'M+Cd'))


# compare chemodiv ----
chem_ASV_div2_HC1.5 %>% subset(compartment == 'Rhizosphere') %>% 
  plot_box(x = trt, y = chem_div)


p1 <- chem_ASV_div2%>% subset(compartment == 'Rhizosphere') %>% 
  plot_box(x = trt, y = chem_div)+
  geom_jitter(alpha = 0.3, size = 1, width = 0.1)+
  ggtitle('Rhizosphere')+
  xlab('Treatment')+
  ylab('Chemodiversity')+
  ylim(3.6, 4.8)+
  annotate("text", x = 1, y = 4.8, label = "a") +
  annotate("text", x = 2, y = 4.8, label = "a") +
  annotate("text", x = 3, y = 4.8, label = "a") +
  annotate("text", x = 4, y = 4.8, label = "a")

# tukey 
anova_rzs_chemdiv <- aov(chem_div~trt, 
                         data = chem_ASV_div2%>% 
                           subset(compartment == 'Rhizosphere'))

summary(anova_rzs_chemdiv)
HSD.test(anova_rzs_chemdiv, 'trt', console=TRUE)
# chem_div groups
# M    4.322831      a
# M+Cd 4.292857      a
# Cd   4.243533      a
# Ctrl 4.204204      a


chem_ASV_div2_HC1.5 %>% subset(compartment == 'Hyphosphere') %>% 
  plot_box(x = trt, y = chem_div)

p2 <- chem_ASV_div2 %>% subset(compartment == 'Hyphosphere') %>% 
  plot_box(x = trt, y = chem_div)+
  geom_jitter(alpha = 0.3, size = 1, width = 0.1)+
  ggtitle('Hyphosphere')+
  xlab('Treatment')+
  ylab('Chemodiversity')+
  ylim(3.6, 4.8)+
  annotate("text", x = 1, y = 4.8, label = "b") +
  annotate("text", x = 2, y = 4.8, label = "b") +
  annotate("text", x = 3, y = 4.8, label = "a") +
  annotate("text", x = 4, y = 4.8, label = "ab")

# tukey 
anova_hps_chemdiv <- aov(chem_div~trt, 
                         data = chem_ASV_div2%>% 
                           subset(compartment == 'Hyphosphere'))

summary(anova_hps_chemdiv)
HSD.test(anova_hps_chemdiv, 'trt', console=TRUE)

# chem_div groups
# M    4.282087      a
# M+Cd 4.165352     ab
# Ctrl 3.984108      b
# Cd   3.951077      b

# ~~~~ Fig. S5. Chemodiversity in the rhizosphere and hyphosphere ----
p_chemodiv <- ggarrange(p1, p2, 
          labels = c('a', 'b'),
          nrow = 1, ncol = 2)
if(save){
  ggsave('out/Fig. S5. chemodiv_compare.pdf', plot = p_chemodiv,
       w = 88, h = 66, units = 'mm')

ggsave('out/Fig. S5. chemodiv_compare.jpg', plot = p_chemodiv,
       w = 88, h = 66, units = 'mm')
}

# all chem div COR asv div ----

plot1 <- chem_ASV_div2_HC1.5 %>% 
  ggplot(aes(chem_div, asv_div))+
  geom_point(aes(color = compartment), alpha = 0.4, size = 3)+
  geom_smooth(aes(color = compartment), method = lm, se = T)+
  ylab('Bacterial diversity')+
  xlab('Chemodiversity')+
  theme(legend.position = c(0.14, 0.18))+
  ggtitle('Overall')




# rzs - all chem div COR asv div ----
plot2 <- chem_ASV_div2_HC1.5 %>% 
  subset(compartment == 'Rhizosphere') %>% 
  ggplot(aes(chem_div, asv_div))+
  geom_point(aes(color = trt), alpha = 0.4, size = 3)+
  geom_smooth(aes(color = trt), method = lm, se = F)+
  ylab('Bacterial diversity')+
  xlab('Chemodiversity')+
  theme(legend.position = c(0.14, 0.25))+
  ggtitle('Rhizosphere')

# rzs - different treatments -- chem_div COR asv_div ----

chem_ASV_div2_HC1.5 %>% 
  subset(compartment == 'Rhizosphere'& trt == 'Ctrl') %>% 
  ggplot(aes(asv_div, chem_div))+
  geom_point()+
  geom_smooth(method = lm, se = T)

chem_ASV_div2_HC1.5 %>% 
  subset(compartment == 'Rhizosphere'& trt == 'Cd') %>% 
  ggplot(aes(asv_div, chem_div))+
  geom_point()+
  geom_smooth(method = lm, se = T)

chem_ASV_div2_HC1.5 %>% 
  subset(compartment == 'Rhizosphere'& trt == 'M') %>% 
  ggplot(aes(asv_div, chem_div))+
  geom_point()+
  geom_smooth(method = lm, se = T)

chem_ASV_div2_HC1.5 %>% 
  subset(compartment == 'Rhizosphere'& trt == 'M+Cd') %>% 
  ggplot(aes(asv_div, chem_div))+
  geom_point()+
  geom_smooth(method = lm, se = T)

# hps - all chem div COR asv div ----

plot3 <- chem_ASV_div2_HC1.5 %>% 
  subset(compartment == 'Hyphosphere') %>% 
  ggplot(aes(chem_div, asv_div))+
  geom_point(aes(color = trt), alpha = 0.4, size = 3)+
  geom_smooth(aes(color = trt), method = lm, se = F)+
  ylab('Bacterial diversity')+
  xlab('Chemodiversity')+
  theme(legend.position = c(0.14, 0.25))+
  ggtitle('Hyphosphere')

# hps - different treatments -- chem_div COR asv_div ----
chem_ASV_div2_HC1.5 %>% 
  subset(compartment == 'Hyphosphere' & trt == 'Ctrl') %>% 
  ggplot(aes(asv_div, chem_div))+
  geom_point()+
  geom_smooth(method = lm, se = T)

chem_ASV_div2_HC1.5 %>% 
  subset(compartment == 'Hyphosphere' & trt == 'Cd') %>% 
  ggplot(aes(asv_div, chem_div))+
  geom_point()+
  geom_smooth(method = lm, se = T)


chem_ASV_div2_HC1.5 %>% 
  subset(compartment == 'Hyphosphere' & trt == 'M') %>% 
  ggplot(aes(asv_div, chem_div))+
  geom_point()+
  geom_smooth(method = lm, se = T)

chem_ASV_div2_HC1.5 %>% 
  subset(compartment == 'Hyphosphere' & trt == 'M+Cd') %>% 
  ggplot(aes(asv_div, chem_div))+
  geom_point()+
  geom_smooth(method = lm, se = T)




# combine plots ----

# ~~~~ Fig. S6. chemodiversity cor alpha diversity ----
p_chemo_cor_alpha_div <- ggarrange(plot1, plot2, plot3, 
          labels = c('a', 'b', 'c'),
          nrow = 2, ncol = 2)

if(save){
  ggsave('out/Fig. S6. chemodiv_COR_bactDiv_HC1.5.pdf', plot = p_chemo_cor_alpha_div,
       width = 88*3, height = 66*3, units = "mm")
  ggsave('out/Fig. S6. chemodiv_COR_bactDiv_HC1.5.jpg', plot = p_chemo_cor_alpha_div,
       width = 88*3, height = 66*3, units = "mm")
}
rm(plot1, plot2, plot3)



# get cor statistics ----

temp_df <- chem_ASV_div2_HC1.5 %>% 
  subset(compartment == 'Hyphosphere' & trt == 'M+Cd' ) 
cor.test(temp_df$chem_div, temp_df$asv_div, method = "spearman")


# overall: hyphosphere, p-value = 0.001271**, rho = -0.5516862 
#     Ctrl: p-value = 0.5821, rho =  -0.2380952
#     Cd: p-value = 0.4279, rho = 0.3333333 
#     M: p-value = 0.9349, rho = 0.04761905 
#     M+Cd: p-value = 0.2431, rho = -0.4761905



# overall: rhizosphere, p-value = 0.0402*, rho = 0.3658358
#     Ctrl: p-value = 0.1966, rho = 0.5238095 
#     Cd: p-value = 1, rho = 0
#     M: p-value = 0.2162, rho = 0.5
#     M+Cd: p-value = 0.2162, rho = 0.5



# chemical relative quantity comparison in the hyphosphere ----



temp_data <- melt(metabo_df %>% select(com_id, com_name, 11:74), id=c("com_id","com_name")) 
write.csv(temp_data, 'data/metabo_data/metabo_df_for_quantity.csv')

# after add compartment and trt 
metabo_df_melt_comp_trt <- read_csv('data/metabo_data/metabo_df_for_quantity_add_comp_trt.csv')

metabo_df_melt_comp_trt$trt <- factor(metabo_df_melt_comp_trt$trt, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))

metabo_df_melt_comp_trt %>% subset(compartment == 'Rhizosphere') %>% 
  plot_box(x = trt, y = log(value))

metabo_df_melt_comp_trt %>% subset(compartment == 'Hyphosphere') %>% 
  plot_box(x = trt, y = log(value))

# comparing the metabolome in Cd vs M+Cd

# plot_box(soil_data, group2, TOC)
# plot_box(soil_data, group2, DOC)
# plot_box(soil_data, group2, moi)
# plot_box(soil_data, group2, `EE-GRSP`)
# plot_box(soil_data, group2, `T-GRSP`)
# plot_box(soil_data, group2, TCd)



# hyphosphere M -- what bacteria are connected to what metabolome? -----
# the following tables were prepared using excel, based on the node file and edge file
hps_M_mother <- read_csv('data/metabo_data/hps_M_metabo_COR_zotu_mother_tab.csv')

hps_M_zotu <- read_xlsx('data/metabo_data/hps_M_zotu_and_taxon.xlsx')

hps_M_metabo <- read_xlsx('data/metabo_data/hps_M_com_IDs_and_com_names.xlsx')

temp <- merge(hps_M_mother, hps_M_zotu, by = 'Zotu', all = T)

hps_M_all_join <- merge(temp, hps_M_metabo, by = 'com_id', all = T)

write.csv(hps_M_all_join, 'data/metabo_data/hps_M_zotu_COR_metabo_all_join.csv')
 
# hyphosphere M+Cd  -- what bacteria are connected to what metabolome? -----
hps_M_Cd_mother <- read_csv('data/metabo_data/hps_M_Cd_metabo_COR_zotu_mother_tab.csv')

hps_M_Cd_zotu <- read_xlsx('data/metabo_data/hps_M_Cd_zotu_and_taxon.xlsx')

hps_M_Cd_metabo <- read_xlsx('data/metabo_data/hps_M_Cd_com_IDs_and_com_names.xlsx')

temp_Cd <- merge(hps_M_Cd_mother, hps_M_Cd_zotu, by = 'Zotu', all = T)

hps_M_Cd_all_join <- merge(temp_Cd, hps_M_Cd_metabo, by = 'com_id', all = T)
write.csv(hps_M_Cd_all_join, 'data/metabo_data/hps_M_Cd_zotu_COR_metabo_all_join.csv')

