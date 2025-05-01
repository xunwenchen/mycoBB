# Heatmap plots differently changed bacterial abundance -----
# subset target groups
# ~~~ M+Cd relative to M (Cd effect) ----
# pl_subset_rzp <- subset_samples(pl_un_rarefied, Compartment == "Rhizoplane" & Mycorrhizal == "M")

pl_subset_rzs_M <- subset_samples(pl_un_rarefied, Compartment == "Rhizosphere" & Mycorrhizal == "M")

pl_subset_hps_M <- subset_samples(pl_un_rarefied, Compartment == "Hyphosphere" & Mycorrhizal == "M")

pl_subset_rzs_Cd <- subset_samples(pl_un_rarefied, Compartment == "Rhizosphere" & Cd == 5)

pl_subset_hps_Cd <- subset_samples(pl_un_rarefied, Compartment == "Hyphosphere" & Cd == 5)



hm1 <- pheatmap_vst_scaled(pl_subset_rzs_M, 
                           trt_group = 'M+Cd', 
                           ref_group = 'M', 
                           low_count_threshold = 20, 
                           top_x = 50)

hm2 <- pheatmap_vst_scaled(pl_subset_hps_M, 
                           trt_group = 'M+Cd', 
                           ref_group = 'M', 
                           low_count_threshold = 20, 
                           top_x = 50)

hm3 <- pheatmap_vst_scaled(pl_subset_rzs_Cd, 
                           trt_group = 'Cd', 
                           ref_group = 'M+Cd', 
                           low_count_threshold = 20, 
                           top_x = 50)

hm4 <- pheatmap_vst_scaled(pl_subset_hps_Cd, 
                           trt_group = 'Cd', 
                           ref_group = 'M+Cd', 
                           low_count_threshold = 20, 
                           top_x = 50)



# Arrange the plots
hm1_grob <- hm1[[4]]
hm2_grob <- hm2[[4]]
hm3_grob <- hm3[[4]]
hm4_grob <- hm4[[4]]

hm_comb <- ggarrange(hm1_grob, hm2_grob, hm3_grob, hm4_grob,
                     labels = c('a', 'b', 'c', 'd'),
                     nrow = 4, ncol = 1)

# ~~~~~~~~ Extended Data Fig. 5. Compartment-wise bacterial abundance changes and KEGG enrichment analysis according to metabolomic profiles (KEGG is from metabo data, not shown in current codes) ----

hm_comb

if(save){
ggsave('out/ExtDataFig. 5. Compt-wise bact ab changes and KEGG.jpg', plot = hm_comb,
       width = 88*2, h = 66*4, units = 'mm', dpi = 600)
ggsave('out/ExtDataFig. 5. Compt-wise bact ab changes and KEGG.pdf', plot = hm_comb,
       width = 88*2, h = 66*4, units = 'mm', dpi = 600)
}


# ~~~~ fungal effects on clean soil ----
pl_subset_hps_clean <- subset_samples(pl_un_rarefied, Compartment == "Hyphosphere" & Cd == "0")

hm5 <- pheatmap_vst_scaled(pl_subset_hps_clean, 
                           trt_group = 'M', 
                           ref_group = 'Ctrl', 
                           low_count_threshold = 10, 
                           top_x = 50)

hm5_un_scaled <- pheatmap_vst(pl_subset_hps_clean, 
                              trt_group = 'M', 
                              ref_group = 'Ctrl', 
                              low_count_threshold = 10, 
                              top_x = 50)



# assign tax to zotu ids of the heatmap manually ----

# check tax of a Zotu
tax_table(pl)['Zotu518', ]

top50zotu <- read_excel('data/Zotu_tax_cd_fungal_effect.xlsx')
# 'data/Zotu_tax_cd_fungal_effect.xlsx' was obtained by capture zotu IDs in the hm_comb manually. It also can be done by DESeq2 and sort out the top tax according to adjusted p-value. Check codes in function.R - function pheatmap_vst()

# take different groups
top50zotu_rzs_cd_eff <- top50zotu$rzs_cd_effect
top50zotu_hps_cd_eff <- top50zotu$hps_cd_effect
top50zotu_rzs_fungal_eff <- top50zotu$rzs_fungal_effect
top50zotu_hps_fungal_eff <- top50zotu$hps_fungal_effect


# zotu and tax table
top50tax_rzs_cd_eff <- (tax_table(pl))[top50zotu_rzs_cd_eff, ]
top50tax_hps_cd_eff <- (tax_table(pl))[top50zotu_hps_cd_eff, ]
top50tax_rzs_fungal_eff <- (tax_table(pl))[top50zotu_rzs_fungal_eff, ]
top50tax_hps_fungal_eff <- (tax_table(pl))[top50zotu_hps_fungal_eff, ]

# the following excel files were further modified in excel manually to identify to the Family level. Some tax not identified to the Family level were noted with in parentheses, e.g., (o) = (Order level). Some tax can reach Genus or Species levels but too many NA, so I go for the Family level as the furthest. 
write_xlsx(as.data.frame(top50tax_rzs_cd_eff), 'out/top50tax_rzs_cd_eff0.xlsx')
write_xlsx(as.data.frame(top50tax_hps_cd_eff), 'out/top50tax_hps_cd_eff0.xlsx')
write_xlsx(as.data.frame(top50tax_rzs_fungal_eff), 'out/top50tax_rzs_fungal_eff0.xlsx')
write_xlsx(as.data.frame(top50tax_hps_fungal_eff), 'out/top50tax_hps_fungal_eff0.xlsx')

# the above excel files were modified and new excel files were named xxx_eff.xlsx. The '0' was deleted to differentiate with original files. The tax info were integrated in heatmap plots, replacing original Zotu IDs. Note the order of Zotu IDs and tax rank should match. Check if the Zotu IDs in head(top50tax_rzs_cd_eff) (look from top to bottom) match in the heatmap hm1 (look from left to right). Check (note the plots may overlap so delete former one before run next one):
head(top50tax_rzs_cd_eff, 10); hm1
head(top50tax_hps_cd_eff, 10); hm2
head(top50tax_rzs_fungal_eff, 10); hm3
head(top50tax_hps_fungal_eff, 10); hm4


# Log2 FoldChange plots ----
# To change tax level, go to 'fun.R', find 'log2fc' function, modify, and source the fun.R again: source('code/fun.R') 

# ~~Cd effects on mycorrhizal soils - compartment-wise ----
# rhizoplane
pl_unra_rzp<- subset_samples(pl_un_rarefied, Compartment == 'Rhizoplane') 
# If DESeq2 package is involved, raw count (z)otu table (i.e., pl_un_rarefied) should be used instead of normalized/rarefied z(otu) table (i.e., pl). 'unra' means 'un-rarefied'.
pl_unra_rzp_AM <- subset_samples(pl_unra_rzp, Mycorrhizal == 'M')
pp1 <- log2fc(pl_unra_rzp_AM, 'Group', 'M+Cd', 'M')+ 
  # note: 'Group','trt_group','ref_group'. This case is M+Cd over M
  ggtitle('Rhizoplane (M+Cd relative to M)')+
  coord_flip()+
  scale_x_discrete(limits=rev)

# rhizosphere
pl_unra_rzs<- subset_samples(pl_un_rarefied, Compartment == 'Rhizosphere')
pl_unra_rzs_AM <- subset_samples(pl_unra_rzs, Mycorrhizal == 'M')
pp2 <- log2fc(pl_unra_rzs_AM, 'Group', 'M+Cd', 'M')+
  ggtitle('Rhizosphere (M+Cd relative to M)')+
  coord_flip()+
  scale_x_discrete(limits=rev)

# hyphosphere
pl_unra_hps<- subset_samples(pl_un_rarefied, Compartment == 'Hyphosphere')
pl_unra_hps_AM <- subset_samples(pl_unra_hps, Mycorrhizal == 'M')
head(sample_data(pl_unra_hps_AM)$Group, 32) # only two levels: M and M+Cd
pp3 <- log2fc(pl_unra_hps_AM, Group = 'Group', 
              trt_group = 'M+Cd', ref_group = 'M')+
  ggtitle('Hyphosphere (M+Cd relative to M)')+
  coord_flip()+
  scale_x_discrete(limits=rev)



# hyphae
pl_unra_hy<- subset_samples(pl_un_rarefied, Compartment == 'Hyphae')
head(sample_data(pl_unra_hy)$Group, 32) # only two levels: M and M+Cd
pp4 <- log2fc(pl_unra_hy, 'Group', 'M+Cd', 'M')+
  ggtitle('Hyphae (M+Cd relative to M)')+
  coord_flip()+
  scale_x_discrete(limits=rev)

# ~~~~~~~~ Fig. S8. Cd_eff_family_log2fc ----
Cd_eff_family_log2fc <- ggarrange(pp1, pp2, pp3, pp4,
                                  labels = c("a", "b", "c", "d"),
                                  ncol = 2, nrow = 2)
Cd_eff_family_log2fc

if(save){
  ggsave('out/Fig. S8. Cd_eff_family_log2fc.jpg', plot = Cd_eff_family_log2fc,
       width = 88*4, height = 66*4, units = "mm")
ggsave('out/Fig. S8. Cd_eff_family_log2fc.pdf', plot = Cd_eff_family_log2fc,
       width = 88*4, height = 66*4, units = "mm")
}

rm(pp1, pp2, pp3, pp4)



# ~~mycorrhizal effects on Cd-polluted soil (Cd vs M+Cd) ----
pl_unra_rzp<- subset_samples(pl_un_rarefied, Compartment == 'Rhizoplane')
pl_unra_rzp_cd <- subset_samples(pl_unra_rzp, Cd == 5)
pp5 <- log2fc(pl_unra_rzp_cd, 'Group', 'M+Cd', 'Cd')+
  ggtitle('Rhizoplane (M+Cd relative to Cd)')+
  coord_flip()+
  scale_x_discrete(limits=rev)


pl_unra_rzs<- subset_samples(pl_un_rarefied, Compartment == 'Rhizosphere')
pl_unra_rzs_cd <- subset_samples(pl_unra_rzs, Cd == 5)
pp6 <- log2fc(pl_unra_rzs_cd, 'Group', 'M+Cd', 'Cd')+
  ggtitle('Rhizosphere (M+Cd relative to Cd)')+
  coord_flip()+
  scale_x_discrete(limits=rev)


pl_unra_hps<- subset_samples(pl_un_rarefied, Compartment == 'Hyphosphere')
pl_unra_hps_cd <- subset_samples(pl_unra_hps, Cd == 5)
head(sample_data(pl_unra_hps_cd)$Group, 32) 
pp7 <- log2fc(pl_unra_hps_cd, Group = 'Group', 
              trt_group = 'M+Cd', ref_group = 'Cd')+
  ggtitle('Hyphosphere (M+Cd relative to Cd)')+
  coord_flip()+
  scale_x_discrete(limits=rev)

# ~~~~~~~~ Fig. S9. myco_eff_family_log2fc ----

myco_eff_family_log2fc <- ggarrange(pp5, pp6, pp7,
                                    labels = c("a", "b", "c"),
                                    ncol = 3, nrow = 1)
myco_eff_family_log2fc

if(save){
ggsave('out/Fig. S9. myco_eff_family_log2fc.jpg', plot = myco_eff_family_log2fc,
       width = 88*7, height = 66*3.5, units = "mm")
ggsave('out/Fig. S9. myco_eff_family_log2fc.pdf', plot = myco_eff_family_log2fc,
       width = 88*7, height = 66*3.5, units = "mm")
}

rm(pp5, pp6, pp7)



# ~~mycorrhizal effects on clean soils ----
# hyphosphere
pl_unra_hps_M_ctrl <- subset_samples(pl_unra_hps, Group == 'M' | Group == 'Ctrl') # subset those without Cd (clean soils)
head(sample_data(pl_unra_hps_M_ctrl)$Group, 32) # only two levels: Ctrl and M

# ~~~~~~~~ Fig. S7. myco_eff_family_log2fc_clean_soil ----

pp_M_vs_ctrl <- log2fc(pl_unra_hps_M_ctrl, Group = 'Group', 
                       trt_group = 'M', ref_group = 'Ctrl')+
  ggtitle('Hyphosphere (M relative to Ctrl)')+
  coord_flip()+
  scale_x_discrete(limits=rev) 

pp_M_vs_ctrl

if(save){
ggsave('out/Fig. S7. myco_eff_family_log2fc_clean_soil.jpg', plot = pp_M_vs_ctrl,
       width = 88*3, height = 66*6.5, units = "mm")
ggsave('out/Fig. S7. myco_eff_family_log2fc_clean_soil.pdf', plot = pp_M_vs_ctrl,
       width = 88*3, height = 66*6.5, units = "mm")
}

# ~~~~~~~~ Fig. S7x. myco_eff_genus_log2fc_clean_soil ----
pp_M_vs_ctrl_hyphosphere_genus <- log2fc_genus(pl_unra_hps_M_ctrl, Group = 'Group', 
                                         trt_group = 'M', ref_group = 'Ctrl')+
  ggtitle('Hyphosphere (M relative to Ctrl)')+
  coord_flip()+
  scale_x_discrete(limits=rev) 
pp_M_vs_ctrl_hyphosphere_genus

if(save){
  ggsave('out/Fig. S7x. myco_eff_genus_log2fc_clean_soil.jpg', plot = pp_M_vs_ctrl_hyphosphere_genus,
       width = 88*3, height = 66*8, units = "mm")
ggsave('out/Fig. S7x. myco_eff_genus_log2fc_clean_soil.pdf', plot = pp_M_vs_ctrl_hyphosphere_genus,
       width = 88*3, height = 66*8, units = "mm")

}
# DONE 
