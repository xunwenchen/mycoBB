# ~ Beta-diversity ----
# ~ ~All compartments ----
# All samples interactive effects of multiple factors
# calculate weighted unifrac distance matrix
w.uni_dist <- UniFrac(physeq = pl, weighted = TRUE, normalized = TRUE)
sample_data <- data.frame(sample_data(pl)) # data.frame() is required

if(is.factor(sample_data$Cd)&&
   is.factor(sample_data$Mycorrhizal)&&
   is.factor(sample_data$Compartment)&&
   is.factor(sample_data$Group)){
  print('Variables are factors, please proceed')
} else {
  print('Not all variables are factors, please check')}

# Permutation ANOVA. Note the order of factors is important when numbers of samples are difference. The most influential factor shall go first. Here I select Compartment. But it is not clear whether Mycorrhizal or Cd should go first.
# ~~~~~~~ three-factor interactive Permanova ----
three_factor_permutation_rs <- adonis2(w.uni_dist ~ Compartment * Mycorrhizal * Cd, data = sample_data, permutations = 999, by = 'terms') # by = "terms" will assess the significance for each term sequentially.


three_factor_permutation_rs


# ~ ~Different compartments - PCoA group by treatment -----------------
# rarefied table 'pl' should be used for PCoA (beta-diversity)
p1 <-  plot_w_unifrac(subset_samples(pl, Compartment == 'Endosphere'), 'Group')+
  ggtitle("Endosphere")+
  stat_ellipse()

p1uw <- plot_uw_unifrac(subset_samples(pl, Compartment == 'Endosphere'), 'Group')+
  ggtitle("Endosphere (U.W. UniFrac)")+
  stat_ellipse()

# adonis
# construct weighted unifrac distance matrix
dist <- phyloseq::distance(subset_samples(pl, Compartment == 'Endosphere'), method ='wunifrac')
# pairwise Adonis
pair_adonis1 <- pairwiseAdonis::pairwise.adonis(dist, sample_data(subset_samples(pl, Compartment == 'Endosphere'))$Group, perm = 999)
pair_adonis1 # adonis results: all n.s.? M vs Cd -> p.adj = 0.024
pair_adonis1$Compartment <- "Endosphere"


# Interactive Adonis using adonis2() to find interactive effects of multiple factors
sample_data_endo <- data.frame(sample_data(subset_samples(pl, Compartment == 'Endosphere'))) # important: data.frame() is needed

library(vegan)
package.version('vegan')# version 2.6-4
inter_adonis_endo <- adonis2(dist ~ Cd * Mycorrhizal, data = sample_data_endo, permutations = 999, by = 'terms')

inter_adonis_endo

rm(dist)

p2 <-  plot_w_unifrac(subset_samples(pl, Compartment == 'Rhizoplane'), 'Group')+
  ggtitle("Rhizoplane")+
  stat_ellipse()

# adonis
dist <- phyloseq::distance(subset_samples(pl, Compartment == 'Rhizoplane'), method ='wunifrac')
pair_adonis2 <- pairwiseAdonis::pairwise.adonis(dist, sample_data(subset_samples(pl, Compartment == 'Rhizoplane'))$Group, perm = 999)
pair_adonis2 # adonis results: <0.01*: Ctrl vs M+Cd, M vs Cd

pair_adonis2$Compartment <- 'Rhizoplane'

# Interactive Adonis using adonis2() to find interactive effects of multiple factors
sample_data_rzp <- data.frame(sample_data(subset_samples(pl, Compartment == 'Rhizoplane'))) # important: data.frame() is needed

inter_adonis_rzp <- adonis2(dist ~ Cd * Mycorrhizal, data = sample_data_rzp, permutations = 999, by = 'terms')

inter_adonis_rzp

rm(dist)


p3 <-  plot_w_unifrac(subset_samples(pl, Compartment == 'Rhizosphere'), 'Group')+
  ggtitle("Rhizosphere")+
  stat_ellipse()

# adonis
dist <- phyloseq::distance(subset_samples(pl, Compartment == 'Rhizosphere'), method ='wunifrac')
pair_adonis3 <- pairwiseAdonis::pairwise.adonis(dist, sample_data(subset_samples(pl, Compartment == 'Rhizosphere'))$Group, perm = 9999)
pair_adonis3
pair_adonis3$Compartment <- 'Rhizosphere'

# Interactive Adonis using adonis2() to find interactive effects of multiple factors
sample_data_rzs <- data.frame(sample_data(subset_samples(pl, Compartment == 'Rhizosphere'))) # important: data.frame() is needed

inter_adonis_rzs <- adonis2(dist ~ Cd * Mycorrhizal, data = sample_data_rzs, permutations = 999, by = 'terms')

inter_adonis_rzs

rm(dist)

p4 <-  plot_w_unifrac(subset_samples(pl, Compartment == 'Hyphosphere'), 'Group')+
  ggtitle("Hyphosphere")+
  stat_ellipse()

p4uw <- plot_uw_unifrac(subset_samples(pl, Compartment == 'Hyphosphere'), 'Group')+
  ggtitle("Hyphosphere (U.W. UniFrac)")+
  stat_ellipse()

p4bc <- plot_PCoA_bray(subset_samples(pl, Compartment == 'Hyphosphere'), 'Group')+
  ggtitle("Hyphosphere (Bray-Curtis)")+
  stat_ellipse()

# adonis
dist <- phyloseq::distance(subset_samples(pl, Compartment == 'Hyphosphere'), method ='wunifrac')
pair_adonis4 <- pairwiseAdonis::pairwise.adonis(dist, sample_data(subset_samples(pl, Compartment == 'Hyphosphere'))$Group, perm = 9999)
pair_adonis4

pair_adonis4$Compartment <- 'Hyphosphere'


# Interactive Adonis using adonis2() to find interactive effects of multiple factors
sample_data_hypho <- data.frame(sample_data(subset_samples(pl, Compartment == 'Hyphosphere'))) # important: data.frame() is needed

inter_adonis_hypho <- adonis2(dist ~ Cd * Mycorrhizal, data = sample_data_hypho, permutations = 999, by = 'terms')

inter_adonis_hypho

rm(dist)


p5 <-  plot_w_unifrac(subset_samples(pl, Compartment == 'Hyphae'), 'Group')+
  scale_shape_manual(values = c(17, 18)) +
  scale_color_manual(values = c("#f4a582", "#ca0020"))+ # this will mask the colors and shapes used in the fun.R
  ggtitle("Hyphae")+
  stat_ellipse()


# adonis
dist <- phyloseq::distance(subset_samples(pl, Compartment == 'Hyphae'), method ='wunifrac')
pair_adonis5 <- pairwiseAdonis::pairwise.adonis(dist, sample_data(subset_samples(pl, Compartment == 'Hyphae'))$Group, perm = 9999)

pair_adonis5$Compartment <- 'Hyphae'

# Interactive Adonis using adonis2() to find interactive effects of multiple factors
sample_data_hyphae <- data.frame(sample_data(subset_samples(pl, Compartment == 'Hyphae'))) # important: data.frame() is needed

inter_adonis_hyphae <- adonis2(dist ~ Cd, data = sample_data_hyphae, permutations = 999, by = 'terms')

inter_adonis_hyphae

rm(dist)

p6 <-  plot_w_unifrac(subset_samples(pl, Compartment == 'Bulk soil'), 'Group')+
  ggtitle("Bulk soil")+
  stat_ellipse()

# adonis
dist <- phyloseq::distance(subset_samples(pl, Compartment == 'Bulk soil'), method ='wunifrac')
pair_adonis6 <- pairwiseAdonis::pairwise.adonis(dist, sample_data(subset_samples(pl, Compartment == 'Bulk soil'))$Group, perm = 9999)
pair_adonis6
pair_adonis6$Compartment <- 'Bulk soil'

# Interactive Adonis using adonis2() to find interactive effects of multiple factors
sample_data_bulk <- data.frame(sample_data(subset_samples(pl, Compartment == 'Bulk soil'))) # important: data.frame() is needed

inter_adonis_bulk <- adonis2(dist ~ Cd * Mycorrhizal, data = sample_data_bulk, permutations = 999, by = 'terms')

inter_adonis_bulk


rm(dist)

# combine all adonis results 
# ~~~~~~~ Table S2. Results of pairwise PERMANOVA (Adonis) of bacterial communities from various treatments using weighted UniFrac distance. ----
PCoA_adonis_results <- rbind(pair_adonis1, pair_adonis2, 
                             pair_adonis3, pair_adonis4,
                             pair_adonis5, pair_adonis6)

if(save){
  write_xlsx(PCoA_adonis_results, "out/Table S2. PCoA_adonis_results_by_group.xlsx")

}

# ~~~~~~~ Table S3. Results of interactive PERMANOVA (Adonis) assessing the effects of specific factors on bacterial community composition using weighted UniFrac distance metrices ----

three_factor_permutation_rs$Compartment <- 'All'
three_factor_permutation_rs <- rownames_to_column(three_factor_permutation_rs, var = "Factor")

inter_adonis_endo$Compartment <- 'Endosphere'
inter_adonis_endo <- rownames_to_column(inter_adonis_endo, var = "Factor")

inter_adonis_rzp$Compartment <- 'Rhizoplane'
inter_adonis_rzp <- rownames_to_column(inter_adonis_rzp, var = "Factor")

inter_adonis_rzs$Compartment <- 'Rhizosphere'
inter_adonis_rzs <- rownames_to_column(inter_adonis_rzs, var = "Factor")

inter_adonis_hypho$Compartment <- 'Hyphosphere'
inter_adonis_hypho <- rownames_to_column(inter_adonis_hypho, var = "Factor")

inter_adonis_hyphae$Compartment <- 'Hyphae'
inter_adonis_hyphae <- rownames_to_column(inter_adonis_hyphae, var = "Factor")

inter_adonis_bulk$Compartment <- 'Bulk soil'
inter_adonis_bulk <- rownames_to_column(inter_adonis_bulk, var = "Factor")

inter_adonis_rs_all <- rbind(inter_adonis_endo,
                             inter_adonis_rzp,
                             inter_adonis_rzs,
                             inter_adonis_hypho,
                             inter_adonis_hyphae,
                             inter_adonis_bulk)

inter_adonis_rs_all

if(save){
  write_xlsx(inter_adonis_rs_all,
             "out/Table S3. interactive_permanova_all.xlsx")
}



# Strong compartmental effect
# no Cd
dist <- phyloseq::distance(subset_samples(pl, Cd == 0), method ='wunifrac')
sample_data_no_cd <- data.frame(sample_data(subset_samples(pl, Cd == 0)))
adonis2(dist ~ Compartment, data = sample_data_no_cd, permutations = 999)

# with Cd
dist <- phyloseq::distance(subset_samples(pl, Cd == 5), method ='wunifrac')
sample_data_cd <- data.frame(sample_data(subset_samples(pl, Cd == 5)))
adonis2(dist ~ Compartment, data = sample_data_cd, permutations = 999)

# ~~~~~ Combine plots and save -------------------------------------
# ~~~~~~~ Fig. 2. Beta div ----
ggarrange(p1, p2, p3,
          p4, p5, p6,
          labels = c("a", "b", "c", 
                     "d", 'e', 'f'),
          ncol = 3, nrow = 2)
if(save){
ggsave('out/Fig. 2. PCoA_diff_compartments_group_by_treatment_R2.pdf',
       width = 88*3, height = 66*2, units = "mm")

ggsave('out/Fig. 2. PCoA_diff_compartments_group_by_treatment_R2.jpg',
       width = 88*3, height = 66*2, units = "mm")
}
rm(p1, p2, p3, p4, p5, p6)
rm(p1uw, p4bc, p4uw) # delete other exploring plots


# within- and between-group Bray-Curtis distance comparison ----
# subsetting different compartments
pl_rzp<- subset_samples(pl, Compartment == 'Rhizoplane')
pl_rzs<- subset_samples(pl, Compartment == 'Rhizosphere')
pl_hps<- subset_samples(pl, Compartment == 'Hyphosphere')
pl_hph<- subset_samples(pl, Compartment == 'Hyphae')

# can use 'distanceMethodList' to check what additional methods for distance calculation then to replace 'bray' here. 

p1 <- beta_boxplot_btw(pl_rzp, 'bray', 'Cd_VS_Ctrl', 'M+Cd_VS_M')+
  ggtitle('Rhizoplane (p < 0.01)')+
  ylab ('Bray-Curtis dissimilarity')+
  xlab ('Grouping')+
  theme(legend.position="none")+
  scale_x_discrete(labels = c('Cd_VS_Ctrl'='Cd0 vs Cd5\n Non-mycorrhizal','M+Cd_VS_M'='Cd0 vs Cd5\n Mycorrhizal'))

p2 <- beta_boxplot_btw(pl_rzs, 'bray', 'Cd_VS_Ctrl', 'M+Cd_VS_M')+
  ggtitle('Rhizosphere (p = 0.5469)')+
  ylab ('Bray-Curtis dissimilarity')+
  xlab ('Grouping')+
  theme(legend.position="none")+
  scale_x_discrete(labels = c('Cd_VS_Ctrl'='Cd0 vs Cd5\n Non-mycorrhizal','M+Cd_VS_M'='Cd0 vs Cd5\n Mycorrhizal'))

p3 <- beta_boxplot_btw(pl_hps, 'bray', 'Cd_VS_Ctrl', 'M+Cd_VS_M')+
  ggtitle('Hyphosphere (p < 0.01)')+
  ylab ('Bray-Curtis dissimilarity')+
  xlab ('Grouping')+
  theme(legend.position="none")+
  scale_x_discrete(labels = c('Cd_VS_Ctrl'='Cd0 vs Cd5\n Non-mycorrhizal','M+Cd_VS_M'='Cd0 vs Cd5\n Mycorrhizal'))


p4 <- beta_boxplot_btw(pl_hph, 'bray', 'M_VS_M', 'M+Cd_VS_M+Cd')+
  ggtitle('Hyphae (p < 0.05)')+
  ylab ('Bray-Curtis dissimilarity')+
  xlab ('Grouping')+
  theme(legend.position="none")+
  scale_x_discrete(labels = c('M_VS_M'='Cd0\nMycorrhizal','M+Cd_VS_M+Cd'='Cd5\nMycorrhizal'))




# for hyphae, the following plot function generate the same plot, so confirm both plotting functions are valid. 
# beta_boxplot(pl_hph, method = "bray", group = 'Group')
# ## Data
# beta_boxplot_hps$data
# 
# ## Plot
# beta_boxplot_hps$plot

# ~~~~~~~ Fig. 5. between-treatment Bray-Curtis dissimilarity ----
beta_comparison_plot <- ggarrange(p1, p2, p3, p4,
                                  labels = c("a", "b", "c", "d"),
                                  ncol = 2, nrow = 2)
beta_comparison_plot


# t test
BC_t.test(pl_rzp, 'NM_VS_NM_5_0', 'M_VS_M_5_0')
BC_t.test(pl_rzs, 'NM_VS_NM_5_0', 'M_VS_M_5_0')
BC_t.test(pl_hps, 'NM_VS_NM_5_0', 'M_VS_M_5_0')
BC_t.test(pl_hph, 'M_VS_M_5_5', 'M_VS_M_0_0') # note all NM should be excluded

if(save){
  ggsave('out/Fig. 5. Between-trt B-C dist.pdf', plot = beta_comparison_plot,
       width = 88*1.6, height = 66*2, units = "mm")
  ggsave('out/Fig. 5. Between-trt B-C dist.jpg', plot = beta_comparison_plot,
       width = 88*1.6, height = 66*2, units = "mm")
}

rm(p1, p2, p3, p4, beta_comparison_plot)


p1 <- beta_boxplot(pl_rzp, method = "bray", group = 'Group')+
  ggtitle('Rhizoplane (within-group)')+
  ylab ('Bray-Curtis dissimilarity')+
  theme(legend.position = 'none')+
  annotate("text", x = 'Ctrl', y = 0.65, label = "b", size = 4) +
  annotate("text", x = 'Cd', y = 0.65, label = "a", size = 4) +
  annotate("text", x = 'M', y = 0.65, label = "b", size = 4) +
  annotate("text", x = 'M+Cd', y = 0.65, label = "b", size = 4)

p2 <- beta_boxplot(pl_rzs, method = "bray", group = 'Group')+
  ggtitle('Rhizosphere (within-group)')+
  ylab ('Bray-Curtis dissimilarity')+
  theme(legend.position = 'none')+
  annotate("text", x = 'Ctrl', y = 0.7, label = "bc", size = 4) +
  annotate("text", x = 'Cd', y = 0.7, label = "b", size = 4) +
  annotate("text", x = 'M', y = 0.7, label = "c", size = 4) +
  annotate("text", x = 'M+Cd', y = 0.7, label = "a", size = 4)

p3 <- beta_boxplot(pl_hps, method = "bray", group = 'Group')+
  ggtitle('Hyphosphere (within-group)')+
  ylab ('Bray-Curtis dissimilarity')+
  theme(legend.position = 'none')+
  annotate("text", x = 'Ctrl', y = 1.05, label = "b", size = 4) +
  annotate("text", x = 'Cd', y = 1.05, label = "b", size = 4) +
  annotate("text", x = 'M', y = 1.05, label = "a", size = 4) +
  annotate("text", x = 'M+Cd', y = 1.05, label = "a", size = 4)

p4 <- beta_boxplot(pl_hph, method = "bray", group = 'Group')+
  ggtitle('Hyphae (within-group)')+
  ylab ('Bray-Curtis dissimilarity')+
  theme(legend.position = 'none')+
    annotate("text", x = 'M', y = 0.68, label = "a", size = 4) +
  annotate("text", x = 'M+Cd', y = 0.68, label = "b", size = 4)



# ~~~~~~~ Fig. 4. within-treatment Bray-Curtis dissimilarity ----
beta_comparison_plot_within_group <- ggarrange(p1, p2, p3, p4,
                                               labels = c("a", "b", "c", "d"),
                                               ncol = 2, nrow = 2)
beta_comparison_plot_within_group

# ANOVA followed by Tukey's test 
BC_tukey(pl_rzp)
BC_tukey(pl_rzs)
BC_tukey(pl_hps)
BC_tukey(pl_hph)

if(save){
ggsave('out/Fig. 4. Within-trt B-C dist.pdf', 
       plot = beta_comparison_plot_within_group,
       width = 88*1.6, height = 66*2, units = "mm")

ggsave('out/Fig. 4. Within-trt B-C dist.jpg', 
       plot = beta_comparison_plot_within_group,
       width = 88*1.6, height = 66*2, units = "mm")
}

rm(p1, p2, p3, p4, beta_comparison_plot_within_group)
# ~~ Between-group beta diversity comparison FINISHED ----


# Compartmental effects in separating bacterial community ----

p1 <- plot_w_unifrac_cpmt(subset_samples(pl, Group == 'Ctrl'), 'Compartment')+
  ggtitle("PCoA on weighted UniFrac distance (Ctrl)")+
  stat_ellipse()+
  theme(legend.position = 'none')

dist <- phyloseq::distance(subset_samples(pl, Group == 'Ctrl'), method ='wunifrac')
pair_adonis_ctrl <- pairwiseAdonis::pairwise.adonis(dist, sample_data(subset_samples(pl, Group == 'Ctrl'))$Compartment, perm = 999)
pair_adonis_ctrl$Group <- 'Ctrl'


p2 <- plot_w_unifrac_cpmt(subset_samples(pl, Group == 'Cd'), 'Compartment')+
  ggtitle("PCoA on weighted UniFrac distance (Cd)")+
  stat_ellipse()+
  theme(legend.position = "none")

dist <- phyloseq::distance(subset_samples(pl, Group == 'Cd'), method ='wunifrac')
pair_adonis_cd <- pairwiseAdonis::pairwise.adonis(dist, sample_data(subset_samples(pl, Group == 'Cd'))$Compartment, perm = 999)
pair_adonis_cd$Group <- 'Cd'

p3 <- plot_w_unifrac_cpmt(subset_samples(pl, Group == 'M' & Compartment != 'Hyphae'), 'Compartment')+
  ggtitle("PCoA on weighted UniFrac distance (M)")+
  stat_ellipse() +
  theme(legend.position = "none")# excluding Hyphae compartment

p4 <- plot_w_unifrac_cpmt(subset_samples(pl, Group == 'M+Cd' & Compartment != 'Hyphae'), 'Compartment')+
  ggtitle("PCoA on weighted UniFrac distance (M+Cd)")+
  stat_ellipse()+
  theme(legend.position = "none") # excluding Hyphae compartment


# ~~~~~~~ Additional Figure. PCoA of bacteria from diff compts (no hyphae) ----
PCoA_by_compartment_plot_no_hy <- ggarrange(p1, p2, p3, p4,
                                            labels = c("a", "b", "c", "d"),
                                            ncol = 2, nrow = 2)
PCoA_by_compartment_plot_no_hy

# ggsave('out/PCoA_by_compartment_plot_no_hyphae.pdf',
#        width = 88*3, height = 66*2, units = "mm")
# 
# ggsave('out/PCoA_by_compartment_plot_no_hyphae.jpg',
#        width = 88*3, height = 66*2, units = "mm")

# include hyphae compartment - only two groups have hyphae: M and M+Cd
p5<- plot_w_unifrac_cpmt_hy(subset_samples(pl, Group == 'M'), 'Compartment')+
  ggtitle("PCoA on weighted UniFrac distance (M)")+
  stat_ellipse() +
  theme(legend.position = c(0.15, 0.72),
        legend.title = element_blank(),
        legend.spacing.y = unit(0.0, "cm"))

dist <- phyloseq::distance(subset_samples(pl, Group == 'M'), method ='wunifrac')
pair_adonis_M <- pairwiseAdonis::pairwise.adonis(dist, sample_data(subset_samples(pl, Group == 'M'))$Compartment, perm = 999)
pair_adonis_M$Group <- 'M'

p6 <- plot_w_unifrac_cpmt_hy(subset_samples(pl, Group == 'M+Cd'), 'Compartment')+
  ggtitle("PCoA on weighted UniFrac distance (M+Cd)")+
  stat_ellipse() +
  theme(legend.position = "none")

dist <- phyloseq::distance(subset_samples(pl, Group == 'M+Cd'), method ='wunifrac')
pair_adonis_M_cd <- pairwiseAdonis::pairwise.adonis(dist, sample_data(subset_samples(pl, Group == 'M+Cd'))$Compartment, perm = 999)
pair_adonis_M_cd$Group <- 'M+Cd'

PCoA_by_compartment_plot_with_hy <- ggarrange(p5, p6,
                                              labels = c("a", "b"),
                                              ncol = 2, nrow = 1)

PCoA_by_compartment_plot_with_hy

# ggsave('out/PCoA_by_compartment_plot_with_hyphae_R1.pdf',
#        width = 88*2, height = 66, units = "mm")
# 
# ggsave('out/PCoA_by_compartment_plot_with_hyphae_R1.jpg',
#        width = 88*2, height = 66, units = "mm")


# ~~~~~~~ Fig. S3. PCoA_by_compartment_plot_with_hyphae ----
PCoA_by_compartment_plot <- ggarrange(p1, p2, p5, p6,
                                      labels = c("a", "b", "c", "d"),
                                      ncol = 2, nrow = 2)
PCoA_by_compartment_plot

if(save){
ggsave('out/Fig.S3. PCoA_by_compartment_plot.pdf',
       width = 88*3, height = 66*3, units = "mm")

ggsave('out/Fig.S3. PCoA_by_compartment_plot.jpg',
       width = 88*3, height = 66*3, units = "mm")
}
rm(p1, p2, p3, p4, p5, p6)

rm(PCoA_by_compartment_plot,
   PCoA_by_compartment_plot_with_hy,
   PCoA_by_compartment_plot_no_hy)


PCoA_pairwise_adonis_by_cpmt <- rbind(pair_adonis_ctrl,
                                      pair_adonis_cd,
                                      pair_adonis_M,
                                      pair_adonis_M_cd)


PCoA_pairwise_adonis_by_cpmt

# ~~~~~~~ Table S1. Results of pairwise PERMANOVA (Adonis) of bacterial communities across various compartments using weighted UniFrac distance matrices. ----

if(save){
  write_xlsx(PCoA_pairwise_adonis_by_cpmt, "out/Table S1. PCoA_adonis_results_by_compartment.xlsx")
}
# Priming effects? The hyphosphere bacterial community drastically changed (phylogenetic diversity increased, as indicated by weighted UniFrac distance) by mycorrhizal symbiosis, and its phylogenetic diversity (UniFrac distance) decreased after receiving Cd. For the endosphere bacterial community seems becoming more phylogenetically diverse after receiving Cd when the plants colonized by AM fungus.

# Resulting plots were adjusted in AI/CorelDraw. 
