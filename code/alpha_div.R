# Alpha diversity  ----
# ~~~~ (use un-rarefied data, that is the 'pl_un_rarefied' phyloseq object) ----
# Ref: 10.1371/journal.pcbi.1003531
# One can consider the SRS method (ref.: 10.7717/peerj.9593) for scaling the ASV table for alpha div

# ~~ all samples by mycorrhizal ----
alpha_myco <- plot_richness(pl_un_rarefied, x="Compartment", 
                            measures=c("Observed", "Shannon", "Chao1"), 
                            color="Mycorrhizal")+
  scale_color_manual(values = c("darkgrey", "indianred4"))+
  geom_boxplot(outlier.color = NA)+
  # geom_jitter(alpha = 0.3, position = position_dodge(width = 0.75))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = 'none',
        # legend.position = c(0.1, 0.4), 
        legend.box.background = element_rect(color = "black", size = 0.5)) +
  ggtitle('Overall mycorrhizal effect')+
  labs(color = "Myco")

alpha_myco$layers <- alpha_myco$layers[-1]
alpha_myco

# ~~ all samples by Cd ----

alpha_cd <- plot_richness(pl_un_rarefied, x="Compartment", 
                          measures=c("Observed", "Shannon", "Chao1"), 
                          color="Cd")+
  scale_color_manual(values = c("darkgrey", "indianred4"))+
  geom_boxplot(outlier.color = NA)+
  # geom_jitter(alpha = 0.3, position = position_dodge(width = 0.75))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = 'none',
        # legend.position = c(0.08, 0.4), 
        legend.box.background = element_rect(color = "black", size = 0.5))+
  ggtitle('Overall Cd effect')


alpha_cd$layers <- alpha_cd$layers[-1]
alpha_cd

# ~~~~~ Combine plots and save -------------------------------------
# ~~~~~~~ Fig. S4. Alpha div various compartments ----
p_alpha_div_overall <- ggarrange(alpha_myco, alpha_cd,
          labels = c("a", "b"),
          ncol = 1, nrow = 2)
p_alpha_div_overall

if(save){
ggsave('out/Fig. S4. overall_alpha_div_by_M_Cd.pdf', 
       plot = p_alpha_div_overall,
       width = 88*2, height = 66*4, units = "mm")

ggsave('out/Fig. S4. overall_alpha_div_by_M_Cd.jpg',
       plot = p_alpha_div_overall,
       width = 88*2, height = 66*4, units = "mm")
}
rm(alpha_myco, alpha_cd)


# ~~ subset those with AM fungal inoculation ----
# compare Cd and non-Cd
AM <- subset_samples(pl_un_rarefied, Mycorrhizal == 'M')
sample_names(AM)
class(AM)
alpha_AM_cd <- plot_richness(AM, x="Compartment", 
                             measures=c("Observed", "Shannon", "Chao1"), 
                             color="Cd")+
  scale_color_manual(values = c("darkgrey", "indianred4"))+
  geom_boxplot(outlier.color = NA)+
  # geom_jitter(alpha = 0.3, position = position_dodge(width = 0.75))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = 'none',
        legend.box.background = element_rect(color = "black", size = 0.5)) +
  ggtitle('Mycorrhizal samples')

alpha_AM_cd$layers <- alpha_AM_cd$layers[-1]
alpha_AM_cd


# ~~ subset those mock ----
# compare Cd and non-Cd

NM <- subset_samples(pl_un_rarefied, Mycorrhizal == 'NM')
sample_names(NM)
alpha_NM_cd <- plot_richness(NM, x="Compartment", 
                             measures=c("Observed", "Shannon", "Chao1"), 
                             color="Cd")+
  scale_color_manual(values = c("darkgrey", "indianred4"))+
  geom_boxplot(outlier.color = NA)+
  # geom_jitter(alpha = 0.3, position = position_dodge(width = 0.75))+
  # the follow code is added to show Hyphae group (non dected actually, but it is better to show it in the plot)
  scale_x_discrete(limits = c('Endosphere', 'Rhizoplane', 'Rhizosphere',
                              'Hyphosphere', 'Hyphae','Bulk soil'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = 'none', 
        legend.box.background = element_rect(color = "black", size = 0.5))+
  ggtitle('Non-mycorrhizal samples')

alpha_NM_cd$layers <- alpha_NM_cd$layers[-1]
alpha_NM_cd



# ~~ subset all with Cd ----
# compare NM and M

Cd <- subset_samples(pl_un_rarefied, Cd == 5)
sample_names(Cd)
alpha_Cd_by_AM <- plot_richness(Cd, x="Compartment", 
                                measures=c("Observed", "Shannon", "Chao1"), 
                                color="Mycorrhizal")+
  scale_color_manual(values = c("darkgrey", "indianred4"))+
  geom_boxplot(outlier.color = NA)+
  # geom_jitter(alpha = 0.3, position = position_dodge(width = 0.75))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = 'none', 
        legend.box.background = element_rect(color = "black", size = 0.5))+
  ggtitle('Cd samples')+
  labs(color = "Myco")

alpha_Cd_by_AM$layers <- alpha_Cd_by_AM$layers[-1]
alpha_Cd_by_AM

# ~~ subset all without Cd ----
non_Cd <- subset_samples(pl_un_rarefied, Cd == 0)
sample_names(non_Cd)
alpha_non_Cd_by_AM <- plot_richness(non_Cd, x="Compartment", 
                                    measures=c("Observed", "Shannon", "Chao1"), 
                                    color="Mycorrhizal")+
  scale_color_manual(values = c("darkgrey", "indianred4"))+
  geom_boxplot(outlier.color = NA)+
  # geom_jitter(alpha = 0.3, position = position_dodge(width = 0.75))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = 'none',
        legend.box.background = element_rect(color = "black", size = 0.5))+
  ggtitle('Non-Cd samples')+
  labs(color = "Myco")

alpha_non_Cd_by_AM$layers <- alpha_non_Cd_by_AM$layers[-1]
alpha_non_Cd_by_AM


# ~~~~~ Combine plots and save -------------------------------------
# ~~~~~~~ Fig. S. Alpha div (subsetted) ----
ggarrange(alpha_non_Cd_by_AM, alpha_Cd_by_AM,
          alpha_NM_cd, alpha_AM_cd,
          labels = c("a", "b", 'c', 'd'),
          ncol = 2, nrow = 2)

if(save){
ggsave('out/Fig. 3. alpha_div_by_AM_Cd.pdf',
       width = 88*3, height = 66*3, units = "mm")

ggsave('out/Fig. 3. alpha_div_by_AM_Cd.jpg',
       width = 88*3, height = 66*3, units = "mm")
}

rm(NM, AM, Cd, non_Cd, alpha_non_Cd_by_AM, alpha_Cd_by_AM)

# plot Shannon index using plot_richness(). Separate Compartment and the group, and show the name of the group on the x-axis. set color for the group: values = c("#bababa","#404040", "#f4a582", "#ca0020")), remove data points except the jitters
p_a_div <- plot_richness(pl_un_rarefied, x="Compartment", measures=c("Shannon", "Chao1"), color="Group") +
  geom_boxplot(aes(group=interaction(Compartment, Group)), position=position_dodge(width=0.8), outliers = F)+ 
  #geom_jitter(alpha=0.2, position=position_dodge(width=0.8)) +
  scale_color_manual(values=c("#bababa","#404040","#f4a582","#ca0020")) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust = 1.05)) +
  # remove grid lines
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # add gray dashed vertical lines between compartments, should be 5 in total for each index
  geom_vline(xintercept=c(1.5, 2.5, 3.5, 4.5, 5.5), linetype="dashed", color="gray") +
  # Change legend title to "Treatment"
  labs(color = "Treatment") 
# theme(legend.title = element_blank(), legend.position = "none")
  

# Remove the original data points, keep the jitters
p_a_div$layers <- p_a_div$layers[-c(1, 2)]
p_a_div

# ANOVA then tukey test for alpha diversity 
# calculate alpha diversity measures
alpha_div <- estimate_richness(pl_un_rarefied, measures = c("Observed", "Shannon", "Chao1"))

# take Group and rownames of sample_data(pl_un_rarefied)
# add Group and rownames to alpha_div
alpha_div <- cbind(alpha_div, sample_data(pl_un_rarefied)[,c('Group')])
head(alpha_div)
# subset rownames that contain 'ed' then perform ANOVA and Tukey test
alpha_div_ed <- alpha_div[grepl('ed', rownames(alpha_div)),]
HSD.test(aov(Observed ~ Group, data = alpha_div_ed), 'Group')$groups
HSD.test(aov(Chao1 ~ Group, data = alpha_div_ed), 'Group')$groups
HSD.test(aov(Shannon ~ Group, data = alpha_div_ed), 'Group')$groups

# do the same for rows that contain 'rp'
alpha_div_rp <- alpha_div[grepl('rp', rownames(alpha_div)),]
HSD.test(aov(Observed ~ Group, data = alpha_div_rp), 'Group')$groups
HSD.test(aov(Chao1 ~ Group, data = alpha_div_rp), 'Group')$groups
HSD.test(aov(Shannon ~ Group, data = alpha_div_rp), 'Group')$groups

# do the same for rows that contain 'rs'
alpha_div_rs <- alpha_div[grepl('rs', rownames(alpha_div)),]
HSD.test(aov(Observed ~ Group, data = alpha_div_rs), 'Group')$groups
HSD.test(aov(Chao1 ~ Group, data = alpha_div_rs), 'Group')$groups
HSD.test(aov(Shannon ~ Group, data = alpha_div_rs), 'Group')$groups

# do the same for rows that contain 'hps'
alpha_div_hps <- alpha_div[grepl('hps', rownames(alpha_div)),]
HSD.test(aov(Observed ~ Group, data = alpha_div_hps), 'Group')$groups
HSD.test(aov(Chao1 ~ Group, data = alpha_div_hps), 'Group')$groups
HSD.test(aov(Shannon ~ Group, data = alpha_div_hps), 'Group')$groups

# do the same for rows that contain 'hp_'
alpha_div_hp <- alpha_div[grepl('hp_', rownames(alpha_div)),]
HSD.test(aov(Observed ~ Group, data = alpha_div_hp), 'Group')$groups
HSD.test(aov(Chao1 ~ Group, data = alpha_div_hp), 'Group')$groups
HSD.test(aov(Shannon ~ Group, data = alpha_div_hp), 'Group')$groups

# do the same for rows that contain 'bs'
alpha_div_bs <- alpha_div[grepl('bs', rownames(alpha_div)),]
HSD.test(aov(Observed ~ Group, data = alpha_div_bs), 'Group')$groups
HSD.test(aov(Chao1 ~ Group, data = alpha_div_bs), 'Group')$groups
HSD.test(aov(Shannon ~ Group, data = alpha_div_bs), 'Group')$groups





# combine with rel_ab_avg and p_a_div
# ~~~~~~~ Fig. 1b. Relative abundance and alpha diversity ----
p_a_div
rel_ab_avg
p_rel_a <- ggarrange(rel_ab_avg, p_a_div,
          labels = c("a", "b"),
          ncol = 1, nrow = 2)
p_rel_a
#save the plot if save is TRUE
if(save){ggsave('out/Fig. 1R. rel_ab_and_a-div.jpg', 
       plot = p_rel_a,
       width = 88*2.3, height = 66*3, units = "mm", dpi = 300)
  
  ggsave('out/Fig. 1R. rel_ab_and_a-div.pdf', 
         plot = p_rel_a,
         width = 88*2.3, height = 66*3, units = "mm", dpi = 300)
}

rm(p_a_div, rel_ab_avg)
