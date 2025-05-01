# Strong and weak correlation analysis ----

# Strong or weak interactions are defined using correlation coefficient and slope of linear regression.

# Strategy for analysis ----
# 1. Agglomerate taxa of the same type using un-rarefied OTU table.
# 1. Take OTU table.
# 2. Subset target group.
# 3. Analyze sub-data set one by one:
#   a) For each sub-data set, there are 8 replicates. For each OTU, if the count sum of the 8 replicates is smaller than a threshold, the OTU will be discarded. Take a same threshold for all sub-data sets first. 
#   b) Conduct OTU pairwise Pearson correlation, linear regression, obtain slopes of regression lines, using the function str_wk()
#   c) Categorize interaction type and strength based on correlation statistics and slopes.


# START the analysis ----

# subset those with p-value =< 0.05. These are corr with significance

# weak interaction contribute to system stability

# we need to define weak interaction here. In complex ecosystems, strong interaction suggests species significantly influence each other's abundance and behavior (such as predator-prey, competition, and mutualism). On the other hand, weak interaction suggests less impact and may not strongly affect population dynamics.

# So we look at the histogram of the p-value and cor. coeff. of the linear regression of the otu pairwise correlation. 

# a higher frequency of low p-value and high R2, suggest more strong interactions in the system. 


# ~~ strong and weak interaction/corr ----
# use the function str_wk() ----


# ~ Bulk soil ----
sw.rs.bs <- str_wk(pl_un_rarefied, gl_tax = 'Family',
                      target_group = 'bs_ck_[1-8]$|bs_f_[1-8]$', c.cutoff = 5) 

sw.rs.bs.cd <- str_wk(pl_un_rarefied, gl_tax = 'Family',
                      target_group = 'bs_cd_[1-8]$|bs_f_cd_[1-8]$', c.cutoff = 5) 

# ~~~~ Fig. S16 ----
plot_sw(sw.rs.bs)
plot_sw(sw.rs.bs.cd)

plot_sw_sig(sw.rs.bs)
plot_sw_sig(sw.rs.bs.cd)

# ~ Hyphosphere ----

sw.rs.hps.cd <- str_wk(pl_un_rarefied, gl_tax = 'Family',
       target_group = 'hps_cd_[1-8]$', c.cutoff = 5) 

sw.rs.hps.f.cd <- str_wk(pl_un_rarefied, gl_tax = 'Family',
       target_group = 'hps_f_cd_[1-8]$', c.cutoff = 5) 


# ~~~~ Fig. 6d-f ----
plot_sw(sw.rs.hps.cd)
plot_sw(sw.rs.hps.f.cd)

plot_sw_sig(sw.rs.hps.cd)
plot_sw_sig(sw.rs.hps.f.cd)


sw.rs.hps.ck <- str_wk(pl_un_rarefied, gl_tax = 'Family',
       target_group = 'hps_ck_[1-8]$', c.cutoff = 5) 

sw.rs.hps.f <- str_wk(pl_un_rarefied, gl_tax = 'Family',
       target_group = 'hps_f_[1-8]$', c.cutoff = 5) 

# ~~~~ Extended Data Fig. 4. ----
plot_sw(sw.rs.hps.ck)
plot_sw(sw.rs.hps.f)

plot_sw_sig(sw.rs.hps.ck)
plot_sw_sig(sw.rs.hps.f)

# ~ Rhizosphere ----

sw.rs.rhs.cd <- str_wk(pl_un_rarefied, gl_tax = 'Family',
                       target_group = 'rs_cd_[1-8]$', c.cutoff = 5) 
sw.rs.rhs.f.cd <- str_wk(pl_un_rarefied, gl_tax = 'Family',
                         target_group = 'rs_f_cd_[1-8]$', c.cutoff = 5) 

# ~~~~ Fig. 6a-c ----

plot_sw(sw.rs.rhs.cd)
plot_sw(sw.rs.rhs.f.cd)

plot_sw_sig(sw.rs.rhs.cd)
plot_sw_sig(sw.rs.rhs.f.cd)

sw.rs.rhs.ck <- str_wk(pl_un_rarefied, gl_tax = 'Family',
                       target_group = 'rs_ck_[1-8]$', c.cutoff = 5) 
sw.rs.rhs.f <- str_wk(pl_un_rarefied, gl_tax = 'Family',
                      target_group = 'rs_f_[1-8]$', c.cutoff = 5) 

# ~~~~ Fig. S17 ----
plot_sw(sw.rs.rhs.ck)
plot_sw(sw.rs.rhs.cd)

plot_sw_sig(sw.rs.rhs.ck)
plot_sw_sig(sw.rs.rhs.cd)

# ~~~~ Fig. S18 ----
plot_sw(sw.rs.rhs.ck)
plot_sw(sw.rs.rhs.f)

plot_sw_sig(sw.rs.rhs.ck)
plot_sw_sig(sw.rs.rhs.f)


# ~ Hyphae ----
sw.rs.hph.f <- str_wk(pl_un_rarefied, gl_tax = 'Family',
                      target_group = 'hp_f_[1-8]$', c.cutoff = 5) 
sw.rs.hph.f.cd <- str_wk(pl_un_rarefied, gl_tax = 'Family',
                     target_group = 'hp_f_cd_[1-8]$', c.cutoff = 5) 

# ~~~~ Fig. S19 ----
plot_sw(sw.rs.hph.f)
plot_sw(sw.rs.hph.f.cd)

plot_sw_sig(sw.rs.hph.f)
plot_sw_sig(sw.rs.hph.f.cd)

# save results ----
# List of rs
rs.list <- list(sw.rs.bs,
                sw.rs.bs.cd,
                
                sw.rs.hps.ck,
                sw.rs.hps.cd,
                sw.rs.hps.f,
                sw.rs.hps.f.cd,
                
                sw.rs.rhs.ck,
                sw.rs.rhs.cd,
                sw.rs.rhs.f,
                sw.rs.rhs.f.cd,
                
                sw.rs.hph.f,
                sw.rs.hph.f.cd)

rs.names <- c("sw.rs.bs",
              "sw.rs.bs.cd",
              "sw.rs.hps.ck",
              "sw.rs.hps.cd",
              "sw.rs.hps.f",
              "sw.rs.hps.f.cd",
              "sw.rs.rhs.ck",
              "sw.rs.rhs.cd",
              "sw.rs.rhs.f",
              "sw.rs.rhs.f.cd",
              "sw.rs.hph.f",
              "sw.rs.hph.f.cd")

# Corresponding file names
sw.rs <- c('out/strong_weak_interaction_rs.bs.xlsx',
           'out/strong_weak_interaction_rs.bs.cd.xlsx',
           
           'out/strong_weak_interaction_hps.ck.xlsx',
           'out/strong_weak_interaction_hps.cd.xlsx',
           'out/strong_weak_interaction_hps.f.xlsx',
           'out/strong_weak_interaction_hps.f.cd.xlsx',
           
           'out/strong_weak_interaction_rhs.ck.xlsx',
           'out/strong_weak_interaction_rhs.cd.xlsx',
           'out/strong_weak_interaction_rhs.f.xlsx',
           'out/strong_weak_interaction_rhs.f.cd.xlsx',
           
           'out/strong_weak_interaction_hph.f.xlsx',
           'out/strong_weak_interaction_hph.f.cd.xlsx'
           
           )

# Loop to save each
for (i in seq_along(rs.list)) {
  write_xlsx(rs.list[[i]], sw.rs[i])
}

# significant percentage calculation ----
for (i in seq_along(rs.list)) {
  print(rs.names[i])
  sig.p.perc(rs.list[[i]])
}

# slope angle categorization ----
# Bulk soil 
a.bs.ck <- angle.perc(sw.rs.bs)
a.bs.cd <- angle.perc(sw.rs.bs.cd)

# Hyphosphere 
a.hps.cd <- angle.perc(sw.rs.hps.cd)
a.hps.f.cd <- angle.perc(sw.rs.hps.f.cd)
a.hps.ck <- angle.perc(sw.rs.hps.ck)
a.hps.f <- angle.perc(sw.rs.hps.f)

# Rhizosphere        
a.rhs.ck<- angle.perc(sw.rs.rhs.ck)
a.rhs.cd<- angle.perc(sw.rs.rhs.cd)
a.rhs.f<- angle.perc(sw.rs.rhs.f)
a.rhs.f.cd<- angle.perc(sw.rs.rhs.f.cd)

# Hyphae             
a.hph.f<- angle.perc(sw.rs.hph.f)
a.hph.f.cd<- angle.perc(sw.rs.hph.f.cd)

# plot slope angle proportion ----

p1 <- plot_pie(a.bs.ck)+labs(title= 'Bulk soil - Ctrl')
p2 <- plot_pie(a.bs.cd)+labs(title= 'Bulk soil - Cd')


p3 <- plot_pie(a.hps.ck)+labs(title= 'Hyphoshere - Ctrl')
p4 <- plot_pie(a.hps.cd)+labs(title= 'Hyphoshere - Cd')
p5 <- plot_pie(a.hps.f)+labs(title= 'Hyphoshere - M')
p6 <- plot_pie(a.hps.f.cd)+labs(title= 'Hyphoshere - M+Cd')


p7 <- plot_pie(a.rhs.ck)+labs(title= 'Rhizosphere - Ctrl')
p8 <- plot_pie(a.rhs.cd)+labs(title= 'Rhizosphere - Cd')
p9 <- plot_pie(a.rhs.f)+labs(title= 'Rhizosphere - M')
p10 <- plot_pie(a.rhs.f.cd)+labs(title= 'Rhizosphere - M+Cd')

p11 <- plot_pie(a.hph.f)+labs(title= 'Hyphae - M')
p12 <- plot_pie(a.hph.f.cd)+labs(title= 'Hyphae - M+Cd')


# ~~~~ Fig. 6a, S16a, S17a, S18a, and S19a ----
a.p.A <- ggarrange(p1, p2, p3,
          p4, p5, p6,
          labels = c('a', 'b', 'c',
                     'd', 'e','f'),
          nrow = 3, ncol = 2)

if(save){
  ggsave('out/angle_combined_A.pdf', plot = a.p.A,
       w = 88*2, h = 66*3, unit = 'mm')
}

a.p.B <- ggarrange(p7, p8, p9,
                   p10, p11, p12,
                   labels = c('a', 'b', 'c',
                              'd', 'e','f'),
                   nrow = 3, ncol = 2)

if(save){
ggsave('out/angle_combined_B.pdf', plot = a.p.B,
       w = 88*2, h = 66*3, unit = 'mm')
}


rm(list = paste0("p", 1:12), a.p.A, a.p.B)


# correlation coefficient grouping ----
# Bulk soil 
cc.perc.rs.bs <- cor_coeff_perc(sw.rs.bs)
cc.perc.rs.bs.cd <- cor_coeff_perc(sw.rs.bs.cd)

# Hyphosphere 
cc.perc.hps.cd <- cor_coeff_perc(sw.rs.hps.cd)
cc.perc.hps.f.cd <- cor_coeff_perc(sw.rs.hps.f.cd)
cc.perc.hps.ck <- cor_coeff_perc(sw.rs.hps.ck)
cc.perc.hps.f <- cor_coeff_perc(sw.rs.hps.f)

# Rhizosphere        
cc.perc.rhs.ck <- cor_coeff_perc(sw.rs.rhs.ck)
cc.perc.rhs.cd <- cor_coeff_perc(sw.rs.rhs.cd)
cc.perc.rhs.f <- cor_coeff_perc(sw.rs.rhs.f)
cc.perc.rhs.f.cd <- cor_coeff_perc(sw.rs.rhs.f.cd)

# Hyphae             
cc.perc.hph.f <- cor_coeff_perc(sw.rs.hph.f)
cc.perc.hph.f.cd <- cor_coeff_perc(sw.rs.hph.f.cd)

# combine
cc.list <- list(cc.perc.rs.bs, 
                cc.perc.rs.bs.cd, 
                   
                cc.perc.hps.ck, 
                cc.perc.hps.cd, 
                cc.perc.hps.f.cd, 
                cc.perc.hps.f, 
                   
                cc.perc.rhs.ck, 
                cc.perc.rhs.cd, 
                cc.perc.rhs.f, 
                cc.perc.rhs.f.cd, 
                   
                cc.perc.hph.f, 
                cc.perc.hph.f.cd)

# Merge all dataframes by the column 'cor_coeff_group'
cc.perc.tab <- Reduce(function(x, y) merge(x, y, by = "cor.coeff.group"), cc.list)

cc.perc.tab$cor.coeff.group <- factor(cc.perc.tab$cor.coeff.group, levels = c('-1.0 ~ -0.7', '-0.7 ~ -0.3', '-0.3 ~ 0.0', '0.0 ~ 0.3', '0.3 ~ 0.7', '0.7 ~ 1.0'))

library(reshape2)
cc.perc.tab.melt <- melt(cc.perc.tab, id.vars = "cor.coeff.group", variable.name = "group", value.name = "perc.")

cc.perc.tab.melt$group <- factor(cc.perc.tab.melt$group, 
                                 levels = c('sw.rs.bs',
                                            'sw.rs.bs.cd',
                                            
                                            'sw.rs.hps.ck',
                                            'sw.rs.hps.cd',
                                            'sw.rs.hps.f',
                                            'sw.rs.hps.f.cd',
                                            
                                            'sw.rs.rhs.ck',
                                            'sw.rs.rhs.cd',
                                            'sw.rs.rhs.f',
                                            'sw.rs.rhs.f.cd',
                                            
                                            'sw.rs.hph.f',
                                            'sw.rs.hph.f.cd'))



# Define the list of arguments
g.list <- c('bs', 'rhs', 'hps', 'hph')
plot_list <- list()

for (g in g.list) {
  plot <- cc.compare(cc.perc.tab.melt, g)
  plot <- plot + theme(legend.position = c(0.1, 0.8),
                       legend.title = element_blank())
  plot_list[[g]] <- plot
}
# Print each plot
# for (plot in plot_list) {
#   print(plot)
# }

# ~~ Fig. Summary of the percentages of all Pearson's rho ----
rho_perc <- ggarrange(plotlist = plot_list, ncol = 2, nrow = 2,
          labels = c('a', 'b', 'c', 'd'))

if(save){
ggsave('out/Summary rho percentage.jpg',
       w = 88*4.5, h = 66*4, unit = 'mm')
ggsave('out/Summary rho percentage.pdf',
       w = 88*4.5, h = 66*4, unit = 'mm')
}

# Diversity of interaction type ----
# ~~ Rough esmimation ----
# ~~~~ Evenness of interaction type ----
even_rs <- data.frame(
  Group = c("Bulk soil - Ctrl", 
            "Bulk soil - Cd",
            "Hyphosphere - Ctrl",
            "Hyphosphere - Cd",
            "Hyphosphere - M",
            "Hyphosphere - M+Cd",
            
            "Rhizosphere - Ctrl",
            "Rhizosphere - Cd",
            "Rhizosphere - M",
            "Rhizosphere - M+Cd",
            
            "Hyphae - M",
            "Hyphae - M+Cd"),
  
  Evenness = c(calculate_evenness(sw.rs.bs),
               calculate_evenness(sw.rs.bs.cd),
               
               calculate_evenness(sw.rs.hps.ck),
               calculate_evenness(sw.rs.hps.cd),
               calculate_evenness(sw.rs.hps.f),
               calculate_evenness(sw.rs.hps.f.cd),
               
               calculate_evenness(sw.rs.rhs.ck),
               calculate_evenness(sw.rs.rhs.cd),
               calculate_evenness(sw.rs.rhs.f),
               calculate_evenness(sw.rs.rhs.f.cd),
               
               calculate_evenness(sw.rs.hph.f),
               calculate_evenness(sw.rs.hph.f.cd))
)

# Print the results
print(even_rs)

# ~~~~ Chao1 of interaction type ----
chao1_rs <- data.frame(
  Group = c("Bulk soil - Ctrl", 
            "Bulk soil - Cd",
            "Hyphosphere - Ctrl",
            "Hyphosphere - Cd",
            "Hyphosphere - M",
            "Hyphosphere - M+Cd",
            
            "Rhizosphere - Ctrl",
            "Rhizosphere - Cd",
            "Rhizosphere - M",
            "Rhizosphere - M+Cd",
            
            "Hyphae - M",
            "Hyphae - M+Cd"),
  
  Chao1 = c(cal.chao1(sw.rs.bs),
                cal.chao1(sw.rs.bs.cd),
               
                cal.chao1(sw.rs.hps.ck),
                cal.chao1(sw.rs.hps.cd),
                cal.chao1(sw.rs.hps.f),
                cal.chao1(sw.rs.hps.f.cd),
               
                cal.chao1(sw.rs.rhs.ck),
                cal.chao1(sw.rs.rhs.cd),
                cal.chao1(sw.rs.rhs.f),
                cal.chao1(sw.rs.rhs.f.cd),
               
                cal.chao1(sw.rs.hph.f),
                cal.chao1(sw.rs.hph.f.cd))
)

# Print the results
print(chao1_rs)



# ~~~~ ACE index of interaction type ----

ace_rs <- data.frame(
  Group = c("Bulk soil - Ctrl", 
            "Bulk soil - Cd",
            "Hyphosphere - Ctrl",
            "Hyphosphere - Cd",
            "Hyphosphere - M",
            "Hyphosphere - M+Cd",
            
            "Rhizosphere - Ctrl",
            "Rhizosphere - Cd",
            "Rhizosphere - M",
            "Rhizosphere - M+Cd",
            
            "Hyphae - M",
            "Hyphae - M+Cd"),
  
  ACE = c(cal.ace(sw.rs.bs),
          cal.ace(sw.rs.bs.cd),
            
          cal.ace(sw.rs.hps.ck),
          cal.ace(sw.rs.hps.cd),
          cal.ace(sw.rs.hps.f),
          cal.ace(sw.rs.hps.f.cd),
            
          cal.ace(sw.rs.rhs.ck),
          cal.ace(sw.rs.rhs.cd),
          cal.ace(sw.rs.rhs.f),
          cal.ace(sw.rs.rhs.f.cd),
            
          cal.ace(sw.rs.hph.f),
          cal.ace(sw.rs.hph.f.cd))
)

# Print the results
print(ace_rs)


temp_merge <- merge(even_rs, chao1_rs, by = 'Group')


merge(temp_merge, ace_rs, by = 'Group')

# ~~ Bootstrapping ----

# Hyphosphere ----
# ~ Bootstrapping Chao1 ----
# Generate bootstrap replicates
set.seed(123)
boot_hps_ctrl <- boot(group.slope(sw.rs.hps.ck$slope_a), boot_chao1, R = 300)
boot_hps_cd <- boot(group.slope(sw.rs.hps.cd$slope_a), boot_chao1, R = 300)
boot_hps_f <- boot(group.slope(sw.rs.hps.f$slope_a), boot_chao1, R = 300)
boot_hps_f_cd <- boot(group.slope(sw.rs.hps.f.cd$slope_a), boot_chao1, R = 300)

boot_hps_chao1_rs <- data.frame(
  Chao1 = c(boot_hps_ctrl$t, boot_hps_cd$t, boot_hps_f$t, boot_hps_f_cd$t),
  Group = factor(rep(c("Ctrl", "Cd", "M", 'M+Cd'), each = 300))
)
boot_hps_chao1_rs$Group <- factor(boot_hps_chao1_rs$Group, levels = c("Ctrl", "Cd", "M", "M+Cd"))


hps_chao1_aov_rs <- aov(Chao1 ~ Group, data = boot_hps_chao1_rs)
summary(hps_chao1_aov_rs)

HSD.test(hps_chao1_aov_rs, 'Group', console=TRUE)

p_boot_chao1 <- ggplot(boot_hps_chao1_rs, aes(x = Group, y = Chao1))+
  geom_boxplot(outlier.size = .1)+
  geom_jitter(alpha = 0.1, size = .1)+
  annotate("text", x = 'Cd', y = 630, label = "a") +
  annotate("text", x = 'Ctrl', y = 630, label = "b") +
  annotate("text", x = 'M+Cd', y = 630, label = "b") +
  annotate("text", x = 'M', y = 630, label = "c")+
  ylab('Diversity of interaction type (Chao1)')+
  labs(title = "hps - Chao1")

# ~ Bootstrapping evenness ----
set.seed(123)
boot_ev_hps_ctrl <- boot(group.slope(sw.rs.hps.ck$slope_a), boot_even, R = 300)
boot_ev_hps_cd <- boot(group.slope(sw.rs.hps.cd$slope_a), boot_even, R = 300)
boot_ev_hps_f <- boot(group.slope(sw.rs.hps.f$slope_a), boot_even, R = 300)
boot_ev_hps_f_cd <- boot(group.slope(sw.rs.hps.f.cd$slope_a), boot_even, R = 300)

boot_hps_ev_rs <- data.frame(
  Evenness = c(boot_ev_hps_ctrl$t, boot_ev_hps_cd$t, boot_ev_hps_f$t, boot_ev_hps_f_cd$t),
  Group = factor(rep(c("Ctrl", "Cd", "M", 'M+Cd'), each = 300))
)
boot_hps_ev_rs$Group <- factor(boot_hps_ev_rs$Group, levels = c("Ctrl", "Cd", "M", "M+Cd"))

hps_ev_aov_rs <- aov(Evenness ~ Group, data = boot_hps_ev_rs)
summary(hps_ev_aov_rs)

HSD.test(hps_ev_aov_rs, 'Group', console=TRUE)

p_boot_even <-ggplot(boot_hps_ev_rs, aes(x = Group, y = Evenness))+
  geom_boxplot(outlier.size = .1)+
  geom_jitter(alpha = 0.1, size = .1)+
  annotate("text", x = 'Cd', y = 0.65, label = "a") +
  annotate("text", x = 'Ctrl', y = 0.65, label = "b") +
  annotate("text", x = 'M+Cd', y = 0.65, label = "b") +
  annotate("text", x = 'M', y = 0.65, label = "c")+
  ylab('Diversity of interaction type (Evenness)')+
  labs(title = "hps - Evenness")


# ~ Bootstrapping ACE ----
set.seed(123)
boot_ace_hps_ctrl <- boot(group.slope(sw.rs.hps.ck$slope_a), boot_ace, R = 300)
boot_ace_hps_cd <- boot(group.slope(sw.rs.hps.cd$slope_a), boot_ace, R = 300)
boot_ace_hps_f <- boot(group.slope(sw.rs.hps.f$slope_a), boot_ace, R = 300)
boot_ace_hps_f_cd <- boot(group.slope(sw.rs.hps.f.cd$slope_a), boot_ace, R = 300)

boot_hps_ace_rs <- data.frame(
  ACE = c(boot_ace_hps_ctrl$t, boot_ace_hps_cd$t, boot_ace_hps_f$t, boot_ace_hps_f_cd$t),
  Group = factor(rep(c("Ctrl", "Cd", "M", 'M+Cd'), each = 300))
)
boot_hps_ace_rs$Group <- factor(boot_hps_ace_rs$Group, levels = c("Ctrl", "Cd", "M", "M+Cd"))

hps_ace_aov_rs <- aov(ACE ~ Group, data = boot_hps_ace_rs)
summary(hps_ace_aov_rs)

HSD.test(hps_ace_aov_rs, 'Group', console=TRUE)

p_boot_ace <-ggplot(boot_hps_ace_rs, aes(x = Group, y = ACE))+
  geom_boxplot(outlier.size = .1)+
  geom_jitter(alpha = 0.1, size = .1)+
  annotate("text", x = 'Cd', y = 600, label = "a") +
  annotate("text", x = 'Ctrl', y = 600, label = "b") +
  annotate("text", x = 'M+Cd', y = 600, label = "b") +
  annotate("text", x = 'M', y = 600, label = "c")+
  ylab('Diversity of interaction type (ACE)')+
  labs(title = "hps - ACE")


hps_inter_div <- ggarrange(p_boot_even, p_boot_chao1, p_boot_ace,
                           labels = c('d', 'e', 'f'), 
                           ncol = 3, nrow = 1)
if(save){
  ggsave('out/hps_diversity_interaction_type.jpg', plot = hps_inter_div,
       w = 88*2, h = 88, units = 'mm', dpi = 300)
ggsave('out/hps_diversity_interaction_type.pdf', plot = hps_inter_div,
       w = 88*2, h = 88, units = 'mm', dpi = 300)
}
# Rhizosphere ----
# ~ Bootstrapping Chao1 ----

# Generate bootstrap replicates
set.seed(123)
boot_rhs_ctrl <- boot(group.slope(sw.rs.rhs.ck$slope_a), boot_chao1, R = 300)
boot_rhs_cd <- boot(group.slope(sw.rs.rhs.cd$slope_a), boot_chao1, R = 300)
boot_rhs_f <- boot(group.slope(sw.rs.rhs.f$slope_a), boot_chao1, R = 300)
boot_rhs_f_cd <- boot(group.slope(sw.rs.rhs.f.cd$slope_a), boot_chao1, R = 300)

boot_rhs_chao1_rs <- data.frame(
  Chao1 = c(boot_rhs_ctrl$t, boot_rhs_cd$t, boot_rhs_f$t, boot_rhs_f_cd$t),
  Group = factor(rep(c("Ctrl", "Cd", "M", 'M+Cd'), each = 300))
)
boot_rhs_chao1_rs$Group <- factor(boot_rhs_chao1_rs$Group, levels = c("Ctrl", "Cd", "M", "M+Cd"))


rhs_chao1_aov_rs <- aov(Chao1 ~ Group, data = boot_rhs_chao1_rs)
summary(rhs_chao1_aov_rs)

HSD.test(rhs_chao1_aov_rs, 'Group', console=TRUE)

p_rhs_boot_chao1 <- ggplot(boot_rhs_chao1_rs, aes(x = Group, y = Chao1))+
  geom_boxplot(outlier.size = .1)+
  geom_jitter(alpha = 0.1, size = .1)+
  annotate("text", x = 'Cd', y = 750, label = "b") +
  annotate("text", x = 'Ctrl', y = 750, label = "b") +
  annotate("text", x = 'M+Cd', y = 750, label = "c") +
  annotate("text", x = 'M', y = 750, label = "a")+
  ylab('Diversity of interaction type (Chao1)')+
  labs(title = "rhs - Chao1")

# ~ Bootstrapping evenness ----
set.seed(123)
boot_ev_rhs_ctrl <- boot(group.slope(sw.rs.rhs.ck$slope_a), boot_even, R = 300)
boot_ev_rhs_cd <- boot(group.slope(sw.rs.rhs.cd$slope_a), boot_even, R = 300)
boot_ev_rhs_f <- boot(group.slope(sw.rs.rhs.f$slope_a), boot_even, R = 300)
boot_ev_rhs_f_cd <- boot(group.slope(sw.rs.rhs.f.cd$slope_a), boot_even, R = 300)

boot_rhs_ev_rs <- data.frame(
  Evenness = c(boot_ev_rhs_ctrl$t, boot_ev_rhs_cd$t, boot_ev_rhs_f$t, boot_ev_rhs_f_cd$t),
  Group = factor(rep(c("Ctrl", "Cd", "M", 'M+Cd'), each = 300))
)
boot_rhs_ev_rs$Group <- factor(boot_rhs_ev_rs$Group, levels = c("Ctrl", "Cd", "M", "M+Cd"))

rhs_ev_aov_rs <- aov(Evenness ~ Group, data = boot_rhs_ev_rs)
summary(rhs_ev_aov_rs)

HSD.test(rhs_ev_aov_rs, 'Group', console=TRUE)

p_rhs_boot_even <-ggplot(boot_rhs_ev_rs, aes(x = Group, y = Evenness))+
  geom_boxplot(outlier.size = .1)+
  geom_jitter(alpha = 0.1, size = .1)+
  annotate("text", x = 'Cd', y = 0.65, label = "a") +
  annotate("text", x = 'Ctrl', y = 0.65, label = "c") +
  annotate("text", x = 'M+Cd', y = 0.65, label = "b") +
  annotate("text", x = 'M', y = 0.65, label = "c")+
  ylab('Diversity of interaction type (Evenness)')+
  labs(title = "rhs - Evenness")


# ~ Bootstrapping ACE ----
set.seed(123)
boot_ace_rhs_ctrl <- boot(group.slope(sw.rs.rhs.ck$slope_a), boot_ace, R = 300)
boot_ace_rhs_cd <- boot(group.slope(sw.rs.rhs.cd$slope_a), boot_ace, R = 300)
boot_ace_rhs_f <- boot(group.slope(sw.rs.rhs.f$slope_a), boot_ace, R = 300)
boot_ace_rhs_f_cd <- boot(group.slope(sw.rs.rhs.f.cd$slope_a), boot_ace, R = 300)

boot_rhs_ace_rs <- data.frame(
  ACE = c(boot_ace_rhs_ctrl$t, boot_ace_rhs_cd$t, boot_ace_rhs_f$t, boot_ace_rhs_f_cd$t),
  Group = factor(rep(c("Ctrl", "Cd", "M", 'M+Cd'), each = 300))
)
boot_rhs_ace_rs$Group <- factor(boot_rhs_ace_rs$Group, levels = c("Ctrl", "Cd", "M", "M+Cd"))

rhs_ace_aov_rs <- aov(ACE ~ Group, data = boot_rhs_ace_rs)
summary(rhs_ace_aov_rs)

HSD.test(rhs_ace_aov_rs, 'Group', console=TRUE)

p_rhs_boot_ace <-ggplot(boot_rhs_ace_rs, aes(x = Group, y = ACE))+
  geom_boxplot(outlier.size = .1)+
  geom_jitter(alpha = 0.1, size = .1)+
  annotate("text", x = 'Cd', y = 600, label = "c") +
  annotate("text", x = 'Ctrl', y = 600, label = "a") +
  annotate("text", x = 'M+Cd', y = 600, label = "d") +
  annotate("text", x = 'M', y = 600, label = "b")+
  ylab('Diversity of interaction type (ACE)')+
  labs(title = "rhs - ACE")


rhs_inter_div <- ggarrange(p_rhs_boot_even, p_rhs_boot_chao1, p_rhs_boot_ace,
                           labels = c('a', 'b', 'c'), 
                           ncol = 3, nrow = 1)
rhs_inter_div
if(save){
  ggsave('out/rhs_diversity_interaction_type.jpg', plot = rhs_inter_div,
       w = 88*2, h = 88, units = 'mm', dpi = 300)

ggsave('out/rhs_diversity_interaction_type.pdf', plot = rhs_inter_div,
       w = 88*2, h = 88, units = 'mm', dpi = 300)
}

# ~~~~ Fig. 7. Diversity of B-B interaction types ----
rhs_hps_inter_div <- ggarrange(rhs_inter_div, hps_inter_div,
          ncol = 1, nrow = 2)
rhs_hps_inter_div

if(save){
ggsave('out/Fig. 7. rhs_hps_diversity_of_interaction_type.jpg', plot = rhs_hps_inter_div,
       w = 88*2, h= 88*2, units = 'mm', dpi = 300)

ggsave('out/Fig. 7. rhs_hps_diversity_of_interaction_type.pdf', plot = rhs_hps_inter_div,
       w = 88*2, h= 88*2, units = 'mm', dpi = 300)
}

# END ----