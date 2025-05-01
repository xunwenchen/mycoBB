# load soil and plant data ----
soil_data <- read_xlsx('data/soil_data.xlsx')
soil_data$group2 <- factor(soil_data$group2, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))

plant_data <- read_xlsx('data/plant_data.xlsx')
plant_data$group <- factor(plant_data$group, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))
# plant_data$group <- factor(plant_data$group, levels = c('Ctrl', 'M', 'M+Cd', 'Cd'))

#################################################

# Plant traits, AMF colonization ----
# use own function(s) to plot
p1 <- plot_box(plant_data, group, colonization_rate, HM)+
  xlab('Treatment')+
  ylab('Colonization rate (%)')+geom_jitter(alpha = 0.4)

p2 <- plot_box(plant_data, group, plant_height)+
  xlab('Treatment')+
  ylab('Plant height (cm)')+geom_jitter(alpha = 0.4)

p3 <- plot_box(plant_data, group, sht_dry_mass)+
  xlab('Treatment')+
  ylab('Shoot dry mass (g/pot)')+geom_jitter(alpha = 0.4)

p4 <- plot_box(plant_data, group, grain_dry_mass)+
  xlab('Treatment')+
  ylab('Grain yield (g/pot DW)')+geom_jitter(alpha = 0.4)

# ~~~~~~ Fig. S2. AM fungal colonization rate and plant biomass ----

plant_trait_plot <- ggarrange(p1, p2, p3, p4,
                              labels = c("a", "b", "c", "d"),
                              ncol = 2, nrow = 2)
plant_trait_plot

if(save){
ggsave('out/Fig. S2. AM col rate and plant biomass.pdf', plot = plant_trait_plot, width = 88*2, height = 66*2, units = "mm")

ggsave('out/Fig. S2. AM col rate and plant biomass.jpg', plot = plant_trait_plot, width = 88*2, height = 66*2, units = "mm")
}
rm(p1, p2, p3, p4)


# Soil properties ----

soil_data_no_bulk <- soil_data %>% subset(compartment != 'Bulk soil') 

soil_data_rhz <- soil_data %>% subset(compartment == 'Rhizosphere') 
soil_data_hps <- soil_data %>% subset(compartment == 'Hyphosphere')
soil_data_bulk <- soil_data %>% subset(compartment == 'Bulk soil') 

#  ~~~ soil TOC ----

toc1 <- plot_box(soil_data_rhz, group2, `TOC`)+
  ggtitle('Rhizosphere')+
  geom_jitter(alpha = 0.3, width = 0.1)+
  ylim(5, 25)+
  xlab('Treatment')+
  ylab('TOC (%)')+
  annotate('text', y = 25, x = 'Ctrl', label = 'a')+
  annotate('text', y = 25, x = 'Cd', label = 'a')+
  annotate('text', y = 25, x = 'M', label = 'a')+
  annotate('text', y = 25, x = 'M+Cd', label = 'a')
toc1
# ANOVA and Tukey test
agricolae::HSD.test(aov(`TOC` ~ group2, data = soil_data_rhz), # this is aov model
                    "group2", console = TRUE)
# non-parametric test
kruskal.test(`TOC`~group2, data = soil_data_rhz)
pairwise.wilcox.test(soil_data_rhz$`TOC`, soil_data_rhz$group2,
                     p.adjust.method = "BH")


  
toc2 <- plot_box(soil_data_hps, group2, `TOC`)+
  ggtitle('Hyphosphere')+
  geom_jitter(alpha = 0.3, width = 0.1)+
  ylim(5, 25)+
  xlab('Treatment')+
  ylab('TOC (%)')+
  annotate('text', y = 25, x = 'Ctrl', label = 'a')+
  annotate('text', y = 25, x = 'Cd', label = 'b')+
  annotate('text', y = 25, x = 'M', label = 'a')+
  annotate('text', y = 25, x = 'M+Cd', label = 'ab')
toc2
# ANOVA and Tukey test
agricolae::HSD.test(aov(`TOC` ~ group2, data = soil_data_hps), # this is aov model
                    "group2", console = TRUE)
# non-parametric test
kruskal.test(`TOC`~group2, data = soil_data_hps)
pairwise.wilcox.test(soil_data_hps$`TOC`, soil_data_hps$group2,
                     p.adjust.method = "BH")


toc3 <- plot_box(soil_data_bulk, group2, `TOC`)+
  ggtitle('Bulk soil')+
  geom_jitter(alpha = 0.3, width = 0.1)+
  ylim(5, 25)+
  xlab('Treatment')+
  ylab('TOC (%)')+
  annotate('text', y = 25, x = 'Ctrl', label = 'a')+
  annotate('text', y = 25, x = 'Cd', label = 'a')+
  annotate('text', y = 25, x = 'M', label = 'a')+
  annotate('text', y = 25, x = 'M+Cd', label = 'a')
toc3
# ANOVA and Tukey test  
agricolae::HSD.test(aov(`TOC` ~ group2, data = soil_data_bulk), # this is aov model
                    "group2", console = TRUE)
# non-parametric test
kruskal.test(`TOC`~group2, data = soil_data_bulk)
pairwise.wilcox.test(soil_data_bulk$`TOC`, soil_data_bulk$group2,
                     p.adjust.method = "BH")

# ~~~~~~ Fig. S10. Total organic carbon in soils ----
toc_by_trt <- ggarrange(toc1, toc2, toc3,
                                    labels = c("a", "b", "c"),
                                    ncol =3, nrow = 1)
toc_by_trt

if(save){
ggsave('out/Fig. S10. toc_by_trt.jpg', plot = toc_by_trt,
       width = 88*2, height = 66*2, units = "mm")
ggsave('out/Fig. S10. toc_by_trt.pdf', plot = toc_by_trt,
       width = 88*2, height = 66*2, units = "mm")
}
rm(toc1, toc2, toc3)


# ~~~ soil DOC  ----

plot_box(soil_data_rhz, group2, `DOC`)+
  ggtitle('Rhizosphere')

plot_box(soil_data_hps, group2, `DOC`)+
  ggtitle('Hyphosphere')

plot_box(soil_data_bulk, group2, `DOC`)+
  ggtitle('Bulk soil')

# ~~~~~ DOC, TOC percentage calculation ----

soil_data %>% subset(compartment == 'Rhizosphere') %>% 
  plot_box(x = group2, y = DOC)

soil_data %>% subset(compartment == 'Rhizosphere') %>% 
  group_by(group2) %>% 
  summarise(doc_mean=mean(DOC))

(0.74-0.7)/0.7
# = 6%

soil_data %>% subset(compartment == 'Rhizosphere') %>% 
  plot_box(x = group2, y = TOC)

soil_data %>% subset(compartment == 'Rhizosphere') %>% 
  group_by(group2) %>% 
  summarise(toc_mean=mean(TOC))
(15.2-13.9)/13.9
# = 9%

soil_data %>% subset(compartment == 'Hyphosphere') %>% 
  plot_box(x = group2, y = DOC)
soil_data %>% subset(compartment == 'Hyphosphere') %>% 
  group_by(group2) %>% 
  summarise(doc_mean=mean(DOC))


soil_data %>% subset(compartment == 'Hyphosphere') %>% 
  plot_box(x = group2, y = TOC)



# ~~~ soil total Cd ----
plot_box(soil_data_hps, group2, TCd)+geom_jitter()

temp_model <- aov(TCd ~ group2, data = soil_data_hps)

# Perform Tukey's test
TukeyHSD(temp_model)


# ~~~ soil P ----
p1 <- soil_data %>% 
  subset(compartment == 'Rhizosphere') %>% 
  plot_box(x = group2, y = EP)+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)+
  xlab('Treatment')+
  ylab('Extractable P (mg/kg)')+
  ggtitle('Rhizosphere soil')
# tukey

soil_data_rzs <- soil_data %>% 
  subset(compartment == 'Rhizosphere')
anova_rzs_EP <- aov(EP~group2, data = soil_data_rzs)
summary(anova_rzs_EP)
HSD.test(anova_rzs_EP, 'group2', console=TRUE)

# EP groups
# Cd   185.8287      a
# M+Cd 178.1150     ab
# Ctrl 176.6025     ab
# M    169.3288      b

#Kruskal-Wallis test

kruskal.test(EP~group2, data = soil_data_rzs)
pairwise.wilcox.test(soil_data_rzs$EP, soil_data_rzs$group2,
                     p.adjust.method = "BH")
# data:  soil_data_rzs$EP and soil_data_rzs$group2 
# 
# Ctrl   Cd     M     
# Cd   0.0413 -      -     
#   M    0.1347 0.0096 -     
#   M+Cd 0.6353 0.0600 0.0413
# 
# P value adjustment method: BH 

p2 <- soil_data %>% 
  subset(compartment == 'Hyphosphere') %>% 
  plot_box(x = group2, y = EP)+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)+
  xlab('Treatment')+
  ylab('Extractable P (mg/kg)')+
  ggtitle('Hyphosphere soil')
# tukey
soil_data_hps <- soil_data %>% 
  subset(compartment == 'Hyphosphere')
anova_hps_EP <- aov(EP~group2, data = soil_data_hps)
summary(anova_hps_EP)
HSD.test(anova_hps_EP, 'group2', console=TRUE)
# EP groups
# Cd   209.4988      a
# Ctrl 206.5700     ab
# M+Cd 200.5675     ab
# M    199.3950      b


p3 <- soil_data %>% 
  subset(compartment == 'Bulk soil') %>% 
  plot_box(x = group2, y = EP)+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)+
  xlab('Treatment')+
  ylab('Extractable P (mg/kg)')+
  ggtitle('Bulk soil')
# tukey
soil_data_bulk <- soil_data %>% 
  subset(compartment == 'Bulk soil')
anova_bulk_EP <- aov(EP~group2, data = soil_data_bulk)
summary(anova_bulk_EP)
HSD.test(anova_bulk_EP, 'group2', console=TRUE)
# EP groups
# Ctrl 215.9437      a
# M    213.5025     ab
# Cd   207.0600     bc
# M+Cd 204.4263      c



# ggarrange(p1, p2, p3,
#           labels = c('a', 'b', 'c'),
#           nrow = 2, ncol = 2)
# ggsave('out/soil_EP.pdf',
#        w = 88*2, h = 66*2, units = 'mm')
# ggsave('out/soil_EP.jpg',
#        w = 88*2, h = 66*2, units = 'mm')


# Plant traits ----
# ~~~~ plant P ----

p4 <- plant_data %>% 
  plot_box(x = group, y = rt_TP)+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)+
  xlab('Treatment')+
  ylab('Total P (mg/kg)')+
  ggtitle('Root')
# tukey

anova_rt_TP <- aov(rt_TP~group, data = plant_data)
summary(anova_rt_TP)
HSD.test(anova_rt_TP, 'group', console=TRUE)

# rt_TP groups
# Ctrl 1.5050      a
# M    1.4300      a
# M+Cd 1.3325      a
# Cd   1.2875      a

kruskal.test(rt_TP~group, data = plant_data)

# data:  rt_TP by group
# Kruskal-Wallis chi-squared = 0.44118, df = 3, p-value =
#   0.9316


p5 <- plant_data %>% 
  plot_box(x = group, y = sht_TP)+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)+
  xlab('Treatment')+
  ylab('Total P (mg/kg)')+
  ggtitle('Shoot')

# tukey
anova_sht_TP <- aov(sht_TP~group, data = plant_data)
summary(anova_sht_TP)
HSD.test(anova_sht_TP, 'group', console=TRUE)

# sht_TP groups
# M    5.1925      a
# Cd   4.9425      a
# Ctrl 4.7875      a
# M+Cd 4.7400      a

kruskal.test(sht_TP~group, data = plant_data)
# data:  sht_TP by group
# Kruskal-Wallis chi-squared = 0.68382, df = 3, p-value =
#   0.877

p6 <- plant_data %>% 
  plot_box(x = group, y = grain_TP)+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)+
  xlab('Treatment')+
  ylab('Total P (mg/kg)')+
  ggtitle('Grain')
# tukey
anova_grain_TP <- aov(grain_TP~group, data = plant_data)
summary(anova_grain_TP)
HSD.test(anova_grain_TP, 'group', console=TRUE)
# grain_TP groups
# M      3.8275      a
# Cd     2.9350      a
# M+Cd   2.6875      a
# Ctrl   2.2550      a
kruskal.test(grain_TP~group, data = plant_data)
# data:  grain_TP by group
# Kruskal-Wallis chi-squared = 5.8676, df = 3, p-value =
#   0.1182

p7 <- plant_data %>% 
  mutate(sht_P_content = sht_TP*sht_dry_mass, grain_P_content = grain_TP*grain_dry_mass) %>% 
  plot_box(x = group, y = sht_P_content)+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)+
  xlab('Treatment')+
  ylab('Total P (micro gram)')+
  ggtitle('Shoot P content')
# tukey
anova_sht_P_cont <- aov(sht_P_content~group, data = plant_data %>% 
                          mutate(sht_P_content = sht_TP*sht_dry_mass, 
                                 grain_P_content = grain_TP*grain_dry_mass))
summary(anova_sht_P_cont)
HSD.test(anova_sht_P_cont, 'group', console=TRUE)
# sht_P_content groups
# M         210.3773      a
# Ctrl      175.5398     ab
# Cd        171.8195     ab
# M+Cd      159.8545      b
kruskal.test(sht_P_content~group, data = plant_data %>% 
               mutate(sht_P_content = sht_TP*sht_dry_mass, 
                      grain_P_content = grain_TP*grain_dry_mass))
# data:  sht_P_content by group
# Kruskal-Wallis chi-squared = 7.8309, df = 3, p-value =
#   0.04964*
plant_data_add_P_cont <- plant_data %>% 
  mutate(sht_P_content = sht_TP*sht_dry_mass, 
         grain_P_content = grain_TP*grain_dry_mass)
pairwise.wilcox.test(plant_data_add_P_cont$sht_P_content,
                     plant_data_add_P_cont$group,
                     p.adjust.method = "BH")
# data:  plant_data_add_P_cont$sht_P_content and plant_data_add_P_cont$group 
# 
# Ctrl Cd   M   
# Cd   0.89 -    -   
#   M    0.23 0.17 -   
#   M+Cd 0.41 0.41 0.17
# 
# P value adjustment method: BH 

p8 <- plant_data %>% 
  mutate(sht_P_content = sht_TP*sht_dry_mass, grain_P_content = grain_TP*grain_dry_mass) %>% 
  plot_box(x = group, y = grain_P_content)+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)+
  xlab('Treatment')+
  ylab('Total P (micro gram)')+
  ggtitle('Grain P content')
# tukey
anova_grain_P_cont <- aov(grain_P_content~group, data = plant_data %>% 
                            mutate(sht_P_content = sht_TP*sht_dry_mass, 
                                   grain_P_content = grain_TP*grain_dry_mass))
summary(anova_grain_P_cont)
HSD.test(anova_grain_P_cont, 'group', console=TRUE)
# grain_P_content groups
# M           55.72222      a
# M+Cd        31.02708      a
# Cd          26.63030      a
# Ctrl        23.48920      a
kruskal.test(grain_P_content~group, data = plant_data %>% 
               mutate(sht_P_content = sht_TP*sht_dry_mass, 
                      grain_P_content = grain_TP*grain_dry_mass))
# Kruskal-Wallis chi-squared = 8.2721, df = 3, p-value =
#   0.04071
pairwise.wilcox.test(plant_data_add_P_cont$grain_P_content,
                     plant_data_add_P_cont$group,
                     p.adjust.method = "BH")
# data:  plant_data_add_P_cont$grain_P_content and plant_data_add_P_cont$group 
# 
# Ctrl Cd   M   
# Cd   0.49 -    -   
#   M    0.17 0.17 -   
#   M+Cd 0.30 0.49 0.23
# 
# P value adjustment method: BH 

# ~~~~ extractable and total P in soil and plant ----
ggarrange(p1, p2, p3,
          p4, p5, p6,         
          labels = c('a', 'b', 'c',
                     'd', 'e', 'f'),
          nrow = 2, ncol = 3)
ggsave('out/soil_plant_P.pdf',
       w = 88*2, h = 66*3, units = 'mm')
ggsave('out/soil_plant_P.jpg',
       w = 88*2, h = 66*3, units = 'mm')


# ~~~~ Shoot and grain P content (unit: mg, not conc.)
ggarrange(p7, p8,         
          labels = c('a', 'b'),
          nrow = 1, ncol = 2)
ggsave('out/plant_P_content.pdf',
       w = 88*2, h = 66*1, units = 'mm')
ggsave('out/plant_P_content.jpg',
       w = 88*2, h = 66*1, units = 'mm')

rm(p1, p2, p3, p4, p5, p6, p7, p8)



# N conc in soil and plant ----
# ~~~ soil N ----

soil_data %>% subset(compartment == 'Rhizosphere') %>% 
  plot_box(x = group2, y = TN)+
  ggtitle('Rhizosphere soil')+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)
soil_data %>% subset(compartment == 'Hyphosphere') %>% 
  plot_box(x = group2, y = TN)+
  ggtitle('Hyphosphere soil')+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)
soil_data %>% subset(compartment == 'Bulk soil') %>% 
  plot_box(x = group2, y = TN)+
  ggtitle('Bulk soil')+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)

soil_data %>% subset(compartment == 'Rhizosphere') %>% 
  plot_box(x = group2, y = amm_N)+
  ggtitle('Rhizosphere soil')+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)
soil_data %>% subset(compartment == 'Hyphosphere') %>% 
  plot_box(x = group2, y = amm_N)+
  ggtitle('Hyphosphere soil')+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)
soil_data %>% subset(compartment == 'Bulk soil') %>% 
  plot_box(x = group2, y = amm_N)+
  ggtitle('Bulk soil')+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)

soil_data %>% subset(compartment == 'Rhizosphere') %>% 
  plot_box(x = group2, y = nitrate)+
  ggtitle('Rhizosphere soil')+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)
soil_data %>% subset(compartment == 'Hyphosphere') %>% 
  plot_box(x = group2, y = nitrate)+
  ggtitle('Hyphosphere soil')+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)
soil_data %>% subset(compartment == 'Bulk soil') %>% 
  plot_box(x = group2, y = nitrate)+
  ggtitle('Bulk soil')+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)

# ~~~ plant N ----
plant_data_add_P_cont %>% 
  plot_box(x = group, y = rt_TN)+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)
plant_data_add_P_cont %>% 
  plot_box(x = group, y = sht_TN)+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)
plant_data_add_P_cont %>% 
  plot_box(x = group, y = grain_TN)+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)

plant_data_add_P_cont %>% 
  plot_box(x = group, y = Gs)+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)

plot_box(plant_data, group, y = apop_cd/(apop_cd+symp_cd))+
  geom_jitter(alpha = 0.3, size = 2, width = 0.1)




# soil NH4+
soil_data %>% subset(compartment == 'Rhizosphere') %>% 
  plot_box(x = group2, y = amm_N)
soil_data %>% subset(compartment == 'Hyphosphere') %>% 
  plot_box(x = group2, y = amm_N)

soil_data %>% subset(compartment == 'Rhizosphere') %>% 
  group_by(group2) %>% 
  summarise(mean_amm_N=mean(amm_N))

# percentage increase -- amm_N in rhizosphere
(10.2-7.54)/7.54 
# = 35%

soil_data %>% subset(compartment == 'Hyphosphere') %>% 
  group_by(group2) %>% 
  summarise(mean_amm_N=mean(amm_N))
# percentage increase -- amm_N in hyphosphere
(10.3-8.92)/8.92 
# = 15%
