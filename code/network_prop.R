# plot number of nodes and node of edges for bacterial co-occurrence network ####

# load data
nw_prop_data <- read_xlsx('data/network_node_degree_and_edge.xlsx')

nw_prop_data$trt <- factor(nw_prop_data$trt, levels=c('Ctrl', 'Cd', 'M', 'M+Cd'))

nw_prop_data$no_node <- as.numeric(nw_prop_data$no_node)
nw_prop_data$no_edge <- as.numeric(nw_prop_data$no_edge)
nw_prop_data$pp <- as.numeric(nw_prop_data$pp)
nw_prop_data$np <- as.numeric(nw_prop_data$np)


# ~~ numbers of nodes and edges plots ------------------------------------------
colnames(nw_prop_data)
tp1 <- nw_prop_data %>% 
  subset(compartment == 'rzs') %>% 
  ggplot(aes(x= trt, y = no_node, group=1))+ 
  # add group = 1 here is needed to create line connecting the points
  geom_line(color = 'grey', size = 1.5)+
  geom_point(shape=21, color="black", fill="gray", size=5)+
  ggtitle('Rhizosphere')+
  xlab('Treatment')+
  ylab('Number of nodes')

tp2 <- nw_prop_data %>% 
  subset(compartment == 'hps') %>% 
  ggplot(aes(x= trt, y = no_node, group=1))+ 
  # add group = 1 here is needed to create line connecting the points
  geom_line(color = 'grey', size = 1.5)+
  geom_point(shape=21, color="black", fill="gray", size=5)+
  ggtitle('Hyphosphere')+
  xlab('Treatment')+
  ylab('Number of nodes')
  
  
tp3 <- nw_prop_data %>% 
  subset(compartment == 'hph') %>% 
  ggplot(aes(x= trt, y = no_node, group=1))+ 
  # add group = 1 here is needed to create line connecting the points
  geom_line(color = 'grey', size = 1.5)+
  geom_point(shape=21, color="black", fill="gray", size=5)+
  ggtitle('Hyphae-associated')+
  xlab('Treatment')+
  ylab('Number of nodes')



tp4 <- nw_prop_data %>% 
  subset(compartment == 'rzs') %>% 
  ggplot(aes(x= trt, y = no_edge, group=1))+ 
  # add group = 1 here is needed to create line connecting the points
  geom_line(color = 'grey', size = 1.5)+
  geom_point(shape=21, color="black", fill="gray", size=5)+
  ggtitle('Rhizosphere')+
  xlab('Treatment')+
  ylab('Number of edges')

tp5 <- nw_prop_data %>% 
  subset(compartment == 'hps') %>% 
  ggplot(aes(x= trt, y = no_edge, group=1))+ 
  # add group = 1 here is needed to create line connecting the points
  geom_line(color = 'grey', size = 1.5)+
  geom_point(shape=21, color="black", fill="gray", size=5)+
  ggtitle('Hyphosphere')+
  xlab('Treatment')+
  ylab('Number of edges')

tp6 <- nw_prop_data %>% 
  subset(compartment == 'hph') %>% 
  ggplot(aes(x= trt, y = no_edge, group=1))+ 
  # add group = 1 here is needed to create line connecting the points
  geom_line(color = 'grey', size = 1.5)+
  geom_point(shape=21, color="black", fill="gray", size=5)+
  ggtitle('Hyphae-associated')+
  xlab('Treatment')+
  ylab('Number of edges')


# ~~ bacterial network - positive and negative correlation percentages ----------------
nw_pos_neg_perc <- read_xlsx('data/network_positive_negative_percentage.xlsx')

nw_pos_neg_perc$trt <- factor(nw_pos_neg_perc$trt, levels=c('Ctrl', 'Cd', 'M', 'M+Cd'))

# bar chart
t1 <- nw_pos_neg_perc %>% 
  subset(compartment == 'rzs') %>% 
  ggplot(aes(x= trt, y = perc, fill = cor_type))+
  geom_bar(stat = 'identity', width = 0.5)+
  ggtitle('Rhizosphere')+
  xlab('Treatment')+
  ylab('Percentage (%)')+
  theme(legend.position = 'none')

t2 <- nw_pos_neg_perc %>% 
  subset(compartment == 'hps') %>% 
  ggplot(aes(x= trt, y = perc, fill = cor_type))+
  geom_bar(stat = 'identity', width = 0.5)+
  ggtitle('Hyphosphere')+
  xlab('Treatment')+
  ylab('Percentage (%)')+
  theme(legend.position = 'none')

t3 <- nw_pos_neg_perc %>% 
  subset(compartment == 'hph') %>% 
  ggplot(aes(x= trt, y = perc, fill = cor_type))+
  geom_bar(stat = 'identity', width = 0.5)+
  ggtitle('Hyphae-associated')+
  xlab('Treatment')+
  ylab('Percentage (%)')+
  theme(legend.position = 'none')

t4 <- nw_pos_neg_perc %>% 
  subset(compartment == 'bs') %>% 
  ggplot(aes(x= trt, y = perc, fill = cor_type))+
  geom_bar(stat = 'identity', width = 0.5)+
  ggtitle('Bulk soil')+
  xlab('Treatment')+
  ylab('Percentage (%)')+
  theme(legend.position = 'none')



# ~~~~ Fig. S11. Network prop of bact co-occur nw and pos/neg by compt. ----
ggarrange(tp1, tp2, tp3, 
          tp4, tp5, tp6,
          t1, t2, t3,
          labels = c('a', 'b', 'c', 
                     'd', 'e', 'f',
                     'g', 'h','i'),
          ncol = 3, nrow = 3)


if(save){
  ggsave('out/Fig. S11. no_nodes_edges_pos.neg.pdf',
       width = 88*3, height = 66*3, units = "mm")

  ggsave('out/Fig. S11. no_nodes_edges_pos.neg.jpg',
       width = 88*3, height = 66*3, units = "mm")
}


# ~~~~ Fig. 8. Proportion of facilitation of B-B inter by compt ----
pos.neg <- ggarrange(t4, t1, 
          t2, t3,
          labels = c('a', 'b', 'c', 
                     'd'),
          ncol = 2, nrow = 2)


if(save){
  ggsave('out/Fig. 8. perc.pos.neg.pdf', plot = pos.neg,
       width = 88*2, height = 66*2, units = "mm")

ggsave('out/Fig. 8. perc.pos.neg.jpg', plot = pos.neg,
       width = 88*2, height = 66*2, units = "mm")
}

rm(tp1, tp2, tp3, tp4, tp5, tp6)
rm(t1, t2, t3, t4)

# point and line chart
# t1 <- nw_pos_neg_perc %>% 
#   subset(compartment == 'bs' & cor_type == 'Positive') %>% 
#   ggplot(aes(x= trt, y = perc, group = 1))+
#   geom_line()+
#   geom_point(size = 3, color = 'darkgrey')+
#   ggtitle('Bulk soil')+
#   xlab('Treatment')+
#   ylab('Facilitation (%)')+
#   ylim(40, 75)
# 
# t2 <- nw_pos_neg_perc %>% 
#   subset(compartment == 'rzs' & cor_type == 'Positive') %>% 
#   ggplot(aes(x= trt, y = perc, group = 1))+
#   geom_line()+
#   geom_point(size = 3, color = 'darkgrey')+
#   ggtitle('Rhizosphere')+
#   xlab('Treatment')+
#   ylab('Facilitation (%)')+
#   ylim(40, 75)
# 
# t3 <- nw_pos_neg_perc %>% 
#   subset(compartment == 'hps' & cor_type == 'Positive') %>% 
#   ggplot(aes(x= trt, y = perc, group = 1))+
#   geom_line()+
#   geom_point(size = 3, color = 'darkgrey')+
#   ggtitle('Hyphosphere')+
#   xlab('Treatment')+
#   ylab('Facilitation (%)')+
#   ylim(40, 75)
# 
# t4 <- nw_pos_neg_perc %>% 
#   subset(compartment == 'hph' & cor_type == 'Positive') %>% 
#   ggplot(aes(x= trt, y = perc, group = 1))+
#   geom_line()+
#   geom_point(size = 3, color = 'darkgrey')+
#   ggtitle('Hyphae-associated')+
#   xlab('Treatment')+
#   ylab('Facilitation (%)')+
#   ylim(40, 75)
# 
# 
# 
# ggarrange(t1, t2, t3, t4,
#           labels = c('a', 'b', 'c', 'd'),
#           ncol = 2, nrow = 2)
# 
# ggsave('out/cor_type_pec_3.pdf',
#        width = 88*1.5, height = 66*2, units = "mm")
# 
# rm(t1, t2, t3, t4)




# plot node and edge numbers - metabo-ASV network ----
# ~~ plot node and edge numbers ----
# load data
nw_prop_data_metabo_COR_asv <- read_xlsx('data/metabo_ASV_network_no_node_edge_np_pp.xlsx')

nw_prop_data_metabo_COR_asv$trt <- factor(nw_prop_data_metabo_COR_asv$trt, levels=c('Ctrl', 'Cd', 'M', 'M+Cd'))

nw_prop_data_metabo_COR_asv$no_node <- as.numeric(nw_prop_data_metabo_COR_asv$no_node)
nw_prop_data_metabo_COR_asv$no_edge <- as.numeric(nw_prop_data_metabo_COR_asv$no_edge)
nw_prop_data_metabo_COR_asv$pp <- as.numeric(nw_prop_data_metabo_COR_asv$pp)
nw_prop_data_metabo_COR_asv$np <- as.numeric(nw_prop_data_metabo_COR_asv$np)

# plot node and edge numbers 

p1 <- nw_prop_data_metabo_COR_asv %>% 
  subset(compartment == 'rzs') %>% 
  ggplot(aes(x= trt, y = no_node, group=1))+ 
  # add group = 1 here is needed to create line connecting the points
  geom_line(color = 'grey', size = 1.5)+
  geom_point(shape=21, color="black", fill="gray", size=5)+
  ggtitle('Rhizosphere')+
  xlab('Treatment')+
  ylab('Number of nodes')

p2 <- nw_prop_data_metabo_COR_asv %>% 
  subset(compartment == 'hps') %>% 
  ggplot(aes(x= trt, y = no_node, group=1))+ 
  # add group = 1 here is needed to create line connecting the points
  geom_line(color = 'grey', size = 1.5)+
  geom_point(shape=21, color="black", fill="gray", size=5)+
  ggtitle('Hyphosphere')+
  xlab('Treatment')+
  ylab('Number of nodes')

p3 <- nw_prop_data_metabo_COR_asv %>% 
  subset(compartment == 'rzs') %>% 
  ggplot(aes(x= trt, y = no_edge, group=1))+ 
  # add group = 1 here is needed to create line connecting the points
  geom_line(color = 'grey', size = 1.5)+
  geom_point(shape=21, color="black", fill="gray", size=5)+
  ggtitle('Rhizosphere')+
  xlab('Treatment')+
  ylab('Number of edges')

p4 <- nw_prop_data_metabo_COR_asv %>% 
  subset(compartment == 'hps') %>% 
  ggplot(aes(x= trt, y = no_edge, group=1))+ 
  # add group = 1 here is needed to create line connecting the points
  geom_line(color = 'grey', size = 1.5)+
  geom_point(shape=21, color="black", fill="gray", size=5)+
  ggtitle('Hyphosphere')+
  xlab('Treatment')+
  ylab('Number of edges')


# ~~ plot positive and negative pecentagge ----
# load data 
metabo_asv_pos_neg_perc <- read_xlsx('data/metabo_COR_ASV_positive_negative_perc.xlsx')
metabo_asv_pos_neg_perc$trt <- factor(metabo_asv_pos_neg_perc$trt, 
                                      levels=c('Ctrl', 'Cd', 'M', 'M+Cd'))

t1 <- metabo_asv_pos_neg_perc %>% 
  subset(compartment == 'rzs') %>% 
  ggplot(aes(x= trt, y = perc, fill = cor_type))+
  geom_bar(stat = 'identity', width = 0.5)+
  ggtitle('Rhizosphere')+
  xlab('Treatment')+
  ylab('Percentage (%)')+
  theme(legend.position = 'none')

t2 <- metabo_asv_pos_neg_perc %>% 
  subset(compartment == 'hps') %>% 
  ggplot(aes(x= trt, y = perc, fill = cor_type))+
  geom_bar(stat = 'identity', width = 0.5)+
  ggtitle('Hyphosphere')+
  xlab('Treatment')+
  ylab('Percentage (%)')+
  theme(legend.position = 'none')


# ~~~~ Fig. S21. Network properties of the metabolome-ASV correlation networks of the rhizosphere and the hyphosphere compartments. ----
ggarrange(p1, p2, 
          p3, p4,
          t1, t2,
          labels = c('a', 'b', 'c', 'd', 'e', 'f'),
          ncol = 2, nrow = 3)

if(save){
  ggsave('out/Fig. S21. no_nodes_and_edges_metabo_COR_asv_pos.neg.pdf',
       width = 88*2, height = 66*3, units = "mm")

ggsave('out/Fig. S21. no_nodes_and_edges_metabo_COR_asv_pos.neg.jpg',
       width = 88*2, height = 66*3, units = "mm")
}

rm(p1, p2, p3, p4)
rm(t1, t2)
