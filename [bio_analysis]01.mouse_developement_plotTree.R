# libraries
library(ggraph)
library(igraph)
library(tidyverse)
library(tidygraph)
#relations <- read.csv('/storage/xuruihong/repository/pygot_data/development/road_map_edge_E0.5-E8.5.csv')
#actors <- read.csv('/storage/xuruihong/repository/pygot_data/development/road_map_node_E0.5-E8.5.csv')

#relations <- read.csv('/storage/xuruihong/repository/pygot_data/development/road_map_edge_E8.5-E13.5.csv')
#actors <- read.csv('/storage/xuruihong/repository/pygot_data/development/road_map_node_E8.5-E13.5.csv')

relations <- read.csv('/storage/xuruihong/repository/pygot_data/development/road_map_edge_E0.5-E13.5_sde.csv')
actors <- read.csv('/storage/xuruihong/repository/pygot_data/development/road_map_node_E0.5-E13.5_sde.csv')

actors
# 1. 筛选出每个 position 中 y 最大的节点作为叶子节点
leaf_nodes <- actors %>%
  group_by(position) %>%
  filter(stage == min(stage)) %>%
  pull(name)

actors %>%
  group_by(position)

actors

actors$show_label <- ifelse(actors$name %in% leaf_nodes, actors$cell_state, NA)

graph <- graph_from_data_frame(d = relations, 
                               directed=TRUE, 
                               vertices=actors)

layout <- create_layout(graph, layout = 'tree')

row.names(actors) <- actors$name
actors[row.names(actors)[1:5],]$position <- 94
#layout$x = -(actors[layout$name,'position'] * 2)
layout$y = -(actors[layout$name,'position'] * 2)
layout$x = (actors[layout$name,'x'] * 5)
#options(ggrepel.max.overlaps = Inf)






p1<-ggraph(graph, layout = layout)  +
  
  geom_edge_diagonal(mapping = aes(edge_color = weight),
                     flipped=T,
                 ) +
  geom_node_point(size = 1.5,
                  mapping = aes(colour = group),
                  ) +
  
  #geom_node_text(mapping = aes(label = show_label), repel = TRUE, size = 1, 
  #               nudge_y = -1,
  #               angle = 0, max.overlaps=Inf) +
  #scale_edge_width_continuous(range = c(0.,1.0)) +
  scale_edge_color_continuous(low = "grey",high = "black")+
  scale_color_manual(values = c(
    'Other'="#000000", 
    'Neuroectoderm'="#145A32",
    'Surface ectoderm'="#78281F",
    'Endoderm'="#9A7D0A",
    'Mesoderm'="#2471A3",
    'ExE embryo'="#E67E22")) +
   theme_graph() + theme(
     legend.key.size = unit(0.5, "cm"),
     legend.text = element_text(family = "Arial", size=8),  # 设置图例文本的字体为 Arial
     legend.title = element_text(family = "Arial", size=10)  # 设置图例标题的字体为 Arial
   )
p1
install.packages('showtext')
library(showtext)
font_add('Arial','/disk/share/xuruihong/font/Arial.ttf') #加载字体，MAC 中字体库在 /Library/Fonts
showtext_auto() #自动调用showtext，否则无法在ggsave()中使用，因为ggsave会自动打开和关闭图形设备。


#ggsave(p1, file='/disk/share/xuruihong/pygot_fig/development_tree.pdf', width=9, height=7)
