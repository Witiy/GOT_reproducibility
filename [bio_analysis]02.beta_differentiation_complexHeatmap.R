library(ComplexHeatmap)
library(viridis)
library(circlize)
library(grid)
library(RColorBrewer)

# 读取数据
A <- read.csv('/disk/xuruihong/tmp.csv', row.names = 1)
A <- as.matrix(A)

meta <- read.csv('/disk/xuruihong/tmp2.csv', row.names = 1)
samples <- meta$prog_type
pseudotime <- meta$pseudotime  # 伪时间信息

# 生成 Spectral colormap（用于热图的值）
spectral_col_fun <- colorRamp2(
  seq(min(A), max(A), length.out = 11),  
  rev(brewer.pal(11, "Spectral"))  
)

# 颜色映射定义（细胞状态）
cell_state_colors <- c(
  'endo' = '#d62728', 
  'exo' = '#1f77b4', 
  'prog_endo' = '#deb887', 
  'prog_exo' = '#95d0fc'
)

# 生成 viridis colormap 的映射函数（用于 pseudotime）
pseudotime_col_fun <- colorRamp2(
  seq(min(pseudotime), max(pseudotime), length.out = 20),  
  viridis(20)  
)



top_annotation <- HeatmapAnnotation(
  CellState = samples,
  Pseudotime = pseudotime,
  col = list(CellState=cell_state_colors, Pseudotime=pseudotime_col_fun),
  #show_annotation_name = TRUE,  # 启用注释名称
  #annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),  # 调整字体
  height = unit(2.5, "mm")  # 增加高度以便显示
  )

row_split_no = 50
col_split_no = 1236
# 自定义行分割标签
row_split_groups <- c(rep("Prog Exo DEG", row_split_no), rep("Prog Endo DEG", nrow(A)-row_split_no))  # 可改成你的分组名称

row_colors <- c(
  'Prog Exo DEG' = "#95d0fc", 
  'Prog Endo DEG' = "#deb887"
)

row_group_anno <- rowAnnotation(
  Group = row_split_groups,
  col=list(Group=row_colors),
  width = unit(4, "mm"),
  show_legend = FALSE
)

# 需要标注的基因列表
genes <- c(
  'CALB1',
  'GPC3',
  'ANXA2',
  'PTPN13',
  'SOX2',
  'PDLIM1',
  'LIN28A',
  'SEMA3C',
  'PAM',
  'COL9A3',
  'GATM',
  'ALDH1A1',
  'LMO3',
  'TM4SF4',
  'DLL1',
  'PTF1A',
  'ADD3',
  'CPA2',
  'F3',
  'SPOCK1',
  'INSM1'
)

# 创建行标注 (基因名称)
row_anno <- rowAnnotation(
  link = anno_mark(
    at = which(rownames(A) %in% genes), 
    labels = rownames(A)[rownames(A) %in% genes], 
    labels_gp = gpar(fontsize = 10)
  )
)
?Heatmap

# 生成热图对象
ht <- Heatmap(A,
              name = "Scaled Exp",
              col = spectral_col_fun, 
              #row_title_gp = gpar(fontsize = 4),
              show_row_names = F, 
              show_column_names = FALSE,
              cluster_rows = T,   # 允许行聚类
              cluster_columns = FALSE,
              row_split = row_split_groups,  # 自定义行拆分
              top_annotation = top_annotation,  
              right_annotation = row_anno,  # 右侧基因标注
              left_annotation = row_group_anno,  # 左侧分组注释
              column_title = c('To Exo Lineage', 'To Endo Lineage'),
              
              
              column_split = c(rep(1, col_split_no), rep(2, ncol(A) - col_split_no)),  # 仍然按列拆分
              row_gap = unit(1.5, "mm"),
              column_gap = unit(1.5, "mm")
             
)

# 保存绘制的热图
pdf(file="~/Fig5_heatmap.pdf", width=10, height=8)
draw(ht , heatmap_legend_side = "right", annotation_legend_side = "left")
dev.off()







edgelist_grn = read.table('/storage/xuruihong/repository/pygot_data/beta_diff/endo_exo_grn.txt', sep='\t', header = T)
plot_new_grn <- function (edgelist_grn, show_labels = "tophubs", top_n_hubs = 5, 
                          interactive = FALSE, layout = igraph::with_kk, arrow.gap = 0.01, 
                          ranked = TRUE, dim_interactive = c(600, 600)) 
{
  requireNamespace("intergraph", quietly = TRUE)
  h <- get_hubs_grn(edgelist_grn, return_degree = TRUE, ranked = ranked)
  geneIDs <- unique(c(as.character(edgelist_grn[, 1]), as.character(edgelist_grn[, 
                                                                                 2])))
  nod_at <- data.frame(Gene = geneIDs, stringsAsFactors = FALSE)
  nod_at$Class <- ifelse(nod_at$Gene %in% edgelist_grn[, 1], 
                         "Regulator", "Target")
  nod_at$Degree <- h$Degree$Degree[h$Degree$Gene %in% nod_at$Gene]
  nod_at$isHub <- ifelse(nod_at$Gene %in% h$Hubs$Gene, TRUE, 
                         FALSE)
  nod_at <- nod_at[order(nod_at$Class, -nod_at$Degree), ]
  if (interactive) {
    graph <- igraph::graph_from_data_frame(d = edgelist_grn, 
                                           vertices = nod_at, directed = TRUE)
    graph <- igraph::simplify(graph)
    graph_d3 <- networkD3::igraph_to_networkD3(graph, group = nod_at$Class)
    graph_d3$nodes <- merge(graph_d3$nodes, nod_at, by.x = "name", 
                            by.y = "Gene", sort = FALSE)
    my_color <- "d3.scaleOrdinal() .domain([\"Regulator\", \"Target\"]) .range([\"forestgreen\", \"orange\"])"
    d <- dim_interactive
    p <- networkD3::forceNetwork(Links = graph_d3$links, 
                                 Nodes = graph_d3$nodes, Source = "source", Target = "target", 
                                 NodeID = "name", Group = "group", colourScale = my_color, 
                                 Nodesize = "Degree", height = d[2], width = d[1], 
                                 opacity = 1, zoom = TRUE, fontSize = 20, legend = TRUE)
  }
  else {
    add_edges <- ggnetwork::geom_edges(color = "grey60", 
                                       alpha = 0.25, arrow = ggplot2::arrow(length = ggnetwork::unit(0.1, 
                                                                                                  "lines"), type = "closed"), curvature = 0.1, 
                                       show.legend = FALSE)
    
    
    if (show_labels == "all") {
      nod_at$Degree2 <- nod_at$Degree * 0.4
      add_nodelabel <- ggnetwork::geom_nodetext(ggplot2::aes_(label = ~name, 
                                                              size = ~Degree2), show.legend = FALSE)
    }
    else if (show_labels == "allhubs") {
      add_nodelabel <- ggnetwork::geom_nodelabel_repel(ggplot2::aes_(label = ~name), 
                                                       color = "azure4", box.padding = ggnetwork::unit(1, 
                                                                                                       "lines"), data = function(x) {
                                                                                                         x[x$isHub, ]
                                                                                                       }, show.legend = FALSE, max.overlaps = Inf)
    }
    else if (show_labels == "tophubs") {
      tophubs <- nod_at[nod_at$isHub == TRUE, 1][seq_len(top_n_hubs)]
      nod_at$isTopHub <- ifelse(nod_at$Gene %in% tophubs, 
                                TRUE, FALSE)
      add_nodelabel <- ggnetwork::geom_nodelabel_repel(ggplot2::aes_(label = ~name), 
                                                       color = "azure4", box.padding = ggnetwork::unit(1, 
                                                                                                       "lines"), data = function(x) {
                                                                                                         x[x$isTopHub, ]
                                                                                                       }, show.legend = FALSE, max.overlaps = Inf)
    }
    else if (show_labels == "none") {
      add_nodelabel <- NULL
    }
    else {
      stop("Please, specify a valid option for 'show_labels'.")
    }
    graph <- igraph::graph_from_data_frame(d = edgelist_grn, 
                                           vertices = nod_at, directed = TRUE)
    graph <- igraph::simplify(graph)
    n <- ggnetwork::ggnetwork(graph, layout = layout(), 
                              arrow.gap = arrow.gap)
    n$Class <- factor(n$Class, levels = c("Target", "Regulator"))
    p <- ggplot2::ggplot(n, ggplot2::aes_(x = ~x, y = ~y, 
                                          xend = ~xend, yend = ~yend)) + add_edges + ggnewscale::new_scale_color() + 
      ggnetwork::geom_nodes(ggplot2::aes_(color = ~Class, 
                                          size = ~Degree)) + ggplot2::scale_color_manual(values = c("orange", 
                                                                                                                    "darkgreen")) + add_nodelabel + ggnetwork::theme_blank()
  }
  return(p)
}


#TODO 
#用网络聚类

p<-plot_new_grn(edgelist_grn[,1:2], 
            interactive = F,
         layout = igraph::with_kk,
         ranked = T, top_n_hubs = 18, arrow.gap = 0.0001)

p

ggsave("~/Fig5_grn.pdf", p, width = 10,
       height = 10)









# 读取数据
A <- read.csv('/storage/xuruihong/repository/pygot_data/beta_diff/endo_exo_grn_binary.txt', row.names = 1, sep='\t')
A <- as.matrix(A)

meta <- read.csv('/storage/xuruihong/repository/pygot_data/beta_diff/endo_exo_grn_binary_meta.txt', row.names = 1, sep='\t')


# 颜色映射定义（细胞状态）
cell_state_colors <- c(
  'endo' = '#d62728', 
  'exo' = '#1f77b4',
  'ns' = 'grey'
  
)

# 生成 viridis colormap 的映射函数（用于 pseudotime）
col_fun1 <- colorRamp2(
  c(min(meta$Prog.stage),0, max(meta$Prog.stage)),  
  c('#deb887','white', '#95d0fc')  
)
#viridis(20)  
col_fun2 <- colorRamp2(
  c(min(meta$Mature.stage),0, max(meta$Mature.stage)),  
  c('#d62728','white', '#1f77b4')  
)
top_annotation <- HeatmapAnnotation(
  ProgDEG = meta$Prog.stage,
  MatureDEG =meta$Mature.stage,
  col = list(ProgDEG=col_fun1, MatureDEG=col_fun2),
  #show_annotation_name = TRUE,  # 启用注释名称
  #annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),  # 调整字体
  height = unit(2.5, "mm")  # 增加高度以便显示
)


# 生成热图对象
ht <- Heatmap(A,
              name = "GRN",
              col = brewer.pal(9, "Greys"),
              #row_title_gp = gpar(fontsize = 4),
              show_row_names = T, 
              show_column_names = FALSE,
              cluster_rows = T,   # 允许行聚类
              cluster_columns = T,
              top_annotation = top_annotation,  

              
)

pdf(file="~/Fig5_grn_heatmap.pdf", width=10, height=8)
draw(ht)
dev.off()





