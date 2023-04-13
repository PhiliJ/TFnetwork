#' @title Contruct Network for plotting
#' @param mouse use mouse genome or not
#' @return net: TF-Target network
#' @export
#' @import magrittr
#' @import dplyr
Cons_net <- function(mouse = FALSE) {
  if (mouse) {
    data("mmDB")
    result <- mmDB %>% filter(TF %in% tfs)
    net <- result %>% filter(Target %in% DErslt$degs)
  } else {
    data("hsDB")
    result <- hsDB %>% filter(TF %in% tfs)
    net <- result %>% filter(Target %in% DErslt$degs)
  }
  assign("net", net, envir = .GlobalEnv)
  write.table(net, file = "TF_Network.txt", sep = "\t", row.names = FALSE)
}

#' @title Plot TF-Target network
#' @param nodeshape shape of the nodes
#' @return TFnet.pdf: network plot
#' @export
#' @import igraph
#' @import magrittr
#' @import dplyr
#' @import ggplot2
#' @import circlize
plot_tfnetwork <- function(net,
                           nodeshape = "circle",
                           TF_color = "#F0E442",
                           Up_color = "#fb8072",
                           Down_color = "#80b1d3",
                           TF_alpha = 0.75,
                           Up_alpha = 0.3,
                           Down_alpha = 0.3,
                           TF_size = 10,
                           Target_size = 8,
                           TF_labelsize = 0.7,
                           Target_labelsize = 0.7,
                           label_color = "black",
                           edge_color = "grey",
                           edge_alpha = 0.75,
                           netlayout = NULL,
                           pdf.name = "TFnet.pdf",
                           pdf.height = 15,
                           pdf.width = 15,
                           edge.arrow.size = 0.4,
                           vertex.label.family = "sans") {
  all_nodes <- unique(c(net$TF, net$Target))
  
  nodecolor <- ifelse(all_nodes %in% net$TF,
                      adjustcolor(TF_color, alpha.f = TF_alpha),
                      ifelse(all_nodes %in% DErslt$up,
                             adjustcolor(Up_color, alpha.f = Up_alpha),
                             adjustcolor(Down_color, alpha.f = Down_alpha)))
  
  nodesize <- ifelse(all_nodes %in% net$TF,
                     scale(table(net$TF),center = F) * TF_size + 6,
                     Target_size)
  labelsize <- ifelse(all_nodes %in% net$TF,
                      scale(table(net$TF),center = F) * TF_labelsize + 0.6,
                      Target_labelsize)
  
  edgecolor <- adjustcolor(edge_color, alpha.f = edge_alpha)
  
  node <- data.frame(nodes = all_nodes,
                     nodesize =  nodesize,
                     nodecolor = nodecolor)
  
  
  edgecolor <- rep(edgecolor, nrow(net))
  
  #定义network
  link <- net
  colnames(link) <- c("From", "To")
  link$color <- edgecolor
  
  network=graph_from_data_frame(d=link, vertices=node, directed=T)
  
  V(network)$color <- V(network)$nodecolor
  V(network)$size <- V(network)$nodesize 
  V(network)$label.cex <- labelsize
  V(network)$shape <- nodeshape
  V(network)$label.color <- label_color 
  
  edge.start <- ends(network, es = E(network), names = FALSE)[, 1] 
  edge.col <- V(network)$color[edge.start]
  E(network)$color <- E(network)$color
  
  #设置网络图的布局
  if (is.null(netlayout)) {
    netlayout <- layout_(network, 
                         with_kk(kkconst = max(vcount(network) + 50, 1)))
  }
  
  #将网络对象存储在内存中
  assign("network", network, envir = .GlobalEnv)
  
  
  #保存为 PDF 文件
  pdf(file = pdf.name, height = pdf.height, width = pdf.width)
  plot(network, layout = netlayout, edge.arrow.size = edge.arrow.size,
       vertex.label.family = vertex.label.family, res = c(1200, 1200))
  dev.off()
  
}


#' @title Plot highlighted TF network
#' @param highlight_node The key node you want to highlight
#' @return TF_highlight.pdf: highlighted network plot
#' @export
#' @import igraph
#' @import magrittr
#' @import dplyr
#' @import ggplot2
#' @import circlize
plot_highlighted_tfnetwork <- function(network, 
                                       highlight_node, 
                                       highlight_TF_color = "#f19d00", 
                                       #                                     highlight_Target_color = "#0072B2", 
                                       highlight_edge_color = "#e41a1c",
                                       netlayout = NULL,
                                       pdf_height = 15, 
                                       pdf_width = 15, 
                                       vertex_label_family = "sans") {
  # 判断 highlight_node 是否为 TF 列中的节点
  if (highlight_node %in% net$TF) {
    # 如果是，将高亮节点设置为 TF 颜色，与之相连的节点设置为目标颜色
    connected_nodes <- neighbors(network, highlight_node, mode = "all")
    V(network)[highlight_node]$color <- adjustcolor(highlight_TF_color, alpha.f = 1)
    for (node in names(connected_nodes)) {
      if (node %in% net$TF) {
        V(network)[node]$color <- adjustcolor("#f19d00", alpha.f = 0.9)
      } else if (node %in% DErslt$up) {
        V(network)[node]$color <- adjustcolor("#fb8072", alpha.f = 0.9)
      } else {
        V(network)[node]$color <- adjustcolor("#80b1d3", alpha.f = 0.9)
      }
    }
    
    #    V(network)[connected_nodes]$color <- adjustcolor(highlight_Target_color, alpha.f = 0.8)
    
  } else {
    # 如果不是，将高亮节点设置为目标颜色，与之相连的节点设置为 TF 颜色
    connected_nodes <- neighbors(network, highlight_node, mode = "all")
    V(network)[highlight_node]$color <-  ifelse(highlight_node %in% DErslt$up,
                                                adjustcolor("#fb8072", alpha.f = 0.9),
                                                adjustcolor("#80b1d3", alpha.f = 0.9))
    for (node in names(connected_nodes)) {
      if (node %in% net$TF) {
        V(network)[node]$color <- adjustcolor("#f19d00", alpha.f = 0.9)
      } else if (node %in% DErslt$up) {
        V(network)[node]$color <- adjustcolor("#fb8072", alpha.f = 0.9)
      } else {
        V(network)[node]$color <- adjustcolor("#80b1d3", alpha.f = 0.9)
      }
    }
    
    #    V(network)[connected_nodes]$color <- adjustcolor(highlight_TF_color, alpha.f = 0.7)
  }
  
  # 获取所有与 highlight_node 相连的边
  connected_edges <- incident_edges(network, highlight_node, mode = "all")
  
  # 设置相关节点和边的颜色
  E(network)[unlist(connected_edges)]$color <- highlight_edge_color
  
  #设置网络图的布局
  if (is.null(netlayout)) {
    netlayout <- layout_(network, 
                         with_kk(kkconst = max(vcount(network) + 50, 1)))
  }
  
  # 设置 PDF 的输出文件名、高度和宽度
  pdf.name <- paste0("TFnet_highlight_", highlight_node, ".pdf")
  
  pdf(file = pdf.name, height = pdf_height, width = pdf_width)
  
  # 绘制网络图
  plot(network, edge.arrow.size = 0.4,
       layout = netlayout, 
       vertex.label.family = vertex_label_family, 
       res = c(1200, 1200))
  
  # 关闭 PDF 设备
  dev.off()
}

#' @title Plot interactive TF network
#' @param network TF network
#' @return interactive TF network plot
#' @export
#' @import igraph
plot_interactive_tfnetwork <- function(network,
                                       canvas.width = 500, 
                                       canvas.height = 500,
                                       vertex.label.family = "sans")
{
  tkplot(network, canvas.width = canvas.width, canvas.height = canvas.height,
         vertex.label.family = vertex.label.family,
         clean.on.unload=FALSE)
}

#' @title Plot 3D TF network
#' @param network TF network
#' @return 3d TF network plot
#' @export
#' @import igraph
#' @import rgl
plot_3d_tfnetwork <- function(network)
{
  rglplot(network)
}
