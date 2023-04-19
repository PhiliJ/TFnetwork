#' @title Construct Network for plotting
#' @description It constructs the TF-Target network and returns network file ("TF_Network.txt") if you want to plot your network using other software like Cytoscape.
#' @param mouse A logical indicating whether to construct the network for mouse or human. Default is FALSE for human.
#' @return net: TF-Target network
#' @return The constructed network to a file named "TF_Network.txt"
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
#' @description This function plots a transcription factor regulatory network, including transcription factors and their target genes. The network is based on a given data frame that contains information about the transcription factors and their target genes.
#' @param net A data frame that contains information about the transcription factors and their target genes. It should have two columns: "TF" and "Target".
#' @param nodeshape The shape of the nodes. Default is "circle". If you want to draw a 3D plot, use "sphere" instead.
#' @param TF_color The color of the transcription factor nodes. Default is "#F0E442".
#' @param Up_color The color of the upregulated target gene nodes. Default is "#fb8072".
#' @param Down_color The color of the downregulated target gene nodes. Default is "#80b1d3".
#' @param TF_alpha The alpha value of the transcription factor nodes. Default is 0.75.
#' @param Up_alpha The alpha value of the upregulated target gene nodes. Default is 0.3.
#' @param Down_alpha The alpha value of the downregulated target gene nodes. Default is 0.3.
#' @param TF_size The size of the transcription factor nodes. Default is 10.
#' @param Target_size The size of the target gene nodes. Default is 8.
#' @param TF_labelsize The label size of the transcription factor nodes. Default is 0.7.
#' @param Target_labelsize The label size of the target gene nodes. Default is 0.7.
#' @param label_color The color of the node labels. Default is "black".
#' @param edge_color The color of the edges. Default is "grey".
#' @param edge_alpha The alpha value of the edges. Default is 0.75.
#' @param netlayout The layout of the network.  Default is NULL, and the function will calculate the layout using the Kamada-Kawai algorithm.
#' @param pdf.name The name of the PDF file to be saved. Default is "TFnet.pdf".
#' @param pdf.height The height of the PDF file. Default is 15.
#' @param pdf.width The width of the PDF file. Default is 15.
#' @param edge.arrow.size The size of the arrows on the edges. Default is 0.4.
#' @param vertex.label.family The font family of the node labels. Default is "sans".
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


#' @title Plot a TF network with highlighted TF or target gene node
#' @description This function plots a transcription factor (TF) network with a specified TF or target gene node highlighted in a different color. The function takes as input the network and the node to be highlighted, and outputs a PDF file of the network plot with the highlighted node and its connected edges colored differently.
#' @param network A graph object representing the TF network, created using the igraph package.
#' @param highlight_node A character string representing the name of the node to be highlighted in the network plot. This can be either a TF or a target gene.
#' @param highlight_TF_color A character string representing the color to be used for highlighting the TF node. Default is "#f19d00".
#' @param highlight_edge_color A character string representing the color to be used for highlighting the edges connected to the highlighted node. Default is "#e41a1c".
#' @param netlayout A layout object representing the layout of the network. Default is NULL, and the function will calculate the layout using the Kamada-Kawai algorithm.
#' @param pdf_height An integer representing the height of the PDF output file in inches. Default is 15.
#' @param pdf_width An integer representing the width of the PDF output file in inches. Default is 15.
#' @param vertex_label_family A character string representing the font family for the vertex labels in the network plot. Default is "sans".
#' @return This function outputs a PDF file of the network plot with the highlighted node and its connected edges colored differently.
#' @export
#' @import igraph
#' @import magrittr
#' @import dplyr
#' @import ggplot2
#' @import circlize
#' @examples
#' # Load example network and DE results
#' data(net)
#' data(DErslt)
#'
#' # Plot the network with the TF node "FOXA1" highlighted
#' plot_highlighted_tfnetwork(net, highlight_node = "FOXA1")
#'
#' # Plot the network with the target gene node "IGFBP3" highlighted
#' plot_highlighted_tfnetwork(net, highlight_node = "IGFBP3")
plot_highlighted_tfnetwork <- function(network, 
                                       highlight_node, 
                                       highlight_TF_color = "#f19d00", 
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

#' @title Plot an interactive TF network using tkplot from igraph
#' @description This function plots an interactive transcription factor (TF) network using the tkplot function from the igraph package. The resulting plot can be manipulated using mouse gestures and keyboard commands.
#' @param network The igraph object representing the TF network to be plotted
#' @param canvas.width The width of the canvas for the tkplot window
#' @param canvas.height The height of the canvas for the tkplot window
#' @param vertex.label.family The font family to use for vertex labels
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

#' @title Plot a 3D TF network using tkplot from igraph
#' @description This function plots the input network as a 3D network graph using the rgl package.
#' @param network an igraph object representing the network to plot
#' @return A 3D network graph plot.
#' @export
#' @import igraph
#' @import rgl
#' @examples
#' plot_3d_tfnetwork(network)
plot_3d_tfnetwork <- function(network)
{
  rglplot(network)
}
