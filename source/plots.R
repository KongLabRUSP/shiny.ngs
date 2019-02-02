# input: a data table returned by function add.clm(), q_value, fold_change
# output: a MA plot, ggplot object

ma <- function(dt.tb, q_value, fold_change){
  p <- ggplot(dt.tb,
              aes(x = mu,
                  y = `log2(Fold_change) normalized`,
                  colour = clr,
                  shape = pch)) +
    scale_shape_manual(name = "Legend:",
                       labels = c("No significance",
                                  paste("p-Value < ",
                                        as.numeric(q_value)),
                                  paste("p-Value < ",
                                        as.numeric(q_value),
                                        " & abs(log2) >= ",
                                        as.numeric(fold_change))),
                       values = c(46, 3, 4)) +
    scale_color_manual(name = "Legend:",
                       values = c("black",
                                  "purple",
                                  "red"),
                       labels = c("No significance",
                                  paste("p-Value < ",
                                        as.numeric(q_value)),
                                  paste("p-Value < ",
                                        as.numeric(q_value),
                                        " & abs(log2) >= ",
                                        as.numeric(fold_change)))) +
    geom_hline(yintercept = c(-as.numeric(fold_change), as.numeric(fold_change)),
               lty = 2) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          legend.position = "top") +
    geom_point() +
    labs(x = "log2Mean", y = "log2FoldChange")
  
  return(p)
  
}

# two circle heatmap 
# input: a sorted up.dn_dn.up_table, like the one produced by cige function
# output: a two cricle heatmap of the normalized logfoldcgange
two.cl.heat <- function(up.dn_dn.up_table){
  up.dn_dn.up_table_long <- melt(data = up.dn_dn.up_table,
                                 id.vars = 1,
                                 measure.vars = 2:3,
                                 variable.name = "Comparison",
                                 value.name = "Gene Expression Diff")
  up.dn_dn.up_table_long$Comparison <- factor(up.dn_dn.up_table_long$Comparison,
                                              # inner circle trt3-trt2 and outer circle trt2-trt1
                                              levels = c(colnames(up.dn_dn.up_table)[3],
                                                         colnames(up.dn_dn.up_table)[2]))
  lvls <- up.dn_dn.up_table_long[up.dn_dn.up_table_long$Comparison == colnames(up.dn_dn.up_table)[2], ]
  up.dn_dn.up_table_long$gene <- factor(up.dn_dn.up_table_long$gene,
                                        levels = lvls$gene[order(lvls$`Gene Expression Diff`)])
  
  two_circle_heatmap <- ggplot(data = up.dn_dn.up_table_long) +
    coord_polar("y",
                start = 0,
                direction = -1) +
    geom_tile(aes(x = as.numeric(Comparison),
                  y = gene,
                  fill = `Gene Expression Diff`),
              color = "white") +
    geom_text(data = up.dn_dn.up_table_long[Comparison == colnames(up.dn_dn.up_table)[2], ],
              aes(x = rep(1.75,
                          nlevels(gene)),
                  y = gene,
                  label = unique(gene),
                  angle = 90 + seq(from = 0,
                                   to = 360,
                                   length.out = nlevels(gene))[as.numeric(gene)]),
              hjust = 0) +
    scale_fill_gradient2(low = "red",
                         high = "green",
                         mid = "grey",
                         midpoint = 0,
                         name = "Gene Expr Diff") +
    scale_x_continuous(limits = c(0,
                                  max(as.numeric(up.dn_dn.up_table_long$Comparison)) + 0.5),
                       expand = c(0, 0)) +
    scale_y_discrete("",
                     expand = c(0, 0)) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  return(two_circle_heatmap)
  
}