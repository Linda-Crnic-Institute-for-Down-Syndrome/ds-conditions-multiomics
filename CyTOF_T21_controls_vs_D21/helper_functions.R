# Common helper functions
mem_used <- function() lobstr::mem_used() %>% as.numeric() %>% R.utils::hsize()
obj_size <- function(x) object.size(x) %>% print(units = "auto")

## Setting and modifying default theme for plots
theme_set(theme_gray(base_size=12, base_family="Arial") +
            theme(
              panel.border=element_rect(colour="black", fill="transparent"), 
              plot.title=element_text(face="bold", hjust=0),
              axis.text=element_text(color="black", size=14), 
              axis.text.x=element_text(angle=0, hjust=0.5),
              axis.ticks = element_line(color = "black"), # make sure tick marks are black
              panel.background=element_blank(),
              panel.grid=element_blank(),
              plot.background=element_blank(),
              strip.background = element_blank(), # facet label borders
              legend.key=element_blank(), # remove grey bg from legend elements
              legend.background=element_blank() 
            )
)

# labelled Volcano plot function -----
volcano_plot_lab <- function(res, title = "",
                             n_labels = 10,
                                subtitle = "down in Pos.                                                  up in Pos.",
                                y_lim = c(0, NA)){
  theme_set(theme_gray(base_size=12, base_family="Arial") +
              theme(panel.border=element_rect(colour="black", fill="transparent"), 
                    plot.title=element_text(face="bold", hjust=0), # lineheight=.8, size=20,
                    axis.text=element_text(color="black", size=14), 
                    axis.text.x=element_text(angle=0, hjust=0.5),
                    panel.background=element_blank(),
                    panel.grid=element_blank(),
                    plot.background=element_blank(),
                    strip.background = element_blank(), # facet label borders
                    legend.key=element_blank(), legend.background=element_blank() # remove grey bg from legend
              )
  )
  res <- res %>% 
    mutate(
      color = if_else(BHadj_pval < 0.1, "padj < 0.1", "All")
    )
  # get max for x-axis
  x_lim <- res %>% 
    summarize(max = max(log2(FoldChange), na.rm = TRUE), min = min(log2(FoldChange), na.rm = TRUE)) %>% 
    abs() %>% 
    max() %>% 
    # ceiling()
    plyr::round_any(accuracy = 0.025, f = ceiling) # CUSTOMIZE
  res %>% 
    ggplot(aes(log2(FoldChange), -log10(BHadj_pval), color = color)) + 
    geom_hline(yintercept = -log10(0.1), linetype = 2) + 
    geom_vline(xintercept = 0, linetype = 2) + 
    geom_point() + 
    scale_color_manual(values = c("padj < 0.1" = "red", "All" = "black")) + 
    xlim(-x_lim, x_lim) +
    ylim(y_lim) +
    geom_text_repel(data = res %>% filter(!is.na(BHadj_pval) & FoldChange > 1) %>% slice_max(order_by = FoldChange, n = 15), aes(label = cluster), min.segment.length = 0, show.legend = FALSE, nudge_x = 0.01, nudge_y = 0.01, max.overlaps = Inf) +
    geom_text_repel(data = res %>% filter(!is.na(BHadj_pval) & FoldChange < 1) %>% slice_min(order_by = FoldChange, n = 15), aes(label = cluster), min.segment.length = 0, show.legend = FALSE, nudge_x = -0.01, nudge_y = 0.01, max.overlaps = Inf) +
    # geom_text_repel(data = res %>% filter(!is.na(BHadj_pval)) %>% slice_min(order_by = BHadj_pval, n = n_labels), aes(label = uniquePopulationName), min.segment.length = 0, show.legend = FALSE, nudge_x = -0, nudge_y = 0.1) +
    theme(aspect.ratio=1.2) +
    labs(
      title = title,
      subtitle = subtitle
    )
} # end of function

