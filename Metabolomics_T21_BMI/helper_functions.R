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

## Density color function
getDenCols <- function(x, y, transform = TRUE) { # set to TRUE if using log2 transformation of data
  if(transform) {
    df <- data.frame(log2(x), log2(y))
  } else{
    df <- data.frame(x, y)
  }
  z <- grDevices::densCols(df, colramp = grDevices::colorRampPalette(c("black", "white")))
  df$dens <- grDevices::col2rgb(z)[1,] + 1L
  cols <-  grDevices::colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
  df$col <- cols[df$dens]
  return(df$dens)
} # End of function

# labelled Volcano plot function for lms -----
volcano_plot_lab_lm <- function(res, title = "", 
                                subtitle = "down in Pos.                                                  up in Pos.",
                                y_lim = c(0, NA)){
  res <- res %>% 
    mutate(
      color = if_else(BHadj_pval < 0.1, "qval < 0.1", "All")
    )
  # get max for x-axis
  x_lim <- res %>% 
    summarize(max = max(log2(FoldChange), na.rm = TRUE), min = min(log2(FoldChange), na.rm = TRUE)) %>% 
    abs() %>% 
    max() %>% 
    ceiling()
  res %>% 
    ggplot(aes(log2(FoldChange), -log10(BHadj_pval), color = color)) + 
    geom_hline(yintercept = -log10(0.1), linetype = 2) + 
    geom_vline(xintercept = 0, linetype = 2) + 
    geom_point() + 
    scale_color_manual(values = c("qval < 0.1" = "red", "All" = "black")) + 
    xlim(-x_lim, x_lim) +
    ylim(y_lim) +
    geom_text_repel(data = res %>% filter(!is.na(BHadj_pval)) %>% filter(BHadj_pval < 0.1 & FoldChange>1) %>% slice_max(order_by = -BHadj_pval, n=10), aes(label = Analyte), min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE, nudge_x = 0.003, nudge_y = 0.05, family = "Arial") +
    geom_text_repel(data = res %>% filter(!is.na(BHadj_pval)) %>% filter(BHadj_pval < 0.1 & FoldChange<1) %>% slice_max(order_by = -BHadj_pval, n=10), aes(label = Analyte), min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE, nudge_x = -0.003, nudge_y = 0.05) +
    #geom_text_repel(data = res %>% filter(!is.na(BHadj_pval)) %>% slice_min(order_by = BHadj_pval, n = 3), aes(label = Aptamer), min.segment.length = 0, show.legend = FALSE, nudge_x = -0, nudge_y = 0.1) +
    theme(aspect.ratio=1.2) +
    labs(
      title = title,
      subtitle = subtitle
    ) 
} # end of function

# GSEA FUNCTIONS -----
# function to get combined pos and neg GSEA results -----
run_fgsea2 <- function(geneset, ranks, weighted = FALSE) {
  library("fgsea")
  weight = 0
  if(weighted) weight = 1#0.5
  # Run positive enrichment
  fgseaRes_POSITIVE <- fgseaMultilevel(
    geneset, 
    ranks, 
    minSize=15, 
    maxSize=500,
    gseaParam = weight,
    # nperm = 1000,
    eps = 0.0, # fgsea has a default lower bound eps=1e-10 for estimating P-values. If you need to estimate P-value more accurately, you can set the eps argument to zero
    scoreType = "pos"
  )
  # Run negative enrichment
  fgseaRes_NEGATIVE <- fgseaMultilevel(
    geneset,
    ranks,
    minSize=15,
    maxSize=500,
    gseaParam = 0,
    # nperm = 1000,
    eps = 0.0, # fgsea has a default lower bound eps=1e-10 for estimating P-values. If you need to estimate P-value more accurately, you can set the eps argument to zero
    scoreType = "neg"
  )
  # Combine positive and negative results + re-adjust pvals
  fgseaRes_POS_NEG <- inner_join(
    fgseaRes_POSITIVE %>% 
      as_tibble(),
    fgseaRes_NEGATIVE %>% 
      as_tibble(),
    by = c("pathway"),
    suffix = c("_POS", "_NEG")
  )
  fgseaRes_COMBINED <- bind_rows(
    fgseaRes_POS_NEG %>% filter(ES_POS > abs(ES_NEG)) %>% select(pathway) %>% inner_join(fgseaRes_POSITIVE),
    fgseaRes_POS_NEG %>% filter(ES_POS < abs(ES_NEG)) %>% select(pathway) %>% inner_join(fgseaRes_NEGATIVE)
  ) %>% 
    mutate(padj = p.adjust(pval, method = "BH"))%>% 
    arrange(padj, -abs(NES))
  return(fgseaRes_COMBINED)
} # end of function

# Customized version of plotEnrichment -----
plotEnrichment2 <- function (pathway, stats, res, title = "") 
{
  # pathway = hallmarks$HALLMARK_INTERFERON_GAMMA_RESPONSE # for testing
  # stats = ranks # for testing
  # title = "Hallmark Interferon Gamma Response"
  gseaParam = 0
  ticksSize = 0.4
  pathname <- deparse(substitute(pathway)) %>% str_remove("\\w+\\$")
  label = paste0(
    "NES = ", (res %>% filter(pathway == pathname))$NES %>% round(2),
    "\n",
    "Q = ", (res %>% filter(pathway == pathname))$padj %>% rstatix::p_format()
  )
  x_label <- length(stats)*0.99
  y_label <- ((res %>% filter(pathway == pathname))$ES)*0.95
  # Setting and modifying default theme for plots
  theme_set(theme_gray(base_size=12, base_family="Arial") +
              theme(panel.border=element_rect(colour="black", fill="transparent"), 
                    plot.title=element_text(face="bold", hjust=0),
                    axis.text=element_text(color="black", size=12), 
                    axis.text.x=element_text(angle=0, hjust=0.5),
                    # axis.text.x=element_text(angle=90, hjust=0.5),
                    # axis.text.x=element_text(angle=45, hjust=1),
                    panel.background=element_blank(),
                    panel.grid=element_blank(),
                    plot.background=element_blank()
              ) +
              # theme(strip.background=element_rect(colour="black", fill="light grey", size=1))
              theme(strip.background = element_blank()) # adjusts facet label borders
  )
  #
  rnk <- rank(-stats) # rank highest values first
  ord <- order(rnk) # get correct order
  statsAdj <- stats[ord] # ensure ranked list is ordered correctly
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam) # gets sign and multiplies by absolute value ^ gsea param
  statsAdj <- statsAdj/max(abs(statsAdj))
  zero_cross <- statsAdj[statsAdj > 0] %>% length() # New; get Zero crossing point
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(stats = statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  g <- ggplot(toPlot, aes(x = x, y = y)) + 
    # geom_point(color = "green", size = 0.1) + # why bother?
    geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") + 
    geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") + 
    geom_hline(yintercept = 0, colour = "black") + 
    geom_vline(xintercept = zero_cross, linetype = 2, color = "grey50") +
    geom_line(color = "green") + # this is the main running ES score
    geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), size = ticksSize) + # gene set ticks
    theme(aspect.ratio = 0.7) +
    annotate(geom = "text", x = zero_cross + (n / 50), y = 0.025, label = paste0("Zero cross at ", zero_cross), hjust = 0, size=4, family="Arial") +
    annotate(geom = "text", x_label, y_label, label = label, hjust = 1, family="Arial") +
    labs(title = title, x = "Rank", y = "Enrichment score")
  g
} # End of function

