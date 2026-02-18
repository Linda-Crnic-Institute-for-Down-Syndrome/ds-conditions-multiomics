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


# Functions to generate final results tbls in our format (use with for loops to get all comparisons) ------
get_results_tbl <- function(x, contrast = c(predictor, "TRUE", "FALSE"), cooks = TRUE, ind_filt = TRUE, shrink_type = "apeglm", norm_counts) {
  # relevel $has_cond and run nbinomWaldTest() to ensure that coefficient of
  # otherwise need to use contrast argument which does not allow use of apeglm
  # shrinkage
  x$has_cond <- x$has_cond %>% relevel(contrast[3]) 
  
  message("------\nRunning negative binomial Wald test for ", paste(contrast[2], "vs", contrast[3], sep="_"))
  x <- x %>% nbinomWaldTest(quiet=TRUE)
  # Get results without LFC shrinkage (default)
  res <- x %>%
    results(contrast,
            cooksCutoff = cooks,
            independentFiltering = ind_filt
    )
  res_tbl <- res %>% # Could also use biobroom::tidy() on results or dds but has less useful colnames
    as.data.frame() %>%
    as_tibble(rownames = "Geneid")
  # Apply LFC shrinkage
  message("Calculating log2 fold-change shrinkage")
  res_shrink <- x %>%
    lfcShrink(coef = paste(contrast[1], contrast[2], "vs", contrast[3], sep="_"), # REPLACE dashes or this fails
              type = shrink_type, # "apeglm" or "normal" or "ashr"; apeglm and ashr are better at preserving large LFCs
              parallel = TRUE,
              res = res
    )
  res_shrink_tbl <- res_shrink %>%
    as.data.frame() %>%
    as_tibble(rownames="Geneid")
  # Combine and mutate to get final results tbl
  message("Assembling final results table")
  res_shrink_tbl %>%
    inner_join(res_tbl %>% select(Geneid, log2FoldChange), by = "Geneid") %>%
    dplyr::rename(log2FoldChange = log2FoldChange.x,
                  log2FoldChange_adj = log2FoldChange.y) %>%
    mutate(FoldChange = 2^log2FoldChange,
           FoldChange_adj = 2^log2FoldChange_adj) %>%
    inner_join(gene_anno, by="Geneid") %>%
    dplyr::select(
      Geneid:baseMean,
      # matches(denom_label),
      # matches(num_label),
      FoldChange,
      log2FoldChange,
      FoldChange_adj,
      log2FoldChange_adj,
      pvalue,
      padj
    ) %>%
    arrange(padj) %>%
    inner_join(gene_anno, by="Geneid") %>%
    select(Gene_name = gene_name, chr, everything())
}

# GSEA FUNCTIONS -----
# function to get combined pos and neg GSEA results -----
run_fgsea2 <- function(geneset, ranks, weighted = FALSE) {
  library("fgsea")
  # with gseaParam = 0, results are VERY similar to original GSEA # this seems to not be operating as expected as N^0 = 1, so all ranking stats would be 1
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

# labelled Volcano plot function ----
volcano_plot_lab <- function(res, 
                             labels = TRUE,
                             n_labels = 15,
                             title = "", 
                             subtitle = "",
                             y_lim = c(0, NA),
                             raster = FALSE
){res <- res %>% 
  mutate(
    color = if_else(padj < 0.1, "padj < 0.1", "All")
  )
# Get max finite -log10(pval) and replace 0s if needed
max_finite <- res %>%
  filter(padj > 0) %$%
  min(padj) %>%
  -log10(.)
res <- res %>% mutate(
  shape = if_else(padj == 0, "infinite", "finite"),
  padj = if_else(padj == 0, 10^-(max_finite * 1.05), padj)
)
# get max for x-axis
x_lim <- res %>% 
  summarize(max = max(log2(FoldChange_adj), na.rm = TRUE), min = min(log2(FoldChange_adj), na.rm = TRUE)) %>% 
  abs() %>% 
  max() %>% 
  ceiling()
p <- res %>% 
  ggplot(aes(log2(FoldChange_adj), -log10(padj), color = color, shape = shape)) + 
  geom_hline(yintercept = -log10(0.1), linetype = 2) + 
  geom_vline(xintercept = 0, linetype = 2) + 
  geom_point() + 
  scale_color_manual(values = c("padj < 0.1" = "red", "All" = "black")) + 
  scale_shape_manual(values = c("infinite" = 2, "finite" = 16)) +
  guides(shape = "none") +
  # xlim(-x_lim, x_lim) +
  # ylim(y_lim) +
  theme(aspect.ratio=1.2) +
  labs(
    title = title,
    subtitle = subtitle
  )
if(labels == TRUE) {
  p <- p +
    geom_text_repel(data = res %>% slice_min(order_by = padj, n = n_labels), aes(label = Gene_name), min.segment.length = 0, show.legend = FALSE) + # , nudge_y = -max_finite / 5
    geom_text_repel(data = res %>% filter(padj < 0.1) %>% slice_min(order_by = FoldChange_adj, n = n_labels), aes(label = Gene_name), min.segment.length = 0, show.legend = FALSE, nudge_x = -x_lim/4, nudge_y = 0, ylim = c(max_finite / 10, NA)) +
    geom_text_repel(data = res %>% filter(padj < 0.1) %>% slice_max(order_by = FoldChange_adj, n = n_labels), aes(label = Gene_name), min.segment.length = 0, show.legend = FALSE, nudge_x = x_lim/4, nudge_y =  max_finite /5, ylim = c(max_finite / 10, NA))
}
#
if(raster) {
  # to rasterize all points:
  # rasterize(p, layers='Point', dpi = 600)
  p <- rasterize(p, layers='Point', dpi = 600, dev = "ragg_png")
  # otherwize can rasterize individual layers using rasterize(geom_sina())
}
#
return(p)
} # end of function
