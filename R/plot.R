#' Visualize control PSMs versus test PSMs on a log-log plot.
#'
#' Graphs control PSMs versus test PSMs on a log-log plot. A pseudocount of 1 is adding to all PSMs for log-log visualization. Colored by Z-test threshold. See README for file formats and example usage on the \href{https://github.com/rachaelcox/diffprot}{diffprot GitHub Repository}.
#'
#' @param data A data frame produced by `enrich()`.
#' @param ylab Text for y-axis label; will correspond to test PSMs.
#' @param xlab Text for x-axis label; will correspond to control PSMs.
#' @param threshold Optional. FDR threshold to color points by. One of c(90, 95, 99). Default = 5% FDR (fdr_bh <= 0.05).
#' @param num_labs Optional. Number of annotation labels to display. Requires input to have an annotation column named `gene_names_primary`. Default = 0.
#' @param  label_file Optional. Comma-separated file containing names of genes to be labeled on the plot; requires input to have an annotation column named `gene_names_primary` and labels must match format in the column.
#' @param outfile_prefix Optional. File prefix for output `.png` and `.pdf` files. Default = "output."
#' @return PSM log-log plots in `.png` and `.pdf` format.
#' @export
psmplot <- function(data, outfile_prefix, threshold, num_labs,
                    label_file, ylab, xlab, point_color){

  theme_set(cowplot::theme_cowplot())
  palette_pretty <- c("#0072B2","#E69F00","#009E24",
                      "#FF0000", "#979797","#5530AA", "#1E1E1E")

  log10_minor_break <- function (...){
    function(x) {
      minx         = floor(min(log10(x), na.rm=T))-1;
      maxx         = ceiling(max(log10(x), na.rm=T))+1;
      n_major      = maxx-minx+1;
      major_breaks = seq(minx, maxx, by=1)
      minor_breaks =
        rep(log10(seq(1, 9, by=1)), times = n_major)+
        rep(major_breaks, each = 9)
      return(10^(minor_breaks))
    }
  }

  if(missing(outfile_prefix)){
    outfile_prefix = "output"
  }

  if(missing(num_labs)){
    num_labs = 10
  }

  if(missing(threshold)){
    threshold = 99
  }

  if(missing(point_color)){
    point_color = palette_pretty[4]
  }

  xcol <- grep('ctrl_PSMs$', names(data), value = TRUE)
  ycol <- grep('exp_PSMs$', names(data), value = TRUE)
  zcol <- grep('_zscore$', names(data), value = TRUE)

  if(any(grepl('gene_names_primary', names(data)))){
    acol <- grep('gene_names_primary$',
                 names(data), value = TRUE)
  } else {
    acol <- grep('accession',
                 names(data), value = TRUE)
  }

  # calculate thresholds
  data_conf <- data %>%
    dplyr::mutate(conf_90 = dplyr::case_when(fdr_bh <= 0.10 ~ TRUE,
                                      TRUE ~ FALSE),
                  conf_95 = dplyr::case_when(fdr_bh <= 0.05 ~ TRUE,
                                      TRUE ~ FALSE),
                  conf_99 = dplyr::case_when(fdr_bh <= 0.01 ~ TRUE,
                                      TRUE ~ FALSE))

  # add 1 to all PSMs for log-log visualization
  print('Adding pseudocounts to PSMs...')

  data_conf[ycol] <- data_conf[ycol]+1
  data_conf[xcol] <- data_conf[xcol]+1

  data_conf %>%
    dplyr::select(accession, xcol, ycol) %>%
    print()

  if(missing(label_file)){
    label_subset <- data_conf %>%
      dplyr::arrange(desc(zcol)) %>%
      dplyr::slice_head(n = num_labs)
  } else {
    labels <- readr::read_csv(label_file, col_names = FALSE)
    label_subset <- data_conf %>%
      dplyr::filter(grepl(paste(labels$X1, collapse = "|"), get(acol),
                   ignore.case = TRUE))
  }

  print(sprintf('Labeled plot points for "%s_PSMloglog.png"', outfile_prefix))
  label_subset %>%
    dplyr::select(accession, zcol, acol) %>%
    print()

  if(threshold == 90){
    p1 <- ggplot(data_conf, aes(x = get(xcol), y = get(ycol),
                                color = conf_90,  label = get(acol))) +
      geom_point(alpha = 0.75, size = 1.5) +
      annotate(geom = "text", x = Inf, y = 1, vjust = 1, hjust = 1,
               label = "10% FDR", color = point_color, size = 4)

  } else if(threshold == 95) {
    p1 <- ggplot(data_conf, aes(x = get(xcol), y = get(ycol),
                                color = conf_95, label = get(acol))) +
      geom_point(alpha = 0.75, size = 1.5) +
      annotate(geom = "text", x = Inf, y = 1, vjust = 1, hjust = 1,
               label = "5% FDR", color = point_color, size = 4)

  } else {
    p1 <- ggplot(data_conf, aes(x = get(xcol), y = get(ycol),
                                color = conf_99, label = get(acol))) +
      geom_point(alpha = 0.75, size = 1.5) +
      annotate(geom = "text", x = Inf, y = 1, vjust = 1, hjust = 1,
               label = "1% FDR", color = point_color, size = 4)

  }

  p1
  p2 <- p1 + ggrepel::geom_text_repel(data = label_subset,
                             size = 6.5/.pt, # font size
                             fontface = "bold",
                             alpha = 1.0,
                             hjust = 0.2,
                             point.padding = unit(0.5, "lines"),
                             box.padding = unit(0.5, "lines"),
                             nudge_y = 0.3,
                             nudge_x = -0.3,
                             colour = palette_pretty[7],
                             set.seed(3333)) +
    scale_x_log10(minor_breaks = log10_minor_break()) +
    scale_y_log10(minor_breaks = log10_minor_break()) +
    scale_color_manual(values = rev(c(point_color, palette_pretty[5]))) +
    theme(legend.position = "none") +
    theme(panel.grid.major = element_line(colour = palette_pretty[5], size = 0.1, linetype = "dashed"),
          panel.grid.minor = element_line(colour = palette_pretty[5], size = 0.1, linetype = "dashed")) +
    xlab(xlab) +
    ylab(ylab)

  ggsave(sprintf("%s_psmplot.pdf", outfile_prefix), plot = p2,
         device = "pdf", width = 8, height = 6, units = "in")

  ggsave(sprintf("%s_psmplot.png", outfile_prefix), plot = p2,
         device = "png", width = 8, height = 6, units = "in")

  print(p2)
  return(p2)
}

#' Plot Z-scores across biological replicates.
#'
#' Colored by a multiple hypothesis corrected false discovery rate (FDR) computed on a weighted Z-score. See README for file formats and example usage on the \href{https://github.com/rachaelcox/diffprot}{diffprot GitHub Repository}.
#'
#' @param data A data frame produced by `combine_reps()`.
#' @param ycol Name of column containing y-axis data; should correspond to replicate #1 Z-scores.
#' @param xcol Name of column containing x-axis data; should correspond to replicate #1 Z-scores.
#' @param ylab Text for y-axis label; should correspond to replicate #1 Z-scores.
#' @param xlab Text for x-axis label; should correspond to replicate #2 Z-scores.
#' @param threshold Optional. FDR threshold to color points by. One of c(90, 95, 99). Default = 5% FDR (fdr_bh <= 0.05).
#' @param num_labs Optional. Number of annotation labels to display. Requires input to have an annotation column named `gene_names_primary`. Default = 0.
#' @param  label_file Optional. Comma-separated file containing names of genes to be labeled on the plot; requires input to have an annotation column named `gene_names_primary` and labels must match format in the column.
#' @param outfile_prefix Optional. File prefix for output `.png` and `.pdf` files. Default = "output."
#' @return Plot depicting replicate 1 z-score vs replicate 2 z-scores in `.png` and `.pdf` format. Colored by FDR.
#' @export
zplot <- function(data, outfile_prefix, ylab, xlab, ycol, xcol,
                      threshold, label_file, num_labels){

  theme_set(cowplot::theme_cowplot())
  palette_pretty <- c("#1E1E1E", "#FF0000", "#0072B2",
                      "#E69F00", "#009E24", "#979797","#5530AA", "#1E1E1E")

  if(missing(outfile_prefix)){
    outfile_prefix = "output"
  }

  if(missing(num_labels)){
    num_labels = 10
  }

  if(missing(threshold)){
    threshold = 95
  }

  if(any(grepl('gene_names_primary', names(data)))){
    acol <- grep('gene_names_primary$',
                 names(data), value = TRUE)
  } else {
    acol <- grep('accession',
                 names(data), value = TRUE)
  }

  zcol <- grep('joint_zscore', names(data), value = TRUE)

  if(length(zcol) > 1){
    zcol <- zcol[2]
  }

  if(missing(label_file)){   # take top n labels if no label file
    label_subset <- data %>%
      dplyr::arrange(desc(zcol)) %>%
      dplyr::slice_head(n = num_labels)

  } else {
    labels <- readr::read_csv(label_file, col_names = FALSE)
    label_subset <- data %>%
      dplyr::filter(grepl(paste(labels$X1, collapse = "|"), get(acol), ignore.case = TRUE))
  }

  print(sprintf('Labeled plot points for "%s_zscore_comparison.png & %s_zscore_comparison.pdf:"', outfile_prefix, outfile_prefix))
  label_subset %>%
    dplyr::select(accession, matches("zscore"), matches("gene_names_primary")) %>%
    print()

  # z-scores vs z-scores plot

  if(threshold == 90){
    p_base <- ggplot(data, aes(x = get(xcol), y = get(ycol),
                               color = conf_90, label = get(acol))) +
      annotate(geom = "text",
               x = max(data[xcol]),
               y = min(data[ycol]),
               vjust = 1,
               label = "10% FDR",
               color = palette_pretty[2], size = 4) +
      geom_hline(yintercept = 1.282, alpha = 0.6) +
      geom_hline(yintercept = -1.282, alpha = 0.6) +
      geom_vline(xintercept = 1.282, alpha = 0.6) +
      geom_vline(xintercept = -1.282, alpha = 0.6) +
      scale_x_continuous(breaks = c(-5, -1.282, 0, 1.282, 5, 10),
                         labels = c(-5, -1.282, 0, 1.282, 5, 10)) +
      scale_y_continuous(breaks = c(-5, -1.282, 0, 1.282, 5, 10),
                         labels = c(-5, -1.282, 0, 1.282, 5, 10))

  } else if(threshold == 95) {
    p_base <- ggplot(data, aes(x = get(xcol), y = get(ycol),
                                      color = conf_95, label = get(acol))) +
      #labs(color = "95% Confidence\n") +
      annotate(geom = "text",
               x = max(data[xcol]),
               y = min(data[ycol]),
               vjust = 1,
               label = "5% FDR",
               color = palette_pretty[2], size = 4) +
      geom_hline(yintercept = 1.645, alpha = 0.6) +
      geom_hline(yintercept = -1.645, alpha = 0.6) +
      geom_vline(xintercept = 1.645, alpha = 0.6) +
      geom_vline(xintercept = -1.645, alpha = 0.6) +
      scale_x_continuous(breaks = c(-5, -1.645, 0, 1.645, 5, 10),
                         labels = c(-5, -1.645, 0, 1.645, 5, 10)) +
      scale_y_continuous(breaks = c(-5, -1.645, 0, 1.645, 5, 10),
                         labels = c(-5, -1.645, 0, 1.645, 5, 10))

  } else {
    p_base <- ggplot(data, aes(x = get(xcol), y = get(ycol),
                                      color = conf_99, label = get(acol))) +
      annotate(geom = "text",
               x = max(data[xcol]),
               y = min(data[ycol]),
               vjust = 1,
               label = "1% FDR",
               color = palette_pretty[2], size = 4) +
      geom_hline(yintercept = 2.33, alpha = 0.6) +
      geom_hline(yintercept = -2.33, alpha = 0.6) +
      geom_vline(xintercept = 2.33, alpha = 0.6) +
      geom_vline(xintercept = -2.33, alpha = 0.6) +
      scale_x_continuous(breaks = c(-5, -2.33, 0, 2.33, 5, 10),
                         labels = c(-5, -2.33, 0, 2.33, 5, 10)) +
      scale_y_continuous(breaks = c(-5, -2.33, 0, 2.33, 5, 10),
                         labels = c(-5, -2.33, 0, 2.33, 5, 10))
  }


  pz <- p_base +
    geom_point(alpha = 0.75, size = 1.5) +
    ggrepel::geom_text_repel(data = label_subset,
                    size = 6.5/.pt, # font size
                    fontface = "bold",
                    alpha = 1.0,
                    hjust = 0.2,
                    point.padding = unit(0.5, "lines"),
                    box.padding = unit(0.5, "lines"),
                    nudge_y = -0.3,
                    nudge_x = 0.8,
                    colour = palette_pretty[1],
                    set.seed(3333)) +
    scale_color_manual(values = palette_pretty) +
    theme(legend.position = "none") +
    #theme(legend.position = "bottom") +
    #theme(legend.position = c(0.7, 0.15),
    #      legend.title = element_blank()) +
    #guides(color = guide_legend(ncol = 1, label.theme = element_text(size = 9))) +
    theme(panel.grid.major = element_line(colour = palette_pretty[6],
                                          size = 0.1, linetype = "dashed"),
          panel.grid.minor = element_line(colour = palette_pretty[6],
                                          size = 0.1, linetype = "dashed")) +
    xlab(xlab) +
    ylab(ylab)

  ggsave(sprintf("%s_zplot.pdf", outfile_prefix), plot = pz,
         width = 8, height = 6, units = "in")

  ggsave(sprintf("%s_zplot.png", outfile_prefix), plot = pz,
         width = 8, height = 6, units = "in")

  print(pz)
  return(pz)
}
