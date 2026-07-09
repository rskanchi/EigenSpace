# =============================================================================
# pca_plots.R  --  ggplot2 builders. Each returns a ggplot object so the caller
# decides how to render/save (the export layer saves them at 600 dpi).
# =============================================================================

# Okabe-Ito colourblind-safe palette (from the original doPCA).
PCA_PALETTE <- c("#009E73", "#CC79A7", "#0072B2", "#E69F00",
                 "#661100", "#332288", "#999933", "#56B4E9", "#D55E00", "#000000")

.pca_base_theme <- function(base = 14) {
  ggplot2::theme_bw(base_size = base) +
    ggplot2::theme(legend.position = "top")
}

#' Scree plot: percent variance explained per PC (bars), with a line overlay.
plot_scree <- function(variance, n_dim = 10) {
  v <- utils::head(variance, n_dim)
  v$PC <- factor(v$PC, levels = v$PC)
  ggplot2::ggplot(v, ggplot2::aes(PC, percent_variance, group = 1)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::geom_line(colour = "red", linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_point(colour = "red", size = 2.5) +
    ggplot2::labs(title = "Scree plot: variance explained by each PC",
                  x = "Principal component", y = "Proportion of variance (%)") +
    .pca_base_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Cumulative variance explained.
plot_cumulative <- function(variance, n_dim = 10) {
  v <- utils::head(variance, n_dim)
  v$PC <- factor(v$PC, levels = v$PC)
  ggplot2::ggplot(v, ggplot2::aes(PC, cumulative_variance, group = 1)) +
    ggplot2::geom_line(colour = "grey50", linetype = "dashed") +
    ggplot2::geom_point(colour = "red", size = 2.5) +
    ggplot2::geom_text(ggplot2::aes(label = round(cumulative_variance, 1)),
                       vjust = -0.8, size = 3.5) +
    ggplot2::ylim(0, 105) +
    ggplot2::labs(title = "Cumulative variance explained",
                  x = "Principal component", y = "Cumulative variance (%)") +
    .pca_base_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

.axis_lab <- function(variance, pc) {
  i <- as.integer(sub("PC", "", pc))
  sprintf("%s (%.1f%%)", pc, variance$percent_variance[i])
}

#' Score plot with no grouping (unsupervised view).
plot_scores_plain <- function(scores, variance, pcs = c("PC1", "PC2")) {
  d <- scores; d$Sample <- rownames(scores)
  ggplot2::ggplot(d, ggplot2::aes(.data[[pcs[1]]], .data[[pcs[2]]], label = Sample)) +
    ggplot2::geom_point(size = 3, colour = "black") +
    ggplot2::geom_text(vjust = -1, size = 3) +
    ggplot2::labs(title = "PCA score plot (samples)",
                  x = .axis_lab(variance, pcs[1]), y = .axis_lab(variance, pcs[2])) +
    .pca_base_theme()
}

#' Score plot coloured by a categorical trait. ellipses = none/stat.
plot_scores_categorical <- function(scores, variance, group, group_name,
                                     pcs = c("PC1", "PC2"), ellipse = FALSE) {
  d <- scores; d$Sample <- rownames(scores); d$Group <- factor(group)
  p <- ggplot2::ggplot(d, ggplot2::aes(.data[[pcs[1]]], .data[[pcs[2]]], colour = Group)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = Sample), vjust = -1, size = 3, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = PCA_PALETTE, name = group_name) +
    ggplot2::labs(
      title = sprintf("PCA score plot coloured by %s%s", group_name,
                      if (ellipse) " (95% ellipses)" else ""),
      x = .axis_lab(variance, pcs[1]), y = .axis_lab(variance, pcs[2])) +
    .pca_base_theme()
  if (ellipse)
    p <- p + ggplot2::stat_ellipse(ggplot2::aes(group = Group), level = 0.95,
                                   type = "norm", linetype = "dashed")
  p
}


#' Score plot coloured by a quantitative trait (continuous gradient).
plot_scores_gradient <- function(scores, variance, trait_values, trait_name,
                                 pcs = c("PC1", "PC2")) {
  d <- scores; d$Sample <- rownames(scores); d$Trait <- as.numeric(trait_values)
  ggplot2::ggplot(d, ggplot2::aes(.data[[pcs[1]]], .data[[pcs[2]]],
                                  colour = Trait, label = Sample)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(vjust = -1, size = 3, colour = "grey30") +
    ggplot2::scale_colour_gradient(low = "#2C7BB6", high = "#D7191C",
                                   name = trait_name) +
    ggplot2::labs(title = sprintf("PCA score plot coloured by %s", trait_name),
                  x = .axis_lab(variance, pcs[1]), y = .axis_lab(variance, pcs[2])) +
    .pca_base_theme()
}

#' Loadings bar plot for one PC (genes ordered by contribution).
plot_loadings <- function(loadings, pc = "PC1") {
  d <- data.frame(Feature = rownames(loadings), Loading = loadings[[pc]],
                  stringsAsFactors = FALSE)
  d <- d[order(d$Loading), ]
  d$Feature <- factor(d$Feature, levels = d$Feature)
  d$Direction <- ifelse(d$Loading > 0, "positive", "negative")
  ggplot2::ggplot(d, ggplot2::aes(Feature, Loading, fill = Direction)) +
    ggplot2::geom_col(width = 0.6, alpha = 0.85) +
    ggplot2::scale_fill_manual(values = c(positive = "#D7191C",
                                          negative = "#2C7BB6")) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = sprintf("Feature contributions to %s", pc),
                  x = "Feature", y = sprintf("%s loading", pc)) +
    .pca_base_theme() +
    ggplot2::theme(legend.position = "none")
}

#' Heatmap of trait x PC correlations (quantitative traits).
plot_trait_pc_heatmap <- function(corr_long, max_pc = 5) {
  corr_long <- corr_long[corr_long$PC %in% paste0("PC", seq_len(max_pc)), ]
  corr_long$PC <- factor(corr_long$PC, levels = paste0("PC", seq_len(max_pc)))
  ggplot2::ggplot(corr_long, ggplot2::aes(PC, trait, fill = r)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", r)), size = 3.5) +
    ggplot2::scale_fill_gradient2(low = "#2C7BB6", mid = "white",
                                  high = "#D7191C", midpoint = 0,
                                  limits = c(-1, 1), name = "Pearson r") +
    ggplot2::labs(title = "Correlation of quantitative traits with PCs",
                  x = "", y = "") +
    .pca_base_theme()
}

#' Violin + box of a single PC's scores across a categorical trait,
#' with the one-way ANOVA p-value in the subtitle.
plot_violin_pc <- function(scores, group, group_name, pc = "PC1") {
  d <- data.frame(Group = factor(group), Score = scores[[pc]])
  fit <- stats::aov(Score ~ Group, data = d)
  p   <- summary(fit)[[1]][["Pr(>F)"]][1]
  ggplot2::ggplot(d, ggplot2::aes(Group, Score, fill = Group)) +
    ggplot2::geom_violin(alpha = 0.4) +
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.05, size = 1.5) +
    ggplot2::scale_fill_manual(values = PCA_PALETTE) +
    ggplot2::labs(title = sprintf("%s scores by %s", pc, group_name),
                  subtitle = sprintf("One-way ANOVA p = %s", formatC(p, format = "e", digits = 2)),
                  x = group_name, y = sprintf("%s score", pc)) +
    .pca_base_theme() +
    ggplot2::theme(legend.position = "none")
}
