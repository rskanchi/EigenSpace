# pca_summary.R:  summary as a starting point for a writeup and the bullet text for the pptx slides

# Short methods sentence
methods_sentence <- function(settings, n_samples, n_features) {
  pre <- c()
  if (isTRUE(settings$log2))   pre <- c(pre, "log2(x + 1)-transformed")
  if (isTRUE(settings$center)) pre <- c(pre, "mean-centered")
  if (isTRUE(settings$scale))  pre <- c(pre, "scaled to unit variance")
  pre_txt <- if (length(pre)) paste0(" Data were ", paste(pre, collapse = ", "),
                                     " prior to decomposition.") else ""
  sprintf(paste0(
    "Principal component analysis was performed on %d samples and %d features ",
    "using the prcomp() function from the stats package in %s.%s ",
    "Visualizations were produced with the ggplot2 package."),
    n_samples, n_features, R.version.string, pre_txt)
} # end of methods_sentence

# Results sentence summarizing variance captured by the leading PCs
results_sentence <- function(variance) {
  pc1 <- variance$percent_variance[1]
  pc2 <- if (nrow(variance) >= 2) variance$percent_variance[2] else NA
  cum2 <- if (nrow(variance) >= 2) variance$cumulative_variance[2] else pc1
  if (is.na(pc2))
    sprintf("PC1 explained %.1f%% of the total variance.", pc1)
  else
    sprintf(paste0("PC1 and PC2 explained %.1f%% and %.1f%% of the total ",
                   "variance, respectively (cumulatively %.1f%%)."),
            pc1, pc2, cum2)
} # end of results_sentence

# Optional trait-association sentences
trait_sentences <- function(anova_tabs = NULL, corr_long = NULL) {
  s <- character(0)
  alpha <- 0.05
  if (!is.null(anova_tabs)) {
    for (nm in names(anova_tabs)) {
      t <- anova_tabs[[nm]]
      r <- t[which.min(t$p_value), ]
      verb <- if (r$p_value < alpha) "differed significantly"
      else "did not differ significantly"
      s <- c(s, sprintf(
        "Scores on %s %s by %s (one-way ANOVA p = %s).",
        r$PC, verb, nm, formatC(r$p_value, format = "e", digits = 2)))
    }
  } # end of anova_tabs
  
  if (!is.null(corr_long) && nrow(corr_long) > 0) {
    best <- corr_long[which.min(corr_long$p_value), ]
    assoc <- if (best$p_value < alpha) "was significantly associated with"
    else "showed no significant association with"
    s <- c(s, sprintf(
      "Among quantitative traits, %s %s %s (Pearson r = %.2f, p = %s).",
      best$trait, assoc, best$PC, best$r, formatC(best$p_value, format = "e", digits = 2)))
  } # end of corr
  s
} # end of trait_sentences

# Assemble the full manuscript paragraph (methods + results + associations).
manuscript_paragraph <- function(settings, report, variance, anova_tabs = NULL, corr_long = NULL) {
  paste(c(
    methods_sentence(settings, report$n_matched, report$n_features),
    results_sentence(variance),
    trait_sentences(anova_tabs, corr_long)
  ), collapse = " ")
} # end of manuscript_paragraph

# Bullet lists for the slides (named list of character vectors).
slide_bullets <- function(settings, report, variance, anova_tabs = NULL, corr_long = NULL) {
  pre <- c(
    sprintf("log2(x + 1) transform: %s", if (settings$log2) "yes" else "no"),
    sprintf("Centering: %s", if (settings$center) "yes" else "no"),
    sprintf("Scaling to unit variance: %s", if (settings$scale) "yes" else "no"))
  methods <- c(
    sprintf("Samples analyzed: %d", report$n_matched),
    sprintf("Features (variables): %d", report$n_features),
    if (!is.na(report$n_meta_samples))
      sprintf("Samples matched between expression and metadata: %d of %d",
              report$n_matched, report$n_expr_samples) else NULL,
    pre,
    "Core computation: stats::prcomp() (base R)",
    "Visualization: ggplot2")
  results <- c(
    sprintf("PC1 explains %.1f%% of variance", variance$percent_variance[1]),
    if (nrow(variance) >= 2)
      sprintf("PC2 explains %.1f%% (cumulative %.1f%%)",
              variance$percent_variance[2], variance$cumulative_variance[2]) else NULL,
    if (nrow(variance) >= 3)
      sprintf("PC1-PC3 cumulative: %.1f%%", variance$cumulative_variance[3]) else NULL)
  list(methods = methods, results = results,
       associations = trait_sentences(anova_tabs, corr_long))
}
