# pca_setup.R includes PCA computation (no Shiny, no plotting)

# 1. Read & align inputs
# Read an expression file and return a numeric matrix in genes x samples form
# Accepts row names OR a first ID column; orientation informs the way file is organized
# the output matrix returned is genes(features) x samples
read_expression <- function(path, orientation = c("features_x_samples","samples_x_features")) {
  orientation <- match.arg(orientation)
  df <- utils::read.csv(path, header = TRUE, check.names = FALSE,
                        stringsAsFactors = FALSE)
  # First column is treated as IDs (row labels).
  ids <- as.character(df[[1]])
  mat <- as.matrix(df[, -1, drop = FALSE])
  rownames(mat) <- ids
  storage.mode(mat) <- "double"
  if (orientation == "samples_x_features") mat <- t(mat)  # -> features x samples
  mat
} # end of read_expression

# Read a metadata file. First column (or row names) must be sample IDs
read_metadata <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path)) return(NULL)
  df <- utils::read.csv(path, header = TRUE, check.names = FALSE,
                        stringsAsFactors = FALSE)
  ids <- as.character(df[[1]])
  df <- df[, -1, drop = FALSE]
  rownames(df) <- ids
  df
}

# Align expression (genes x samples) with metadata by shared sample IDs
# Returns the objects for analysis and a report
align_inputs <- function(expr_genes_x_samples, metadata = NULL) {
  expr_samples <- colnames(expr_genes_x_samples)
  if (is.null(metadata)) {
    common <- expr_samples
  } else {
    common <- intersect(expr_samples, rownames(metadata))
  }
  if (length(common) < 3)
    stop("Fewer than 3 samples are shared between the files (or present in ",
         "the expression data). PCA needs at least 3 samples. Check that the ",
         "sample IDs match exactly (watch for stray spaces or naming like ",
         "'H1' vs 'H_1').")

  expr_aligned <- expr_genes_x_samples[, common, drop = FALSE]
  meta_aligned <- if (is.null(metadata)) NULL else metadata[common, , drop = FALSE]

  report <- list(
    n_expr_samples = length(expr_samples),
    n_meta_samples = if (is.null(metadata)) NA_integer_ else nrow(metadata),
    n_matched      = length(common),
    n_features     = nrow(expr_genes_x_samples),
    dropped_from_expr = setdiff(expr_samples, common),
    dropped_from_meta = if (is.null(metadata)) character(0)
                        else setdiff(rownames(metadata), common)
  )
  list(expr = expr_aligned, metadata = meta_aligned, report = report)
} # end of align_inputs

# Preprocess & run PCA

# Optional log2(x + 1) to a genes x samples matrix
log2_transform <- function(expr_genes_x_samples, do_log2 = FALSE) {
  if (!do_log2) return(expr_genes_x_samples)
  if (any(expr_genes_x_samples < 0, na.rm = TRUE))
    warning("Negative values present; log2(x + 1) of a negative number is NaN. ",
            "Check whether these data are already log-transformed.")
  log2(expr_genes_x_samples + 1)
} # end of log2_transform

# Run PCA
# Input = genes x samples expr matrix; prcomp requires samples as rows
# center/scale are passed as arguemnts
# Returns the prcomp object plus tidy result tables
run_pca <- function(expr_genes_x_samples, center = TRUE, scale = FALSE, do_log2 = FALSE) {
  expr <- log2_transform(expr_genes_x_samples, do_log2)
  X <- t(expr)                                  # samples x features
  # Drop zero-variance features when scaling (prcomp errors otherwise).
  if (scale) {
    keep <- apply(X, 2, stats::sd) > 0
    if (any(!keep))
      warning(sum(!keep), " feature(s) with zero variance dropped before ",
              "scaling: ", paste(colnames(X)[!keep], collapse = ", "))
    X <- X[, keep, drop = FALSE]
  }
  pca <- stats::prcomp(X, center = center, scale. = scale)

  # Variance table
  sdev <- pca$sdev
  eig  <- sdev^2
  pct  <- eig / sum(eig) * 100
  variance <- data.frame(
    PC                = paste0("PC", seq_along(sdev)),
    stdev             = round(sdev, 4),
    eigenvalue        = round(eig, 4),
    percent_variance  = round(pct, 2),
    cumulative_variance = round(cumsum(pct), 2),
    stringsAsFactors  = FALSE
  )

  scores   <- as.data.frame(pca$x)                 # samples x PC
  loadings <- as.data.frame(pca$rotation)          # features x PC

  list(prcomp = pca, variance = variance,
       scores = scores, loadings = loadings,
       settings = list(center = center, scale = scale, log2 = do_log2))
} # end of run_pca

# Metadata trait classification & tests

# Auto-classify each metadata column as "categorical" or "quantitative"
# Numeric columns with few distinct values are treated as categorical
classify_traits <- function(metadata, max_levels_numeric = 5L) {
  if (is.null(metadata) || ncol(metadata) == 0)
    return(data.frame(trait = character(0), type = character(0),
                      n_levels = integer(0), stringsAsFactors = FALSE))
  out <- lapply(names(metadata), function(nm) {
    x <- metadata[[nm]]
    n_levels <- length(unique(stats::na.omit(x)))
    type <- if (is.numeric(x) && n_levels > max_levels_numeric)
      "quantitative" else "categorical"
    data.frame(trait = nm, type = type, n_levels = n_levels,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
} # end of classify_traits

# One-way ANOVA of each PC against a categorical trait
# Returns PC, F, p for the requested PCs
anova_pc_by_group <- function(scores, group, pcs = c("PC1", "PC2")) {
  pcs <- intersect(pcs, colnames(scores))
  group <- factor(group)
  res <- lapply(pcs, function(pc) {
    fit <- stats::aov(scores[[pc]] ~ group)
    s <- summary(fit)[[1]]
    data.frame(PC = pc, F_value = round(s[["F value"]][1], 3),
               p_value = signif(s[["Pr(>F)"]][1], 3),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
} # end of anova_pc_by_group

# Pearson correlation of each PC with a quantitative trait
cor_pc_with_trait <- function(scores, trait_values, pcs = colnames(scores)) {
  res <- lapply(pcs, function(pc) {
    ct <- suppressWarnings(stats::cor.test(scores[[pc]], trait_values,
                                           method = "pearson"))
    data.frame(PC = pc, r = round(unname(ct$estimate), 3),
               p_value = signif(ct$p.value, 3), stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}

# Full trait x PC correlation table (quantitative traits) in long format for a heatmap and for export
trait_pc_correlations <- function(scores, metadata, quant_traits, pcs = colnames(scores)) {
  if (length(quant_traits) == 0) return(NULL)
  rows <- lapply(quant_traits, function(tr) {
    tab <- cor_pc_with_trait(scores, as.numeric(metadata[[tr]]), pcs)
    tab$trait <- tr
    tab
  })
  out <- do.call(rbind, rows)
  out[, c("trait", "PC", "r", "p_value")]
}
