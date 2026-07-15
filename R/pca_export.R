# pca_export.R: export results, 600-dpi figures, CSV tables, 
# a bulleted PowerPoint (officer), a manuscript paragraph, all zipped.
# Depends on the engine, plots, and summary files being sourced first.

DPI <- 600

# Decide which traits to use and of what type, honoring user overrides.
# overrides: optional named character vector, e.g.
#   c(Group = "categorical", RIN = "quantitative", Batch = "ignore")
resolve_traits <- function(metadata, overrides = NULL) {
  if (is.null(metadata)) return(list(categorical = character(0),
                                     quantitative = character(0)))
  cls <- classify_traits(metadata)
  types <- setNames(cls$type, cls$trait)
  if (!is.null(overrides)) {
    for (nm in names(overrides)) if (nm %in% names(types)) types[nm] <- overrides[[nm]]
  }
  list(categorical  = names(types)[types == "categorical"],
       quantitative = names(types)[types == "quantitative"])
}

# Run the whole analysis and return every object needed for display & export.
build_run <- function(expr_genes_x_samples, metadata = NULL,
                      center = TRUE, scale = FALSE, do_log2 = FALSE,
                      trait_overrides = NULL, n_dim = 10) {

  aligned <- align_inputs(expr_genes_x_samples, metadata)
  pca <- run_pca(aligned$expr, center = center, scale = scale, do_log2 = do_log2)
  traits <- resolve_traits(aligned$metadata, trait_overrides)

  # Plots (named list of ggplots)
  plots <- list(
    scree      = plot_scree(pca$variance, n_dim),
    cumulative = plot_cumulative(pca$variance, n_dim),
    scores     = plot_scores_plain(pca$scores, pca$variance),
    loadings_PC1 = plot_loadings(pca$loadings, "PC1"),
    loadings_PC2 = plot_loadings(pca$loadings, "PC2")
  )

  anova_tabs <- list()
  for (g in traits$categorical) {
    grp <- aligned$metadata[[g]]
    plots[[paste0("scores_", g)]]          <- plot_scores_categorical(
      pca$scores, pca$variance, grp, g, ellipse = FALSE)
    plots[[paste0("scores_", g, "_ellipse")]] <- plot_scores_categorical(
      pca$scores, pca$variance, grp, g, ellipse = TRUE)
    plots[[paste0("violin_", g, "_PC1")]]  <- plot_violin_pc(pca$scores, grp, g, "PC1")
    plots[[paste0("violin_", g, "_PC2")]]  <- plot_violin_pc(pca$scores, grp, g, "PC2")
    anova_tabs[[g]] <- anova_pc_by_group(pca$scores, grp)
  }

  corr_long <- NULL
  if (length(traits$quantitative) > 0) {
    corr_long <- trait_pc_correlations(pca$scores, aligned$metadata,
                                       traits$quantitative)
    plots[["trait_pc_heatmap"]] <- plot_trait_pc_heatmap(corr_long)
    for (q in traits$quantitative)
      plots[[paste0("scores_", q, "_gradient")]] <- plot_scores_gradient(
        pca$scores, pca$variance, aligned$metadata[[q]], q)
  }

  bullets   <- slide_bullets(pca$settings, aligned$report, pca$variance,
                             anova_tabs, corr_long)
  paragraph <- manuscript_paragraph(pca$settings, aligned$report, pca$variance,
                                    anova_tabs, corr_long)

  list(pca = pca, aligned = aligned, traits = traits, plots = plots,
       anova_tabs = anova_tabs, corr_long = corr_long,
       bullets = bullets, paragraph = paragraph)
}

# ---- Outputs

write_tables <- function(run, dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(run$pca$variance, file.path(dir, "variance_explained.csv"),
                   row.names = FALSE)
  utils::write.csv(data.frame(Sample = rownames(run$pca$scores), run$pca$scores),
                   file.path(dir, "pc_scores.csv"), row.names = FALSE)
  utils::write.csv(data.frame(Feature = rownames(run$pca$loadings), run$pca$loadings),
                   file.path(dir, "loadings.csv"), row.names = FALSE)
  for (nm in names(run$anova_tabs))
    utils::write.csv(run$anova_tabs[[nm]],
                     file.path(dir, sprintf("anova_%s.csv", nm)), row.names = FALSE)
  if (!is.null(run$corr_long))
    utils::write.csv(run$corr_long, file.path(dir, "trait_pc_correlations.csv"),
                     row.names = FALSE)
}

write_figures <- function(run, dir, width = 7, height = 5, dpi = DPI) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  for (nm in names(run$plots)) {
    ggplot2::ggsave(file.path(dir, paste0(nm, ".png")), run$plots[[nm]],
                    width = width, height = height, dpi = dpi, bg = "white")
  }
}

write_summary <- function(run, path) {
  writeLines(c("PCA ANALYSIS SUMMARY", "====================", "",
               "Manuscript paragraph:", run$paragraph, "",
               "Methods bullets:", paste0("  - ", run$bullets$methods), "",
               "Results bullets:", paste0("  - ", run$bullets$results),
               if (length(run$bullets$associations))
                 c("", "Trait associations:",
                   paste0("  - ", run$bullets$associations))),
             path)
}

# Build the PowerPoint with officer; fig_dir must already contain the PNGs.
write_pptx <- function(run, path, fig_dir, project_name = "PCA analysis") {
  if (!requireNamespace("officer", quietly = TRUE))
    stop("Package 'officer' is required to write the PowerPoint.")
  library(officer)

  add_fig <- function(doc, title, png, bullets = NULL) {
    doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
    doc <- ph_with(doc, value = title, location = ph_location_type(type = "title"))
    doc <- ph_with(doc, value = external_img(png, width = 6, height = 4.3),
                   location = ph_location(left = 0.4, top = 1.3))
    if (!is.null(bullets)) {
      props <- lapply(seq_along(bullets), function(i) fp_text(font.size = 16, color = "#222222"))
      bl <- unordered_list(str_list = bullets, level_list = rep(1, length(bullets)), style = props)
      doc <- ph_with(doc, value = bl,
                     location = ph_location(left = 6.6, top = 1.3,
                                            width = 3.2, height = 5))
    }
    doc
  }
  add_text <- function(doc, title, bullets) {
    doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
    doc <- ph_with(doc, value = title, location = ph_location_type(type = "title"))
    bl <- unordered_list(str_list = bullets, level_list = rep(1, length(bullets)))
    doc <- ph_with(doc, value = bl, location = ph_location_type(type = "body"))
    doc
  }
  fp <- function(nm) file.path(fig_dir, paste0(nm, ".png"))

  doc <- read_pptx()
  # Title slide
  doc <- add_slide(doc, layout = "Title Slide", master = "Office Theme")
  doc <- ph_with(doc, value = project_name, location = ph_location_type(type = "ctrTitle"))
  doc <- ph_with(doc, value = paste("Principal Component Analysis -",
                                    format(Sys.Date())),
                 location = ph_location_type(type = "subTitle"))
  # Methods + results text slides
  doc <- add_text(doc, "Methods", run$bullets$methods)
  doc <- add_fig(doc, "Variance explained", fp("scree"), run$bullets$results)
  doc <- add_fig(doc, "Cumulative variance", fp("cumulative"))
  doc <- add_fig(doc, "PCA score plot", fp("scores"))
  doc <- add_fig(doc, "Feature contributions to PC1", fp("loadings_PC1"))
  doc <- add_fig(doc, "Feature contributions to PC2", fp("loadings_PC2"))
  # Categorical trait slides
  for (g in run$traits$categorical) {
    doc <- add_fig(doc, paste("Scores by", g), fp(paste0("scores_", g)))
    doc <- add_fig(doc, paste("Scores by", g, "(95% ellipses)"),
                   fp(paste0("scores_", g, "_ellipse")))
    doc <- add_fig(doc, paste(g, "- PC1 by group"), fp(paste0("violin_", g, "_PC1")))
  }
  # Quantitative trait slides
  if (length(run$traits$quantitative) > 0) {
    doc <- add_fig(doc, "Trait-PC correlations", fp("trait_pc_heatmap"))
    for (q in run$traits$quantitative)
      doc <- add_fig(doc, paste("Scores coloured by", q),
                     fp(paste0("scores_", q, "_gradient")))
  }
  # Manuscript paragraph slide
  doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
  doc <- ph_with(doc, value = "Summary",
                 location = ph_location_type(type = "title"))
  para <- block_list(fpar(ftext(run$paragraph, fp_text(font.size = 14))))
  doc <- ph_with(doc, value = para, location = ph_location_type(type = "body"))
  print(doc, target = path)
  invisible(path)
} # end of write_pptx

# Top-level: write everything for a run into <outdir>/<project_name>/... and zip it
# Returns the path to the zip file
export_run <- function(run, outdir, project_name = "pca_run") {
  safe <- gsub("[^A-Za-z0-9_.-]+", "_", project_name)
  root <- file.path(outdir, safe)
  unlink(root, recursive = TRUE); dir.create(root, recursive = TRUE)

  fig_dir <- file.path(root, "figures")
  tab_dir <- file.path(root, "tables")
  write_figures(run, fig_dir)
  write_tables(run, tab_dir)
  write_summary(run, file.path(root, "summary.txt"))
  write_pptx(run, file.path(root, paste0(safe, "_PCA.pptx")), fig_dir, project_name)

  zip_path <- file.path(outdir, paste0(safe, ".zip"))
  old <- setwd(outdir); on.exit(setwd(old))
  utils::zip(zipfile = basename(zip_path), files = safe, flags = "-r9Xq")
  normalizePath(zip_path)
}
