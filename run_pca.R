# =============================================================================
# run_pca.R  --  Run the full pipeline from a script (no Shiny needed).
#
# Edit the values below and run:  Rscript run_pca.R
# =============================================================================

# --- Load the setup ---------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(officer)
})
src <- file.path("R", c("pca_setup.R", "pca_plots.R",
                        "pca_summary.R", "pca_export.R"))
invisible(lapply(src, source))

# --- Settings (edit these) ---------------------------------------------------
expression_file <- "data/demo_expression.csv"
metadata_file   <- "data/demo_metadata.csv"   # or NULL for none
orientation     <- "features_x_samples"        # or "samples_x_features"
do_log2 <- FALSE
center  <- TRUE
scale   <- TRUE                                  # demo data is raw + un-logged
project_name <- "demo_run"
outdir   <- "output"

# Optional: override auto-detected trait types, e.g.
# trait_overrides <- c(Batch = "categorical", RIN = "quantitative")
trait_overrides <- NULL

# --- Run ---------------------------------------------------------------------
expr <- read_expression(expression_file, orientation)
meta <- read_metadata(metadata_file)

run <- build_run(expr, meta, center = center, scale = scale, do_log2 = do_log2,
                 trait_overrides = trait_overrides)

cat(run$paragraph, "\n\n")
zip_path <- export_run(run, outdir, project_name)
cat("All artifacts written. Zip:", zip_path, "\n")
