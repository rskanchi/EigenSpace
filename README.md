# EigenSpace

A standalone R/Shiny app that takes expression data and optional metadata, runs PCA, and produces a PowerPoint deck, 600-dpi
figures, CSV tables, and a manuscript-ready paragraph.

## Structure

```
EigenSpace/
├── app.R              # Shiny app
├── run_pca.R          # R script to run PCA without GUI/Shiny
├── R/
│   ├── pca_setup.R   # read & align inputs, preprocess, prcomp, result tables, trait tests
│   ├── pca_plots.R    # generate plots using ggplot2
│   ├── pca_summary.R  # manuscript paragraph & slide bullets
│   └── pca_export.R   # build run, save 600 dpi PNGs + CSVs + officer pptx, zip
└── data/
    ├── demo_expression.csv   # 6 genes x 8 samples (example, demo)
    └── demo_metadata.csv     # Group (categorical), Batch (categorical), RIN (quantitative)
```

The setup, plots, summary, and export functions have no Shiny dependency, so you can `source()` them and call `build_run()` / `export_run()` from any script.

## Install

```r
install.packages(c("shiny", "ggplot2", "officer"))
```

`prcomp` (the core computation) is part of base R `stats`; nothing to install.

## Run the app

From inside the `EigenSpace/` folder:

```r
shiny::runApp()
```

Click **Load demo data** to try it immediately, or upload your own files.

## Run from a script (no GUI required)

Edit the settings block at the top of `run_pca.R`, then:

```bash
Rscript run_pca.R
```

## Input format

- **Expression CSV**: first column = IDs. Use the orientation toggle to say whether rows are features or samples. (Demo is features × samples.)
- **Metadata CSV** (optional): first column = sample IDs matching the expression sample IDs exactly. Remaining columns are traits; the app auto-detects categorical vs quantitative and lets you override.

## Outputs (one zip, named by your run label)

```
<project_name>/
├── <project_name>_PCA.pptx   # title, methods (bullets), results, per-trait slides
├── figures/              # every plot as a 600-dpi PNG
├── tables/               # variance, scores, loadings, ANOVA, trait–PC correlations (CSV)
└── summary.txt           # manuscript paragraph + bullets
```

If you would like an overview of PCA based on the demo dataset of 8 samples and 6 genes, see [Understanding PCA](docs/UnderstandingPCA.md).
 
## Notes

- Core PCA: `stats::prcomp()`; Visualization: `ggplot2`; Deck: `officer`. 
- `log2` uses `log2(x + 1)` so zeros stay finite
- Categorical traits: grouped score plots (with/without 95% ellipses) + per-PC violin/box with one-way ANOVA. 
- Quantitative traits: gradient score plots + Pearson trait×PC correlation table and heatmap

