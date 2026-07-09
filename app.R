# =============================================================================
# app.R  --  Standalone Shiny app for PCA. Thin shell over the engine in R/.
# Run with:  shiny::runApp()   (from this folder)
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
})
# officer is needed only when exporting; loaded lazily in pca_export.R.
src <- file.path("R", c("pca_engine.R", "pca_plots.R",
                        "pca_summary.R", "pca_export.R"))
invisible(lapply(src, source))

DEMO_EXPR <- "data/demo_expression.csv"
DEMO_META <- "data/demo_metadata.csv"

# ---- UI ---------------------------------------------------------------------
ui <- fluidPage(
  tags$head(tags$style(HTML(
    "#paragraph { white-space: pre-wrap; word-wrap: break-word; }"
  ))),
  titlePanel("Principal Component Analysis"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("1. Data"),
      fileInput("expr_file", "Expression file (CSV)", accept = ".csv"),
      fileInput("meta_file", "Metadata file (CSV) - optional", accept = ".csv"),
      helpText(
        "Sample IDs must be present in both files - as row names (first column) ",
        "in the expression file, and as the first column of the metadata file - ",
        "and must match exactly. Only samples present in both are analyzed."),
      actionButton("use_demo", "Load demo data", class = "btn-info btn-sm"),
      tags$hr(),
      h4("2. Orientation & preprocessing"),
      radioButtons("orientation", "Expression layout",
                   c("Features (rows) x Samples (columns)" = "features_x_samples",
                     "Samples (rows) x Features (columns)" = "samples_x_features")),
      checkboxInput("do_log2", "Apply log2(x + 1) transform", FALSE),
      checkboxInput("center", "Center (subtract mean)", TRUE),
      checkboxInput("scale", "Scale to unit variance", FALSE),
      helpText(
        "For expression data a common choice is log2 + center without scaling; ",
        "scaling is more appropriate when features remain on very different ",
        "scales. See the 'Help / Learn' tab for the full explanation."),
      tags$hr(),
      h4("3. Run"),
      textInput("run_name", "project name (folder name to save output)", "pca_run"),
      uiOutput("trait_ui"),
      actionButton("run", "Run PCA", class = "btn-primary"),
      tags$br(), tags$br(),
      downloadButton("download_zip", "Download results (.zip)")
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        id = "tabs",
        tabPanel("Status", verbatimTextOutput("status")),
        tabPanel("Variance", tableOutput("variance_tbl"),
                 plotOutput("scree_plot", height = 350, width = "60%")),
        tabPanel("Scores", plotOutput("score_plot", height = 450, width = "60%")),
        tabPanel("Loadings", plotOutput("loading_plot", height = 450, width = "75%")),
        tabPanel("Trait plots", uiOutput("trait_plot_ui")),
        tabPanel("Summary", h4("Summary, use as a starting point..."),
                 verbatimTextOutput("paragraph")),
        tabPanel("Help / Learn", uiOutput("help_ui"))
      )
    ) # end of mainPanel
  )
)

# ---- Server -----------------------------------------------------------------
server <- function(input, output, session) {

  `%||%` <- function(a, b) if (is.null(a)) b else a

  rv <- reactiveValues(expr_path = NULL, meta_path = NULL, run = NULL,
                       zip = NULL)

  observeEvent(input$expr_file, rv$expr_path <- input$expr_file$datapath)
  observeEvent(input$meta_file, rv$meta_path <- input$meta_file$datapath)
  observeEvent(input$use_demo, {
    rv$expr_path <- DEMO_EXPR; rv$meta_path <- DEMO_META
    updateCheckboxInput(session, "scale", value = TRUE)  # demo is raw + un-logged
    showNotification("Demo data loaded (genes x samples; Group, Batch, RIN).",
                     type = "message")
  })

  metadata <- reactive({
    if (is.null(rv$meta_path)) return(NULL)
    read_metadata(rv$meta_path)
  })

  # Dynamic trait-type override UI, one selector per metadata column.
  output$trait_ui <- renderUI({
    md <- metadata(); if (is.null(md)) return(NULL)
    cls <- classify_traits(md)
    tagList(
      tags$hr(), strong("Metadata trait types (auto-detected; override here):"),
      lapply(seq_len(nrow(cls)), function(i) {
        nm <- cls$trait[i]
        selectInput(paste0("trait_", nm), nm,
                    choices = c("categorical", "quantitative", "ignore"),
                    selected = cls$type[i])
      })
    )
  })

  trait_overrides <- reactive({
    md <- metadata(); if (is.null(md)) return(NULL)
    nms <- names(md)
    vals <- vapply(nms, function(nm) input[[paste0("trait_", nm)]] %||% NA_character_,
                   character(1))
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) NULL else vals
  })

  observeEvent(input$run, {
    req(rv$expr_path)
    withProgress(message = "Running PCA...", {
      expr <- read_expression(rv$expr_path, input$orientation)
      meta <- metadata()
      rv$run <- build_run(expr, meta, center = input$center, scale = input$scale,
                          do_log2 = input$do_log2,
                          trait_overrides = trait_overrides())
    })
    updateTabsetPanel(session, "tabs", selected = "Variance")
  })

  output$status <- renderText({
    if (is.null(rv$run)) {
      if (is.null(rv$expr_path))
        "Upload an expression file (and optionally metadata), or click 'Load demo data', then 'Run PCA'."
      else "Data loaded. Set options and click 'Run PCA'."
    } else {
      rep <- rv$run$aligned$report
      paste0(
        "Run complete.\n",
        "Samples analyzed: ", rep$n_matched, "\n",
        "Features: ", rep$n_features, "\n",
        if (!is.na(rep$n_meta_samples))
          paste0("Matched ", rep$n_matched, " of ", rep$n_expr_samples,
                 " expression samples with metadata.\n") else "",
        if (length(rep$dropped_from_expr))
          paste0("Dropped (no metadata match): ",
                 paste(rep$dropped_from_expr, collapse = ", "), "\n") else "")
    }
  })

  output$variance_tbl <- renderTable({ req(rv$run); rv$run$pca$variance })
  output$scree_plot    <- renderPlot({ req(rv$run); rv$run$plots$scree })
  output$score_plot    <- renderPlot({ req(rv$run); rv$run$plots$scores })
  output$loading_plot  <- renderPlot({ req(rv$run); rv$run$plots$loadings_PC1 })
  output$paragraph     <- renderText({ req(rv$run); rv$run$paragraph })

  # Trait plots: render every trait-specific plot stacked.
  output$trait_plot_ui <- renderUI({
    req(rv$run)
    nms <- grep("^(scores_|violin_|trait_pc)", names(rv$run$plots), value = TRUE)
    if (length(nms) == 0) return(helpText("No metadata traits to plot."))
    lapply(nms, function(nm) {
      w <- if (grepl("^scores_", nm)) "55%" else "75%"
      plotOutput(paste0("tp_", nm), height = 380, width = w)
    })
  })
  observe({
    req(rv$run)
    nms <- grep("^(scores_|violin_|trait_pc)", names(rv$run$plots), value = TRUE)
    lapply(nms, function(nm) {
      output[[paste0("tp_", nm)]] <- renderPlot(rv$run$plots[[nm]])
    })
  })

  output$download_zip <- downloadHandler(
    filename = function() paste0(gsub("[^A-Za-z0-9_.-]+", "_", input$run_name), ".zip"),
    content = function(file) {
      req(rv$run)
      tmp <- tempfile("pca_out"); dir.create(tmp)
      zip_path <- export_run(rv$run, tmp, input$run_name)
      file.copy(zip_path, file, overwrite = TRUE)
    }
  )

  output$help_ui <- renderUI({
    HTML(paste0(
      "<h3>Choosing log2, centering, and scaling</h3>",
      "<p>These are three independent options that do different jobs, so set ",
      "them deliberately rather than as a package. <b>Centering</b> (subtracting ",
      "each feature's mean) should almost always be on: PCA describes variation ",
      "<i>around the mean</i>, and without centering the first component can ",
      "simply capture the overall offset of the data rather than its structure. ",
      "<b>log2 transformation</b> compresses dynamic range - it turns large ",
      "multiplicative differences into smaller additive ones and reduces ",
      "right-skew - so highly expressed features stop dominating purely because ",
      "of their magnitude. <b>Scaling</b> (dividing each feature by its standard ",
      "deviation) forces every feature to contribute equally regardless of ",
      "variability.</p>",
      "<p>For expression data a common, well-justified choice is <b>log2 + ",
      "centering, without scaling</b>: once data are on the log scale, much of ",
      "the magnitude-dominance problem is already handled, and the features that ",
      "remain highly variable are usually the biologically informative ones you ",
      "<i>want</i> to weight more. Scaling here would lift flat, near-noise ",
      "features (for example, genes near the detection floor) up to the same ",
      "weight as strongly differential ones, amplifying noise. Scaling is the ",
      "better choice in the opposite situation - when features remain on very ",
      "different scales and have <i>not</i> been brought together by a transform ",
      "(the textbook case for scaling raw, untransformed data).</p>",
      "<h3>About the demo data</h3>",
      "<p>The built-in demo is a small synthetic gene-expression dataset: 6 genes ",
      "x 8 samples (4 healthy, 4 cancer), with metadata columns for Group ",
      "(categorical), Batch (categorical), and RIN (quantitative). Gene1 and ",
      "Gene4 are high in cancer; Gene2 and Gene6 are high in healthy; Gene3 and ",
      "Gene5 differ little. PC1 separates the two groups and explains the ",
      "majority of the variance. It is raw and un-logged, so the demo runs with ",
      "centering + scaling.</p>"))
  })
}

shinyApp(ui, server)
