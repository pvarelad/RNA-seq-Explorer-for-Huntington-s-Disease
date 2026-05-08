library(shiny)
library(tidyverse)
library(DT)
library(bslib)
library(pheatmap)   # clustered heatmap with color bar legend

# ─────────────────────────────────────────────
#  Helper: reusable file-type validator
#  Call inside any reactive() that loads a file.
#  Throws a shiny::validate() error if the file
#  extension is not in the allowed `types` vector.
# ─────────────────────────────────────────────

# handle input file size 
options(shiny.maxRequestSize = 500 * 1024^2)


validate_file <- function(file_input, types = c("csv", "tsv")) {
  ext <- tools::file_ext(file_input$name)
  validate(
    need(
      tolower(ext) %in% types,
      paste0("Invalid file type '.", ext, "'. Please upload a ",
             paste(toupper(types), collapse = " or "), " file.")
    )
  )
}

# ─────────────────────────────────────────────
#  UI
# ─────────────────────────────────────────────
ui <- fluidPage(
  theme = bslib::bs_theme(
    bootswatch   = "sandstone",
    primary      = "#2c7bb6",   # matches your existing plot blue
    danger       = "#c0392b",   # matches your existing plot red
    base_font    = bslib::font_google("Source Sans Pro"),
    heading_font = bslib::font_google("Source Sans Pro"),
    `font-size-base` = "0.95rem"
  ),
  
  p("Final Project - R for Biological Sciences (BF530) - Penelope Varela"),
  # App title
  titlePanel("Huntington's Disease Prefrontal Cortex — RNA-seq Explorer"),
  
  # App description shown at the top of every tab
  p(
    "This app explores post-mortem RNA-seq data from Huntington's Disease (HD)",
    "prefrontal cortex samples compared with neurologically healthy controls.",
    "Use the tabs below to explore sample metadata, filtered gene counts,",
    "differential expression results, and GSEA pathway enrichment."
  ),
  
  hr(),  # horizontal rule to separate header from tabs
  
  # ── Top-level tab set ──────────────────────
  tabsetPanel(
    
    # ══════════════════════════════════════════
    #  TAB 1 — SAMPLES
    #  Purpose: upload & explore sample metadata
    # ══════════════════════════════════════════
    tabPanel(
      "Sample Information",
      br(),
      sidebarLayout(
        
        sidebarPanel(
          h4("Sample Metadata"),
          p("Upload a CSV file containing sample metadata. The file should include
             columns for sample ID, condition (HD vs. Control), age, tissue,
             and RIN score."),
          
          # File upload for metadata CSV
          fileInput(
            inputId = "samples_file",
            label   = "Upload Metadata CSV",
            accept  = c(".csv", ".tsv")
          ),
          
          hr(),
          
          # Let the user pick which column to use for grouping in plots
          selectInput(
            inputId  = "samples_group_col",
            label    = "Color plots by:",
            choices  = c("condition", "tissue"),  # updated after file loads
            selected = "condition"
          )
        ),
        
        mainPanel(
          tabsetPanel(
            
            # Summary: per-column type + mean(sd) or distinct values
            tabPanel(
              "Summary",
              br(),
              p("A column-by-column summary of the metadata table. Numeric columns
                 show mean ± SD; categorical columns list their distinct values.
                 Row and column counts are shown above the table."),
              # Display total row/col counts as plain text above the table
              verbatimTextOutput("samples_dim_text"),
              br(),
              DTOutput("samples_summary_table")
            ),
            
            # Full metadata table
            tabPanel(
              "Table",
              br(),
              p("Full sample metadata table. Click column headers to sort."),
              DTOutput("samples_full_table")
            ),
            
            # Histograms for continuous variables + bar chart for condition
            tabPanel(
              "Plots",
              br(),
              p("Histograms of continuous metadata variables (Age, RIN, PMI) colored
                 by condition, and a bar chart showing sample counts per group.
                 Use the 'Color plots by' selector in the sidebar to change grouping."),
              h5("Age Distribution"),
              plotOutput("samples_hist_age"),
              br(),
              h5("RIN Distribution"),
              plotOutput("samples_hist_rin"),
              br(),
              h5("PMI Distribution"),
              plotOutput("samples_hist_pmi"),
              br(),
              h5("Sample Counts by Condition"),
              plotOutput("samples_bar_plot")
            )
          )
        )
      )
    ),  # end Samples tab
    
    # ══════════════════════════════════════════
    #  TAB 2 — COUNTS
    #  Purpose: upload normalized counts matrix,
    #           filter by variance, PCA, heatmap
    # ══════════════════════════════════════════
    tabPanel(
      "Counts Matrix",
      br(),
      sidebarLayout(
        
        sidebarPanel(
          h4("Normalized Counts Matrix"),
          p("Upload the normalized counts matrix (genes as rows, samples as columns).
             Use the controls below to filter genes by variance before plotting."),
          
          # File upload for normalized counts
          fileInput(
            inputId = "counts_file",
            label   = "Upload Normalized Counts CSV",
            accept  = c(".csv", ".tsv")
          ),
          
          hr(),
          h5("Gene Filtering"),
          
          # Variance percentile slider: keep only the top X% most variable genes.
          # E.g. 50% keeps only the half of genes with the highest variance.
          sliderInput(
            inputId = "var_percentile",
            label   = "Keep top X% of genes by variance",
            min     = 0, max = 100, value = 50, step = 5,
            post    = "%"
          ),
          
          hr(),
          h5("Heatmap Options"),
          
          # How many top-variance genes to include in the clustered heatmap
          numericInput(
            inputId = "top_genes_heatmap",
            label   = "Number of top genes for heatmap",
            value   = 50, min = 10, max = 500, step = 10
          ),
          
          # Toggle log2 transformation of counts before heatmap rendering
          checkboxInput(
            inputId = "heatmap_log",
            label   = "Log2-transform counts for heatmap",
            value   = TRUE
          ),
          
          hr(),
          h5("PCA Options"),
          
          # Let user choose which two PCs to plot on X and Y axes
          selectInput(
            inputId  = "pca_x",
            label    = "X axis principal component",
            choices  = paste0("PC", 1:10),
            selected = "PC1"
          ),
          selectInput(
            inputId  = "pca_y",
            label    = "Y axis principal component",
            choices  = paste0("PC", 1:10),
            selected = "PC2"
          )
        ),
        
        mainPanel(
          tabsetPanel(
            
            # ── Subtab 1: Filtering summary table ──
            tabPanel(
              "Filter Summary",
              br(),
              p("Summary of the effect of the current variance percentile filter.
                 Adjust the slider in the sidebar to see counts update in real time."),
              DTOutput("counts_summary_table")
            ),
            
            # ── Subtab 2: Diagnostic scatter plots ──
            # Genes passing the filter are plotted in a darker color;
            # genes that did not pass are shown in a lighter color.
            tabPanel(
              "Diagnostic Plots",
              br(),
              p("Each point is a gene. Dark color = passes the current variance filter;
                 light color = filtered out. Both plots use log10 scales to handle the
                 wide dynamic range typical of RNA-seq data."),
              h5("Median Count vs. Variance"),
              plotOutput("counts_scatter_var"),
              br(),
              h5("Median Count vs. Number of Zero-Count Samples"),
              plotOutput("counts_scatter_zeros")
            ),
            
            # ── Subtab 3: Clustered heatmap ──
            tabPanel(
              "Heatmap",
              br(),
              p("Clustered heatmap of the top N highest-variance genes after filtering
                 (set N in the sidebar). Rows and columns are both hierarchically
                 clustered. Toggle log2 transformation in the sidebar to improve
                 visualization of lowly expressed genes. The color bar legend shows
                 expression level."),
              plotOutput("counts_heatmap", height = "600px")
            ),
            
            # ── Subtab 4: PCA scatter plot ──
            tabPanel(
              "PCA",
              br(),
              p("Principal Component Analysis (PCA) projection of samples using the
                 variance-filtered gene set. Select which principal components to plot
                 on each axis using the sidebar dropdowns. The percentage of total
                 variance explained by each component is shown in the axis labels.
                 Points are colored by condition (HD vs. Control) if metadata is loaded
                 in the Samples tab; otherwise colored by sample name."),
              plotOutput("counts_pca_plot")
            )
          )
        )
      )
    ),  # end Counts tab
    
    # ══════════════════════════════════════════
    #  TAB 3 — DIFFERENTIAL EXPRESSION (DE)
    #  Purpose: upload pre-computed DESeq2
    #           results and explore them
    # ══════════════════════════════════════════
    tabPanel(
      "Differential Expression",
      br(),
      sidebarLayout(
        
        sidebarPanel(
          h4("Differential Expression Results"),
          p("Upload a CSV of pre-computed DESeq2 results. The file must contain
             columns: gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj."),
          
          # File upload for DESeq2 results
          fileInput(
            inputId = "de_file",
            label   = "Upload DESeq2 Results CSV",
            accept  = c(".csv", ".tsv")
          ),
          
          hr(),
          h5("Significance Thresholds"),
          
          # Adjusted p-value cutoff slider
          sliderInput(
            inputId = "padj_thresh",
            label   = "Adjusted p-value (padj) threshold",
            min = 0.001, max = 0.25, value = 0.05, step = 0.005
          ),
          
          # Log2 fold change cutoff slider
          sliderInput(
            inputId = "lfc_thresh",
            label   = "| log2 Fold Change | threshold",
            min = 0, max = 5, value = 1, step = 0.1
          ),
          
          hr(),
          h5("Volcano Plot Axes"),
          
          # Radio buttons let the user swap axes on the volcano plot
          radioButtons(
            inputId  = "volcano_x",
            label    = "X axis",
            choices  = c("log2FoldChange", "stat"),
            selected = "log2FoldChange"
          ),
          radioButtons(
            inputId  = "volcano_y",
            label    = "Y axis",
            choices  = c("padj", "pvalue"),
            selected = "padj"
          )
        ),
        
        mainPanel(
          tabsetPanel(
            
            # Summary: counts of up/down/NS genes
            tabPanel(
              "Summary",
              br(),
              p("Summary of significantly differentially expressed genes at the
                 selected padj and log2FC thresholds."),
              DTOutput("de_summary_table")
            ),
            
            # Full DE results table
            tabPanel(
              "Table",
              br(),
              p("Full DESeq2 results table. Significant genes are highlighted.
                 Click column headers to sort; use the search box to find a gene."),
              DTOutput("de_results_table")
            ),
            
            # Volcano and MA plots
            tabPanel(
              "Plots",
              br(),
              p("Volcano plot: each point is a gene. Red = significantly up-regulated,
                 blue = significantly down-regulated, grey = not significant at the
                 chosen thresholds."),
              plotOutput("de_volcano_plot"),
              br(),
              p("MA plot: log2 fold change versus mean expression (baseMean)."),
              plotOutput("de_ma_plot")
            )
          )
        )
      )
    ),  # end DE tab
    
    # ══════════════════════════════════════════
    #  TAB 4 — GSEA
    #  Purpose: upload or auto-load fgsea results
    #           and explore pathway enrichment
    # ══════════════════════════════════════════
    tabPanel(
      "GSEA",
      br(),
      sidebarLayout(
        sidebarPanel(
          h4("GSEA / Pathway Enrichment"),
          p("Upload a CSV of pre-computed fgsea results."),
          
          fileInput(
            inputId = "gsea_file",
            label   = "Upload GSEA Results CSV",
            accept  = c(".csv", ".tsv")
          )
        ),
        
        mainPanel(
          tabsetPanel(
            
            # ── Tab 1: Bar plot ──────────────────────────────
            tabPanel(
              "Top Pathways",
              br(),
              sidebarLayout(
                sidebarPanel(
                  sliderInput(
                    inputId = "top_n_pathways",
                    label   = "Number of top pathways (by adjusted p-value)",
                    min = 5, max = 50, value = 10, step = 5
                  ),
                  br(),
                  p("Click a bar in the plot to see the full table entry for that pathway."),
                  DTOutput("gsea_clicked_row")
                ),
                mainPanel(
                  plotOutput("gsea_bar_plot",
                             click = "gsea_bar_click",
                             height = "500px")
                )
              )
            ),
            
            # ── Tab 2: Filtered table ────────────────────────
            tabPanel(
              "Table",
              br(),
              sidebarLayout(
                sidebarPanel(
                  sliderInput(
                    inputId = "gsea_padj_thresh",
                    label   = "Adjusted p-value threshold",
                    min = 0.001, max = 1, value = 0.25, step = 0.005
                  ),
                  radioButtons(
                    inputId  = "nes_direction",
                    label    = "NES direction",
                    choices  = c("All", "Positive (activated)", "Negative (suppressed)"),
                    selected = "All"
                  ),
                  br(),
                  downloadButton("gsea_download", "Download Filtered Results")
                ),
                mainPanel(
                  DTOutput("gsea_results_table")
                )
              )
            ),
            
            # ── Tab 3: NES vs -log10(padj) scatter ──────────
            tabPanel(
              "Scatter Plot",
              br(),
              sidebarLayout(
                sidebarPanel(
                  sliderInput(
                    inputId = "gsea_scatter_padj",
                    label   = "Adjusted p-value threshold",
                    min = 0.001, max = 1, value = 0.25, step = 0.005
                  )
                ),
                mainPanel(
                  plotOutput("gsea_scatter_plot", height = "500px")
                )
              )
            )
            
          )
        )
      )
    )
    
  )  # end tabsetPanel
)  # end fluidPage


# ─────────────────────────────────────────────
#  SERVER
#  Each reactive() block corresponds to one
#  file input or derived dataset. Outputs are
#  wired to UI output IDs defined above.
# ─────────────────────────────────────────────
server <- function(input, output, session) {
  
  # ── SAMPLES ───────────────────────────────
  
  # Load and validate the metadata CSV
  metadata <- reactive({
    req(input$samples_file)                     # wait until file is uploaded
    validate_file(input$samples_file)           # check extension
    read_csv(input$samples_file$datapath,
             show_col_types = FALSE)
  })
  
  # Dimension text: "Number of rows: X   Number of columns: Y"
  output$samples_dim_text <- renderText({
    df <- metadata()
    paste0("Number of rows: ", nrow(df), "\nNumber of columns: ", ncol(df))
  })
  
  # Column-by-column summary table.
  # For numeric columns  → "mean (± sd)"
  # For character/factor → comma-separated list of distinct values (capped at 5)
  output$samples_summary_table <- renderDT({
    df <- metadata()
    
    summary_rows <- lapply(names(df), function(col) {
      x <- df[[col]]
      
      if (is.numeric(x)) {
        # Format as "mean (± sd)" rounded to 2 decimal places
        col_type   <- "numeric"
        col_values <- paste0(round(mean(x, na.rm = TRUE), 2),
                             " (\u00b1 ",                  # ± symbol
                             round(sd(x,   na.rm = TRUE), 2), ")")
      } else {
        # Treat everything else as categorical; list up to 5 distinct values
        col_type   <- "factor/character"
        uniq       <- unique(na.omit(as.character(x)))
        if (length(uniq) > 5) {
          # Truncate and note how many more values exist
          col_values <- paste0(paste(head(uniq, 5), collapse = ", "),
                               " ... (", length(uniq), " total)")
        } else {
          col_values <- paste(uniq, collapse = ", ")
        }
      }
      
      tibble(`Column Name` = col,
             `Type`        = col_type,
             `Mean (sd) or Distinct Values` = col_values)
    })
    
    bind_rows(summary_rows) %>%
      datatable(options = list(dom = "t", ordering = FALSE,
                               pageLength = 30),
                rownames = FALSE)
  })
  
  # Full metadata table
  output$samples_full_table <- renderDT({
    datatable(metadata(), rownames = FALSE,
              options = list(pageLength = 15, scrollX = TRUE))
  })
  
  # Helper: build a histogram for any continuous metadata column.
  # Filled and colored by the grouping variable chosen in the sidebar.
  make_hist <- function(col_name, x_label) {
    renderPlot({
      df <- metadata()
      validate(need(col_name %in% names(df),
                    paste0("Column '", col_name, "' not found in metadata.")))
      group_col <- input$samples_group_col
      ggplot(df, aes(x    = .data[[col_name]],
                     fill = .data[[group_col]])) +
        geom_histogram(bins = 15, alpha = 0.7, position = "identity",
                       color = "white") +
        labs(title = paste(x_label, "by", group_col),
             x     = x_label,
             y     = "Count",
             fill  = group_col) +
        theme_bw()
    })
  }
  
  # Wire each continuous variable to its own histogram output
  output$samples_hist_age <- make_hist("age", "Age (years)")
  output$samples_hist_rin <- make_hist("rin", "RIN Score")
  output$samples_hist_pmi <- make_hist("pmi", "PMI (hours)")
  
  # Bar chart: sample counts per group (categorical overview)
  output$samples_bar_plot <- renderPlot({
    group_col <- input$samples_group_col
    ggplot(metadata(), aes(x    = .data[[group_col]],
                           fill = .data[[group_col]])) +
      geom_bar(show.legend = FALSE) +
      labs(title = paste("Sample counts by", group_col),
           x     = group_col,
           y     = "Number of samples") +
      theme_bw()
  })
  
  # ── COUNTS ────────────────────────────────
  
  # Load and validate the normalized counts CSV
  # Expected format: first column = gene names, remaining = one column per sample
  norm_counts <- reactive({
    req(input$counts_file)
    validate_file(input$counts_file)
    
    mat <- read_csv(input$counts_file$datapath, show_col_types = FALSE)
    
    # Drop any unnamed index column (auto-named "...1" by readr when the
    # CSV was saved with row numbers but no column header for that column)
    mat <- mat %>% select(-any_of(c("...1")))
    
    # Now the first column should be GeneID — set it as rownames
    # and confirm the remaining columns are all numeric sample counts
    gene_col <- names(mat)[1]   # should be "GeneID"
    mat <- column_to_rownames(mat, var = gene_col)
    
    # Validate that what remains is numeric (catches future format issues)
    validate(
      need(all(sapply(mat, is.numeric)),
           paste("Non-numeric columns remain after removing gene ID column:",
                 paste(names(mat)[!sapply(mat, is.numeric)], collapse = ", ")))
    )
    
    mat
  })
  
  # Step 2: Filter genes by variance percentile.
  # This reactive depends on norm_counts() and the slider input, so it
  # re-runs automatically whenever either changes.
  filtered_counts <- reactive({
    counts    <- norm_counts()
    gene_vars <- apply(counts, 1, var)
    threshold <- quantile(gene_vars,
                          probs = input$var_percentile / 100,
                          na.rm = TRUE)
    counts[gene_vars >= threshold, ]
  })
  
  
  # ── COUNTS: per-gene statistics ───────────────────────────────────────────
  # Compute per-gene summary stats (median, variance, zero count) across ALL
  # genes once, then label each gene as passing or failing the current filter.
  # This single reactive feeds both diagnostic scatter plots and the summary table,
  # avoiding redundant computation when the slider changes.
  gene_stats <- reactive({
    req(norm_counts())
    counts    <- norm_counts()
    gene_vars <- apply(counts, 1, var)                   # variance per gene
    medians   <- apply(counts, 1, median)                # median count per gene
    n_zeros   <- apply(counts, 1, function(x) sum(x == 0))  # zeros per gene
    threshold <- quantile(gene_vars,
                          probs = input$var_percentile / 100,
                          na.rm = TRUE)
    tibble(
      gene     = rownames(counts),
      median   = medians,
      variance = gene_vars,
      n_zeros  = n_zeros,
      # Label each gene so plots can color passing vs. filtered-out genes
      passes   = gene_vars >= threshold
    )
  })
  
  # Summary table: number/% of genes passing and failing the variance filter
  output$counts_summary_table <- renderDT({
    req(gene_stats())
    gs          <- gene_stats()
    n_total     <- nrow(gs)
    n_pass      <- sum(gs$passes)
    n_fail      <- n_total - n_pass
    
    tibble(
      Metric = c(
        "Number of samples",
        "Total number of genes",
        "Genes passing filter",
        "% passing filter",
        "Genes NOT passing filter",
        "% not passing filter",
        "Variance percentile cutoff"
      ),
      Value = c(
        ncol(norm_counts()),
        n_total,
        n_pass,
        paste0(round(100 * n_pass / n_total, 1), "%"),
        n_fail,
        paste0(round(100 * n_fail / n_total, 1), "%"),
        paste0(input$var_percentile, "%")
      )
    ) %>%
      datatable(options = list(dom = "t", ordering = FALSE), rownames = FALSE)
  })
  
  # Diagnostic scatter 1: median count vs. variance (both on log10 scale).
  # Genes passing the filter are dark blue; filtered-out genes are light grey.
  output$counts_scatter_var <- renderPlot({
    req(gene_stats())
    ggplot(gene_stats(),
           aes(x     = log10(median + 1),
               y     = log10(variance + 1),
               color = passes)) +
      geom_point(alpha = 0.5, size = 0.8) +
      # Dark color for passing genes, light for filtered-out genes
      scale_color_manual(
        values = c("TRUE" = "#1a5276", "FALSE" = "#aed6f1"),
        labels = c("TRUE" = "Passes filter", "FALSE" = "Filtered out"),
        name   = ""
      ) +
      labs(
        title = "Median count vs. variance (log10 scale)",
        x     = "log10(median count + 1)",
        y     = "log10(variance + 1)"
      ) +
      theme_bw() +
      theme(legend.position = "top")
  })
  
  # Diagnostic scatter 2: median count vs. number of zero-count samples.
  # Same color scheme as scatter 1 for consistency.
  output$counts_scatter_zeros <- renderPlot({
    req(gene_stats())
    ggplot(gene_stats(),
           aes(x     = log10(median + 1),
               y     = n_zeros,
               color = passes)) +
      geom_point(alpha = 0.5, size = 0.8) +
      scale_color_manual(
        values = c("TRUE" = "#1a5276", "FALSE" = "#aed6f1"),
        labels = c("TRUE" = "Passes filter", "FALSE" = "Filtered out"),
        name   = ""
      ) +
      labs(
        title = "Median count vs. number of zero-count samples",
        x     = "log10(median count + 1)",
        y     = "Number of samples with zero counts"
      ) +
      theme_bw() +
      theme(legend.position = "top")
  })
  
  # Clustered heatmap of the top N highest-variance genes after filtering.
  # Uses pheatmap which draws row/col dendrograms and a color bar legend.
  output$counts_heatmap <- renderPlot({
    req(filtered_counts())
    counts <- filtered_counts()
    
    # Select top N genes by variance for the heatmap (subset if needed)
    gene_vars  <- apply(counts, 1, var)
    top_genes  <- names(sort(gene_vars, decreasing = TRUE))[
      seq_len(min(input$top_genes_heatmap, nrow(counts)))
    ]
    mat <- as.matrix(counts[top_genes, ])
    
    # Optionally log2-transform for better color scaling on lowly expressed genes
    if (input$heatmap_log) {
      mat <- log2(mat + 1)   # +1 avoids log(0)
      legend_title <- "log2(count + 1)"
    } else {
      legend_title <- "Normalized count"
    }
    
    # pheatmap handles clustering, color bar legend, and dendrograms automatically.
    # scale = "row" z-scores each gene so color reflects relative expression.
    pheatmap::pheatmap(
      mat,
      scale            = "row",          # z-score per gene row
      cluster_rows     = TRUE,           # hierarchical clustering of genes
      cluster_cols     = TRUE,           # hierarchical clustering of samples
      show_rownames    = FALSE,          # hide gene names (too many to display)
      show_colnames    = TRUE,
      color            = colorRampPalette(c("#2c7bb6", "white", "#d7191c"))(100),
      main             = paste0("Top ", input$top_genes_heatmap,
                                " variance genes (", legend_title, ", row z-scored)"),
      fontsize_col     = 8,
      legend           = TRUE            # show color bar legend
    )
  })
  
  # PCA scatter plot — user selects which two PCs to show via sidebar dropdowns.
  # Colors by condition if metadata is loaded; falls back to sample name coloring.
  output$counts_pca_plot <- renderPlot({
    req(filtered_counts())
    counts <- filtered_counts()
    
    validate(
      need(nrow(counts) >= 2, "Not enough genes to run PCA after filtering."),
      need(ncol(counts) >= 2, "Not enough samples to run PCA.")
    )
    
    # prcomp expects samples as rows → transpose the genes × samples matrix
    pca <- tryCatch(
      prcomp(t(counts), scale. = TRUE),
      error = function(e) validate(need(FALSE, paste("PCA error:", e$message)))
    )
    
    # Extract % variance explained for every PC (used in axis labels)
    pct_var <- round(summary(pca)$importance[2, ] * 100, 1)
    
    # Build a data frame of PC scores with sample IDs as a column
    pca_df <- as.data.frame(pca$x) %>%
      rownames_to_column("sample_id")
    
    # Determine the color variable: use condition from metadata if available,
    # otherwise fall back to coloring each sample individually by name
    if (!is.null(input$samples_file)) {
      meta <- tryCatch(metadata(), error = function(e) NULL)
      if (!is.null(meta)) {
        overlap <- intersect(pca_df$sample_id, meta[[1]])
        validate(
          need(length(overlap) > 0,
               "No matching sample IDs between counts columns and metadata first column.")
        )
        pca_df <- left_join(pca_df, meta,
                            by = c("sample_id" = "geo_accession"))
        
        color_col <- "condition"
      } else {
        color_col <- "sample_id"
      }
    } else {
      color_col <- "sample_id"
    }
    
    # Retrieve the user-selected PC names (e.g. "PC1", "PC3")
    x_pc <- input$pca_x
    y_pc <- input$pca_y
    
    # Look up % variance for the chosen PCs to annotate axes
    x_pct <- pct_var[x_pc]
    y_pct <- pct_var[y_pc]
    
    ggplot(pca_df, aes(x     = .data[[x_pc]],
                       y     = .data[[y_pc]],
                       color = .data[[color_col]],
                       label = sample_id)) +
      geom_point(size = 3) +
      labs(
        title = paste0("PCA — ", x_pc, " vs. ", y_pc,
                       " (variance-filtered genes)"),
        x     = paste0(x_pc, " (", x_pct, "% variance)"),
        y     = paste0(y_pc, " (", y_pct, "% variance)"),
        color = color_col
      ) +
      theme_bw()
  })
  
  # ── DE ────────────────────────────────────
  
  # Load and validate DESeq2 results CSV
  de_results <- reactive({
    req(input$de_file)
    validate_file(input$de_file)
    df <- read_csv(input$de_file$datapath, show_col_types = FALSE)
    # Confirm required columns are present
    required_cols <- c("baseMean", "log2FoldChange", "pvalue", "padj")
    validate(
      need(all(required_cols %in% names(df)),
           paste("Missing required columns:",
                 paste(setdiff(required_cols, names(df)), collapse = ", ")))
    )
    df
  })
  
  # Classify each gene as Up, Down, or NS based on threshold sliders
  de_classified <- reactive({
    de_results() %>%
      mutate(
        significance = case_when(
          padj < input$padj_thresh & log2FoldChange >  input$lfc_thresh ~ "Up",
          padj < input$padj_thresh & log2FoldChange < -input$lfc_thresh ~ "Down",
          TRUE ~ "NS"
        )
      )
  })
  
  # Summary table: how many genes in each category
  output$de_summary_table <- renderDT({
    de_classified() %>%
      count(significance, name = "n_genes") %>%
      datatable(options = list(dom = "t"), rownames = FALSE)
  })
  
  # Full DE results table with row highlighting
  output$de_results_table <- renderDT({
    datatable(
      de_classified(),
      rownames  = FALSE,
      filter    = "top",
      options   = list(pageLength = 15, scrollX = TRUE)
    ) %>%
      formatRound(columns = c("baseMean", "log2FoldChange",
                              "lfcSE", "stat"), digits = 3) %>%
      formatSignif(columns = c("pvalue", "padj"), digits = 3)
  })
  
  # Volcano plot
  output$de_volcano_plot <- renderPlot({
    df   <- de_classified()
    x_col <- input$volcano_x
    y_col <- input$volcano_y
    
    ggplot(df, aes(x = .data[[x_col]],
                   y = -log10(.data[[y_col]]),
                   color = significance)) +
      geom_point(alpha = 0.5, size = 1.2) +
      scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
      geom_hline(yintercept = -log10(input$padj_thresh),
                 linetype = "dashed", color = "black") +
      geom_vline(xintercept = c(-input$lfc_thresh, input$lfc_thresh),
                 linetype = "dashed", color = "black") +
      labs(
        title = "Volcano Plot — HD vs. Control",
        x     = x_col,
        y     = paste0("-log10(", y_col, ")")
      ) +
      theme_bw()
  })
  
  # MA plot: fold change vs. mean expression
  output$de_ma_plot <- renderPlot({
    df <- de_classified()
    ggplot(df, aes(x = log10(baseMean + 1),
                   y = log2FoldChange,
                   color = significance)) +
      geom_point(alpha = 0.4, size = 1) +
      scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(
        title = "MA Plot — HD vs. Control",
        x     = "log10(baseMean + 1)",
        y     = "log2 Fold Change"
      ) +
      theme_bw()
  })
  
  # ── GSEA ──────────────────────────────────────────────────────────────────
  
  # Load and validate GSEA results CSV
  gsea_results <- reactive({
    req(input$gsea_file)
    validate_file(input$gsea_file)
    df <- read_csv(input$gsea_file$datapath, show_col_types = FALSE)
    required_cols <- c("pathway", "NES", "padj")
    validate(
      need(all(required_cols %in% names(df)),
           paste("Missing required columns:",
                 paste(setdiff(required_cols, names(df)), collapse = ", ")))
    )
    df
  })
  
  # ── Tab 1: Bar plot data — top N pathways by adjusted p-value ─────────────
  gsea_top <- reactive({
    gsea_results() %>%
      arrange(padj) %>%
      slice_head(n = input$top_n_pathways) %>%
      mutate(
        pathway_label = str_trunc(pathway, 50),
        direction     = ifelse(NES > 0, "Activated", "Suppressed")
      )
  })
  
  output$gsea_bar_plot <- renderPlot({
    df <- gsea_top()
    ggplot(df, aes(x    = reorder(pathway_label, NES),
                   y    = NES,
                   fill = direction)) +
      geom_col() +
      scale_fill_manual(values = c("Activated" = "#c0392b",
                                   "Suppressed" = "#2980b9")) +
      coord_flip() +
      labs(
        title = paste("Top", input$top_n_pathways,
                      "GSEA Pathways by Adjusted P-value"),
        x     = "Pathway",
        y     = "Normalized Enrichment Score (NES)",
        fill  = "Direction"
      ) +
      theme_bw() +
      theme(axis.text.y = element_text(size = 9))
  })
  
  # Click on bar → show that pathway's table row below the plot
  output$gsea_clicked_row <- renderDT({
    req(input$gsea_bar_click)
    
    # nearPoints works on the flipped coordinates — match on y (the NES axis)
    # We identify the clicked pathway by finding the closest NES value
    df     <- gsea_top()
    click  <- input$gsea_bar_click
    
    # After coord_flip, x in plot space = NES (numeric), y = pathway (factor level)
    # Use the y click coordinate to identify the pathway rank position
    clicked_rank <- round(click$y)
    
    # Reorder to match plot order (reorder by NES), then index by rank
    df_ordered <- df %>% arrange(NES)
    
    if (clicked_rank >= 1 && clicked_rank <= nrow(df_ordered)) {
      selected <- df_ordered[clicked_rank, ] %>%
        select(-pathway_label, -direction) %>%
        mutate(across(where(is.list), ~ sapply(., paste, collapse = "; ")))
      
      datatable(selected, rownames = FALSE,
                options = list(dom = "t", scrollX = TRUE)) %>%
        formatRound(columns  = intersect(c("ES", "NES", "log2err"), names(selected)),
                    digits   = 3) %>%
        formatSignif(columns = intersect(c("pval", "padj"), names(selected)),
                     digits  = 3)
    }
  })
  
  # ── Tab 2: Filtered table ─────────────────────────────────────────────────
  gsea_table_filtered <- reactive({
    df <- gsea_results() %>%
      filter(padj <= input$gsea_padj_thresh)
    
    switch(input$nes_direction,
           "Positive (activated)"  = filter(df, NES > 0),
           "Negative (suppressed)" = filter(df, NES < 0),
           df  # "All"
    )
  })
  
  output$gsea_results_table <- renderDT({
    df <- gsea_table_filtered() %>%
      mutate(across(where(is.list), ~ sapply(., paste, collapse = "; ")))
    
    datatable(
      df,
      rownames  = FALSE,
      filter    = "top",
      options   = list(pageLength = 15, scrollX = TRUE)
    ) %>%
      formatRound(columns  = intersect(c("ES", "NES", "log2err"), names(df)),
                  digits   = 3) %>%
      formatSignif(columns = intersect(c("pval", "padj"), names(df)),
                   digits  = 3)
  })
  
  output$gsea_download <- downloadHandler(
    filename = function() {
      paste0("gsea_filtered_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- gsea_table_filtered() %>%
        mutate(across(where(is.list), ~ sapply(., paste, collapse = "; ")))
      write_csv(df, file)
    }
  )
  
  # ── Tab 3: NES vs -log10(padj) scatter ───────────────────────────────────
  output$gsea_scatter_plot <- renderPlot({
    df <- gsea_results() %>%
      mutate(
        sig       = padj <= input$gsea_scatter_padj,
        direction = case_when(
          sig & NES > 0 ~ "Activated",
          sig & NES < 0 ~ "Suppressed",
          TRUE          ~ "Not significant"
        ),
        label = ifelse(sig, str_trunc(pathway, 40), NA_character_)
      )
    
    ggplot(df, aes(x     = NES,
                   y     = -log10(padj),
                   color = direction,
                   label = label)) +
      geom_point(alpha = 0.7, size = 2) +
      geom_hline(yintercept = -log10(input$gsea_scatter_padj),
                 linetype = "dashed", color = "black") +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
      scale_color_manual(
        values = c("Activated"       = "#c0392b",
                   "Suppressed"      = "#2980b9",
                   "Not significant" = "grey70")
      ) +
      labs(
        title = "NES vs. Significance",
        x     = "Normalized Enrichment Score (NES)",
        y     = "-log10(adjusted p-value)",
        color = ""
      ) +
      theme_bw() +
      theme(legend.position = "top")
  })
}
 
# Run the application
shinyApp(ui = ui, server = server)