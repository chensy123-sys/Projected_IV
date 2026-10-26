# library(dplyr); i <- 1
# robust_mean <- function(x) {
#   x <- x[!is.na(x)]
#   if (length(x) == 0) return(NA)
#   mean(x, trim = 0.0)
# }
# for(n in c(2000,5000,10000)){
#   for (i in 1:7) {
#     path <- paste0("results/",n,"/RES",i,"_5.rds")
#     RES <- readRDS(path)
#     RES %>% group_by(method) %>%
#       summarise(
#         Bias = abs(robust_mean(est)),
#         SD = sd(est),
#         SE = robust_mean(sd),
#         CR95 = mean(0 >= (est - 1.96 * sd) & 0 <= (est + 1.96 * sd))
#       ) %>% ungroup() -> rec
#     print(i)
#     print(rec)
#   }
# }



for(n in c(2000,5000,10000)){
  args <- commandArgs(trailingOnly = TRUE)
  results_dir <- if (length(args) >= 1L) args[[1L]] else paste0("results/",n)
  output_tex <- if (length(args) >= 2L) args[[2L]] else paste0("results/tex/simulation_table_",n,".tex")
  
  dgp_ids <- 1:6
  method_order <- c(
    "PIV1", "PIV2", "PIV3",
    "GE1", "GE2", "GE3",
    "NUC", "TSHS", "sisVIVE"
  )
  
  method_headers <- c(
    "PIV$_1$", "PIV$_2$", "PIV$_3$",
    "GE$_1$", "GE$_2$", "GE$_3$",
    "NUC", "TSHS", "sisVIVE"
  )
  
  metric_order <- c("BIAS", "SD", "SE", "CR")
  
  safe_mean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else mean(x,trim = 0.0)
  }
  
  robust_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2L) return(NA_real_)
    iqr <- unname(quantile(x, 0.75) - quantile(x, 0.25))
    iqr / (qnorm(0.75) - qnorm(0.25))
  }
  
  coverage_rate <- function(est, se) {
    keep <- is.finite(est) & is.finite(se)
    if (!any(keep)) return(NA_real_)
    mean(
      0 >= est[keep] - 1.96 * se[keep] &
        0 <= est[keep] + 1.96 * se[keep]
    )
  }
  
  format_value <- function(x) {
    if (!is.finite(x)) return("--")
    value <- sprintf("%.4f", x)
    value <- sub("^0\\.", ".", value)
    sub("^-0\\.", "-.", value)
  }
  
  summarise_method <- function(data, method) {
    part <- data[data$method == method, , drop = FALSE]
    if (nrow(part) == 0L) {
      stop(sprintf("Missing method '%s'.", method), call. = FALSE)
    }
    
    c(
      BIAS = abs(safe_mean(part$est)),
      SD = sd(part$est),
      SE = safe_mean(part$sd),
      CR = coverage_rate(part$est, part$sd)
    )
  }
  
  result_files <- file.path(results_dir, sprintf("RES%d_5.rds", dgp_ids))
  missing_files <- result_files[!file.exists(result_files)]
  if (length(missing_files) > 0L) {
    stop(
      paste("Missing result file(s):", paste(missing_files, collapse = ", ")),
      call. = FALSE
    )
  }
  
  all_summaries <- vector("list", length(dgp_ids))
  replication_counts <- integer(length(dgp_ids))
  
  for (index in seq_along(dgp_ids)) {
    dgp <- dgp_ids[[index]]
    data <- readRDS(result_files[[index]])
    
    required_columns <- c("method", "est", "sd")
    missing_columns <- setdiff(required_columns, names(data))
    if (length(missing_columns) > 0L) {
      stop(
        sprintf(
          "%s is missing column(s): %s",
          result_files[[index]],
          paste(missing_columns, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    
    if ("loop" %in% names(data)) {
      replication_counts[[index]] <- length(unique(data$loop))
    } else {
      method_counts <- table(data$method)
      replication_counts[[index]] <- as.integer(method_counts[[method_order[[1L]]]])
    }
    
    summary_matrix <- vapply(
      method_order,
      function(method) summarise_method(data, method),
      numeric(length(metric_order))
    )
    rownames(summary_matrix) <- metric_order
    colnames(summary_matrix) <- method_order
    all_summaries[[index]] <- summary_matrix
  }
  
  if (length(unique(replication_counts)) != 1L) {
    warning(
      "DGP files have different replication counts: ",
      paste(replication_counts, collapse = ", ")
    )
  }
  replications <- replication_counts[[1L]]
  
  dir_name <- basename(normalizePath(results_dir, mustWork = TRUE))
  sample_size_text <- strsplit(dir_name, "_", fixed = TRUE)[[1L]][[1L]]
  sample_size <- suppressWarnings(as.integer(sample_size_text))
  if (is.na(sample_size)) {
    sample_size_caption <- sample_size_text
    sample_size_label <- gsub("[^[:alnum:]]+", "-", sample_size_text)
  } else {
    sample_size_caption <- format(sample_size, big.mark = ",", scientific = FALSE)
    sample_size_label <- as.character(sample_size)
  }
  
  lines <- c(
    "\\begin{table}[t]",
    "  \\centering",
    sprintf(
      "  \\caption{Comparison of simulation results for $n=%s$ over %d replications.}",
      sample_size_caption,
      replications
    ),
    sprintf("  \\label{tbl: simu %s}", sample_size_label),
    "  \\renewcommand{\\arraystretch}{1}",
    "  \\resizebox{0.8\\textwidth}{!}{%",
    "    \\begin{tabular}{cc|ccc|ccc|cccc}",
    "      \\toprule",
    paste0(
      "      DGP & Metrics & ",
      paste(method_headers, collapse = " & "),
      " \\\\"
    ),
    "      \\midrule"
  )
  
  for (index in seq_along(dgp_ids)) {
    dgp <- dgp_ids[[index]]
    summary_matrix <- all_summaries[[index]]
    lines <- c(lines, sprintf("      \\multirow{4}{*}{%d}", dgp))
    
    for (metric_index in seq_along(metric_order)) {
      metric <- metric_order[[metric_index]]
      values <- vapply(summary_matrix[metric, ], format_value, character(1L))
      row_end <- if (metric_index == length(metric_order) && index < length(dgp_ids)) {
        " \\\\[10pt]"
      } else {
        " \\\\"
      }
      lines <- c(
        lines,
        paste0(
          "      & ",
          sprintf("%-4s", metric),
          " & ",
          paste(values, collapse = " & "),
          row_end
        )
      )
    }
    
    if (index < length(dgp_ids)) lines <- c(lines, "")
  }
  
  lines <- c(
    lines,
    "      \\bottomrule",
    "    \\end{tabular}",
    "  }",
    "\\end{table}",
    ""
  )
  
  output_dir <- dirname(output_tex)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  writeLines(lines, output_tex, useBytes = TRUE)
  
  message(sprintf("Wrote %s", normalizePath(output_tex, mustWork = TRUE)))
  
}
