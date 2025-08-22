pr_files <- list.files("pr-benchmark-results", pattern = "*.bench.Rds",
                       full.names = TRUE)

output <- paste0(
  "report<<EOF\n### Performance benchmark results\n\n",
  "| Call     | Status | Change | Time (ms) |\n",
  "|----------|--------|--------|-----------|\n"
)

regressions <- FALSE

for (pr_file in pr_files) {
  file_name <- basename(pr_file)
  replicate_file <- file.path("pr2-benchmark-results", file_name)
  main_file <- file.path("main-benchmark-results", file_name)
  if (!file.exists(main_file)) next;
  
  # Load the results
  rep_exists <- file.exists(replicate_file)
  pr1 <- readRDS(pr_file)
  pr2 <- if (rep_exists) readRDS(replicate_file) else pr1
  main <- readRDS(main_file)
  
  # Prepare a report
  report <- list()
  
  # Iterate over each function benchmarked
  for (fn_name in unique(as.character(unlist(pr1[["expression"]])))) {
    pr1_times <-  as.numeric(pr1[["time"]][[1]])
    pr2_times <-  as.numeric(pr2[["time"]][[1]])
    pr_times <- if (rep_exists) c(pr1_times, pr2_times) else pr1_times
    main_times <- as.numeric(main[["time"]][[1]])
    matched <- if (length(main_times)) {
      TRUE
    } else {
      main_times <- main[, "time"]
      FALSE
    }
    
    median_pr <- median(pr_times)
    median_main <- median(main_times)
    mad_main <- mad(main_times)
    percentage_change <- ((median_main - median_pr) / median_main) * 100
    
    q <- 0.2
    main_iqr <- quantile(main_times, c(q, 1 - q))
    interrun_var <- median(pr1_times) - median(pr2_times)
    noise_range <- range(
      median(main_times) + c(1, -1) * interrun_var,
      main_iqr)
    
    threshold_percent <- 4
    
    is_faster <- matched &&
      median_pr < median_main - 2 * mad_main &&
      median_pr < noise_range[[1]]
    
    is_slower <- matched &&
      median_pr > median_main + 2 * mad_main &&
      median_pr > noise_range[[2]]
    
    report[[fn_name]] <- list(
      matched = matched,
      slower = is_slower,
      faster = is_faster,
      median_pr = median(pr1_times),
      median_cf = median(pr2_times),
      median_main = median_main,
      change = percentage_change
    )
  }
  
  # Create a markdown-formatted message
  has_significant_regression <- FALSE
  
  for (fn_name in names(report)) {
    res <- report[[fn_name]]
    status <- if (res$matched) {
      if (res$slower) {
        if (abs(percentage_change) > threshold_percent) {
          has_significant_regression <- TRUE
          "\U1F7E0 Slower \U1F641"
        } else {
          "\U1F7E3 ~same"
        }
      } else if (res$faster) {
        if (abs(percentage_change) > threshold_percent) {
          "\U1F7E2 Faster!"
        } else {
          "\U1F7E3 ~same"
        }
      } else {
        "\U26AA NSD"
      }
    } else {
      "\U1F7E4 ?Mismatch"
    }
    
    bold <- ifelse(res$faster | res$slower, "**", "")
    
    message <- paste0(
      "| `", fn_name, "` | ", status, " | ", 
      bold, round(res$change, 2), "%", bold, " | ", 
      signif(res$median_main * 1e3, 3), " \u2192<br />",
      signif(res$median_pr   * 1e3, 3), ",  ",
      signif(res$median_cf   * 1e3, 3), " |\n"
    )
  }
  
  if (has_significant_regression) {
    regressions <- TRUE
  }
  
  cat(message)
  output <- paste0(output, message)
}

cat(paste0(output, "\nEOF"), file = Sys.getenv("GITHUB_OUTPUT"), append = TRUE)

# Fail the build if there is a significant regression
if (any(regressions)) {
  stop("Significant performance regression detected.")
}
