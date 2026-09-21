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
  main_replicate_file <- file.path("main2-benchmark-results", file_name)
  rep_exists <- file.exists(replicate_file)
  main_rep_exists <- file.exists(main_replicate_file)
  pr1 <- readRDS(pr_file)
  pr2 <- if (rep_exists) readRDS(replicate_file) else pr1
  main <- readRDS(main_file)
  main2 <- if (main_rep_exists) readRDS(main_replicate_file) else main
  
  # Prepare a report
  report <- list()
  
  # Use deparse1 for reliable expression-to-string conversion;
  # as.character(unlist()) decomposes call objects into their components.
  expr_names <- vapply(pr1[["expression"]], deparse1, "")
  
  # Iterate over each function benchmarked
  for (fn_name in unique(expr_names)) {
    pr1_times <-  as.numeric(pr1[["time"]][[1]])
    pr2_times <-  as.numeric(pr2[["time"]][[1]])
    pr_times <- if (rep_exists) c(pr1_times, pr2_times) else pr1_times
    main1_times <- as.numeric(main[["time"]][[1]])
    matched <- if (length(main1_times)) {
      TRUE
    } else {
      main1_times <- main[, "time"]
      FALSE
    }
    main2_times <- if (main_rep_exists && matched) {
      as.numeric(main2[["time"]][[1]])
    } else {
      main1_times
    }
    main_times <- if (main_rep_exists) c(main1_times, main2_times) else main1_times

    median_pr <- median(pr_times)
    median_main <- median(main_times)
    percentage_change <- ((median_main - median_pr) / median_main) * 100

    # Two sources of spread, and the wider one governs. Iteration-to-iteration
    # jitter within a run is the tighter of the two; the gap between two runs of
    # the SAME tree carries everything that changes between steps of the job
    # (frequency, page placement, a noisy neighbour), and on this runner that
    # has reached 25% on a 7 ms call. Measuring both trees twice is what makes
    # the second estimate available for each of them.
    within_run <- max(mad(main_times), mad(pr_times))
    between_run <- max(abs(median(pr1_times) - median(pr2_times)),
                       abs(median(main1_times) - median(main2_times)))
    noise <- max(within_run, between_run)

    threshold_percent <- 6 #  Changes of ~5% are frequent
    # Sub-millisecond benchmarks are dominated by system jitter on CI runners;
    # require an absolute difference floor before flagging.
    min_meaningful_diff <- 2e-4 # 0.2 ms (times are in seconds)
    abs_diff <- abs(median_pr - median_main)

    is_faster <- matched &&
      abs_diff > min_meaningful_diff &&
      median_pr < median_main - 2 * noise

    is_slower <- matched &&
      abs_diff > min_meaningful_diff &&
      median_pr > median_main + 2 * noise
    
    report[[fn_name]] <- list(
      matched = matched,
      slower = is_slower,
      faster = is_faster,
      median_pr = median(pr1_times),
      median_cf = median(pr2_times),
      median_main = median_main,
      median_main1 = median(main1_times),
      median_main2 = median(main2_times),
      change = percentage_change
    )
  }
  
  # Create a markdown-formatted message
  has_significant_regression <- FALSE
  
  for (i in seq_along(report)) {
    fn_name <- names(report)[[i]]
    res <- report[[i]]
    if (is.null(res) || is.null(res$matched)) next
    
    status <- if (isTRUE(res$matched)) {
      if (res$slower) {
        if (abs(res$change) > threshold_percent) {
          has_significant_regression <- TRUE
          "\U1F7E0 Slower \U1F641"
        } else {
          "\U1F7E3 ~same"
        }
      } else if (res$faster) {
        if (abs(res$change) > threshold_percent) {
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
      signif(res$median_main1 * 1e3, 3), ",  ",
      signif(res$median_main2 * 1e3, 3), " \u2192<br />",
      signif(res$median_pr   * 1e3, 3), ",  ",
      signif(res$median_cf   * 1e3, 3), " |\n"
    )
    
    cat(message)
    output <- paste0(output, message)
  }
  
  if (has_significant_regression) {
    regressions <- TRUE
  }
}

cat(paste0(output, "\nEOF"), file = Sys.getenv("GITHUB_OUTPUT"), append = TRUE)

# Fail the build if there is a significant regression
if (any(regressions)) {
  stop("Significant performance regression detected.")
}
