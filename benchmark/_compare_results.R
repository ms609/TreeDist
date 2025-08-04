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
  pr_results <- readRDS(pr_file)
  pr_replicate <- if (file.exists(replicate_file)) readRDS(replicate_file) else
    pr_results
  main_results <- readRDS(main_file)
  
  
  # Get the data frames of timings
  pr_df <- as.data.frame(pr_results)
  cf_df <- as.data.frame(pr_replicate)
  main_df <- as.data.frame(main_results)
  
  # Prepare a report
  report <- list()
  
  # Iterate over each function benchmarked
  for (fn_name in unique(pr_df$expr)) {
    pr_times <- pr_df[pr_df$expr == fn_name, "time"]
    cf_times <- cf_df[cf_df$expr == fn_name, "time"]
    main_times <- main_df[main_df$expr == fn_name, "time"]
    
    better_result <- t.test(pr_times, main_times, alternative = "less")
    worse_result <- t.test(pr_times, main_times, alternative = "greater")
    
    # The p-value tells us if the PR's performance is significantly slower
    # A small p-value (e.g., < 0.05) suggests it is.
    median_pr <- median(pr_times)
    median_cf <- median(cf_times)
    median_main <- median(main_times)
    percentage_change <- ((median_main - median_pr) / median_main) * 100
    
    delta <- abs(median_pr - median_main)
    df_delta <- abs(median_pr - median_cf)
    not_noise <- delta > (df_delta * 2)
    
    is_faster <- better_result$p.value < 0.01 && not_noise
    is_slower <- worse_result$p.value < 0.01 && not_noise
    
    report[[fn_name]] <- list(
      slower = is_slower,
      faster = is_faster,
      p_value = worse_result$p.value,
      median_pr = median_pr,
      median_cf = median_cf,
      median_main = median_main,
      change = percentage_change
    )
  }
  
  # Create a markdown-formatted message
  has_significant_regression <- FALSE
  
  for (fn_name in names(report)) {
    res <- report[[fn_name]]
    status <- ifelse(res$slower, "\U1F7E0 Slower \U1F641",
                     ifelse(res$faster, "\U1F7E2 Faster!",
                            ifelse(res$p_value < 0.05,
                                   "\U1F7E1 ?Slower",
                                   "\U26AA NSD")
                     )
    )
    if (res$slower) {
      has_significant_regression <- TRUE
    }
    
    bold <- ifelse(res$faster | res$slower, "**", "")
    
    message <- paste0(
      "| `", fn_name, "` | ", status, " | ", 
      bold, round(res$change, 2), "%", bold, "<br />(p: ", 
      format.pval(res$p_value), ") | ",
      signif(res$median_main / 1e6, 3), " \u2192<br />",
      signif(res$median_pr / 1e6, 3), ",  ",
      signif(res$median_cf / 1e6, 3), " |\n"
    )
  }
  
  if (has_significant_regression) {
    regressions <- TRUE
  }
  
  cat(message)
  output <- paste0(output, message)
}

cat(paste0(output, "\nEOF"), file = Sys.getenv("GITHUB_OUTPUT"), append = TRUE)

message(readLines(Sys.getenv("GITHUB_OUTPUT")))

# Fail the build if there is a significant regression
if (any(regressions)) {
  stop("Significant performance regression detected.")
}
