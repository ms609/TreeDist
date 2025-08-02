pr_files <- list.files("pr-benchmark-results", pattern = "*.bench.Rds",
                       full.names = TRUE)

regressions <- vapply(pr_files, function(pr_file) {
  file_name <- basename(pr_file)
  main_file <- file.path("main-benchmark-results", file_name)
  if (!file.exists(main_file)) return(NA);
  
  # Load the results
  pr_results <- readRDS(pr_file)
  main_results <- readRDS(main_file)
  
  
  # A simple comparison function using t-test
  # You can make this much more sophisticated.
  compare_timings <- function(pr_results, main_results) {
    # Get the data frames of timings
    pr_df <- as.data.frame(pr_results)
    main_df <- as.data.frame(main_results)
    
    # Prepare a report
    report <- list()
    
    # Iterate over each function benchmarked
    for (fn_name in unique(pr_df$expr)) {
      pr_times <- pr_df[pr_df$expr == fn_name, "time"]
      main_times <- main_df[main_df$expr == fn_name, "time"]
      
      better_result <- t.test(pr_times, main_times, alternative = "less")
      worse_result <- t.test(pr_times, main_times, alternative = "greater")
      
      # The p-value tells us if the PR's performance is significantly slower
      # A small p-value (e.g., < 0.05) suggests it is.
      is_faster <- better_result$p.value < 0.01
      is_slower <- worse_result$p.value < 0.01
      mean_pr <- mean(pr_times)
      mean_main <- mean(main_times)
      percentage_change <- ((mean_pr - mean_main) / mean_main) * 100
      
      report[[fn_name]] <- list(
        slower = is_slower,
        faster = is_faster,
        p_value = worse_result$p.value,
        mean_pr = mean_pr,
        mean_main = mean_main,
        change = percentage_change
      )
    }
    
    return(report)
  }
  
  # Generate the report
  report <- compare_timings(pr_results, main_results)
  
  # Create a markdown-formatted message
  message <- "### Performance Benchmark Results\n\n"
  has_significant_regression <- FALSE
  
  for (fn_name in names(report)) {
    res <- report[[fn_name]]
    status <- ifelse(res$slower, "\U1F7E0 Slower", 
                     ifelse(res$faster, "\U1F7E2 Faster!",
                            ifelse(res$p_value < 0.05,
                                   "\U1F7E1 A little slower (0.01 < p < 0.05)",
                                   "\U26AA No significant change")
                     )
    )
    if (res$slower) {
      has_significant_regression <- TRUE
    }
    
    message <- paste0(
      message,
      "#### `", fn_name, "`\n",
      "- Status: ", status, "\n",
      "- Mean time (PR): ", round(res$mean_pr / 1e6, 2), " ms\n",
      "- Mean time (Main): ", round(res$mean_main / 1e6, 2), " ms\n",
      "- Change: ", round(res$change, 2), "%\n",
      "- p-value: ", format.pval(res$p_value), "\n\n"
    )
  }
  if (has_significant_regression) {
    message <- paste0(message, "**Performance regression detected!**\n\n\n\n")
  }
  cat(message)
  has_significant_regression
}, FALSE)

# Fail the build if there is a significant regression
if (any(regressions)) {
  stop("Significant performance regression detected.")
} else {
  cat(message)
}

install.packages("gh")
gh::gh("POST /repos/{owner}/{repo}/issues/{issue_number}/comments",
  owner = "ms609",
  repo = "TreeDist",
  issue_number = Sys.getenv("PR_NUMBER"),
  body = message
)
