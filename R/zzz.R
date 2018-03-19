## Reminders when releasing for CRAN
release_questions <- function() {
  c(
    "Have you cleared GitHub issues for this release milestone?",
    "Have you checked the Vignettes for sanity?"
  )
}

#rhub::check_on_windows()
#check_with_rdevel() # redundifies check_on_debian()
#check_on_ubuntu()
#check_on_fedora()
#check_on_centos() 
#check_with_valgrind() # runs the build and check on Linux, in valgrind to find memory leaks and pointer errors.
#check_with_sanitizers() # runs all package package tests, examples and vignettes with Address Sanitizer and Undefined Behavior Sanitizer.
