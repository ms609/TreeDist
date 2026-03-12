# Redirect all graphics to a null device for the test session.
# This prevents bare plot() calls and vdiffr snapshot rendering from
# appearing in the interactive graphics device (e.g. the RStudio viewer).
# vdiffr opens its own svglite device on top of this one, so snapshot
# tests are unaffected.
pdf(NULL)
teardown(dev.off())
