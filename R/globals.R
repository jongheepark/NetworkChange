# Global variables declaration for R CMD check
# These are column names used in ggplot2 aes() calls with non-standard evaluation

utils::globalVariables(c(
  # plotV.R
  "Time", "Value", "Dimension",
  # plotU.R
  "First", "Second", "Name", "Names",
  # drawpostanalysis.R
  "Cluster",
  # BreakDiagnostic.R
  "model", "value", "Metric",
  # Utils.R
  "Size"
))
