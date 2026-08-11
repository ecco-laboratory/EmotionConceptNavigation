library(BayesFactor)

wide <- read.csv(
  "./outputs/6fold_by_2consistency_subavg_OFC2016ConstantinescuR5.csv",
  check.names = FALSE
)

head(wide)
names(wide)

# Keep participants with both conditions
wide_complete <- wide[complete.cases(wide[, c("0", "1")]), ]

# Paired Bayesian t-test:
# "1" = consistent
# "0" = inconsistent
bf_result <- ttestBF(
  x = wide_complete[["1"]],
  y = wide_complete[["0"]],
  paired = TRUE
)

bf_result
