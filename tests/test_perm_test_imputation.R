# Run from the repository root: Rscript tests/test_perm_test_imputation.R
source("randomization inference.R")

assignments = t(apply(combn(8, 4), 2, function(treated) {
  allocation = rep(2, 8)
  allocation[treated] = 1
  allocation
}))
observed = assignments[12, ]
outcomes = c(1, 4, 2, 8, 3, 9, 6, 7)

# Construct one potential-outcome schedule from the observed assignment, then
# enumerate every balanced assignment independently of perm_test's loops.
for(contrast in list(c(1, -1), c(2, -3), c(1, 1)))
{
  for(null in c(0, 3))
  {
    potential_1 = ifelse(observed == 1, outcomes,
                        (null - contrast[2] * outcomes) / contrast[1])
    potential_2 = ifelse(observed == 2, outcomes,
                        (null - contrast[1] * outcomes) / contrast[2])
    stopifnot(isTRUE(all.equal(contrast[1] * potential_1 + contrast[2] * potential_2,
                              rep(null, 8), tolerance = 1e-14)))
    expected = apply(assignments, 1, function(allocation) {
      contrast[1] * mean(potential_1[allocation == 1]) +
        contrast[2] * mean(potential_2[allocation == 2]) - null
    })
    expected_observed = contrast[1] * mean(outcomes[observed == 1]) +
      contrast[2] * mean(outcomes[observed == 2]) - null
    for(alternative in c("greater", "less", "not equal"))
    {
      actual = perm_test(outcomes, observed, assignments, contrast, null,
                         alternative = alternative)
      stopifnot(isTRUE(all.equal(actual$test_stat_iter, expected, tolerance = 1e-14)),
                isTRUE(all.equal(actual$test_stat_obs, expected_observed, tolerance = 1e-14)))
      expected_p = switch(alternative,
        greater = mean(expected >= expected_observed),
        less = mean(expected <= expected_observed),
        "not equal" = mean(abs(expected) >= abs(expected_observed)))
      stopifnot(actual$p_value == expected_p)
    }
    if(identical(contrast, c(1, -1)) && null == 0)
    {
      unchanged = apply(assignments, 1, function(allocation) {
        mean(outcomes[allocation == 1]) - mean(outcomes[allocation == 2])
      })
      stopifnot(identical(actual$test_stat_iter, unchanged))
    }
  }
}
cat("perm_test imputation: all 70 assignments, three contrasts and alternatives passed.\n")
