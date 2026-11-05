# Run from the repository root: Rscript tests/test_ess_rand.R
source("design_assessment.R")

set.seed(1)
n_units = 24
covariate = rnorm(n_units)
outcomes = cbind(covariate + rnorm(n_units), covariate + rnorm(n_units) + 1)
make_assignments = function(n_iterations)
{
  assignments = array(0, c(n_iterations, n_units, 2))
  assignments[, , 1] = t(replicate(n_iterations, sample(rep(1:2, n_units / 2))))
  assignments[, , 2] = t(replicate(n_iterations, {
    allocation = integer(n_units)
    ordered = order(covariate)
    for(pair in seq(1, n_units, 2))
    {
      allocation[ordered[pair:(pair + 1)]] = sample(1:2)
    }
    allocation
  }))
  assignments
}

# Independent contrasts over the rows of each design, including R < N and R > N.
for(n_iterations in c(200, 10))
{
  assignments = make_assignments(n_iterations)
  for(contrast in list(c(1, -1), c(2, -3)))
  {
    variances = vapply(1:2, function(design) {
      statistics = apply(assignments[, , design], 1, function(allocation) {
        contrast[1] * mean(outcomes[allocation == 1, 1]) +
          contrast[2] * mean(outcomes[allocation == 2, 2])
      })
      var(statistics)
    }, numeric(1))
    expected = n_units * min(variances) / variances
    actual = ess_rand(assignments, outcomes, contrast)
    stopifnot(isTRUE(all.equal(actual, expected, tolerance = 1e-12)))
    reversed = ess_rand(assignments[n_iterations:1, , ], outcomes, contrast)
    stopifnot(isTRUE(all.equal(reversed, expected, tolerance = 1e-12)))
  }
}
cat("ess_rand: all iterations, both dimension orders and contrasts passed.\n")
