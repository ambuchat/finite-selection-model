# Run from the repository root: Rscript tests/test_perm_test_ties.R
source("randomization inference.R")

alternatives = c("greater", "less", "not equal")
extreme = function(values, observed, alternative)
{
  switch(alternative,
         greater = values >= observed,
         less = values <= observed,
         "not equal" = abs(values) >= abs(observed))
}

# Twice the ranks and their group sums are integers. Cross multiplication
# compares rank contrasts exactly without reproducing floating-point means.
check_ranks = function(outcomes, assignments, observed)
{
  ranks = 2 * rank(outcomes, ties.method = "average")
  numerator = function(assignment)
  {
    n1 = sum(assignment == 1)
    n2 = sum(assignment == 2)
    n2 * sum(ranks[assignment == 1]) - n1 * sum(ranks[assignment == 2])
  }
  denominator = function(assignment)
  {
    2 * sum(assignment == 1) * sum(assignment == 2)
  }
  numerators = apply(assignments, 1, numerator)
  denominators = apply(assignments, 1, denominator)
  observed_numerator = numerator(observed)
  observed_denominator = denominator(observed)
  for(alternative in alternatives)
  {
    expected = mean(extreme(numerators * observed_denominator,
                           observed_numerator * denominators, alternative))
    for(scale in c(2^-40, 1, 2^40))
    {
      fit = perm_test(outcomes, observed, assignments,
                      contrast = scale * c(1, -1),
                      method = "marginal rank", alternative = alternative)
      stopifnot(fit$p_value == expected)
    }
  }
}

set.seed(2)
outcomes = round(rnorm(24, 10, 3), 2)
assignments = t(replicate(500, sample(rep(1:2, 12))))
observed = assignments[1, ]
check_ranks(outcomes, assignments, observed)
fit = perm_test(outcomes, observed, assignments, method = "marginal rank")
stopifnot(fit$p_value == 0.836)

# Exhaust every allocation, with repeated outcomes and a zero contrast.
assignments = t(combn(8, 4, FUN = function(treated)
{
  assignment = rep(2, 8)
  assignment[treated] = 1
  assignment
}))
outcomes = c(0, 0, 1, 2, 2, 3, 4, 4)
for(row in seq_len(nrow(assignments)))
{
  check_ranks(outcomes, assignments, assignments[row, ])
}
check_ranks(rep(1, 8), assignments, assignments[1, ])
check_ranks(-outcomes, assignments, assignments[1, ])
check_ranks(outcomes * 1e-10, assignments, assignments[1, ])
check_ranks(outcomes * 1e10, assignments, assignments[1, ])

# Unequal groups exercise different denominators in the rank means.
unequal = t(combn(9, 3, FUN = function(treated)
{
  assignment = rep(2, 9)
  assignment[treated] = 1
  assignment
}))
for(row in c(1, 42, 84))
{
  check_ranks(c(1, 1, 2, 3, 3, 4, 5, 6, 9), unequal, unequal[row, ])
}

# Distinct mean statistics must not become ties when outcome units shrink.
# The integer sums provide an independent oracle at power-of-two scales.
outcomes = 1:8
observed = assignments[1, ]
integers = apply(assignments, 1, function(assignment)
{
  sum(outcomes[assignment == 1]) - sum(outcomes[assignment == 2])
})
for(scale in c(2^-40, 2^-20, 1, 2^20, 2^40))
{
  for(alternative in alternatives)
  {
    expected = mean(extreme(integers, integers[1], alternative))
    fit = perm_test(outcomes * scale, observed, assignments,
                    alternative = alternative)
    stopifnot(fit$p_value == expected)
  }
}
fit = perm_test(outcomes * 1e-10, observed, assignments)
stopifnot(fit$p_value == 2 / 70)
cat("perm_test ties: integer rank oracle, all alternatives, cancellation and scaling passed.\n")
