# Run from the repository root: Rscript tests/test_criteria.R
source("fsm_without_strata.R")
source("fsm_batch_without_strata.R")

new_units = data.frame(Index = 1:4, x = 1:4)
selection_order = data.frame(Treat = c(1, 2, 1, 2))
past_units = data.frame(Index = 5:8, x = 5:8, Treat = c(1, 2, 1, 2))
expected_allocation = data.frame(Index = 1:4, x = 1:4, Treat = c(1, 2, 1, 2))
expected_criteria = matrix(rep(-(1:4), each = 4), nrow = 4)
expected_criteria[lower.tri(expected_criteria)] = NA_real_

fits = list(
  fsm(new_units, selection_order, s_function = "min average", standardize = FALSE,
      ties = "smallest", units_print = FALSE),
  fsm_batch(new_units, past_units, t_ind = "Treat", SOM = selection_order,
            s_function = "min average", standardize = FALSE,
            ties = "smallest", units_print = FALSE)
)
for(fit in fits)
{
  stopifnot(identical(fit$criteria, expected_criteria),
            isTRUE(all.equal(fit$data_frame_allocated, expected_allocation)))
}
cat("fsm and fsm_batch: genuine -1 criteria, unused cells and allocations passed.\n")
