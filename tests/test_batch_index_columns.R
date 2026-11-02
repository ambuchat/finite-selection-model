# Run from the repository root: Rscript tests/test_batch_index_columns.R
source("fsm_without_strata.R")
source("fsm_batch_without_strata.R")

set.seed(72)
new = data.frame(Index = 101:112, x = rnorm(12), y = rnorm(12))
past = data.frame(Index = 201:212, x = rnorm(12), y = rnorm(12), Treat = rep(1:3, 4))
order = data.frame(Treat = rep(1:3, 4))
checks = 0
for(standardize in c(FALSE, TRUE))
{
  expected = fsm_batch(new, past, "Treat", order, standardize = standardize,
                       ties = "smallest", units_print = FALSE)
  for(new_index in c(FALSE, TRUE))
  {
    for(past_index in c(FALSE, TRUE))
    {
      fit = fsm_batch(if(new_index) new else new[, -1],
                      if(past_index) past else past[, -1], "Treat", order,
                      index_col = new_index, index_col_past = past_index,
                      standardize = standardize, ties = "smallest", units_print = FALSE)
      stopifnot(identical(fit$data_frame_allocated$Treat, expected$data_frame_allocated$Treat),
                identical(fit$criteria, expected$criteria),
                identical(fit$som_appended[c("x", "y")], expected$som_appended[c("x", "y")]),
                identical(fit$data_frame_allocated_augmented[c("x", "y", "Treat")],
                          expected$data_frame_allocated_augmented[c("x", "y", "Treat")]))
      checks = checks + 1
    }
  }
}
cat("Batch index columns:", checks, "representations preserve allocation, criteria and combined data.\n")
