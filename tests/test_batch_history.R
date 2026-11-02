# Run from the repository root: Rscript tests/test_batch_history.R
source("fsm_without_strata.R")
source("fsm_batch_without_strata.R")

criteria = c("Dopt average", "Dopt pc", "marginal var sum")
requested = commandArgs(trailingOnly = TRUE)
if(length(requested)) criteria = match.arg(requested, criteria)
set.seed(412)
new = matrix(rnorm(16), 8, 2)
past = matrix(rnorm(12, 5, 4), 6, 2)
colnames(new) = colnames(past) = c("x", "y")
past_group = rep(1:2, 3)
order = rep(1:2, 4)
checks = 0
for(n_covariates in c(1, 2))
{
  for(standardize in c(FALSE, TRUE))
  {
    values = rbind(past, new)[, seq_len(n_covariates), drop = FALSE]
    if(standardize) values = scale(values)
    averages = rowMeans(values)
    # Independent eigendecomposition gives the shared principal direction.
    # Its arbitrary sign and the centering constant cancel from dispersion.
    direction = eigen(cov(values), symmetric = TRUE)$vectors[, 1]
    pc = drop(values %*% direction)
    for(criterion in criteria)
    {
      fit = fsm_batch(data.frame(Index = 1:8, new[, seq_len(n_covariates), drop = FALSE]),
        data.frame(Index = 9:14, past[, seq_len(n_covariates), drop = FALSE], Treat = past_group),
        "Treat", data.frame(Treat = order), s_function = criterion,
        standardize = standardize, ties = "smallest", units_print = FALSE)
      selected = fit$som_appended$Index
      stopifnot(setequal(selected, 1:8), !anyDuplicated(selected),
                identical(as.integer(table(fit$data_frame_allocated$Treat)), c(4L, 4L)))
      for(step in 1:8)
      {
        previous = seq_len(step - 1)
        group = c(which(past_group == order[step]),
                  nrow(past) + selected[previous[order[previous] == order[step]]])
        available = setdiff(1:8, selected[previous])
        expected = vapply(available, function(candidate) {
          rows = c(group, nrow(past) + candidate)
          if(criterion == "marginal var sum" && n_covariates > 1)
            return(sum(apply(values[rows, , drop = FALSE], 2, var)))
          scores = switch(criterion,
            "Dopt average" = averages[rows],
            "Dopt pc" = pc[rows],
            "marginal var sum" = values[rows, 1])
          sum((scores - mean(scores))^2)
        }, numeric(1))
        expected = round(expected, 7)
        stopifnot(isTRUE(all.equal(unname(fit$criteria[step, available]), expected, tolerance = 1e-8)),
                  selected[step] == available[which.max(expected)])
      }
      checks = checks + 1
    }
  }
}

# A transparent example: prior values -20,-10 require selecting new value 10.
if("Dopt average" %in% criteria)
{
  fit = fsm_batch(data.frame(Index = 1:4, x = c(0, 1, 8, 10)),
    data.frame(Index = 5:8, x = c(-20, -10, 10, 20), Treat = c(1, 1, 2, 2)),
    "Treat", data.frame(Treat = c(1, 2, 1, 2)), s_function = "Dopt average",
    standardize = FALSE, ties = "smallest", units_print = FALSE)
  stopifnot(fit$som_appended$Index[1] == 4,
            isTRUE(all.equal(unname(fit$criteria[1, ]), round(c(200, 662 / 3, 1208 / 3, 1400 / 3), 7))))
}
cat("Batch history:", checks, "allocations and", checks * 8, "independent dispersion choices passed.\n")
