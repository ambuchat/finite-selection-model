# Run from the repository root: Rscript tests/test_batch_dimensions.R
source("fsm_without_strata.R")
source("fsm_batch_without_strata.R")

set.seed(48)
new = matrix(rnorm(36), 12, 3)
colnames(new) = paste0("x", 1:3)
order = rep(1:3, 4)
checks = 0
for(n_past in c(1, 2, 4))
{
  past = matrix(rnorm(9 * n_past), 3 * n_past, 3)
  colnames(past) = colnames(new)
  past_group = rep(1:3, n_past)
  for(n_covariates in c(1, 3))
  {
    for(intercept in c(FALSE, TRUE))
    {
      for(standardize in c(FALSE, TRUE))
      {
        all_covariates = rbind(past, new)[, seq_len(n_covariates), drop = FALSE]
        if(standardize) all_covariates = scale(all_covariates)
        design = if(intercept) cbind(1, all_covariates) else all_covariates
        prior = design[seq_len(nrow(past)), , drop = FALSE]
        current = design[nrow(past) + seq_len(nrow(new)), , drop = FALSE]
        fit = fsm_batch(data.frame(Index = 1:12, new[, seq_len(n_covariates), drop = FALSE]),
          data.frame(Index = 12 + seq_len(nrow(past)), past[, seq_len(n_covariates), drop = FALSE],
                     Treat = past_group), "Treat", data.frame(Treat = order),
          intercept = intercept, standardize = standardize, ties = "smallest", units_print = FALSE)
        selected = fit$som_appended$Index
        stopifnot(setequal(selected, 1:12), !anyDuplicated(selected),
                  identical(as.integer(table(fit$data_frame_allocated$Treat)), rep(4L, 3)))
        for(step in 1:12)
        {
          previous = seq_len(step - 1)
          group = rbind(prior[past_group == order[step], , drop = FALSE],
                        current[selected[previous[order[previous] == order[step]]], , drop = FALSE])
          information = crossprod(group)
          if(qr(group)$rank < ncol(group))
            information = information / nrow(group) + .001 * crossprod(design) / nrow(design)
          available = setdiff(1:12, selected[previous])
          # Recompute determinants directly instead of the quadratic-form shortcut.
          score = vapply(available, function(i) {
            as.numeric(determinant(information + tcrossprod(current[i, ]), logarithm = TRUE)$modulus)
          }, numeric(1))
          stopifnot(max(score) - score[match(selected[step], available)] < 1e-8)
        }
        checks = checks + 1
      }
    }
  }
}
cat("Batch dimensions:", checks, "allocations and", 12 * checks, "direct determinant choices passed.\n")
