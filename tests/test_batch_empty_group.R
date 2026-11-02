# Run from the repository root: Rscript tests/test_batch_empty_group.R
source("fsm_without_strata.R")
source("fsm_batch_without_strata.R")

set.seed(191)
new = data.frame(Index = 1:12, x = rnorm(12), y = rnorm(12))
past = data.frame(Index = 13:18, x = rnorm(6), y = rnorm(6), Treat = rep(1:2, each = 3))
order = rep(c(3, 1, 2), 4)
checks = 0
for(empty_history in c(FALSE, TRUE))
{
  prior_data = if(empty_history) past[FALSE, ] else past
  n_past = nrow(prior_data)
  for(standardize in c(FALSE, TRUE))
  {
    for(intercept in c(FALSE, TRUE))
    {
      values = as.matrix(rbind(prior_data[c("x", "y")], new[c("x", "y")]))
      if(standardize) values = scale(values)
      design = if(intercept) cbind(1, values) else values
      current = design[n_past + 1:12, , drop = FALSE]
      target = diag(seq_len(ncol(design)))
      for(criterion in c("Dopt", "Aopt"))
      {
        fit = fsm_batch(new, prior_data, "Treat", data.frame(Treat = order),
          s_function = criterion, standardize = standardize, intercept = intercept,
          Pol_mat = diag(ncol(design)), w_pol = seq_len(ncol(design)),
          ties = "smallest", units_print = FALSE)
        selected = fit$som_appended$Index
        stopifnot(setequal(selected, 1:12), !anyDuplicated(selected),
                  identical(as.integer(table(fit$data_frame_allocated$Treat)), rep(4L, 3)))
        for(step in 1:12)
        {
          previous = seq_len(step - 1)
          rows = c(which(prior_data$Treat == order[step]),
                   n_past + selected[previous[order[previous] == order[step]]])
          group = design[rows, , drop = FALSE]
          information = crossprod(group)
          if(!nrow(group)) information = .001 * crossprod(design) / nrow(design)
          else if(qr(group)$rank < ncol(group))
            information = information / nrow(group) + .001 * crossprod(design) / nrow(design)
          available = setdiff(1:12, selected[previous])
          # Rebuild each candidate's information matrix and evaluate its objective.
          scores = vapply(available, function(i) {
            augmented = information + tcrossprod(current[i, ])
            if(criterion == "Dopt") as.numeric(determinant(augmented, logarithm = TRUE)$modulus)
            else -sum(diag(target %*% solve(augmented)))
          }, numeric(1))
          stopifnot(max(scores) - scores[match(selected[step], available)] < 1e-7)
        }
        if(empty_history)
        {
          ordinary = fsm(new, data.frame(Treat = order), s_function = criterion,
            standardize = standardize, intercept = intercept,
            Pol_mat = diag(ncol(design)), w_pol = seq_len(ncol(design)),
            ties = "smallest", units_print = FALSE)
          stopifnot(identical(fit$criteria, ordinary$criteria),
                    identical(fit$som_appended, ordinary$som_appended))
        }
        checks = checks + 1
      }
    }
  }
}
cat("Empty batch groups:", checks, "allocations and", 12 * checks, "direct matrix choices passed.\n")
