# Run from the repository root: Rscript tests/test_scalar_policy_weight.R
source("fsm_without_strata.R")
source("fsm_batch_without_strata.R")

set.seed(19)
new = matrix(rnorm(24), 12, 2)
past = matrix(rnorm(12), 6, 2)
past_group = rep(1:2, 3)
order = rep(1:2, 6)
checks = 0
for(batch in c(FALSE, TRUE))
{
  for(weight in c(.5, 1, 2))
  {
    all = if(batch) rbind(past, new) else new
    design = cbind(1, all)
    n_past = if(batch) nrow(past) else 0
    current = design[n_past + 1:12, , drop = FALSE]
    policy = matrix(c(0, 1, -2), nrow = 1)
    # One policy's weighted covariance is w * p p', independent of diag().
    target = weight * crossprod(policy)
    args = list(data_frame = data.frame(Index = 1:12, new), SOM = data.frame(Treat = order),
                s_function = "Aopt", Pol_mat = policy, w_pol = weight,
                standardize = FALSE, ties = "smallest", units_print = FALSE)
    if(batch)
    {
      args$data_frame_past = data.frame(Index = 13:18, past, Treat = past_group)
      args$t_ind = "Treat"
      fit = do.call(fsm_batch, args)
    }
    else fit = do.call(fsm, args)
    selected = fit$som_appended$Index
    for(step in 1:12)
    {
      previous = seq_len(step - 1)
      group = current[selected[previous[order[previous] == order[step]]], , drop = FALSE]
      if(batch) group = rbind(design[which(past_group == order[step]), , drop = FALSE], group)
      information = crossprod(group)
      if(!nrow(group)) information = .001 * crossprod(design) / nrow(design)
      else if(qr(group)$rank < ncol(group))
        information = information / nrow(group) + .001 * crossprod(design) / nrow(design)
      available = setdiff(1:12, selected[previous])
      # Direct covariance reduction, not the implementation's rank-one formula.
      baseline = sum(diag(target %*% solve(information)))
      reduction = vapply(available, function(i) {
        baseline - sum(diag(target %*% solve(information + tcrossprod(current[i, ]))))
      }, numeric(1))
      stopifnot(isTRUE(all.equal(unname(fit$criteria[step, available]), round(reduction, 7),
                                tolerance = 1e-7)),
                max(reduction) - reduction[match(selected[step], available)] < 1e-7)
    }
    stopifnot(setequal(selected, 1:12), !anyDuplicated(selected))
    checks = checks + 1
  }
}
cat("Scalar policy weights:", checks, "allocations and", 12 * checks, "direct weighted-covariance choices passed.\n")
