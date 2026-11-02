
#' Fisher's randomization test for sharp null hypothesis. 
#'
#' @description 
#' Performs Fisher's randomization test for sharp null hypotheses of the form 
#' \eqn{H_0: c_1 Y_i(1) + c_2 Y_i(2) - \tau = 0}, for a vector of contrasts \eqn{(c_1, c_2)}. 
#' @param Y_obs Vector of observed outcome.
#' @param alloc_obs Vector of observed treatment assignment.
#' @param alloc A matrix of treatment assignments over which the randomization distribution of the test statistic
#' is computed. Each row of \code{alloc} should correspond to an assignment vector. 
#' @param treat_lab Label of the treatment group in which the TASMD is computed. 
#' @param contrast A vector of the coefficients of the treatment contrast of interest. For example, for estimating the
#' average treatment effect of treatment 1 versus treatment 2, \code{contrast = c(1,-1)}.
#' @param tau The value of the treatment contrast specified by the sharp null hypothesis.
#' @param method The method of computing the test statistic. If \code{method = 'marginal mean'}, the test statistic
#' is \eqn{c_1 \hat{Y}_i(1) + c_2 \hat{Y}_i(2)}, where \eqn{\hat{Y}(z)} 
#' is the mean of the observed outcome in the group \eqn{Z = z}, for \eqn{z = 0,1}. If \code{method = 'marginal rank'}, 
#' the test statistic is \eqn{c_1 \hat{Y}_i(1) + c_2 \hat{Y}_i(2)}, 
#' where \eqn{\hat{Y}(z)} is the mean of the rank of the observed outcome in the group \eqn{Z = z}, for \eqn{z = 0,1}   
#' @param alternative The type of alternative hypothesis used. For right-sided test, \code{alternative = 'greater'}.
#' For left-sided test, \code{alternative = 'less'}. For both-sided test, \code{alternative = 'not equal'}.
#' @return A list containing the following items.
#' 
#' \code{test_stat_obs}: The observed value of the test statistic.
#' 
#' \code{test_stat_iter}: A vector of values of the test statistic across repeated randomizations.
#' 
#' \code{p_value}: p-value of the test. 
#' @export
#' @author Ambarish Chattopadhyay, Carl N. Morris and Jose R. Zubizarreta.
#' @references 
#' Chattopadhyay, A., Morris, C. N., and Zubizarreta, J. R. (2020), ``Randomized and Balanced Allocation 
#' of Units into Treatment Groups Using the Finite Selection Model for \code{R}".
#' @examples
#' # Consider N = 12, n1 = n2 = 6. 
#' # We test the sharp null of no treatment effect under CRD.
#' df_sample = data.frame(index = 1:12, x = c(20,30,40,40,50,60,20,30,40,40,50,60))
#' # True potential outcomes.
#' Y_1_true = 100 + (df_sample$x - mean(df_sample$x)) + rnorm(12, 0, 4)
#' Y_2_true = Y_1_true + 50
#' # Generate the realized assignment under CRD.
#' fc = crd(data_frame = df_sample, n_treat = 2, treat_sizes = c(6,6), control = FALSE)
#' Z_crd_obs = fc$Treat
#' # Get the observed outcomes
#' Y_obs = Y_1_true
#' Y_obs[Z_crd_obs == 2] = Y_2_true[Z_crd_obs == 2]
#' # Generate 5000 assignments under CRD.
#' Z_crd_iter = matrix(rep(0, 5000 * 12), nrow = 5000)
#' for(i in 1:5000)
#' {
#' fc = crd(data_frame = df_sample, n_treat = 2, treat_sizes = c(6,6), control = FALSE)
#' Z_crd_iter[i,] = fc$Treat
#' }
#' # Test for the sharp null H0: Y_i(1) = Y_i(0) for all i.
#' # Alternative: not H0 (two-sided test).
#' perm = perm_test(Y_obs = Y_obs, alloc_obs = Z_crd_obs, alloc = Z_crd_iter, 
#' contrast = c(1,-1), tau = 0, method = "marginal mean", alternative = 'not equal')
#' # Obtain the p-value.
#' perm$p_value




## Test for the sharp null hypothesis with two treatment groups
# H0: c1 * Y1 + c2 * Y2 - tau = 0

perm_test = function(Y_obs, alloc_obs, alloc, contrast = c(1,-1), tau = 0, 
                     method = "marginal mean", alternative = 'not equal')
{
  
  N = length(alloc_obs)
  level = sort(unique(alloc_obs))
  n.iter = nrow(alloc)
  
  Z_obs = alloc_obs
  
  S_perm = rep(0, n.iter)
  error_perm = rep(0, n.iter)
  error_obs = 0
  # Rank accumulation, division, weighting, summation and null subtraction.
  # gamma bounds their rounding error; inflate it for the computed scale.
  u = .Machine$double.eps / 2
  gamma = (N + 4) * u / (1 - (N + 4) * u)
  rank_error = gamma / (1 - gamma)
  
  for(iter in 1:n.iter)
  {
    Y1 = Y_obs
    Y2 = Y_obs
    
    Z_iter = alloc[iter,]
    # Fill in the missing potential outcomes under sharp null
    for(i in 1:N)
    {
      if(Z_iter[i] == level[1])
      {
        Y2[i] = (tau - contrast[1]*Y_obs[i])/contrast[2]
      }
      if(Z_iter[i] == level[2])
      {
        Y1[i] = (tau - contrast[2]*Y_obs[i])/contrast[1]
      }
    }
    
    # Y1 and Y2 now give the full set of potential outcomes   
    
    # Get the vector of realized outcomes under this assignment
    Y_real = Y1
    Y_real[Z_iter == level[2]] = Y2[Z_iter == level[2]]
    
    if(method == "marginal mean")
    {
      
      # Compute the mean response in each treatment group
      T_mean = c(mean(Y_real[Z_iter == level[1]]), mean(Y_real[Z_iter == level[2]]))
      S_perm[iter] = sum(T_mean * contrast) - tau
      
    }
    
    if(method == "marginal rank")
    {
      # Ranking the outcomes in ascending order
      Rj = rank(Y_real, ties.method = "average")
      
      # Compute the means of ranks in each group
      T_rank = c(mean(Rj[Z_iter == level[1]]), mean(Rj[Z_iter == level[2]]))
      S_perm[iter] = sum(T_rank * contrast) - tau
      error_perm[iter] = rank_error * (sum(abs(T_rank * contrast)) + abs(tau))
      
    }
    
  }
  
  ## Compute the observed test statistic
  if(method == "marginal mean")
  {
    # Compute the mean response in each treatment group
    T_mean = c(mean(Y_obs[Z_obs == level[1]]), mean(Y_obs[Z_obs == level[2]]))
    S_obs = sum(T_mean * contrast) - tau
  }
  
  if(method == "marginal rank")
  {
    Rj = rank(Y_obs, ties.method = "average")
    # Compute the means of ranks in each group
    T_rank = c(mean(Rj[Z_obs == level[1]]), mean(Rj[Z_obs == level[2]]))
    S_obs = sum(T_rank * contrast) - tau
    error_obs = rank_error * (sum(abs(T_rank * contrast)) + abs(tau))
  }
  
  
  ## Compute the p value
  # Count rank ties within the combined rounding bounds. Marginal-mean
  # comparisons remain exact, so outcome units do not set a tie threshold.
  tol = error_perm + error_obs
  if(alternative == "greater")
  {
    p_value = length(S_perm[S_perm >= S_obs - tol])/length(S_perm)  
  }
  if(alternative == "less")
  {
    p_value = length(S_perm[S_perm <= S_obs + tol])/length(S_perm)  
  }
  if(alternative == "not equal")
  {
    p_value = length(S_perm[abs(S_perm) >= abs(S_obs) - tol])/length(S_perm)  
  }
  
  return(list(test_stat_obs = S_obs, test_stat_iter = S_perm, p_value = p_value))
  
}
