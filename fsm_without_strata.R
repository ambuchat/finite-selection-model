## Basic FSM without stratification

library(MASS)
library(matrixcalc)

## Finding hcf/gcd of a vector of numbers
hcf = function(x)
{
  h = 1
  for(i in 1:min(x))
  {
    if(max(x %% i) == 0)
    {
      h = i
    }
    
  }
  return(h)
}


########################################################################################################

#########################################################################################################

### Function to construct a Selection Order Matrix



#' Selection Order Matrix (SOM)
#'
#' @description 
#' Generates a Selection Order Matrix (SOM) in a deterministic/random manner. 
#' @param data_frame A (optional) data frame corresponding to the full sample of units. 
#' Required if \code{include_discard = TRUE}.
#' @param n_treat Number of treatment groups.
#' @param treat_sizes A vector of treatment group sizes. If \code{control = TRUE}, 
#' the first element of \code{treat_sizes} should be the control group size.
#' @param include_discard \code{TRUE} if a discard group is considered.
#' @param method Specifies the selection strategy used among \code{'global percentage'}, 
#' \code{'randomized chunk'}, \code{'SCOMARS'}. \code{'SCOMARS'} is applicable only if \code{n_treat = 2}. 
#' @param control If \code{TRUE}, treatments are labeled as 0,1,...,g-1 (0 representing the control group). 
#' If \code{FALSE}, they are labeled as 1,2,...,g.
#' @param marginal_treat A vector of marginal probabilities, the jth element being the probability that treatment group
#'  (or treatment group 2 in case \code{control = FALSE}) gets to choose at the jth stage given
#' the total number of choices made by treatment group upto the (j-1)th stage. 
#' Only applicable when \code{method = 'SCOMARS'}. 
#' @return A data frame containing the selection order of treatments, i.e. the labels of treatment groups 
#' at each stage of selection. If \code{method = 'SCOMARS'}, the data frame contains an additional column of 
#' the conditional selection probabilities. 
#' @export
#' @author Ambarish Chattopadhyay, Carl N. Morris and Jose R. Zubizarreta.
#' @references 
#' Chattopadhyay, A., Morris, C. N., and Zubizarreta, J. R. (2020), 
#' ``Randomized and Balanced Allocation of Units into Treatment Groups Using the Finite Selection Model for \code{R}''.
#' 
#' Morris, C. (1983), ``Sequentially controlled Markovian random sampling (SCOMARS)'', Institute of 
#' Mathematical Statistics Bulletin,12(5), 237.
#' @examples
#' Generate an SOM with N = 12, n1 = n2 = 6.
#' som_sample = som(data_frame = NULL, n_treat = 2, treat_sizes = c(6,6), include_discard = FALSE, 
#' method = 'SCOMARS', control = FALSE, marginal_treat = rep(6/12,12))


som = function(data_frame = NULL, n_treat, treat_sizes, include_discard = FALSE, method = 'SCOMARS', control = FALSE, marginal_treat = NULL)
{
  if(include_discard == TRUE)
  {
    # treat the discard group as an additional treatment group
    g = n_treat + 1 
    N = nrow(data_frame)
    discard.size = N - sum(treat_sizes)
    n = c(treat_sizes, discard.size)
  }
  if(include_discard == FALSE)
  {
    # leave certain individuals unallocated, if any
    g = n_treat
    N = sum(treat_sizes)
    n = treat_sizes
  }
  
  sel.order = data.frame(Treat = rep(0,N)) # initialize SOM
  
  if(method == 'global percentage')
  {
    current.sizes = rep(0,g)
    for(i in 1:N)
    {
      crit = (current.sizes + 0.5)/n # evaluate the relative completeness of each group
      candidates = which(crit == min(crit))
      # resolve ties randomly
      candidate.opt = candidates[sample(x = 1:length(candidates),size = 1)]
      # 'candidate.opt' chooses a unit at this stage
      sel.order[i,] = candidate.opt
      # re-adjusting the group sizes
      current.sizes[candidate.opt] = current.sizes[candidate.opt] + 1 
    }
  }
  
  
  if(method == 'randomized chunk') # applies for proportional group sizes
  {
    n.hcf = hcf(n)
    chunk.prop = n/n.hcf
    chunk.size = sum(chunk.prop)
    
    # randomly permute 1st chunk.size selections, (chunk.size+1)th to (2*chunk.size)th selection and so on
    for(chunk in 1:n.hcf)
    {
      tr = rep(1:g,chunk.prop)
      sel.order[((chunk-1)*chunk.size+1): ((chunk-1)*chunk.size+chunk.size),] = sample(tr,chunk.size,replace = FALSE)
    }
    
  }
  
  if(method == 'SCOMARS') 
  {
    # set the marginal probabilities of receiving treatment
    p = marginal_treat
    # vector whose ith element gives the number of times treatment gets to choose upto ith stage
    S = rep(0,N)
    # FF = E(S)
    FF = cumsum(p)
    
    # initialize the probability that treatment gets ti choose at the ith stage
    prob.treat = p
    
    # initial step: treatment gets to choose at the 1st step with prob p1
    sel.order[1,] = rbinom(1,1,prob.treat[1])
    S[1] = as.vector(sel.order[1,])
    
    for(i in 1:(N-1))
    {
      prob.treat[i+1] = punif((p[i+1] - max(c(0,S[i] - FF[i])))/(1-abs(S[i] - FF[i])))
      sel.order[(i+1),] = rbinom(1,1,prob.treat[i+1])
      S[i+1] = S[i] + as.vector(sel.order[(i+1),])

    }
    
    sel.order[['Treat']] = sel.order[['Treat']] + 1
  }
  
  
  if(control == TRUE)
  {
    sel.order[['Treat']] = sel.order[['Treat']] - 1 
  }
  
  if(method == 'SCOMARS')
  {
    sel.order2 = data.frame(Probs = round(prob.treat,4), Treat = sel.order[['Treat']])
    return(sel.order2)
  }
  else
  {
    return(sel.order)
  }
  
}


#################################################################################################


### Function for obtaining an assignment under the FSM given an SOM


#' Finite Selection Model (FSM)
#'
#' @description 
#' Generates a randomized assignment of a group of units to multiple groups of pre-determined 
#' sizes using the Finite Selection Model (FSM).
#' @param data_frame A data frame containing a column of unit indices (optional) and covariates (or transformations thereof).
#' @param SOM A selection order matrix.
#' @param s_function Specifies a selection function, a string among \code{'constant'}, \code{'Dopt'}, 
#' \code{'Aopt'}, \code{'max pc'}, \code{'min pc'}, \code{'Dopt pc'}, \code{'max average'}, \code{'min average'},
#' \code{'Dopt average'}. \code{'constant'} selection function puts a constant value on every unselected unit. 
#' \code{'Dopt'} use the D-optimality criteria based on the full set of covariates to select units. 
#' \code{'Aopt'} uses the A-optimality criteria. \code{'max pc'} (respectively, \code{'min pc'}) selects that 
#' unit that has the maximum (respectively, minimum) value of the first principal component. 
#' \code{'Dopt pc'} uses the D-optimality criteria on the first principal component, \code{'max average'} 
#' (respectively, \code{'min average'}) selects that unit that has the maximum (respectively, minimum) 
#' value of the simple average of the covariates. \code{'Dopt average'} uses the D-optimality criteria on the 
#' simple average of the covariates.
#' @param Q_initial A (optional) non-singular matrix (called 'initial matrix') that is added the \eqn{(X^T X)} 
#' matrix of the choosing treatment group at any stage, when the \eqn{(X^T X)} matrix of that treatment group
#' at that stage is non-invertible. If \code{FALSE}, the \eqn{(X^T X)} matrix for the full set of observations is used
#' as the non-singular matrix. Applicable if \code{s_function = 'Dopt'} or \code{'Aopt'}.
#' @param eps Proportionality constant for \code{Q_initial}, the default value is 0.001.
#' @param ties Specifies how to deal with ties in the values of the selection function. If \code{ties = 'random'},
#'  a unit is selected randomly from the set of candidate units. If \code{ties = 'smallest'}, the unit 
#'  that appears earlier in the data frame, i.e. the unit with the smallest index gets selected.
#' @param intercept if \code{TRUE}, the design matrix of each treatment group includes a column of intercepts.
#' @param standardize if \code{TRUE}, the columns of the \eqn{X} matrix other than the column for the intercept (if any), 
#' are standardized.
#' @param units_print if \code{TRUE}, the function automatically prints the candidate units at each step of selection.
#' @param index_col if \code{TRUE}, data_frame contains a column of unit indices.
#' @param Pol_mat Policy matrix. Applicable only when \code{s_function = 'Aopt'}.
#' @param w_pol A vector of policy weights. Applicable only when \code{s_function = 'Aopt'}.
#' @export
#' @return A list containing the following items.
#'
#' \code{data_frame_allocated}:  The original data frame augmented with the column of the treatment indicator.
#' 
#' \code{som_appended}:  The SOM with augmented columns for the indices and covariate values for units selected.
#' 
#' \code{som_split}:  som_appended, split by the levels of the treatment.
#' 
#' \code{crit_print}:  The value of the objective function, at each stage of build up process. At each stage, 
#' the unit that maximizes the objective function is selected.
#' @author Ambarish Chattopadhyay, Carl N. Morris and Jose R. Zubizarreta
#' @references 
#' Chattopadhyay, A., Morris, C. N., and Zubizarreta, J. R. (2020), ``Randomized and Balanced Allocation of 
#' Units into Treatment Groups Using the Finite Selection Model for \code{R}''.
#' 
#' Morris, C. (1979), ``A finite selection model for experimental design of the health insurance study'',
#' Journal of Econometrics, 11(1), 43–61.
#' 
#' Morris, C., Hill, J. (2000), ``The health insurance experiment: design using the finite selection model'', 
#' Public policy and statistics: case studies from RAND, Springer Science & Business Media, 29–53.
#' @examples
#' # Load the data.
#' df_sample = data.frame(index = 1:12, x = c(20,30,40,40,50,60,20,30,40,40,50,60))
#' Generate an SOM with N = 12, n1 = n2 = 6.
#' som_sample = som(n_treat = 2, treat_sizes = c(6,6), method = 'SCOMARS', control = TRUE, 
#' marginal_treat = rep(6/12,12))
#' Assign units given the SOM.
#' f = fsm(data_frame = df_sample, SOM = som_sample, s_function = 'Dopt', 
#' eps = 0.001, ties = 'random', intercept = TRUE, standardize = TRUE, units_print = TRUE, 
#' index_col = TRUE)



fsm = function(data_frame, SOM, s_function = 'Dopt', 
               Q_initial = NULL, eps = 0.001, ties = 'random', intercept = TRUE, 
               standardize = TRUE, units_print = TRUE, index_col = TRUE, Pol_mat = NULL, w_pol = NULL)
{
  # names of all possible selection functions
  sf.names = c('constant', 'Dopt', 'Aopt', 'negative Dopt',
               'max pc', 'min pc', 'Dopt pc', 'max average', 'min average', 'Dopt average',
               'marginal var sum')
  
  
  if(ncol(SOM)>1)
  {
    som_order = SOM[['Treat']] # treatments should be labelled 1,2,...,g or 0,1,...,g-1
  }
  if(ncol(SOM)==1)
  {
    som_order = SOM[,1] # treatments should be labelled 1,2,...,g or 0,1,...,g-1
  }
  
  
  if(index_col == TRUE)
  {
    unit.identity = data_frame[['Index']]
  }
  
  unit.index = 1:nrow(data_frame)
  g = length(table(som_order)) # no. of treatments
  n = as.vector(table(som_order)) # vector of treatment group sizes
  N = sum(n) # total no. of units in the sample
  
  ## build-up phase
  units.selected = rep(0,N)
  if(index_col == TRUE)
  {
    X_cov = as.matrix(data_frame[,-1]) # matrix of all covariates
    
  }
  if(index_col == FALSE)
  {
    X_cov = as.matrix(data_frame) # matrix of all covariates
  }
  
  k = ncol(X_cov) # no. of covariates'
  
  if(standardize == TRUE)
  {
    # standardize the columns of X_cov
    for(j in 1:ncol(X_cov))
    {
      X_cov[,j] = (X_cov[,j] - mean(X_cov[,j]))/sd(X_cov[,j])
    }
  }
  
  if(intercept == TRUE)
  {
    X_N = as.matrix(cbind(rep(1,N), X_cov)) # n x (k+1) design matrix
    colnames(X_N) = c('intercept', sprintf('x%d', 1:k))
    
  }
  
  if(intercept == FALSE)
  {
    X_N = as.matrix(X_cov) # n x k design matrix
    colnames(X_N) = sprintf('x%d', 1:k)
    
  }
  
  ## set initital nonsingular matrix
  if(is.null(Q_initial) == TRUE)
  {
    Q0 = t(X_N)%*%X_N #invertible(?) matrix
  } else{
    Q0 = Q_initial
  }
  
  
  # Treatments take turn in selecting the units
  Z = rep(-1,N) # treatment indicator initialized at all -1
  
  crit_print = matrix(NA_real_, nrow = N, ncol = N)
  
  for(i in 1:N)
  {
    t.index = som_order[i] # treatment group that will pick a unit at this stage
    units.current = unit.index[Z == t.index] # units already in that treatment group
    # {
    if(length(units.current) == 0)
    {
      Sn = solve(eps * Q0/N) # starting matrix 
    }
    else
    {
      X_n = X_N[units.current,] 
      dim(X_n) = c(length(units.current),ncol(X_N)) 
      # treats X_n as a row vector when units.current has length 1
      
      # reciprocal 'Condition number'
      if((kappa((t(X_n)%*% X_n)/i))^{-1} < 1e-6)
      {
        #print(t(X_n)%*% X_n + eps * Q0)
        Sn = solve((t(X_n)%*% X_n)/nrow(X_n) + eps * (Q0/N))
      }
      else
      {
        #print(t(X_n)%*% X_n)
        #Sn = solve((t(X_n)%*% X_n)/i)
        Sn = solve((t(X_n)%*% X_n))
      }
    }
    #  }
    units.search = unit.index[Z == -1] # units yet to be allocated
    if(units_print == TRUE)
    {
      print(units.search)
    }
    # evaluate the criterion function on each of these units
    crit.func = rep(0, length(units.search))
    for(u in 1:length(units.search))
    {
      
      if(s_function == 'constant')
      {
        crit.func[u] = 10 #a constant
      }
      
      
      if(s_function == 'Dopt')
      {
        crit.func[u] = as.vector(t(X_N[units.search[u],]) %*% Sn %*% X_N[units.search[u],] )
      }
      
      
      if(s_function == 'negative Dopt')
      {
        crit.func[u] = (-1) * as.vector(t(X_N[units.search[u],]) %*% Sn %*% X_N[units.search[u],] )
      }
      
      if(s_function == 'Aopt')
      {
        # Policy matrix - Pol_mat, Matrix of weights = w_pol
        # I am assuming cost to be constant
        Pol_mat = as.matrix(Pol_mat)
        W = diag(w_pol)
        T.mat = t(Pol_mat) %*% W %*% Pol_mat
        
        crit.func[u] = as.vector(t(X_N[units.search[u],]) %*% (Sn %*% T.mat %*% Sn) %*% 
                                   X_N[units.search[u],])/(1 + as.vector(t(X_N[units.search[u],]) %*%  
                                                                           Sn %*% X_N[units.search[u],]) ) 
      }
      
      if(s_function %in% c('max pc', 'min pc', 'Dopt pc'))
      {              
        # Extract the 1st principal component from the full set of covariates
        pca = prcomp(X_cov)
        x.pc = pca$x[,1] #1st principal component
        
        if(s_function == 'max pc')
        {
          
          # include the unit which maximizes the 1st principal component
          crit.func[u] = x.pc[units.search[u]]
        }
        
        if(s_function == 'min pc')
        {
          # include the unit which minimizes the 1st principal component
          crit.func[u] = -x.pc[units.search[u]]
        }
        
        if(s_function == 'Dopt pc')
        {
          # include the unit which maximizes the dispersion of the 1st principal component
          # within the chooser group
          x.pc.append = x.pc[c(units.current,units.search[u])]
          crit.func[u] = sum((x.pc.append - mean(x.pc.append))^2)
        }
      }
      
      if(s_function == 'max average') # takes simple average of all covariates
      {
        # first check if there are more than one covariates
        if(k==1)
        {
          crit.func[u] =  as.vector(X_cov)[units.search[u]] 
        }else
        {
          crit.func[u] =  rowMeans(X_cov)[units.search[u]]
        }
        
      }
      
      if(s_function == 'min average')
      {
        # first check if there are more than one covariates
        if(k==1)
        {
          crit.func[u] =  -as.vector(X_cov)[units.search[u]] 
        }
        else
        {
          crit.func[u] =  -rowMeans(X_cov)[units.search[u]]
        }
      }
      
      if(s_function == 'Dopt average')
      {
        # first check if there are more than one covariate
        if(k==1)
        {
          x.avg.append = as.vector(X_cov)[c(units.current,units.search[u])]
          crit.func[u] = sum((x.avg.append - mean(x.avg.append))^2) 
        }
        else
        {
          x.avg = rowMeans(X_cov)
          x.avg.append = x.avg[c(units.current,units.search[u])]
          crit.func[u] = sum((x.avg.append - mean(x.avg.append))^2) 
        }
      }
      
      if(s_function == 'marginal var sum')
      {
        # first check if there are more than one covariate
        if(k==1)
        {
          x.append = as.vector(X_cov)[c(units.current,units.search[u])]
          crit.func[u] = sum((x.append - mean(x.append))^2) 
        }
        else
        {
          X_treat = X_cov[c(units.current,units.search[u]),]
          if(is.null(nrow(X_treat)) == TRUE)
          {
            crit.func[u] = 0 
          }
          if(is.null(nrow(X_treat)) == FALSE)
          {
            crit.func[u] = matrix.trace(cov(X_treat)) 
          }
          
        } 
      }
      
      if((s_function %in% sf.names) == FALSE)
      {
        stop('Invalid selection function')
      }
      
    }
    
    crit.func = round(crit.func,7) # rounding off the values (?)
    
    crit_print[i,units.search] = crit.func
    
    # all candidate units
    units.opt = units.search[which(crit.func == max(crit.func))]
    
    # resolve ties
    if(ties == 'random')
    {
      unit.opt = units.opt[sample(x = 1:length(units.opt),size = 1)]
      Z[unit.opt] = t.index
      
    }
    if(ties == 'smallest')
    {
      unit.opt = units.opt[1]
      Z[unit.opt] = t.index
      
    }
    # unit number that is selected at this stage
    units.selected[i] = unit.opt
  }
  
  
  data_frame_allocated = cbind(data_frame,Z)
  colnames(data_frame_allocated)[ncol(data_frame_allocated)] = 'Treat'
  som_appended = cbind(as.vector(som_order), data_frame[units.selected,])
  colnames(som_appended)[1] = 'Treat'
  rownames(som_appended) = 1:N
  som_appended = as.data.frame(som_appended)
  som_split  = split(som_appended, som_appended$Treat)
  
  return(list(data_frame_allocated = data_frame_allocated, som_appended = som_appended, 
              som_split = som_split, criteria = crit_print))
  
}


