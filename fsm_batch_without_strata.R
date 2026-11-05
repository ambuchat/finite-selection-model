#' Batched FSM for sequential experiments
#'
#' @description 
#' Extension of the FSM to cases where units arrive sequentially in batches.
#' @param data_frame Data frame containing a column of unit indices (optional) and covariates (or transformations thereof).
#' @param data_frame_past A data frame of units already allocated to treatment groups. 
#' Data frame contains a column of unit indices (optional), columns of covariates (or transformations thereof), 
#' and a column for treatment indicator. 
#' @param t_ind column name containing the treatment indicator in \code{data_frame_past}.
#' @param SOM Selection Order Matrix.

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
#' @param index_col_past \code{TRUE} if column of unit indices is present in \code{data_frame_past}.
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
#' \code{data_frame_allocated_augmented}: data frame combining \code{data_frame_allocated} and \code{data_frame_past}. 
#' @author Ambarish Chattopadhyay, Carl N. Morris and Jose R. Zubizarreta
#' @references 
#' Chattopadhyay, A., Morris, C. N., and Zubizarreta, J. R. (2020), ``Randomized and Balanced Allocation of Units into Treatment Groups Using the Finite Selection Model for \code{R}'.
#' @examples
#' # Consider N=18, number of treatments = 2, n1 = n2 = 9, batch sizes = 6,6,6.
#' # Get data frame for the first batch.
#' df_sample_1 = data.frame(index = 1:6, age = c(20,30,40,40,50,60))
#' # Obtain SOM for all the 12 units.
#' som_gen = som(data_frame = NULL, n_treat = 2, treat_sizes = c(9,9), 
#' include_discard = FALSE, method = 'SCOMARS', marginal_treat = rep((9/18),18), control = F)
#' # Assign the first batch.
#' f1 = fsm(data_frame = df_sample_1, SOM = som_gen[1:6,], s_function = 'Dopt', 
#' eps = 0.0001, ties = 'random', intercept = TRUE, standardize = TRUE, units_print = TRUE)
#' f1_app = f1$data_frame_allocated
#' # Get data frame for the second batch.
#' df_sample_2 = data.frame(index = 7:12, age = c(20,30,40,40,50,60))
#' # Assign the second batch.
#' f2 = fsm_batch(data_frame = df_sample_2, SOM = som_gen[7:12,], s_function = 'Dopt', 
#' eps = 0.0001, ties = 'random', intercept = TRUE, standardize = TRUE, units_print = TRUE,
#' data_frame_past = f1_app, t_ind = 'Treat', index_col_past = TRUE)
#' f2_app = f2$data_frame_allocated_augmented
#' # Get data frame for the third batch.
#' df_sample_3 = data.frame(index = 13:18, age = c(20,30,40,40,50,60))
#' # Assign the third batch.
#' f3 = fsm_batch(data_frame = df_sample_3, SOM = som_gen[13:18,], s_function = 'Dopt', 
#' eps = 0.0001, ties = 'random', intercept = TRUE, standardize = TRUE, units_print = TRUE,
#' data_frame_past = f2_app, t_ind = 'Treat', index_col_past = TRUE)
#' f3_app = f3$data_frame_allocated_augmented




fsm_batch = function(data_frame, data_frame_past, t_ind, SOM, s_function = 'Dopt', 
                     Q_initial = NULL, eps = 0.001, ties = 'random', intercept = TRUE, 
                     index_col_past = TRUE, standardize = TRUE, units_print = TRUE, 
                     index_col = TRUE, Pol_mat = NULL, w_pol = NULL)
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
  
  # if a column contains the same values, remove it to avoid singularity
  
  
  
  
  k = ncol(X_cov) # no. of covariates'
  
  # total size of past units
  N_past = nrow(data_frame_past)
  # treatment indicator for past units
  Z_past = data_frame_past[,t_ind]
  
  
  if(index_col_past == TRUE)
  {
    X_cov_past = as.matrix(data_frame_past[,-c(1, which(colnames(data_frame_past) == t_ind))]) # matrix of all covariates
    
  }
  if(index_col == FALSE)
  {
    X_cov_past = as.matrix(data_frame_past[,-which(colnames(data_frame_past) == t_ind)]) # matrix of all covariates
  }
  
  
  
  
  if(standardize == TRUE)
  {
    # standardize the columns of X_cov
    for(j in 1:ncol(X_cov))
    {
      #X_cov[,j] = (X_cov[,j] - mean(X_cov[,j]))/sd(X_cov[,j])
      #X_cov_past[,j] = (X_cov_past[,j] - mean(X_cov_past[,j]))/sd(X_cov_past[,j])
      # can also standardize based on the augmented matrix
      loc = mean(c(X_cov_past[,j], X_cov[,j]))
      scale = sd(c(X_cov_past[,j], X_cov[,j]))
      X_cov[,j] = (X_cov[,j] - loc)/scale
      X_cov_past[,j] = (X_cov_past[,j] - loc)/scale
      
      
    }
  }
  
  if(intercept == TRUE)
  {
    X_N = as.matrix(cbind(rep(1,N), X_cov)) # n x (k+1) design matrix for the new units
    colnames(X_N) = c('intercept', sprintf('x%d', 1:k))
    
    X_N_past = as.matrix(cbind(rep(1,N_past), X_cov_past)) # n x (k+1) design matrix for the past units
    colnames(X_N_past) = c('intercept', sprintf('x%d', 1:k))
    
  }
  
  if(intercept == FALSE)
  {
    X_N = as.matrix(X_cov) # n x k design matrix for the new units
    colnames(X_N) = sprintf('x%d', 1:k)
    
    X_N_past = as.matrix(X_cov_past) # n x k design matrix for the past units
    colnames(X_N) = sprintf('x%d', 1:k)
    
  }
  
  ## set initital nonsingular matrix
  if(is.null(Q_initial) == TRUE)
  {
    X_N_combined = rbind(X_N_past, X_N)
    
    Q0 = t(X_N_combined)%*%X_N_combined #invertible(?) matrix
  } else{
    Q0 = Q_initial
  }
  
  
  # Treatments take turn in selecting the units
  Z = rep(-1,N) # treatment indicator initialized at all -1
  
  crit_print = matrix(NA_real_, nrow = N, ncol = N)
  
  for(i in 1:N)
  {
    t.index = som_order[i] # treatment group that will pick a unit at this stage
    
    units.current = unit.index[Z == t.index] # new units already in that treatment group
    
    X_N_group_past = X_N_past[Z_past == t.index, ] # design matrix for the past units in the choosing treatment group
    
    # augment the past design matrix with the new design matrix
    if(length(units.current) == 0)
    {
      X_n = X_N_group_past
    }
    if(length(units.current)>0)
    {
      X_n = rbind(X_N_group_past, X_N[units.current,]) 
      
    }
    
    
    # reciprocal 'Condition number'
    if((kappa((t(X_n)%*% X_n)/i))^{-1} < 1e-6)
    {
      #print(t(X_n)%*% X_n + eps * Q0)
      Sn = solve((t(X_n)%*% X_n)/nrow(X_n) + eps * (Q0/(N_past + N)))
    } else{
      #print(t(X_n)%*% X_n)
      #Sn = solve((t(X_n)%*% X_n)/i)
      Sn = solve((t(X_n)%*% X_n))
    }
    
    
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
  #stage = 1:N
  som_appended = cbind(as.vector(som_order), data_frame[units.selected,])
  #colnames(som_appended)[1] = 'Sel order'
  colnames(som_appended)[1] = 'Treat'
  rownames(som_appended) = 1:N
  som_appended = as.data.frame(som_appended)
  som_split  = split(som_appended,som_appended$Treat)
  
  # augmented data frame after allocation
  if(index_col_past == TRUE && index_col == TRUE)
  {
    data_frame_allocated_augmented = as.data.frame(rbind(data_frame_past, data_frame_allocated))
  } 
  if(index_col_past == TRUE && index_col == FALSE)
  {
    temp1 = as.data.frame(rbind(data_frame_past[,-1], data_frame_allocated))
    data_frame_allocated_augmented = data.frame(Index = 1:(N_past + N), temp1)
  }
  if(index_col_past == FALSE && index_col == TRUE)
  {
    temp1 = as.data.frame(rbind(data_frame_past, data_frame_allocated[,-1]))
    data_frame_allocated_augmented = data.frame(Index = 1:(N_past + N), temp1)
  }
  if(index_col_past == FALSE && index_col == FALSE)
  {
    temp1 = as.data.frame(rbind(data_frame_past, data_frame_allocated))
    data_frame_allocated_augmented = data.frame(Index = 1:(N_past + N), temp1)
  }
  
  
  return(list(data_frame_allocated = data_frame_allocated, som_appended = som_appended, 
              som_split = som_split,
              data_frame_allocated_augmented = data_frame_allocated_augmented,
              criteria = crit_print))
  
  
}




