## Love plot for experimental design with multiple treatment groups


# Function for love plot of ASMD
love_plot_asmd = function(data_frame, alloc1, alloc2 = NULL, treat_lab = c(0,1), 
                          vline = "", index_col = TRUE, xupper = 1.0,
                          legend_text= "FSM",legend_position = 'topright') 
{
  
  Z_data = alloc1
  
  # get the matrix of covariates X_mat
  if(index_col == TRUE)
  {
    X_mat = data_frame[,-1]
  }
  
  if(index_col == FALSE)
  {
    X_mat = data_frame
  }
  
  # separating the matrix of covariates for the two treatment groups
  X_mat_1_1 = X_mat[Z_data == treat_lab[1],]
  X_mat_2_1 = X_mat[Z_data == treat_lab[2],]
  
  # column wise means and variances of covariates in the two treatment groups
  means_tr1_1 = apply(X_mat_1_1, 2, mean)
  var_tr1_1 = apply(X_mat_1_1, 2, var)
  means_tr2_1 = apply(X_mat_2_1, 2, mean)
  var_tr2_1 = apply(X_mat_2_1, 2, var)
  
  # calculate asmd
  asmd_1 = abs(means_tr2_1 - means_tr1_1)/sqrt((var_tr1_1 + var_tr2_1)/2)
  
  if(is.null(alloc2) == TRUE)
  {
    
    n_aux = length(asmd_1)
    
    dotchart(asmd_1[n_aux:1], labels=colnames(X_mat)[n_aux:1], 
             cex = 0.7, pch="", color= , main="", xlim=c(0,xupper),
             xlab = "ASMD")
    points(asmd_1[n_aux:1], y=1:ncol(X_mat), cex=1.2, 
           pch=19)
    
    legend(legend_position, legend_text, 
           cex=0.8, bty="n", pch=c(16,0), col=c("black","red"))
    abline(v = 0, lty = 2)
    abline(v=vline, lty=2)
    
  }
  
  if(is.null(alloc2) == FALSE)
  {
    Z_data = alloc2
    # separating the matrix of covariates for the two treatment groups
    X_mat_1_2 = X_mat[Z_data == treat_lab[1],]
    X_mat_2_2 = X_mat[Z_data == treat_lab[2],]
    
    # column wise means and variances of covariates in the two treatment groups
    means_tr1_2 = apply(X_mat_1_2, 2, mean)
    var_tr1_2 = apply(X_mat_1_2, 2, var)
    means_tr2_2 = apply(X_mat_2_2, 2, mean)
    var_tr2_2 = apply(X_mat_2_2, 2, var)
    
    # calculate asmd
    asmd_2 = abs(means_tr2_2 - means_tr1_2)/sqrt((var_tr1_2 + var_tr2_2)/2)
    
    
    n_aux = length(asmd_2)
    
    dotchart(asmd_1[n_aux:1], labels=colnames(X_mat)[n_aux:1], 
             cex = 0.7, pch="", color= , main="", xlim=c(0,xupper), 
             xlab = "ASMD")
    points(asmd_1[n_aux:1], y=1:ncol(X_mat), cex=1.2, 
           pch=19, col = "black")
    points(asmd_2[n_aux:1], y=1:ncol(X_mat), cex=1.2, 
           pch=4, col = "black")
    legend(legend_position, legend_text, 
           cex=0.8, bty="n", pch=c(19,4), col=c("black","black"))
    abline(v = 0, lty = 2)
    abline(v=vline, lty=2)
    
  }
  
}


################################################################

# Function for love plot of TASMD
love_plot_tasmd = function(data_frame, alloc1, alloc2 = NULL, treat_lab = 1, 
                           vline = "", index_col = TRUE, xupper = 1.0, 
                           mean_tar = NULL, sd_tar = NULL, denom = 'target',
                           legend_text= "FSM",legend_position = 'topright') 
{
  
  Z_data = alloc1
  
  # get the matrix of covariates X_mat
  if(index_col == TRUE)
  {
    X_mat = data_frame[,-1]
  }
  
  if(index_col == FALSE)
  {
    X_mat = data_frame
  }
  
  # separating the matrix of covariates for the given treatment group
  X_mat_t = X_mat[Z_data == treat_lab,]
  
  # column wise means of covariates in the two treatment groups
  means_t = apply(X_mat_t, 2, mean)
  
  if(is.null(mean_tar) == TRUE)
  {
    means_tar = apply(X_mat, 2, mean)
    if(denom == 'target')
    {
      sd_denom = apply(X_mat, 2, sd)
    }
    if(denom == 'group')
    {
      sd_denom = apply(X_mat_t, 2, sd)
    }
    
  }
  
  if(is.null(mean_tar) == FALSE)
  {
    means_tar = mean_tar
    if(denom == 'target')
    {
      sd_denom = sd_tar
    }
    if(denom == 'group')
    {
      sd_denom = apply(X_mat_t, 2, sd)
    }
    
  }
  
  # calculate tasmd
  tasmd = abs(means_t - means_tar)/sd_denom
  
  if(is.null(alloc2) == TRUE)
  {
    
    n_aux = length(tasmd)
    
    dotchart(tasmd[n_aux:1], labels=colnames(X_mat)[n_aux:1], 
             cex = 0.7, pch="", color= , main="", xlim=c(0,xupper),
             xlab = "TASMD")
    points(tasmd[n_aux:1], y=1:ncol(X_mat), cex=1.2, 
           pch=19,)
    
    legend(legend_position, legend_text, 
           cex=0.8, bty="n", pch=c(16,0), col=c("black","red"))
    abline(v = 0, lty = 2)
    abline(v=vline, lty=2)
    
  }
  
  if(is.null(alloc2) == FALSE)
  {
    Z_data = alloc2
    # separating the matrix of covariates for the two treatment groups
    X_mat_t2 = X_mat[Z_data == treat_lab,]
    
    # column wise means of covariates in the two treatment groups
    means_t2 = apply(X_mat_t2, 2, mean)
    
    if(denom == 'group')
    {
      sd_denom = apply(X_mat_t2, 2, sd)
    }
    
    # calculate asmd
    tasmd2 = abs(means_t2 - means_tar)/sd_denom
    
    
    n_aux = length(tasmd2)
    
    dotchart(tasmd[n_aux:1], labels = colnames(X_mat)[n_aux:1], 
             cex = 0.7, pch = "", color = , main = "", xlim = c(0,xupper), 
             xlab = "TASMD")
    points(tasmd[n_aux:1], y = 1:ncol(X_mat), cex = 1.2, 
           pch = 19, col = "black")
    points(tasmd2[n_aux:1], y = 1:ncol(X_mat), cex = 1.2, 
           pch = 4, col = "black")
    legend(legend_position, legend_text, 
           cex = 0.8, bty = "n", pch=c(19,4), col = c("black","black"))
    abline(v = 0, lty = 2)
    abline(v = vline, lty = 2)
    
  }
  
}



#' Love plot 
#'
#' @description 
#' Generates a Love plot of Absolute Standardized Mean Differences (ASMD) or Target Absolute Standardized Differences (TASMD) 
#' between two groups under one or two designs. 
#' @param data_frame Data frame containing a column of unit indices (optional) and covariates (or transformations thereof).
#' @param index_col if \code{TRUE}, \code{data_frame} contains a column of unit indices.
#' @param alloc1 A vector of treatment assignment.
#' @param alloc2 A (optional) vector of treatment assignment. 
#' @param imbalance Measure of imbalance used. If \code{imbalance = 'TASMD'}, imbalance is computed using 
#' the Target Absolute Standardized Mean Differences (TASMD). If \code{imbalance = 'ASMD'}, 
#' imbalance is computed using the Absolute Standardized Mean Differences (ASMD)
#' @param treat_lab Label of the treatment group in which the TASMD is computed. Applicable only when \code{imbalance = 'TASMD'}. 
#' @param vline A (optional) x-coordinate at which a vertical line is drawn.
#' @param xupper Upper limit of the x-axis.
#' @param mean_tar A (optional) vector of target profile of the covariates under consideration, 
#' e.g., mean of the covariates in the target population. Applicable only when \code{imbalance = 'TASMD'}. 
#' If \code{mean_tar = NULL}, the full-sample average of the covariates is considered as the target profile. 
#' @param sd_tar A optional vector of the standard deviation of the covariates in the target population. 
#' Applicable only when \code{imbalance = 'TASMD'}.
#' @param denom Specifies the denominator for the computation of TASMD. If \code{denom = 'target'}, 
#' the standard deviations of the covariates in the target population are used. If \code{denom = 'group'}, 
#' the standard deviations of the covariates in the treatment group given by \code{treat_lab} are used. 
#' Applicable only when \code{imbalance = 'TASMD'}. 
#' @param legend_text Legend of the two designs under consideration.
#' @param legend_position = Position of the legend in the plot. The default is \code{'topright'}.
#' @return Love plot of the ASMD/TASMD of the covariates.
#' @export
#' @author Ambarish Chattopadhyay, Carl N. Morris and Jose R. Zubizarreta.
#' @references 
#' Chattopadhyay, A., Morris, C. N., and Zubizarreta, J. R. (2020), ``Randomized and Balanced Allocation 
#' of Units into Treatment Groups Using the Finite Selection Model for \code{R}".
#' 
#' Love, T. (2004), “Graphical display of covariate balance”, Presentation, 
#' See http://chrp.org/love/JSM2004RoundTableHandout.pdf, 1364.
#' @examples
#' # Consider the Lalonde dataset.
#' # Get the full sample size.
#' N = nrow(Lalonde)
#' # Get the treatment group sizes.
#' n1 = floor(N/2)
#' n2 = N-n1
#' # Generate an SOM.
#' som_obs = som(n_treat = 2, treat_sizes = c(n1,n2),include_discard = FALSE,
#' method = 'SCOMARS', marginal_treat = rep((n2/N),N), control = FALSE)
#' # Generate a treatment assignment given som_obs.
#' f = fsm(data_frame = Lalonde, SOM = som_obs, s_function = 'Dopt', eps = 0.0001, 
#' ties = 'random', intercept = TRUE, standardize = TRUE, units_print = FALSE)
#' # Get assignment vector under the FSM.
#' Z_fsm_obs = f$data_frame_allocated$Treat
#' # Draw a random CRD.
#' Z_crd_obs = crd(data_frame = Lalonde, n_treat = 2, treat_sizes = c(n1, n2), 
#' control = FALSE)$Treat
#' # Draw Love plot.
#' love_plot(data_frame = Lalonde, index_col = TRUE, alloc1 = Z_fsm_obs, alloc2 = Z_crd_obs, 
#' imbalance = 'TASMD', treat_lab = 1, mean_tar = NULL, sd_tar = NULL, denom = 'target',
#' vline = "", legend_text = c("FSM","CRD"), xupper = 0.15, legend_position = 'topright') 



## General loveplot function

love_plot = function(data_frame, index_col = TRUE, alloc1, alloc2 = NULL, imbalance = 'TASMD',
                     treat_lab = 1, vline = "", xupper = 1.0,
                     mean_tar = NULL, sd_tar = NULL, denom = 'target',
                     legend_text= "FSM",legend_position = 'topright'){
  if(imbalance == 'ASMD')
  {
    love_plot_asmd(data_frame = data_frame, alloc1 = alloc1, alloc2 = alloc2, treat_lab = treat_lab, 
                   vline = vline, index_col = index_col, xupper = xupper,
                   legend_text = legend_text, legend_position = legend_position)
  }
  
  if(imbalance == 'TASMD')
  {
    love_plot_tasmd(data_frame = data_frame, alloc1 = alloc1, alloc2 = alloc2, treat_lab = treat_lab, 
                    vline = vline, index_col = index_col, xupper = xupper, 
                    mean_tar = mean_tar, sd_tar = sd_tar, denom = denom,
                    legend_text = legend_text, legend_position = legend_position)
  }
  
  
}  


#####################################################################################################
#####################################################################################################
#####################################################################################################

## Create squares and pairwise products of columns of a data frame.

#' Squares and two-way interactions of variables
#'
#' @description 
#' Generates squares and/or two-way interactions (pairwise products) of the columns of a data frame.
#' @param data_frame Data frame containing the variables whose squares and interactions are to be created.
#' @param is_square If \code{TRUE}, square of each column of \code{data_frame} is created.
#' @param is_inter If \code{TRUE}, product of every pair of columns of \code{data_frame} is created.
#' @param keep_marginal If \code{TRUE}, the original columns of \code{data_frame} are retained 
#' in the resulting data frame.
#' @return A data frame containing the squares and/or pairwise products of \code{data_frame}.
#' @export
#' @author Ambarish Chattopadhyay, Carl N. Morris and Jose R. Zubizarreta.
#' @examples
#' # Consider a data frame with N = 12 units and 2 covariates.
#' data_frame_sample = data.frame(male = c(rep(1,6),rep(0,6)), 
#' age = c(20,30,40,40,50,60,20,30,40,40,50,60))
#' # Get a data frame with all possible squares and first order interactions.
#' make_sq_inter(data_frame = data_frame_sample, is_square = TRUE, is_inter = TRUE, keep_marginal = FALSE)


make_sq_inter = function(data_frame, is_square = TRUE, is_inter = TRUE, keep_marginal = TRUE)
{
  
  if(keep_marginal == TRUE)
  {
    data_frame_appended = data_frame
  }
  
  
  if(is_inter == TRUE)
  {
    data_frame_inter = matrix(rep(0, nrow(data_frame)* choose(ncol(data_frame),2)), nrow = nrow(data_frame))
    
    # introduce a pseudo outcome and run a linear model
    
    y.pseudo = data_frame[,1]
    data_frame_with.outcome = data.frame(y.pseudo = y.pseudo, data_frame)
    
    data_frame_inter = data.frame(model.matrix(y.pseudo~.*.-1, data = data_frame_with.outcome))
    names.inter = colnames(data_frame_inter)
    
    if(keep_marginal == TRUE)
    {
      data_frame_appended = data_frame_inter
    }
    if(keep_marginal == FALSE)
    {
      data_frame_appended = as.data.frame(data_frame_inter[,-(1:ncol(data_frame))])
      colnames(data_frame_appended) <- names.inter[-(1:ncol(data_frame))]
    }
    
  }
  
  
  
  if(is_square == TRUE)
  {
    data_frame_square  = data_frame 
    for(j in 1:ncol(data_frame))
    {
      data_frame_square[,j] = (data_frame[,j])^2
      colnames(data_frame_square)[j] = paste(colnames(data_frame)[j],"^2", sep = "")
    }
    
    
    if(is_inter == TRUE)
    {
      data_frame_appended = cbind(data_frame_appended, data_frame_square)
    }
    if(is_inter == FALSE)
    {
      if(keep_marginal == TRUE)
      {
        data_frame_appended = cbind(data_frame_appended, data_frame_square)
      }
      if(keep_marginal == FALSE)
      {
        data_frame_appended = data_frame_square
      }
      
    }
    
    
  }
  
  
  return(as.data.frame(data_frame_appended))
}

#####################################################################################################
#####################################################################################################
#####################################################################################################

# Function to compute TASMD

tasmd_func = function(data_frame, index_col = FALSE, Z_1, Z_2, treat_lab = 1,
                      mean_tar = NULL, sd_tar = NULL, denom = 'target')
{
  # get the matrix of covariates X_mat
  if(index_col == TRUE)
  {
    X_mat = data_frame[,-1]
  }
  
  if(index_col == FALSE)
  {
    X_mat = data_frame
  }
  
  # separating the matrix of covariates for the given treatment group
  X_mat_t1 = X_mat[Z_1 == treat_lab,]
  X_mat_t2 = X_mat[Z_2 == treat_lab,]
  
  # column wise means of covariates in the two treatment groups
  means_t1 = apply(X_mat_t1, 2, mean)
  means_t2 = apply(X_mat_t2, 2, mean)
  
  
  if(is.null(mean_tar) == TRUE)
  {
    means_tar = apply(X_mat, 2, mean)
    if(denom == 'target')
    {
      sd_denom1 = apply(X_mat, 2, sd)
      sd_denom2 = apply(X_mat, 2, sd)
    }
    if(denom == 'group')
    {
      sd_denom1 = apply(X_mat_t1, 2, sd)
      sd_denom2 = apply(X_mat_t2, 2, sd)
    }
    
  }
  
  if(is.null(mean_tar) == FALSE)
  {
    means_tar = mean_tar
    if(denom == 'target')
    {
      sd_denom1 = sd_tar
      sd_denom2 = sd_tar
    }
    if(denom == 'group')
    {
      sd_denom1 = apply(X_mat_t1, 2, sd)
      sd_denom2 = apply(X_mat_t2, 2, sd)
    }
    
  }
  
  # calculate tasmd
  tasmd1 = abs(means_t1 - means_tar)/sd_denom1
  tasmd2 = abs(means_t2 - means_tar)/sd_denom2
  
  A = cbind(tasmd1, tasmd2)
  rownames(A) = colnames(X_mat)
  return(A)
}


#' Target Absolute Standardized Mean Differences (TASMD) 
#'
#' @description 
#' Computes the mean and standard deviation of Target Absolute Standardized Mean Differences (TASMD) of 
#' multiple covariates (or transformations thereof) in a treatment group relative to a target population 
#' or a target individual for a set of assignments under one or two designs. 
#' @param data_frame Data frame containing a column of unit indices (optional) and covariates 
#' (or transformations thereof).
#' @param index_col if \code{TRUE}, \code{data_frame} contains a column of unit indices.
#' @param alloc1 A matrix or vector of treatment assignments. If \code{alloc1} is a matrix, then each row
#' should correspond to an assignment vector. 
#' @param alloc2 A (optional) matrix or vector of treatment assignment. If \code{alloc2} is a matrix, then each row
#' should correspond to an assignment vector. 
#' @param treat_lab Label of the treatment group in which the TASMD is computed. 
#' @param mean_tar A (optional) vector of target profile of the covariates under consideration, 
#' e.g., mean of the covariates in the target population. Applicable only when \code{imbalance = 'TASMD'}. 
#' If \code{mean_tar = NULL}, the full-sample average of the covariates is considered as the target profile. 
#' @param sd_tar A optional vector of the standard deviation of the covariates in the target population. 
#' Applicable only when \code{imbalance = 'TASMD'}.
#' @param denom Specifies the denominator for the computation of TASMD. If \code{denom = 'target'}, 
#' the standard deviations of the covariates in the target population are used. If \code{denom = 'group'}, 
#' the standard deviations of the covariates in the treatment group given by \code{treat_lab} are used. 
#' Applicable only when \code{imbalance = 'TASMD'}. 
#' @param legend_text Legend of the two designs under consideration.
#' @param rounoff A number indicating the number of decimal places to be used for rounding off the TASMDs.
#' @return A list containing the following items (if \code{alloc1} and \code{alloc2} are matrices)
#' 
#' \code{tasmd_table}: A matrix containing the means (standard deviations in parenthesis) of the TASMDs
#' for the designs under consideration. If \code{alloc1} or \code{alloc2} is a vector, the
#' TASMD of the corresponding assignment is returned.
#' 
#' \code{tasmd_mean}: A matrix containing the means of the TASMDs for the designs under consideration.
#' 
#' \code{tasmd_sd}: A matrix containing the standard deviations of the TASMDs for the designs under consideration.
#' 
#' If \code{alloc1} and \code{alloc2} are vectors, \code{tasmd_rand} produces a data frame of the corresponding TASMDs.
#' @export
#' @author Ambarish Chattopadhyay, Carl N. Morris and Jose R. Zubizarreta.
#' @references 
#' Chattopadhyay, A., Morris, C. N., and Zubizarreta, J. R. (2020), ``Randomized and Balanced Allocation 
#' of Units into Treatment Groups Using the Finite Selection Model for \code{R}".
#' @examples
#' # Consider the Lalonde dataset.
#' # Get the full sample size.
#' N = nrow(Lalonde)
#' # Get the treatment group sizes.
#' n1 = floor(N/2)
#' n2 = N-n1
#' # Generate an SOM.
#' som_obs = som(n_treat = 2, treat_sizes = c(n1,n2),include_discard = FALSE,
#' method = 'SCOMARS', marginal_treat = rep((n2/N),N), control = FALSE)
#' # Generate a treatment assignment given som_obs.
#' f = fsm(data_frame = Lalonde, SOM = som_obs, s_function = 'Dopt', eps = 0.0001, 
#' ties = 'random', intercept = TRUE, standardize = TRUE, units_print = FALSE)
#' # Get assignment vector under the FSM.
#' Z_fsm_obs = f$data_frame_allocated$Treat
#' # Draw a random CRD.
#' Z_crd = crd(data_frame = Lalonde, n_treat = 2, treat_sizes = c(n1, n2), 
#' control = FALSE)$Treat
#' # Calculate the TASMD.
#' TASMD = tasmd_rand(data_frame = Lalonde, index_col = TRUE, alloc1 = Z_crd, 
#' alloc2 = Z_fsm, treat_lab = 1, mean_tar = NULL, sd_tar = NULL, 
#' denom = 'target', legend = c('CRD','FSM'), roundoff = 3)


tasmd_rand = function(data_frame, index_col = FALSE, alloc1, alloc2,
                      treat_lab = 1, legend = c('CRD', 'FSM'),
                      mean_tar = NULL, sd_tar = NULL, denom = 'target', roundoff = 3)
{
  # get number of iterations
  r = nrow(alloc1)
  
  # get the no. of covariates
  if(index_col == TRUE)
  {
    n.cov = ncol(data_frame[,-1])
    var.names = colnames(data_frame[,-1])
  }
  
  if(index_col == FALSE)
  {
    n.cov = ncol(data_frame)
    var.names = colnames(data_frame)
  }
  
  if(is.null(nrow(alloc1)) == TRUE)
  {
    imb = tasmd_func(data_frame = data_frame, index_col = index_col, Z_1 = alloc1, Z_2 = alloc2, 
                     treat_lab = treat_lab,
                     mean_tar = mean_tar, sd_tar = sd_tar, denom = denom)
    colnames(imb) = legend
    return(round(imb,roundoff))
  }
  
  if(is.null(nrow(alloc1)) == FALSE)
  {
    tasmd_iter1 = matrix(rep(0,r*n.cov),ncol = r)
    tasmd_iter2 = matrix(rep(0,r*n.cov),ncol = r)
    
    for(i in 1:r)
    {
      imb = tasmd_func(data_frame = data_frame, index_col = index_col, Z_1 = alloc1[i,], Z_2 = alloc2[i,], 
                       treat_lab = treat_lab,
                       mean_tar = mean_tar, sd_tar = sd_tar, denom = denom)
      tasmd_iter1[,i] = as.vector(imb[,1])
      tasmd_iter2[,i] = as.vector(imb[,2])
      
    }
    
    # Calculate mean and sd of the tasmd values
    tasmd_mean = round(cbind(rowMeans(tasmd_iter1), rowMeans(tasmd_iter2)), roundoff)
    tasmd_sd = round(cbind(apply(tasmd_iter1,1,sd), apply(tasmd_iter2,1,sd)), roundoff)
    
    tasmd_table = matrix(rep(0,2*n.cov),ncol = 2)
    
    for(i in 1:n.cov)
    {
      for(j in 1:2)
      {
        tasmd_table[i,j] = paste(tasmd_mean[i,j],"(",tasmd_sd[i,j],")", sep = "")
      }
    }
    
    rownames(tasmd_table) = var.names
    colnames(tasmd_table) = legend
    rownames(tasmd_mean) = var.names
    colnames(tasmd_mean) = legend
    rownames(tasmd_sd) = var.names
    colnames(tasmd_sd) = legend
    
    return(list(tasmd_table = tasmd_table, tasmd_mean = tasmd_mean, tasmd_sd = tasmd_sd))
    
  }
  
}

#####################################################################################################
#####################################################################################################
#####################################################################################################


# Function to compute Model-based ESS 


#' Model-based Effective Sample Size (ESS)
#'
#' @description 
#' Computes the model-based effective sample size (ESS) of a collection of assignments 
#' under a given set of potential outcomes.
#' @param X_cov  A matrix of covariates or transformations thereof that will be used
#'  as explanatory variables in the linear outcome models within each treatment group.
#' @param assign_matrix A matrix containing a collection of treatment assignment vectors, each column
#' containing a particular assignment vector. 
#' @param Y_mat A matrix of potential outcomes, where rows represent units and columns represent treatment
#' levels (ordered). 
#' @param contrast A vector of the coefficients of the treatment contrast of interest. For example, for estimating the
#' average treatment effect of treatment 1 versus treatment 2, \code{contrast = c(1,-1)}.
#' @return A vector of effective sample sizes for the given collection of assignments.
#' @export
#' @author Ambarish Chattopadhyay, Carl N. Morris and Jose R. Zubizarreta.
#' @references 
#' Chattopadhyay, A., Morris, C. N., and Zubizarreta, J. R. (2020), ``Randomized and Balanced Allocation 
#' of Units into Treatment Groups Using the Finite Selection Model for \code{R}".
#' @examples
#' # Consider the Lalonde dataset.
#' # Get the full sample size.
#' N = nrow(Lalonde)
#' # Get the treatment group sizes.
#' n1 = floor(N/2)
#' n2 = N-n1
#' # Generate an SOM.
#' som_obs = som(n_treat = 2, treat_sizes = c(n1,n2),include_discard = FALSE,
#' method = 'SCOMARS', marginal_treat = rep((n2/N),N), control = FALSE)
#' # Generate a treatment assignment given som_obs.
#' f = fsm(data_frame = Lalonde, SOM = som_obs, s_function = 'Dopt', eps = 0.0001, 
#' ties = 'random', intercept = TRUE, standardize = TRUE, units_print = FALSE)
#' # Get assignment vector under the FSM.
#' Z_fsm_obs = f$data_frame_allocated$Treat
#' # Draw a random CRD.
#' Z_crd_obs = crd(data_frame = Lalonde, n_treat = 2, treat_sizes = c(n1, n2), 
#' control = FALSE)$Treat
#' Z_big = cbind(Z_crd_obs, Z_fsm_obs)
#' # Generate the potential outcomes.
#' Y_1 = 100 - Lalonde$Age + 6 * Lalonde$Education - 20 * Lalonde$Black + 
#' 20 * Lalonde$Hispanic + 0.003 * Lalonde$Re75 + rnorm(N,0,4)
#' Y_1 = round(Y_1,2)
#' # Set unit level causal effect = tau = 0.
#' tau = 0
#' Y_2 = Y_1 + tau
#' # Get the matrix of potential outcomes.
#' Y_appended = cbind(Y_1, Y_2)
#' # Get the matrix of covariates.
#' X_cov = Lalonde[,-1]
#' ess = ess_model(X_cov = X_cov, assign_matrix = Z_big, Y_mat = Y_appended, contrast = c(1,-1))


ess_model = function(X_cov, assign_matrix, Y_mat, contrast = c(1,-1))
{
  # Extract treatment labels
  N = nrow(X_cov)
  treat.labs = sort(unique(assign_matrix[,1]))
  
  # de-meaned columns of X
  X_cov_demean = X_cov
  for(j in 1:ncol(X_cov))
  {
    X_cov_demean[,j] = X_cov[,j] - mean(X_cov[,j])
  }
  
  X_cov_demean = as.matrix(X_cov_demean)
  
  #### Matrix of assignments
  Z.big = assign_matrix
  
  ## Fit the model after demeaning
  V_est = rep(0,ncol(Z.big))
  N_eff = rep(0,ncol(Z.big))
  
  ## Calculate the model-based variances for the ith design
  for(i in 1:ncol(Z.big))
  {
    var_coef = rep(0, length(treat.labs))
    
    # calculate the variance of E^{hat}[Y(s)]
    for(s in 1:length(treat.labs))
    {
      Y_pot = Y_mat[,s]
      fit_s =  lm(Y_pot[Z.big[,i] == treat.labs[s]] ~ X_cov_demean[Z.big[,i] == treat.labs[s],])
      var_coef[s] = (coef(summary(fit_s))[1, "Std. Error"])^2
    }
    
    # variance of the estimator
    V_est[i] = sum((contrast^2)*var_coef)
    
  }
  
  N_eff = N/(V_est/min(V_est))
  
  return(N_eff)  
}    

######################################################################################################
######################################################################################################
######################################################################################################


#' Randomization-based Effective Sample Size (ESS)
#'
#' @description 
#' Computes the randomization-based effective sample size (ESS) of a collection of assignments 
#' under a given set of potential outcomes.
#' @param assign_array A three dimensional array containing a set of independent realizations of a
#' collection the designs. The first coordinate of the array represents the iterations for
#' each design. The second coordinate represents the units. The third coordinate represents the design.
#' @param Y_mat A matrix of potential outcomes, where rows represent units and columns represent treatment
#' levels (ordered). 
#' @param contrast A vector of the coefficients of the treatment contrast of interest. For example, for estimating the
#' average treatment effect of treatment 1 versus treatment 2, \code{contrast = c(1,-1)}.
#' @return A vector of effective sample sizes for the given collection of assignments.
#' @export
#' @author Ambarish Chattopadhyay, Carl N. Morris and Jose R. Zubizarreta.
#' @references 
#' Chattopadhyay, A., Morris, C. N., and Zubizarreta, J. R. (2020), ``Randomized and Balanced Allocation 
#' of Units into Treatment Groups Using the Finite Selection Model for \code{R}".
#' @examples
#' # Consider N = 12, n1 = n2 = 6.
#' df_sample = data.frame(index = 1:12, x = c(20,30,40,40,50,60,20,30,40,40,50,60))
#' # Generate the potential outcomes.
#' Y_1 = 100 + (df_sample$x - mean(df_sample$x)) + rnorm(12, 0, 4)
#' Y_2 = Y_1_true + 50
#' # Create matrix of potential outcomes.
#' Y_appended = cbind(Y_1, Y_2)
#' # Generate 1000 assignments under CRD and the FSM.
#' Z_crd_iter = matrix(rep(0, 1000 * 12), nrow = 1000)
#' Z_fsm_iter = matrix(rep(0, 1000 * 12), nrow = 1000)
#' for(i in 1:1000)
#' {
#' # Generate an assignment vector under CRD.
#' fc = crd(data_frame = df_sample, n_treat = 2, treat_sizes = c(6,6), control = FALSE)
#' Z_crd_iter[i,] = fc$Treat
#' # Generate an assignment vector under the FSM.
#' som_iter = som(data_frame = NULL, n_treat = 2, 
#' treat_sizes = c(6, 6),include_discard = FALSE,
#' method = 'SCOMARS', marginal_treat = rep((6/12), 12), control = FALSE)
#' f = fsm(data_frame = df_sample, SOM = som_iter, s_function = 'Dopt',eps = 0.0001, 
#' ties = 'random', intercept = TRUE, standardize = TRUE, units_print = FALSE)
#' Z_fsm_iter[i,] = f$data_frame_allocated$Treat
#' }
#' # Create a 3-dim array of assignments.
#' Z_array = array(0,dim=c(1000, 12, 2))
#' Z_array[,,1] = Z_crd_iter
#' Z_array[,,2] = Z_fsm_iter
#' # Calculate the ESS.
#' ess_rand(assign_array = Z_array, Y_mat = Y_appended, contrast = c(1,-1))



## Function for calculating randomization-based ESS 

# assign_array = 3rd coordinate indicates design, 1st coordinate indicates iterations, 2nd coordinate indicates unit 

ess_rand = function(assign_array, Y_mat, contrast = c(1,-1))
{
  # Extract treatment labels.
  N = nrow(Y_mat)
  n_designs = dim(assign_array)[3]
  level = sort(unique(assign_array[1,,1]))  
  
  Y1 = Y_mat[,1]
  Y2 = Y_mat[,2]
  
  n.iter = dim(assign_array)[1]
  # Initialize matrix of estimators - rows indicate iterations, columns indicate designs
  T_est = matrix(rep(0,n.iter*n_designs), nrow = n.iter)
  
  for(i in 1:n.iter)
  {
    for(j in 1: n_designs)
    {
      Z_iter = assign_array[i,,j]
      # Compute the mean response in each treatment group
      T_mean = c(mean(Y1[Z_iter == level[1]]), mean(Y2[Z_iter == level[2]]))
      T_est[i,j] = sum(T_mean * contrast)
      
    }
  }    
  # Calculate the randomization-based variance of the estimator under each design.
  V_est = apply(T_est, 2, var)
  # Calculate the ESS.
  N_eff = N/(V_est/min(V_est))
  
  return(N_eff)  
  
}












