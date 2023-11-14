#'  Function to calculate the Rsq function as a total effect size measure for mediation effect using cross-fitted estimation
#'
#' @param Y vector of the outcome of interest; outcome has to follow a Gaussian distribution.
#' @param M matrix of putative mediators
#' @param Covar covariates matrix
#' @param X vector of the independent variable of interest, e.g. environmental variable
#' @param iter.max Maximum number of iterations used in iSIS, default = 3 (details see the SIS package).
#' @param nsis Number of predictors recruited by iSIS, default = NULL
#' @param first.half TRUE: split sample into two halves by the order in the dataset. FALSE: randomly split samples into halves, default = TRUE.
#' @param seed Random seed used for sample splitting, default = 2022.
#' @param tune Method for tuning the regularization parameter of the penalized likelihood subproblems and of the final model selected by (i)SIS. Options include tune = 'bic' and tune = 'aic'.
#' @param penalty The penalty to be applied in the regularized likelihood subproblems. 'MCP', and 'lasso' are provided. 'MCP' is recommended.
#' @return Output Vector consisting of Rsq mediated(Rsq.mediated), Lower confidence bound constructed by the asymptotic variance (CI_asym_low), Upper confidence bound constructed by the asymptotic variance (CI_asym_up), Lower confidence bound constructed by the conservative variance (CI_cons_low), Upper confidence bound constructed by the conservative variance (CI_cons_up), number of selected mediators in subsample 1 (pab1), number of selected mediators in subsample 2 (pab2), and the Rsq that used to calculate the Rsq measure: variance of outcome explained by mediator (Rsq.YM), variance of outcome explained by the independent variable (Rsq.YX), and variance of outcome explained by mediator and independent variable (Rsq.YMX); Sample Size in analysis (Sample Size)
#' @return Name of selected mediators in subsample 1 (select1)
#' @return Name of selected mediators in subsample 2 (select2)
#' @export
#' @examples{
#'\donttest{
#' data(example)
#' attach(example)
#'CF_Rsq.measure(Y=Y, M=M, X=X, tune = "bic", penalty = "MCP")
#'}
#'}

CF_Rsq.measure <- function(Y, M, Covar = NULL, X, iter.max = 3, nsis = NULL, first.half = TRUE, seed = 2022,
                           tune = c("aic", "bic"), penalty = c("MCP", "lasso")){


   if (length(Y)!=nrow(M)) stop('Sample sizes do not match.') 
   if (length(Y)!=length(X)) stop('Sample sizes do not match.') 
   if (nrow(M)!=length(X)) stop('Sample sizes do not match.')
   if (!is.null(Covar)) {Covar <- as.matrix(Covar)
   if (nrow(Covar)!=length(X)) stop('Sample sizes do not match.')}
   if (length(unique(Y)) < 3) message('Warning: Current algorithm can only deal with continuous outcome.')
   if (is.vector(M)) stop('This algorithm does not support single mediator model.')
   if (is.null(colnames(M))) {
    	message("Column name of Mediators are not specified. Naming the M by its position in the dataset. ") 
    	colnames(M) <- paste0('M', 1:ncol(M))
    	}
   if (sum(is.na(Y)) + sum(is.na(M)) + sum(is.na(Covar)) + sum(is.na(X)) > 0) stop('Algorithm cannot deal with missing data! Impute or subset the data first. ')
  	
  # Whether use the first half or random sample split for variable selection
  if (first.half == TRUE) {
    idx1 <- 1:(nrow(M) * 1/2)
  } else {
    set.seed(seed)
    idx1 <- sample(1:nrow(M), ceiling(nrow(M)/2), replace = FALSE)
  }
  
  if (max(idx1) < 20){stop("Please specifiy a larger dataset.")}
  if (is.null(tune)){stop("Please specifiy the tuning parameter.")}
 
  d <- ncol(M)
  n <- nrow(M)
  
  # ---- iSIS variable selection ----

  # Scale M matrix
  M_res_1 <- apply(M, 2, scale)[idx1, ]
  M_res_2 <- apply(M, 2, scale)[-idx1, ]

  # Regress covariates and X out
  if (!is.null(Covar)){
    tdat_1 <- data.frame(y = Y[idx1], 
                         envir = scale(X[idx1]),
                         cov = Covar[idx1, ])
    f_1 <- stats::lm(y ~ ., data = tdat_1)
  } else{
    tdat_1 <- data.frame(y = Y[idx1], 
                         envir = scale(X[idx1]))
    f_1 <- stats::lm(y ~ ., data = tdat_1)
  }
  Y_res_1 <- stats::residuals(f_1)  
  
  if (!is.null(Covar)){
    tdat_2 <- data.frame(y = Y[-idx1], 
                         envir = scale(X[-idx1]),
                         cov = Covar[-idx1, ])
    f_2 <- stats::lm(y ~ ., data = tdat_2)
  } else{
    tdat_2 <- data.frame(y = Y[-idx1], 
                         envir = scale(X[-idx1]))
    f_2 <- stats::lm(y ~ ., data = tdat_2)
  }
  Y_res_2 <- stats::residuals(f_2) 

  # Regress Covariates out
  if (!is.null(Covar)){
    tdat_3 <- data.frame(y = Y[idx1], 
                         cov = Covar[idx1, ])
    f_3 <- stats::lm(y ~ ., data = tdat_3)
    Y_res_3 <- stats::residuals(f_3) 
  } else{
    Y_res_3 <- Y[idx1]
  }
  
  if (!is.null(Covar)){
    tdat_4 <- data.frame(y = Y[-idx1], 
                         cov = Covar[-idx1, ])
    f_4 <- stats::lm(y ~ ., data = tdat_4)
    Y_res_4 <- stats::residuals(f_4)  
  } else{
    Y_res_4 <- Y[-idx1]
  }
 
  # Scale the independent variable
  X_res_1 <- as.numeric(scale(X[idx1]))
  X_res_2 <- as.numeric(scale(X[-idx1]))

  # ---- iSIS ----
  model1 <- invisible(SIS::SIS(x = M_res_1, y = Y_res_1, 
                               family = 'gaussian', tune = tune,  seed = 1234, 
                               penalty = penalty,
                               nsis = nsis,
                               iter.max = iter.max))
  pab_1 <- length(model1$ix)
  m1 <- model1$ix

# ------ Estimation Procedure ------
  if (length(m1) > 0){
    ols1.YW <- stats::lm(Y_res_4 ~ cbind(M_res_2, X_res_2)[, c(m1, d + 1)])   # Y ~ M + X 
    err_yw1 <- ols1.YW$residuals
    
    ols1.YZ <- stats::lm(Y_res_4 ~ M_res_2[, m1])                             # Y ~ M
    err_yz1 <- ols1.YZ$residuals
    
    ols1.YX <- stats::lm(Y_res_4 ~ X_res_2)                                   # Y ~ X
    err_yx1 <- ols1.YX$residuals
    
    RYX_1 <- summary(ols1.YX)$adj.r.squared
    V.YW1 <- sum(ols1.YW$residuals^2)/ols1.YW$df.residual                       
    V.YZ1 <- sum(ols1.YZ$residuals^2)/ols1.YZ$df.residual
    V.YX1 <- sum(ols1.YX$residuals^2)/ols1.YX$df.residual
    
    err_y1 <- Y_res_4 - mean(Y_res_4)
    
    v_yx1 <- mean(err_yx1^2)
    v_yw1 <- mean(err_yw1^2)
    v_yz1 <- mean(err_yz1^2)
    v_y1 <- mean(err_y1^2)
    
    err1 <- cbind(err_yx1^2, err_yz1^2, err_yw1^2, err_y1^2)
    A1 <- stats::cov(err1)
    
  }else{
    message("There is no mediators selected in the 1st half")
  }

# ------ Subset 2 Estimation ------

# ------ iSIS
  model2 <- invisible(SIS::SIS(x = M_res_2, y = Y_res_2, 
                             family = 'gaussian', tune = tune,  seed = 1234, 
                             penalty = penalty,
                             nsis = nsis,
                             iter.max = iter.max))
  pab_2 <- length(model2$ix)
  m2 <- model2$ix


# -----Estimation Process 
  if (length(m2) > 0) {
    ols2.YW <- stats::lm(Y_res_3 ~ cbind(M_res_1, X_res_1)[, c(m2, d + 1)])
    err_yw2 <- ols2.YW$residuals
    
    ols2.YZ <- stats::lm(Y_res_3 ~ M_res_1[, m2])
    err_yz2 <- ols2.YZ$residuals
    
    ols2.YX <- stats::lm(Y_res_3 ~ X_res_1)
    err_yx2 <- ols2.YX$residuals
    
    RYX_2 <- summary(ols2.YX)$adj.r.squared
    V.YW2 <- sum(ols2.YW$residuals^2)/ols2.YW$df.residual
    V.YZ2 <- sum(ols2.YZ$residuals^2)/ols2.YZ$df.residual
    V.YX2 <- sum(ols2.YX$residuals^2)/ols2.YX$df.residual

    err_y2 <- Y_res_3 - mean(Y_res_3)
    
    v_yx2 <- mean(err_yx2^2)
    v_yw2 <- mean(err_yw2^2)
    v_yz2 <- mean(err_yz2^2)
    v_y2 <- mean(err_y2^2)
    
    err2 <- cbind(err_yx2^2, err_yz2^2, err_yw2^2, err_y2^2)
    A2 <- stats::cov(err2)
    
  }else{
    message("There is no mediators selected in the 2nd half")
  }
  
  if(length(m1) > 0 & length(m2) > 0){
    A <- 0.5 * (A1 + A2)
    v_yw <- 0.5 * (v_yw1 + v_yw2)
    v_yz <- 0.5 * (v_yz1 + v_yz2)
    v_yx <- 0.5 * (v_yx1 + v_yx2)
    v_y <- stats::var(c(Y_res_3, Y_res_4))
    
    a <- c(-1/v_y, -1/v_y, 1/v_y, (v_yx + v_yz - v_yw)/v_y^2)
    v <- t(a) %*% A %*% a
    
    V.YW <- 0.5 * (V.YW1 + V.YW2)
    V.YZ <- 0.5 * (V.YZ1 + V.YZ2)
    
    # Compute the R2
    Rsq.mediated <- 1.0 - (v_yx + v_yz - v_yw) / v_y
    CI_width_asym <- stats::qnorm(0.975) * sqrt(v) / sqrt(n)
    v_asym <- sqrt(v) / sqrt(n)
    
    
    # V.YX <- 0.5 * (V.YX1 + V.YX2) # 1
    ols.YX <- stats::lm(c(Y_res_3, Y_res_4) ~ c(X_res_1, X_res_2)) 
    V.YX <- sum(ols.YX$residuals^2)/ols.YX$df.residual
    RYX_12 <- summary(ols.YX)$adj.r.squared
    V.Y  <- stats::var(c(Y_res_3, Y_res_4))
  } else{
    stop("No mediators selected.")
  }
    
  output <- c(Rsq.mediated = Rsq.mediated, 
              CI_asym_low = Rsq.mediated - CI_width_asym, CI_asym_up = Rsq.mediated + CI_width_asym, 
              pab1 = round(pab_1, 0), pab2 = round(pab_2, 0), 
              Rsq.YX = RYX_12, Rsq.YMX = 1 - V.YW/V.Y, Rsq.YM = 1 - V.YZ/V.Y,
              SampleSize = n)
  
  Med1 <- M[-idx1, m1, drop = F]
  Med2 <- M[idx1, m2, drop = F]
  
  if (is.null(colnames(Med1)) | is.null(colnames(Med2))){
    select1 <- "Column name of Mediators are not specified."
    select2 <- "Column name of Mediators are not specified."
    } else{
    select1 <- colnames(Med1)
    select2 <- colnames(Med2)
    }
  
  return(list(output = output,
              select1 = select1,
              select2 = select2))

}




