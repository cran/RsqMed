#' Functions to generate the confidence interval of the Rsq measure using nonparametric bootstrap. 
#'
#' @param p Proportion of the training dataset for selecting mediators regarding to the whole dataset, default is set as 1/2
#' @param Y Vector of outcome type of interest; Only Gaussian distributed outcome is accepted.
#' @param M Matrix of putative mediators
#' @param Covar Covariate matrix
#' @param X Vector of the exposure or independent variable of interest, e.g. environmental exposure
#' @param method Method used to screen out non-mediators. When no variable selection is required, method='ALL'; otherwise, iterative sure independence screening (SIS) is used for variable selection, i.e., method='iSIS'.
#' @param B Number of bootstrap samples, default is 100
#' @param iter.max Maximum number of iteration used in iSIS, default=10 (details see the SIS package)
#' @param nsis Number of pedictors recruited by iSIS
#' @param filtering T: filtering mediators based on the strength of indp and mediators as a preprocessing step; F: all putative mediators are included, default=F.
#' @param init.FDR.cutoff FDR threshold for the filtering method on M1 type of mediators.
#' @return CI: The 95 percent confidence intervals of Rsq measure (Rsq.mediated), shared over simple effects (SOS), number of mediators selected (pab), variance of outcome explained by mediator (Rsq.YM), variance of outcome explained by the independent variable (Rsq.YX), and variance of outcome explained by mediator and independent variable (Rsq.YMX). The estimates for each bootstrap are also returned.
#' @export
#' @examples {
#'  \dontrun{
#' data(example)
#' attach(example)
#' CI.Rsq.measure( p=1/2, Y=Y,M=M,X=X,method='ALL', B=1, iter.max=1)
#' }
#'}
CI.Rsq.measure <- function( p=1/2, Y,M,Covar=NULL,X, method=c('iSIS','ALL'),B=200, iter.max=3, nsis=NULL, init.FDR.cutoff=0.1, filtering=F){
if (is.null(nsis)) nsis=length(Y)/log(length(Y))
output <- matrix(NA, nrow=B, ncol=10)
for (kk in 1:B){
set.seed(kk)
ID <- sample(1:length(Y),length(Y), replace=T)
#ID <- 1:length(Y)
outcome_cl <- Y[ID]
med_cl <- M[ID,]
if (!is.null(Covar)) covar_cl <- Covar[ID,] else covar_cl <- NULL
indp_cl <- X[ID]
output[kk,] <- Rsq.measure(p=p,  Y=outcome_cl ,M=med_cl ,Covar=covar_cl, X=indp_cl , method=method, iter.max=iter.max, nsis=nsis, init.FDR.cutoff=init.FDR.cutoff, filtering=filtering)$output
}
output <- as.data.frame(output)
colnames(output) <- c('Rsq.mediated','SOS','pab','Rsq.YM','Rsq.YX','Rsq.YMX','total','abprod','n.train','n.estimate')
CI <- apply(output,2,stats::quantile,probs=c(0.025,0.975),na.rm=T)
rownames(CI) <- c('2.5%','97.5%')
return(list(CI=CI, output=output))
}


