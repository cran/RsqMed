#' Functions to generate the confidence interval of the Rsq measure using nonparametric bootstrap. 
#'
#' @param p Proportion of the training dataset for selecting mediators regarding to the whole dataset, default is set as 1/2
#' @param outcome Vector of outcome type of interest; Only Gaussian distributed outcome is accepted.
#' @param med Matrix of putative mediators
#' @param covar Covariate matrix
#' @param indp Vector of the independent variable of interest, e.g. environmental exposure
#' @param method Method used to screen out non-mediators. When no variable selection is required, method='ALL'; otherwise, iterative sure independence screening (SIS) is used for variable selection, i.e., method='iSIS'.
#' @param B Number of bootstrap samples, default is 100
#' @param iter.max Maximum number of iteration used in iSIS, default=10 (details see the SIS package)
#' @param nsis Number of pedictors recruited by iSIS
#' @param screening T: filtering mediators based on the strength of indp and mediators as a preprocessing step; F: all putative mediators are included, default=F.
#' @param init.cutoff The percentage of mediators remaining after the screening step.
#' @return CI: The 95 percent confidence intervals of Rsq measure (Rsq.mediated), shared over simple effects (SOS), number of mediators selected (pab), variance of outcome explained by mediator (Rsq.YM), variance of outcome explained by the independent variable (Rsq.YX), and variance of outcome explained by mediator and independent variable (Rsq.YMX). The estimates for each bootstrap are also returned.
#' @export
#' @examples {
#'  \donttest{
#' data(example)
#' attach(example)
#' CI.Rsq.measure( p=1/2, outcome=Y,med=M[,1:2],covar=NULL,indp=X,method='ALL', B=1)
#' }
#'}
CI.Rsq.measure <- function( p=1/2, outcome,med,covar,indp, method=c('iSIS','All'),B=200, iter.max=10, nsis=NULL, init.cutoff=0.1, screening=F){
if (is.null(nsis)) nsis=length(outcome)/log(length(outcome))
output <- matrix(NA, nrow=B, ncol=9)
for (i in 1:B){
set.seed(i)
ID <- sort(sample(1:length(outcome),length(outcome), replace=T))
outcome_cl <- outcome[ID]
med_cl <- med[ID,]
if (!is.null(covar)) covar_cl <- covar[ID,] else covar_cl <- NULL
indp_cl <- indp[ID]
output[i,] <- Rsq.measure(p=p,  outcome=outcome_cl ,med=med_cl ,covar=covar_cl, indp=indp_cl , method=method, iter.max=iter.max, nsis=nsis, init.cutoff=init.cutoff, screening=screening)$output
}
output <- as.data.frame(output)
colnames(output) <- c('Rsq.mediated','SOS','pab','Rsq.YM','Rsq.YX','Rsq.YMX','total','abprod','n')
CI <- apply(output,2,stats::quantile,probs=c(0.025,0.975),na.rm=T)
rownames(CI) <- c('2.5%','97.5%')
return(list(CI=CI, output=output))

}

