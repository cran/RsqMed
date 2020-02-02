#'  Function to calculate the Rsq function as a total mediation effect size measure (Gaussian outcome only). If method='iSIS', a two-step procedure is performed, where the first step filters the non-mediators based on part of the data and the second step calculates the point estimates for Rsq using random-effect models on the remaining data. If method='ALL', Rsq is calculated based on all subjects and variables. 
#'
#' @param p Proportion of the training dataset for selecting mediators regarding the whole dataset, default is set as 1/2. If method='ALL', keep p at default. 
#' @param outcome Vector of outcome type of interest; Only Gaussian distributed outcome is accepted.
#' @param med Matrix of putative mediators
#' @param covar Covariate matrix
#' @param indp Vector of the independent variable of interest, e.g. environmental variable
#' @param method Method used to screen out non-mediators. When no variable selection is required, method='ALL'; otherwise, iterative sure independence screening (SIS) is used for variable selection, i.e., method='iSIS'. Note that when method='ALL', no screening is performed, i.e., the Rsq measure is calculated on all data and all variables included.  
#' @param iter.max Maximum number of iteration used in iSIS, default=10 (please refer the SIS package for detail explanation)
#' @param nsis Number of variables recruited by iterative SIS
#' @param screening T if filtering mediators based on the strength of independent variable and mediators as a preprocessing step; F if all putative mediators are included, default=F.
#' @param init.cutoff The percentage of mediators remaining after the screening step.
#' @return Output vector consist of Rsq mediated(Rsq.mediated), shared over simple effects (SOS), number of selected mediators (pab), and the Rsq that used to calculate the Rsq measure: variance of outcome explained by mediator (Rsq.YM), variance of outcome explained by the independent variable (Rsq.YX), and variance of outcome explained by mediator and independent variable (Rsq.YMX), n is the sample size based on which the random effect models are fitted. 
#' @return Name of selected mediators (select)
#' @export
#' @examples{
#'\donttest{
#' data(example)
#' attach(example)
#'Rsq.measure(p=1/2, outcome=Y, med=M,covar=Cov,indp=X,method='iSIS', iter.max=1)
#'}
#'}
Rsq.measure<-function(p=1/2, outcome,med,covar,indp,method=c('iSIS','ALL'), iter.max=10, nsis=NULL, init.cutoff=0.1, screening=FALSE){

   if (is.null(nsis)) nsis=length(outcome)/log(length(outcome))
   # run variable selection
   # take the first p proportion as the training dataset
   train <- 1:(nrow(med)*p)
   if (min(train) < 1) stop('Please specifiy larger training dataset!')
   if (method=='ALL' & screening==T) message('Note: the analysis is performed on all data and all included variables. No screening is performed.')
   # standardized Med and independent variable
   Med<-apply(med,2,scale,center=T, scale=T)
   indp.std <- scale(indp[train],center=T, scale=T)
   if (!is.null(covar)) {covar <- as.matrix(covar); covar.train <- covar[train,]; covar.test <- covar[-train,]} 


   # screen by FDR
   cal_alpha_simple<-function(x){
   if (!is.null(covar)) data2<-data.frame(Med=x,envir=indp.std, cov=covar.train ) else data2 <- data.frame(Med=x,envir=indp.std)
   l<-summary(stats::lm('Med ~.', data=data2))
   invalphap<-(l$coefficients)['envir','Pr(>|t|)']
   return(invalphap)
   }
   # keep the top candidates for ultra high dimensional data
   if(screening) {
   inv.alpha.p<-apply(Med[train,],2, cal_alpha_simple)
   order <- order(stats::p.adjust(inv.alpha.p, method='fdr'))[1:round(init.cutoff *ncol(Med)) ]
   } else order <- 1:ncol(Med)

   Med <- Med[,order, drop=F]

   if (!is.null(covar)) tdat <- data.frame(y=outcome[train], cov=cbind(indp.std, covar.train)) else tdat <- data.frame(y=outcome[train], cov=indp.std)
   f0 <- stats::lm(y~., data=tdat)
   res <- stats::residuals(f0)

   # iterative SIS is suggested for mediation analysis,
   X <- Med[train, ]
   if (method=='iSIS'){
   model1 <- invisible(SIS::SIS(x=X, y=res, family='gaussian',tune='bic',  seed=1234, penalty='MCP',nsis=nsis, iter.max=iter.max))
   pab <- length(model1$ix)
   select <- (model1$ix)
   outcome.test <- outcome[-train]
   indp.test <- indp[-train]
   } else {
   select <- 1:ncol(Med)
   pab <- ncol(Med)
   covar.test <- covar
   outcome.test <- outcome
   indp.test <- indp
   }

   #in the testing data
   if (!is.null(covar)){
   tdat <- data.frame(y=outcome.test, cov=covar.test)
   f0 <- stats::lm(y~. , data=tdat)
   res <- stats::residuals(f0)} else res <- outcome.test

   indp.std <- scale(indp.test,center=T, scale=T)

   testing <- data.frame(pheno=res,  envir=indp.std)
   f3<-tryCatch(stats::lm(pheno~envir, data=testing), error=function(c) NA, warning=function(c) NA)
   if (all(is.na(f3))) total<-Rsq.YX<-NA else {
   Rsq.YX=summary(f3)$adj.r.squared #******
   total<-f3$coefficients[2]
   names(total) <-NULL }
   sy<-stats::sd(res,na.rm=T)
   SST=sum((res-mean(res))^2)
   MST=SST/(nrow(testing)-1)


   if (pab==0) {
   output <- c(Rsq.mediated=0,SOS=0, pab=0, Rsq.YM=0, Rsq.YX=Rsq.YX, Rsq.YMX=Rsq.YX, total=total, abproduct=0, n=nrow(testing) )
   return(list(output=output, select=NULL))
   } else {

   if (method=='iSIS') Med <- Med[-train,select,drop=F] else Med <- Med

   #Y~M+X
   Kins <- as.matrix(Med)%*%t(as.matrix(Med))
   colnames(Kins) <- rownames(Kins) <- 1:nrow(testing)
   testing <- data.frame(cbind(res=res, indp.std=indp.std))
   testing$id <- 1:nrow(testing)
   fit1 <- tryCatch(GMMAT::glmmkin(res~indp.std, id='id', data=testing, kins=Kins, family=stats::gaussian(link='identity')), error=function(c) NA)
   if (all(is.na(fit1))) { tau<-phi<-direct<-NA } else { if (fit1$converged) {tau<-fit1$theta[2] ; phi<-fit1$theta[1];  direct<-fit1$coefficients[2]; } else tau<-phi<-direct<-NA }
   Rsq.YMX <- 1-phi/MST

   #Y~M
   fit2 <- tryCatch(GMMAT::glmmkin(res~1, data=testing, kins=Kins,id='id', family=stats::gaussian(link='identity')), error=function(c) NA)
   if (all(is.na(fit2))) { tau2<-phi2<-NA } else {
   if (fit2$converged) {tau2<-fit2$theta[2] ; phi2=fit2$theta[1];} else tau2<-phi2<-NA }
   Rsq.YM <- 1-phi2/MST

   Rsq.mediated <- Rsq.YM+Rsq.YX-Rsq.YMX
   SOS <-  Rsq.mediated/Rsq.YX
   abproduct= as.numeric(total-direct)
   output <- c(Rsq.mediated=Rsq.mediated,  SOS=SOS, pab=pab, Rsq.YM=Rsq.YM, Rsq.YX=Rsq.YX, Rsq.YMX=Rsq.YMX, total=total,  abprod=abproduct, n=nrow(testing) )

    if (!is.null(colnames(Med))) select <- 'column name of Mediators are not specified' else select <- colnames(Med)

   #select is the selected mediators, use med[,select] to look at the selected mediators (for future pathway analysis)
   return(list(output=output, select=select))
}
}


