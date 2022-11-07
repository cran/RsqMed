#'  Function to calculate the Rsq function as a total mediation effect size measure (Gaussian outcome only). If method = 'iSIS', a two-step procedure is performed, where the first step filters the non-mediators based on the first p proportion of the data and the second step calculates the point estimates for Rsq using random-effect models on the remaining data. If method = 'ALL', Rsq is calculated based on all subjects and mediators (assuming all mediators are the true mediators). It is optional to adding filtering step on putative mediators to exclude M1 type of non-mediators (See Yang et al, BMC bioinformatics, 2021).
#'
#' @param p Proportion of the training dataset for selecting mediators regarding the whole dataset, default is set as 1/2. If method = 'ALL', all subjects are included. 
#' @param Y Vector of outcome; Only Gaussian distributed outcome is accepted.
#' @param M Matrix of putative mediators
#' @param Covar Covariate matrix or vector, default = NULL 
#' @param X Vector of the independent variable of interest, e.g. environmental variable
#' @param method Method used to screen out M2 type of non-mediators. When method = 'ALL', no variable selection is performed and Rsq is calculated on all data; otherwise, iterative sure independence screening (SIS) is performed on training dataset, i.e., method = 'iSIS'. 
#' @param iter.max Maximum number of iteration used in iSIS, default = 3 
#' @param nsis Number of variables recruited by iterative SIS (please refer the SIS package for detailed explanations)
#' @param filtering When filtering = T, filtering based on the strength of independent variable and mediators is performed; When filtering = F, no preprocessing is performed before variable selection. By default filtering = F.
#' @param init.FDR.cutoff FDR threshold for the filtering. 
#' @return Output Vector consists of Rsq mediated(Rsq.mediated), shared over simple effects (SOS), number of selected mediators (pab), and the Rsq that used to calculate the Rsq measure: variance of outcome explained by mediator (Rsq.YM), variance of outcome explained by the independent variable (Rsq.YX), and variance of outcome explained by mediator and independent variable (Rsq.YMX), n.train is the sample size on which variable selection is performed, n.estimate is the sample size based on which the mediation effect is estimated. 
#' @return Name of selected putative mediators (select). Note that M1 type of non-mediators may still be included in the model, but it would not impact the estimation of total mediation effect under certain assumptions. 
#' @export
#' @examples{
#'\dontrun{
#' data(example)
#' attach(example)
#'Rsq.measure(p=1/2, Y=Y, M=M,Covar=Cov,X=X,method='iSIS', iter.max=1, filtering=TRUE)
#'}
#'}
Rsq.measure<-function(p=1/2, Y,M,Covar=NULL,X,method=c('iSIS','ALL'), iter.max=3, nsis=NULL, init.FDR.cutoff=0.1, filtering=FALSE){

   if (length(Y)!=nrow(M)) stop('Sample sizes do not match.') 
   if (length(Y)!=length(X)) stop('Sample sizes do not match.') 
   if (nrow(M)!=length(X)) stop('Sample sizes do not match.')
   if (is.vector(M)) stop('This algorithm does not support single mediator model.')
   if (!is.null(Covar)) {Covar <- as.matrix(Covar)
   if (nrow(Covar)!=length(X)) stop('Sample sizes do not match.')}
   if (sum(is.na(Y)) + sum(is.na(M)) + sum(is.na(Covar)) + sum(is.na(X)) > 0) stop('Algorithm cannot deal with missing data! Impute or subset the data first. ')
   if (length(unique(Y)) < 3) message('Warning: Current algorithm can only deal with continuous outcome.')
     
   if (is.null(nsis)) nsis=length(Y)/log(length(Y))
   # run variable selection
   # take the first p proportion as the training dataset
   train <- 1:(nrow(M)*p)
   test <- (nrow(M)*p + 1):length(X)
   if (max(train) < 20) stop('Please specifiy a larger training dataset!')
   if (method == 'ALL' & filtering == F) message('Note: the analysis is performed on all included variables and data.')
   if (method =='ALL' & filtering == T) message('Note: Filtering on M1 type of non-mediators will not be performed. ')

   if (is.null(colnames(M))) {
    	message("Column name of putative Mediators are not specified. Naming the M by its position in the dataset. ") 
    	colnames(M) <- paste0('M', 1:ncol(M))
    	}
    	
  #adjusting for covariate variable
  #get the residual of y
   if (!is.null(Covar)) {tdat <- data.frame(y=Y, cov=Covar) 
   f0 <- stats::lm(y~., data=tdat)
   res.y <- stats::residuals(f0)  
   
   tdat <- data.frame(y=X, cov=Covar) 
   f0 <- stats::lm(y~., data=tdat)
   res.x <- stats::residuals(f0)   
     
   #get the residual of med
   res.med <- NULL
   for (i in 1:ncol(M)){
   m <- M[,i]
   tdat <- data.frame(y=m, cov=Covar)
   f1 <- stats::lm(y~., data=tdat)
   res.m <- stats::residuals(f1) 
   res.med <- cbind(res.med, res.m)
   }
   res.med <- apply(res.med, 2,scale,center=T, scale=T)   
   colnames(res.med) <- colnames(M)
   
   y.train <- res.y[train]
   y.test <- res.y[test]
   Med.train <- res.med[train,]
   Med.test <- res.med[test,]
   expo.train <- res.x[train]  
   expo.test <- res.x[test]
   } else {
   y.train <- Y[train]
   y.test <- Y[test]
   Med.train <- M[train,]
   Med.test <- M[test,]
   expo.train <- X[train]  
   expo.test <- X[test]
   }

   # optional filtering step by FDR on alpha
   cal_alpha_simple<-function(x){
   data2 <- data.frame(Med=x,envir=expo.train)
   l<-summary(stats::lm('Med ~.', data=data2))
   invalphap<-(l$coefficients)['envir','Pr(>|t|)']
   return(invalphap)
   }
   # keep the top candidates for ultra high dimensional data
   if(filtering & method !='ALL') {
   inv.alpha.p<-apply(Med.train,2, cal_alpha_simple)
   alpha.fdr <- stats::p.adjust(inv.alpha.p, method='fdr')
   order <- which(alpha.fdr < init.FDR.cutoff) 
   } else order <- 1:ncol(Med.train)

   #med <- med[,order, drop=F]
   Med.train <- Med.train[,order, drop=F]
   Med.test <- Med.test[,order, drop=F]

   #regress Y on X
   tdat <- data.frame(y=y.train, x=expo.train) 
   f0 <- stats::lm(y~., data=tdat)
   res.y.train <- stats::residuals(f0)  
   
   #regress M on X 
   res.med.train <- NULL
   for (i in 1:ncol(Med.train)){
   m <- Med.train[,i]
   tdat <- data.frame(y=m, x=expo.train)
   f1 <- stats::lm(y~., data=tdat)
   res.m.train <- stats::residuals(f1) 
   res.med.train <- cbind(res.med.train, res.m.train)
   }
   res.med.train <- apply(res.med.train, 2,scale,center=T, scale=T)   
  
   #perform variable selection
   if (method == 'iSIS'){
   model1 <- invisible(SIS::SIS(x=res.med.train, y=res.y.train, family='gaussian',tune='bic',  seed=1234, penalty='MCP',nsis=nsis, iter.max=iter.max))
   pab <- length(model1$ix)
   select <- (model1$ix)
   } 
   
   #all subjects and mediators are included
   if (method == 'ALL'){
   select <- 1:ncol(M)
   pab <- ncol(M)
   Med.test <- M
   y.test <- Y
   expo.test <- X
   }

   #in the testing data  
   sy<-stats::sd(y.test,na.rm=T)
 
   #Y~X
   testing <- data.frame(pheno=y.test,  envir=expo.test)
   f3<-tryCatch(stats::lm(pheno~envir, data=testing), error=function(c) NA, warning=function(c) NA)
   if (all(is.na(f3))) total<-Rsq.YX<-NA else {
   Rsq.YX=summary(f3)$r.squared #******
   total<-f3$coefficients[2]
   names(total) <-NULL }

   if (pab==0) {
   output <- c(Rsq.mediated=0,SOS=0, pab=0, Rsq.YM=0, Rsq.YX=Rsq.YX, Rsq.YMX=Rsq.YX, total=total, abproduct=0, n.train=length(y.train), n.estimate=length(y.test))
   paste0("There is no mediators selected.")
   return(list(output=output, select=NULL))
   } else {

   if (method=='iSIS') Med <- Med.test[,select,drop=F] else Med <- Med.test
   Med <- apply(Med, 2,scale,center=T, scale=T)

   #Y~M+X
   Kins <- as.matrix(Med)%*%t(as.matrix(Med))
   expo.std = scale(expo.test,center=T, scale=T)
   testing <- data.frame(res=y.test, expo=expo.std)
   colnames(Kins) <- rownames(Kins) <- 1:nrow(testing)
  
   testing$id <- 1:nrow(testing)
   fit11 <- tryCatch(GMMAT::glmmkin(res~expo, id='id', data=testing, kins=Kins, family=stats::gaussian(link='identity')), error=function(c) NA)
   if (all(is.na(fit11))) { tau<-phi<-direct<-NA } else { if (fit11$converged) {tau<-fit11$theta[2] ; phi<-fit11$theta[1];  direct<-fit11$coefficients[2]; } else tau<-phi<-direct<-NA }
   Rsq.YMX <- 1-phi/sy^2

   #Y~M
   fit22 <- tryCatch(GMMAT::glmmkin(res~1, data=testing, kins=Kins,id='id', family=stats::gaussian(link='identity')), error=function(c) NA)
   if (all(is.na(fit22))) { tau2<-phi2<-NA } else {
   if (fit22$converged) {tau2<-fit22$theta[2] ; phi2=fit22$theta[1];} else tau2<-phi2<-NA }
   Rsq.YM <- 1-phi2/sy^2

   Rsq.mediated <- Rsq.YM+Rsq.YX-Rsq.YMX
   SOS <-  Rsq.mediated/Rsq.YX
   abproduct= as.numeric(total-direct)
   output <- c(Rsq.mediated=Rsq.mediated,  SOS=SOS, pab=pab, Rsq.YM=Rsq.YM, Rsq.YX=Rsq.YX, Rsq.YMX=Rsq.YMX, total=total,  abprod=abproduct, n.train=length(y.train), n.estimate=length(y.test) )

   select <- colnames(Med)

   #select is the selected mediators, use med[,select] to look at the selected mediators (for future pathway analysis)
   return(list(output=output, select=select))
}
}
