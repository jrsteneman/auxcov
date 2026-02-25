# CITE: Steneman, J. and Vinci, G. (2024). Covariance matrix completion via auxiliary information. arXiv:2402.05767




#######################################
# Function "positivize" (Algorithm 5) #
#######################################
# INPUT
#- A: pxp symmetric matrix with positive diagonals. This is the matrix to be positivized.
#- delta: positive real number that is added to the diagonals of the rescaled version of A 
#         at each iteration.
#- verbose: logical value. If TRUE, the function prints the total value added to the diagonals
#           required to achieve positivization.
#
# OUTPUT
#- A.tilde: positive definite version of matrix A.

positivize = function(A,delta=.001,verbose=FALSE){

  SD = sqrt(diag(A))
  B = cov2cor(A)
  count = 0
  
  while(min(eigen(B)$values) <= 0){
    diag(B) = diag(B) + delta
    count=count+1
  }
  
  if(count>0 & verbose==TRUE){
    print(paste0("positive def. correction (+",count*delta,")"))
  }
  
  A.tilde = diag(SD)%*%cov2cor(B)%*%diag(SD)
  
 return(A.tilde)  
}






#############################################
# Function "Psi.gauss" (Corollary 3.1(iii)) #
#############################################
# INPUT
#- X: nxp data matrix. It may contain NAs.
#- Sigma: pxp covariance matrix about X (either true matrix or estimate). It may contain NAs 
#         with pairwise missingness consistent with missingness in X.
#- positivize.delta: A real number. If positive, the computed Psi matrix is positivized via 
#                    the function "positivize" with delta=positivize.delta.
#
# OUTPUT
#- Psi: Gaussian based Psi matrix estimate.

Psi.gauss = function(X,Sigma,positivize.delta=0.001){ 
  
  Us.fun = function(X){
    p=ncol(X)
    Sigma.hat.O = cov(X,use='pairwise.complete.obs')
    U = cbind(rep(1:p,each=p),rep(1:p,p))
    U = U[!is.na(c(Sigma.hat.O)),]
    U.BAR = U[(U[,2]-U[,1])>=0,]
    U = U[(U[,2]-U[,1])>0,]
    return(list(U=U,U.BAR=U.BAR))
  }
  
  H.fun = function(X,Sigma,U.BAR){

    H = matrix(0,nc=nrow(U.BAR),nr=nrow(U.BAR))
    OBS = !is.na(X)
    
    for(s in 1:nrow(U.BAR)){
      i = U.BAR[s,1]       
      j = U.BAR[s,2]   
      vec.ij = OBS[,i]*OBS[,j]
      nij = sum(vec.ij)
      
      for(z in s:nrow(U.BAR)){
        k = U.BAR[z,1]              
        l = U.BAR[z,2]              
        vec.kl = OBS[,k]*OBS[,l]              
        nijkl = sum(vec.ij*vec.kl)              
        nkl = sum(vec.kl)
        
        if(nijkl>0){
          H[s,z] = H[z,s] = ((Sigma[i,k]*Sigma[j,l]) + (Sigma[i,l] * Sigma[j,k]))*nijkl/(nij*nkl)
        }

      }
            
    }
    
    H = nrow(X)*H
    
    return(H)
  }
  

  J.fun = function(X,Sigma,U,U.BAR){

    Jacobian = matrix(NA,nc=nrow(U.BAR),nr=nrow(U))
    
    for(s in 1:nrow(U)){
      i = U[s,1]            
      j = U[s,2]
      
      for(z in 1:nrow(U.BAR)){
        k = U.BAR[z,1]
        l = U.BAR[z,2]
        
        if(i == j){
          
          Jacobian[s,z] = 0
          
        }else if(i == k && j == l && i != j){
          
          Jacobian[s,z]  = 1/(sqrt(Sigma[i,i]*Sigma[j,j]))
          
        }else if(i == l && j == k && i != j){
          
          Jacobian[s,z] = 1/(sqrt(Sigma[i,i]*Sigma[j,j]))
          
        }else if(i == l && i == k && i != j){
          
          Jacobian[s,z] = -Sigma[i,j]/(2*(sqrt(Sigma[i,i]))^3 * sqrt(Sigma[j,j]))
          
        }else if(j == l && j == k && i != j){
          
          Jacobian[s,z] = -Sigma[i,j]/(2*(sqrt(Sigma[j,j]))^3 * sqrt(Sigma[i,i]))
          
        }else{
          
          Jacobian[s,z] = 0
          
        }        
      }
    }

    return(Jacobian)
  }
  

  F.fun = function(cors){
    return(diag(1/(1-(cors^2))))
  }
 
  
  US = Us.fun(X)
  H = H.fun(X=X,Sigma=Sigma,U.BAR=US$U.BAR)
  J = J.fun(X=X,Sigma=Sigma,U=US$U,U.BAR=US$U.BAR)
  R = cov2cor(Sigma)
  CORS = sapply(1:nrow(US$U),function(i) R[US$U[i,1],US$U[i,2]])
  F = F.fun(CORS)
  
  Psi = F%*%(J%*%H%*%t(J))%*%F

  if(positivize.delta>0){
    Psi = positivize(Psi,delta=positivize.delta)
  }

  return(Psi)

}





#########################################
# Function "Psi.emp" (Corollary 3.1(ii) #
#########################################
# INPUT
#- X: an nxp data matrix. It may contain NAs.
#- positivize.delta: A real number. If positive, the computed Psi matrix is positivized via the 
#                    function "positivize" with delta=positivize.delta.
#
# OUTPUT
#- Psi: Empirical Psi matrix.

Psi.emp = function(X,positivize.delta=0.001){ 
  
  
  Us.fun = function(X){
    p=ncol(X)
    Sigma.hat.O = cov(X,use='pairwise.complete.obs')
    U = cbind(rep(1:p,each=p),rep(1:p,p))
    U = U[!is.na(c(Sigma.hat.O)),]
    U.BAR = U[(U[,2]-U[,1])>=0,]
    U = U[(U[,2]-U[,1])>0,]
    return(list(U=U,U.BAR=U.BAR))
  }
  
  
  H.fun = function(X,U.BAR){
    
    H = matrix(0,nc=nrow(U.BAR),nr=nrow(U.BAR))
    OBS = !is.na(X)
    
    for(s in 1:nrow(U.BAR)){
      i = U.BAR[s,1]       
      j = U.BAR[s,2]   
      vec.ij = OBS[,i]*OBS[,j]
      nij = sum(vec.ij)
      
      for(z in s:nrow(U.BAR)){
        k = U.BAR[z,1]              
        l = U.BAR[z,2]              
        vec.kl = OBS[,k]*OBS[,l]              
        nijkl = sum(vec.ij*vec.kl)              
        nkl = sum(vec.kl)
        
        if(nijkl>0){

          Z.i = X[,i]-mean(X[,i],na.rm = TRUE)
          Z.j = X[,j]-mean(X[,j],na.rm = TRUE)
          Z.k = X[,k]-mean(X[,k],na.rm = TRUE)
          Z.l = X[,l]-mean(X[,l],na.rm = TRUE)

          Xijkl = mean(Z.i*Z.j*Z.k*Z.l,na.rm = TRUE)
          Xij = mean(Z.i*Z.j,na.rm = TRUE)
          Xkl = mean(Z.k*Z.l,na.rm = TRUE)
          
          H[s,z] = H[z,s] = (Xijkl - (Xij*Xkl))*nijkl/(nij*nkl)

        }
        
      }
      
    }
    
    H = nrow(X)*H
    
    return(H)
  }
  
  
  J.fun = function(X,Sigma,U,U.BAR){
    
    Jacobian = matrix(NA,nc=nrow(U.BAR),nr=nrow(U))
    
    for(s in 1:nrow(U)){
      i = U[s,1]            
      j = U[s,2]
      
      for(z in 1:nrow(U.BAR)){
        k = U.BAR[z,1]
        l = U.BAR[z,2]
        
        if(i == j){
          
          Jacobian[s,z] = 0
          
        }else if(i == k && j == l && i != j){
          
          Jacobian[s,z]  = 1/(sqrt(Sigma[i,i]*Sigma[j,j]))
          
        }else if(i == l && j == k && i != j){
          
          Jacobian[s,z] = 1/(sqrt(Sigma[i,i]*Sigma[j,j]))
          
        }else if(i == l && i == k && i != j){
          
          Jacobian[s,z] = -Sigma[i,j]/(2*(sqrt(Sigma[i,i]))^3 * sqrt(Sigma[j,j]))
          
        }else if(j == l && j == k && i != j){
          
          Jacobian[s,z] = -Sigma[i,j]/(2*(sqrt(Sigma[j,j]))^3 * sqrt(Sigma[i,i]))
          
        }else{
          
          Jacobian[s,z] = 0
          
        }        
      }
    }
    
    return(Jacobian)
  }
  
  
  F.fun = function(cors){
    return(diag(1/(1-(cors^2))))
  }
  
  Sigma.hat = cov(X,use='pairwise.complete.obs')
  R = cov2cor(Sigma.hat)
  US = Us.fun(X)
  H = H.fun(X=X,U.BAR=US$U.BAR)
  J = J.fun(X=X,Sigma=Sigma.hat,U=US$U,U.BAR=US$U.BAR)
  CORS = sapply(1:nrow(US$U),function(i) R[US$U[i,1],US$U[i,2]])
  F = F.fun(CORS)
  
  Psi = F%*%(J%*%H%*%t(J))%*%F
  
  if(positivize.delta>0){
    Psi = positivize(Psi,delta=positivize.delta)
  }

  return(Psi)
  
}



###################################
# Function "auxcov" (Algorithm 1) #
###################################
# INPUT
#- Sigma0: pxp observed sample covariance matrix estimate that may contain NAs.
#- X: nxp data matrix that may contain NAs. It is used to compute Sigma0 if not provided.
#- W: A list containing q pxp symmetric matrices of auxiliary variables. Only the upper triangle 
#     of each matrix is used.
#- alpha: A constant between 0 and 1 inclusive. Tuning parameter which determines shrinkage 
#         towards auxiliary baseline values.
#- regmod: A string. It specifies the type of baseline regression model used. It can be set equal
#          to 'ols', 'gls', or 'splines'.
#- intercept: A logical value. If TRUE, an intercept term is included in the regression.
#- Psi: A |U| x |U| symmetric covariance matrix of the measurement error. Usually the output of
#      Psi.gauss or Psi.emp. It is used only if regmod = 'gls'. If left NULL, Psi is computed 
#      automatically (see argument gls.Psi).
#- splines.knots: A non-negative integer specifying the number of knots used in spline regression.
#                 It is used only if regmod = 'splines'.
#- g: A function that maps (-1,1) to R. Default is the Fisher transformation.
#- g.inv: The inverse of g. Defaul is the inverse Fisher transformation.
#- positivize.delta: A real number. If positive, the function "positivize" with delta=positivize.delta
#                    is applied to various matrices involved in the algorithm.
#- gls.min.iters: A positive integer. The minimum number of iterations for the gradient ascent 
#                 algorithm used when regmod = 'gls'.
#- gls.max.iters: A positive integer. The maximum number of iterations for the gradient ascent 
#                 algorithm used when regmod = 'gls'.
#- b: A positive constant. The step size for the gradient ascent algorithm used when regmod = 'gls'
#- gls.stoppct: A positive constant. When regmod = 'gls', it determines the stop threshold for the
#               gradient ascent as a percent change threshold in the objective function.
#- gls.Psi: A string. When regmod = 'gls', will compute Psi based either on a Gaussian assumption
#           (gls.Psi = 'gaussian') or nonparametrically (gls.Psi = 'empirical').
#
# OUTPUT
#- Sigma.alpha: A symmetric pxp matrix. The completed covariance matrix, computed as a convex combination
#               of Cor.bar and Cor.tilde, weighted by alpha.
#- Cor.bar: A symmetric pxp matrix. The auxiliary baseline correlation matrix predicted by the auxiliary 
#           baseline regression.
#- Cor.tilde: A symmetric pxp matrix. The sample correlation matrix, with its missing values filled in 
#             with values predicted by the auxiliary baseline regression.
#- mod.ols: An lm object containing the fitted regression model if regmod = 'ols'.
#- mod.gls: A list containing the fitted regression coefficients ('Beta'), the estimated variance of 
#           the error term epsilon ('sigma2'), and the computed log likelihoods in the gradient ascent 
#           algorithm ('logliks') if regmod = 'gls'.
#- mod.splines: The fitted splines regression model if regmod = 'splines'. 



auxcov = function(SigmaO = NULL, X=NULL, W, alpha = 0.5, regmod = 'ols', intercept=TRUE, Psi=NULL, 
          splines.knots = 0, g=atanh, g.inv=tanh, positivize.delta=0.001,gls.min.iters=30,
          gls.max.iters=500,gls.b=0.001,gls.stoppct = 1e-7, gls.Psi='gaussian'){
  

  # VARIABLES PREPARATION

  if(is.null(SigmaO)){
    SigmaO=cov(X,use = 'pairwise.complete.obs')
  }
  SD = diag(sqrt(diag(SigmaO)))
  UP = upper.tri(SigmaO)
  CorO = cov2cor(SigmaO)
  y = g(CorO[UP])
  W.mat = matrix(data = 0,nrow = length(y),ncol = length(W))
    for(i in 1:length(W)){
      W.mat[,i] = W[[i]][UP]
    }  
  YW= as.data.frame(cbind(y,W.mat))


  # REGRESSION

  mod.ols = NA
  mod.gls = NA
  mod.splines = NA

  if(regmod=='ols'){
    if(intercept==FALSE){
      mod = lm(y ~. -1, data = YW)
    }
    if(intercept==TRUE){
      mod = lm(y ~., data = YW)
    }

    mod.ols=mod

    y.pred = predict(mod,newdata=YW)

  }
  

  if(regmod=='gls'){

    GLS = function(Y,X,Psi,min.iters=10,max.iters=100,b=.001,stoppct = 1e-7,b.acc=1.4,Plot=FALSE){

      X = as.matrix(X)
      if(min(eigen(Psi)$val)<0){
        print('Warning: Psi is not positive definite')
      }
      n = length(Y)

     loglik = function(Beta,logsigma2) {
        SIG = exp(logsigma2)*diag(n)+Psi
        IMAT = solve(SIG)
        ERR = Y-X%*%Beta
        A = as.numeric(determinant(SIG,logarithm=TRUE)$mod)
        B = sum(diag(t(ERR)%*%IMAT%*%ERR))
       -(A+B)/2
      }

      gradf = function(Beta,logsigma2,Y,X, Psi){
        THETA = solve(exp(logsigma2)*diag(n)+Psi)
        Z = Y-X%*%Beta
        dBeta = t(X)%*%THETA%*%Z
        dlogsigma2 = -(1/2)*exp(logsigma2)*sum(diag(THETA))+(1/2)*exp(logsigma2)*t(Z)%*%THETA%*%THETA%*%Z
        return(c(dlogsigma2,dBeta))
      }

      # start values
      LM.start = lm(Y~X-1)
      Beta = as.numeric(LM.start$coeff)
      logsigma2 = 2*log(summary(LM.start)$sigma)
      
      # iterations

            if(Plot==TRUE){
              dev.new()
            }

        logliks = numeric()
        loglik.opt = -1e99
        it = 1

        repeat{
          GRAD = gradf(Beta,logsigma2,Y,X,Psi)
          Beta.it = Beta+b*GRAD[-1]
          logsigma2.it = logsigma2+b*GRAD[1]
          loglik.it = loglik(Beta.it,logsigma2.it)

          if(it > max.iters){
            break
          }
          if(loglik.it > loglik.opt && abs((loglik.it - loglik.opt)/loglik.opt) < stoppct && it>min.iters){
            break
          }
          if(is.nan(loglik.it)){
            loglik.it = loglik.opt
          }
          if(loglik.it<loglik.opt){
            b=b/b.acc
          }
          if(loglik.it>loglik.opt){
            loglik.opt = loglik.it
            Beta=Beta.it
            logsigma2=logsigma2.it
            b = b*b.acc
          }
          logliks=c(logliks,loglik.opt)

          if(Plot==TRUE){
            plot(1:it,logliks,xlab='iterations',ylab='log-likelihood')
          }

          it = it + 1
        }

      return(list(Beta=Beta,sigma2=exp(logsigma2),logliks=logliks))
    }

    if(is.null(Psi)){
      if(gls.Psi=='gaussian'){
        Psi = Psi.gauss(X=X,Sigma = SigmaO,positivize.delta=positivize.delta)/nrow(X)
      }
      if(gls.Psi=='empirical'){
        Psi = Psi.emp(X=X,K.pairwise = FALSE,positivize.delta=positivize.delta)/nrow(X)
      }
    }
    
    W_gls = W.mat

    if(intercept==TRUE){
      W_gls = cbind(rep(1,nrow(W_gls)),W_gls)
    }

    GLS_out = GLS(Y=y[!is.na(y)],X=W_gls[!is.na(y),],Psi=Psi,min.iters=gls.min.iters,
                  max.iters=gls.max.iters,b=gls.b,stoppct=gls.stoppct)
    Beta = GLS_out$Beta

    mod.gls = GLS_out

    y.pred = as.vector(W_gls%*%Beta)

  }
  

  if(regmod=='splines'){
    
    df = 1 + splines.knots
    
    if(df == 1){
      
      splineframe = matrix(data = 0,nrow = length(y),ncol = length(W))
      splineframe = as.data.frame(splineframe)
      
      for(i in 1:length(W)){
        splineframe[,i] = splines::bs(W[[i]][UP])
      }

      mod = lm(y ~. ,data = splineframe)  
      
    }
    
    if(df > 1){
      
      splineframe = matrix(data = 0,nrow = length(y),ncol = length(W))
      splineframe = as.data.frame(splineframe)
      
      for(i in 1:length(W)){
        s = splines::ns(W[[i]][UP], df = df)
        splineframe[,i] = splines::bs(W[[i]][UP],knots = as.numeric(attr(s,"knots")))
      }
      
      mod = lm(y ~. ,data = splineframe)
    }
    
    mod.splines=mod
    y.pred = predict(mod,newdata=splineframe)

  }
  

  # auxiliary baseline "Cor.bar"
  Cor.bar = diag(ncol(CorO))*0
  Cor.bar[UP] = g.inv(y.pred)
  Cor.bar = Cor.bar + t(Cor.bar)
  diag(Cor.bar) = 1
  
  
  # complete matrix "Cor.tilde"
  Cor.tilde = CorO
  Cor.tilde[is.na(CorO)] = Cor.bar[is.na(CorO)]
  
  
  if(positivize.delta>0){
    Cor.bar = positivize(Cor.bar,delta=positivize.delta)
    Cor.tilde = positivize(Cor.tilde,delta=positivize.delta)
  }
  
  
  # final covariance matrix "Sigma.alpha"
  if(length(alpha)==1){
    Cor.alpha = alpha*Cor.bar+(1-alpha)*Cor.tilde 
    Sigma.alpha = SD%*%Cor.alpha%*%SD
  }
  

  if(length(alpha)>1){  
    Sigma.alpha = lapply(alpha,function(aa){
      Cor.aa = aa*Cor.bar+(1-aa)*Cor.tilde     
      Sigma.aa = SD%*%Cor.aa%*%SD
     return(Sigma.aa) 
    })
  }
  
  return(list(Sigma.alpha=Sigma.alpha,Cor.bar=Cor.bar,Cor.tilde=Cor.tilde,mod.ols=mod.ols,
              mod.gls=mod.gls,mod.splines=mod.splines))
  
}





######################################
# Function "auxcov.cv" (Algorithm 2) #
######################################
# INPUT
#- X: nxp data matrix that may contain NAs.
#- W: A list containing q pxp symmetric matrices of auxiliary variables. Only the upper triangle 
#     of each matrix is used.
#- alpha: A vector of constants between 0 and 1 inclusive (see function "auxcov"). Cross-Validation 
#         will select one of these values.
#- folds: Number of folds to be used in cross-validation.
#- cores: Number of CPU cores for parallel computing on compatible machines. Setting cores = 1 
#         disables this feature.
#- regmod: A string. It specifies the type of baseline regression model used. It can be set equal
#          to 'ols', 'gls', or 'splines'.
#- splines.knots.seq: Vector of nonnegative integers specifying the number of knots in splines regression.
#                     Cross-Validation will select one of these values if regmod = 'splines'.
#- Plot.cv: A logical value. If TRUE, it plots cross validation risk versus alpha (and spline knots, if
#           applicable)
#- positivize.delta: A real number. If positive, the function "positivize" with delta=positivize.delta
#                    is applied to various matrices involved in the algorithm.
#- ... : Additional arguments valid for the auxcov function. See auxcov function.
#
# OUTPUT
#- RISK: An vector containing the CV risk corresponding to the vector alpha, or a matrix containing the 
#        CV risk corresponding to the vector alpha and vector splines.knots.seq if regmod='splines'.
#- alpha.cv: Optimal alpha selected via CV.
#- splines.knots.cv: Optimal number of splines knots selected via CV (if regmod = 'splines').
#- Sigma.alpha.cv: Complete pxp covariance matrix computed via auxcov using the alpha.cv and splines.knots.cv.
#- mod.ols: An lm object containing the fitted regression model if regmod = 'ols'.
#- mod.gls: A list containing the fitted regression coefficients ('Beta'), the estimated variance of 
#           the error term epsilon ('sigma2'), and the computed log likelihoods in the gradient ascent 
#           algorithm ('logliks') if regmod = 'gls'.
#- mod.splines: The fitted splines regression model if regmod = 'splines', with number of knots splines.knots.cv.



auxcov.cv = function(X,W,alpha,folds=10,cores=1,regmod = 'ols',splines.knots.seq=0,Plot.cv=FALSE,positivize.delta=0.001,...){
  
  LAPPLY = lapply
  if(cores>1){
    LAPPLY = function(...) parallel::mclapply(...,mc.cores=cores)
  }
  

  loss.fun = function(Sigma1,Sigma2){
    O = !is.na(Sigma1+Sigma2)    
    Sigma1 = cov2cor(Sigma1)
    Sigma2 = cov2cor(Sigma2)
    UP = upper.tri(diag(ncol(Sigma1)))
    UP.O = UP & O
    sum((Sigma1[UP.O]-Sigma2[UP.O])^2)
  }
  
  folds.maker = function(X,N=10){
    
    block.identity = function(X){
      Vr = lapply(1:nrow(X),function(i)which(!is.na(X[i,])))
      Vs=unique(lapply(1:nrow(X),function(i)which(!is.na(X[i,]))))
      
      blocks = numeric()
      for(i in 1:length(Vs)){
        for(r in 1:length(Vr)){
          if(setequal(Vr[[r]],Vs[[i]])){
            blocks[r] = i
          }
          
        }
      }
      
      return(blocks)
    }
    
    BLOCKS = block.identity(X)
    SAMPLES = 1:nrow(X)
    
    for(i in 1:max(BLOCKS)){
      SAMPLES[BLOCKS==i] = sample(SAMPLES[BLOCKS==i])
    }
    
    FOLDS = list()
    
    for(k in 1:N){
      
      FOLDS[[k]] = c(unlist(lapply(1:max(BLOCKS),function(i){
        sam.i = SAMPLES[BLOCKS==i]
        n.i.N = round(length(sam.i)/N)
        low = ((k-1)*n.i.N)+1
        up = (k*n.i.N)
        if(k==N){
          up = sum(BLOCKS==i)
        }
        sel = low:up
        return(sam.i[sel])
      })))
      
    }
    
    if(var(sapply(1:N,function(i) length(FOLDS[[i]])))>0){
      print('FOLDS HAVE DIFFERENT SIZES')
    }

   return(FOLDS)
  }


  FOLDS = folds.maker(X, N = folds)


    if(regmod%in%c('ols','gls')){

      LOSS = matrix(NA,nr=folds,nc=length(alpha))

      losses.alpha = LAPPLY(1:folds,function(i){
        
        testing = X[FOLDS[[i]],]   
        training = X[-(FOLDS[[i]]),]
        
        Sigma_training = auxcov(SigmaO = cov(training, use='pairwise.complete.obs'), W = W, alpha = alpha, 
          regmod=regmod, X=training, positivize.delta=positivize.delta,...)$Sigma.alpha    
        Sigma_testing = cov(testing, use = "pairwise.complete.obs")
        
        return(sapply(1:length(Sigma_training),function(j) loss.fun(Sigma_training[[j]],Sigma_testing)))
        
      })

      for(i in 1:folds){
        LOSS[i,] = losses.alpha[[i]]
      }
      
      RISK = colMeans(LOSS)
      alpha.cv = alpha[which.min(RISK)]
      
      if(Plot.cv == TRUE){  
        plot(alpha,RISK,type='b',pch=16,xlab=expression(alpha))
        abline(v=alpha.cv)
      }
      
      OPT.cv = auxcov(SigmaO = cov(X, use='pairwise.complete.obs'), W=W, alpha=alpha.cv, regmod=regmod, X=X, 
        positivize.delta=positivize.delta,...)
      splines.knots.cv = FALSE

    }


    if(regmod=='splines'){

      LOSS = array(NA,dim=c(folds,length(alpha),length(splines.knots.seq)))
      losslist = list()
      
      for(k in 1:length(splines.knots.seq)){
        losses.alpha.k = LAPPLY(1:folds,function(i){

          testing = X[FOLDS[[i]],]
          training = X[-(FOLDS[[i]]),]
          
          Sigma_training = auxcov(cov(training, use='pairwise.complete.obs'),W = W,alpha = alpha,regmod=regmod, 
            X=training, splines.knots = splines.knots.seq[k],...)$Sigma.alpha
          Sigma_testing = cov(testing, use='pairwise.complete.obs')

          return(sapply(1:length(Sigma_training),function(j) loss.fun(Sigma_training[[j]],Sigma_testing)))

        })

        losslist[[k]] = losses.alpha.k
      }
      
      for(i in 1:folds){
        for(k in 1:length(splines.knots.seq)){
          LOSS[i,,k] = losslist[[k]][[i]]
        }
      }
      
      RISK = apply(LOSS,c(2,3),mean)
      
      min.coord=function(A){
        row.min = which.min(apply(A,1,min))
        col.min = which.min(apply(A,2,min))
       return(c(row.min,col.min))
      }

      opt.coord = min.coord(RISK)
      alpha.cv = alpha[opt.coord[1]]
      splines.knots.cv = splines.knots.seq[opt.coord[2]]
      
      OPT.cv = auxcov(SigmaO = cov(X, use='pairwise.complete.obs'), W=W, alpha=alpha.cv, regmod=regmod, X=X, 
        splines.knots=splines.knots.cv, positivize.delta=positivize.delta,...)

      if(Plot.cv == TRUE){
        par(mfrow=c(1,2))
        image(alpha,splines.knots.seq,RISK)
        points(alpha.cv,splines.knots.cv,pch=16, cex=3,col='green')
        persp(alpha,splines.knots.seq,RISK)
      }
      
  }

 return(list(RISK = RISK, alpha.cv = alpha.cv, splines.knots.cv=splines.knots.cv, Sigma.alpha.cv=OPT.cv$Sigma.alpha, 
            mod.ols=OPT.cv$mod.ols, mod.gls=OPT.cv$mod.gls, mod.splines=OPT.cv$mod.splines))
  
}







###########################################
# Function "auxcov.boot" (Algorithms 3-4) #
###########################################
# INPUT:
#- X: X: nxp data matrix. It may contain NAs.
#- W: A list containing q pxp symmetric matrices of auxiliary variables. Only the upper triangle 
#     of each matrix is used.
#- alpha: A vector of constants between 0 and 1 inclusive (see function "auxcov"). Cross-Validation 
#         will select one of these values.
#- folds: Number of folds to be used in cross-validation.
#- mod: A string. Bootstrap samples are drawn from a Normal distribution if mod = 'par.gauss', or 
#       from the empirical distribution if mod = 'nonpar'.
#- boots: Number of bootstrap samples to be generated.
#- cores: Number of CPU cores for parallel computing on compatible machines. Setting cores = 1 
#         disables this feature.
#- transf: function of covariance matrix estimator of which we want to approximate the standard deviation.
#- ... : Additional arguments valid for the auxcov or auxcov.cv functions. See auxcov or auxcov.cv functions.
#
# OUTPUT:
#- Sigma.boots: A list containing the pxp covariance matrix estimates computed by auxcov.cv for each boostrap sample.
#- SD: A pxp matrix. The estimated standard errors for the entries of the completed covariance matrix, computed from Sigma.boots.


auxcov.boot = function(X,W,alpha,folds=10,mod='nonpar',boots=200,cores=1,transf=identity,...){

  LAPPLY = lapply
    if(cores>1){
      LAPPLY = function(...) parallel::mclapply(...,mc.cores=cores)
    }


  SD.mat = function(SIGMAS,transf){ 

   N = length(SIGMAS)

   S1 = lapply(1:N,function(b) transf(SIGMAS[[b]]))
   S2 = lapply(1:N,function(b) transf(SIGMAS[[b]])^2)

   M1 = Reduce("+",S1)/N
   M2 = Reduce("+",S2)/N
   SD = sqrt(M2-M1^2)

   return(SD)
  }



  if(mod=='nonpar'){

      block.identity = function(X){

        Vr = lapply(1:nrow(X),function(i)which(!is.na(X[i,])))
        Vs = unique(lapply(1:nrow(X),function(i)which(!is.na(X[i,]))))
    
        blocks = numeric()
        for(i in 1:length(Vs)){
          for(r in 1:length(Vr)){
            if(setequal(Vr[[r]],Vs[[i]])){
              blocks[r] = i
            }       
         }
        } 
       return(blocks)
      }
  

      boot.maker = function(X,boots=100,cores=1){

        BLOCKS = block.identity(X)

        X.boots = LAPPLY(1:boots,function(b){
          X.boot = X
            for(i in 1:max(BLOCKS)){
             X.boot[BLOCKS==i,] = X[sample(which(BLOCKS==i),replace=TRUE),]
            }
           return(X.boot)
          })
       return(X.boots)
      }


   X.boots = boot.maker(X=X,boots=boots,cores=cores)
  
  }



  if(mod=='par.gauss'){

    SIGMA.est = auxcov.cv(X=X,W=W,alpha=alpha,folds=folds,cores=cores,Plot.cv=FALSE,...)$Sigma.alpha.cv

    MASK = X*0+1
    n=nrow(X); p=ncol(X)

    X.boots = LAPPLY(1:boots,function(b){      
      X.b = MASS::mvrnorm(n,mu=rep(0,p),Sigma=SIGMA.est)*MASK
       return(X.b)
    })

  }


  Sigma.boots = LAPPLY(1:boots,function(b){

      auxcov.cv(X=X.boots[[b]],W=W,alpha=alpha,folds=folds,cores=1,Plot.cv=FALSE,...)$Sigma.alpha.cv

    })

  SD = SD.mat(SIGMAS=Sigma.boots,transf=transf)


 return(list(Sigma.boots=Sigma.boots,SD=SD))

}











###########
# Example #
###########


# Simulation settings 
p = 30
n = 500
gamma = 0.8
vars = rep(1,p)
alpha = seq(0,1,by = 0.1)

set.seed(1)

# Generate auxiliary variable matrix
W = matrix(0,nc=p,nr=p)
UP = upper.tri(W)
W[UP] = runif(p*(p-1)/2,-1,1)
W = W+t(W)
diag(W)=1
W.mat = list(W)

# Generate error term matrix
E = matrix(0,nc=p,nr=p)
UP = upper.tri(E)
E[UP] = runif(p*(p-1)/2,-1,1)
E = E+t(E)
diag(E)=1

# Generate true Sigma as in Section 4, ensuring it is positive definite
Sigma = (sqrt(gamma)*W+sqrt(1-gamma)*E)/sqrt(2)
Sigma = positivize(Sigma,delta=0.001)
Sigma = diag(vars)%*%Sigma%*%diag(vars)

# Generate data matrix, adding missingness
X = MASS::mvrnorm(n,mu=rep(0,p),Sigma = Sigma)
X[1:(n/2),1:(p/3)] = NA
X[((n/2)+1):n,((2*p/3)+1):p] = NA

# Performing AUXCOV-CV
CV.out = auxcov.cv(X = X,W = W.mat,alpha = alpha,folds=10,cores=1,regmod = 'ols',splines.knots.seq=0,
        Plot.cv=TRUE,positivize.delta=0.001)

alpha.CV = CV.out$alpha.cv
Sigma.CV = CV.out$Sigma.alpha.cv

dev.new()
plot(cov2cor(Sigma)[UP],cov2cor(Sigma.CV)[UP],pch=16,col=rgb(0,0,0,.4),xlab='true correlation', 
      ylab='correlation estimate')
abline(a=0,b=1)


# Performing AUXCOV nonparametric Bootstrap
Boot.out = auxcov.boot(X = X,W = W.mat,alpha = alpha,folds=10,mod='nonpar',boots=200,cores=1,transf=identity)
SD = Boot.out$SD
