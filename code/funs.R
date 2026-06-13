# rm(list=ls())
library(mvnfast)
library(mccr)
library(glmnet)


##### Functions #####

### l2 error 
l2 = function(hat,est){
  return( sqrt( sum( (hat-est)^2 ) ) )
}

### summary data
mysummary = function(hat, real){
  p=length(hat)
  ressum = rep(0,5)
  names(ressum) = c( "l_2", "l_sigma", "l_inf", "mcc", "FDRFNR" )
  delta =  hat - real
  
  ressum[1] = sqrt( sum( delta^2 ) )
  ressum[2] = sqrt( t(delta) %*% toeplitz( 0.5^(0:(p-1)) ) %*% delta )
  ressum[3] = max( abs(delta) )
  
  estloc = ifelse(hat==0, 0, 1)
  trueloc = ifelse(real==0, 0, 1)
  ressum[4] = mccr( estloc, trueloc )
  
  TP <- sum( hat &  real) # True Positive
  FN <- sum(!hat &  real) # False Negative
  TN <- sum(!hat & !real) # True Negative
  FP <- sum( hat & !real) # False Positive
  ressum[5] = FP/( max(1, TP+FP ) ) + FN/ ( TP+FN )
  
  return(ressum)
}


### ADIHT
# X             current batch X
# Y             current batch Y
# g_prime       g'
# g_second      g''
# warm          start with beta_ini (usually the last batch output), T; or 0_p, F
# beta_ini      initialization beta_0
# Inter Hess    summary stats
# N             cumulative sample size (up to the last batch)
# lambda_0      initial threshold
# lambda_inf    stat thresholds
# kappa         decay rate of lambda
# eta           learn rate
# roundcoef     learning round, learn additional roundcoef*log(Nnew) round
# 
ADIHT.online = function(X, Y, g_prime, g_second, warm = FALSE, 
                        beta_ini = NULL, Inter = NULL, Hess = NULL, N = 0,
                        lambda_0 = 5, lambda_inf = 0.5, etaloc= FALSE,
                        kappa = 0.95, eta=0.5, method="IHT", roundcoef = 20) {
  
  # Batch size and dimension
  n = dim(X)[1]; p = dim(X)[2]
  
  # Initialization
  if (is.null(Inter)) {  Inter = rep(0, p) }
  if (is.null(Hess)) {  Hess = matrix(0, p, p) }
  if (warm == F | is.null(beta_ini) ) {  beta_ini = rep(0,p) }
  
  # Cumulative sample
  Nnew = N + n
  D = 1e3; eta = eta/0.9
  
  # Decrease learn rate eta if the iteration not converge
  while(D >= 1e3 | is.na(D) ){
    
    # Initialization
    eta = eta*0.9
    t=0; betat = beta_ini; lamt = lambda_0
    # Do dynamic IHT
    while( t <= log(lambda_0/lambda_inf)/log(1/kappa) + roundcoef*log(Nnew) ){
      tempb = betat
      
      gradt = t(X) %*% ( g_prime(X%*%betat) - Y )
      H = betat - eta/Nnew * ( gradt + Inter + Hess %*% betat)
      
      if(method=="IHT"){
        lamt = max(kappa*lamt, lambda_inf)
        betat = ifelse( abs(H) >= lamt*( eta*(etaloc==TRUE) +  (etaloc== FALSE) )  , H, 0 )
      }
      if(method=="lasso"){
        if(etaloc ==F){
          betat = sign(H) * ifelse( abs(H) >= lambda_inf, 
                                    abs(H) -  lambda_inf, 0 )
        }else{
          betat = sign(H) * ifelse( abs(H) >= eta*lambda_inf, 
                                    abs(H) - eta* lambda_inf, 0 )
        }
      }
      t=t+1
      
      D = sqrt( sum( betat^2 ) )
      if(D >= 1e3 | is.na(D) ){break}
      
      delt = sum( (tempb - betat)^2 ) 
      
      if( (delt<= 1e-4) & (t > log(lambda_0/lambda_inf)/log(1/kappa) + 2*log(Nnew) ) ){
        break}
      
    }
    
  }
  
  gradb = t(X) %*% ( g_prime(X%*%betat) - Y )
  Hessb = t(X) %*% diag( c(g_second(X%*%betat)) ) %*% X
  
  #Update 
  Internew <- Inter + gradb - (Hessb %*% betat)
  Hessnew <- Hess + Hessb 
  
  return(list(beta_old = beta_ini,  beta_new = betat,
              #Inter_old= Inter,     
              Inter_new= Internew,
              #Hess_old = Hess,  
              Hess_new = Hessnew, 
              N_old    = N,         N_new    = Nnew,
              eta      = eta,       laminf   = lambda_inf ))
}



### Renewable from Luo EJS (implement by PGD)
# X             current batch X
# Y             current batch Y
# g_prime       g'
# g_second      g''
# warm          start with beta_ini (usually the last batch output), T; or 0_p, F
# beta_ini      initialization beta_0
# Hess          summary Hessian matrix
# N             cumulative sample size (up to the last batch)
# lambda_inf    regularization  
# eta           learn rate
# method        get lasso/scad estimation
#
Renew.online = function(X, Y, g_prime, g_second, warm = FALSE, 
                        beta_ini = NULL, Hess = NULL, N = 0, etaloc =T,
                        lambda_0 = 5, kappa = 0.95, roundcoef = 20,
                        lambda_inf = 0.5, eta=0.5, method = "lasso" ) {

  # Batch size and dimension
  n = dim(X)[1]; p = dim(X)[2]
  
  # Initialization
  if (is.null(Hess)) {  Hess <- matrix(0, p, p) }
  
  # Take beta_0 as the initialization
  if (warm == F | is.null(beta_ini) ) {
    beta_0 = rep(0,p) 
  }else{
    beta_0 = beta_ini
  }
  
  # Cumulative sample
  Nnew <- N + n
  D = 1e3;  eta = eta/0.9;  
  
  # Decrease learn rate eta if the iteration not converge
  while(D >= 1e3 | is.na(D) ){
    
    # Initialization
    eta = eta*0.9;  delta = 1
    t=0;  betat = beta_0;  lamt = lambda_0
    # Do PGD
    while( t <= log(lambda_0/lambda_inf)/log(1/kappa) + roundcoef *log(Nnew) ){
      temp = betat
      gradt = t(X) %*% ( g_prime(X%*%betat) - Y )
      H = betat - eta/Nnew * ( gradt + Hess %*% (betat-beta_ini) )
      if(method == "lasso"){
        if(etaloc ==F){
          betat = sign(H) * ifelse( abs(H) >= lambda_inf, 
                                    abs(H) -  lambda_inf, 0 )
        }else{
          betat = sign(H) * ifelse( abs(H) >= eta*lambda_inf, 
                                    abs(H) - eta* lambda_inf, 0 )
        }
      }
      if(method=="IHT"){
        lamt = max(kappa*lamt, lambda_inf)
        betat = ifelse( abs(H) >= lamt, H, 0 )
      }
      if(method == "scad"){
        betat = H
        for(j in 1:p){
          if( abs(H[j]) <= 2*lambda_inf){
            betat[j] = sign(H[j])* max( 0, abs(H[j]) - lambda_inf )
          }else if( abs(H[j]) <= 3.7*lambda_inf ) {
            betat[j] = (2.7*H[j] - 3.7*sign(H[j])*lambda_inf  )/1.7
          }
        }
      }
      
      t=t+1
      delta = sqrt( sum( (temp - betat)^2 ) )
      betat[1:15]
    }
    
    D = sqrt( sum( betat^2 ) )
  }
  
  #Update 
  Hessb = t(X) %*% diag( c(g_second(X%*%betat)) ) %*% X
  Hessnew <- Hess + Hessb 
  
  return(list(beta_old = beta_ini,  beta_new = betat,
              # Hess_old = Hess,      
              Hess_new = Hessnew, 
              N_old    = N,         N_new    = Nnew,
              eta      = eta,       laminf   = lambda_inf ))
}


### Renewable from Huang JMLR (implement by PGD)
# X             current batch X
# Y             current batch Y
# g_prime       g'
# g_second      g''
# warm          start with beta_ini (usually the last batch output), T; or 0_p, F
# beta_ini1     last batch output 1, used in the second part
# Hess1         summary Hessian matrix, first part data, second part output
# beta_ini2     last batch output 2, used in the first part
# Hess2         summary Hessian matrix, second part data, first part output
# N             cumulative sample size (up to the last batch)
# lambda_inf    regularization  
# eta           learn rate
# method        get lasso/scad estimation
#
OSIM = function(X, Y, g_prime, g_second, warm = FALSE, 
                beta_ini1 = NULL, Hess1 = NULL,
                beta_ini2 = NULL, Hess2 = NULL, etaloc = T,
                N = 0, lambda_inf = 0.5, roundcoef = 20, 
                eta=0.5, method = "lasso"){
  
  # Batch size and dimension
  n = dim(X)[1];        p = dim(X)[2]
  part1 = 1:(n/2);  part2 = (n/2+1):n
  X1 = X[part1, ];     Y1 = Y[part1]
  X2 = X[part2, ];     Y2 = Y[part2]
  
  # Initialization
  if (is.null(Hess1)) {  Hess1 <- matrix(0, p, p) }
  if (is.null(Hess2)) {  Hess2 <- matrix(0, p, p) }
  
  # Take beta_01,2 as the initialization
  if (warm == F | is.null(beta_ini1) ) {
    beta_01 = beta_02 = rep(0,p);  
  }else{
    beta_01 = beta_ini1; beta_02 = beta_ini2
  }
  
  # Cumulative sample
  Nnew <- N + n
  D = 1e3; eta = eta/0.9
  
  # Decrease learn rate eta if the iteration not converge
  while(D >= 1e3 | is.na(D) ){
    
    # Initialization
    eta = eta*0.9
    
    # Do ISTA1
    t1=0; betat1 = beta_01; delta1 = 1
    while( t1 <= roundcoef *log(Nnew) & (delta1 > 1e-3) ){
      temp1 = betat1
      gradt1 = t(X1) %*% ( g_prime(X1%*%betat1) - Y1 )
      H1 = betat1 - eta/Nnew * ( 2*gradt1 + Hess1 %*% (betat1-beta_ini2) )
      if(method == "lasso"){
        if(etaloc==F){
          betat1 = sign(H1) * ifelse( abs(H1) >= lambda_inf, abs(H1)- lambda_inf, 0 )
        }else{
          betat1 = sign(H1) * ifelse( abs(H1) >=eta*lambda_inf, abs(H1)-eta*lambda_inf, 0 )
        }
      }
      if(method == "scad"){
        betat1 = H1
        for(j in 1:p){
          if( abs(H1[j]) <= 2*lambda_inf){
            betat1[j] = sign(H1[j])* max( 0, abs(H1[j]) - lambda_inf )
          }else if( abs(H1[j]) <= 3.7*lambda_inf ) {
            betat1[j] = (2.7*H1[j] - 3.7*sign(H1[j])*lambda_inf  )/1.7
          }
        }
      }
      
      t1=t1+1 
      delta1 = sqrt( sum( (temp1 - betat1)^2 ) )
      if( is.na(delta1) ){delta1=0}  #if NA: stop iteration
    }
    
    
    # Do ISTA2
    t2=0; betat2 = beta_02; delta2 = 1
    while( t2 <= roundcoef*log(Nnew) & (delta2 > 1e-3) ){
      temp2 = betat2
      gradt2 = t(X2) %*% ( g_prime(X2%*%betat2) - Y2 )
      H2 = betat2 - eta/Nnew * ( 2*gradt2 + Hess2 %*% (betat2-beta_ini1) )
      if(method == "lasso"){
        if(etaloc==F){
          betat2 = sign(H2) * ifelse( abs(H2) >= lambda_inf, abs(H2)- lambda_inf, 0 )
        }else{
          betat2 = sign(H2) * ifelse( abs(H2) >=eta*lambda_inf, abs(H2)-eta*lambda_inf, 0 )
        }
        
      }
      if(method == "scad"){
        betat2 = H2
        for(j in 1:p){
          if( abs(H2[j]) <= 2*lambda_inf){
            betat2[j] = sign(H2[j])* max( 0, abs(H2[j]) - lambda_inf )
          }else if( abs(H2[j]) <= 3.7*lambda_inf ) {
            betat2[j] = (2.7*H2[j] - 3.7*sign(H2[j])*lambda_inf  )/1.7
          }
        }
      }
      
      t2=t2+1 
      delta2 = sqrt( sum( (temp2 - betat2)^2 ) )
      if( is.na(delta2) ){delta2=0}  #if NA: stop iteration
    }
    
    D = sqrt( max( sum( betat1^2 ), sum( betat2^2 ) )  )
  }
  
  #Update 
  Hessb1 = 2* t(X1) %*% diag( c(g_second(X1%*%betat2)) ) %*% X1
  Hessb2 = 2* t(X2) %*% diag( c(g_second(X2%*%betat1)) ) %*% X2
  
  return(list(beta_old1 = beta_ini1,  beta_new1 = betat1,
              beta_old2 = beta_ini2,  beta_new2 = betat2,
              # Hess_old1 = Hess1,      
              Hess_new1 = Hess1 + Hessb1 , 
              # Hess_old2 = Hess2,      
              Hess_new2 = Hess2 + Hessb2 , 
              N_old    = N,         N_new    = Nnew,
              eta      = eta,       laminf   = lambda_inf,
              beta_ave = ( betat1 + betat2 ) / 2,
              iter = c(t1,t2) ) )
}


### RADAR (tuning guided by Han2026JASA, online GLM)
# X              current batch X
# Y              current batch Y
# g_fun          function g
# g_prime        g' 
# betapast       last batch (epoch) output
# lastlambdalist lambda list used in the last batch learning
# betapastlist   the auxiliary beta list used for check the first sample's lambda
# mylambdalist   the regularization parameter set in current batch learning
# Rk             radius
# alphak         learn rate
#
RADAR = function(X, Y, g_fun, g_prime, betapast,
                 lastlambdalist = 10^seq(-4, 0, l=9),
                 betapastlist = matrix(0, 1000, 9),
                 mylambdalist = 10^seq(-4, 0, l=9)/sqrt(2), 
                 Rk=6, alphak=10 ){
  
  # lambdalist = 10^seq(-4, 0, l=9); alphak=10; Rk=6
  
  Len = length(mylambdalist)
  n = length(Y)                
  p = dim(X)[2]
  
  # prepare [Len] paths with different lambda in mylambdalist
  mut = matrix(0, p, Len)      ;  bt = matrix(0, p, Len);  bt[,] = betapast
  
  # collect real-time estimates through different lambda in mylambdalist
  temp = matrix(0, p , Len) 
  
  # prepare optimal path via forward cross-validation
  muopt = rep(0,p)             ;  bopt = betapast
  
  # collect optimal path by adaptively tuning lambda
  tempopt = rep(0, p)  
  
  
  q1= 2*log(p)/ ( 2*log(p) -1 );  q2 = q1/ ( q1 -1 )
  
  for( t in 1:n){
    
    # get optimal regularization para lambda by using last estimate
    if( t ==1 ){
      # if t==1, use the last batch estimates in betapastlist
      xbpast = X[t,] %*% betapastlist
      lamidx = which.min( -Y[t]* t(xbpast) + g_fun( t(xbpast) ) )[[1]]
      
      lamopt = lastlambdalist[lamidx]
    }else{
      #if t>=2, use last last-time estimates in temp
      xbpast = X[t,] %*% temp
      lamidx = which.min( -Y[t]* t(xbpast) + g_fun( t(xbpast) ) )[[1]]
      
      lamopt = mylambdalist[lamidx]
    } 
    
    # optimal learning 
    muopt = muopt + c( g_prime( t(X[t,])%*%bopt ) - Y[t] )*X[t,] + 
      lamopt*sign(bopt)
    
    at = alphak/sqrt(t)
    muoptq2norm = ( sum( abs(muopt)^q2 ) )^{1/q2}
    xi = max( 0, (q1-1)*at*muoptq2norm*Rk -1 )
    
    bopt = betapast - Rk^2*(q1-1)*at/(1+xi)*muoptq2norm^(2-q2)* 
      sign(muopt)*( abs(muopt)^{q2-1} )
    tempopt = tempopt + bopt
    
    
    # learn one sample with different lambda in mylambdalist
    for( lamidx in 1: Len ){
      
      mut[, lamidx] = mut[, lamidx] + c( g_prime( t(X[t,])%*%bt[, lamidx] ) - Y[t] )*X[t,] +
        mylambdalist[lamidx] * sign(bt[, lamidx])
      
      at = alphak/sqrt(t)
      mutq2norm = ( sum( abs(mut[, lamidx])^q2 ) )^{1/q2}
      xi = max( 0, (q1-1)*at*mutq2norm*Rk -1 )
      
      bt[, lamidx] = betapast- Rk^2*(q1-1)*at/(1+xi)*mutq2norm^(2-q2)*
        sign(mut[, lamidx])*( abs(mut[, lamidx])^{q2-1} )
      
      # Polyak-Ruppert averaging
      temp[, lamidx] = temp[, lamidx]*(t-1)/t + bt[, lamidx]/t
    }
    
  }
  
  betaoptout = tempopt/n
  betaoptout = ifelse( abs(betaoptout)> 1e-2, betaoptout, 0) 
  
  return(list( betahat = betaoptout,
               usedRk = Rk,
               usedalphak = alphak,
               usedlambdalist = mylambdalist,
               usedbetaseq = temp ) )
}


