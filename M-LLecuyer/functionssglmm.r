##########################################################################################
## functions to fit GLGM #################################################################
##########################################################################################

##########################################################################################
## Link functions ########################################################################
##########################################################################################

## Inverse by logit
inv.logit <- function(x){1/(1+exp(-x))}

## Inverse by log
inv.log <- function(x){exp(x)}


##########################################################################################
## Função para simular dos modelos #######################################################
##########################################################################################

rglgm <-function(n.sample, family, beta, extra.prec = NA, X , cov.pars, nugget, kappa, cov.model, ntrial = 1, offset = 1){
  s <- grf(n=n.sample, cov.pars=cov.pars, cov.model=cov.model, kappa = kappa, nugget = nugget,messages=FALSE)
  Xbeta <- X%*%beta
  if(family == "binomial"){
    p <- inv.logit(Xbeta + s$data)
    y <- rbinom(n.sample, size = ntrial, prob = p)}
  if(family == "poisson"){
    lambda <- offset*inv.log(Xbeta + s$data)
    y <- rpois(n.sample, lambda = lambda)}
  if(family == "negative.binomial"){
    mu <- offset*inv.log(Xbeta + s$data)
    y <- rnbinom(n.sample,mu=mu,size=extra.prec)}
  if(family == "gamma"){
    mu <- inv.log(Xbeta + s$data)
    y <- rgamma(n.sample, shape = extra.prec, scale = mu/extra.prec)}
  if(family == "beta"){
    mu <- inv.logit(Xbeta + s$data)
    y <- rbeta(n.sample, mu*extra.prec, (1-mu)*extra.prec)}
  saida <- data.frame("y"=y, "coord.X" = s$coords[,1], "coord.Y" = s$coords[,2], "efeito" = s$data)
  return(saida)
}

##########################################################################################
## Function to building the covariance matrix ############################################
##########################################################################################

monta.sigma <- function(cov.pars,cov.model,nugget,kappa,mat.dist){
  Sigma <- varcov.spatial(dists.lowertri = mat.dist, cov.model = cov.model, kappa=kappa,
                          nugget=nugget, cov.pars=cov.pars)
  return(Sigma)
}

##########################################################################################
## Function to evaluate the gaussian multivariate distribution ###########################
##########################################################################################

my.gauss <- function(b, det.Sigma, inv.Sigma){
 n <- length(b)
 dens <- (-n/2)*log(2*pi) - 0.5*det.Sigma - 0.5*t(b)%*%inv.Sigma%*%b
 return(dens)
}

##########################################################################################
## Integrand at log-scale Q(b) ###########################################################
##########################################################################################

## Case Poisson
Q.b.poisson <- function(b,Xbeta,Y,det.Sigma,inv.Sigma){
  eta = Xbeta + b
  dens = sum(dpois(Y,lambda = exp(eta),log=TRUE)) + 
         my.gauss(b,det.Sigma = det.Sigma, inv.Sigma = inv.Sigma)
  return(dens)
}

## Case Binomial
Q.b.binomial <- function(b,Xbeta,Y,det.Sigma,inv.Sigma,ntrial){
  eta = Xbeta + b
  p <- inv.logit(eta)
  dens = sum(dbinom(Y,size = ntrial, prob = p,log=TRUE)) + my.gauss(b,det.Sigma = det.Sigma, inv.Sigma = inv.Sigma)
  return(dens)
}

## Case Negative binomial
Q.b.negbin <- function(b,Xbeta,prec,Y,det.Sigma,inv.Sigma){
  eta = Xbeta + b
  dens = sum(dnbinom(Y, size=prec, mu = exp(eta), log=TRUE)) + my.gauss(b,det.Sigma = det.Sigma, inv.Sigma = inv.Sigma)
  return(as.numeric(dens))
}

## Case Gamma
Q.b.gamma <- function(b,Xbeta,prec,Y,det.Sigma,inv.Sigma){
  eta = Xbeta + b
  mu = exp(eta)
  dens = sum(dgamma(Y, shape=prec, scale = mu/prec,log=TRUE)) + my.gauss(b,det.Sigma = det.Sigma, inv.Sigma = inv.Sigma)
  return(as.numeric(dens))
}

## Case Beta
Q.b.beta <- function(b,Xbeta,prec,Y,det.Sigma,inv.Sigma){
  eta = Xbeta + b
  mu = inv.logit(eta)
  dens = sum(dbeta(Y, mu*prec, (1-mu)*prec, log=TRUE)) + my.gauss(b,det.Sigma = det.Sigma, inv.Sigma = inv.Sigma)
  return(as.numeric(dens))}


##########################################################################################
## First derived the integrand - Q'(b)  ##################################################
##########################################################################################

## Case Poisson
Q1.b.poisson <- function(b,Xbeta,Y,det.Sigma, inv.Sigma){
  grad <- t((Y - exp(Xbeta + b))) - t(b)%*%inv.Sigma
  return(grad)
}

## Case Binomial
Q1.b.binomial <- function(b,Xbeta,Y,det.Sigma, inv.Sigma,ntrial){
  b1 <- ntrial*inv.logit(Xbeta +b)
  grad <- t(Y - b1) - t(b)%*%inv.Sigma
  return(grad)
}

## Case Binomial negativa
Q1.b.negbin <- function(b,Xbeta, prec,Y,det.Sigma, inv.Sigma){
  et1 = exp(Xbeta)
  et2 = exp(b)
  et12 = et1*et2
  grad <- t((Y - (prec+Y)*(et12*(et12+prec)^-1))) - t(b)%*%inv.Sigma
  return(grad)
}

## Case Gamma
Q1.b.gamma <- function(b,Xbeta, prec,Y,det.Sigma, inv.Sigma){
  p1 = -prec + prec*(exp(-Xbeta -b))*Y 
  grad <- t(p1) - t(b)%*%inv.Sigma
  return(grad)
}

## Case Beta
Q1.b.beta <- function(b,Xbeta, prec,Y,det.Sigma, inv.Sigma){
  eta <- Xbeta + b
  mu <- inv.logit(eta)
  D.mu.b <- exp(-eta)/ ((1+exp(-eta))^2)
  dg1 <- digamma((1-mu)*prec)
  dg2 <- digamma(mu*prec)
  logY <- log(Y/(1-Y))
  part1 <- D.mu.b*prec*(dg1 - dg2 + logY)
  grad <- t(part1) - t(b)%*%inv.Sigma
  return(grad)
}

##########################################################################################
## Second derived the integrand ##########################################################
##########################################################################################

## Case Poisson
Q2.b.poisson <- function(b, Xbeta,Y,det.Sigma,inv.Sigma){
  eta = exp(Xbeta + b)
  diag(inv.Sigma) = eta + diag(inv.Sigma)
  Hess <- -inv.Sigma
  return(Hess)
}

## Case Binomial
Q2.b.binomial <- function(b, Xbeta,Y,det.Sigma,inv.Sigma,ntrial){
  b1 <- inv.logit(Xbeta + b)
  b2 <- exp(2*Xbeta + 2*b)/ ((1+exp(Xbeta + b))^2 )
  D <- b1-b2
  diag(inv.Sigma) = ntrial*D + diag(inv.Sigma)
  Hess <- -inv.Sigma
  return(Hess)
}

## Case Binomial negativa
Q2.b.negbin <- function(b, Xbeta,prec,Y,det.Sigma,inv.Sigma,ntrial){
  et1 = exp(Xbeta)
  et2 = exp(b)
  et12 = et1*et2
  p1 <- et12*((et12+prec)^-1)
  p2 <- p1^2
  D <- p1-p2
  diag(inv.Sigma) = (prec+Y)*D + diag(inv.Sigma)
  Hess <- -inv.Sigma
  return(Hess)
}

## Caso Gama
Q2.b.gamma <- function(b, Xbeta,prec,Y,det.Sigma,inv.Sigma,ntrial){
  p2 <- prec*Y*exp(-Xbeta - b)
  diag(inv.Sigma) = p2 + diag(inv.Sigma)
  Hess <- -inv.Sigma
  return(Hess)
}

## Caso Beta simplificada
Q2.b.beta <- function(b, Xbeta,prec, Y, det.Sigma,inv.Sigma){
  eta <- Xbeta + b
  mu <- inv.logit(eta)
  p1 <- mu*prec
  p2 <- (1-mu)*prec
  med <- mu*(1-mu)
  logY <- log(Y/(1-Y))
  d2 <- (1-mu)^2 - mu^2
  part1 <- -prec^2*(trigamma(p1) + trigamma(p2))*med
  part2 <- prec*(digamma(p2) - digamma(p1) + logY)*d2
  part3 <- (part1 + part2)*med
  diag(inv.Sigma) <- -part3 + diag(inv.Sigma)
  Hess <- - inv.Sigma
  return(Hess)
}

##########################################################################################
## Newton-Raphson generic ################################################################
##########################################################################################
newton_raphson <- function(initial, escore, hessiano, tol=0.0001, max.iter, n.dim,...){
  solucao <- matrix(NA, max.iter,n.dim)
  solucao[1,] <- initial
  for(i in 2:max.iter){
    HSS <- hessiano(initial, ...)
    ESC <- t(escore(initial, ...))
    solucao[i,] <- initial - solve(HSS,ESC)
    initial <- solucao[i,]
    tolera <- abs(solucao[i,] - solucao[i-1,])
    if(all(tolera < tol) == TRUE)break
  }
  saida <- list()
  saida[1][[1]] <- HSS
  saida[2][[1]] <- initial
  saida <<- initial
  return(saida)
	}

#newton_raphson <- function(initial, escore, hessiano, tol=0.0001, max.iter, n.dim,...){
#  solucao <- matrix(NA, max.iter,n.dim)
#  solucao[1,] <- initial
#print("Aqui")
#  for(i in 2:max.iter){
#    #print(i)
#    HSS <- hessiano(initial, ...)
#    ESC <- t(escore(initial, ...))
#    s.barra <- solve(HSS, -ESC)
#    alpha <- 0.5
#    h <- 10^-6
#    Ce <- (hessiano(initial+h*s.barra,...) - HSS)*(h^-1)
#    s2.barra <- solve(HSS+alpha*Ce, -0.5*Ce%*%s.barra)
#    #solucao[i,] <- initial - solve(HSS,ESC)
#    #print(c(initial+s.barra+s2.barra))
#    solucao[i,] <- initial + s.barra + s2.barra
#    initial <- solucao[i,]
#    tolera <- abs(solucao[i,] - solucao[i-1,])
#    if(all(tolera < tol) == TRUE)break
#  }
#  saida <- list()
#  saida[1][[1]] <- HSS
#  saida[2][[1]] <- initial
#  saida <<- initial
#  return(saida)
#}

##########################################################################################
## Numerical integration by Laplace method multidimensional ##############################
##########################################################################################

laplace <- function(Q.b, gr, hess, metodo, otimizador, n.dim, ...){
  log.integral <- -sqrt(.Machine$double.xmax)
  inicial <- rep(0,n.dim)
  if(metodo == "BFGS"){
  temp <- try(optim(inicial,Q.b, gr = gr, ...,
                    method=otimizador,hessian=TRUE,control=list(fnscale=-1)), silent=TRUE)
  }
  if(metodo == "NR"){
    temp <- try(newton_raphson(initial = inicial, escore = gr, hessiano = hess, n.dim = n.dim, max.iter = 100, ...), silent=TRUE)}
  if(metodo == "QNR"){
    temp <- try(qq.newton_raphson(initial = inicial, escore = gr, hessiano = hess, n.dim = n.dim, max.iter = 100, ...), silent=TRUE)}
  
  if(class(temp) != "try-error" & metodo == "BFGS"){
     log.integral <- temp$value + ((n.dim/2)*log(2*pi) - 0.5*determinant(-temp$hessian)$modulus)}

  if(class(temp) != "try-error" & metodo == "NR"){
    value <- Q.b(b = temp[2][[1]], ...)
    log.integral <- value + ((n.dim/2)*log(2*pi) - 0.5*determinant(-temp[1][[1]])$modulus)}
  .preditos <<- new.env(parent = .GlobalEnv)
  if(metodo == "NR"){
    assign("pred", value = temp[2][[1]], envir = .preditos)}
  if(metodo == "BFGS"){
    assign("pred", value = temp$par, envir = .preditos)}
  return(log.integral)
}

##########################################################################################
## Marginal likelihood ###################################################################
##########################################################################################

nlikGLGM <- function(par, Y, X, kappa, nugget, mat.dist, cov.model, family,metodo,ntrial=1, offset=NA){
  #print(par)
  I = -sqrt(.Machine$double.xmax)
  n <- length(Y)
  n.beta <- dim(X)[2]
  beta <- as.numeric(par[1:n.beta])
  Xbeta <- X%*%beta
  if(is.na(offset)[1] != TRUE){Xbeta <- cbind(X,log(offset))%*%c(beta,1)}
  sigma <- exp(as.numeric(par[c(n.beta+1)]))
  phi <- exp(as.numeric(par[c(n.beta+2)]))
  if(nugget == TRUE){tau2 <- exp(as.numeric(par[c(n.beta+3)]))}
  if(nugget == FALSE){tau2 <- 0}
  if(family == "negative.binomial"   & nugget == TRUE){prec <- exp(as.numeric(par[c(n.beta+4)]))}
  if(family == "gamma" & nugget == TRUE){prec <- exp(as.numeric(par[c(n.beta+4)]))}
  if(family == "beta" & nugget == TRUE){prec <- exp(as.numeric(par[c(n.beta+4)]))}
                                                     
  if(family == "negative.binomial" & nugget == FALSE){prec <- exp(as.numeric(par[c(n.beta+3)]))}
  if(family == "gamma" & nugget == FALSE){prec <- exp(as.numeric(par[c(n.beta+3)]))}
  if(family == "beta" & nugget == FALSE){prec <- exp(as.numeric(par[c(n.beta+3)]))}
                                                                                        
  if(kappa != "NULL"){kappa=exp(as.numeric(kappa))}
  #print(kappa)
  # Montando a matriz de covariancia
  Sigma <- as.matrix(forceSymmetric(monta.sigma(cov.pars = c(sigma,phi), 
                    cov.model=cov.model, nugget = tau2, kappa = kappa, mat.dist = mat.dist)$varcov))
  chol.Sigma <- try(chol(Sigma),silent=TRUE)
  det.Sigma <- try(sum(log(diag(chol.Sigma)))*2,silent=TRUE)
  inv.Sigma <- try(chol2inv(chol.Sigma), silent=TRUE)
  if(class(chol.Sigma) != "try-error"){
    if(class(inv.Sigma)[1] != "try-error"){
  #### Avaliando a verossimilhança
  ## Caso Poisson
    if(family == "poisson"){
      I <- laplace(Q.b.poisson, gr = Q1.b.poisson, hess = Q2.b.poisson, metodo= metodo, otimizador="BFGS", 
                   n.dim = n, Xbeta = Xbeta, Y = Y, det.Sigma = det.Sigma, inv.Sigma = inv.Sigma)}
  ## Caso Binomial
    if(family == "binomial"){
      I <- laplace(Q.b.binomial, gr = Q1.b.binomial, hess = Q2.b.binomial,
                   metodo= metodo, otimizador="BFGS", n.dim = n, Xbeta = Xbeta, Y = Y, det.Sigma = det.Sigma, 
                   inv.Sigma = inv.Sigma,ntrial=ntrial)}
  if(family == "negative.binomial"){
    I <- laplace(Q.b.negbin, gr = Q1.b.negbin, hess = Q2.b.negbin, metodo= metodo, otimizador="BFGS", n.dim = n,
                   Xbeta = Xbeta, prec = prec,Y = Y, det.Sigma = det.Sigma, inv.Sigma = inv.Sigma)}
   if(family == "gamma"){
    I <- laplace(Q.b.gamma, gr = Q1.b.gamma, hess = Q2.b.gamma, metodo= metodo, otimizador="BFGS", n.dim = n,
                   Xbeta = Xbeta, prec = prec,Y = Y, det.Sigma = det.Sigma, inv.Sigma = inv.Sigma)}
  if(family == "beta"){
    I <- laplace(Q.b.beta, gr = Q1.b.beta, hess = Q2.b.beta, metodo= metodo, otimizador="BFGS", n.dim = n,
                   Xbeta = Xbeta, prec = prec,Y = Y, det.Sigma = det.Sigma, inv.Sigma = inv.Sigma)}
  }
  }
 return(-I)
}

##########################################################################################
## Fit function  #########################################################################
##########################################################################################

glgm <- function(formula, cov.model, kappa, inits, data, coords, nugget, family, ntrial=1, offset=1, method.optim, method.integrate, predict=TRUE){
  formula <- as.formula(formula)
  mf <- model.frame(formula,data=data)
  Y <- model.response(mf)
  X <- model.matrix(formula, data=data)
  mat.dist <- dist(coords)
  names <- c(colnames(X),"logsigma2", "logphi")
  n.beta <- dim(X)[2]
  
  if(family == "negative.binomial" & nugget == TRUE){names <- c(names,"logtau2","logprec")}
  if(family == "gamma" & nugget == TRUE){names <- c(names,"logtau2","logprec")}
  if(family == "beta" & nugget == TRUE){names <- c(names,"logtau2","logprec")}

  if(family == "negative.binomial"   & nugget == FALSE){names <- c(names,"logprec")}
  if(family == "gamma" & nugget == FALSE){names <- c(names,"logprec")}
  if(family == "beta" & nugget == FALSE){names <- c(names,"logprec")}
  
  if(nugget == TRUE & family == "poisson"){names <- c(names,"logtau2")}
  if(nugget == TRUE & family == "binomial"){names <- c(names,"logtau2")}
   names(inits) <- parnames(nlikGLGM) <- names
  estimativas <- mle2(nlikGLGM, start=inits,
                     vecpar=TRUE,
                     method=method.optim,
                     control=list(maxit=1000),
                     skip.hessian=FALSE,
                     data=list(Y = Y, X = X, mat.dist = mat.dist,cov.model= cov.model,ntrial=ntrial, nugget=nugget,
                     family=family,metodo=method.integrate,kappa=kappa,offset=offset))
  n.pars <- length(coef(estimativas))
  summary.estimativas <- summary(estimativas)
  summary.estimativas@coef[,1][c(n.beta+1):n.pars] <- exp(summary.estimativas@coef[,1][c(n.beta+1):n.pars])
  std.error = sqrt(exp(summary.estimativas@coef[,1][c(n.beta+1):n.pars])^2*(summary.estimativas@coef[,2][c(n.beta+1):n.pars]^2))
  summary.estimativas@coef[,2] <- c(summary.estimativas@coef[,2][1:n.beta],std.error)
  summary.estimativas@coef[,3] <- summary.estimativas@coef[,1]/summary.estimativas@coef[,2]
  summary.estimativas@coef[,4] <- NA
  #if(nugget == FALSE){nomes <- c("sigma2","phi")}
  #if(nugget == TRUE){nomes <- c("sigma2","phi","tau2")}
  #rownames(summary.estimativas@coef) <- c(colnames(X),nomes)
  if(predict == TRUE) {
  preditos <- .preditos$pred
  saida <- list()
  saida[1][[1]] <- summary.estimativas
  saida[7][[1]] <- logLik(estimativas)
  saida[2][[1]] <- preditos
  saida[3][[1]] <- coords
  saida[4][[1]] <- cov.model
  saida[5][[1]] <- family
  saida[6][[1]] <- exp(coef(estimativas)[c(n.beta+1):n.pars])
  saida[8][[1]] <- ifelse(cov.model == "matern", exp(kappa), "NULL")
  saida[9][[1]] <- estimativas
  return(saida)}
  if(predict == FALSE){return(estimativas)}
}

##########################################################################################
## Function to obtain start values #######################################################
##########################################################################################

start.values.glgm <- function(formula, data, coords, nugget, family, ntrial=1, offset=1){
  mf <- model.frame(formula,data)
  Y <- model.response(mf)
  X <- model.matrix(formula ,data=data)

  if( family == "binomial"){
    response <- cbind(Y,ntrial -Y)
    fit <- glm(response ~ -1 + X, data=data, family="binomial")
    print(logLik(fit))
    esp <- predict(fit, type="response")
    #sigma = var(Y/ntrial - esp)
    res1 <<- Y/ntrial - esp 
    sigma = sd(Y/ntrial - esp)
    phi <- 0.1*max(dist(coords))
    saida <- c(coef(fit),log(sigma),log(phi))
    names(saida) <- c(colnames(X),"logsigma","logphi") # changed colnames(X) to names(coef(fit))
    if(nugget == TRUE){ nugget = 0.1*sigma
                        saida <- c(saida,"logtau" = log(nugget))}}
  if(family == "poisson"){
    fit <- glm(Y ~ -1 + X, family= "poisson", data=data, offset=log(offset))
    print(logLik(fit))
    esp <- predict(fit)
    sigma <- var( esp - log(Y+1))
    phi <- 0.1*max(dist(coords))
    saida <- c(coef(fit),log(sigma),log(phi))
    names(saida) <- c(colnames(X),"logsigma","logphi")
    if(nugget == TRUE){ nugget = 0.1*sigma
                        saida <- c(saida,"logtau" = log(nugget))}}

   if(family == "negative.binomial"){
    fit <- glm.nb(Y ~ -1 + X + offset(log(offset)), data=data)
    print(logLik(fit))
    esp <- predict(fit)
    sigma <- var( esp - log(Y+1))
    phi <- 0.1*max(dist(coords))
    saida <- c(coef(fit),log(sigma),log(phi))
    names(saida) <- c(colnames(X),"logsigma","logphi")
    if(nugget == TRUE){ nugget = 0.1*sigma
                        saida <- c(saida,"logtau" = log(nugget))}
    saida <- c(saida, "logtheta" = log(fit$theta))
  }
 if(family == "gamma"){
    fit <- glm(Y ~ -1 + X, family=Gamma(link="log"),data=data)
    te <- summary(fit)
    print(logLik(fit))
    esp <- predict(fit)
    sigma <- var(esp - log(Y))
    phi <- 0.1*max(dist(coords))
    saida <- c(coef(fit),log(sigma),log(phi))
    names(saida) <- c(colnames(X),"logsigma","logphi")
    if(nugget == TRUE){ nugget = 0.1*sigma
                        saida <- c(saida,"logtau" = log(nugget))}
    saida <- c(saida, "logtheta" = te$dispersion)}
  if(family == "beta"){
    fit <- betareg(Y ~ -1 + X,data=data)
    te <- summary(fit)
    n.beta <- dim(X)[2]
    print(logLik(fit))
    esp <- predict(fit)
    sigma <- var(esp - Y)
    phi <- 0.1*max(dist(coords))
    saida <- c(coef(fit)[1:n.beta],log(sigma),log(phi))
    names(saida) <- c(colnames(X),"logsigma","logphi")
    if(nugget == TRUE){ nugget = 0.1*sigma
                        saida <- c(saida,"logtau" = log(nugget))}
    saida <- c(saida, "logtheta" = log(as.numeric(coef(fit)[n.beta+1])))}


  
    return(saida)
}

##########################################################################################
## summary.glgm ##########################################################################
##########################################################################################

summary.glgm <- function(obj.glgm){
  obj.glgm[[1]]}

logLik.glgm <- function(obj.glgm){
  obj.glgm[[7]]}

##########################################################################################
## Function to predict random effect #####################################################
##########################################################################################

prediction <- function(model,locations){
  names(locations) <- c("coord.X", "coord.Y")
  n.locations <- dim(locations)[1]
  coords <- model[[3]]
  names(coords) <- c("coord.X", "coord.Y")
  gridi <- rbind(locations,coords)
  pars <- model[6][[1]]
  if(length(pars) == 2){pars <- c(pars,0)}
  mat.dist <- dist(gridi)
  n.total <- dim(gridi)[1]
  Sigma <- monta.sigma(cov.pars = c(pars[1:2]),cov.model = model[4][[1]],
                       nugget = pars[3], kappa = model[8][[1]], mat.dist = mat.dist)$varcov
  S11 <- Sigma[1:n.locations,1:n.locations]
  S22 <- Sigma[c(n.locations+1):n.total,c(n.locations+1):n.total]
  S12 <- Sigma[1:n.locations,c(n.locations+1):n.total]
  S21 <- Sigma[c(n.locations+1):n.total,1:n.locations]
  preditos <- S12%*%solve(S22,model[2][[1]])
  #var.preditos <- Sigma11 - Sigma12%*%solve(Sigma22)%*%Sigma21
  return(preditos)
}
  
  

