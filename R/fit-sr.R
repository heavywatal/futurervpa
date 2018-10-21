#' Fit SR
#'
#' @description
#' 西嶋加筆 # 2018/06/07
#' 加入の残差の自己相関を考慮した再生産関係の推定
#' L1ノルム（最小絶対値）も推定できる (sigmaはSD)
#' TMB = TRUEでmarginal likelihood (.cppファイルが必要)
#' @rdname fit-sr
#' @export
fit.SR <- function(SRdata,SR="HS",method="L2",AR=1,TMB=FALSE,hessian=FALSE,w=rep(1,length(SRdata$year)),length=20){

  argname <- ls()
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname

  rec <- SRdata$R
  ssb <- SRdata$SSB

  N <- length(rec)

  #  if (SR=="HS") SRF <- function(x,a,b) a*(x+sqrt(b^2+gamma^2/4)-sqrt((x-b)^2+gamma^2/4))
  if (SR=="HS") SRF <- function(x,a,b) ifelse(x>b,b*a,x*a)
  if (SR=="BH") SRF <- function(x,a,b) a*x/(1+b*x)
  if (SR=="RI") SRF <- function(x,a,b) a*x*exp(-b*x)

  obj.f <- function(a,b,rho){
    resid <- sapply(1:N,function(i) log(rec[i]) - log(SRF(ssb[i],a,b)))
    resid2 <- NULL
    for (i in 1:N) {
      resid2[i] <- ifelse(i==1,resid[i], resid[i]-rho*resid2[i-1])
    }

    if (method == "L2") {
      sd <- sqrt(sum(resid2^2)/(N-rho^2))
      sd2 <- c(sd/sqrt(1-rho^2), rep(sd,N-1))
      obj <- -sum(w*dnorm(resid2,0,sd2,log=TRUE))
    } else {
      sd <- sum(abs(resid2))/(N-rho^2)
      sd2 <- c(sd/sqrt(1-rho^2), rep(sd,N-1))
      obj <- -sum(w*sapply(1:N, function(i){-log(2*sd2[i])-abs(resid2[i]/sd2[i])}))
    }
    return(obj)
  }

  a.range <- range(rec/ssb)
  b.range <- range(1/ssb)
  if (SR == "HS") b.range <- range(ssb)
  grids <- as.matrix(expand.grid(
    seq(a.range[1],a.range[2],len=length),
    seq(b.range[1],b.range[2],len=length)
  ))
  init <- as.numeric(grids[which.min(sapply(1:nrow(grids),function(i) obj.f(grids[i,1],grids[i,2],0))),])
  init[1] <- log(init[1])
  init[2] <- ifelse (SR == "HS",-log(max(0.000001,(max(ssb)-min(ssb))/max(init[2]-min(ssb),0.000001)-1)),log(init[2]))
  if (AR != 0) init[3] <- 0

  if (SR == "HS") {
    if (AR == 0) {
      obj.f2 <- function(x) obj.f(exp(x[1]),min(ssb)+(max(ssb)-min(ssb))/(1+exp(-x[2])),0)
    } else {
      obj.f2 <-  function(x) obj.f(exp(x[1]),min(ssb)+(max(ssb)-min(ssb))/(1+exp(-x[2])),1/(1+exp(-x[3])))
    }
  } else {
    if (AR == 0) {
      obj.f2 <- function(x) obj.f(exp(x[1]),exp(x[2]),0)
    } else {
      obj.f2 <-  function(x) obj.f(exp(x[1]),exp(x[2]),1/(1+exp(-x[3])))
    }
  }

  opt <- optim(init,obj.f2)
  opt <- optim(opt$par,obj.f2,method="BFGS",hessian=hessian)

  Res <- list()
  Res$input <- arglist
  Res$opt <- opt

  a <- exp(opt$par[1])
  b <- ifelse(SR=="HS",min(ssb)+(max(ssb)-min(ssb))/(1+exp(-opt$par[2])),exp(opt$par[2]))
  rho <- ifelse(AR==0,0,1/(1+exp(-opt$par[3])))
  resid <- sapply(1:N,function(i) log(rec[i]) - log(SRF(ssb[i],a,b)))
  resid2 <- NULL
  for (i in 1:N) {
    resid2[i] <- ifelse(i == 1,resid[i], resid[i]-rho*resid2[i-1])
  }
  sd <- ifelse(method=="L2",sqrt(sum(resid2^2)/(N-rho^2)),sqrt(2)*sum(abs(resid2))/(N-rho^2))

  Res$resid <- resid
  Res$resid2 <- resid2

  Res$pars <- c(a,b,sd,rho)

  if (method!="L2") {
    if (AR!=0) {
      arres <- ar(resid,aic=FALSE,order.max=1)
      Res$pars[3] <- sqrt(arres$var.pred)
      Res$pars[4] <- arres$ar
    }
  }

  Res$loglik <- loglik <- -opt$value

  if (method=="L2") {
    if (TMB) {
      data <- list()
      data$rec <- rec
      data$ssb <- ssb
      if (SR=="HS") data$SR <- 0
      if (SR=="BH") data$SR <- 1
      if (SR=="RI") data$SR <- 2
      #      data$gamma <- gamma

      params <- list()
      params$rec_loga <- opt$par[1]
      params$rec_logb <- ifelse(SR=="HS",-log(Res$pars[2]),opt$par[2])
      params$log_sd <- log(Res$pars[3]/(sqrt(1-Res$pars[4]^2)))
      params$logit_rho <- ifelse(AR==0,-20,opt$par[3])

      map <- list()
      if (AR==0) map$logit_rho<-factor(NA)
      obj <- MakeADFun(data, params, map=map,DLL="autoregressiveSR2",silent=TRUE)
      lower <- obj$par*0-Inf
      upper <- obj$par*0+Inf

      if (SR == "HS") {
        lower["rec_logb"] <- -log(max(ssb))
        upper["rec_logb"] <- -log(min(ssb))
      }
      opt <- nlminb(obj$par, obj$fn, obj$gr, lower=lower, upper=upper)
      rep <- sdreport(obj)

      # grid search
      if (SR != "HS") {
        grids <- expand.grid(seq(opt$par[1]-2,opt$par[1]+2,length=5),
                             seq(opt$par[2]-2,opt$par[2]+2,length=5))
      } else {
        grids <- expand.grid(seq(opt$par[1]-2,opt$par[1]+2,length=5),
                             seq(-log(max(ssb)),-log(min(ssb)),length=5))
      }
      params2 <- params
      params2$log_sd <- opt$par[3]
      if (AR == 1) params2$logit_rho <- opt$par[4]
      for (j in 1:nrow(grids)) {
        params2$rec_loga <- grids[j,1]
        params2$rec_logb <- grids[j,2]
        obj2 <- MakeADFun(data, params2, map=map, DLL="autoregressiveSR2",silent=TRUE)
        opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr, lower=lower, upper=upper)
        if (opt2$objective < opt$objective) {
          opt <- opt2
          obj <- obj2
          rep <- sdreport(obj2)
        }
      }

      Res$opt <- opt
      Res$rep <- rep
      if (SR=="HS") {
        Res$pars <- c(exp(rep$par.fixed[1]),1/exp(rep$par.fixed[2]),exp(rep$par.fixed[3]),ifelse(AR==0,0,1/(1+exp(-rep$par.fixed[4]))))
      } else {
        Res$pars <- c(exp(rep$par.fixed[1]),exp(rep$par.fixed[2]),exp(rep$par.fixed[3]),ifelse(AR==0,0,1/(1+exp(-rep$par.fixed[4]))))
      }
      Res$pars[3] <- sqrt(1-Res$pars[4]^2)*Res$pars[3]
      Res$loglik <- loglik <- -opt$objective

      a <- Res$pars[1]
      b <- Res$pars[2]
      rho <- Res$pars[4]
      resid <- sapply(1:N,function(i) log(rec[i]) - log(SRF(ssb[i],a,b)))
      resid2 <- NULL
      for (i in 1:N) {
        resid2[i] <- ifelse(i == 1,resid[i], resid[i]-rho*resid2[i-1])
      }
      Res$resid <- as.numeric(resid)
      Res$resid2 <- as.numeric(resid2)
    }
  }
  names(Res$pars) <- c("a","b","sd","rho")
  Res$pars <- data.frame(t(Res$pars))
  #  Res$gamma <- gamma

  ssb.tmp <- seq(from=0,to=max(ssb)*1.3,length=100)
  R.tmp <- sapply(1:length(ssb.tmp), function(i) SRF(ssb.tmp[i],a,b))
  pred.data <- data.frame(SSB=ssb.tmp,R=R.tmp)
  Res$pred <- pred.data

  Res$k <- k <- sum(Res$pars>0)
  Res$AIC <- -2*loglik+2*k
  Res$AICc <- Res$AIC+2*k*(k+1)/(N-k-1)
  Res$BIC <- -2*loglik+k*log(N)
  return(Res)
}
# Hockey-stick
