##
## ABC Calculation
##

calc.beta <- function(res,mY=5,prob.beta=c(0.5,0.9),prob.delta=c(0.9,0.95),beta=1,delta=1,beta.est=TRUE,delta.est=FALSE,beta.range=c(0,1),delta.range=c(0.1,5),Fm2.max=5,thin=1,step1=0.2,tol=0.0001,
                      Btar=res$Btar, # いちおう、各種管理基準値は外からでも与えられるようにした
                      Blow=res$Blow,
                      Blim=res$Blim,
                      Bban=res$Bban,
                      Fmsy=res$Fmsy)
{
  stockid <- res$stockid
  farg <- res$farg
  N <- res$N
  GT <- res$GT
  sim0 <- res$sim0
  nyear <- res$nyear
  lag <- res$lag
  future.vpa1 <- get(res$future.function.name)
  fm <- res$fm
  w.recent <- res$w.recent

  nY <- nyear
  eyear <- res$eyear
  SRdata <- res$SRdata

  seed <- res$seed

  B.cur <- SRdata$SSB[length(SRdata$SSB)]

  # alpha & abc calculation

  target.func <- res$target.func

  targ <- res$TARres

  Nlast <- targ$naa
  error.last <- targ$eaa+w.recent

  Dim.targ <- dim(Nlast)

  if (lag == 0) Blast <- NULL else Blast <- targ$ssb[,Dim.targ[2]-(lag-1),]

  targ.calc <- function(x) colSums(target.func(Fmsy,farg,mY=mY,eyear=1,seed=seed,N=N,naa0=Nlast[,Dim.targ[2],],eaa0=error.last[Dim.targ[2],],ssb0=Blast,beta=x,delta=0,Blim=Blim,Bban=Bban)$ssb[,1,])

  if (beta.est){
    beta.f1 <- function(x) {
        ssb.tmp <- targ.calc(x)
        Prob1 <- mean(ssb.tmp > Btar)
        x1 <- Prob1-prob.beta[1]
        dist1 <- x1^2
        dist1
    }

    beta.f2 <- function(x) {
        ssb.tmp <- targ.calc(x)
        Prob2 <- mean(ssb.tmp > Blim)

        x2 <- Prob2-prob.beta[2]
        dist1 <- x2^2
        dist1
    }

    res.beta1 <- optimize(beta.f1,beta.range)
    res.beta2 <- optimize(beta.f2,beta.range)

    beta <- min(res.beta1$minimum,res.beta2$minimum)
#    beta <- floor(beta * 100)/100
  }

  ## グラフによる図示
#  ssb.msy <- apply(targ$ssb,c(2,3),sum)[5,]
#  plot(density(ssb.msy),type="l",title="SSB")
#  abline(v=Blim,lty=2,col=2)
#  abline(v=Btar,lty=2)
#  mean(ssb.msy>Blim)
#  mean(ssb.msy>Btar)

  ssb.tmp <- targ.calc(beta)
  Prob.b1 <- mean(ssb.tmp > Btar)
  Prob.b2 <- mean(ssb.tmp > Blim)


  ## Limit

  if (B.cur >= Blim) lim <- res$LIMres else {
      N0 <- sim0$naa
      e0 <- sim0$eaa
      if (lag==0) SSB0 <- NULL else SSB0 <- sim0$ssb[,nY-(lag-1),]

      GT2 <- round(GT*2)

      Fm2 <- seq(1,Fm2.max,by=step1)

      FSYm <- lapply(Fm2, function(x) target.func(x*Fmsy,farg,seed=seed,mY=GT2,N=round(N/thin),eyear=1,naa0=N0[,nY,1:round(N/thin)],eaa0=e0[nY,1:round(N/thin)],ssb0=SSB0[,1:round(N/thin)]))

      FSYmest.s <- sapply(1:length(Fm2), function(x) mean(colSums(FSYm[[x]]$ssb[,1,])))

      num.l <- which.min((FSYmest.s - B.cur)^2)

      F.l <- Fm2[num.l]

      lim.calc <- function(x) target.func(x*Fmsy,farg,mY=GT2,seed=seed,N=N,eyear=1,naa0=N0[,nY,],eaa0=e0[nY,],ssb0=SSB0)

      obj.l <- function(x) (mean(colSums((lim.calc(x))$ssb[,1,]))-B.cur)^2
      res.obj.l <- optimize(obj.l, pmin(pmax(c(F.l-step1,F.l+step1),1),max(Fm2)),tol=tol)

      F.l <- res.obj.l$minimum

      lim <- target.func(F.l*Fmsy,farg,mY=GT2,seed=seed,N=N,eyear=eyear,naa0=N0[,nY,],eaa0=e0[nY,],ssb0=SSB0)
    }

    Nlast <- lim$naa
    error.last <- lim$eaa+w.recent

    Dim.lim <- dim(Nlast)

    if (lag == 0) Blast <- NULL else Blast <- lim$ssb[,Dim.lim[2]-(lag-1),]

    ssb0 <- colSums(target.func(0,farg,mY=mY,seed=seed,N=N,eyear=1,naa0=Nlast[,Dim.lim[2],],eaa0=error.last[Dim.lim[2],],ssb0=Blast,beta=1,delta=0,Blim=Blim,Bban=Bban)$ssb[,1,])


    Prob01 <- mean(ssb0 > Blim)
    Prob02 <- mean(ssb0 > Bban)

    # delta calculation

    lim.calc <- function(x) colSums(target.func(Fmsy,farg,mY=mY,seed=seed,N=N,eyear=1,naa0=Nlast[,Dim.lim[2],],eaa0=error.last[Dim.lim[2],],ssb0=Blast,beta=beta,delta=x,Blim=Blim,Bban=Bban)$ssb[,1,])

  if (delta.est){

    delta.f <- function(x) {
      Prob1 <- mean(lim.calc(x) > Blim)
      Prob2 <- mean(lim.calc(x) > Bban)

      dist1 <- (Prob1-Prob01*prob.delta[1])^2+(Prob2-Prob02*prob.delta[2])^2
      dist1
    }

    res.delta <- optimize(delta.f,delta.range,tol=tol)

    delta <- res.delta$minimum
  }

  Prob.d1 <- mean(lim.calc(delta) > Blim)
  Prob.d2 <- mean(lim.calc(delta) > Bban)

  farg <- res$farg

  farg$beta <- beta
  farg$delta <- delta
  farg$multi <- Fmsy
  farg$Blim <- Blim
  farg$Bban <- Bban
  farg$N <- N
  farg$nyear <- mY+1
  farg$naa0 <- NULL
  farg$eaa0 <- NULL
  farg$ssb0 <- NULL
#  if(!is.null(farg$vpares)) farg$vpares$Fc.at.age <- farg$vpares$Fc.at.age/fm
#  if(!is.null(farg$res0)) farg$res0$Fc.at.age <- farg$res0$Fc.at.age/fm

  future.pred <- do.call(future.vpa1,farg)
  ABC <- mean(future.pred$vwcaa[3,])  # calculation of ABC
  Fabc <- mean(colMeans(future.pred$faa[,3,]))

  future.B <- mean(colSums(future.pred$ssb[,mY,]))
  future.C <- mean(future.pred$vwcaa[mY,])

  # output

  Ccur <- res$Ccur
  Blim.cur <- res$Blim.cur

  Bref <- c(Btar, Blim, Bban)
  names(Bref) <- c("Target","Limit","Ban")
#  Fc.at.age <- farg$res0$Fc.at.age
#  {if(!is.null(farg$vpares)){
#       Fc.at.age <- farg$vpares$Fc.at.age
#   }
#n   else{
#
#   }
#  }

  cat("beta=",round(beta,2),"\n")

  out <- list(stockid=stockid, Bref=Bref, beta=beta, delta=delta, w.recent=w.recent,P0=c(Prob01, Prob02), P.beta=c(Prob.b1, Prob.b2), P.delta=c(Prob.d1,Prob.d2), Fmsy=Fmsy, B.cur=B.cur, ssb.tmp=ssb.tmp,Ccur=Ccur, Level=c(B.cur/Btar,B.cur/Blim), Blim.ratio=(Blim/10^6)/Blim.cur, future.B=future.B, future.C=future.C, Fabc=Fabc, abc=ABC, abc.ratio=ABC/Ccur)

  invisible(list(out,future.pred))
}
