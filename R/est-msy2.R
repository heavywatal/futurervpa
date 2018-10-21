
est.MSY2 <- function(vpares,N=1000,res1=NULL,sim0=NULL,nyear=NULL,pgy=0.9,lim=0.6,ban=0.1,mY=5,long.term=20,
                     Fmsy.max=3, # current FがFmsyに比べて小さすぎる場合、うまく収束しない場合があります。そのときはこのオプションでFmsy.max=10とかしてください。
                     Fmsy.step=0.1,thin=1,inc=1,SRtype="L2",fm=5,tol=NULL,
                     AutoCor=FALSE,# 関数内部で自己相関係数を推定するか "future.vpa"を使う場合はどちらでも良い
                     AutoCorOut=FALSE, # フィットさせたあと残差の自己相関を計算する場合
                     current.resid=0, # 最近年何年分の自己相関を平均するか
                     future.function.name="future.vpa1",seed=1){

  if (is.null(tol)) tol <- .Machine$double.eps^0.25

  Ccur <- sum(tail(t(vpares$input$dat$caa*vpares$input$dat$waa),1),na.rm=TRUE)
  Blim.cur <- vpares$Blim

  A <- nrow(vpares$input$dat$caa)

  if (is.na(vpares$Fc.at.age[length(vpares$Fc.at.age)])){
    vpares$input$dat$caa <- vpares$input$dat$caa[-A,]
    vpares$input$dat$waa <- vpares$input$dat$waa[-A,]
    vpares$input$dat$maa <- vpares$input$dat$maa[-A,]
    vpares$input$dat$M <- vpares$input$dat$M[-A,]
    vpares$naa <- vpares$naa[-A,]
    vpares$faa <- vpares$faa[-A,]
    vpares$Fc.at.age <- vpares$Fc.at.age[-A]
    vpares$ssb[is.na(vpares$ssb)] <- 0
    A <- A-1
  }

  vpares$Fc.at.age <- fm*vpares$Fc.at.age

  SRdata <- get.SRdata(vpares)

# fit SR

  if (is.null(res1) && is.null(sim0)){
      res0 <- estSR(SRdata,SR="HS",type=SRtype,Length=20,rho.range=0,AutoCor=FALSE)
      if (AutoCor){
          res1 <- estSR(SRdata,SR="HS",type=SRtype,Length=20,rho.range=0,AutoCor=AutoCor)
          if (res0$AICc <= res1$AICc) res1 <- res0
#          if(!is.null(res0$aic)) if(res0$aic <= res1$aic) res1 <- res0
#          if(!is.null(res0$AICc)) if(res0$AICc <= res1$AICc) res1 <- res0
      } else{
          res1 <- res0
      }
  }

  if (class(res1$pars)!="numeric") res1$pars <- as.numeric(res1$pars)

  if (AutoCorOut){
    ar1 <- ar(ts(res1$resid),order.max=1)
    rho <- ar1$ar
    if (length(rho)==0) rho <- 0
    if (abs(rho) >= 0.99) rho <- sign(rho)*0.99
    res1$pars[3] <- sqrt(ar1$var.pred)
    res1$pars[4] <- rho
  }

  if (current.resid > 0) w.recent <- mean(rev(res1$resid)[1:current.resid]) else w.recent <- 0

# Initial Setting

    lag <- as.numeric(rownames(vpares$input$dat$caa))[1]

    Pope <- vpares$input$Pope
    future.vpa1 <- get(future.function.name)

    years <- sort(as.numeric(rev(names(vpares$naa))[1:5]))

    Surv <- exp(-tail(t(vpares$input$dat$M),1))

    L <- cumprod(Surv)
    L[A-1] <- L[A-1]/(1-Surv[A])
    L <- c(1,L[1:(A-1)])

    GT <- Generation.Time(vpares,maa.year=years,M.year=years)  # Generation Time
    GT2 <- round(GT*2)

    det.naa0 <- res1$pars[1]*res1$pars[2]*L

    waa <- vpares$input$dat$waa
    waa <- rowMeans(waa[,colnames(waa) %in% years])
    maa <- vpares$input$dat$maa
    maa <- rowMeans(maa[,colnames(maa) %in% years])

    det.B0 <-  sum(det.naa0*waa*maa)

    if(is.null(nyear)) nyear <- round(GT*long.term)

    if (is.null(sim0)){
       sim0 <- future.vpa1(vpares,
                   multi=0,
                   nyear=nyear, # 将来予測の年数
                   N=N, # 確率的計算の繰り返し回数
                   ABC.year=max(years)+1, # ABCを計算する年
                   waa.year=years, # 生物パラメータの参照年
                   maa.year=years,
                   M.year=years,
                   seed=seed,
                   naa0=det.naa0,
                   Pope=Pope,
                   # recfuncに対する引数
                   rec.arg=list(a=res1$pars[1],b=res1$pars[2],gamma=res1$gamma,sd=res1$pars[3],bias.correction=TRUE,rho=res1$pars[4],resid=res1$resid)
      )

    sim1 <- future.vpa1(vpares,
                   multi=1,
                   nyear=nyear, # 将来予測の年数
                   N=1, # 確率的計算の繰り返し回数
                   ABC.year=max(years)+1, # ABCを計算する年
                   waa.year=years, # 生物パラメータの参照年
                   maa.year=years,
                   M.year=years,
                   seed=seed,
                   naa0=det.naa0,
                   Pope=Pope,
                   # recfuncに対する引数
                   rec.arg=list(a=res1$pars[1],b=res1$pars[2],gamma=res1$gamma,sd=res1$pars[3],bias.correction=TRUE,rho=res1$pars[4],resid=res1$resid)
                   )
    } else{
        farg <- sim0$input
        farg$N <- N
        farg$nyear <- nyear
        farg$multi <- 0
        farg$ABC.year <- max(years)+1
        farg$naa0 <- det.naa0
        if(!is.null(farg$pre.catch)){
            farg$pre.catch <- NULL # pre.catchオプションがあるとうまくいかないのでなかったことにする
            cat("notice: option \"pre.catch\" is turned off in estimating MSY.\n")
        }
        if(!is.null(farg$new.rec)){
            farg$rec.new <- NULL # rec.newプションがあるとうまくいかないのでなかったことにする
            cat("notice: option \"rec.new\" is turned off in estimating MSY.\n")
        }
        farg$add.year <- 1
        farg$is.plot <- FALSE
        farg$silent <- TRUE
        farg$det.run <- FALSE
        sim0 <- do.call(future.vpa1,farg)

        farg$N <- 2
        farg$multi <- 1
        sim1 <- do.call(future.vpa1,farg)
    }

##    MSY推定

    farg <- sim1$input
    nY <- nyear+1
    eyear <- mY+(lag > 0)*(lag-1)

    syfunc <- function(x,farg,nyear=50,N=100,eyear=4,naa0=NULL,eaa0=NULL,ssb0=NULL,faa0=NULL,sd=NULL){
      farg$multi <- x
      farg$N <- N
      farg$nyear <- nyear
      farg$naa0 <- naa0
      farg$eaa0 <- eaa0
      farg$ssb0 <- ssb0
      farg$faa0 <- faa0
      if (!is.null(sd)) farg$rec.arg$sd <- sd
      fout <- do.call(future.vpa1,farg)

      nY <- nyear+1

      out <- list(catch=fout$vwcaa[(nY-(eyear-1)):nY,,drop=FALSE],ssb=fout$ssb[,(nY-(eyear-1)):nY,,drop=FALSE],naa=fout$naa[,(nY-(eyear-1)):nY,,drop=FALSE],baa=fout$baa[,(nY-(eyear-1)):nY,,drop=FALSE],eaa=fout$eaa[(nY-(eyear-1)):nY,,drop=FALSE])
      return(out)
    }

    F.multi <- seq(0,Fmsy.max,by=Fmsy.step)

    N0 <- sim0$naa[,nY,]
    e0 <- sim0$eaa[nY,]

    if(lag==0) SSB0 <- NULL else SSB0 <- sim0$ssb[,nY-(lag-1),]

    FSYest <- lapply(F.multi, function(x) syfunc(x,farg,N=round(N/thin),nyear=nyear,eyear=eyear,naa0=N0[,1:round(N/thin)],eaa=e0[1:round(N/thin)],ssb0=SSB0[,1:round(N/thin)]))

    FSYest.c <- sapply(1:length(F.multi), function(i) mean(FSYest[[i]]$catch))

    num.msy <- which.max(FSYest.c)

    num.msy0 <- num.msy

    Fmsy.multi <- F.multi[num.msy]

    obj.msy <- function(x) -mean(syfunc(x,farg,N=N,nyear=nyear,eyear=eyear,naa0=N0,eaa0=e0,ssb0=SSB0)$catch)
    res.msy <- optimize(obj.msy, pmin(pmax(c(Fmsy.multi-inc*Fmsy.step,Fmsy.multi+inc*Fmsy.step),min(F.multi)),max(F.multi)),tol=tol)

    Fmsy.multi <- res.msy$minimum
    MSY <- -res.msy$objective
    MSYres <- syfunc(Fmsy.multi,farg,N=N,nyear=nyear,eyear=eyear,naa0=N0,eaa0=e0,ssb0=SSB0)
    ssb.msy <- mean(apply(MSYres$ssb,c(2,3),sum,na.rm=TRUE))

##    PGY推定

    id.pgy0 <- num.msy0:length(F.multi)
    id.pgy <- which.min((FSYest.c[id.pgy0] - pgy*MSY)^2)

    pgy.low <- id.pgy0[id.pgy]

    Flow <- F.multi[pgy.low]

    N.m <- MSYres$naa[,eyear,]
    e.m <- MSYres$eaa[eyear,]
    if(lag==0) SSB.m <- NULL else SSB.m <- MSYres$ssb[,eyear-(lag-1),]

    obj.pgy <- function(x) (mean(syfunc(x,farg,N=N,nyear=nyear,eyear=eyear,naa0=N.m,eaa0=e.m,ssb0=SSB.m)$catch)-pgy*MSY)^2
    res.pgy <- optimize(obj.pgy, pmin(pmax(c(Flow-inc*Fmsy.step,Flow+inc*Fmsy.step),Fmsy.multi),max(F.multi)),tol=tol)

    Flow.multi <- res.pgy$minimum
    PGYlow.res <- syfunc(Flow.multi,farg,N=N,nyear=nyear,eyear=eyear,naa0=N.m,eaa0=e.m,ssb0=SSB.m)
    PGYlow <- mean(PGYlow.res$catch)
    ssb.low <- mean(apply(PGYlow.res$ssb,c(2,3),sum,na.rm=T))

    id.pgy0 <- 1:num.msy0
    id.pgy <- which.min((FSYest.c[id.pgy0] - pgy*MSY)^2)

    pgy.high <- id.pgy0[id.pgy]

    Fhigh <- F.multi[pgy.high]

    res.pgy <- optimize(obj.pgy, pmin(pmax(c(Fhigh-inc*Fmsy.step,Fhigh+inc*Fmsy.step),min(F.multi)),Fmsy.multi),tol=tol)

    Fhigh.multi <- res.pgy$minimum
    PGYhigh.res <- syfunc(Fhigh.multi,farg,N=N,nyear=nyear,eyear=eyear,naa0=N.m,eaa0=e.m,ssb0=SSB.m)
    PGYhigh <- mean(PGYhigh.res$catch)
    ssb.high <- mean(apply(PGYhigh.res$ssb,c(2,3),sum,na.rm=T))

##  Bhs推定

    if(res1$input$SR=="HS"){
        det.Bhs <- res1$pars[2]

        if (!is.null(SSB0)) SSB0.HS <- as.matrix(rowMeans(SSB0)) else SSB0.HS <- NULL

        obj.HS <- function(x) sum(syfunc(x,farg,N=2,nyear=nyear,eyear=eyear,naa0=as.matrix(rowMeans(N0)),eaa0=NULL,ssb0=SSB0.HS,sd=0)$ssb[,eyear,1])-det.Bhs
    res.HS <- uniroot(obj.HS, c(0,2*max(F.multi)),tol=tol)

        FHS.multi <- res.HS$root

        HSres <- syfunc(FHS.multi,farg,N=N,nyear=nyear,eyear=eyear,naa0=N0,eaa0=e0,ssb0=SSB0)
        HScat <- mean(HSres$catch)
        ssb.hs <- mean(apply(HSres$ssb,c(2,3),sum,na.rm=T))
    }
    else{
        det.Bhs <- SSB0.HS <- obj.HS <- res.HS <- FHS.multi <- HSres <- HScat <- ssb.hs <- NULL
    }

##  target function

    # 実際には最後の再生産関係の残差も入れてやる必要がある（自己相関を入れるときに必要）

    ##    target.func <- function(x,farg,naa0=NULL,eaa0=NULL,ssb0=NULL,faa0=NULL,mY=5,N=1,seed=1,eyear=4,p=1,beta=1,delta=0,Blim=0,Bban=0,sd0=NULL){
    target.func <- function(x,farg,naa0=NULL,eaa0=NULL,ssb0=NULL,faa0=NULL,mY=5,N=2,seed=1,eyear=4,p=1,beta=1,delta=0,Blim=0,Bban=0,sd0=NULL){
      farg$multi <- x
      farg$seed <- seed
      farg$N <- N
      farg$nyear <- mY
      farg$naa0 <- p*naa0
      farg$eaa0 <- eaa0
      farg$ssb0 <- p*ssb0
      farg$faa0 <- faa0
      farg$beta <- beta
      farg$delta <- delta
      farg$Blim <- Blim
      farg$Bban <- Bban
      if(!is.null(farg$ABC.year)) farg$ABC.year <- farg$start.year
      if (!is.null(sd0)) farg$rec.arg$sd <- sd0
      fout <- do.call(future.vpa1,farg)

      nY <- mY+1

      out <- list(catch=fout$vwcaa[(nY-(eyear-1)):nY,,drop=FALSE],ssb=fout$ssb[,(nY-(eyear-1)):nY,,drop=FALSE],naa=fout$naa[,(nY-(eyear-1)):nY,,drop=FALSE],baa=fout$baa[,(nY-(eyear-1)):nY,,drop=FALSE],eaa=fout$eaa[(nY-(eyear-1)):nY,,drop=FALSE])
      return(out)
    }

##  Blim0推定

    Lim0.res <- LIMtoLOW <- list()
    Flim.multi <- Lim0 <- ssb.lim0 <- PRT.lim <- numeric(length(lim))

    PRT.range <- 1:(4*GT2)

    N.m <- MSYres$naa[,eyear,]
    e.m <- MSYres$eaa[eyear,]
    if(lag==0) SSB.m <- NULL else SSB.m <- MSYres$ssb[,eyear-(lag-1),]

    for (j in 1:length(lim)){
      id.lim0 <- num.msy0:length(F.multi)
      id.lim <- which.min((FSYest.c[id.lim0] - lim[j]*MSY)^2)

      lim.num <- id.lim0[id.lim]

      Flim <- F.multi[lim.num]

      obj.lim <- function(x) (mean(syfunc(x,farg,N=N,nyear=nyear,eyear=eyear,naa0=N.m,eaa0=e.m,ssb0=SSB.m)$catch)-lim[j]*MSY)^2
      res.lim <- optimize(obj.lim, pmin(pmax(c(Flim-inc*Fmsy.step,Flim+inc*Fmsy.step),Fmsy.multi),max(F.multi)),tol=tol)

      Flim.multi[j] <- res.lim$minimum
      Lim0.res[[j]] <- syfunc(Flim.multi[[j]],farg,N=N,nyear=nyear,eyear=eyear,naa0=N.m,eaa0=e.m,ssb0=SSB.m)
      Lim0[j] <- mean(Lim0.res[[j]]$catch)
      ssb.lim0[j] <- mean(apply(Lim0.res[[j]]$ssb,c(2,3),sum,na.rm=T))

      ## PRT

        N.p <- Lim0.res[[j]]$naa[,eyear,]
        e.p <- Lim0.res[[j]]$eaa[eyear,]

        if(lag==0) SSB.p <- NULL else SSB.p <- Lim0.res[[j]]$ssb[,eyear-(lag-1),]

        LIMtoLOW[[j]] <- target.func(Fmsy.multi,farg,mY=4*GT2,seed=seed,N=N,eyear=4*GT2,naa0=N.p,eaa0=e.p,ssb0=SSB.p)

        LIMtoLOW.ssb <- sapply(PRT.range, function(x) mean(colSums(LIMtoLOW[[j]]$ssb[,x,])))

        PRT.lim[j] <- min(which(LIMtoLOW.ssb >= ssb.low))
    }

    ## PRT.lim <= mYを満たすものが一個もなくwarningを返すことがあるが、放置しています（使っていない計算結果なので）
    nlim.est <- min(which(PRT.lim <= mY))
    if(is.na(nlim.est) | nlim.est == Inf) nlim.est <- length(lim)

    lim1 <- lim[nlim.est]
    Lim1.res <- Lim0.res[[nlim.est]]
    Flim1.multi <- Flim.multi[nlim.est]
    Lim1 <- Lim0[nlim.est]
    ssb.lim1 <- ssb.lim0[nlim.est]
    PRT.lim1 <- PRT.lim[nlim.est]

##  Bban0推定

    Ban0.res <- BANtoLIM <- list()
    Fban.multi <- Ban0 <- ssb.ban0 <- PRT.ban <- numeric(length(ban))

    for (j in 1:length(ban)){
      id.ban0 <- num.msy0:length(F.multi)
      id.ban <- which.min((FSYest.c[id.ban0] - ban[j]*MSY)^2)

      ban.num <- id.ban0[id.ban]

      Fban <- F.multi[ban.num]

      obj.ban <- function(x) (mean(syfunc(x,farg,N=N,nyear=nyear,eyear=eyear,naa0=N.m,eaa0=e.m,ssb0=SSB.m)$catch)-ban[j]*MSY)^2
      res.ban <- optimize(obj.ban, pmin(pmax(c(Fban-inc*Fmsy.step,Fban+inc*Fmsy.step),Fmsy.multi),max(F.multi)),tol=tol)

      Fban.multi[j] <- res.ban$minimum
      Ban0.res[[j]] <- syfunc(Fban.multi[j],farg,N=N,nyear=nyear,eyear=eyear,naa0=N.m,eaa0=e.m,ssb0=SSB.m)
      Ban0[j] <- mean(Ban0.res[[j]]$catch)
      ssb.ban0[j] <- mean(apply(Ban0.res[[j]]$ssb,c(2,3),sum,na.rm=T))

      ## PRT

        N.p <- Ban0.res[[j]]$naa[,eyear,]
        e.p <- Ban0.res[[j]]$eaa[eyear,]

        if(lag==0) SSB.p <- NULL else SSB.p <- Ban0.res[[j]]$ssb[,eyear-(lag-1),]

        BANtoLIM[[j]] <- target.func(Fmsy.multi,farg,mY=4*GT2,seed=seed,N=N,eyear=4*GT2,naa0=N.p,eaa0=e.p,ssb0=SSB.p,delta=1,Blim=ssb.lim1,Bban=ssb.ban0[j])

        BANtoLIM.ssb <- sapply(PRT.range, function(x) mean(colSums(BANtoLIM[[j]]$ssb[,x,])))

        PRT.ban[j] <- min(which(BANtoLIM.ssb >= ssb.lim1))
    }

    nban.est <- min(which(PRT.ban <= mY))
    if(is.na(nban.est) | nban.est == Inf) nban.est <- length(ban)

    ban1 <- ban[nban.est]
    Ban1.res <- Ban0.res[[nban.est]]
    Fban1.multi <- Fban.multi[nban.est]
    Ban1 <- Ban0[nban.est]
    ssb.ban1 <- ssb.ban0[nban.est]
    PRT.ban1 <- PRT.ban[nban.est]

## Btarget推定

    ## Btarget

    Ftar.multi <- Fmsy.multi

    TARres <- target.func(Ftar.multi,farg,mY=mY,seed=seed,N=N,eyear=mY,naa0=N.m,eaa0=e.m+w.recent,ssb0=SSB.m)
    Btar <- mean(colSums(TARres$ssb[,mY,]))

    # Blow 推定

    N.low <- PGYlow.res$naa[,1+(lag>0)*(lag-1),]
    e.low <- PGYlow.res$eaa[1+(lag>0)*(lag-1),]
    if(lag==0) SSB.low <- NULL else SSB.low <- PGYlow.res$ssb[,1,]

    LOWres <- target.func(Flow.multi,farg,mY=mY,seed=seed,N=N,eyear=mY,naa0=N.low,eaa0=e.low+w.recent,ssb0=SSB.low)

    Blow <- mean(colSums(LOWres$ssb[,mY,]))

    P.low <- Blow/Btar

    # Blim 推定

    N.lim <- Lim1.res$naa[,1+(lag>0)*(lag-1),]
    e.lim <- Lim1.res$eaa[1+(lag>0)*(lag-1),]
    if(lag==0) SSB.lim <- NULL else SSB.lim <- Lim1.res$ssb[,1,]

    LIMres <- target.func(Flim1.multi,farg,mY=mY,seed=seed,N=N,eyear=mY,naa0=N.lim,eaa0=e.lim+w.recent,ssb0=SSB.lim)

    Blim <- mean(colSums(LIMres$ssb[,mY,]))

    P.lim <- Blim/Btar

    # Bban 推定

    N.ban <- Ban1.res$naa[,1+(lag>0)*(lag-1),]
    e.ban <- Ban1.res$eaa[1+(lag>0)*(lag-1),]
    if(lag==0) SSB.ban <- NULL else SSB.ban <- Ban1.res$ssb[,1,]

    BANres <- target.func(Fban1.multi,farg,mY=mY,seed=seed,N=N,eyear=mY,naa0=N.ban,eaa0=e.ban+w.recent,ssb0=SSB.ban)

    Bban <- mean(colSums(BANres$ssb[,mY,]))

    P.ban <- Bban/Btar

    Pref <- c(P.low, P.lim, P.ban)
    names(Pref) <- c("Low","Lim","Ban")

    if(0){
        # plotはしない
        x.range <- range(SRdata$SSB,Btar,Bban)
        plot(SRdata$SSB, SRdata$R,xlab="SSB",ylab="R",pch=16,col="blue",cex=1.2,xlim=x.range)
        x.SSB <- seq(0,x.range[2],len=100)
        lines(x.SSB,HS(x.SSB,res1$pars[1],res1$pars[2]),col="pink",lwd=3)
        abline(v=Btar,col="green",lwd=3,lty=2)
        abline(v=Blim,col="orange",lwd=3,lty=2)
        abline(v=Bban,col="red",lwd=3,lty=2)
    }

    if(0){
      argname <- ls()
      arglist <- lapply(argname,function(x) eval(parse(text=x)))
      names(arglist) <- argname
      arglist$MSYres <- NULL
      arglist$PGYlow.res <- NULL
      invisible(arglist)
    }

    out <- list(
      stockid=vpares$stockid,
      seed=seed,
      SRdata=SRdata,
      res1=res1,
      det.B0=det.B0,
      det.Bhs=det.Bhs,
      B0=mean(colSums(sim0$ssb[,nY,])),
      Bmsy=ssb.msy,
      Bpgy.low=ssb.low,
      Bpgy.high=ssb.high,
      Blim1=ssb.lim1,
      Bban1=ssb.ban1,
      Btar=Btar,
      Blow=Blow,
      Blim=Blim,
      Bban=Bban,
      Bhs=ssb.hs,
      pgy=pgy,
      lim=lim,
      ban=ban,
      fm=fm,
      lag=lag,
      eyear=eyear,
      Fmsy=Fmsy.multi,
      Flow=Flow.multi,
      Fhigh=Fhigh.multi,
      Flim=Flim1.multi,
      Fban=Fban1.multi,
      Fhs=FHS.multi,
      farg=farg,
      N=N,
      GT=GT,
      sim0=sim0,
      w.recent=w.recent,
      nyear=nyear,
      MSYres=MSYres,
      lim.est=lim1,
      ban.est=ban1,
      PRT.lim0=PRT.lim,
      PRT.ban0=PRT.ban,
      PRT.lim=PRT.lim1,
      PRT.ban=PRT.ban1,
      PGYlow.res=PGYlow.res,
      PGYhigh.res=PGYhigh.res,
      TARres=TARres,
      LOWres=LOWres,
      Lim0.res=Lim0.res,
      Ban0.res=Ban0.res,
      Lim1.res=Lim1.res,
      Ban1.res=Ban1.res,
      LIMres=LIMres,
      BANres=BANres,
      HSres=HSres,
      Pref=Pref,
      syfunc=syfunc,
      target.func=target.func,
      future.function.name=future.function.name,
      Blim.cur=Blim.cur,
      Ccur=Ccur
    )

    b.table <- unlist(out[c("Bmsy","Btar","Fmsy","Bpgy.low","Blow","Flow",
                             "Blim1","Blim","Flim","Bban1","Bban","Fban")])
#    names(b.table) <- c("Bmsy\n(Equiribrium)","Bmsy\n(with AR)=Btarget","Fmsy",
#                         "Bpgy90%-low\n(Equiribrium)","Bpgy90%-low\n(with AR)=Blow","Flow")
#    b.limit <- c(out$Blim0,out$Blim,out$Flim,out$Bban0,out$Bban,out$Fban)
    b.table <- t(matrix(b.table,3,4))
#    b.table[,1:2] <- b.table[,1:2]/1000
    b.table <- cbind(apply(b.table[,1:2],2, function(x) round(as.numeric(x),0)),
                     round(b.table[,3],2))
    w.recent2 <- ifelse(is.null(sim0$input$rec.arg$rho),NA,w.recent)
    b.table <- rbind(b.table,c(NA,w.recent2,NA))
    colnames(b.table) <- c("Equiribrium","with AR","Fref/Fcurrent")
    rownames(b.table) <- c("Bmsy","B_pgy_90%_L","B_limit (B_pgy_60%_L)","B_ban (B_pgy_10%_L)","Recent residual")

    out$summary <- b.table

    return(out)
}
