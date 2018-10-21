#---------------- 管理基準値計算のための関数 ------------------------
# ref.F
ref.F <- function(
  res, # VPAの結果のオブジェクト
  sel=NULL, # 仮定する選択率．NULLの場合，res$Fc.at.ageが使われる
  waa=NULL, # 仮定する生物パラメータ．直接の値を入れるか，年を指定するやり方のどちらでも動く。直接指定するほうが優先。
  maa=NULL,
  M=NULL,
  waa.catch=NULL,
  M.year=NULL,
  waa.year=NULL, # 年を指定して生物パラメータを仮定する場合．年の範囲の平均値が用いられる．NULLの場合，VPA最終年の値が使われる
  maa.year=NULL,
  rps.year = NULL, # Fmedの計算に使うRPSの年の範囲．NULLの場合，全範囲が用いられる
  max.age = Inf, # 加入年齢を０歳としたときに、SPR計算で考慮される最大の年齢（年齢の数ではないことに注意）。加入年齢が１歳以上のときは、SPR計算で考慮したい年齢-加入年齢を入力する、またはmin.ageの引数に加入年齢を設定する。
  min.age = 0, # 加入年齢が0歳でないときに指定できる
  d = 0.001,
  Fspr.init = 0.5, # F%SPRの初期値
  Fmax.init = 1.5, # Fmaxの初期値
  F0.1.init = 0.7, # F0.1の初期値
  pSPR = seq(10,90,by=10), # F%SPRを計算するときの％SPR
  iterlim=1000,
  plot=TRUE,
  Pope=FALSE, # 2014.7.4追加
  F.range = seq(from=0,to=2,length=101)  # YPR, SPR曲線を書くときのFの範囲
){

    argname <- ls()
    arglist <- lapply(argname,function(x) eval(parse(text=x)))
    names(arglist) <- argname

#  caa <- res$input$dat$caa
  naa <- res$naa
  ssb <- res$ssb
  ny <- ncol(naa)
  years <- dimnames(naa)[[2]]
  ages <- dimnames(naa)[[1]]

  if(is.null(sel)){
    Fc.at.age <- res$Fc.at.age
    sel <- Fc.at.age/max(Fc.at.age,na.rm=TRUE)
  }
  else{
    Fc.at.age <- sel
  }
  sel <- sel/max(sel,na.rm=T)

  na <- sum(!is.na(sel))

  if(is.null(waa.year)) waa.year <- rev(years)[1]
  if(is.null(maa.year)) maa.year <- rev(years)[1]
  if(is.null(M.year)) M.year <- rev(years)[1]
  if(is.null(rps.year)) rps.year <- as.numeric(colnames(res$naa))

  if(is.null(waa))  waa <- apply(as.matrix(as.data.frame(res$input$dat$waa)[as.character(waa.year)]),1,mean)
  if(is.null(M))  M <- apply(as.matrix(as.data.frame(res$input$dat$M)[as.character(M.year)]),1,mean)
  if(is.null(maa))  maa <- apply(as.matrix(as.data.frame(res$input$dat$maa)[as.character(maa.year)]),1,mean)

  if(is.null(waa.catch)){
      if(is.null(res$input$dat$waa.catch)){
          waa.catch <- waa
      }
      else{
          waa.catch <- apply(as.matrix(as.data.frame(res$input$dat$waa.catch)[as.character(waa.year)]),1,mean)
          }
  }

  ssb.coef <- ifelse(is.null(res$ssb.coef),0,res$ssb.coef)

  min.age <- min(as.numeric(rownames(res$naa)))
  if(min.age==0) slide.tmp <- TRUE else slide.tmp <- -1:-min.age
  rps.data <- data.frame(year=as.numeric(names(colSums(ssb,na.rm=T))),
                         ssb=as.numeric(colSums(ssb,na.rm=T)),
                         recruit=as.numeric(c(naa[1,slide.tmp],rep(NA,min.age))))
  if (sum(is.na(rps.data$year))>0) rps.data <- rps.data[-which(is.na(rps.data$year)),]
  rps.data$rps <- rps <- rps.data$recruit/rps.data$ssb
#  rps <- as.numeric(naa[1,]/colSums(ssb,na.rm=TRUE))

#  if (is.null(rps.year)) rps.year <- years

  tmp <- rps.data$year %in% rps.year
  rps.q <- quantile(rps[tmp], na.rm=TRUE, probs=c(0.1,0.5,0.9))
  rps.q <- c(rps.q,mean(rps[tmp], na.rm=TRUE))
  names(rps.q)[4] <- "mean"
  spr.q <- 1/rps.q

#  browser()
  # F.spr

  spr.f.est <- function(log.p, out=FALSE, sub="med", spr0=NULL){
    Fr <- exp(log.p)

    tmp <- calc.rel.abund(sel,Fr,na,M,waa,waa.catch,maa,min.age=min.age,max.age=max.age,Pope=Pope,ssb.coef=ssb.coef)
    rel.abund <- tmp$rel.abund
    spr <- sum(tmp$spr)
#    rel.abund <- rep(NA, na)
#    rel.abund[1] <- 1
#    for (i in 2:(na-1)) {
#      rel.abund[i] <- rel.abund[i-1]*exp(-M[i-1]-sel[i-1]*Fr)
#    }
#    rel.abund[na] <- rel.abund[na-1]*exp(-M[na-1]-sel[na-1]*Fr)*(1-exp(-(max.age-(na-2))*(M[na]+sel[na]*Fr)))/(1-exp(-M[na]-sel[na]*Fr))

#    spr <- sum(rel.abund*waa*maa)

    if (isTRUE(out)) obj <- spr
    else{
#      browser()
     if(sub=="mean") obj <- (spr-spr.q[4])^2
     if(sub=="low") obj <- (spr-spr.q[3])^2
     if(sub=="med") obj <- (spr-spr.q[2])^2
     if(sub=="high") obj <- (spr-spr.q[1])^2
     if(is.numeric(sub)) obj <- (spr/spr0-sub/100)^2

    }

    return(obj)
  }

  spr0 <- spr.f.est(-Inf, out=TRUE)

  Fmed.res <- nlm(spr.f.est, Fspr.init, out=FALSE, sub="med", iterlim = iterlim)
  Fmean.res <- nlm(spr.f.est, Fspr.init, out=FALSE, sub="mean", iterlim = iterlim)
#  browser()
  Flow.res <- nlm(spr.f.est, Fspr.init, out=FALSE, sub="low", iterlim = iterlim)
  Fhigh.res <- nlm(spr.f.est, Fspr.init, out=FALSE, sub="high", iterlim = iterlim)

  Fmean <- exp(Fmean.res$estimate)
  Fmed <- exp(Fmed.res$estimate)
  Flow <- exp(Flow.res$estimate)
  Fhigh <- exp(Fhigh.res$estimate)

  if (!is.null(pSPR)){
    FpSPR <- NULL

    for (i in pSPR){
      FpSPR.res <- nlm(spr.f.est, Fspr.init, out=FALSE, sub=i, spr0=spr0, iterlim=iterlim)
#      print(FpSPR.res)
#       cat("i", FpSPR.res$code," ")
      FpSPR <- c(FpSPR, exp(FpSPR.res$estimate))
    }
    names(FpSPR) <- paste(pSPR,"%SPR",sep="")
  }

  # Fmax

  ypr.f.est <- function(log.p, out=FALSE){
    Fr <- exp(log.p)

#    rel.abund <- rep(NA, na)
#    rel.abund[1] <- 1
#    for (i in 2:(na-1)) {
#      rel.abund[i] <- rel.abund[i-1]*exp(-M[i-1]-sel[i-1]*Fr)
#    }
#    rel.abund[na] <- rel.abund[na-1]*exp(-M[na-1]-sel[na-1]*Fr)*(1-exp(-(max.age-(na-2))*(M[na]+sel[na]*Fr)))/(1-exp(-M[na]-sel[na]*Fr))
#    ypr <- sum(rel.abund*waa*(1-exp(-sel*Fr))*exp(-M/2))
    tmp <- calc.rel.abund(sel,Fr,na,M,waa,waa.catch,maa,max.age=max.age,Pope=Pope,ssb.coef=ssb.coef)
    rel.abund <- tmp$rel.abund
    ypr <- sum(tmp$ypr)

    if (isTRUE(out)) obj <- ypr else obj <- -ypr

   return(obj)
  }

  Fmax.res <- nlm(ypr.f.est, log(Fmax.init), out=FALSE)

  Fmax <- exp(Fmax.res$estimate)

  # F0.1

  Fp <- function(log.p, out=FALSE){
    Fr <- exp(log.p)

    tmp <- calc.rel.abund(sel,Fr,na,M,waa,waa.catch,maa,max.age=max.age,Pope=Pope,ssb.coef=ssb.coef)
    rel.abund <- tmp$rel.abund
    ypr <- sum(tmp$ypr)

#    rel.abund <- rep(NA, na)
#    rel.abund[1] <- 1
#    for (i in 2:(na-1)) {
#      rel.abund[i] <- rel.abund[i-1]*exp(-M[i-1]-sel[i-1]*Fr)
#    }
#   rel.abund[na] <- rel.abund[na-1]*exp(-M[na-1]-sel[na-1]*Fr)*(1-exp(-(max.age-(na-2))*(M[na]+sel[na]*Fr)))/(1-exp(-M[na]-sel[na]*Fr))
#   ypr <- sum(rel.abund*waa*(1-exp(-sel*Fr))*exp(-M/2))

   if (isTRUE(out)) obj <- ypr else obj <- -ypr

   return(obj)
  }

  F0.1.est <- function(log.p){
    p <- exp(log.p)
    ref.trend <- (Fp(log(d))-Fp(log(0)))/d
    trend <- (Fp(log(p+d)) - Fp(log(p)))/d

    obj <- (ref.trend/10 - trend)^2

    obj
  }

  F0.1.res <- nlm(F0.1.est,log(F0.1.init))

  F0.1 <- exp(F0.1.res$estimate)

  # Fcurrent
  Fcurrent <- c(max(Fc.at.age,na.rm=T), mean(Fc.at.age,na.rm=T))

  # output
  f.mean <- function(x) mean(x*sel, na.rm=T)

  Fmean <- c(Fmean, f.mean(Fmean))
  Fmed <- c(Fmed, f.mean(Fmed))
  Flow <- c(Flow, f.mean(Flow))
  Fhigh <- c(Fhigh, f.mean(Fhigh))
  Fmax <- c(Fmax, f.mean(Fmax))
  F0.1 <- c(F0.1, f.mean(F0.1))

  names(Fcurrent) <- names(Fmed) <- names(Fmean) <- names(Flow) <- names(Fhigh) <- names(Fmax) <- names(F0.1) <- c("max","mean")

  Res <- list(sel=sel, min.age=min.age, max.age=max.age, rps.q=rps.q, spr.q=spr.q, Fcurrent=Fcurrent, Fmed=Fmed, Flow=Flow, Fhigh=Fhigh, Fmax=Fmax, F0.1=F0.1, Fmean=Fmean,rps.data=rps.data)

  if (!is.null(pSPR)){
    FpSPR <- rbind(FpSPR, sapply(FpSPR, f.mean))
    rownames(FpSPR) <- c("max","mean")
    Res$FpSPR <- FpSPR
  }

  #-----  YPR & SPR figure -----

  spr0 <- sum(calc.rel.abund(sel,0,na,M,waa,waa.catch,maa,min.age=min.age,max.age=max.age,Pope=Pope,ssb.coef=ssb.coef)$spr)
  tmp <- lapply(F.range, function(x) calc.rel.abund(sel,x,na,M,waa,waa.catch,maa,min.age=min.age,max.age=max.age,Pope=Pope,ssb.coef=ssb.coef))
  ypr <- sapply(tmp,function(x) sum(x$ypr))
  spr <- sapply(tmp,function(x) sum(x$spr))/spr0*100

  if (isTRUE(plot)){
  plot(F.range,spr,xlab="F at selectivity=1",ylab="%SPR",type="l",ylim=c(0,max(spr)))
  par(new=T)
  plot(F.range,ypr,axes=F,xlab="",ylab="",lty=2,type="l",ylim=c(0,max(ypr)))
  axis(side=4)
  mtext("YPR",side=4,line=2)
  abline(v=xx <- c(Res$Fmax[1],Res$Fcurrent[1],Res$F0.1[1],Res$Fmed[1]))
  text(xx,rep(max(ypr)*c(0.4,0.3,0.2,0.1),length(xx)),c("Fmax","Fcur","F0.1","Fmed"))
  }

  ypr.spr <- data.frame(F.range=F.range,ypr=ypr,spr=spr)
  Res$ypr.spr  <- ypr.spr
  Res$waa <- waa
  Res$waa.catch <- waa.catch
  Res$maa <- maa
  #------------------------------

  Res$summary <- as.data.frame(Res[substr(names(Res),1,1)=="F"])
  Res$summary <- rbind(Res$summary,Res$summary[1,]/Res$summary[1,1])
  dimnames(Res$summary)[[1]][3] <- "Fref/Fcur"
  Res$arglist <- arglist
  Res$spr0 <- spr0
  class(Res) <- "ref"
#  print(Res)
  return(Res)
}


calc.rel.abund <- function(sel,Fr,na,M,waa,waa.catch=NULL,maa,min.age=0,max.age=Inf,Pope=TRUE,ssb.coef=0){
    if(is.null(waa.catch)) waa.catch <- waa
    rel.abund <- rep(NA, na)
    rel.abund[1] <- 1
    for (i in 2:(na-1)) {
        rel.abund[i] <- rel.abund[i-1]*exp(-M[i-1]-sel[i-1]*Fr)
    }
    rel.abund[na] <- rel.abund[na-1]*exp(-M[na-1]-sel[na-1]*Fr)*(1-exp(-((max.age-min.age)-(na-2))*(M[na]+sel[na]*Fr)))/(1-exp(-M[na]-sel[na]*Fr))

    if(isTRUE(Pope)){
        ypr1 <- rel.abund*waa.catch[1:na]*(1-exp(-sel[1:na]*Fr))*exp(-M[1:na]/2)
    }
    else{
                                        # use Baranov catch equation
        ypr1 <- rel.abund*(1-exp(-sel[1:na]*Fr-M[1:na]))*sel[1:na]*Fr/
            (sel[1:na]*Fr+M[1:na])*waa.catch[1:na]
    }
    spr <- rel.abund*waa[1:na]*maa[1:na]*exp(-ssb.coef*(sel[1:na]*Fr+M[1:na]))
  return(list(rel.abund=rel.abund,ypr=ypr1,spr=spr))
}
