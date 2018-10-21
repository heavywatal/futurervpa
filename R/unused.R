# Functions not appeared elsewhere

sref.F2 <- function(res0,target.year=c(2018,2023),current.year=2011,Blim,
                   interval=c(0,3),...){
  ssb <- apply(res0$ssb,2,sum)
  Frec <- numeric()
  Frec[1] <- ssb[current.year]/Blim

  for(i in 1:length(target.year)){
    tmpfunc <- function(x,res0,Blim,...){
      fres <- future.vpa(res0=res0,multi=x,...)
      cat(x," ")
      return((fres$vssb[rownames(fres$vssb)==target.year[i]]-Blim)^2)
    }
    Frec[i+1] <- optimize(tmpfunc,interval=interval,res0=res0,Blim=Blim,...)$minimum
  }
  return(Frec)
}

# 2012. 8. 3 -- 管理基準値計算は外に出す
getABC <- function(res.vpa, # VPAの結果
                   res.ref, # 管理基準値計算の結果
                   res.future, # 将来予測計算の結果
                   ref.case="all",
                   multi=NULL,
                   N=NULL,
                   SSBcur=1000,
                   Blim=1000,Bban=0,
                   target.year=NULL, # NULLの場合，ABC.year+4
                   catch.year=NULL, # 2013:2017など、漁獲量の平均を出したい期間、NULLの場合、ABC.year:ABC.year+4
                   is.plot=TRUE){
  if(all(ref.case=="all")) ref.case <- names(res.ref$summary)
  if(all(is.null(multi))) multi <- rep(1,length(ref.case))

  nref <- length(ref.case)

  ABC.year <- res.future$input$ABC.year
  if(is.null(target.year)) target.year <- ABC.year+4
  ABC <- wariai <- aveF <- catch5u <- catch5l <- upperSSBlim <- upperSSBcur <- SSBlim <- SSBcur.tmp <- rep(NA,nref)
  names(ABC) <- names(wariai) <- names(aveF) <- paste(ref.case,"x",round(multi,3))
  wcatch <- matrix(NA,5,nref,dimnames=list(((min(ABC.year)):(min(ABC.year)+4)),names(aveF)))

  fres <- list()
  i.tmp <- match(ref.case,names(res.ref$summary))

  if(any(is.na(i.tmp)))
    stop(paste("ref.case specification of is wrong!"))

  years <- res.future$year
  currentF <- res.ref$Fcurrent["max"] * res.ref$sel
  N <- ifelse(is.null(N),dim(res.future$naa)[[3]],N)

  for(i in 1:nref){
    tmp <- res.ref$summary[i.tmp[i]][1,1] * res.ref$sel
    tmp <- max(tmp,na.rm=T)/max(currentF,na.rm=T)*multi[i]
    tmpF <- tmp * currentF
    aveF[i] <- mean(tmpF,na.rm=T)
    input.tmp <- res.future$input
    input.tmp$multi <- tmp
    input.tmp$is.plot <- FALSE
    input.tmp$N <- N

    # Frecで使われたシードはとっておかないといけない=> seedはFrecの引数の外に出すこと！
    input.tmp$Frec <- NULL

    fres[[i]] <- do.call(future.vpa, input.tmp)
    ABC[i] <- fres[[i]]$ABC[1]
#    browser()
    if(res.future$input$ts>1){ # ts>2の場合、漁獲量などの計算は暦年を使う
      input.tmp <- res.future$input
      input.tmp$multi <- tmp
      input.tmp$ts <- 1
      input.tmp$is.plot <- FALSE
      input.tmp$ABC.year <- ABC.year <- floor(min(input.tmp$ABC.year))
      input.tmp$waa <- input.tmp$maa <- input.tmp$M <- input.tmp$waa.catch <- NULL
      input.tmp$N <- N
      fres[[i]] <- do.call(future.vpa, input.tmp)
      years <- fres[[i]]$year
    }
    wariai[i] <- sum(fres[[i]]$wcaa[,years==ABC.year,1],na.rm=T)/
            sum(fres[[i]]$biom[,years==ABC.year,1],na.rm=T)
    catch.year <- (ABC.year):(ABC.year+4)
    wcatch[,i] <- apply(fres[[i]]$vwcaa[years %in% (catch.year),-1],1,mean,na.rm=T)
    catch5u[i] <- quantile(fres[[i]]$vwcaa[years==max(catch.year),-1],probs=0.9) # catchは2017年
    catch5l[i] <- quantile(fres[[i]]$vwcaa[years==max(catch.year),-1],probs=0.1)

    tmp.year <- years %in% target.year
    if(is.null(SSBcur)) SSBcur <- fres[[i]]$vssb[years==(ABC.year),1]

    SSBcur.tmp[i] <- SSBcur
    upperSSBlim[i] <- sum(fres[[i]]$vssb[tmp.year,-1]>Blim)/N*100 # SSBは2018年当初まで
    upperSSBcur[i] <- sum(fres[[i]]$vssb[tmp.year,-1]>SSBcur)/N*100
    SSBlim[i] <- Blim
  }

  if(is.plot){
    par(mfrow=c(1,2),mar=c(4,4,2,1))
    vssb <- apply(res.vpa$ssb,2,sum,na.rm=T)/1000
    x <- sapply(fres,function(x) x$vssb[,1])/1000
    plot(range(c(as.numeric(names(vssb)),years)),
         c(0,max(x)*1.1),type="n",xlab="Year",ylab="SSB (x1000)")
    matpoints(years,x,col=1:nref,type="l",lty=1,
            ylim=c(0,max(x)))
    points(as.numeric(names(vssb)),vssb,type="b")
    abline(h=c(SSBlim/1000,SSBcur/1000),col="gray")
    title("SSB in deterministic runs")
    plot(0,axes=F,xlab="",ylab="")
    legend("topleft",col=1:nref,lty=1,legend=names(ABC))
  }
  average <- apply(wcatch,2,mean)
  res.ref$ABC <- rbind(aveF,wariai,catch5l,catch5u,average,
                         upperSSBcur,SSBcur.tmp,upperSSBlim,SSBlim,ABC)
  rownames(res.ref$ABC)[3] <- paste("catch5l during ",min(catch.year),"-",max(catch.year),sep="")
  rownames(res.ref$ABC)[4] <- paste("catch5u during ",min(catch.year),"-",max(catch.year),sep="")
  rownames(res.ref$ABC)[5] <- paste("average catch during ",min(catch.year),"-",max(catch.year),sep="")
  rownames(res.ref$ABC)[6] <- paste("upperSSBcur at",target.year)
  rownames(res.ref$ABC)[8] <- paste("upperSSBlim at",target.year)
  fres0 <- fres
  write.table(round(res.ref$ABC,2),sep="\t")
  save(fres0,file="fres0.R") # 将来予測の全結果はfres0.Rにてセーブされている

  # Kobe chartの作成
  kobe.array <- array(NA,dim=c(length(fres),nrow(fres[[1]]$vssb),5))
  dimnames(kobe.array) <- list(names(ABC),rownames(fres[[1]]$vssb),
                               c("catch","Biomass","SSB","upperBlimit","upperBban"))
  for(i in 1:length(fres)){
      kobe.array[i,,] <- as.matrix(get.kobematrix(fres[[i]],
                                   Blim=Blim,Bban=Bban,ssb=TRUE))
  }
  return(list(ABC=res.ref$ABC,kobe.matrix=kobe.array))
}


solv.Feq <- function(cvec,nvec,mvec){
  Fres <- rep(0,length(cvec))
 # cat(nvec," ")
  for(i in 1:length(cvec)){
    F0 <- cvec[i]/nvec[i]
    F1 <- cvec[i]*(F0+mvec[i])/nvec[i]/(1-exp(-F0-mvec[i]))
    if(round(cvec[i],6)<round(nvec[i],6)){
      while(abs(F0-F1)>0.0001 ){
        F0 <- F1
        F1 <- cvec[i]*(F0+mvec[i])/nvec[i]/(1-exp(-F0-mvec[i]))
        if(F0-F1==-Inf) cat("\n",cvec[i]," ",nvec[i]," \n")
      }
      Fres[i] <- F1
    }
    else{
      Fres[i] <- 10
      cat("Warning: catch exceeded tot_num at: ",i," ",
          round(cvec[i],6)," ",round(nvec[i],6),"\n")
      }
  }
  Fres
}


## 管理基準値を取り出す関数
get.Bref <- function(res,SRfunc="hs",B0=c(0.3),SPR0=c(0.3),HS=c(1,1.3),PGY=c("PGY_0.9_upper_hs","PGY_0.9_lower_hs")){
    sumref <- res$summary[rownames(res$summary)==SRfunc,]
    refpoints <- list()
    ## MSY管理基準値をピックアップ
    refpoints$BMSY <- sumref$"SSB_MSY"

    ## B0基準の管理基準値はB0×％
    ## B0の値はmout$summary$"B0(SSB)"にある。１番目がHSの結果
    refpoints$B0per <- sumref$"B0(SSB)"[1] * B0 # B0_10,20,30,35,40%の値
    names(refpoints$B0per) <- paste(B0*100,"%",sep="")

    ## B_HS関連の管理基準値
    refpoints$BHS <- sumref$b[1] *  HS
    names(refpoints$BHS) <- paste("B_HSx",HS,sep="")

    ## B_PGY関連の管理基準値(HSをもとにしたものはPGY.biom.hsにあります)
    x <- res$PGY.biom.hs["ssb.mean"]
    refpoints$BPGY <- x[match(PGY,rownames(x)),1]
    names(refpoints$BPGY) <- PGY

    ## SSB_current
    refpoints$SSBcur <- rev(as.numeric(res$vpares$ssb))[1]

    ## SSB_max
    refpoints$SSBmax <- max(as.numeric(res$vpares$ssb))
    return(unlist(refpoints))
}


plot.waa <- function(vres){
    lm.list <- list()
    nage <- nrow(vres$naa)
    col.tmp <- rainbow(nage)
    logx <- log(unlist(vres$naa))
    logy <- log(unlist(vres$input$dat$waa))
    ages <- as.numeric(rep(rownames(vres$naa),ncol(vres$naa)))
    u.age <- unique(ages)
    plot(logx,logy,col=col.tmp[1+ages],xlab="log(N)",ylab="log(weight)")
    for(i in 1:length(u.age)){
        tmp <- ages==u.age[i] & logy>-Inf & logx>-Inf
        if(sum(tmp,na.rm=TRUE)>0){
            lm.list[[i]] <- lm(logy[tmp]~logx[tmp])
            l.type <- ifelse(summary(lm.list[[i]])$coeff[2,4]<0.05,1,2)
            if(!is.na(l.type)) abline(lm.list[[i]],col=col.tmp[1+ages[i]],lty=l.type)
        }
    }
    title(vres$stockid,line=0.2)
    legend("bottomleft",lty=c(1:2,rep(1,nage)),
           col=c(1,1,col.tmp),
           legend=c("p<0.05","p>0.05",paste("Age",u.age)))
    return(lm.list)
}


### dynamics MSYを計算してみる
dyn.msy <- function(naa.past,naa.init=NULL,fmsy,a,b,resid,resid.year,waa,maa,M,SR=TRUE){
    nyear <- length(resid)
    if(is.null(naa.init)) nage <- nrow(naa.past) else nage <- length(naa.init)
    naa <- matrix(0,nage,nyear)
    ssb <- numeric()
    if(is.null(naa.init)) naa[,1] <- naa.past[,colnames(naa.past)==min(resid.year)]
    else naa[,1] <- naa.init
    colnames(naa) <- resid.year
    if(is.null(naa.init)){
        waa <- waa[,colnames(naa.past)%in%resid.year]
        maa <- maa[,colnames(naa.past)%in%resid.year]
        M <- M[,colnames(naa.past)%in%resid.year]
    }
    for(i in 2:nyear){
        ssb[i-1] <- sum(naa[,i-1]*waa[,i-1]*maa[,i-1],na.rm=T)
        if(SR==TRUE){
            naa[1,i] <- HS(ssb[i-1],a,b)*exp(resid[i])
        }
        else{
            naa[1,i] <- naa.past[1,i]
            }
        for(j in 2:(nage-1)) naa[j,i] <- naa[j-1,i-1] * exp(-fmsy[j-1]-M[j-1,i-1])
        naa[nage,i] <- naa[nage-1,i-1] * exp(-fmsy[j-1]-M[j-1,i-1]) + naa[nage,i-1] * exp(-fmsy[nage]-M[nage,i-1])
    }
    i <- nyear ; ssb[i] <- sum(naa[,i]*waa[,i]*maa[,i])
    list(naa=naa,ssb=ssb)
}


plot.HCR <- function(alpha=1,bban=0,blimit=1,btarget=2,add=FALSE,yscale=1.3,xlim=NULL,
                     Fmsy=1,scale=1,
                     ssb.cur=NULL,...) {
    if(is.null(xlim)) xlim <- c(0,max(c(bban,balimit,btarget)))/scale
    b.tmp <- seq(from=0,to=max(xlim),length=300)
    y <- (b.tmp-bban)/(blimit-bban)*alpha
    y <- ifelse(b.tmp>blimit,alpha,y)
    y <- ifelse(y<0,0,y)
    if(!isTRUE(add)) plot(b.tmp/scale,Fmsy*y,type="n",ylim=c(0,yscale),xlab="SSB",ylab="multiplier to current F",xlim=xlim/scale)
    points(b.tmp/scale,Fmsy*y,type="l",...)
    abline(h=Fmsy,col="gray")
    text(0,Fmsy+0.02,paste("Fmsy=",Fmsy,"Fcurrent"),adj=0)
    text(0,Fmsy*alpha+0.02,paste("Ftarget=",round(alpha*Fmsy,2),"Fcurrent (",round(alpha,2),"*Fmsy)",sep=""),adj=0)

    if(!is.null(ssb.cur)){
        Frec <- (ssb.cur-bban)/(blimit-bban)
        Frec <- ifelse(Frec>1,1,Frec)
        lines(c(0,ssb.cur/scale,ssb.cur/scale),c(Frec*alpha*Fmsy,Frec*alpha*Fmsy,0),lty=2)
        points(ssb.cur/scale,Frec*alpha*Fmsy,lty=2,pch=4)
        text(0,0.8*alpha*Fmsy+0.02,
             paste("F2018=",round(Frec*alpha*Fmsy,2),"","Fcurrent (",
                   round(Frec,2),"*",round(alpha,2),"*",round(Fmsy,2),"*Fmsy)",sep=""),adj=0)
    }
}

draw.refline <- function(reftable,horiz=TRUE,scale=1000,lwd=3){
    if(horiz==FALSE){
        abline(v=reftable[1:4,2]/scale,
               col=c("darkgreen","darkblue","darkred","red"),lwd=lwd,lty="22")
        abline(v=reftable[1:4,1]/scale,
               col=c("darkgreen","darkblue","darkred","red"),lwd=lwd,lty=1)
    }
    else{
        abline(h=reftable[1:4,2]/scale,
               col=c("darkgreen","darkblue","darkred","red"),lwd=lwd,lty="22")
        abline(h=reftable[1:4,1]/scale,
               col=c("darkgreen","darkblue","darkred","red"),lwd=lwd,lty=1)
        }
}
