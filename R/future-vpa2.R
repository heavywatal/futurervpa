## est.SR用の将来予測関数。ベクトル化されているので多少速い
future.vpa2 <- function(res0,
                       currentF=NULL, # 管理前のF
                       multi=1, # 管理後（ABC.yearから）のF (current F x multi)
                       nyear=10,Pope=res0$input$Pope,
                        seed=NULL,
                       multi.year=1,#ある特定の年だけFを変えたい場合。デフォルトは1。変える場合は、指定した年またはタイムステップの要素数のベクトルで指定。
                       # 年数の指定
                       start.year=NULL, # 将来予測の開始年，NULLの場合はVPA計算の最終年の次の年
                       ABC.year=NULL, # ABC yearを計算する年。NULLの場合はVPA計算の最終年の次の次の年
                       waa.year=NULL, # VPA結果から生物パラメータをもってきて平均する期間
                                      # NULLの場合，VPAの最終年のパラメータを持ってくる
                       maa.year=NULL, # VPA結果から生物パラメータをもってきて平均する期間
                       M.year=NULL, # VPA結果から生物パラメータをもってきて平均する期間

                       plus.group=res0$input$plus.group,
                       N=1000,# 確率的なシミュレーションをする場合の繰り返し回数。
                              # N+1の結果が返され、1列目に決定論的な結果が
                              # 0を与えると決定論的な結果のみを出力
                        silent=FALSE, is.plot=TRUE, # 計算条件を出力、プロットするか

                        pre.catch=NULL, #list(year=2012,wcatch=13000), 漁獲重量をgivenで与える場合
                        outtype="FULL", # 結果の出力を小さくするか。FULL=しない。それ以外＝する。
                       #-------- 加入に関する設定 -----------------
                       rec.new=NULL, # 指定した年の加入量
                                     # 年を指定しないで与える場合は、自動的にスタート年の加入になる。
                                     # list(year=, rec=)で与える場合は、対応する年の加入を置き換える。
                       #--- 加入関数
                       recfunc=RPS.simple.rec, # 太平洋マサバ、ゴマサバ以外はRPS.simple.recを使う
                       rec.arg=list(upper.ssb=Inf,upper.recruit=Inf), # 加入の各種設定
                       #--- Frecオプション；Frec計算のための設定リストを与えると、指定された設定でのFrecに対応するFで将来予測を行う
                       Frec=NULL, # list(stochastic=TRUE, # TRUEの場合、stochastic simulationで50%の確率でBlimitを越す(PMS, TMI)
                                                          # FALSEの場合、RPS固定のprojectionがBilmitと一致する(NSK)
                                 #      future.year=2018, # 何年の資源量を見るか？
                                 #      Blimit=450*1000,  # Blimit (xトン)
                            #      seed=100) # 乱数のシード
                       # 対馬サバに対応するオプション。ts=2のとき、1年を2つの季節に分けて将来予測する
                       ts=1, # 時間ステップ。1年間の季節の数。普通は１。対馬サバの場合2。ts=1の場合、以下の引数はすべて無視される。
                #---- 以下、ts>2のときに必要な引数 -----
                       waa=NULL,waa.catch=NULL,maa=NULL,M=NULL, # 季節毎の生物パラメータ
                       rec.season=1, # 加入がおこる季節
                       waa.multi="opt", # waa.optyearに対応する年について、暦年漁獲量と一致するようにwaaを最適化するか？ "opt"の場合、内部で最適化。waa.optyearの長さ分のベクトルを与えて指定することもできる
                       waa.optyear=2011:2013, # waa.optyearをするときに、置き換えるwaaの年
                       replace.rec.year=2012, # 加入量を暦年の将来予測での加入量に置き換えるか？
                       partial.F=NULL         # 季節毎のpartial F
                       ){

  if(!is.null(seed)) set.seed(seed)
  argname <- ls()
  arglist <- lapply(argname,function(x) eval(parse(text=x)))
  names(arglist) <- argname

  if(is.null(res0$input$unit.waa)) res0$input$unit.waa <- 1
  if(is.null(res0$input$unit.caa)) res0$input$unit.caa <- 1
  if(is.null(res0$input$unit.biom)) res0$input$unit.biom <- 1

  if(ts>1 && is.null(partial.F)){
    stop("When ts>1, partial.F should be given")
  }
  #--------------------------------------------------
  N <- N + 1
  years <- as.numeric(dimnames(res0$naa)[[2]])

  #------------- set default options
  if(is.null(currentF)) currentF <- res0$Fc.at.age
  if(is.null(waa.year)) waa.year <- rev(years)[1]
  if(is.null(maa.year)) maa.year <- rev(years)[1]
  if(is.null(M.year)) M.year <- rev(years)[1]
  if(is.null(start.year)) start.year <- rev(years)[1]+1
  if(is.null(ABC.year)) ABC.year <- rev(years)[1]+1
  arglist$ABC.year <- ABC.year
  #-------------

#  fyears <- start.year:(start.year+nyear-1)
  fyears <- seq(from=start.year,to=start.year+nyear-1,by=1/ts)
  fyear.year <- floor(fyears)
  fyear.season <- #fyears-fyear.year
                  rep(1:ts,nyear)
  fyear.season <- fyear.season[1:length(fyears)]
  ntime <- length(fyears)
#  if(is.null(multi.year)) multi.year <- rep(1,nyear)*ts
  ages <- as.numeric(dimnames(res0$naa)[[1]])
#  nage <- length(ages)
  min.age <- min(as.numeric(ages))
  if(ts>1){
    ages <- seq(from=min(ages),to=max(ages)+1/ts,by=1/ts)
    nage <- length(ages) # naaにNAが入っていて、かつ、半年毎の将来予測をする場合対応できない可能性がある
  }
  if(any(is.na(res0$naa[,ncol(res0$naa)]))){
    nage <- sum(!is.na(res0$naa[,ncol(res0$naa)])) # naaにNAが入っている対馬マイワシ対応
  }
  else{
    nage <- length(ages)
  }

  if(!silent)  cat("F multiplier= ", multi,"seed=",seed,"\n")
  #------------Frecオプションの場合 -------------
  if(!is.null(Frec)){
    if(is.null(Frec$stochastic)) Frec$stochastice <- TRUE
    if(is.null(Frec$method)) Frec$method <- "nibun"
    if(is.null(Frec$seed)) Frec$seed <- as.numeric(Sys.time())

    getFrec <- function(x,arglist){
      set.seed(Frec$seed)
      arglist.tmp <- arglist
      arglist.tmp$multi <- x
      arglist.tmp$Frec <- NULL
      arglist.tmp$is.plot <- FALSE
      if(Frec$stochastic==FALSE){
        arglist.tmp$N <- 0
      }
      fres.tmp <- do.call(future.vpa,arglist.tmp)
      tmp <- rownames(fres.tmp$vssb)==Frec$future.year
      if(all(tmp==FALSE)) stop("nyear should be longer than Frec$future.year.")
      if(Frec$stochastic==TRUE){
        is.lower.ssb <- fres.tmp$vssb<Frec$Blimit
        probs <- (sum(is.lower.ssb[tmp,],na.rm=T)-1)/
          (length(is.lower.ssb[tmp,])-1)*100
        return.obj <- probs-50
      }
      else{
        return.obj <- Frec$Blimit-fres.tmp$vssb[tmp,1]
      }
      return(ifelse(Frec$method=="nibun",return.obj,return.obj^2))
    }

    if(Frec$method=="nibun"){
      # 二分法
      eps <- ifelse(Frec$stochastic==TRUE,0.5,0.001)
      x.high <- 2 ; x.low <- 0.01;  fx <- Inf
      max.count <- 1000
      s <- 1
      while(abs(fx)>eps && s<max.count){
        x <- (x.high+x.low)/2
        fx <- getFrec(x,arglist)
        if(fx>0) x.high <- x
        if(fx<0) x.low <- x
        cat("fx =",fx,"\n")
        s <- s+1
      }
      multi <- x
    }
    else{
      # optimizeを使う場合=>収束基準が厳しいので時間がかかる
      res <- optimize(getFrec,interval=c(0.01,2),arglist=arglist)
      multi <- res$minimum
    }
  }

  #-------------- main function ---------------------
  # ts>1 (半年毎の将来予測の場合、半年毎のwaa, maa, Mを別に与える必要がある)
  if(ts>1 && ((any(sapply(list(waa,maa,M),is.null))) || (any(sapply(list(waa,maa,M),length)!=length(ages))))){
    stop("Appropriate biological paramters of waa, maa, M should be given when ts>1.")
  }
  else{
    waa.org <- waa
    waa.catch.org <- waa.catch
    maa.org <- maa
    M.org <- M
  }

  faa <- naa <- waa <- waa.catch <- maa <- M <- caa <-
          array(NA,dim=c(length(ages),ntime,N),dimnames=list(ages,fyears,1:N))
  # future biological patameter
  if(!is.null(M.org))  M[] <- M.org  else M[] <- apply(as.matrix(res0$input$dat$M[,years %in% M.year]),1,mean)
  if(!is.null(waa.org))  waa[] <- waa.org  else waa[] <- apply(as.matrix(res0$input$dat$waa[,years %in% waa.year]),1,mean)
  if(!is.null(maa.org))  maa[] <- maa.org  else maa[] <- apply(as.matrix(res0$input$dat$maa[,years %in% maa.year]),1,mean)
  if(!is.null(waa.catch.org)){
      waa.catch[] <- waa.catch.org
  }
  else{
      if(!is.null(res0$input$dat$waa.catch)) waa.catch[] <- apply(as.matrix(res0$input$dat$waa.catch[,years %in% waa.year]),1,mean)
      else waa.catch <- waa
  }

  # time step definition (change F and M)
  M <- M/ts
  if(ts>1){
    currentF <- as.numeric(sweep(matrix(partial.F,ts,nage/ts),2,currentF,FUN="*"))
  }

  # future F matrix
  faa[] <- currentF*multi
  faa[,fyears<min(ABC.year),] <- currentF
#  browser()
  if(length(tmp <- which(fyear.year %in% years))>0){
      tmp0 <- which(years %in% fyear.year)
#      tmp1 <- which(fyear.year %in% years)
      for(jj in 1:length(tmp0)){
        for(j in 1:length(tmp)){
          if(ts>1){
            # VPAデータを季節で展開
            faa[,tmp[j],] <-
              as.numeric(sweep(matrix(partial.F,ts,nage/ts),2,res0$faa[,tmp0[jj]],FUN="*"))
          }
          else{
            if(any(res0$faa[,tmp0[jj]]>0)){ # もしfaaがゼロでないなら（PMIの場合、2012までデータが入っているが、faaはゼロになっているので
              faa[,tmp[j],] <- res0$faa[,tmp0[jj]]
              waa[,tmp[j],] <- res0$input$dat$waa[,tmp0[jj]]
              if(!is.null(res0$input$dat$waa.catch)){
                  waa.catch[,tmp[j],] <- res0$input$dat$waa.catch[,tmp0[jj]]
              }
              else{
                  waa.catch[,tmp[j],] <- res0$input$dat$waa[,tmp0[jj]]
                  }
            }
          }
        }}}

  if(ts>1){
    for(j in 1:ts){
      for(kk in 1:N){
#          faa[max(floor(ages))==floor(ages),fyear.season==j,][j,,kk]
        # plus goupのFやwaaは季節によって変えないといけない
        # (plus groupに限らない？。1年に複数回の加入がある場合、季節によるFの違いなのか、加入群に対するFの違いなのかによって仕様を変える必要がある)
        tmp <- t(faa[max(floor(ages))==floor(ages),fyear.season==j,kk])
        tmp[] <- faa[max(floor(ages))==floor(ages),fyear.season==j,,drop=F][j,,kk]
        faa[max(floor(ages))==floor(ages),fyear.season==j,kk] <- t(tmp)

        tmp <- t(waa[max(floor(ages))==floor(ages),fyear.season==j,kk])
        tmp[] <- waa[max(floor(ages))==floor(ages),fyear.season==j,,drop=F][j,,kk]
        waa[max(floor(ages))==floor(ages),fyear.season==j,kk] <- t(tmp)
      }
    }

    #waaは歴年の漁獲量と同じになるように最適化する

    arglist.tmp <- arglist
    arglist.tmp$ts <- 1
    arglist.tmp$N <- 0
    arglist.tmp$silent <- TRUE
    arglist.tmp$is.plot <- FALSE
    arglist.tmp$waa <- arglist.tmp$maa <- arglist.tmp$M <- NULL
    # SSB用
    fres.cyear <- do.call(future.vpa,arglist.tmp)
    # waaの補正用
    arglist.tmp2 <- arglist.tmp
    arglist.tmp2$multi <- 1
    a <- do.call(future.vpa,arglist.tmp2)
    if(!is.numeric(waa.multi)){ # if waa="opt"
      optfunc <- function(x,arglist,a,waa.optyear,replace.rec.year){
        opt.catch <- a$vwcaa[names(a$vwcaa[,1])%in%waa.optyear]
        arglist.tmp <- arglist
        arglist.tmp$N <- 0
        arglist.tmp$silent <- TRUE
        arglist.tmp$is.plot <- FALSE
        arglist.tmp$waa.multi <- x
          #        browser()
        arglist.tmp$rec.new <- list(year=replace.rec.year,rec=a$naa[1,a$year==replace.rec.year,1])
#        cat(arglist.tmp$rec.new$rec,"\n")
        a.tmp <- do.call(future.vpa,arglist.tmp)
        pre.catch <- tapply(a.tmp$vwcaa[,1],a.tmp$fyear.year,sum)

        xx <- sum((pre.catch[names(pre.catch)%in%waa.optyear]-opt.catch)^2)
#        cat(xx,"\n")
        return(xx)
      }
#      browser()
#      tmp <- optfunc(c(1,1,1),arglist=arglist,opt.catch=opt.catch,waa.optyear=waa.optyear)
#      debug(future.vpa)
      est <- optim(rep(1,length(waa.optyear)),optfunc,
                   arglist=arglist,a=a,waa.optyear=waa.optyear,replace.rec.year=replace.rec.year)
      waa.multi <- est$par
      cat(waa.multi,"\n")
      rec.new <- list(year=replace.rec.year,rec=a$naa[1,a$year==replace.rec.year,1])
    }
    for(kk in 1:length(waa.optyear)){
      waa[,fyear.year==waa.optyear[kk],] <- waa[,fyear.year==waa.optyear[kk],] * waa.multi[kk]
    }
  }
  tmp <- aperm(faa,c(2,1,3))
  tmp <- tmp*multi.year
  faa <- aperm(tmp,c(2,1,3))

  #  vpa.multi <- ifelse(is.null(vpa.mode),1,vpa.mode$multi)
  # rps assumption
  rps.mat <- array(NA,dim=c(ntime,N),dimnames=list(fyears,1:N))
#  rps.mat[] <- sample(rps.range2,nyear*N,replace=TRUE)  # 平均を揃えたもの
#  rps.mat[,1] <- rps.med
  rec.tmp <- list(rec.resample=NULL,tmparg=NULL)

  if(!is.null(Frec$seed)) set.seed(Frec$seed)

#  for(k in 1:N){  #k loopを消す
    # future N matrix
    if(sum(start.year==years)==0){
      # VPA結果が2011年まで、将来予測が2012年からだったら、VPA結果を使って2011年まで1年前進計算を行う
      if(start.year==(max(years)+1)){
#        tmp <- forward.calc(res0$faa,res0$naa,res0$input$dat$M,
#                            rep(nage,length(years)+1),length(years)+1)
        tmp <- forward.calc.simple(res0$faa[,length(years)],
                                     res0$naa[,length(years)],
                                     res0$input$dat$M[,length(years)],
                                   plus.group=plus.group)
        if(ts==1){
          naa[1:nage,1,] <- tmp
        }
        else{
          naa[1:nage,1,] <- 0
          naa[(ages-floor(ages))==0,1,] <- tmp
        }
        if(all(is.na(naa[1,1,]))){
          if(fyears[1]-min.age < start.year){
            thisyear.ssb <- rep(sum(res0$ssb[,as.character(fyears[1]-min.age)],na.rm=T),N)
          }
          else{
              thisyear.ssb <- colSums(naa[,1,]*waa[,1,]*maa[,1,],na.rm=T)*res0$input$unit.waa/res0$input$unit.biom
          }
          rec.tmp0 <- recfunc(thisyear.ssb[1],res0,
                             rec.resample=rec.tmp$rec.resample,
                             rec.arg=rec.arg,
                             deterministic=TRUE)
          rec.tmp1 <- recfunc(thisyear.ssb[-1],res0,
                             rec.resample=rec.tmp$rec.resample,
                             rec.arg=rec.arg,
                             deterministic=FALSE)
          if(!is.null(rec.tmp$rec.arg)) rec.arg <- rec.tmp$rec.arg
          naa[1,1,] <- c(rec.tmp0$rec,rec.tmp1$rec)
          rps.mat[1,] <- naa[1,1,]/thisyear.ssb
        }
      }
      else{
        stop("ERROR Set appropriate year to start projection\n")
      }
    }
    else{
      if(any(ts==rec.season)){
        naa[,1,] <- res0$naa[,start.year==years]
      }
      else{
        naa[,1,] <- 0
        naa[(ages-floor(ages))==0,1,] <- res0$naa[,start.year==years]
      }
    }

    if(!is.null(rec.new)){
      if(!is.list(rec.new)){
        naa[1,1,] <- rec.new
      }
      else{ # rec.newがlistの場合
        naa[1,fyears==rec.new$year,] <- rec.new$rec
      }}

    for(i in 1:(ntime-1)){
      if(Pope){
       caa[,i,] <- naa[,i,]*(1-exp(-faa[,i,]))*exp(-M[,i,]/2)
     }
      else{
        caa[,i,] <- naa[,i,]*(1-exp(-faa[,i,]-M[,i,]))*faa[,i,]/(faa[,i,]+M[,i,])
      }

      #漁獲量がgivenの場合
      if(!is.null(pre.catch) && fyears[i]==pre.catch$year){
        for(k in 1:N){
          tmp <- caa.est(naa[,i,k],faa[,i,k]/max(faa[,i,k]),
                         waa.catch[,i,k],M[,i,k],pre.catch$wcatch*1000,Pope=Pope)
          faa.new <- tmp$x * faa[,i,k]/max(faa[,i,k])
          caa[,i,k] <- tmp$caa
          faa[,i,k] <- faa.new
        }}

      tmp <- forward.calc.mat(faa[,i,],naa[,i,],M[,i,],plus.group=plus.group)
      naa[,i+1,][is.na(naa[,i+1,])] <- tmp[is.na(naa[,i+1,])]

      # 当年の加入の定義
#      if(ifelse(is.null(vpa.mode),TRUE, sum(years==fyears[i+1])==0|vpa.mode$rec=="recfun")){
      if(fyears[i+1]-min.age < start.year){
        thisyear.ssb <- rep(sum(res0$ssb[,as.character(fyears[i+1]-min.age)],na.rm=T),N)
      }
      else{
        if(ts==1){
          thisyear.ssb <- colSums(naa[,i+1-min.age,]*waa[,i+1-min.age,]*
                              maa[,i+1-min.age,],na.rm=T)*res0$input$unit.waa/res0$input$unit.biom
        }
        else{
          # 暦年の将来予測での親魚資源量から加入量を推定する（もしかしたらreplace.rec.yearはこれがあれば必要ないのかもしれない）
          # min.ageが0才より大きく、半年毎の計算をする場合対応できない
           # stochasticのときはどうする？
          cssb <- fres.cyear$vssb[,1]
          thisyear.ssb <- cssb[as.numeric(names(cssb))==fyears[i+1]]
          if(length(thisyear.ssb)==0) thisyear.ssb <- 0
        }
      }
      rec.tmp0 <- recfunc(thisyear.ssb[1],res0,
                         rec.resample=rec.tmp$rec.resample,
                         rec.arg=rec.arg,
                         deterministic=TRUE)
      rec.tmp <- recfunc(thisyear.ssb[-1],res0,
                         rec.resample=rec.tmp$rec.resample,
                         rec.arg=rec.arg,
                         deterministic=FALSE)
      if(!is.null(rec.tmp$rec.arg)) rec.arg <- rec.tmp$rec.arg
      if(all(is.na(naa[1,i+1,])) && (floor(fyears[i+1])-fyears[i+1])==0){ # 加入は最初の季節にのみおこる
        naa[1,i+1,] <- c(rec.tmp0$rec,rec.tmp$rec)
      }
      else{
        if(all(is.na(naa[1,i+1,]))) naa[1,i+1,] <- 0
      }
      rps.mat[i+1,] <- naa[1,i+1,]/thisyear.ssb
    }
    if(Pope){
      caa[,ntime,] <- naa[,ntime,]*(1-exp(-faa[,ntime,]))*exp(-M[,ntime,]/2)
    }
    else{
      caa[,ntime,] <- naa[,ntime,]*(1-exp(-faa[,ntime,]-M[,ntime,]))*faa[,ntime,]/(faa[,ntime,]+M[,ntime,])
    }

  biom <- naa*waa*res0$input$unit.waa/res0$input$unit.biom
  ssb <- naa*waa*maa*res0$input$unit.waa/res0$input$unit.biom

  wcaa <- caa*waa.catch*res0$input$unit.waa/res0$input$unit.biom
  vwcaa <- apply(wcaa,c(2,3),sum,na.rm=T)

  ABC <- apply(as.matrix(vwcaa[fyears%in%ABC.year,,drop=F]),2,sum)

  if(outtype=="FULL"){
      fres <- list(faa=faa,naa=naa,biom=biom,ssb=ssb,wcaa=wcaa,caa=caa,M=M,rps=rps.mat,
                   maa=maa,vbiom=apply(biom,c(2,3),sum,na.rm=T),
                   waa=waa,waa.catch=waa.catch,currentF=currentF,
                   vssb=apply(ssb,c(2,3),sum,na.rm=T),vwcaa=vwcaa,
                   years=fyears,fyear.year=fyear.year,ABC=ABC,recfunc=recfunc,rec.arg=rec.arg,
                   waa.year=waa.year,maa.year=maa.year,multi=multi,multi.year=multi.year,
                   Frec=Frec,rec.new=rec.new,pre.catch=pre.catch,input=arglist)
  }
  else{
      fres <- list(faa=faa[,,1],M=M[,,1],recruit=naa[1,,],
                   maa=maa[,,1],vbiom=apply(biom,c(2,3),sum,na.rm=T),
                   waa=waa[,,1],waa.catch=waa.catch[,,1],currentF=currentF,
                   vssb=apply(ssb,c(2,3),sum,na.rm=T),vwcaa=vwcaa,
                   years=fyears,fyear.year=fyear.year,ABC=ABC,recfunc=recfunc,rec.arg=rec.arg,
                   waa.year=waa.year,maa.year=maa.year,multi=multi,multi.year=multi.year,
                   Frec=Frec,rec.new=rec.new,pre.catch=pre.catch,input=arglist)
  }
  class(fres) <- "future"
  if(is.plot){
    par(mfrow=c(2,2))
    plot.future(fres)
  }
  invisible(fres)
}



#----------------------------------------------------------------------
#----------   加入に関する関数。魚種specific        -------------------
#----------------------------------------------------------------------

#-------------- VPA mode 用関数 -------------------
caa.est <- function(naa,saa,waa,M,catch.obs,Pope){
  saa <- saa/max(saa)
  tmpfunc <- function(x,catch.obs=catch.obs,naa=naa,saa=saa,waa=waa,M=M,out=FALSE,Pope=Pope){
    if(isTRUE(Pope)){
      caa <- naa*(1-exp(-saa*x))*exp(-M/2)
    }
    else{
      caa <- naa*(1-exp(-saa*x-M))*saa*x/(saa*x+M)
    }
    wcaa <- caa*waa
    if(out==FALSE){
      return((sum(wcaa,na.rm=T)-catch.obs)^2)
    }
    else{
      return(caa)
    }
  }
  tmp <- optimize(tmpfunc,c(0,5),catch.obs=catch.obs,naa=naa,saa=saa,waa=waa,M=M,Pope=Pope,out=FALSE)
  tmp2 <- tmpfunc(x=tmp$minimum,catch.obs=catch.obs,naa=naa,saa=saa,waa=waa,M=M,Pope=Pope,out=TRUE)
  return(list(x=tmp$minimum,caa=tmp2))
}


#type="TorF" # true or false
#type="diff" # excel-RVPA
#type="%" # (excel-RVPA)/excel
check.res <- function(res,fres,tdata,digits=3,type="%"){

  check.twomats <- function(mat1,mat2,digits=3,type="%"){
    if(!is.null(colnames(mat1))){
      tmp1 <- mat1[,colnames(mat1)%in%colnames(mat2)]
      tmp2 <- mat2[,colnames(mat2)%in%colnames(mat1)]
    }
    else{
      tmp1 <- mat1
      tmp2 <- mat2
    }
    if(type=="TorF"){
      tmp <- round(tmp1,digits) == round(tmp2,digits)
    }
    if(type=="diff"){
      tmp <- round(tmp1-tmp2,digits)
    }
    if(type=="%"){
      tmp <- round((tmp1-tmp2)/tmp1*100,digits)
    }
    return(tmp)
  }

  naa.res <- check.twomats(tdata$naa,res$naa,digits=digits,type=type)
  faa.res <- check.twomats(tdata$faa,res$faa,digits=digits,type=type)
  fcaa.res <- check.twomats(tdata$Fc.at.age,res$Fc.at.age,digits=digits,type=type)

  tmp.list <- list(naa=naa.res,faa=faa.res,Fc.at.age=fcaa.res)
  return(tmp.list)
}


#---------------- 結果の確かめ用関数 ---------------------
# --------USAGE-------
# tdata <- get.tdata("vpa_results.csv")
# check.res(res.pms,list(fres,fres),tdata,digits=2,type="%")

get.data <- function(tfile){
  tmpdata <- read.csv(tfile,header=F,as.is=F,colClasses="character")
  flags <- which(substr(tmpdata[,1],1,1)=="#")
  tlist <- list()
  for(i in 1:(length(flags)-1)){
      tmp <- tmpdata[(flags[i]+1):(flags[i+1]-1),]
      if(dim(tmp)[[1]]>1){
        dimnames(tmp) <- list(tmp[,1],tmp[1,])
        tmp <- tmp[,!apply(tmp=="",2,all)]
        tlist[[i]] <- sapply((tmp[-1,-1]),as.numeric)
      }
     else{
        tlist[[i]] <- as.numeric(tmp[tmp!=""])
      }
  }
  names(tlist)[1:4] <- c("naa","faa","Biomass","Fc.at.age")
  dimnames(tlist[[3]])[[1]] <- c("SSB","Biomass")
  for(i in 1:tlist[[5]]){
    names(tlist)[(4+(i-1)*4+1):(4+(i*4))] <- c("fnaa","ffaa","fwcaa","ABC")
  }
  return(tlist)
}
