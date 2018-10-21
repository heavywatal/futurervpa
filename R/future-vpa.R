#' 将来予測関数
#'
#' @description
#' multiのオプションは管理後のFのmultiplier（管理前後でselectivityが同じ）
#' 年数がNULLの場合，VPAの最終年のパラメータを持ってくる。
#' @param currentF 管理前のF
#' @param multi 管理後（ABC.yearから）のF (current F x multi)
#' @param multi.year ある特定の年だけFを変えたい場合。デフォルトは1。変える場合は、指定した年またはタイムステップの要素数のベクトルで指定。
#' @param start.year 将来予測の開始年，NULLの場合はVPA計算の最終年の次の年
#' @param ABC.year ABC yearを計算する年。NULLの場合はVPA計算の最終年の次の次の年
#' @param waa.year VPA結果から生物パラメータをもってきて平均する期間
#' @param maa.year VPA結果から生物パラメータをもってきて平均する期間
#' @param M.year VPA結果から生物パラメータをもってきて平均する期間
#' @param strategy F: 漁獲係数一定, E: 漁獲割合一定、C: 漁獲量一定（pre.catchで漁獲量を指定）
#' @param HCR HCRを使う場合、list(Blim=154500, Bban=49400,beta=1)のように指定するか、以下の引数をセットする,
#' @param N 確率的なシミュレーションをする場合の繰り返し回数。
#'        N+1の結果が返され、1列目に決定論的な結果が
#'        0を与えると決定論的な結果のみを出力
#' @param silent,is.plot 計算条件を出力、プロットするか
#' @param random.select 選択率をランダムリサンプリングする場合、ランダムリサンプリングする年を入れる
#'        strategy="C"または"E"のときのみ有効
#' @param pre.catch list(year=2012,wcatch=13000), 漁獲重量をgivenで与える場合
#'        list(year=2012:2017,E=rep(0.5,6)), 漁獲割合をgivenで与える場合
#' @param rec.new 指定した年の加入量
#'        年を指定しないで与える場合は、自動的にスタート年の加入になる。
#'        list(year=, rec=)で与える場合は、対応する年の加入を置き換える。
#' @param recfunc 再生産関係の関数
#' @param rec.arg 加入関数の各種設定
#' @param Frec Frec計算のための設定リストを与えると、指定された設定でのFrecに対応するFで将来予測を行う
#'        list(stochastic=TRUE, # TRUEの場合、stochastic simulationで50%の確率でBlimitを越す(PMS, TMI)
#'        FALSEの場合、RPS固定のprojectionがBilmitと一致する(NSK)
#'        future.year=2018, # 何年の資源量を見るか？
#'        Blimit=450*1000,  # Blimit (xトン)
#'        scenario="catch.mean" or "blimit" (デフォルトはblimit; "catch.mean"とするとstochastic simulationにおける平均漁獲量がBlimitで指定した値と一致するようになる)
#'        Frange=c(0.01,2*mult)) # Fの探索範囲
#' @param waa,waa.catch,maa,M 季節毎の生物パラメータ、または、生物パラメータを外から与える場合
#' @param waa.multi waa.optyearに対応する年について、暦年漁獲量と一致するようにwaaを最適化するか？ "opt"の場合、内部で最適化。waa.optyearの長さ分のベクトルを与えて指定することもできる（ts>1のときのオプション）
#' @param waa.optyear waa.optyearをするときに、置き換えるwaaの年
#' @param replace.rec.year 加入量を暦年の将来予測での加入量に置き換えるか？
#' @param waa.fun waaをnaaのfunctionとするか
#' @param add.year 岡村オプションに対応。=1で1年分余計に計算する
#' @param det.run 1回めのランは決定論的将来予測をする（完璧には対応していない）
#' @rdname future-vpa
#' @export
future.vpa <- function(
    res0,
    currentF=NULL,
    multi=1,
    nyear=10,Pope=res0$input$Pope,
    outtype="FULL",
    multi.year=1,
    start.year=NULL,
    ABC.year=NULL,
    waa.year=NULL,
    maa.year=NULL,
    M.year=NULL,
    seed=NULL,
    strategy="F",
    HCR=NULL,
    beta=NULL,delta=NULL,Blim=0,Bban=0,
    plus.group=res0$input$plus.group,
    N=1000,
    silent=FALSE, is.plot=TRUE,
    random.select=NULL,
    pre.catch=NULL,
    rec.new=NULL,
    recfunc=HS.recAR,
    rec.arg=list(upper.ssb=Inf,upper.recruit=Inf),
    Frec=NULL,
    waa=NULL,waa.catch=NULL,maa=NULL,M=NULL,
    waa.multi="opt",
    waa.optyear=2011:2013,
    replace.rec.year=2012,
    F.sigma=0,
    waa.fun=FALSE,
    naa0=NULL,eaa0=NULL,ssb0=NULL,faa0=NULL,
    add.year=0,
    det.run=TRUE
  ){

      argname <- ls()
      arglist <- lapply(argname,function(x) eval(parse(text=x)))
      names(arglist) <- argname

      if(is.null(res0$input$unit.waa)) res0$input$unit.waa <- 1
      if(is.null(res0$input$unit.caa)) res0$input$unit.caa <- 1
      if(is.null(res0$input$unit.biom)) res0$input$unit.biom <- 1
      if(is.null(plus.group)) plus.group <- TRUE
      if(is.null(Pope)) Pope <- FALSE

      ##--------------------------------------------------
      if(isTRUE(det.run)) N <- N + 1
      years <- as.numeric(dimnames(res0$naa)[[2]])

      ##------------- set default options
      if(is.null(currentF)) currentF <- res0$Fc.at.age
      if(is.null(waa.year)) waa.year <- rev(years)[1]
      if(is.null(maa.year)) maa.year <- rev(years)[1]
      if(is.null(M.year)) M.year <- rev(years)[1]
      if(is.null(start.year)) start.year <- rev(years)[1]+1
      if(is.null(ABC.year)) ABC.year <- rev(years)[1]+1
      ##    if(!is.null(Bban)) Bban$is.Bban <- rep(FALSE,N)
      arglist$ABC.year <- ABC.year
      ##-------------

      ##---- set S-R functin option -----
      ## 使う関数によっては必要ないオプションもあるが、使わないオプションを入れてもエラーは出ないので、
      # rec.arg$resampleがNULLかどうかで、パラメトリックな誤差分布かそうでないか（残差リサンプリング）を判別する
      if(is.null(rec.arg$sd2)) rec.arg$sd2 <- sqrt(rec.arg$sd^2/(1-rec.arg$rho^2)) #rho込み平均補正用SD # HS.recAR

      if(is.null(rec.arg$resample)|!isTRUE(rec.arg$resample)){
          if(is.null(rec.arg$bias.correction)) rec.arg$bias.correction <- TRUE # HS.recAR, HS.rec0
#          if(is.null(rec.arg$resample)) rec.arg$resample <- FALSE      # HS.rec0オプション
          if(is.null(rec.arg$rho)){
              rec.arg$rho <- 0 # HS.recAR, HS.rec0
              rec.arg$resid <- 0
          }
          if(!is.null(rec.arg$rho)){
              if(rec.arg$rho>0){
                  if(is.null(eaa0)) rec.arg$resid <- rep(rev(rec.arg$resid)[1],N)
                  else{ rec.arg$resid <- eaa0 }
              }
              else{
                  rec.arg$resid <- rep(0,N)
              }
          }
      }
      else{
         rec.arg$rho <- 0 # resamplingの場合に自己相関は考慮できないのでrhoは強制的にゼロ
      }

      if(!is.null(rec.arg$sd)) rec.arg$sd <- c(0,rep(rec.arg$sd,N-1))
      if(!is.null(rec.arg$sd2)) rec.arg$sd2 <- c(0,rep(rec.arg$sd2,N-1))
      ##---------------------------------

      if(!is.null(beta)){
          HCR$beta <- beta
          HCR$Blim <- Blim
          HCR$Bban <- Bban
      }

    #  fyears <- seq(from=start.year,to=start.year+nyear-1,by=1/ts)
    fyears <- seq(from=start.year,to=start.year+nyear+add.year,by=1)

    fyear.year <- floor(fyears)
    ntime <- length(fyears)
    ages <- as.numeric(dimnames(res0$naa)[[1]])
    min.age <- min(as.numeric(ages))
#    if(ts>1){
#      ages <- seq(from=min(ages),to=max(ages)+1/ts,by=1/ts)
#      nage <- length(ages) # naaにNAが入っていて、かつ、半年毎の将来予測をする場合対応できない可能性がある
#    }
    if(any(is.na(res0$naa[,ncol(res0$naa)]))){
      nage <- sum(!is.na(res0$naa[,ncol(res0$naa)])) # naaにNAが入っている対馬マイワシ対応
    }
    else{
      nage <- length(ages)
    }

    if(!silent)  cat("F multiplier= ", multi,"seed=",seed,"\n")

    # シードの設定
    if(is.null(seed)) arglist$seed <- as.numeric(Sys.time())

    #------------Frecオプションの場合 -------------
    if(!is.null(Frec)){
      multi.org <- multi
      if(is.null(Frec$stochastic)) Frec$stochastice <- TRUE
#      if(is.null(Frec$method)) Frec$method <- "optimize"
      if(is.null(Frec$target.probs)) Frec$target.probs <- 50
      if(is.null(Frec$scenario)) Frec$scenario <- "blimit" # 2017/12/25追記
      if(is.null(Frec$Frange)) Frec$Frange <- c(0.01,multi.org*2)   # 2017/12/25追記(探索するFの範囲の指定)
      #      arglist$Frec <- Frec

      getFrec <- function(x,arglist){
        set.seed(arglist$seed)
        arglist.tmp <- arglist
        arglist.tmp$multi <- x
        arglist.tmp$silent <- TRUE
        arglist.tmp$Frec <- NULL
        arglist.tmp$is.plot <- FALSE
        if(Frec$stochastic==FALSE){
          arglist.tmp$N <- 0
        }
        fres.tmp <- do.call(future.vpa,arglist.tmp)
        tmp <- rownames(fres.tmp$vssb)==Frec$future.year
        if(all(tmp==FALSE)) stop("nyear should be longer than Frec$future.year.")
        if(Frec$stochastic==TRUE){
          if(Frec$scenario=="blimit"){
            is.lower.ssb <- fres.tmp$vssb<Frec$Blimit
            probs <- (sum(is.lower.ssb[tmp,-1],na.rm=T)-1)/
              (length(is.lower.ssb[tmp,-1])-1)*100
            return.obj <- probs-Frec$target.probs
          }
          # stochastic projectionにおける平均漁獲量を目的の値に一致させる
          if(Frec$scenario=="catch.mean"){
            return.obj <- (log(Frec$Blimit)-log(mean(fres.tmp$vwcaa[tmp,-1])))^2
          }
          # stochastic projectionにおける平均親魚資源量を目的の値に一致させる
          if(Frec$scenario=="ssb.mean"){
            return.obj <- (log(Frec$Blimit)-log(mean(fres.tmp$vssb[tmp,-1])))^2
          }
        }
        else{
          return.obj <- Frec$Blimit-fres.tmp$vssb[tmp,1]
        }
#        return(ifelse(Frec$method=="nibun",return.obj,return.obj^2))
        return(return.obj^2)
      }

      res <- optimize(getFrec,interval=Frec$Frange,arglist=arglist)
      multi <- res$minimum
      cat("F multiplier=",multi,"\n")
    }

    #-------------- main function ---------------------
    waa.org <- waa
    waa.catch.org <- waa.catch
    maa.org <- maa
    M.org <- M

    if(strategy=="C"|strategy=="E") multi.catch <- multi else multi.catch <- 1

    faa <- naa <- waa <- waa.catch <- maa <- M <- caa <-
      array(NA,dim=c(length(ages),ntime,N),dimnames=list(ages,fyears,1:N))
    # future biological patameter
    if(!is.null(M.org))  M[] <- M.org  else M[] <- apply(as.matrix(res0$input$dat$M[,years %in% M.year]),1,mean)
    if(!is.null(waa.org))  waa[] <- waa.org  else waa[] <- apply(as.matrix(res0$input$dat$waa[,years %in% waa.year]),1,mean)
    if(!is.null(maa.org))  maa[] <- maa.org  else maa[] <- apply(as.matrix(res0$input$dat$maa[,years %in% maa.year]),1,mean)
    if(!is.null(waa.catch.org))  waa.catch[] <- waa.catch.org
    else{
      if(!is.null(res0$input$dat$waa.catch)) waa.catch[] <- apply(as.matrix(res0$input$dat$waa.catch[,years %in% waa.year]),1,mean)
      else waa.catch <- waa
    }


    # future F matrix
    faa[] <- currentF*multi # *exp(rnorm(length(faa),0,F.sigma))
    # ABCyear以前はcurrent Fを使う。
    faa[,fyears<min(ABC.year),] <- currentF*exp(rnorm(length(faa[,fyears<min(ABC.year),]),0,F.sigma))

    ## VPA期間と将来予測期間が被っている場合、VPA期間のFはVPAの結果を使う
    if(length(tmp <- which(fyear.year %in% years))>0){
      tmp0 <- which(years %in% fyear.year)
      for(jj in 1:length(tmp0)){
        for(j in 1:length(tmp)){
            if(any(res0$faa[,tmp0[jj]]>0) && !is.null(res0$input$dat$waa[,tmp0[jj]])){ # もしfaaがゼロでないなら（PMIの場合、2012までデータが入っているが、faaはゼロになっているので
              faa[,tmp[j],] <- res0$faa[,tmp0[jj]]
              waa[,tmp[j],] <- res0$input$dat$waa[,tmp0[jj]]
              if(!is.null(res0$input$dat$waa.catch)){
                waa.catch[,tmp[j],] <- res0$input$dat$waa.catch[,tmp0[jj]]
              }
              else{
                waa.catch[,tmp[j],] <- res0$input$dat$waa[,tmp0[jj]]
              }
            }
        }}}

    tmp <- aperm(faa,c(2,1,3))
    tmp <- tmp*multi.year
    faa <- aperm(tmp,c(2,1,3))

    #  vpa.multi <- ifelse(is.null(vpa.mode),1,vpa.mode$multi)
    # rps assumption
      rps.mat <- array(NA,dim=c(ntime,N),dimnames=list(fyears,1:N))
      eaa <- matrix(0,ntime,N)
      rec.tmp <- list(rec.resample=NULL,tmparg=NULL)

    if (waa.fun){ #年齢別体重の予測関数
      WAA <- res0$input$dat$waa
      NAA <- res0$naa
#      nage <- nrow(WAA)
      WAA.res <- lapply(1:nage, function(i) {
        log.w <- as.numeric(log(WAA[i,]))
        log.n <- as.numeric(log(NAA[i,]))
        lm(log.w~log.n)
      })
      WAA.cv <- sapply(1:nage, function(i) sqrt(mean(WAA.res[[i]]$residuals^2)))
      WAA.b0 <- sapply(1:nage, function(i) as.numeric(WAA.res[[i]]$coef[1]))
      WAA.b1 <- sapply(1:nage, function(i) as.numeric(WAA.res[[i]]$coef[2]))
      ##      waa.rand <- array(0,dim=c(al,nyear+1-min.age,N))
      set.seed(0)
      cv.vec <- rep(WAA.cv,N*ntime)
      waa.rand <- array(rnorm(length(cv.vec),-0.5*cv.vec^2,cv.vec),dim=c(nage,ntime,N))
      waa.rand[,,1] <- 0
#      for (ii in 1:N) {
#        if (ii==1) {
#          waa.rand[,,ii] <- t(sapply(1:nage, function (j) rnorm(ntime,0,0)))
#        } else {
#          waa.rand[,,ii] <- t(sapply(1:nage, function (j) rnorm(ntime,-0.5*WAA.cv[j]^2,WAA.cv[j])))
#        }
                                       #      }
    }

      set.seed(arglist$seed)

      # 1年目の年齢組成を入れる
      if(!start.year%in%years){
        # VPA結果が2011年まで、将来予測が2012年の場合
        if(start.year==(max(years)+1)){
            {if(is.null(res0$input$dat$M)){
                M.lastyear <- M.org
            }
            else{
                M.lastyear <- res0$input$dat$M[,length(years)]
            }}
            tmp <- forward.calc.simple(res0$faa[,length(years)],
                                     res0$naa[,length(years)],
#                                     res0$input$dat$M[,length(years)],
                                     M.lastyear,
                                     plus.group=plus.group)
            naa[1:nage,1,] <- tmp
#            if(is.na(naa[1,1,])){
            if(fyears[1]-min.age < start.year){
                thisyear.ssb <- sum(res0$ssb[,as.character(fyears[1]-min.age)],na.rm=T)
            }
            else{
                if(waa.fun){
                    waa[2:nage,1,] <- t(sapply(2:nage, function(ii) as.numeric(exp(WAA.b0[ii]+WAA.b1[ii]*log(naa[ii,1,])+waa.rand[ii,1,]))))
                }
                thisyear.ssb <- colSums(naa[,1,]*waa[,1,]*maa[,1,],na.rm=T)*res0$input$unit.waa/res0$input$unit.biom                           }

            thisyear.ssb <- thisyear.ssb+(1e-10)
            rec.tmp <- recfunc(thisyear.ssb,res0,
                               rec.resample=rec.tmp$rec.resample,
                               rec.arg=rec.arg)
            eaa[1,] <- rec.tmp$rec.resample[1:N]
            rec.arg$resid <- rec.tmp$rec.resample # ARオプションに対応

            if(!is.null(rec.tmp$rec.arg)) rec.arg <- rec.tmp$rec.arg
            naa[1,1,] <- rec.tmp$rec
            if (waa.fun) {
              waa[1,1,] <- as.numeric(exp(WAA.b0[1]+WAA.b1[1]*log(naa[1,1,])+waa.rand[1,1,]))
            }
            rps.mat[1,] <- naa[1,1,]/thisyear.ssb
        }
        else{
          stop("ERROR Set appropriate year to start projection\n")
        }
      }
      else{
          naa[,1,] <- res0$naa[,start.year==years]
      }

      ### 任意の年齢組成からスタートさせたい場合###
      if(!is.null(naa0)){
          naa[,1,] <- naa0
          if(is.null(faa0)) faa0 <- res0$Fc.at.age
          faa[] <- faa0*multi
      }

      if(!is.null(rec.new)){
        if(!is.list(rec.new)){
          naa[1,1,] <- rec.new
        }
        else{ # rec.newがlistの場合
          naa[1,fyears%in%rec.new$year,] <- rec.new$rec
        }}

      for(i in 1:(ntime-1)){

        #漁獲量がgivenの場合
        if(!is.null(pre.catch) && fyears[i]%in%pre.catch$year){
          if(!is.null(pre.catch$wcatch)){
            if(fyears[i]<ABC.year){
              tmpcatch <- as.numeric(pre.catch$wcatch[pre.catch$year==fyears[i]])
            }
            else{
              tmpcatch <- as.numeric(pre.catch$wcatch[pre.catch$year==fyears[i]]) * multi.catch
            }
          }
          if(!is.null(pre.catch$E)){
            biom <- sum(naa[,i,]*waa[,i,]*res0$input$unit.waa/res0$input$unit.biom)
            if(fyears[i]<ABC.year){
              tmpcatch <- as.numeric(pre.catch$E[pre.catch$year==fyears[i]])  * biom
            }
            else{
              tmpcatch <- as.numeric(pre.catch$E[pre.catch$year==fyears[i]]) * biom * multi.catch
            }
          }

          # 選択率をランダムサンプリングする場合
#          if(!is.null(random.select)) saa.tmp <- as.numeric(res0$saa[,colnames(res0$saa)==sample(random.select,1)])
          saa.tmp <- sweep(faa[,i,],2,apply(faa[,i,],2,max),FUN="/")
          tmp <- lapply(1:dim(naa)[[3]],function(x) caa.est.mat(naa[,i,x],saa.tmp[,x],
                                                                waa.catch[,i,x],M[,i,x],tmpcatch,Pope=Pope))
          faa.new <- sapply(tmp,function(x) x$x) * saa.tmp
          caa[,i,] <- sapply(tmp,function(x) x$caa)
          faa[,i,] <- faa.new
        }

          ## HCRを使う場合(当年の資源量から当年のFを変更する)
          if(!is.null(HCR)&&fyears[i]>=ABC.year){
              ssb.tmp <- colSums(naa[,i,]*waa[,i,]*maa[,i,],na.rm=T)*
                                               res0$input$unit.waa/res0$input$unit.biom
              alpha <- ifelse(ssb.tmp<HCR$Blim,HCR$beta*(ssb.tmp-HCR$Bban)/(HCR$Blim-HCR$Bban),HCR$beta)
              faa[,i,] <- sweep(faa[,i,],2,alpha,FUN="*")
              faa[,i,] <- ifelse(faa[,i,]<0,0,faa[,i,])
          }

          ## 漁獲して１年分前進（加入はまだいれていない）
          tmp <- forward.calc.mat2(faa[,i,],naa[,i,],M[,i,],plus.group=plus.group)
          # 既に値が入っているところ（１年目の加入量）は除いて翌年のNAAを入れる
          naa.tmp <- naa[,i+1,]
          naa.tmp[is.na(naa.tmp)] <- tmp[is.na(naa.tmp)]
          naa[,i+1, ] <- naa.tmp

          ## 当年の加入の計算
          if(fyears[i+1]-min.age < start.year){
              # 参照する親魚資源量がVPA期間である場合、VPA期間のSSBをとってくる
              thisyear.ssb <- sum(res0$ssb[,as.character(fyears[i+1]-min.age)],na.rm=T)*res0$input$unit.waa/res0$input$unit.biom
              if(!is.null(ssb0)) thisyear.ssb <- colSums(ssb0)
          }
          else{
              # そうでない場合
            if(waa.fun){
                # 動的なwaaは対応する年のwaaを書き換えた上で使う？
                waa[2:nage,i+1-min.age,] <- t(sapply(2:nage, function(ii) as.numeric(exp(WAA.b0[ii]+WAA.b1[ii]*log(naa[ii,i+1-min.age,])+waa.rand[ii,i+1-min.age,]))))

            }
            thisyear.ssb <- colSums(naa[,i+1-min.age,]*waa[,i+1-min.age,]*maa[,i+1-min.age,],na.rm=T)*res0$input$unit.waa/res0$input$unit.biom
          }

          thisyear.ssb <- thisyear.ssb+(1e-10)
          rec.tmp <- recfunc(thisyear.ssb,res0,
                             rec.resample=rec.tmp$rec.resample,
                             rec.arg=rec.arg)
          if(is.na(naa[1,i+1,1]))  naa[1,i+1,] <- rec.tmp$rec
#          if(!is.null(rec.tmp$rec.arg)) rec.arg <- rec.tmp$rec.arg
          rps.mat[i+1,] <- naa[1,i+1,]/thisyear.ssb
          eaa[i,] <- rec.tmp$rec.resample[1:N]
          rec.arg$resid <- rec.tmp$rec.resample # ARオプションに対応
      }

#      if(Pope){
#        caa[,ntime,] <- naa[,ntime,]*(1-exp(-faa[,ntime,]))*exp(-M[,ntime,]/2)
#      }
#      else{
#        caa[,ntime,] <- naa[,ntime,]*(1-exp(-faa[,ntime,]-M[,ntime,]))*faa[,ntime,]/(faa[,ntime,]+M[,ntime,])
#      }
      if (!is.null(rec.arg$rho)) rec.tmp$rec.resample <- NULL


#      if(Pope){
#          caa[,i,] <- naa[,i,]*(1-exp(-faa[,i,]))*exp(-M[,i,]/2)
#      }
#      else{
#          caa[,i,] <- naa[,i,]*(1-exp(-faa[,i,]-M[,i,]))*faa[,i,]/(faa[,i,]+M[,i,])
#      }

      if(Pope){
          caa[] <- naa*(1-exp(-faa))*exp(-M/2)
      }
      else{
          caa[] <- naa*(1-exp(-faa-M))*faa/(faa+M)
      }

      caa <- caa[,-ntime,,drop=F]
      waa.catch <- waa.catch[,-ntime,,drop=F]
      waa <- waa[,-ntime,,drop=F]
      maa <- maa[,-ntime,,drop=F]
      naa <- naa[,-ntime,,drop=F]
      faa <- faa[,-ntime,,drop=F]
      M <- M[,-ntime,,drop=F]
      fyears <- fyears[-ntime]

      biom <- naa*waa*res0$input$unit.waa/res0$input$unit.biom
      ssb <- naa*waa*maa*res0$input$unit.waa/res0$input$unit.biom

      wcaa <- caa*waa.catch*res0$input$unit.waa/res0$input$unit.biom
      vwcaa <- apply(wcaa,c(2,3),sum,na.rm=T)

      ABC <- apply(as.matrix(vwcaa[fyears%in%ABC.year,,drop=F]),2,sum)

      if(!is.null(rec.arg$resample)) if(rec.arg$resample==TRUE) eaa[] <- NA # resamplingする場合にはeaaにはなにも入れない

      if(outtype=="FULL"){
          fres <- list(faa=faa,naa=naa,biom=biom,baa=biom,ssb=ssb,wcaa=wcaa,caa=caa,M=M,rps=rps.mat,
                       maa=maa,vbiom=apply(biom,c(2,3),sum,na.rm=T),
                       eaa=eaa,
                       waa=waa,waa.catch=waa.catch,currentF=currentF,
                       vssb=apply(ssb,c(2,3),sum,na.rm=T),vwcaa=vwcaa,
                       years=fyears,fyear.year=fyear.year,ABC=ABC,recfunc=recfunc,rec.arg=rec.arg,
                       waa.year=waa.year,maa.year=maa.year,multi=multi,multi.year=multi.year,
                       Frec=Frec,rec.new=rec.new,pre.catch=pre.catch,input=arglist)
    }
      else{
          fres <- list(faa=faa[,,1],M=M[,,1],recruit=naa[1,,],eaa=eaa,baa=biom,
                       maa=maa[,,1],vbiom=apply(biom,c(2,3),sum,na.rm=T),
                       waa=waa[,,1],waa.catch=waa.catch[,,1],currentF=currentF,
                       vssb=apply(ssb,c(2,3),sum,na.rm=T),vwcaa=vwcaa,
                       years=fyears,fyear.year=fyear.year,ABC=ABC,recfunc=recfunc,rec.arg=rec.arg,
                       waa.year=waa.year,maa.year=maa.year,multi=multi,multi.year=multi.year,
                       Frec=Frec,rec.new=rec.new,pre.catch=pre.catch,input=arglist)
      }

      ## if(non.det==TRUE){
      ##     fres <- list(faa=faa[,,-1,drop=F],naa=naa[,,-1,drop=F],biom=biom[,,-1,drop=F],
      ##                  ssb=ssb[,,-1,drop=F],wcaa=wcaa[,,-1,drop=F],caa=caa[,,-1,drop=F],
      ##                  M=M[,,-1,drop=F],rps=rps.mat[,-1,drop=F],
      ##                  maa=maa[,,-1,drop=F],vbiom=apply(biom[,,-1,drop=F],c(2,3),sum,na.rm=T),
      ##                  eaa=eaa[,-1,drop=F],
      ##                  waa=waa[,,-1,drop=F],waa.catch=waa.catch[,,-1,drop=F],currentF=currentF,
      ##                  vssb=apply(ssb[,,-1,drop=F],c(2,3),sum,na.rm=T),vwcaa=vwcaa[,-1,drop=F],
      ##                  years=fyears,fyear.year=fyear.year,ABC=ABC,recfunc=recfunc,rec.arg=rec.arg,
      ##                  waa.year=waa.year,maa.year=maa.year,multi=multi,multi.year=multi.year,
      ##                  Frec=Frec,rec.new=rec.new,pre.catch=pre.catch,input=arglist)
      ## }

      class(fres) <- "future"
      if(is.plot){
          par(mfrow=c(2,2))
          plot.future(fres)
      }
      if(waa.fun) fres$waa.reg <- WAA.res
      invisible(fres)
  }



forward.calc.mat2 <- function(fav,nav,Mv,plus.group=TRUE){
  nage <- max(which(!is.na(nav[,1])))#length(fav)
  na.age <- which(is.na(nav[-1,1]))
#  naa <- matrix(NA,nage,dim(nav)[[2]])
  naa <- matrix(NA,dim(nav)[[1]],dim(nav)[[2]])
#  for(a in 2:(nage-1)){
  naa[-c(1,nage,na.age),] <- nav[-c(nage,nage-1,na.age),]*
      exp(-fav[-c(nage,nage-1,na.age),]-Mv[-c(nage,nage-1,na.age),])
#  }
  naa[nage,] <- nav[nage-1,]*exp(-fav[nage-1,]-Mv[nage-1,])
  pg <- nav[nage,]*exp(-fav[nage,]-Mv[nage,])
  if(plus.group) naa[nage,] <- naa[nage,] + pg
  return(naa)
}

caa.est.mat <- function(naa,saa,waa,M,catch.obs,Pope){
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

