############
# RVPAの結果からMSYを計算する関数
# 主に使うのはSR.est(再生産関係をフィットし、MSYを計算)とSR.plot（フィットした結果をプロット）
############

############
# 使い方
############
if(0){
                                        # マサバ太平洋のデータを読み込み; modelAはvpaの帰り値
    modelA <- readRDS("modelA_res.Rdata")
                                        # MSY計算
    res1 <- SR.est(modelA,
                   what.est=c(TRUE,TRUE,TRUE), # HS,BH,RIのどれをフィットするか。
                   bref.year=2013:2015, # 生物パラメータを用いる期間
                   years=c(1970:2013), # 観測されたSR関係を用いる期間
                   er.log=TRUE, # 誤差。TRUEで対数正規誤差
                   fc.year=2013:2015, # MSY計算のさいに選択率を平均する期間
                   seed=1 # 乱数の種。この値を変えると乱数が変わるので結果も変わる
                   )

    res1$summary # 推定パラメータ、管理基準値の確認
                                        # 再生産パラメータa,bはエクセルとほぼ一致するはずだが、管理基準値は確率的シミュレーションをもとに計算しているので、エクセルとは必ずしも一致しない。±５％くらいの違いはあるみたい

                                        # 結果のプロット(HSのみ)
    res1.pred <- plot.SR(res1,what.plot=c("hs"))
                                        # 結果のプロット(HS,BH,RIを全て)
    res1.pred <- plot.SR(res1,what.plot=c("hs","bh","ri"))
    allplot(res1) # 要約表・グラフの出力

}

############
# fit to S-R relationship & MSE estimation -- old version?
############

SR.est <- function(vpares,SSB.dat=NULL,R.dat=NULL,gamma1=0.0001,er.log=TRUE,
                   years=as.numeric(colnames(vpares$naa)), # 親子関係に推定に使う年のベクトル
                   bref.year=2011:2013,# B0やMSYを計算するさいの生物パラメータの範囲(2011-2013に変更、2016-06-06）
                   fc.year=bref.year, # 将来予測をするときに仮定する選択率をとる年の範囲
                   seed=1,n.imputation=1,
                   nyear=100,
                   bias.correction=TRUE, # 確率的な将来予測を行うときにbias correctionをするかどうか
                   eyear=0, # 将来予測の最後のeyear+1年分を平衡状態とする
#                   FUN=median, # 漁獲量の何を最大化するか？
                   FUN=mean, # 漁獲量の何を最大化するか？
                   sigma=-1, #加入変動のCV。-1の場合にはobservedの値を使う
                   N=1000, # stochastic計算するときの繰り返し回数
                   is.small=FALSE, # 将来予測の結果を返さない。
                   is.boot=1000,# 正の値であれば、SRフィットのノンパラメトリックブートストラップを行う
                   is.Kobe=c(FALSE,FALSE,FALSE), # Kobeの計算をするかどうか。順番に、HS, BH, RIの順
                   is.5perlower=FALSE, # HSの折れ点を5%の確率で下回るときの親魚資源量
                   PGY=NULL, # PGY管理基準値を計算するかどうか。計算しない場合はNULLを、計算する場合はc(0.8,0.9,0.95)のように割合を入れる
                   what.est=c(TRUE,TRUE,TRUE) # MSY等を推定するか。順番に、HS, BH, RIの順
                   ){

    #####-------- 内部で使う関数の定義
    HS <- function(p,R,SSB,gamma=gamma1,er.log=er.log,MLE=FALSE,a=NULL,b=NULL){
        if(!is.null(a)) p[1] <- a
        if(!is.null(b)) p[2] <- b

        a <- exp(p[1])
        b <- max(SSB)/(1+exp(-p[2]))
        if(isTRUE(MLE)) sigma <- exp(p[3])
        Pred <- function(SSB) a*(SSB+sqrt(b^2+gamma^2/4)-sqrt((SSB-b)^2+gamma^2/4))
        if(!isTRUE(MLE)){
            if(er.log==FALSE)   return(sum((R-Pred(SSB))^2))
            else return(sum((log(R)-log(Pred(SSB)))^2))
        }
        else{
            if(er.log==FALSE){
                obj <- length(R)*log(1/(sqrt(2*pi)*sigma))-1/2/sigma^2*sum( (R-Pred(SSB))^2 )
                return(-obj)
            }
            else{
                obj <- length(R)*log(1/(sqrt(2*pi)*sigma))-1/2/sigma^2*sum( (log(R)-log(Pred(SSB)))^2 )
                return(-obj)
            }
        }
    }

    BH <- function(p,R,SSB,er.log=er.log,MLE=FALSE){
        a <- exp(p[1])
        b <- exp(p[2])
        if(isTRUE(MLE)) sigma <- exp(p[3])

        Pred <- function(SSB) a*SSB/(1+b*SSB)

        if(!isTRUE(MLE)){
            if(er.log==FALSE) return(sum((R-Pred(SSB))^2))
            else return(sum((log(R)-log(Pred(SSB)))^2))
        }
        else{
            if(er.log==FALSE){
                obj <- length(R)*log(1/(sqrt(2*pi)*sigma))-1/2/sigma^2*sum( (R-Pred(SSB))^2 )
                return(-obj)
            }
            else{
                obj <- length(R)*log(1/(sqrt(2*pi)*sigma))-1/2/sigma^2*sum( (log(R)-log(Pred(SSB)))^2 )
                return(-obj)
            }
        }
    }


    SL <- function(p,R,SSB,er.log=er.log,MLE=FALSE){
        a <- exp(p[1])
#        b <- exp(p[2])
        if(isTRUE(MLE)) sigma <- exp(p[2])

        Pred <- function(SSB) a*SSB

        if(!isTRUE(MLE)){
            if(er.log==FALSE) return(sum((R-Pred(SSB))^2))
            else return(sum((log(R)-log(Pred(SSB)))^2))
        }
        else{
            if(er.log==FALSE){
                obj <- length(R)*log(1/(sqrt(2*pi)*sigma))-1/2/sigma^2*sum( (R-Pred(SSB))^2 )
                return(-obj)
            }
            else{
                obj <- length(R)*log(1/(sqrt(2*pi)*sigma))-1/2/sigma^2*sum( (log(R)-log(Pred(SSB)))^2 )
                return(-obj)
            }
        }
    }


    RI <- function(p,R,SSB,er.log=er.log,MLE=FALSE){
        a <- exp(p[1])
        b <- exp(p[2])
        if(isTRUE(MLE)) sigma <- exp(p[3])

        Pred <- function(SSB) a*SSB*exp(-b*SSB)

        if(!isTRUE(MLE)){
            if(er.log==FALSE) return(sum((R-Pred(SSB))^2))
            else return(sum((log(R)-log(Pred(SSB)))^2))
        }
        else{
            if(er.log==FALSE){
                obj <- length(R)*log(1/(sqrt(2*pi)*sigma))-1/2/sigma^2*sum( (R-Pred(SSB))^2 )
                return(-obj)
            }
            else{
                obj <- length(R)*log(1/(sqrt(2*pi)*sigma))-1/2/sigma^2*sum( (log(R)-log(Pred(SSB)))^2 )
                return(-obj)
            }
        }
    }

    # HSを推定するための関数
    get.HS <- function(R,SSB,er.log,gamma1,do.profile=TRUE){
        reg0 <- lm(R~SSB-1)
        a0 <- reg0$coef
        b0 <- 0.9
        # hockey-stick
        res.HS <-  optim(c(log(a0),logit(b0)),HS,R=R,SSB=SSB,er.log=er.log,gamma=gamma1)
        s <- 1
        for (j in seq(0.95,0.1,by=-0.05)){
            res.HS0 <-  optim(c(log(a0),logit(j)),HS,R=R,SSB=SSB,er.log=er.log,gamma=gamma1)
            if (res.HS0$value < res.HS$value) res.HS <- res.HS0
        }
        res.HS <-  optim(res.HS$par,HS,R=R,SSB=SSB,method="BFGS",er.log=er.log,gamma=gamma1)
        ofv.least.square <- res.HS$value

        res.HS <-  optim(c(res.HS$par,log(sqrt(res.HS$value/length(R)))),HS,R=R,SSB=SSB,method="BFGS",er.log=er.log,gamma=gamma1,MLE=TRUE)

        a.HS <- exp(res.HS$par[1])
        names(a.HS) <- NULL
        b.HS <- max(SSB)/(1+exp(-res.HS$par[2])) # 曲がるところのx軸(ssb_hs)
        # r0の計算
        r0.HS <- pred.HS(b.HS,a=a.HS,b=b.HS,gamma=gamma1)
        # もし、b.HSが最大・最小SSBよりも大きい・小さかったら
        if(b.HS>max(SSB)|b.HS<min(SSB)){
            b.HS <- ifelse(b.HS>max(SSB),max(SSB),b.HS)
            b.HS <- ifelse(b.HS<min(SSB),min(SSB),b.HS)
            tmpfunc <- function(x,r0,...) (pred.HS(a=x,...)-r0)^2
            tmp <- optimize(tmpfunc,c(0,a.HS*10),b=b.HS,gamma=gamma1,SSB=b.HS,r0=r0.HS)
            a.HS <- tmp$minimum
        }

        # 尤度surfaceの計算
        if(isTRUE(do.profile)){
            a.grid <- c(seq(from=0.1,to=0.9,by=0.1),seq(from=0.91,to=1.09,by=0.02),seq(from=1.1,to=1.5,by=0.1)) * a.HS
            b.grid <- c(seq(from=0.1,to=0.9,by=0.1),seq(from=0.91,to=1.09,by=0.02),seq(from=1.1,to=1.5,by=0.1)) * b.HS
            b.grid <- b.grid[b.grid<max(SSB)]
            obj.data <- expand.grid(a=a.grid,b=b.grid)
            obj.data$obj <- NA
            obj.data$log.a <- log(obj.data$a)
            obj.data$conv.b <- -log(max(SSB)/obj.data$b-1)
            for(i in 1:nrow(obj.data))
            {
                obj.data$obj[i] <- HS(c(obj.data$log.a[i],obj.data$conv.b[i]),R=R,SSB=SSB,
                                      MLE=FALSE,er.log=er.log,gamma=gamma1)

            }
        }
        else{
            obj.data <- NA
        }

        return(list(a=a.HS,b=b.HS,r0=r0.HS,res=res.HS,obj.data=obj.data,ofv.least.square=ofv.least.square))
    }


    # Beverton-Holt
    get.BH <- function(R,SSB,er.log){
        reg0 <- lm(R~SSB-1)
        a0 <- reg0$coef
        b0 <- max(SSB)
        res.BH <-  optim(c(log(a0),log(1/b0)),BH,R=R,SSB=SSB,method="BFGS",er.log=er.log)
        for (j in seq(0.9,0.1,by=-0.1)){
            res.BH0 <-  optim(c(log(a0),log(j/b0)),BH,R=R,SSB=SSB,er.log=er.log)
            if (res.BH0$value < res.BH$value) res.BH <- res.BH0
        }

        # 最尤法で計算しなおしたもので上書き
        res.BH <-  optim(c(res.BH$par,log(sqrt(res.BH$value/length(R)))),BH,R=R,SSB=SSB,method="BFGS",er.log=er.log,MLE=TRUE)

        a.BH <- exp(res.BH$par[1])
        b.BH <- exp(res.BH$par[2])
        return(list(a=a.BH,b=b.BH,res=res.BH))
    }


    get.RI <- function(R,SSB,er.log){
        reg0 <- lm(R~SSB-1)
        a0 <- reg0$coef
        b0 <- max(SSB)
        # Ricker
        res.RI <- optim(c(log(a0),log(1/b0)),RI,R=R,SSB=SSB,method="BFGS",er.log=er.log)
        for (j in seq(0.9,0.1,by=-0.1)){
            res.RI0 <-  optim(c(log(a0),log(j/b0)),RI,R=R,SSB=SSB,er.log=er.log)
            if (res.RI0$value < res.RI$value) res.RI <- res.RI0
        }
        #　最尤法
        res.RI <- optim(c(res.RI$par,log(sqrt(res.RI$value/length(R)))),
                        RI,R=R,SSB=SSB,method="BFGS",er.log=er.log,MLE=TRUE)
        a.RI <- exp(res.RI$par[1])
        b.RI <- exp(res.RI$par[2])
        return(list(a=a.RI,b=b.RI,res=res.RI))
    }
    ##### 関数定義終わり

    # R.datとSSB.datだけが与えられた場合、それを使ってシンプルにフィットする
    if(!is.null(R.dat) & !is.null(SSB.dat)){
        dat <- data.frame(R=R.dat,SSB=SSB.dat,years=1:length(R.dat))
    }
    else{
        vpares$Fc.at.age <- rowMeans(vpares$faa[as.character(fc.year)])

    # データの整形
        n <- ncol(vpares$naa)
        L <- as.numeric(rownames(vpares$naa)[1])

        dat <- list()
        dat$R <- as.numeric(vpares$naa[1,])
        dat$SSB <- as.numeric(colSums(vpares$ssb))
        dat$year <- as.numeric(colnames(vpares$ssb))
    # 加入年齢分だけずらす
        dat$R <- dat$R[(L+1):n]
        dat$SSB <- dat$SSB[1:(n-L)]
        dat$year <- dat$year[(L+1):n]

                                        # データの抽出
        dat <- as.data.frame(dat)
        dat <- dat[dat$year%in%years,]
    }

    R <- dat$R
    SSB <- dat$SSB

    # HS推定
#    if(what.est[1]==TRUE){
        tmp <- get.HS(R,SSB,er.log,gamma1)
        a.HS <- tmp$a; b.HS <- tmp$b ; r0.HS <- tmp$r0
        sd.HS <- exp(tmp$res$par[3])
        surface.HS <- tmp$obj.data
        ofv.least.square <- tmp$ofv.least.square
        res.HS <- tmp$res
        boot.HS <- matrix(NA,is.boot,3)
        jack.HS <- matrix(NA,length(R),3)
        colnames(boot.HS) <- colnames(jack.HS) <-  c("a","b","r0")
        if(what.est[1]==TRUE&&is.boot>0){ # ブートストラップ
            for(i in 1:is.boot){
                rand <- sample(length(R),size=length(R)*n.imputation,replace=TRUE)
                tmp <- get.HS(R[rand],SSB[rand],er.log,gamma1,do.profile=FALSE)
                boot.HS[i,] <- unlist(tmp[c("a","b","r0")])
            }
            for(i in 1:length(R)){
                tmp <- get.HS(R[-i],SSB[-i],er.log,gamma1,do.profile=FALSE)
                jack.HS[i,] <- unlist(tmp[c("a","b","r0")])
            }
#            rownames(jack.HS) <- years
        }
        # 予測値
        dat$pred.HS <- pred.HS(dat$SSB,a=a.HS,b=b.HS,gamma=gamma1)
        dat$log.resid.HS <- log(dat$R) - log(dat$pred.HS)
#    }

    if(0){
        # 直線回帰
        reg0 <- lm(R~SSB-1)
        a0 <- reg0$coef
        res.SL <-  optimize(SL,c(0,log(a0)*10),R=R,SSB=SSB,er.log=er.log)
        res.SL <-  optim(c(res.SL$minimum,log(sqrt(res.SL$objective/length(R)))),
                         SL,R=R,SSB=SSB,er.log=er.log,MLE=TRUE)
        #    res.SL$value <- res.SL$objective
        a.SL <- exp(res.SL$par[1])
        boot.SL <- rep(NA,is.boot)
        jack.SL <- rep(NA,length(R))
        if(is.boot>0){
            for(i in 1:is.boot){
                rand <- sample(length(R),replace=TRUE)
                tmp <-  optimize(SL,c(0,log(a0)*10),R=R[rand],SSB=SSB[rand],er.log=er.log)
                boot.SL[i] <- exp(tmp$minimum[1])
            }
            for(i in 1:length(R)){
                tmp <-  optimize(SL,c(0,log(a0)*10),R=R[-i],SSB=SSB[-i],er.log=er.log)
                jack.SL[i] <- exp(tmp$minimum[1])
            }
            rownames(jack.SL) <- years
        }
    }

    if(what.est[2]==TRUE){
        # BH推定
        tmp <- get.BH(R,SSB,er.log)
        a.BH <- tmp$a; b.BH <- tmp$b; res.BH <- tmp$res
        sd.BH <- exp(tmp$res$par[3])
        boot.BH <- matrix(NA,is.boot,2)
        jack.BH <- matrix(NA,length(R),2)
        colnames(boot.BH) <- colnames(jack.BH) <-  c("a","b")
        if(is.boot>0){ # ブートストラップ
            for(i in 1:is.boot){
                rand <- sample(length(R),replace=TRUE)
                tmp <- get.BH(R[rand],SSB[rand],er.log)
                boot.BH[i,] <- unlist(tmp[c("a","b")])
            }
            # ジャックナイフも
            for(i in 1:length(R)){
                tmp <- get.BH(R[-i],SSB[-i],er.log)
                jack.BH[i,] <- unlist(tmp[c("a","b")])
            }
            rownames(jack.BH) <- years
        }
        ##
        dat$pred.BH <- pred.BH(dat$SSB,a=a.BH,b=b.BH)
        dat$log.resid.BH <- log(dat$R) - log(dat$pred.BH)
    }

    if(what.est[3]==TRUE){
        # RI推定
        tmp <- get.RI(R,SSB,er.log)
        a.RI <- tmp$a ; b.RI <- tmp$b ; res.RI <- tmp$res
        sd.RI <- exp(tmp$res$par[3])
        boot.RI <- matrix(NA,is.boot,2)
        jack.RI <- matrix(NA,length(R),2)
        colnames(boot.RI) <- colnames(jack.RI) <-  c("a","b")
        if(is.boot>0){ # ブートストラップ
            for(i in 1:is.boot){
                rand <- sample(length(R),replace=TRUE)
                tmp <- get.RI(R[rand],SSB[rand],er.log)
                boot.RI[i,] <- unlist(tmp[c("a","b")])
            }
            # ジャックナイフも
            for(i in 1:length(R)){
                tmp <- get.RI(R[-i],SSB[-i],er.log)
                jack.RI[i,] <- unlist(tmp[c("a","b")])
            }
            rownames(jack.RI) <- years
        }
        ##
        dat$pred.RI <- pred.RI(dat$SSB,a=a.RI,b=b.RI)
        dat$log.resid.RI <- log(dat$R) - log(dat$pred.RI)
    }

    # 単に回帰だけする場合
    if(!is.null(R.dat) & !is.null(SSB.dat)){
        res <- list()
        paste2 <- function(x,...) paste(x,...,sep="")
        for(j in which(what.est)){
            SR <- c("HS","BH","RI")
            xx <- c(get(paste2("a.",SR[j])),
                    get(paste2("b.",SR[j])),
                    get(paste2("sd.",SR[j])),
                    get(paste2("res.",SR[j]))$value)
            names(xx) <- c("a","b","sd","value")
            res[[j]] <- list(parameter=xx,
                             boot=get(paste2("boot.",SR[j])),
                             jack=get(paste2("jack.",SR[j])))
            }
        names(res) <- SR[what.est]
        return(res)
    }

    #-------------------- B0 & MSY for HS --------------------
    # function to minimize

    # シミュレーション回数ぶんの漁獲量のFUN（mean, geomean, median）を最大化するFを選ぶ
    tmpfunc <- function(x,f.arg,FUN=FUN,eyear=eyear){
      f.arg$multi <- x
      fout <- do.call(future.vpa2,f.arg)
      return(-FUN(fout$vwcaa[(nrow(fout$vwcaa)-eyear):nrow(fout$vwcaa),-1]))
    }

    tmpfunc2 <- function(x,f.arg,FUN=FUN,eyear=eyear,hsp=0){
      f.arg$multi <- x
      fout <- do.call(future.vpa2,f.arg)
      tmp <- as.numeric(fout$vssb[(nrow(fout$vssb)-eyear):nrow(fout$vssb),-1])
      lhs <- sum(tmp<hsp)/length(tmp)
      return( (lhs-0.05)^2 + as.numeric(lhs==0) + as.numeric(lhs==1)  )
    }

    get.Fhist <- function(farg,vpares,eyear,trace,hsp=0){
        Fhist <- NULL
        original.sel <- farg$res0$Fc.at.age # original F
        for(j in 1:ncol(vpares$faa)){
            farg$res0$Fc.at.age <- vpares$faa[,j] # change the selectivity
            farg$multi <- 1
            tmp <- do.call(future.vpa2,farg)
            tmp2 <- get.stat(tmp,eyear=eyear,hsp=hsp)
#            browser()
            xx <- which.min(abs(trace$ssb.median-tmp2$ssb.median))+c(-1,1)
            range.tmp <- trace$fmulti[xx]
            if(is.na(range.tmp[2])) range.tmp[2] <- max(trace$fmulti)*2
            if(xx[1]==0) range.tmp <- c(0,range.tmp[1])
            tmpfunc <- function(x,farg,ssb.target,eyear){
                farg$multi <- x
                return((get.stat(do.call(future.vpa2,farg),eyear=eyear,hsp=hsp)$ssb.mean-ssb.target)^2)
            }
            farg$res0$Fc.at.age <- original.sel # current Fをもとにもどす
            # originalな選択率のもとで、それを何倍にすればi年目のFで漁獲した時の親魚資源量と同じになるか
            ores <- optimize(tmpfunc,range.tmp,farg=farg,ssb.target=tmp2$ssb.mean,eyear=eyear)
#            farg$multi <- ores$minimum
#            tmp3 <- do.call(future.vpa2,farg)
            tmp2$fmulti <- ores$minimum
            Fhist <- rbind(Fhist,tmp2)
        }
        return(as.data.frame(Fhist))
    }

    trace.func <- function(farg,eyear,hsp=0,
                           fmulti=c(seq(from=0,to=0.9,by=0.1),1,seq(from=1.1,to=2,by=0.1),3:5,7,20,100)){
        trace.res <- NULL
        farg$outtype <- "FULL"
        for(i in 1:length(fmulti)){
            farg$multi <- fmulti[i]
            tmp <- do.call(future.vpa2,farg)
#            trace.res <- rbind(trace.res,get.stat(tmp,eyear=eyear,hsp=hsp))
            tmp2 <- cbind(get.stat(tmp,eyear=eyear,hsp=hsp),
                          get.stat2(tmp,eyear=eyear,hsp=hsp))
            trace.res <- rbind(trace.res,tmp2)
            if(tmp2$"ssb.mean"<trace.res$"ssb.mean"[1]/1000){
                fmulti <- fmulti[1:i]
                break()
            }
          }
        trace.res <- as.data.frame(trace.res)
        trace.res$fmulti <- fmulti
        return(trace.res)
    }

    b0.HS <- b0.BH <- b0.RI <- numeric() # B0
    fout.HS <- fout.BH <- fout.RI <- list()
    fout0.HS <- fout0.BH <- fout0.RI <- list()
    trace.HS <- trace.BH <- trace.RI <- list()
    Fhist.HS <- Fhist.BH <- Fhist.RI <- list()
    fout.HS.5per <- list()

    for(kk in 1:length(sigma)){
      ref.year <- as.numeric(rev(colnames(vpares$naa))[1])
      if(sigma[kk]==-1){
          if(isTRUE(what.est[1])){
              sigma.tmp <- exp(res.HS$par[3])
          }
          else{
              if(isTRUE(what.est[2]))  sigma.tmp <- exp(res.BH$par[3])
              if(isTRUE(what.est[3]))  sigma.tmp <- exp(res.RI$par[3])
          }
      }
      else{
          sigma.tmp <- sigma[kk]
      }

      #--------- Hockey stick
      fout0.HS[[kk]] <- future.vpa2(vpares,multi=0,nyear=nyear,start.year=ref.year,
                          N=ifelse(sigma[kk]==0,1,N),
                          ABC.year=ref.year+1,waa.year=bref.year,maa.year=bref.year,
                          M.year=bref.year,is.plot=FALSE,
                          recfunc=HS.rec,seed=seed,outtype="simple",
                          rec.arg=list(a=a.HS,b=b.HS,gamma=gamma1,sd=sigma.tmp,bias.correction=bias.correction))

#      b0.HS[kk] <- fout0$vssb[nrow(fout0$vssb),1] #static b0
      farg.HS <- fout0.HS[[kk]]$input

      which.min2 <- function(x){
          max(which(min(x)==x))
      }

      if(isTRUE(what.est[1])){
          trace.HS[[kk]] <- trace.func(farg.HS,eyear,hsp=b.HS)

          xx <- which.max(trace.HS[[kk]]$catch.median)+c(-1,1)
          range.tmp <- trace.HS[[kk]]$fmulti[xx]
          if(xx[1]==0) range.tmp <- c(0,range.tmp)
#          if(is.na(xx[2])) range.tmp[2] <- max(trace.HS[[kk]]$fmulti)*10
          if(is.na(range.tmp[2])) range.tmp[2] <- max(trace.HS[[kk]]$fmulti)*10

          tmp <- optimize(tmpfunc,range.tmp,f.arg=farg.HS,eyear=eyear,FUN=FUN)
          farg.HS$multi <- tmp$minimum # Fc.at.a * multiがFmsy
          fout.HS[[kk]] <- do.call(future.vpa2,farg.HS)
          Fmsy.HS <- tmp$minimum * farg.HS$res0$Fc.at.age

          ## ここでtraceを追加
          trace.HS[[kk]] <- rbind(trace.HS[[kk]],trace.func(farg.HS,eyear,hsp=b.HS,
                                  fmulti=tmp$minimum+c(-0.025,-0.05,-0.075,0,0.025,0.05,0.075)))
          trace.HS[[kk]] <- trace.HS[[kk]][order(trace.HS[[kk]]$fmulti),]
          ###

          if(is.Kobe[1]) Fhist.HS[[kk]] <- get.Fhist(farg.HS,vpares,eyear=eyear,trace=trace.HS[[kk]])
          if(is.5perlower){
              xx <- which.min2((trace.HS[[kk]]$lower-0.05)^2)+c(-1,1)
              range.tmp <- trace.HS[[kk]]$fmulti[xx]
              if(xx[1]==0) range.tmp <- c(0,range.tmp)
              if(is.na(xx[2])) range.tmp[2] <- max(trace.HS[[kk]]$fmulti)*10
              tmp <- optimize(tmpfunc2,range.tmp,f.arg=farg.HS,eyear=eyear,FUN=FUN,hsp=b.HS)
              farg.HS$multi <- tmp$minimum
              fout.HS.5per[[kk]] <- do.call(future.vpa2,farg.HS)
          }
      }

      #---------------------- calculation of MSY for BH
      if(isTRUE(what.est[2])){
          if(sigma[kk]==-1){
              sigma.tmp <- exp(res.BH$par[3])
          }
          else{
              sigma.tmp <- sigma[kk]
          }
          farg.BH <- farg.HS
          farg.BH$recfunc <- BH.rec
          farg.BH$rec.arg <- list(a=a.BH,b=b.BH,sd=sigma.tmp,bias.correction=bias.correction)
          farg.BH$multi <- 0
          fout0.BH[[kk]] <- do.call(future.vpa2,farg.BH)
          #      b0.BH[kk] <- fout0.BH$vssb[nrow(fout0$vssb),1] #static b0

          trace.BH[[kk]] <- trace.func(farg.BH,eyear)
          #      tmp <- optimize(tmpfunc,c(0,10),f.arg=farg.BH,eyear=eyear,FUN=FUN)
          xx <- which.max(trace.BH[[kk]]$catch.median)+c(-1,1)
          range.tmp <- trace.BH[[kk]]$fmulti[xx]
          if(xx[1]==0) range.tmp <- c(0,range.tmp)
#          if(is.na(xx[2])) range.tmp[2] <- max(trace.BH[[kk]]$fmulti)*10
          if(is.na(range.tmp[2])) range.tmp[2] <- max(trace.BH[[kk]]$fmulti)*10
          tmp <- optimize(tmpfunc,range.tmp,f.arg=farg.BH,eyear=eyear,FUN=FUN)


          farg.BH$multi <- tmp$minimum
          fout.BH[[kk]] <- do.call(future.vpa2,farg.BH)
          Fmsy.BH <- tmp$minimum * farg.BH$res0$Fc.at.age
          if(is.Kobe[2])  Fhist.BH[[kk]] <- get.Fhist(farg.BH,vpares,eyear=eyear,trace.BH[[kk]])

          ## ここでtraceを追加
          trace.BH[[kk]] <- rbind(trace.BH[[kk]],
                                  trace.func(farg.BH,eyear,hsp=b.BH,fmulti=tmp$minimum+c(-0.025,-0.05,-0.075,0,0.025,0.05,0.075)))
          trace.BH[[kk]] <- trace.BH[[kk]][order(trace.BH[[kk]]$fmulti),]
          ###
      }

      #------------------- calculation of MSY for RI
      if(isTRUE(what.est[3])){
          if(sigma[kk]==-1){
              sigma.tmp <- exp(res.RI$par[3])
          }
          else{
              sigma.tmp <- sigma[kk]
          }
          farg.RI <- farg.HS
          farg.RI$recfunc <- RI.rec
          farg.RI$rec.arg <- list(a=a.RI,b=b.RI,sd=sigma.tmp,bias.correction=bias.correction)
          farg.RI$multi <- 0
          fout0.RI[[kk]] <- do.call(future.vpa2,farg.RI)
          #      b0.RI[kk] <- fout0$vssb[nrow(fout0$vssb),1] #static b0

          trace.RI[[kk]] <- trace.func(farg.RI,eyear)

          xx <- which.max(trace.RI[[kk]]$catch.median)+c(-1,1)
          range.tmp <- trace.RI[[kk]]$fmulti[xx]
          if(xx[1]==0) range.tmp <- c(0,range.tmp)
#          if(is.na(xx[2])) range.tmp[2] <- max(trace.RI[[kk]]$fmulti)*10
          if(is.na(range.tmp[2])) range.tmp[2] <- max(trace.RI[[kk]]$fmulti)*10

          tmp <- optimize(tmpfunc,range.tmp,f.arg=farg.RI,eyear=eyear,FUN=FUN)

          farg.RI$multi <- tmp$minimum
          fout.RI[[kk]] <- do.call(future.vpa2,farg.RI)
          Fmsy.RI <- tmp$minimum * farg.RI$res0$Fc.at.age

          if(is.Kobe[3])  Fhist.RI[[kk]] <- get.Fhist(farg.RI,vpares,eyear=eyear,trace.RI[[kk]])

          ## ここでtraceを追加
          trace.RI[[kk]] <- rbind(trace.RI[[kk]],
                                  trace.func(farg.RI,eyear,hsp=b.RI,
                                             fmulti=tmp$minimum+c(-0.025,-0.05,-0.075,0,0.025,0.05,0.075)))
          trace.RI[[kk]] <- trace.RI[[kk]][order(trace.RI[[kk]]$fmulti),]
          ###
      }
    }
    #--------------------------------------

    if(isTRUE(is.5perlower)){
        tmp <- as.data.frame(t(sapply(fout.HS.5per,get.stat,eyear=eyear,hsp=b.HS)))
        tmp$f <- sapply(fout.HS.5per,function(x)x$multi)
    }
    else{
        tmp <- NA
    }

    # 関数を返すとsaveしたときに異常にファイルサイズが大きくなる。原因は不明。
    # とりあえず、関数を返すのをやめる
    output <- list(dat=dat,sigma=sigma,vpares=vpares)
    if(what.est[1]==TRUE)
        output$hs <- list(a=a.HS,b=b.HS,sd=sd.HS,gamma=gamma1,ofv=res.HS$value,ofv.least.square=ofv.least.square,
                           res=res.HS,r0=r0.HS,Fhist=Fhist.HS,
                           trace=trace.HS,boot=as.data.frame(boot.HS),
                           jack=as.data.frame(jack.HS),farg=farg.HS,
                           f.msy=sapply(fout.HS,function(x)x$multi),
                          Fmsy=Fmsy.HS,surface=surface.HS,
                          fout=fout.HS,
                         # 最大化したときのFを使って将来予測したときのサマリーをMSYのreference pointとする
                           MSY=as.data.frame(t(sapply(fout.HS,get.stat,eyear=eyear,hsp=b.HS))),
                           B0=as.data.frame(t(sapply(fout0.HS,get.stat,eyear=eyear,hsp=b.HS))),
                           per5=tmp)

#                   sl=list(a=a.SL,
#                       res=res.SL,jack=jack.SL,boot=boot.SL),
    if(what.est[2]==TRUE)
        output$bh <- list(a=a.BH,b=b.BH,sd=sd.BH,
                       res=res.BH,r0=NA,#R0を入れないといけない
                       Fhist=Fhist.BH,ofv=res.BH$value,
                       trace=trace.BH,b0=b0.BH,boot=as.data.frame(boot.BH),jack=as.data.frame(jack.BH),
                       f.msy=sapply(fout.BH,function(x)x$multi),
                       fout=fout.BH,
                       Fmsy=Fmsy.BH,farg=farg.BH,
                       MSY=as.data.frame(t(sapply(fout.BH,get.stat,eyear=eyear))),
                          B0=as.data.frame(t(sapply(fout0.BH,get.stat,eyear=eyear))))

    if(what.est[3]==TRUE)
        output$ri <- list(a=a.RI,b=b.RI,sd=sd.RI,ofv=res.RI$value,
                           res=res.RI,r0=NA,#R0を入れないといけない,
                           Fhist=Fhist.RI,farg=farg.RI,
                           trace=trace.RI,b0=b0.RI,boot=as.data.frame(boot.RI),
                          jack=as.data.frame(jack.RI),
                          fout=fout.RI,
                           f.msy=sapply(fout.RI,function(x)x$multi),
                       Fmsy=Fmsy.RI,
                     MSY=as.data.frame(t(sapply(fout.RI,get.stat,eyear=eyear))),
                          B0=as.data.frame(t(sapply(fout0.RI,get.stat,eyear=eyear))))
    index <- c("a","b","R0","sd","MSY","B0","f.msy","Fmsy")
    tmp <- NULL
    if(what.est[1]==TRUE) tmp <- rbind(tmp,unlist(output$hs[index]))
    if(what.est[2]==TRUE) tmp <- rbind(tmp,unlist(output$bh[index]))
    if(what.est[3]==TRUE) tmp <- rbind(tmp,unlist(output$ri[index]))

    tmp <- as.data.frame(tmp)
    rownames(tmp) <- c("hs","bh","ri")[what.est]
#    tmp$nLL <- output$ofv
    output$summary0 <- tmp
    colnames(output$summary0)[1] <- "a"
    output$summary <- output$summary0[c("a","b","sd","MSY.ssb.mean.ssb.mean",
                                        "MSY.biom.mean.biom.mean",
                                        "MSY.U.mean.U.mean",
                                        "MSY.catch.mean.catch.mean",
                                        "B0.ssb.mean.ssb.mean",
                                        "B0.biom.mean.biom.mean","f.msy")]
    colnames(output$summary) <- c("a","b","sd","SSB_MSY","B_MSY","U_MSY","MSY","B0(SSB)","B0(Biomass)","FMSY/Fcurrent")
    output$summary <- cbind(output$summary,output$summary0[,substr(colnames(output$summary0),1,4)=="Fmsy"])
    class(output) <- "SR"

    ##--- PGY管理基準値を計算する
    if(!is.null(PGY)){
        k.tmp <- which(what.est)
        for(k in 1:length(k.tmp)){
            fout.list2 <- list()
            s <- 1
            for(j in 1:length(PGY)){
                outtmp <- output[[which(names(output)==c("hs","bh","ri")[k.tmp[k]])[1]]]
#                outtmp$trace
#                frange.list <- list(c(output[[which(names(output)==c("hs","bh","ri")[1])[1]]]$f.msy,2),
#                                    c(0.01,output[[which(names(output)==c("hs","bh","ri")[1])[1]]]$f.msy))
                ttmp <- outtmp$trace[[1]]$catch.mean-PGY[j]*output$summary$MSY[k]
                ttmp <- which(diff(sign(ttmp))!=0)
                frange.list <- list(outtmp$trace[[1]]$fmulti[ttmp[1]+0:1],
                                   outtmp$trace[[1]]$fmulti[ttmp[2]+0:1])
#                browser()
                for(i in 1:2){
                    if(k.tmp[k]==1) farg.tmp <- farg.HS
                    if(k.tmp[k]==2) farg.tmp <- farg.BH
                    if(k.tmp[k]==3) farg.tmp <- farg.RI
                    farg.tmp$outtype <- NULL
                    farg.tmp$Frec <- list(stochastic=TRUE,
                                          future.year=rev(rownames(outtmp$fout[[1]]$vssb))[1],
                                          Blimit=PGY[j]*output$summary$MSY[k],
                                          scenario="catch.mean",Frange=frange.list[[i]])
                    fout.list2[[s]] <- do.call(future.vpa,farg.tmp)
                    s <- s+1
                }}
            PGY.biom <- as.data.frame(t(sapply(fout.list2,get.stat,eyear=eyear)))
            rownames(PGY.biom) <- paste("PGY",rep(PGY,each=2),rep(c("upper","lower"),length(PGY)),c("hs","bh","ri")[k.tmp[k]],sep="_")
            PGY.biom$target.catch <- rep(PGY*output$summary$MSY[k],each=2)
            if(k.tmp[k]==1) output$PGY.biom.hs <- PGY.biom
            if(k.tmp[k]==2) output$PGY.biom.bh <- PGY.biom
            if(k.tmp[k]==3) output$PGY.biom.ri <- PGY.biom
        }
    }
    ##---

    if(isTRUE(is.small)){
        output$hs$fout <- NULL
        output$bh$fout <- NULL
        output$ri$fout <- NULL
        }


    return(output)
}

pred.RI <- function(SSB,a,b) a*SSB*exp(-b*SSB)
pred.BH <- function(SSB,a,b) a*SSB/(1+b*SSB)
pred.HS <- function(SSB,a,b,gamma) a*(SSB+sqrt(b^2+gamma^2/4)-sqrt((SSB-b)^2+gamma^2/4))
pred.SL <- function(SSB,a) a*SSB
