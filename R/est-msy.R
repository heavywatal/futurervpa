#' Estimate MSY
#'
#' @rdname est-msy
#' @export
est.MSY <- function(vpares,farg,
                   seed=1,n.imputation=1,
                   nyear=50,
                   eyear=0, # 将来予測の最後のeyear+1年分を平衡状態とする
#                   FUN=median, # 漁獲量の何を最大化するか？
                   FUN=mean, # 漁獲量の何を最大化するか？
                   N=1000, # stochastic計算するときの繰り返し回数
                   onlylower.pgy=FALSE,# PGY計算するとき下限のみ計算する（計算時間省略のため）
                   is.small=FALSE, # 将来予測の結果を返さない。
                   is.Kobe=FALSE, # Kobeの計算をするかどうか。順番に、HS, BH, RIの順
                   is.5perlower=FALSE, # HSの折れ点を5%の確率で下回るときの親魚資源量
                   optim.method="optimize",
                   max.target="catch.mean", # method="optimize"以外を使うとき、どの指標を最大化するか。他のオプションとしては"catch.median" (漁獲量のmedianの最大化)
                   calc.yieldcurve=TRUE, # yield curveを正確に計算するかどうか。TRUEだと計算時間が余計にかかる。FALSEだと、yield curveは正確ではない
                   Blimit=0,
                   trace.multi=c(seq(from=0,to=0.9,by=0.1),1,seq(from=1.1,to=2,by=0.1),3:5,7,20,100), # Fmsyを探索したり、Yield curveを書くときにグリッドサーチをするときのFの刻み。Fcurrentに対する乗数。Fが異常に大きい場合、親魚=0になって加入＝NA
                   is.plot=TRUE,
                   PGY=NULL, # PGY管理基準値を計算するかどうか。計算しない場合はNULLを、計算する場合はc(0.8,0.9,0.95)のように割合を入れる
                   B0percent=NULL # B0_XX%の管理基準値を計算するかどうか
                   ){


    #-------------------- B0 & MSY for HS --------------------
    # 最小化のための関数
    # シミュレーション回数ぶんの漁獲量のFUN（mean, geomean, median）を最大化するFを選ぶ
    tmpfunc <- function(x,f.arg,FUN=FUN,eyear=eyear){
      f.arg$multi <- x
      fout <- do.call(future.vpa,f.arg)
      return(-FUN(fout$vwcaa[(nrow(fout$vwcaa)-eyear):nrow(fout$vwcaa),-1]))
    }

    tmpfunc2 <- function(x,f.arg,FUN=FUN,eyear=eyear,hsp=0){
      f.arg$multi <- x
      fout <- do.call(future.vpa2,f.arg)
      tmp <- as.numeric(fout$vssb[(nrow(fout$vssb)-eyear):nrow(fout$vssb),-1])
      lhs <- sum(tmp<hsp)/length(tmp)
      return( (lhs-0.05)^2 + as.numeric(lhs==0) + as.numeric(lhs==1)  )
    }

    #
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

    trace.func <- function(farg,eyear,hsp=0,trace.N=farg$N,
                           fmulti=c(seq(from=0,to=0.9,by=0.1),1,seq(from=1.1,to=2,by=0.1),3:5,7,20,100)){
        trace.res <- NULL
#        ssb.array <- array(0,dim=c(farg$nyear,farg$N+1,length(fmulti)))
        farg$outtype <- "FULL"
        farg$N <- trace.N
        for(i in 1:length(fmulti)){
            farg$multi <- fmulti[i]
            tmp <- do.call(future.vpa,farg)
#            ssb.array[,,i] <- tmp$vssb
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
        return(list(table=trace.res))
    }

    ######## 関数定義おわり
    b0 <- numeric() # B0
    fout <- fout0 <- trace <- Fhist <- fout.HS.5par <- list()

    farg.org <- farg.tmp <- farg
    farg.tmp$outtype <- "FULL"
    farg.tmp$nyear <- nyear
    farg.tmp$N <- N
#    trace.N <- ifelse(isTRUE(calc.yieldcurve),N,3)
    trace.N <- N
    farg.tmp$silent <- TRUE
    farg.tmp$is.plot <- FALSE

    # B0の計算
    farg.tmp$multi <- 0
    fout0 <- do.call(future.vpa,farg.tmp)
    B0 <- get.stat(fout0,eyear=eyear,hsp=Blimit)
    B0 <- cbind(B0,get.stat2(fout0,eyear=eyear,hsp=Blimit))
    rownames(B0) <- "B0"

    which.min2 <- function(x){
        max(which(min(x)==x))
    }

    trace <- trace.func(farg.tmp,eyear,hsp=Blimit,fmulti=trace.multi,trace.N=trace.N)

    xx <- which.max(trace$table$catch.mean)+c(-1,1)
    range.tmp <- trace$table$fmulti[xx]
    if(xx[1]==0) range.tmp <- c(0,range.tmp)
    if(is.na(range.tmp[2])) range.tmp[2] <- max(trace$table$fmulti)*10

    farg.tmp$multi <- 1
    cat("Estimating MSY\n")
    if(optim.method=="optimize"){
        tmp <- optimize(tmpfunc,range.tmp,f.arg=farg.tmp,eyear=eyear,FUN=FUN)
        # 壁にあたっている限り続ける
        while(sum(round(tmp$minimum,3)==range.tmp)>0){
            tmp0 <- round(tmp$minimum,3)==range.tmp
            range.tmp <- sort(c(range.tmp[tmp0],
                                range.tmp[tmp0] -2*(mean(range.tmp) - range.tmp[tmp0])))
            range.tmp <- ifelse(range.tmp<0,0,range.tmp)
            tmp <- optimize(tmpfunc,range.tmp,f.arg=farg.tmp,eyear=eyear,FUN=FUN)
        }
        farg.msy <- farg.tmp
        farg.msy$multi <- tmp$minimum # Fc.at.a * multiがFmsy
        cat("F multiplier=",tmp$minimum,"\n")
        fout.msy <- do.call(future.vpa,farg.msy)
        if(calc.yieldcurve){
            trace$table <- rbind(trace$table,trace.func(farg.msy,eyear,hsp=Blimit,trace.N=trace.N,
                                                    fmulti=tmp$minimum+c(-0.025,-0.05,-0.075,0,0.025,0.05,0.075))$table)
            trace$table <- trace$table[order(trace$table$fmulti),]
        }
    }
    # optimizeでなくgridでやる場合
    else{
        Fmulti <- seq(from=min(range.tmp),to=max(range.tmp),by=0.01)
        trace.tmp <- trace.func(farg.tmp,eyear,hsp=Blimit,fmulti=Fmulti,trace.N=trace.N)
        farg.msy <- farg.tmp
        farg.msy$multi <- trace.tmp$table$fmulti[which.max(unlist(trace.tmp$table[max.target]))]
        cat("F multiplier=",farg.msy$multi,"\n")
        fout.msy <- do.call(future.vpa,farg.msy)
        trace$table <- rbind(trace$table,trace.tmp$table)
        trace$table <- trace$table[order(trace$table$fmulti),]
    }

    MSY <- get.stat(fout.msy,eyear=eyear)
    MSY <- cbind(MSY,get.stat2(fout.msy,eyear=eyear))
    rownames(MSY) <- "MSY"
#    cat(" SSB=",MSY$"ssb.mean","\n")

    if(is.Kobe){
        cat("Estimating Kobe plot\n")
        Fhist <- get.Fhist(farg.org,vpares,eyear=eyear,trace=trace$table)
    }
    if(is.5perlower){
              xx <- which.min2((trace$table$lower-0.05)^2)+c(-1,1)
              range.tmp <- trace$table$fmulti[xx]
              if(xx[1]==0) range.tmp <- c(0,range.tmp)
              if(is.na(xx[2])) range.tmp[2] <- max(trace$table$fmulti)*10
              tmp <- optimize(tmpfunc2,range.tmp,f.arg=farg.org,eyear=eyear,FUN=FUN,hsp=Blimit)
              farg.tmp$multi <- tmp$minimum
              fout.HS.5per <- do.call(future.vpa,farg.tmp)
    }


    if(isTRUE(is.5perlower)){
        tmp <- as.data.frame(t(sapply(fout.HS.5per,get.stat,eyear=eyear,hsp=Blimit)))
        tmp$f <- sapply(fout.HS.5per,function(x)x$multi)
    }
    else{
        tmp <- NA
    }

    ## PGYの計算
    fout.PGY <- list()
    if(!is.null(PGY)){
        s <- 1
        for(j in 1:length(PGY)){
            cat("Estimating PGY ",PGY[j]*100,"%\n")
            ttmp <- trace$table$catch.mean-PGY[j]*MSY$catch.mean
            ttmp <- which(diff(sign(ttmp))!=0)
            frange.list <- list(trace$table$fmulti[ttmp[1]+0:1],
                                trace$table$fmulti[ttmp[2]+0:1])
            if(isTRUE(onlylower.pgy)) i.tmp <- 2  else i.tmp <- 1:2
            for(i in i.tmp){
                farg.pgy <- farg.tmp
                if(sum(is.na(frange.list[[i]]))>0) frange.list[[i]] <- c(0,300)
                farg.pgy$Frec <- list(stochastic=TRUE,
                                      future.year=rev(rownames(fout0$vssb))[1],
                                      Blimit=PGY[j]*MSY$catch.mean,
                                      scenario="catch.mean",Frange=frange.list[[i]])
                fout.PGY[[s]] <- do.call(future.vpa,farg.pgy)
                s <- s+1
            }
        }
        PGYstat <- as.data.frame(t(sapply(fout.PGY,get.stat,eyear=eyear,hsp=Blimit)))
        PGYstat <- cbind(PGYstat,as.data.frame(t(sapply(fout.PGY,get.stat2,eyear=eyear,hsp=Blimit))))
        rownames(PGYstat) <- names(fout.PGY) <- paste("PGY",rep(PGY,each=length(i.tmp)),
                                                      rep(c("upper","lower")[i.tmp],length(PGY)),sep="_")
    }
    else{
        PGYstat <-  NULL
        }
    ###

    ## B0_%の計算
    fout.list3 <- list()
    if(!is.null(B0percent)){
        for(j in 1:length(B0percent)){
            cat("Estimating B0 ",B0percent[j]*100,"%\n")
            ttmp <- trace$table$ssb.mean-B0percent[j]*B0$ssb.mean
            ttmp <- which(diff(sign(ttmp))!=0)
            frange.list <- trace$table$fmulti[ttmp[1]+0:1]
            farg.b0 <- farg.tmp
            farg.b0$Frec <- list(stochastic=TRUE,
                                 future.year=rev(rownames(fout0$vssb))[1],
                                 Blimit=B0percent[j]*B0$ssb.mean,
                                 scenario="ssb.mean",Frange=frange.list)
            fout.list3[[j]] <- do.call(future.vpa,farg.b0)
        }
        B0stat <- as.data.frame(t(sapply(fout.list3,get.stat,eyear=eyear,hsp=Blimit)))
        B0stat <- cbind(B0stat,as.data.frame(t(sapply(fout.list3,get.stat2,eyear=eyear,hsp=Blimit))))
        rownames(B0stat) <- names(fout.list3) <- paste("B0-",B0percent*100,"%",sep="")
    }
    else{
        B0stat <-  NULL
        }
    ###

    refvalue <- rbind(MSY,B0,PGYstat,B0stat)
    sumvalue <- refvalue[,c("ssb.mean","biom.mean","U.mean","catch.mean","Fref2Fcurrent")]
    colnames(sumvalue) <- c("SSB","B","U","Catch","Fref/Fcur")
    sumvalue <- cbind(sumvalue,refvalue[,substr(colnames(refvalue),1,1)=="F"])

    if(isTRUE(is.plot)){
        par(mfrow=c(1,2),mar=c(4,4,1,1))
        plot(trace$table$fmulti,trace$table$"ssb.mean"*1.2,type="n",xlab="Fref/Fcurrent",ylab="SSB")
        abline(v=sumvalue$Fref2Fcurrent,col="gray")
        text(sumvalue$Fref2Fcurrent,max(trace$table$"ssb.mean")*seq(from=1.1,to=0.8,length=nrow(sumvalue)),rownames(sumvalue))
        menplot(trace$table$fmulti,cbind(0,trace$table$"ssb.mean"),col="skyblue",line.col="darkblue")

        plot(trace$table$fmulti,trace$table$"catch.mean",type="n",xlab="Fref/Fcurrent",ylab="Catch")
        abline(v=sumvalue$Fref2Fcurrent,col="gray")
        menplot(trace$table$fmulti,cbind(0,trace$table$"catch.mean"),col="lightgreen",line.col="darkgreen")
    }

    ## kobe II matrix
    #kobe2 <- array(0,dim=c(dim(trace$array)[[1]],dim(trace$array)[[2]],length(sumvalue$SSb)))
    #for(i in 1:length(sumvalue$SSB)){
    #tmp <- trace$array > sumvalue$SSB[i]
    #kobe2[,,i] <- cbind(kobe2,apply(tmp,c(1,2),mean))
    #  }
    #dimnames(kobe2) <- list()

    invisible(list(all.stat=refvalue,summary=sumvalue,trace=trace$table,fout0=fout0,fout.msy=fout.msy,fout.B0percent=fout.list3,fout.PGY=fout.PGY))
}
