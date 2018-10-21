## kobe.matrixの計算
# Pr(B<Btarget)のみを返す単純なやつ

get.kobemat <- function(fout,N=fout$input$N,nyear=fout$input$nyear,Btarget=0,
                      fmulti=seq(from=0.3,to=1,by=0.1)){
    multi.org <- 1
    fres.short <- list()
    farg <- fout$input
    farg$Frec <- NULL
    farg$N <- N
    farg$nyear <- nyear
    for(i in 1:length(fmulti)){
        farg$multi <- multi.org * fmulti[i]
        fres.short[[i]] <- do.call(future.vpa,farg)
    }
  	prob.btarget <- sapply(fres.short,function(x) apply(x$vssb>Btarget,1,mean))
	colnames(prob.btarget) <- fmulti

    invisible(prob.btarget)
}

# Btargetはベクトルで入力、平均親魚量なども出力
get.kobemat2 <- function(fout,N=fout$input$N,nyear=fout$input$nyear,Btarget=0,
                      fmulti=seq(from=0.3,to=1,by=0.1),target.name=1:length(Btarget)){
    multi.org <- 1
    fres.short <- list()
    farg <- fout$input
    farg$Frec <- NULL
    farg$N <- N
    farg$nyear <- nyear
    for(i in 1:length(fmulti)){
        farg$multi <- multi.org * fmulti[i]
        fres.short[[i]] <- do.call(future.vpa,farg)
    }

    # 結果の取り出し
    prob.btarget <- list()
    for(i in 1:length(Btarget)){
		prob.btarget[[i]] <- sapply(fres.short,function(x) apply(x$vssb>Btarget,1,mean))
		colnames(prob.btarget[[i]]) <- fmulti
      }
	names(prob.btarget) <- target.name

	# SSB, biomass, catch
    ssb <- sapply(fres.short,function(x) apply(x$vssb,1,mean))
	biom <- sapply(fres.short,function(x) apply(x$vbiom,1,mean))
	catch <- sapply(fres.short,function(x) apply(x$vwcaa, 1,mean))
	colnames(ssb) <- colnames(biom) <- colnames(catch) <- fmulti

    invisible(list(prob.btarget=prob.btarget,ssb=ssb,biom=biom,catch=catch))
}

## ちょっと複雑なkobe.plot
# foutsが複数の将来予測の結果。brefsは複数の管理基準値
#' @importFrom foreach foreach %do%
#' @importFrom RColorBrewer brewer.pal
get.kobemat2 <- function(fouts,brefs,xlim=NULL,target.prob=0.5){
#    brefs <- sort(brefs)
    years <- as.numeric(rownames(fouts[[1]]$vssb))
    probs <- matrix(0,length(years),length(fouts))
    for(j in 1:ncol(brefs)){
        probs <- probs + foreach(i=1:length(fouts),.combine=cbind) %do%
            as.numeric(rowMeans(fouts[[i]]$vssb > brefs[i,j])>target.prob)
        }
    if(is.null(xlim)) xlim <- range(years)
    plot(range(years),
         range(0.5,nrow(brefs)+1.5),
         type="n",xlab="Years",cex=3,xlim=xlim,
         ylab="Strategies",yaxt="n")
    abline(h=1:ncol(brefs),v=years,col="gray")
    axis(side=2,at=1:nrow(brefs),label=rownames(brefs))

    cols <- brewer.pal(ncol(brefs), "Paired")

    for(i in 1:length(fouts)){
        points(years,rep(i,length(years)),
               col=cols[probs[,i]],pch=20,cex=3)
    }
    legend("topright",pch=20,cex=1,col=cols,ncol=ceiling(ncol(brefs)/2),
           legend=paste("Prob(B>",colnames(brefs),")>",round(target.prob*100),"%"))
}

get.kobematrix <- function(fres,Blim=0,Bban=0,ssb=TRUE){
    if(isTRUE(ssb))  tmp <- fres$vssb[,-1]
    else  tmp <- fres$vbiom[,-1]

    res <- data.frame(
        # 漁獲量
        catch.deterministic=fres$vwcaa[,1],
        # 資源量
        biom.deterministic=fres$vbiom[,1],
        # 親魚量
        ssb.deterministic=fres$vssb[,1],
        # Blim回復確率
        probability.upper.Blim=apply(tmp>Blim,1,mean)*100,
        # Bban以上確率
        probability.upper.Bban=apply(tmp>Bban,1,mean)*100)

    return(res)
}
