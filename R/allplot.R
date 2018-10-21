#' Plot all
#'
#' @description
#' 単位はcatch at ageの尾数が100万尾、waaがgの場合、重量の単位がちょうどトンになるようになっている。
#' @rdname allplot
#' @importFrom png readPNG
#' @export
allplot <- function(res0,target="hs",biomass.scale=1000,
                    pngfile="../buri.png",detail.plot=1:3){
    dmodel <- res0$vpares
    summary <- res0$summary[rownames(res0$summary)==target,]
    res00 <- res0[names(res0)==target][[1]]
    if(is.null(dmodel$Bban))  dmodel$Bban <- NA
    if(is.null(dmodel$Blimit))  dmodel$Blimit <- NA
    if(is.null(dmodel$ABC))  dmodel$ABC <- NA
#    summary$Bban <- dmodel$Bban
    summary[1] <- res00$r0
    colnames(summary)[1] <- "R0"
    summary$Blimit <- ifelse(length(dmodel$Blimit)>0,dmodel$Blimit,NA)
    summary$"Fmsy/Fcurrent" <- res00$f.msy
    summary$ABC2017 <- dmodel$ABC2017
    summary$ABC2015 <- dmodel$ABC2015
    summary$"newABC2015" <- rev(colSums(dmodel$baa))[1] * summary$"U_MSY"

    summary$"Pr(SSB_2021>SSB_MSY) (%)" <- mean(res00$fout[[1]]$vssb[6,]>summary$"SSB_MSY")*100
    summary$"Pr(SSB_2021>SSB_hs) (%)" <- mean(res00$fout[[1]]$vssb[6,]>summary$"b")*100

    colnames(summary)[2] <- "SSB_HS"
    colnames(summary)[8:9] <- c("SSB0","B0")


    summary1 <- data.frame(parameter=names(summary),value=as.numeric(summary))

    ## currentF projection
    arg.tmp <- res0[names(res0)==target][[1]]$farg
    arg.tmp$multi <- 1
    fout0 <- do.call(future.vpa2,arg.tmp) # currentFでの将来予測の結果
    fout1 <- res0[names(res0)==target][[1]]$fout[[1]] # Fmsyでの将来予測の結果

    layout(t(matrix(c(1,2,3,4,5,6,7,8),2,4)),heights=c(0.7,1,1,1))

    if(1%in%detail.plot){
        ## plot talbe
        par(mar=c(1,4.3,3,1))
        n <- floor(nrow(summary1)/2)
        plot.info(summary1[1:n,])
        title(dmodel$jstockid)
        plot.info(summary1[(n+1):nrow(summary1),])

        ## SR plot(all)
        par(mar=c(4.3,4.3,3,1))
        aa <- summary[c("SSB_MSY","Blimit","SSB_HS")]
        fit.tmp <- plot.SR(res0,what.plot=rownames(res0$summary),pick="SSB_MSY",xyscale=c(1.8,1.3))
#        abline(v=aa,col=c("chartreuse3","orange","red"))
        title("R vs SSB (HS, BH, RI)")

        ## SR plot(HS)
        par(mar=c(4.3,4.3,3,1))
        aa <- summary[c("SSB_MSY","Blimit","SSB_HS")]
        fit.tmp <- plot.SR(res0,what.plot=target,pick="SSB_MSY",xyscale=c(1.8,1.3),is.legend=FALSE)
#        abline(v=aa,col=c("chartreuse3","orange","red"))
        title(paste("R vs SSB (only ",target,")",sep=""))

        ## Residual plot
        plot(x <- res0$dat$year,
             y <- fit.tmp$resid[[1]],
             ylim=c(-1.5,1.5),type="p",pch=20,xlab="Year",ylab="log(Obs)-log(Pred)")
        abline(h=0,lty=2)
        xx <- loess(y~x)
        points(x,xx$fitted,type="l",col=2,lwd=2)
        title("Residual to HS prediction")

        ## plot diagnostics
        if(target=="hs"){
            show.likeprof(res0)
            lines(c(res0$hs$b,res0$hs$b),quantile(res0$hs$boot$a,probs=c(0.05,0.95)))
            lines(quantile(res0$hs$boot$b,probs=c(0.05,0.95)),c(res0$hs$a,res0$hs$a))
            points(res0$hs$jack$b,res0$hs$jack$a)
            legend("topleft",lty=c(1,NA),pch=c(NA,1),legend=c("Bootstrap 90%", "Jackknife"),ncol=2,bg="white")
        }

        ## yield curve
        par(mar=c(4.3,4.3,3,4.3))
        plotyield(res00)

    }

    if(2%in%detail.plot){
        layout(t(matrix(c(1,2,3,4,5,6),2,3)),heights=c(1.2,1,1))
                                        # Kobe
        par(mar=c(4.3,4.3,5,1))
        a <- plot.Kobe2(res0,SR=target,pickU="U_MSY",pickB="SSB_MSY",is.plot=T)
        abline(v=summary$"SSB_HS"/summary$"SSB_MSY",lty=2)
        title("Kobe chart",line=0.5)
        title(dmodel$jstockid,line=2)

        # plot selectivity
        matplot(rownames(res0$vpares$faa),
                res0$vpares$faa,pch=20,col="gray",xlab="Age",ylab="F",ylim=c(0,max(res0$vpares$faa,na.rm=T)))
        points(rownames(res0$vpares$faa),res00$fout[[1]]$faa[,1],type="b")
        legend("topleft",pch=c(1,20),legend=c("Current F","Past Fs"),col=c(1,"gray"))
        title("Current F")

        ## plot SSB
        par(mar=c(4.3,4.3,3,1))
        years <- as.numeric(colnames(dmodel$ssb))
        y <- as.numeric(colSums(dmodel$ssb))
        plot(years,y,ylim=c(0,1.1*max(c(y,unlist(aa)),na.rm=T)),xlab="Years",ylab="SSB",type="o",
             xlim=c(min(years),max(years)+10))

        ## projection
        menplot(rownames(fout0$vssb),t(apply(fout0$vssb,1,quantile,probs=c(0.1,0.9))),
                col=rgb(40/255,96/255,163/255,0.2),border="blue",lty=2)
        menplot(rownames(fout1$vssb),t(apply(fout1$vssb,1,quantile,probs=c(0.1,0.9))),
                col=rgb(40/255,96/255,163/255,0.2),border="red",lty=2)
        points(rownames(fout0$vssb),apply(fout0$vssb,1,mean),type="o",pch=20,
               col="blue")
        points(rownames(fout1$vssb),apply(fout1$vssb,1,mean),type="o",pch=20,
               col=2)

        abline(h=aa,col=c("chartreuse3","orange","red"))
        legend("topleft",col=c("chartreuse3","orange","red","blue","red"),lty=1,pch=c(NA,NA,NA,20,20),
               legend=c("SSB_MSY","SSB_limit","SSB_HS","Projection (Fcur)","Projection (Fmsy)"),cex=0.8,
               bg="white",ncol=2)
        title("SSB",line=0.5)

                                        # plot catch
#        y <- as.numeric(colSums(dmodel$input$dat$caa * dmodel$input$dat$waa,na.rm=T)) # * N.unit * waa.unit
        y <- as.numeric(colSums(dmodel$wcaa))
        aa <- summary[c("MSY")]
        plot(years,y,ylim=c(0,1.1*max(c(y,unlist(aa)))),xlab="Years",ylab="Total catch",type="o",
             xlim=c(min(years),max(years)+10))

        ## projection
        menplot(rownames(fout0$vssb),t(apply(fout0$vwcaa,1,quantile,probs=c(0.1,0.9))),
                col=rgb(210/255,94/255,44/255,0.3),border="blue",lty=2)
        menplot(rownames(fout1$vssb),t(apply(fout1$vwcaa,1,quantile,probs=c(0.1,0.9))),
                col=rgb(210/255,94/255,44/255,0.3),border="red",lty=2)
        points(rownames(fout0$vssb),apply(fout0$vwcaa,1,mean),type="o",pch=20,
               col="blue")
        points(rownames(fout1$vssb),apply(fout1$vwcaa,1,mean),type="o",pch=20,
               col=2)

        abline(h=aa,col="chartreuse3")
        tmp <- names(summary)%in%c("ABC2017","ABC2015","newABC2015")
        tmp2 <- c("ABC2017","ABC2015","newABC2015")%in%names(summary)
        points(c(2017,2015,2015)[tmp2],summary[tmp],pch=c(3,2,1)[tmp2],col=c(1,1,2)[tmp2])
        points(c(2017,2015,2015)[tmp2],summary[tmp],pch=c(3,2,1)[tmp2],col=c(1,1,2)[tmp2])
        tmp3 <- c(TRUE,tmp2,rep(TRUE,3))
        legend("topleft",col=c("chartreuse3",1,1,2,"blue","red")[tmp3],
               bg="white",
               pch=c(NA,3,2,1,20,20)[tmp3],legend=c("MSY","ABC2017","ABC2015","newABC2015","Projection (Fcur)","Projection (Fmsy)")[tmp3],ncol=2,cex=0.8,
               lty=c(1,NA,NA,NA,1,1)[tmp3])



        title("Catch")

                                        # plot exploitation rates
        y <- y/as.numeric(colSums(dmodel$baa))
        aa <- summary[c("U_MSY")]
        plot(years,y,ylim=c(0,1.1*max(c(y,unlist(aa)))),xlab="Years",ylab="Exploitation rates (Catch/B)",type="o",
             xlim=c(min(years),max(years)+10))

        # projection
        points(rownames(fout0$vssb),apply(fout0$vwcaa,1,mean)/apply(fout0$vbiom,1,mean),type="o",pch=20,
               col="blue")
        points(rownames(fout1$vssb),apply(fout1$vwcaa,1,mean)/apply(fout1$vbiom,1,mean),type="o",pch=20,
               col=2)


        abline(h=aa,col=c("chartreuse3"))
        legend("topright",col=c("chartreuse3"),lty=1,legend=c("U_MSY"))
        title("Exploitation rates (U)")

        ## plot %SPR
        #y <- dmodel$SPR$ysdata$perSPR
        y <- get.SPR(res0)$ysdata$perSPR
        plot(years,y,ylim=c(0,100),xlab="Years",ylab="%SPR",type="o",
             xlim=c(min(years),max(years)))

        aa <- summary[c("SSB_MSY","SSB0")]
        abline(h=aa[1]/aa[2]*100,col=c("chartreuse3"))
        legend("topright",col=c("chartreuse3"),lty=1,legend=c("SSB_MSY/SSB0"))
        title("%SPR")

    }

    if(3%in%detail.plot){
        layout(t(matrix(c(1,2,3,3,4,4,5,5),2,4)),heights=c(1,2,1,1))
        tres0 <- res00$trace[[1]]
        ssb <- res00$trace[[1]]$ssb.mean/biomass.scale
        tmp <- substr(colnames(tres0),1,5)=="TB-MA"
        tb <- tres0[,tmp]/biomass.scale
        tb2 <- sapply(1:ncol(tb),function(x) apply(tb[,1:x,drop=F],1,sum,na.rm=T))
        tmp <- substr(colnames(tres0),1,5)=="TC-MA"
        tc <- tres0[,tmp]/biomass.scale
        tc2 <- sapply(1:ncol(tc),function(x) apply(tc[,1:x,drop=F],1,sum,na.rm=T))
        if(file.exists(pngfile)) image <- readPNG(pngfile)
        else image <- NULL

#        if(ncol(tc)<10)    col.tmp <- brewer.pal(ncol(tc),"Greens")
#        else col.tmp <- c(brewer.pal(9,"Greens"),rev(brewer.pal(ifelse(ncol(tc)-9<3,3,ncol(tc)-9),"GnBu")))
        col.tmp1 <- rgb(40/255,96/255,163/255,seq(from=0.1,to=0.9,length=ncol(tc)))
        col.tmp2 <- rgb(210/255,94/255,44/255,seq(from=0.1,to=0.9,length=ncol(tc)))

        ## plot table
        par(mar=c(1,4.3,3,2))
        plot.info(summary1[c(2,4,5,6,10,11,12,14),])
        title(dmodel$jstockid)

        plot.SR(res0,what.plot=target,pick="SSB_MSY",xyscale=c(1.8,1.3),
                is.legend=FALSE)

        par(mar=c(2,4.3,3,4.3))
        scale <- max(ssb,na.rm=T) * 0.1
                                        #    ysize <- max(ssb,na.rm=T)/max(tb,na.rm=T)
        year.tmp <- rev(colnames(res0$vpares$ssb))[1:5]
        range1 <- range(res0$vpares$ssb)/biomass.scale
        range2 <- range(as.data.frame(res0$vpares$ssb)[as.character(year.tmp)])/biomass.scale


        ### plot of SSB
        ssb.max <- max(c(range1,summary$"SSB_MSY"),na.rm=T) *1.5 /biomass.scale
        tb3 <- tb2[which(ssb<ssb.max),]
        matplot(ssb,tb2,type="n",xlab="",ylab=paste("Biomass (",biomass.scale,")",sep=""),xaxs="i",yaxs="i",
                ylim=c(0,max(tb2[which(ssb<ssb.max),])*1.2),xlim=c(0,ssb.max))
                                        #    menplot(range1,cbind(c(-100,-100),rep(max(tb2)*1.5,2)),col=gray(0.9),border=NA)
                                        #    menplot(range2,cbind(c(-100,-100),rep(max(tb2)*1.5,2)),col=gray(0.7),border=NA)
        non.na <- !is.na(ssb)
        for(i in 1:ncol(tb2)) menplot(ssb[non.na], cbind(0,tb2)[non.na,i:(i+1)],col=col.tmp1[i],border=NA)
        waa.tmp <- (apply(res0$vpares$input$dat$waa,1,mean))^(1/3)*10
        waa.tmp <- waa.tmp/max(waa.tmp) * 0.9
        x <- tb3[1,]
        if(!is.null(image)){
            plotfish(image,x=rep(ssb.max*0.88,ncol(tb)),y=x-diff(c(0,x))/2,
                     size=waa.tmp,scale=scale,ysize=1)
        }
        text(rep(ssb.max*0.93,ncol(tb)),x-diff(c(0,x))/2,
             paste(0:(ncol(tb2)-1),"y/o: ",round(apply(res0$vpares$input$dat$waa,1,mean),0),""),cex=1)

        matpoints(ssb,tb2[,1],type="l",lwd=2,col="gray",lty=3)
        points(x <- colSums(res0$vpares$ssb)/biomass.scale,
               y <- colSums(res0$vpares$baa)/biomass.scale,type="o",
               col=gray(c(seq(from=0.7,to=0,length=length(x)))),pch=20,cex=1.2,
               lwd=3)
        text(x[1],y[1],colnames(x)[1],adj=0)
        text(rev(x[1]),rev(y)[1],rev(colnames(x))[1],adj=0)
                                        #    browser()
        abline(v=summary$"SSB_MSY"/biomass.scale,lty=2)
        abline(v=summary$"Blimit"/biomass.scale,lty=2)
        abline(v=summary$"SSB_HS"/biomass.scale,lty=2)
        text(summary$"SSB_MSY"/biomass.scale,max(tb3)*1.1,
             paste("SSB_MSY=",format(round(summary$"SSB_MSY"/biomass.scale),big.mark=","),"",sep=""),adj=0)
        text(summary$"Blimit"/biomass.scale,max(tb3)*1.0,
             paste("SSB_limit=",format(round(summary$"Blimit"/biomass.scale),big.mark=","),
                   "",sep=""),adj=0)
        text(summary$"SSB_HS"/biomass.scale,max(tb3)*1.05,
             paste("SSB_HS=",format(round(summary$"SSB_HS"/biomass.scale),big.mark=","),
                   "",sep=""),adj=0)
                                        #    text(max(ssb)*0.8,max(tb[,1],na.rm=T)*1.05,"予測加入量",col=2)
        ## plot of SSB CV
        par(new=T)
        y <- res00$trace[[1]]$ssb.CV
        plot(ssb,y,type="l",lwd=3,col=rgb(0.8,0.8,0,0.6),axes=F,xlab="",ylab="",
             ylim=c(0,ifelse(max(y,na.rm=T)>1.5,1.5,max(y,na.rm=T))))
        axis(side=4)
        mtext(side=4,"CV of SSB",line=3,col=rgb(0.8,0.8,0,0.6),cex=0.8)

        ##  catch
        par(mar=c(2,4.3,1,4.3))
        if(!is.null(res0$vpares$wcaa)) wcatch <- as.numeric(colSums(res0$vpares$wcaa))
        else{
            wcatch <- as.numeric(colSums(res0$vpares$input$dat$caa * res0$vpares$input$dat$waa/biomass.scale,na.rm=T))
        }
        matplot(ssb,tc2,type="n",,xaxs="i",yaxs="i",ylab=paste("Catch (",biomass.scale,")",sep=""),
                                        #            ylim=c(0,max(tc2,wcatch)*1.2),xlim=c(0,ssb.max))
                ylim=c(0,max(tc2)*1.2),xlim=c(0,ssb.max))

                                        #    menplot(range1,cbind(c(-100,-100),rep(max(tb2)*1.5,2)),col=gray(0.9),border=NA)
                                        #    menplot(range2,cbind(c(-100,-100),rep(max(tb2)*1.5,2)),col=gray(0.7),border=NA)
        scale <- max(ssb)/max(tc2)
        for(i in 1:ncol(tc2)) menplot(ssb[non.na], cbind(0,tc2)[non.na,i:(i+1)],col=col.tmp2[i],border=NA)
        abline(v=summary$"SSB_MSY"/biomass.scale,lty=2)
        abline(v=summary$"Blimit"/biomass.scale,lty=2)
        abline(v=summary$"SSB_HS"/biomass.scale,lty=2)
        points(x <- as.numeric(colSums(res0$vpares$ssb))/biomass.scale,
               y <- wcatch,pch=20,lwd=3,
               type="o",col=gray(c(seq(from=0.7,to=0,length=length(x)))))
        text(x[1],y[1],colnames(res0$vpares$ssb)[1],adj=0)
        text(rev(x)[1],rev(y)[1],rev(colnames(res0$vpares$ssb))[1],adj=0)
        points(x <- apply(res00$fout[[1]]$vssb,1,mean)[1:10]/biomass.scale,
               y <- apply(res00$fout[[1]]$vwcaa,1,mean)[1:10]/biomass.scale,col=2,
               type="o",pch=20,lwd=3)
        text(rev(x)[1],rev(y)[1],
             paste("Projection ",rownames(res00$fout[[1]]$vssb)[10],"(F_MSY)",sep=""),adj=-0.1,col=2)

        points(x <- apply(fout0$vssb,1,mean)[1:10]/biomass.scale,
               y <- apply(fout0$vwcaa,1,mean)[1:10]/biomass.scale,col="blue",type="o",pch=20,lwd=3)
        text(rev(x)[1],rev(y)[1],
             paste("Projection ",rownames(res00$fout[[1]]$vssb)[10],"(F_current)",sep=""),adj=-0.1,col="blue")

        ## plot of catch CV
        par(new=T)
        y <- res00$trace[[1]]$catch.CV
        plot(ssb,y,type="l",lwd=2,col=rgb(0.8,0.8,0,0.6),axes=F,xlab="",ylab="",
             ylim=c(0,ifelse(max(y,na.rm=T)>1.5,1.5,max(y,na.rm=T))))
        axis(side=4)
        mtext(side=4,"CV of Catch",line=3,col=rgb(0.8,0.8,0,0.6),cex=0.8)

        ## 努力量 plot
        par(mar=c(4.3,4.3,2,4.3))
        tmp <- round(ssb*biomass.scale)>0 & !is.na(ssb)
        matplot(ssb,tres0$fmulti,type="n",ylab="Fishing efforts (current=1)",xaxs="i",yaxs="i",xlab=paste("SSB(",biomass.scale,")",sep=""),xlim=c(0,ssb.max),
                ylim=c(0,max(tres0$fmulti[tmp]*1.2)))
                                        #    menplot(range1,cbind(c(-100,-100),rep(max(tb2)*1.5,2)),col=gray(0.9),border=NA)
                                        #    menplot(range2,cbind(c(-100,-100),rep(max(tb2)*1.5,2)),col=gray(0.7),border=NA)
        menplot(ssb[tmp],cbind(0,tres0$fmulti[tmp]),col=rgb(221/255,159/255,33/255,0.5),border=NA)
        if(length(res0$hs$Fhist)>0){
            points(x <- unlist(colSums(res0$vpares$ssb))/biomass.scale,
                   y <- res0$hs$Fhist[[1]]$fmulti,pch=20,lwd=3,
                   type="o",col=gray(c(seq(from=0.7,to=0,length=length(x)))))
        }
        abline(h=1,lty=2)
        abline(v=summary$"SSB_MSY"/biomass.scale,lty=2)
        abline(v=summary$"Blimit"/biomass.scale,lty=2)
        abline(v=summary$"SSB_HS"/biomass.scale,lty=2)
                                        #    points(x <- ssb[which.min(abs(tres0$fmulti-1))],y <- tres0$fmulti[which.min(abs(tres0$fmulti-1))],pch=4)
                                        #    text(x,y,"Recent 3 years",adj=-0.1)

    }
    else{
        a <- NULL
    }

    Fcurrent <- list(wcatch=mean(fout0$vwcaa[nrow(fout0$vwcaa),],na.rm=T),
                     ssb=mean(fout0$vssb[nrow(fout0$vssb),],na.rm=T))

    invisible(list(a=a,summary=summary,Fcurrent=Fcurrent))

}



show.likeprof <- function(res){
    x <- tapply(res$hs$surface$obj,list(res$hs$surface$b,res$hs$surface$a),function(x) x)
    image(as.numeric(rownames(x)),as.numeric(colnames(x)),log(x/min(x)),col=rev(heat.colors(12)),ylab="a",xlab="b")
    contour(as.numeric(rownames(x)),as.numeric(colnames(x)),log(x/min(x)),add=T,nlevels=10,zlim=c(0,0.3))
    points(res$hs$b,res$hs$a)
    title("Diagnostics")
}

plot.info <- function(a,xpos=7){
    plot(1:(nrow(a)+2),type="n",ylab="",xlab="",axes=F)
    units <- ceiling(-1*log10(a[,2]))
    units <- units + 2
    units <- ifelse(units<0,0,units)
    for(i in 1:nrow(a)){
      text(1,nrow(a)-i+2,a[i,1],adj=c(0,1),cex=1)
      text(xpos,nrow(a)-i+2,format(round(a[i,2],units[i]),big.mark=",",
                                   scientific=F),adj=c(1,1))
    }
}

plotfish <- function(image,x,y,size,scale=1,ysize=1){
#    image <- readJPEG("../buri.jpg")
    xx <- dim(image)[1]/dim(image)[2]
    rasterImage(image,
                x-size*xinch(1), y-size*yinch(1)*xx*ysize, x+size*xinch(1), y+size*yinch(1)*xx*ysize)
}

menplot <- function(x,y,line.col=1,...){
    polygon(c(x,rev(x)),c(y[,1],rev(y[,2])),...)
    if(dim(y)[[2]]>2) points(x,y[,3],type="l",lwd=2,col=line.col)
}

menplot2 <- function(xy,probs=c(0.1,0.9),new=FALSE,xlab=NULL,ylab=NULL,...){
    xx <- rownames(xy)
    yy <- t(apply(xy,1,quantile,probs=c(0.1,0.9)))
    if(isTRUE(new)) matplot(xx,yy,type="n",xlab=xlab,ylab=ylab)
    menplot(xx,yy,...)
}


get.SPR <- function(dres){
    # Fの歴史的な%SPRを見てみる
    # 毎年異なるFや生物パラメータに対して、YPR,SPR、SPR0がどのくらい変わっているのか見る(Rコード例2)
    dres$ysdata <- matrix(0,ncol(dres$faa),4)
    dimnames(dres$ysdata) <- list(colnames(dres$faa),c("perSPR","YPR","SPR","SPR0"))
    for(i in 1:ncol(dres$faa)){
	dres$Fc.at.age <- dres$faa[,i] # Fc.at.ageに対象年のFAAを入れる
        if(all(dres$Fc.at.age>0)){
            byear <- colnames(dres$faa)[i] # 何年の生物パラメータを使うか
        # RVPAのref.F関数でYPRなどを計算。
        # 配布している1.3から1.4にアップデートしているので、新しいほうの関数を使うこと(返り値がちょっと違う)
            a <- ref.F(dres,waa.year=byear,maa.year=byear,M.year=byear,rps.year=2000:2011,
                       F.range=c(seq(from=0,to=ceiling(max(dres$Fc.at.age,na.rm=T)*2),
                                     length=101),max(dres$Fc.at.age,na.rm=T)),plot=FALSE)
        # YPRと%SPR
            dres$ysdata[i,1:2] <- (as.numeric(rev(a$ypr.spr[nrow(a$ypr.spr),-1])))
        # SPR
        dres$ysdata[i,3] <- a$spr0*dres$ysdata[i,1]/100
        # SPR0
        dres$ysdata[i,4] <- a$spr0
        }
        else{
            break;
            }
    }
    dres$ysdata <- as.data.frame(dres$ysdata)
    return(dres)
}
