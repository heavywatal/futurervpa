#' 資源量の上積みグラフを書く
#'
#' @rdname plot-b-fish
#' @importFrom png readPNG
#' @export
plotBfish <- function(res0, # SR.estの結果
                      Bref,
                      unit.waa=1,ssb.max=Inf,
                      target="hs",biomass.scale=1000,pngfile="fish.png"){


    summary <- res0$summary[rownames(res0$summary)==target,]
    res00 <- res0[names(res0)==target][[1]]
    tres0 <- res00$trace[[1]]
    ssb <- res00$trace[[1]]$ssb.mean/biomass.scale

    tmp <- substr(colnames(tres0),1,5)=="TB-MA"
    tb <- tres0[,tmp]/biomass.scale * unit.waa
    tb2 <- sapply(1:ncol(tb),function(x) apply(tb[,1:x,drop=F],1,sum,na.rm=T))

    tmp <- substr(colnames(tres0),1,5)=="TC-MA"
    tc <- tres0[,tmp]/biomass.scale * unit.waa
    tc2 <- sapply(1:ncol(tc),function(x) apply(tc[,1:x,drop=F],1,sum,na.rm=T))
    if(file.exists(pngfile)) image <- readPNG(pngfile)
    else image <- NULL

    year.tmp <- rev(colnames(res0$vpares$ssb))[1:5]
    range1 <- range(res0$vpares$ssb)/biomass.scale
    range2 <- range(as.data.frame(res0$vpares$ssb)[as.character(year.tmp)])/biomass.scale

    col.tmp1 <- rgb(40/255,96/255,163/255,seq(from=0.1,to=0.9,length=ncol(tc)))
    col.tmp2 <- rgb(210/255,94/255,44/255,seq(from=0.1,to=0.9,length=ncol(tc)))

    ### plot of SSB
    ssb.max <- min(ssb.max,
                   max(c(range1,summary$"SSB_MSY"),na.rm=T)) *1.5 /biomass.scale
    tb3 <- tb2[which(ssb<ssb.max),]
    matplot(ssb,tb2,type="n",ylab=paste("Biomass (",biomass.scale," MT)",sep=""),xaxs="i",yaxs="i",
            xlab="SSB",
            ylim=c(0,max(tb2[which(ssb<ssb.max),])*1.2),xlim=c(0,ssb.max))
#            ylim=c(0,max(tb2)),xlim=c(0,ssb.max))
                                        #    menplot(range1,cbind(c(-100,-100),rep(max(tb2)*1.5,2)),col=gray(0.9),border=NA)
                                        #    menplot(range2,cbind(c(-100,-100),rep(max(tb2)*1.5,2)),col=gray(0.7),border=NA)
    # 管理基準値のプロット
    plot.RP(Bref,biomass.scale=biomass.scale,ymax=max(tb3)*1.1)

    # 過去の時系列
#        matpoints(ssb,tb2[,1],type="l",lwd=2,col="gray",lty=3)
        points(x <- colSums(res0$vpares$ssb)/biomass.scale,
               y <- colSums(res0$vpares$baa)/biomass.scale,type="o",
               col=gray(c(seq(from=0.7,to=0,length=length(x)))),pch=20,cex=1.2,
               lwd=3)
        text(x[1],y[1],colnames(x)[1],adj=0)
        text(rev(x[1]),rev(y)[1],rev(colnames(x))[1],adj=0)

    ## 積み上げグラフ
    non.na <- !is.na(ssb)
    for(i in 1:ncol(tb2)) menplot(ssb[non.na], cbind(0,tb2)[non.na,i:(i+1)],col=col.tmp1[i],border=NA)
    title("Total biomass",line=-1,adj=0.1)

                                        #    browser()
        ## abline(v=summary$"SSB_MSY"/biomass.scale,lty=2)
        ## abline(v=summary$"Blimit"/biomass.scale,lty=2)
        ## abline(v=summary$"SSB_HS"/biomass.scale,lty=2)
        ## text(summary$"SSB_MSY"/biomass.scale,max(tb3)*1.1,
        ##      paste("SSB_MSY=",format(round(summary$"SSB_MSY"/biomass.scale),big.mark=","),"",sep=""),adj=0)
        ## text(summary$"Blimit"/biomass.scale,max(tb3)*1.0,
        ##      paste("SSB_limit=",format(round(summary$"Blimit"/biomass.scale),big.mark=","),
        ##            "",sep=""),adj=0)
        ## text(summary$"SSB_HS"/biomass.scale,max(tb3)*1.05,
        ##      paste("SSB_HS=",format(round(summary$"SSB_HS"/biomass.scale),big.mark=","),
        ##            "",sep=""),adj=0)


        ##  catch
        if(!is.null(res0$vpares$wcaa)) wcatch <- as.numeric(colSums(res0$vpares$wcaa))
        else{
            wcatch <- as.numeric(colSums(res0$vpares$input$dat$caa * res0$vpares$input$dat$waa,na.rm=T))*unit.waa
        }
        matplot(ssb,tc2,type="n",,xaxs="i",yaxs="i",ylab=paste("Catch (",biomass.scale,") MT",sep=""),
                xlab="SSB",
                                        #            ylim=c(0,max(tc2,wcatch)*1.2),xlim=c(0,ssb.max))
                ylim=c(0,max(tc2)*1.2),xlim=c(0,ssb.max))

                                        #    menplot(range1,cbind(c(-100,-100),rep(max(tb2)*1.5,2)),col=gray(0.9),border=NA)
                                        #    menplot(range2,cbind(c(-100,-100),rep(max(tb2)*1.5,2)),col=gray(0.7),border=NA)
    points(x <- as.numeric(colSums(res0$vpares$ssb))/biomass.scale,
           y <- wcatch/biomass.scale,pch=20,lwd=3,
           type="o",col=gray(c(seq(from=0.7,to=0,length=length(x)))))
    plot.RP(Bref,biomass.scale=biomass.scale,ymax=max(tc2)*1.1,is.text=FALSE)
#    scale <- max(ssb)/max(tc2) * 0.8
    for(i in 1:ncol(tc2)) menplot(ssb[non.na], cbind(0,tc2)[non.na,i:(i+1)],col=col.tmp2[i],border=NA)

    ## abline(v=summary$"SSB_MSY"/biomass.scale,lty=2)
    ## abline(v=summary$"Blimit"/biomass.scale,lty=2)
    ##     abline(v=summary$"SSB_HS"/biomass.scale,lty=2)
    ##     text(x[1],y[1],colnames(res0$vpares$ssb)[1],adj=0)
    ##     text(rev(x)[1],rev(y)[1],rev(colnames(res0$vpares$ssb))[1],adj=0)
#        points(x <- apply(res00$fout[[1]]$vssb,1,mean)[1:10]/biomass.scale,
#               y <- apply(res00$fout[[1]]$vwcaa,1,mean)[1:10]/biomass.scale,col=2,
#               type="o",pch=20,lwd=3)
#        text(rev(x)[1],rev(y)[1],
#             paste("Projection ",rownames(res00$fout[[1]]$vssb)[10],"(F_MSY)",sep=""),adj=-0.1,col=2)

#        points(x <- apply(fout0$vssb,1,mean)[1:10]/biomass.scale,
#               y <- apply(fout0$vwcaa,1,mean)[1:10]/biomass.scale,col="blue",type="o",pch=20,lwd=3)
#        text(rev(x)[1],rev(y)[1],paste("現在のFでの10年将来予測"),adj=-0.1,col="blue")

        # 魚のプロット
        waa.tmp <- (apply(res0$vpares$input$dat$waa,1,mean))^(1/3)*10
        waa.tmp <- waa.tmp/max(waa.tmp) * 0.9
        x <- tc2[which.min(abs(ssb-ssb.max*0.88)),]

        if(!is.null(image)){
            plotfish(image,x=rep(ssb.max*0.88,ncol(tc2)),y=x-diff(c(0,x))/2,
                     size=waa.tmp*0.8,scale=scale,ysize=1)
        }
        text(rep(ssb.max*0.9,ncol(tc2)),x-diff(c(0,x))/2,
             paste(0:(ncol(tc2)-1),"y/o: ",round(apply(res0$vpares$input$dat$waa,1,mean),0)," g"),cex=1)

    title("Total catch",line=-1,adj=0.1)

    ## 努力量やCVのプロット
    tmp <- round(ssb*biomass.scale)>0 & !is.na(ssb)
    matplot(ssb,tres0$fmulti,type="l",ylab="Efforts (Current=1)",col=1,lwd=2,
            xaxs="i",yaxs="i",xlab=paste("SSB (",biomass.scale,"MT)",sep=""),xlim=c(0,ssb.max),
            ylim=c(0,max(tres0$fmulti[tmp]*1.2)))
                                        #    menplot(range1,cbind(c(-100,-100),rep(max(tb2)*1.5,2)),col=gray(0.9),border=NA)
                                        #    menplot(range2,cbind(c(-100,-100),rep(max(tb2)*1.5,2)),col=gray(0.7),border=NA)
#        menplot(ssb[tmp],cbind(0,tres0$fmulti[tmp]),col=rgb(221/255,159/255,33/255,0.5),border=NA)
    plot.RP(Bref,biomass.scale=biomass.scale,ymax=max(tc2)*1.1,is.text=FALSE)
    title("Efforts",line=-1.5,font=2,adj=0.1)
    par(new=T)
#        y <- res00$trace[[1]]$ssb.CV
        y <- res00$trace[[1]]$catch.CV
        plot(ssb,y,type="l",lwd=2,col=2,axes=F,xlab="",ylab="",
             ylim=c(0,ifelse(max(y,na.rm=T)>1.5,1.5,max(y,na.rm=T))))
#        points(ssb,y,type="l",lwd=2,col=3)
        axis(side=4)
        mtext(side=4,"Catch CV",line=3,col=2,cex=0.8)

}

###############################
#### 資源量の上積みグラフを書く
###############################

plotBfish <- function(tres0,vpares, # SR.estの結果
                      b.target,ssb.max=Inf,
                      biomass.scale=1000){

    ssb <- tres0$ssb.mean/biomass.scale

    tmp <- substr(colnames(tres0),1,5)=="TB-MA"
    tb <- tres0[,tmp]/biomass.scale
    tb2 <- sapply(1:ncol(tb),function(x) apply(tb[,1:x,drop=F],1,sum,na.rm=T))

    tmp <- substr(colnames(tres0),1,5)=="TC-MA"
    tc <- tres0[,tmp]/biomass.scale
    tc2 <- sapply(1:ncol(tc),function(x) apply(tc[,1:x,drop=F],1,sum,na.rm=T))
#    {if(file.exists(pngfile)) image <- readPNG(pngfile)
#    else image <- NULL}

    year.tmp <- rev(colnames(vpares$ssb))[1:5]
    range1 <- c(0,min(ssb.max,max(ssb)))
    range2 <- range(as.data.frame(vpares$ssb)[as.character(year.tmp)])

    col.tmp1 <- rgb(40/255,96/255,163/255,seq(from=0.1,to=0.9,length=ncol(tc)))
    col.tmp2 <- rgb(100/255,200/255,44/255,seq(from=0.1,to=0.9,length=ncol(tc)))

    ### plot of SSB
#    tb3 <- tb2[which(ssb<ssb.max),]
    matplot(ssb,tb2,type="n",ylab=paste("Total biomass (",biomass.scale," MT)",sep=""),xaxs="i",yaxs="i",
            xlab="SSB (1000MT)", xlim=range1,ylim=c(0,max(tb2)))
    # 過去の時系列
        matpoints(ssb,tb2[,1],type="l",lwd=2,col="gray",lty=3)
#        points(x <- colSums(vpares$ssb)/biomass.scale,
#               y <- colSums(vpares$baa)/biomass.scale,type="o",
#               col=gray(c(seq(from=0.7,to=0,length=length(x)))),pch=20,cex=1.2,
#               lwd=3)
#        text(x[1],y[1],colnames(x)[1],adj=0)
#text(rev(x[1]),rev(y)[1],rev(colnames(x))[1],adj=0)

    text(rep(ssb[10],ncol(tb2)),
         tb2[10,],
         paste(0:(ncol(tb2)-1),"y/o"))#,round(apply(vpares$input$dat$waa,1,mean),0)," g"),cex=1)

    ssb.hist <- range(colSums(vpares$ssb)/biomass.scale)
    polygon(c(ssb.hist,rev(ssb.hist)),c(0,0,max(tb2)*10,max(tb2)*10),col=gray(0.8),border=NA)
    abline(v=b.target,col="gray")


    ## 積み上げグラフ
    non.na <- !is.na(ssb)
    for(i in 1:ncol(tb2)) menplot(ssb[non.na], cbind(0,tb2)[non.na,i:(i+1)],col=col.tmp1[i],border=NA)
                                        #    title("Total biomass",line=-1,adj=0.1)

    ##  catch
    {if(!is.null(vpares$wcaa)) wcatch <- as.numeric(colSums(vpares$wcaa))
        else{
            wcatch <- as.numeric(colSums(vpares$input$dat$caa * vpares$input$dat$waa,na.rm=T))
        }}
        matplot(ssb,tc2,type="n",,xaxs="i",yaxs="i",ylab=paste("Catch (",biomass.scale," MT)",sep=""),
                xlab="SSB (1000MT)",
                ylim=c(0,max(tc2)*1.2),xlim=range1)

    polygon(c(ssb.hist,rev(ssb.hist)),c(0,0,max(tc2)*10,max(tc2)*10),col=gray(0.8),border=NA)
    abline(v=b.target,col="gray")

    for(i in 1:ncol(tc2)) menplot(ssb[non.na], cbind(0,tc2)[non.na,i:(i+1)],col=col.tmp2[i],border=NA)

    text(rep(ssb[10],ncol(tc2)),
         tc2[10,],
         paste(0:(ncol(tc2)-1),"y/o"),col="darkgreen")#,round(apply(vpares$input$dat$waa,1,mean),0)," g"),cex=1)
                                        #    title("Total catch",line=-1,adj=0.1)

    if(0){
        plot(ssb,tres0$fmulti,type="n",lwd=2,xlim=range1,ylab="Fishing efforts",xlab="SSB (1000MT)")
        polygon(c(ssb.hist,rev(ssb.hist)),c(0,0,10,10),col=gray(0.8),border=NA)
        abline(v=b.target,col="gray")
        points(ssb,tres0$fmulti,type="l",lwd=2,xlim=range1,ylab="Fishing efforts")
        par(new=T)
        matplot(ssb,cbind(tres0$ssb.CV,tres0$catch.CV),xlim=range1,ylim=c(0,2),type="l",lty="22",col=3:4,lwd=2,axes=F,ylab="",xlab="")
        axis(side=4)
        mtext(side=4,"CV",cex=1,line=2.3)
    }
}


plot.RP <- function(rdata,RP=NULL,biomass.scale=1,ymax=1,is.text=TRUE){
    n <- length(rdata)
    rdata <- sort(rdata)
    if(is.null(RP)) RP <- names(rdata)
    ymax <- ymax * seq(from=0.5,to=1,length=n)
    for(j in 1:n){
        abline(v=rdata[j]/biomass.scale,lty=1,lwd=2,col=rgb(40/255,96/255,40/255,0.5))
        if(isTRUE(is.text)){
            text(rdata[j]/biomass.scale,ymax[j],
             paste(RP[j],"=\n",format(round(rdata[j]/biomass.scale),big.mark=","),"",sep=""),adj=0)
        }
    }
}
