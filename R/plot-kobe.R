# MSY計算で仮定されている選択率で漁獲したとき何倍になるか？ => 計算時間が、、。
# %SPR？MSYを達成した時のFを%SPR換算
plot.Kobe0 <- function(srres,pickB="",what.plot="hs",plot.history=FALSE){
    tmp <- which(names(srres)==what.plot)-3
    years <- colnames(srres$vpares$ssb)
    y <- srres[what.plot][[1]]$Fhist[[1]]$fmulti/srres[what.plot][[1]]$f.msy
    x <- as.numeric(colSums(srres$vpares$ssb))/
                      srres$summary[pickB][[1]][tmp]
    x <- x[y>0.001]
    years <- years[y>0.001]
    y <- y[y>0.001]
    plot(x,y,type="n",xlim=c(0,ifelse(max(x)<3,3,max(x,na.rm=T))),
#         ylim=c(0,ifelse(max(y)<3,3,max(y))),
         ylim=c(0,4),
         pch=c(3,rep(1,length(y)-2),20),col=c(1,rep(1,length(y)-2),2),
         cex=c(1,rep(1,length(y)-2),2),ylab="F/Fmsy",xlab="B/Bmsy")
    polygon(c(-1,1,1,-1),c(-1,-1,1,1),col="khaki1",border=NA)
    polygon(c(1,4,4,1),c(-1,-1,1,1),col="olivedrab2",border=NA)
    polygon(c(1,4,4,1),c(1,1,6,6),col="khaki1",border=NA)
    polygon(c(-1,1,1,-1),c(1,1,6,6),col="indianred1",border=NA)
    axis(side=1:2)

    points(x,y,type="o",
           pch=c(3,rep(1,length(y)-1),20),
           col=c(1,rep(1,length(y)-1),1),
           cex=c(1,rep(1,length(y)-1),2),ylab="F/Fmsy",xlab="B/Bmsy")
    points(rev(x)[1],rev(y)[1],pch=20)

    if(isTRUE(plot.history)){
      plot(years,y,type="b",ylab="F/Fmsy",ylim=c(0,max(y)))
      abline(h=1)
      plot(years,x,type="b",xlab="B/Bmsy",ylim=c(0,max(y)))
      abline(h=1)
    }

    invisible(data.frame(years=years,F=y,B=x))
}

plot.Kobe2 <- get.trend <- function(srres,UBdata=NULL,SR="hs",plot.history=FALSE,is.plot=FALSE,pickU="",pickB="",ylab.tmp="U/Umsy",xlab.tmp="SSB/SSBmsy"){

    dres <- srres$vpares
    tmp <- which(names(srres)==SR)-3

    if(is.null(dres$TC.MT)) dres$TC.MT <- as.numeric(colSums(dres$wcaa))

    if(is.null(UBdata)){
    U <- data.frame(years=as.numeric(colnames(dres$baa)),
                    U=as.numeric(dres$TC.MT)/as.numeric(colSums(dres$baa,na.rm=T)))
    B <- data.frame(years=as.numeric(colnames(dres$ssb)),
                    B=as.numeric(colSums(dres$ssb)))
    Catch <- data.frame(years=as.numeric(colnames(dres$baa)),
                    C=as.numeric(dres$TC.MT))
    UBdata <- merge(U,B)
    UBdata <- merge(UBdata,Catch)

#    U <- data.frame(years=as.numeric(ts$YEAR),
#                    U=as.numeric(ts$"TC-MT")/as.numeric(ts$"TB-MT"))

#    UBdata$Umsy <- srres$summary$MSY.U.median2[tmp]
#    UBdata$Bmsy <- srres$summary$MSY.ssb.median2[tmp]

    UBdata$Umsy <- srres$summary[pickU][tmp,]
    UBdata$Bmsy <- srres$summary[pickB][tmp,]

    UBdata$Uratio <- UBdata$U/UBdata$Umsy
    UBdata$Bratio <- UBdata$B/UBdata$Bmsy
    }

    if(is.plot){
      plot(x <- UBdata$Bratio,
           y <- UBdata$Uratio,type="n",xlim=c(0,ifelse(max(x)<2,2,max(x,na.rm=T))),
           ylim=c(0,ifelse(max(y,na.rm=T)<3,3,max(y,na.rm=T))),
           cex=c(1,rep(1,length(y)-2),3),ylab=ylab.tmp,xlab=xlab.tmp)
      polygon(c(-1,1,1,-1),c(-1,-1,1,1),col="khaki1",border=NA)
      polygon(c(1,6,6,1),c(-1,-1,1,1),col="olivedrab2",border=NA)
      polygon(c(1,6,6,1),c(1,1,6,6),col="khaki1",border=NA)
      polygon(c(-1,0.5,0.5,-1),c(1,1,6,6),col="indianred1",border=NA)
      polygon(c(0.5,1,1,0.5),c(1,1,6,6),col="tan1",border=NA)
      polygon(c(-1,0.5,0.5,-1),c(-1,-1,1,1),col="khaki2",border=NA)
      polygon(c(0.5,1,1,0.5),c(-1,-1,1,1),col="khaki1",border=NA)
      axis(side=1:2)

#      points(x,y,type="o",pch=c(3,rep(1,length(y)-2),20),col=c(1,rep(1,length(y)-2),1),cex=c(1,rep(1,length(y)-2),1.5))

      points(x,y,type="l",pch=20,col=1,lwd=1)
      points(x,y,type="p",pch=20,col=gray(c(seq(from=0.7,to=0,length=length(x)))),cex=1.2)
      points(rev(x)[1],rev(y)[1],type="p",pch=20,cex=2.5)

    if(isTRUE(plot.history)){
      plot(UBdata$years,y,type="b",ylab="F/Fmsy",ylim=c(0,max(y)))
      abline(h=1)
      plot(UBdata$years,x,type="b",xlab="B/Bmsy",ylim=c(0,max(y)))
      abline(h=1)
    }}

    invisible(UBdata)
}


plot.kobemat <- function(xx,title.name="",line=0){
    yy  <- as.data.frame.table(xx)
    yy$pch <- 20
    yy$color <- "gray"
    yy$color[as.numeric(yy[,3])>0.5] <- "red"
    yy$color[yy[,3]<0.5 & yy[,3]>0.4] <- "pink"
    plot(as.numeric(as.character(yy[,1])),
         as.numeric(as.character(yy[,2])),type="n",col=yy$color,xlab="Years",pch=yy$pch,cex=3,
         ylab="multiplier to F_current")
    abline(h=seq(from=0,to=10,by=0.1),v=2000:2100,col="gray")
    points(as.numeric(as.character(yy[,1])),
         as.numeric(as.character(yy[,2])),col=yy$color,pch=yy$pch,cex=3)
    title(title.name)
    abline(h=line,col="red")
    text(2+min(as.numeric(as.character(yy[,1]))),line,paste("F_",title.name,sep=""))
    legend("topleft",col=c("gray","pink","red"),legend=c("<40%","40-50%",">50%"),
           title="Pr(B>Btarget)",pch=20,pt.cex=3,bg="white")
}

plot.kobemat2 <- function(yy,...){
	xx <- yy$prob.btarget
    for(i in 1:length(xx)){
		plot.kobemat(xx[[i]],title.name=names(xx)[i],line=-1)
	}
	matplot(yy$ssb,type="l",ylim=c(0,max(yy$ssb)),lty=1)
	title("SSB",line=-1)
	matplot(yy$biom,type="l",ylim=c(0,max(yy$biom)),lty=1)
	title("Biomass",line=-1)
	matplot(yy$catch,type="l",ylim=c(0,max(yy$catch)),lty=1)
	title("Catch",line=-1)
	legend("bottomright",col=1:ncol(yy$ssb),legend=colnames(yy$ssb),lty=1,title="Fcurrentx")
}


plot.kobe <- function(vpares,Bref,Uref,plot.history=FALSE,is.plot=FALSE,pickU="",pickB="",ylab.tmp="U/Umsy",xlab.tmp="SSB/SSBmsy",title.tmp=""){

    vpares$TC.MT <- as.numeric(colSums(vpares$wcaa))
    UBdata <- data.frame(years=as.numeric(colnames(vpares$baa)),
                         U=as.numeric(vpares$TC.MT)/as.numeric(colSums(vpares$baa,na.rm=T))/Uref,
                         B=as.numeric(colSums(vpares$ssb))/Bref)

    x <- UBdata$B
    y <- UBdata$U
    tmp <- x>0 & y>0
    x <- x[tmp]
    y <- y[tmp]
    plot(x,
         y,type="n",xlim=c(0,ifelse(max(x)<2,2,max(x,na.rm=T))),
         ylim=c(0,ifelse(max(y,na.rm=T)<3,3,max(y,na.rm=T))),
         cex=c(1,rep(1,length(y)-2),3),ylab=ylab.tmp,xlab=xlab.tmp)
    polygon(c(-1,1,1,-1),c(-1,-1,1,1),col="khaki1",border=NA)
    polygon(c(1,6,6,1),c(-1,-1,1,1),col="olivedrab2",border=NA)
    polygon(c(1,6,6,1),c(1,1,6,6),col="khaki1",border=NA)
    polygon(c(-1,0.5,0.5,-1),c(1,1,6,6),col="indianred1",border=NA)
    polygon(c(0.5,1,1,0.5),c(1,1,6,6),col="tan1",border=NA)
    polygon(c(-1,0.5,0.5,-1),c(-1,-1,1,1),col="khaki2",border=NA)
    polygon(c(0.5,1,1,0.5),c(-1,-1,1,1),col="khaki1",border=NA)
    axis(side=1:2)

#      points(x,y,type="o",pch=c(3,rep(1,length(y)-2),20),col=c(1,rep(1,length(y)-2),1),cex=c(1,r
      points(x,y,type="l",pch=20,col=1,lwd=1)
      points(x,y,type="p",pch=20,col=gray(c(seq(from=0.7,to=0,length=length(x)))),cex=1.2)
      points(rev(x)[1],rev(y)[1],type="p",pch=20,cex=2.5)
    title(title.tmp,adj=0.8,line=-2)

    if(isTRUE(plot.history)){
      plot(UBdata$years,y,type="b",ylab="F/Fmsy",ylim=c(0,max(y)))
      abline(h=1)
      plot(UBdata$years,x,type="b",xlab="B/Bmsy",ylim=c(0,max(y)))
      abline(h=1)
    }

    invisible(UBdata)
}
