plotyield <- function(res00,int.res=NULL,detail.plot=FALSE){
#    par(mfrow=c(2,1))
    arg.tmp <- res00$farg
    arg.tmp$rec.arg$sd <- 0
    arg.tmp$N <- 1
#    fout.tmp <- do.call(future.vpa2,arg.tmp)

    # average
    plot(x <- res00$trace[[1]]$fmulti,y <- res00$trace[[1]]$catch.mean,type="n",xlim=c(0,max(x)),
         xlab="Multiplier to current F",ylab="Catch weight",ylim=c(0,max(res00$trace[[1]]$catch.det,y)))
    menplot(res00$trace[[1]]$fmulti,cbind(res00$trace[[1]]$catch.L10,res00$trace[[1]]$catch.H10),
            col=rgb(210/255,94/255,44/255,0.3),border=NA)

    ## integrate
    if(!is.null(int.res)){
        points(x,y <- int.res$yield,lty=2,type="o",lwd=1,col="gray")
        points(fmax5 <- x[which.max(y)],y[which.max(y)],pch=20,col="gray")
    }

    points(x <- res00$trace[[1]]$fmulti,y <- res00$trace[[1]]$catch.mean,type="l",xlim=c(0,max(x)),
           xlab="Multiplier to current F",ylab="Catch weight",ylim=c(0,max(res00$trace[[1]]$catch.det,y)))
    points(fmax1 <- x[which.max(y)],y[which.max(y)],pch=20,col=1)

    if(isTRUE(detail.plot)){
    # geomean
        points(x <- res00$trace[[1]]$fmulti,y <- res00$trace[[1]]$catch.geomean,col=2,type="l",xlim=c(0,2))
    points(fmax2 <- x[which.max(y)],y[which.max(y)],pch=20,col=2)

    # median
        points(x <- res00$trace[[1]]$fmulti,y <- res00$trace[[1]]$catch.median,col=3,type="l",xlim=c(0,2))
        points(fmax3 <- x[which.max(y)],y[which.max(y)],pch=20,col=3)
    }

    # deteministic
    points(x <- res00$trace[[1]]$fmulti,y <- res00$trace[[1]]$catch.det,col=4,
           type="l",xlim=c(0,2))
    points(fmax4 <- x[which.max(y)],y[which.max(y)],pch=20,col=4)

    title("Yield vs. F")

    ## plot CV of yield
    par(new=T)
    y <- res00$trace[[1]]$catch.CV
    plot(x,y,type="l",lwd=3,
         col=rgb(0.8,0.8,0,0.6),axes=F,xlab="",ylab="",
         ylim=c(0,ifelse(max(y,na.rm=T)>1.5,1.5,max(y,na.rm=T))))
    axis(side=4)
    mtext(side=4,"CV of catch",line=2.5,col=rgb(0.8,0.8,0,0.6),cex=0.8)

    ### plot SSB
    plot(x <- res00$trace[[1]]$fmulti,y <- res00$trace[[1]]$ssb.mean,type="n",xlim=c(0,max(x)),
         xlab="Relative F (to current F)",ylab="SSB")
    menplot(res00$trace[[1]]$fmulti,cbind(res00$trace[[1]]$ssb.L10,res00$trace[[1]]$ssb.H10),
            col=rgb(40/255,96/255,163/255,0.3),border=NA)

    ## integrate
    if(!is.null(int.res)){
        points(x,y <- int.res$ssb,lty=2,type="o",lwd=1,col="gray")
        points(fmax5,y[x==fmax5],pch=20,col="gray")
    }

    points(x <- res00$trace[[1]]$fmulti,y <- res00$trace[[1]]$ssb.mean,type="l",xlim=c(0,max(x)),
         xlab="Relative F (to current F)",ylab="SSB")
    points(fmax1,y[x==fmax1],pch=20,col=1)

    if(isTRUE(detail.plot)){
        points(x <- res00$trace[[1]]$fmulti,y <- res00$trace[[1]]$ssb.geomean,col=2,type="l",xlim=c(0,2))
        points(fmax2,y[x==fmax2],pch=20,col=2)

        points(x <- res00$trace[[1]]$fmulti,y <- res00$trace[[1]]$ssb.median,col=3,type="l",xlim=c(0,2))
        points(fmax3,y[x==fmax3],pch=20,col=3)
    }

    points(x <- res00$trace[[1]]$fmulti,y <- res00$trace[[1]]$ssb.det,
           col=4,type="l")
    points(fmax4,y[x==fmax4],pch=20,col=4)
    title("SSB vs. F")
    if(!is.null(int.res)){
        legend("topright",lty=c(1,1,1,1,2,NA),col=c(1:4,"gray",NA),legend=c("Simple mean","Geometric mean","Median","Deterministic","Integrate","fill: 80% conf"),bty="n")
    }
    else{
        if(isTRUE(detail.plot)){
            legend("topright",lty=c(1,1,1,1,NA),col=c(1:4,NA),
                   legend=c("Simple mean","Geometric mean","Median","Deterministic","fill: 80% conf"))
        }
        else{
            legend("topright",lty=c(1,1,NA),col=c(c(1,4),NA),
                   legend=c("Simple mean","Deterministic","fill: 80% conf"))
        }
    }

    #### CV plot
    par(new=T)
    y <- res00$trace[[1]]$ssb.CV
    plot(x,y,type="l",lwd=3,
         col=rgb(0.8,0.8,0,0.6),axes=F,xlab="",ylab="",
         ylim=c(0,ifelse(max(y,na.rm=T)>1.5,1.5,max(y,na.rm=T))))
    axis(side=4)
    mtext(side=4,"CV of SSB",line=2.5,col=rgb(0.8,0.8,0,0.6),cex=0.8)

#    points(fout.tmp$multi,fout.tmp$vssb[100,1],pch=4)
}
