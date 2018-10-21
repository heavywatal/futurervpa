
plot.SR <- function(srres,what.plot=c("hs","bh","ri","sl"),xyscale=c(1.3,1.3),xscale=FALSE,is.legend=TRUE,what.sigma=1,FUN="mean",is.MSYline=TRUE,
                    pick="SSB_MSY"){

    col.tmp <- c(rgb(0.3,0.8,0.3,alpha=0.8),rgb(0.8,0.3,0.3,alpha=0.8),rgb(0.3,0.3,0.8,alpha=0.8))

    # xscale=TRUEの場合、B0が再生産関係によって異なってくるので、複数の重ね書きはしないこSと！
#    tmp <- which(names(srres)==what.plot)
    tmp <- which(names(srres)%in%what.plot)
    resid <- list()

    if(isTRUE(xscale)){
        ssb0 <- srres$summary$B0.ssb.mean1[tmp]
        xrange <- seq(from=0,to=ssb0,length=100)
    }
    else{
        ssb0 <- 1
        xrange <- seq(from=0,to=xyscale[1]*max(srres$dat$SSB,na.rm=T),length=100)
    }

    plot(x <- srres$dat$SSB/ssb0,y <- srres$dat$R,type="l",pch=20,xlim=range(xrange/ssb0),col="gray",
         ylim=c(0,xyscale[2]*max(y,na.rm=T)),xaxs="i",yaxs="i",xlab=ifelse(!xscale,"Spawning biomass (MT)","SB/SB0"),
         ylab="Number of recruits",lwd=1)
    points(x,y,type="p",pch=20,col=gray(c(seq(from=0.7,to=0,length=length(x)))))
    points(rev(x)[1],rev(y)[1],type="p",pch=20,cex=2.5)
    for(i in 1:length(what.plot)){
        Bmsy <- srres$summary[pick][which(what.plot[i]==rownames(srres$summary)),]
        if(what.plot[i]=="hs"){
            points(xpred <- xrange/ssb0,
                   ypred <- pred.HS(SSB=xrange,
                                       a=srres[what.plot[i]][[1]]$a,b=srres[what.plot[i]][[1]]$b,gamma=srres[what.plot[i]][[1]]$gamma),
                   type="l",lwd=2,col=col.tmp[i],lty=1)
            resid[[i]] <- pred.HS(SSB=x,
                             a=srres[what.plot[i]][[1]]$a,
                             b=srres[what.plot[i]][[1]]$b,
                             gamma=srres[what.plot[i]][[1]]$gamma)
            resid[[i]] <- log(y)-log(resid[[i]])
            Rmsy <- pred.HS(SSB=Bmsy,
                            a=srres[what.plot[i]][[1]]$a,b=srres[what.plot[i]][[1]]$b,gamma=srres[what.plot[i]][[1]]$gamma)
        }
        if(what.plot[i]=="bh"|what.plot[i]=="ri"){
          if(what.plot[i]=="bh") tmpfunc <- pred.BH
          if(what.plot[i]=="ri") tmpfunc <- pred.RI
          points(xpred <- xrange/ssb0,
                 ypred <- tmpfunc(SSB=xrange,
                                     a=srres[what.plot[i]][[1]]$a,b=srres[what.plot[i]][[1]]$b),type="l",lwd=2,col=col.tmp[i],lty=ifelse(what.plot[i]=="bh",2,3))
          resid[[i]] <- tmpfunc(SSB=x,
                                a=srres[what.plot[i]][[1]]$a,
                                b=srres[what.plot[i]][[1]]$b)
          resid[[i]] <- log(y)-log(resid[[i]])
          Rmsy <- tmpfunc(SSB=Bmsy,
                          a=srres[what.plot[i]][[1]]$a,b=srres[what.plot[i]][[1]]$b)
        }
        if(what.plot[i]=="sl")
          points(xrange/ssb0,pred.SL(SSB=xrange,
                    a=srres[what.plot[i]][[1]]$a),type="l",lwd=1,col=col.tmp[i])



        if(is.MSYline){ #abline(v=Bmsy/ssb0,col=i+1,lty=2)
            arrows(Bmsy/ssb0,Rmsy*1.2,Bmsy/ssb0,Rmsy,col=col.tmp[i],lty=1,lwd=2,length=.1)
            if(Bmsy/ssb0>rev(xrange)[1]){
                ymax <- xyscale[2]*max(y,na.rm=T)
                arrows(rev(xrange)[2],ifelse(ymax<Rmsy,ymax*0.8,Rmsy*0.8),
                       rev(xrange)[2],ifelse(ymax<Rmsy,ymax,Rmsy),col=col.tmp[i],lty=1,lwd=2,length=.1)
            }
        }
    }

    neg.LL <- sapply(srres[what.plot],function(x) x$res$value)
    k <- sapply(srres[what.plot],function(x) length(x$res$par))
    n <- length(srres$dat$R)
    AICc <- 2*neg.LL+2*k+2*k*(k+1)/(n-k-1)
    if(is.legend){
        legend("topright",legend=paste(toupper(what.plot[order(AICc)]),round(AICc[order(AICc)],2)),
               col=col.tmp[order(AICc)],lwd=1,title="AICc",bg="white",ncol=3)
    }

    return(list(AICc=AICc,resid=resid,x=xpred,y=ypred))
}
