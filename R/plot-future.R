
plot.futures <- function(fres.list,conf=c(0.1,0.5,0.9),target="SSB",legend.text="",xlim.tmp=NULL,y.scale=1){
    if(target=="SSB")  aa <- lapply(fres.list,function(x) apply(x$vssb[,-1],1,quantile,probs=conf))
    if(target=="Biomass") aa <- lapply(fres.list,function(x) apply(x$vbiom[,-1],1,quantile,probs=conf))
    if(target=="Catch") aa <- lapply(fres.list,function(x) apply(x$vwcaa[,-1],1,quantile,probs=conf))
    if(target=="Recruit") aa <- lapply(fres.list,function(x) apply(x$naa[1,,-1],1,quantile,probs=conf))

    if(is.null(xlim.tmp)) xlim.tmp <- as.numeric(range(unlist(sapply(aa,function(x) colnames(x)))))
    plot(0,max(unlist(aa)),type="n",xlim=xlim.tmp,
         ylim=y.scale*c(0,max(unlist(aa))),xlab="Year",ylab=target)
    lapply(1:length(aa),function(i) matpoints(colnames(aa[[i]]),t(aa[[i]]),col=i,type="l",lty=c(2,1,2)))
    legend("bottomright",col=1:length(aa),legend=legend.text,lty=1)
    invisible(aa)
}

plot.future <- function(fres0,ylim.tmp=NULL,xlim.tmp=NULL,vpares=NULL,what=c(TRUE,TRUE,TRUE),conf=0.1,N.line=0,
                        label=c("Biomass","SSB","Catch"),is.legend=TRUE,add=FALSE,col=NULL,...){
    ## 暗黙に、vssbなどのmatrixの1列目は決定論的なランの結果と仮定している
    if(is.null(col)) col <- 1
    matplot2 <- function(x,add=FALSE,...){
        if(add==FALSE) matplot(rownames(x),x,type="l",lty=c(2,1,2),col=col,xlab="Year",...)
        if(add==TRUE) matpoints(rownames(x),x,type="l",lty=c(2,1,2),col=col,xlab="Year",...)
    }

    if(is.null(xlim.tmp)) xlim.tmp <- range(as.numeric(colnames(fres0$naa)))

    if(what[1]){
        matplot2(x <- t(apply(fres0$vbiom[,-1],1,quantile,probs=c(conf,0.5,1-conf))),
                 ylim=c(0,ifelse(is.null(ylim.tmp),max(x),ylim.tmp[1])),
                 xlim=xlim.tmp,
                 ylab=label[1],main=label[1],add=add,...)
        points(rownames(fres0$vbiom),apply(fres0$vbiom[,-1],1,mean),type="b",pch=1)
        points(rownames(fres0$vbiom),as.numeric(fres0$vbiom[,1]),type="b",pch=3)
        if(!is.null(vpares)){
            points(colnames(vpares$baa),colSums(vpares$baa),type="o",pch=20)
        }
    }

  if(what[2]){
    matplot2(x <- t(apply(fres0$vssb[,-1],1,quantile,probs=c(conf,0.5,1-conf))),
             ylim=c(0,ifelse(is.null(ylim.tmp),max(x),ylim.tmp[2])),
             xlim=xlim.tmp,
             ylab=label[2],main=label[2],add=add,...)
    points(rownames(fres0$vssb),apply(fres0$vssb[,-1],1,mean),type="b",pch=1)
    points(rownames(fres0$vssb),as.numeric(fres0$vssb[,1]),type="b",pch=3)
    if(!is.null(fres0$input$Frec))
        if(!is.null(fres0$input$Frec$scenario))
        if(fres0$input$Frec$scenario!="catch.mean"){
            abline(h=fres0$input$Frec$Blimit,col=2)
            abline(v=fres0$input$Frec$future.year,col=2)
        }
    if(!is.null(vpares)){
      points(colnames(vpares$ssb),colSums(vpares$ssb),type="o",pch=20)
    }
    if(N.line>0) matpoints(rownames(fres0$vssb),fres0$vssb[,2:(N.line+1)],col="gray",type="l",lty=1)
  }

  if(what[3]){
    matplot2(x <- t(apply(fres0$vwcaa[,-1],1,quantile,probs=c(conf,0.5,1-conf))),
             ylim=c(0,ifelse(is.null(ylim.tmp),max(x),ylim.tmp[3])),
             xlim=xlim.tmp,
             ylab=label[3],main=label[3],add=add,...)
    points(rownames(fres0$vwcaa),apply(fres0$vwcaa[,-1],1,mean),type="b",pch=1)
    points(rownames(fres0$vwcaa),as.numeric(fres0$vwcaa[,1]),type="b",pch=3)
    if(!is.null(fres0$input$Frec))
        if(fres0$input$Frec$scenario=="catch.mean"){
        abline(h=fres0$input$Frec$Blimit,col=2)
        abline(v=fres0$input$Frec$future.year,col=2)
        }
    if(!is.null(vpares)){
      points(colnames(vpares$baa),colSums(vpares$input$dat$caa*vpares$input$dat$waa),type="o",pch=20)
    }
    if(N.line>0) matpoints(rownames(fres0$vwcaa),fres0$vwcaa[,2:(N.line+1)],col="gray",type="l",lty=1)
  }
  if(is.legend){
    if(sum(what)>1) plot(1:10,type = "n",ylab = "", xlab = "", axes = F)
    legend("topleft",lty=c(NA,NA,1,2),legend=c("Deterministic","Mean","Median",paste(100-(conf*2)*100,"%conf")),pch=c(3,1,NA,NA))
  }

}

#print.future <- function(fres){ # S3 method を使いたいんですが、まだいまいちわかりません
#  cat(fres$ABC[1])
#}
#
