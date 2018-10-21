##
get.stat <- function(fout,eyear=0,hsp=NULL,tmp.year=NULL){
    col.target <- ifelse(fout$input$N==0,1,-1)
    tmp <- as.numeric(fout$vssb[(nrow(fout$vssb)-eyear):nrow(fout$vssb),col.target])
    lhs <- sum(tmp<hsp)/length(tmp)
    if(is.null(tmp.year)) tmp.year <- (nrow(fout$vwcaa)-eyear):nrow(fout$vwcaa)

    a <- data.frame("catch.mean"=mean(fout$vwcaa[tmp.year,col.target]),
                    "catch.sd"=sd(fout$vwcaa[tmp.year,col.target]),
                    "catch.geomean"=geomean(fout$vwcaa[tmp.year,col.target]),
                    "catch.median"=median(fout$vwcaa[tmp.year,col.target],na.rm=T),
                    "catch.det"=mean(fout$vwcaa[tmp.year,1],na.rm=T),
                    "catch.L10"=quantile(fout$vwcaa[tmp.year,col.target],na.rm=T,probs=0.1),
                    "catch.H10"=quantile(fout$vwcaa[tmp.year,col.target],na.rm=T,probs=0.9),
                    "ssb.mean"=mean(fout$vssb[tmp.year,col.target]),
                    "ssb.sd"=sd(fout$vssb[tmp.year,col.target]),
                        "ssb.geomean"=geomean(fout$vssb[tmp.year,col.target]),
                        "ssb.median"=median(fout$vssb[tmp.year,col.target],na.rm=T),
                        "ssb.det"=mean(fout$vssb[tmp.year,1],na.rm=T),
                        "ssb.L10"=quantile(fout$vssb[tmp.year,col.target],na.rm=T,probs=0.1),
                        "ssb.H10"=quantile(fout$vssb[tmp.year,col.target],na.rm=T,probs=0.9),

                        "biom.mean"=mean(fout$vbiom[tmp.year,col.target]),
                        "biom.sd"=sd(fout$vbiom[tmp.year,col.target]),
                        "biom.geomean"=geomean(fout$vbiom[tmp.year,col.target]),
                        "biom.median"=median(fout$vbiom[tmp.year,col.target],na.rm=T),
                        "biom.det"=mean(fout$vbiom[tmp.year,1],na.rm=T),
                        "biom.L10"=quantile(fout$vbiom[tmp.year,col.target],na.rm=T,probs=0.1),
                        "biom.H10"=quantile(fout$vbiom[tmp.year,col.target],na.rm=T,probs=0.9),
                        "lower.HSpoint"=lhs,
                        "Fref2Fcurrent"=fout$multi
                        )
        a$U.mean <- a$catch.mean/a$biom.mean
        a$U.median <- a$catch.median/a$biom.median
        a$U.geomean <- a$catch.geomean/a$biom.geomean
        a$U.det <- a$catch.det/a$biom.det

        a$catch.CV <- a$catch.sd/a$catch.mean
        a$ssb.CV <- a$ssb.sd/a$ssb.mean
        a$biom.CV <- a$biom.sd/a$biom.mean

        Faa <- as.data.frame(t(fout$multi * fout$input$res0$Fc.at.age))
        colnames(Faa) <- paste("F",dimnames(fout$naa)[[1]],sep="")
        a <- cbind(a,Faa)
        return(a)
    }

get.stat2 <- function(fout,unit.waa=1,eyear=2,hsp=NULL,tmp.year=NULL){
    col.target <- ifelse(fout$input$N==0,1,-1)
    if(is.null(tmp.year)) tmp.year <- (nrow(fout$vwcaa)-eyear):nrow(fout$vwcaa)
        nage <- dim(fout$naa)[[1]]
        tb <- fout$naa * fout$waa * unit.waa
        if(is.null(fout$waa.catch)) fout$waa.catch <- fout$waa
        tc <- fout$caa * fout$waa.catch * unit.waa
        ssb <- fout$naa * fout$waa *fout$maa  * unit.waa
        tb.mat <- tc.mat <- ssb.mat <- matrix(0,nage,6)
        for(i in 1:nage){
            tb.mat[i,1] <- mean(tb[i,tmp.year,col.target])
            tb.mat[i,2] <- median(tb[i,tmp.year,col.target])
            tb.mat[i,3] <- geomean(tb[i,tmp.year,col.target])
            tb.mat[i,4] <- mean(tb[i,tmp.year,1])
            tb.mat[i,5:6] <- quantile(tb[i,tmp.year,col.target],probs=c(0.1,0.9),na.rm=T)

            tc.mat[i,1] <- mean(tc[i,tmp.year,col.target])
            tc.mat[i,2] <- median(tc[i,tmp.year,col.target])
            tc.mat[i,3] <- geomean(tc[i,tmp.year,col.target])
            tc.mat[i,4] <- mean(tc[i,tmp.year,1])
            tc.mat[i,5:6] <- quantile(tc[i,tmp.year,col.target],probs=c(0.1,0.9),na.rm=T)

            ssb.mat[i,1] <- mean(ssb[i,tmp.year,col.target])
            ssb.mat[i,2] <- median(ssb[i,tmp.year,col.target])
            ssb.mat[i,3] <- geomean(ssb[i,tmp.year,col.target])
            ssb.mat[i,4] <- mean(ssb[i,tmp.year,1])
            ssb.mat[i,5:6] <- quantile(ssb[i,tmp.year,col.target],probs=c(0.1,0.9),na.rm=T)
        }
        tc.mat <- as.numeric(tc.mat)
        tb.mat <- as.numeric(tb.mat)
        ssb.mat <- as.numeric(ssb.mat)

        # MA; mean, ME; median, GM; geometric mean
        names(tc.mat) <- c(paste("TC-MA-A",1:nage,sep=""),paste("TC-ME-A",1:nage,sep=""),
                           paste("TC-GM-A",1:nage,sep=""),paste("TC-DE-A",1:nage,sep=""),
                           paste("TC-L10-A",1:nage,sep=""),paste("TC-H10-A",1:nage,sep=""))
        names(tb.mat) <- c(paste("TB-MA-A",1:nage,sep=""),paste("TB-ME-A",1:nage,sep=""),
                           paste("TB-GM-A",1:nage,sep=""),paste("TB-DE-A",1:nage,sep=""),
                           paste("TB-L10-A",1:nage,sep=""),paste("TB-H10-A",1:nage,sep=""))
        names(ssb.mat) <- c(paste("SSB-GA-A",1:nage,sep=""),paste("SSB-ME-A",1:nage,sep=""),
                            paste("SSB-GM-A",1:nage,sep=""),paste("SSB-DE-A",1:nage,sep=""),
                            paste("SSB-L10-A",1:nage,sep=""),paste("SSB-H10-A",1:nage,sep=""))

        return(as.data.frame(t(c(tb.mat,tc.mat,ssb.mat))))
    }





geomean <- function(x)
{
  ifelse(all(x > 0), exp(mean(log(x))), NA)
}

