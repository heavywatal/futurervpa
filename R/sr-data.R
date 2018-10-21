
get.SRdata <- function(vpares,R.dat=NULL,SSB.dat=NULL,years=as.numeric(colnames(vpares$naa))){
    # R.datとSSB.datだけが与えられた場合、それを使ってシンプルにフィットする
    if(!is.null(R.dat) & !is.null(SSB.dat)){
        dat <- data.frame(R=R.dat,SSB=SSB.dat,year=1:length(R.dat))
    }
    else{
    # データの整形
        n <- ncol(vpares$naa)
        L <- as.numeric(rownames(vpares$naa)[1])

        dat <- list()
        dat$R <- as.numeric(vpares$naa[1,])
        dat$SSB <- as.numeric(colSums(vpares$ssb,na.rm=TRUE))
        dat$year <- as.numeric(colnames(vpares$ssb))
    # 加入年齢分だけずらす
        dat$R <- dat$R[(L+1):n]
        dat$SSB <- dat$SSB[1:(n-L)]
        dat$year <- dat$year[(L+1):n]

                                        # データの抽出
        dat <- as.data.frame(dat)
        dat <- dat[dat$year%in%years,]
    }

    class(dat) <- "SRdata"
    return(dat[c("year","SSB","R")])
}

plot.SRdata <- function(SRdata){
    plot(SRdata$SSB,SRdata$R,xlab="SSB",ylab="R",xlim=c(0,max(SRdata$SSB)),ylim=c(0,max(SRdata$R)))
}
