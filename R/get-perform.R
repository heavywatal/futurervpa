#' Get perform
#'
#' @rdname get-perform
#' @export
get.perform <- function(fout0,Blimit=0,longyear=50,smallcatch=0.5,N=NULL,
                        shortyear=c(3,5,10),tmp.year=NULL){
    stat1 <- get.stat(fout0,eyear=0,hsp=Blimit,tmp.year=tmp.year)[c("catch.mean","catch.CV","biom.mean","biom.CV","ssb.mean","lower.HSpoint")]
    stat2 <- get.stat2(fout0,eyear=0,tmp.year=tmp.year)
    stat2 <- data.frame(t(as.data.frame(strsplit(colnames(stat2),"-"))),value=as.numeric(stat2))
    rownames(stat2) <- NULL

    # waaによる加重平均年齢&組成
    xx <- subset(stat2,X1=="TB" & X2=="MA")
    nage <- sum(!is.na(xx$value))
    tmp <- c(rep(2,ceiling(nage/3)),rep(3,ceiling(nage/3)))
    tmp <- c(rep(1,nage-length(tmp)),tmp)
    if(sum(tmp==1)==0 & sum(tmp==2)>1) tmp[1] <- 1

    xx$bvalue <- xx$value * fout0$waa[,1,1]
    xx$waa <- fout0$waa[,1,1]
    large.portion1 <- tapply(xx$bvalue[!is.na(xx$bvalue)],tmp,sum,na.rm=T)
    stat1$largefish.nature <- large.portion1[names(large.portion1)==3]/sum(large.portion1)
    aage.biom <- sum(xx$bvalue * 0:(length(xx$bvalue)-1))/sum(xx$bvalue)

    xx <- subset(stat2,X1=="TC" & X2=="MA")
    xx$bvalue <- xx$value * fout0$waa[,1,1]
    aage.catch <- sum(xx$bvalue * 0:(length(xx$bvalue)-1))/sum(xx$bvalue)
    large.portion2 <- tapply(xx$bvalue[!is.na(xx$bvalue)],tmp,sum,na.rm=T)
    stat1$largefish.catch <- large.portion2[names(large.portion2)==3]/sum(large.portion2)

    # 漁獲量<0.5平均漁獲量の頻度
    if(is.null(tmp.year)) tmp.year <- nrow(fout0$vwcaa)
    stat1$catch.safe <- 1/mean(fout0$vwcaa[tmp.year,]<smallcatch*mean(fout0$vwcaa[tmp.year,]))
    stat1$catch.safe <- ifelse(stat1$catch.safe>longyear,longyear,stat1$catch.safe)

    # 親魚量<Blimitの頻度　→　確率の逆数
    stat1$ssb.safe <- 1/stat1$"lower.HSpoint"
    stat1$ssb.safe <- ifelse(stat1$ssb.safe>longyear,longyear,stat1$ssb.safe)

    # ABC.yearから5年目までの平均累積漁獲量
    short.catch <- numeric()
    for(i in 1:length(shortyear)){
        years <- fout0$input$ABC.year:(fout0$input$ABC.year+shortyear[i])
        short.catch[i] <- mean(apply(fout0$vwcaa[rownames(fout0$vwcaa)%in%years,-1],2,sum))
    }
    names(short.catch) <- paste("short.catch",shortyear,sep="")
    short.catch <- as.data.frame(t(short.catch))

    # 平衡状態になった年
    years <- names(fout0$vssb[,1])[-1]
    heikou.diff <- which(diff(fout0$vssb[,1])/fout0$vssb[-1,1]<0.01)
    if(length(heikou.diff)>0) stat1$eq.year <- years[min(heikou.diff)] else stat1$eq.year <- Inf

    dat <- data.frame(stat1,short.catch,aage.biom=aage.biom,aage.catch=aage.catch,effort=fout0$multi,
                      waa=as.data.frame(t(fout0$waa[,1,1])),meigara=as.data.frame(t(tmp)))
    return(dat)
}
