###

read.vpa <- function(tfile,
                     caa.label="catch at age",
                     maa.label="maturity at age",
                     waa.label="weight at age",
                     waa.biomass.label="weight at age for biomass calculation",
                     waa.catch.label="weight at age for catch calculation",
                     M.label="M at age",
                     faa.label="fishing mortality at age",
                     Fc.label="Current F",
                     naa.label="numbers at age",
                     Blimit=NULL,Pope=TRUE,fc.year=NULL){

    tmpdata <- read.csv(tfile,header=F,as.is=F,colClasses="character")

    tmpfunc <- function(tmpdata,label,type=NULL){
        flags <- which(substr(tmpdata[,1],1,1)=="#")
        flag.name <- tmpdata[flags,1]
        flag.name <- gsub("#","",flag.name)
        flag.name <- gsub(" ","",flag.name)
        get.flag <- which(flag.name==gsub(" ","",label))
        {if(length(get.flag)>0){
            tdat <- tmpdata[(flags[get.flag]+1):(flags[min(which(flags[get.flag]<flags))]-1),]
            if(!is.null(type)){
                tdat <- tdat[,!apply(tdat=="",2,all)]
                tdat <- as.numeric(tdat)
            }
            else{
                tmp.names <- list(tdat[-1,1],tdat[1,-1])
                tdat <- tdat[,!apply(tdat=="",2,all)]
                tdat <- tdat[,!apply(tdat=="",1,all)]
                tdat <- sapply((tdat[-1,-1]),as.numeric)
                dimnames(tdat) <- tmp.names
                tdat <- as.data.frame(tdat)
            }
        }
        else{
                tdat <- NULL
            }}
        tdat
    }

  dres <- list()
  dres$naa <- tmpfunc(tmpdata,naa.label)
  dres$faa <- tmpfunc(tmpdata,faa.label)

  dres$Fc.at.age <- tmpfunc(tmpdata,Fc.label,type="Fc")

  dres$input <- list()
  dres$input$dat <- list()
  dres$input$dat$maa <- tmpfunc(tmpdata,maa.label)
  dres$input$dat$caa <- tmpfunc(tmpdata,caa.label)
  dres$input$dat$M <- tmpfunc(tmpdata,M.label)
  dres$input$dat$waa <- tmpfunc(tmpdata,waa.label)
  dres$input$dat$waa <- tmpfunc(tmpdata,waa.biomass.label)
  dres$input$dat$waa.catch <- tmpfunc(tmpdata,waa.catch.label)

  dres$ssb <- dres$input$dat$waa * dres$input$dat$maa * dres$naa
  dres$ssb <- as.data.frame(dres$ssb)

  dres$baa <- dres$input$dat$waa * dres$naa
  dres$baa <- as.data.frame(dres$baa)

  # setting total catch
  dres$waa.catch <- dres$input$dat$waa * dres$input$dat$caa
  dres$waa.catch <- as.data.frame(dres$waa.catch)

  dres$Blimit <- Blimit

  if(is.null(Pope)){
      caa.pope  <- dres$naa*(1-exp(-dres$faa))*exp(-dres$input$dat$M/2)
        diff.pope <- mean(unlist(dres$input$dat$caa/caa.pope))

        faa <- dres$faa
        M <- dres$input$dat$M
        caa.bara <- dres$naa*faa/(faa+M)*(1-exp(-faa-M))
        diff.bara <- mean(unlist(dres$input$dat$caa/caa.bara))

        if(abs(1-mean(diff.bara))>abs(1-mean(diff.pope))){
            dres$Pope <- TRUE
        }
        else{
            dres$Pope <- FALSE
        }
  }
    if(is.null(dres$Fc.at.age) && !is.null(fc.year)) dres$Fc.at.age <- apply(dres$faa[,colnames(dres$faa)%in%fc.year],1,mean)
  dres
}
