############## 岡村さん作成関数 ##################3

future.vpa1 <- function(
     vpares,
     multi=1,
     nyear=50, # 将来予測の年数
     ABC.year=2018, # ABCを計算する年
     waa.year=2015:2016, # 生物パラメータの参照年
     maa.year=2015:2016,
     M.year=2015:2016,
     seed=1,
     N=100,
     naa0=NULL,
     eaa0=NULL,
     beta=1,
     delta=0,
     Blim=0,
     Bban=0,
     Pope=FALSE,
     ssb0=NULL,
     faa0=NULL,
     # recfuncに対する引数
     rec.arg=list(a=res1$pars[1],b=res1$pars[2],gamma=res1$gamma,
     sd=res1$pars[3],bias.correction=TRUE,rho=res1$pars[4],resid=res1$resid)
){

#    print(multi)
    argname <- ls()
    arglist <- lapply(argname,function(x) eval(parse(text=x)))
    names(arglist) <- argname

set.seed(seed)

lag <- as.numeric(rownames(vpares$input$dat$caa))[1]

nY <- ncol(vpares$input$dat$caa)
final.year <- as.numeric(colnames(vpares$input$dat$caa)[nY])

waa <- vpares$input$dat$waa
waa <- rowMeans(waa[,colnames(waa) %in% waa.year])
maa <- vpares$input$dat$maa
maa <- rowMeans(maa[,colnames(maa) %in% maa.year])
M <- vpares$input$dat$M
A <- nrow(M)
M <- rowMeans(M[,colnames(M) %in% M.year])

waa <- array(waa,dim=c(A,nyear+1,N))
maa <- array(maa,dim=c(A,nyear+1,N))
M <- array(M,dim=c(A,nyear+1,N))

naa <- baa <- ssb <- faa <- caa <- wcaa <- array(NA,dim=c(A,nyear+1,N))

if (is.null(faa0)) Fc.at.age <- vpares$Fc.at.age else Fc.at.age <- faa0

# 1st year

  if (is.null(naa0)) {
    waa[,1,] <- vpares$input$dat$waa[,nY]
    maa[,1,] <- vpares$input$dat$maa[,nY]
    M[,1,] <- vpares$input$dat$M[,nY]
    naa[,1,] <- vpares$naa[,nY]
    faa[,1,] <- vpares$faa[,nY]
    baa[,1,] <- naa[,1,]*waa[,1,]
    ssb[,1,] <- baa[,1,]*maa[,1,]
    faa[,2,] <- Fc.at.age
    faa[,-(1:2),] <- multi*Fc.at.age
  } else {
    naa[,1,] <- naa0
    faa[,,] <- multi*Fc.at.age
    baa[,1,] <- naa[,1,]*waa[,1,]
    ssb[,1,] <- baa[,1,]*maa[,1,]
  }

eaa <- matrix(NA,nrow=nyear+1,ncol=N)

sd2 <- sqrt(rec.arg$sd^2/(1-rec.arg$rho^2))
if (is.null(eaa0)) eaa[1,] <- rec.arg$resid[length(rec.arg$resid)] else eaa[1,] <- eaa0

eaa[-1,] <- rnorm(nyear*N,0,rec.arg$sd)

if(class(ssb0)=="matrix") ssb0 <- array(ssb0,dim=c(A,1,N))

if (is.null(ssb0) & (2 - lag <= 0)) {
  ssb0 <- array(NA, dim=c(A,abs(2-lag)+1,N))
  for (j in 1:(abs(2-lag)+1)){
    ssb0[,j,] <- vpares$ssb[,nY-j]
  }
}

alpha <- array(1,dim=c(A,nyear+1,N))

# 2nd year and onward

for (i in 2:(nyear+1)){
  eaa[i,] <- rec.arg$rho*eaa[i-1,]+eaa[i,]
  naa[2:A,i,] <- naa[1:(A-1),i-1,]*exp(-M[1:(A-1),i-1,]-alpha[1:(A-1),i-1,]*faa[1:(A-1),i-1,])
  naa[A,i,] <- naa[A,i,]+naa[A,i-1,]*exp(-M[A,i-1,]-alpha[A,i-1,]*faa[A,i-1,])
  ssb[,i,] <- naa[,i,]*waa[,i,]*maa[,i,]
  if(i-lag > 0) SSB <- colSums(ssb[,i-lag,,drop=FALSE],na.rm=TRUE) else SSB <- colSums(ssb0[,i-1,,drop=FALSE])
  naa[1,i,] <- HS(SSB,rec.arg$a,rec.arg$b,rec.arg$gamma,HStype="HS")*exp(eaa[i,]-0.5*sd2^2)
  baa[,i,] <- naa[,i,]*waa[,i,]
  ssb[,i,] <- baa[,i,]*maa[,i,]
  Bcur <- colSums(ssb[,i,,drop=FALSE],na.rm=TRUE)
  alpha[,i,] <- beta*matrix(ifelse(Bcur > Blim, 1, ifelse(Bcur > Bban, ((Bcur-Bban)/(Blim-Bban))^delta, 0)),byrow=TRUE,nrow=A,ncol=N)
}

  if (Pope) caa <- (1-exp(-alpha*faa))*exp(-M/2)*naa else caa <- naa*(1-exp(-alpha*faa-M))*alpha*faa/(alpha*faa+M)

  wcaa <- caa*waa
  vwcaa <- apply(wcaa,c(2,3),sum,na.rm=T)

res <- list(beta=beta,delta=delta,alpha=alpha,Blim=Blim,Bban=Bban,waa=waa,maa=maa,M=M,naa=naa,baa=baa,ssb=ssb,faa=faa,caa=caa,wcaa=wcaa,vwcaa=vwcaa,eaa=eaa,multi=multi,input=arglist)
}

#

HS <- function(x,a,b,gamma1=0.001,HStype="HS") if (HStype=="Mesnil") a*(x+sqrt(b^2+(gamma1^2)/4)-sqrt((x-b)^2+(gamma1^2)/4)) else ifelse(x > b, a*b, a*x)
