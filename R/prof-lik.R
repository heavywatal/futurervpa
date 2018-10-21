### profile likelihood
prof.lik <- function(Res,a=Res$pars$a,b=Res$pars$b,sd=Res$pars$sd,rho=Res$pars$rho) {
  SRdata <- Res$input$SRdata
  rec <- SRdata$R
  ssb <- SRdata$SSB
  N <- length(rec)
  SR <- Res$input$SR
  gamma <- Res$gamma
  method <- Res$input$method
  w <- Res$input$w

#  if (SR=="HS") SRF <- function(x,a,b) a*(x+sqrt(b^2+gamma^2/4)-sqrt((x-b)^2+gamma^2/4))
  if (SR=="HS") SRF <- function(x,a,b) ifelse(x>b,b*a,x*a)
  if (SR=="BH") SRF <- function(x,a,b) a*x/(1+b*x)
  if (SR=="RI") SRF <- function(x,a,b) a*x*exp(-b*x)

  resid <- sapply(1:N,function(i) log(rec[i]) - log(SRF(ssb[i],a,b)))
  resid2 <- NULL
  for (i in 1:N) {
    resid2[i] <- ifelse(i==1,resid[i], resid[i]-rho*resid2[i-1])
  }

  obj <- NULL
  if (method == "L2") {
    for (i in 1:N) {
      if (i==1) {
        obj <- c(obj,-0.5*log(2*pi)-log(sd^2/(1-rho^2))-resid2[i]^2/(2*sd^2/(1-rho^2)))
      } else {
        obj <- c(obj, -0.5*log(2*pi)-0.5*log(sd^2)-resid2[i]^2/(2*sd^2))
      }
    }
  } else {
    for (i in 1:N) {
      if (i==1) {
        obj <- c(obj,-log(2*sqrt(sd^2/(1-rho^2)))-abs(resid2[i])/sqrt(sd^2/(1-rho^2)))
      } else {
        obj <- c(obj, -log(2*sd)-abs(resid2[i])/sd)
      }
    }
  }
  obj <- sum(w*obj) # exact likelihood
  return(exp(obj))
}
