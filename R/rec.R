# HS用; ARには対応していないが、残差リサンプリングには対応している
HS.rec <- function(ssb,vpares,#deterministic=FALSE,
                   rec.resample=NULL,
                   rec.arg=list(a=1000,b=1000,sd=0.1, # Mesnil関数のparameter
                                        resample=FALSE,resid=0, # 残差リサンプリングする場合、resample=TRUEにして、residにリサンプリングする残差（対数）を入れる
                                        bias.correction=TRUE)){

    rec0 <- ifelse(ssb>rec.arg$b,rec.arg$a*rec.arg$b,rec.arg$a*ssb)
    if(!isTRUE(rec.arg$resample)){
        if(isTRUE(rec.arg$bias.correction)){
            rec <- rec0*exp(rnorm(length(ssb),-0.5*(rec.arg$sd)^2,rec.arg$sd))
        }
        else{
            rec <- rec0*exp(rnorm(length(ssb),0,rec.arg$sd))
        }
    }
    else{
        if(isTRUE(rec.arg$bias.correction)){
            rec <- c(rec0[1],exp(log(rec0[-1])+sample(rec.arg$resid,length(ssb)-1,replace=TRUE))/mean(exp(rec.arg$resid)))
        }
        else{
            rec <- c(rec0[1],exp(log(rec0[-1])+sample(rec.arg$resid,length(ssb)-1,replace=TRUE)))
        }
    }
  return(list(rec=rec,rec.resample=rec.arg$resid)) # 暫定的変更
}

# RI用; ARには対応していないが、残差リサンプリングには対応している
RI.rec <- function(ssb,vpares,#deterministic=FALSE,
                   rec.resample=NULL,
                   rec.arg=list(a=1000,b=1000,sd=0.1, # Mesnil関数のparameter
                                        resample=FALSE,resid=0, # 残差リサンプリングする場合、resample=TRUEにして、residにリサンプリングする残差（対数）を入れる
                                        bias.correction=TRUE)){

    rec0 <- rec.arg$a*ssb*exp(-rec.arg$b*ssb) # rec.arg$a*ssb/(1+rec.arg$b*ssb)
#    rec0 <- ifelse(ssb>rec.arg$b,rec.arg$a*rec.arg$b,rec.arg$a*ssb)
    if(!isTRUE(rec.arg$resample)){
        if(isTRUE(rec.arg$bias.correction)){
            rec <- rec0*exp(rnorm(length(ssb),-0.5*(rec.arg$sd)^2,rec.arg$sd))
        }
        else{
            rec <- rec0*exp(rnorm(length(ssb),0,rec.arg$sd))
        }
    }
    else{
        if(isTRUE(rec.arg$bias.correction)){
            rec <- c(rec0[1],exp(log(rec0[-1])+sample(rec.arg$resid,length(ssb)-1,replace=TRUE))/mean(exp(rec.arg$resid)))
        }
        else{
            rec <- c(rec0[1],exp(log(rec0[-1])+sample(rec.arg$resid,length(ssb)-1,replace=TRUE)))
        }
    }
  return(list(rec=rec,rec.resample=rec.arg$resid)) # 暫定的変更
}


# RI用; ARには対応していないが、残差リサンプリングには対応している
BH.rec <- function(ssb,vpares,#deterministic=FALSE,
                   rec.resample=NULL,
                   rec.arg=list(a=1000,b=1000,sd=0.1, # Mesnil関数のparameter
                                        resample=FALSE,resid=0, # 残差リサンプリングする場合、resample=TRUEにして、residにリサンプリングする残差（対数）を入れる
                                        bias.correction=TRUE)){
    rec0 <- rec.arg$a*SSB/(1+rec.arg$b*SSB)
#    rec0 <- ifelse(ssb>rec.arg$b,rec.arg$a*rec.arg$b,rec.arg$a*ssb)
    if(!isTRUE(rec.arg$resample)){
        if(isTRUE(rec.arg$bias.correction)){
            rec <- rec0*exp(rnorm(length(ssb),-0.5*(rec.arg$sd)^2,rec.arg$sd))
        }
        else{
            rec <- rec0*exp(rnorm(length(ssb),0,rec.arg$sd))
        }
    }
    else{
        if(isTRUE(rec.arg$bias.correction)){
            rec <- c(rec0[1],exp(log(rec0[-1])+sample(rec.arg$resid,length(ssb)-1,replace=TRUE))/mean(exp(rec.arg$resid)))
        }
        else{
            rec <- c(rec0[1],exp(log(rec0[-1])+sample(rec.arg$resid,length(ssb)-1,replace=TRUE)))
        }
    }
  return(list(rec=rec,rec.resample=rec.arg$resid)) # 暫定的変更
}


# Hockey-stick(bias.correctionのオプションは削除。どうせするので）
HS.recAR <- function(ssb,vpares,#deterministic=FALSE,
                      rec.resample=NULL,
                      rec.arg=list(a=1000,b=1000,#gamma=0.01,
                                   sd=0.1, rho=0,
                                   resid=0)#, bias.correction=TRUE)
                      ){
    ## 再生産関係からの予測値
#    rec0 <- rec.arg$a*(ssb+sqrt(rec.arg$b^2+(rec.arg$gamma^2)/4)-sqrt((ssb-rec.arg$b)^2+(rec.arg$gamma^2)/4))
    rec0 <- ifelse(ssb>rec.arg$b,rec.arg$a*rec.arg$b,rec.arg$a*ssb)
    rec <- rec0*exp(rec.arg$rho*rec.arg$resid) # 自己相関込みの予測値

    rec <- rec*exp(rnorm(length(ssb),-0.5*rec.arg$sd2^2,rec.arg$sd))
    new.resid <- log(rec/rec0)+0.5*rec.arg$sd2^2
    return(list(rec=rec,rec.resample=new.resid))
}


# Beverton-Holt
BH.recAR <- function(ssb,vpares,deterministic=FALSE,rec.resample=NULL,
                   rec.arg=list(a=1000,b=1000,sd=0.1,bias.correction=TRUE)){
  rec0 <- rec.arg$a*ssb/(1+rec.arg$b*ssb)
  rec <- rec0*exp(rec.arg$rho*rec.arg$resid) # 自己相関込みの予測値
  rec <- rec*exp(rnorm(length(ssb),-0.5*rec.arg$sd2^2,rec.arg$sd))
  new.resid <- log(rec/rec0)+0.5*rec.arg$sd2^2
  return(list(rec=rec,rec.resample=new.resid))
}

# Ricker
RI.recAR <- function(ssb,vpares,deterministic=FALSE,rec.resample=NULL,
                   rec.arg=list(a=1000,b=1000,sd=0.1,bias.correction=TRUE)){
    rec0 <- rec.arg$a*ssb*exp(-rec.arg$b*ssb) # rec.arg$a*ssb/(1+rec.arg$b*ssb)
    rec <- rec0*exp(rec.arg$rho*rec.arg$resid) # 自己相関込みの予測値
    rec <- rec*exp(rnorm(length(ssb),-0.5*rec.arg$sd2^2,rec.arg$sd))
    new.resid <- log(rec/rec0)+0.5*rec.arg$sd2^2
    return(list(rec=rec,rec.resample=new.resid))
}
