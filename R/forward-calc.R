forward.calc.simple <- function(fav,nav,Mv,plus.group=TRUE){
  nage <- max(which(!is.na(nav)))#length(fav)
  naa <- rep(NA,nage)
#  for(a in 2:(nage-1)){
    naa[c(-1,-nage)] <- nav[c(-nage,-(nage-1))]*exp(-fav[c(-nage,-(nage-1))]-Mv[c(-nage,-(nage-1))])
#  }
  naa[nage] <- nav[nage-1]*exp(-fav[nage-1]-Mv[nage-1])
  pg <- nav[nage]*exp(-fav[nage]-Mv[nage])
  if(plus.group) naa[nage] <- naa[nage] + pg
  return(naa)
}

forward.calc.mat <- function(fav,nav,Mv,plus.group=TRUE){
  nage <- max(which(!is.na(nav[,1])))#length(fav)
  na.age <- which(is.na(nav[,1]))
#  naa <- matrix(NA,nage,dim(nav)[[2]])
  naa <- matrix(NA,dim(nav)[[1]],dim(nav)[[2]])
#  for(a in 2:(nage-1)){
    naa[c(-1,-nage,-na.age),] <- nav[c(-nage,-(nage-1),-na.age),]*
        exp(-fav[c(-nage,-(nage-1),-na.age),]-Mv[c(-nage,-(nage-1),-na.age),])
#  }
  naa[nage,] <- nav[nage-1,]*exp(-fav[nage-1,]-Mv[nage-1,])
  pg <- nav[nage,]*exp(-fav[nage,]-Mv[nage,])
  if(plus.group) naa[nage,] <- naa[nage,] + pg
  return(naa)
}
