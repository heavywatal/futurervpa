
Generation.Time <- function(vpares,
  maa.year=2014:2015,
  M.year=2014:2015,
  Plus = 19
){

  maa <- vpares$input$dat$maa
  maa <- rowMeans(maa[,colnames(maa) %in% maa.year])
  M <- vpares$input$dat$M
  M <- rowMeans(M[,colnames(M) %in% M.year])

  age <- as.numeric(names(maa))

  maa <- c(maa, rep(1,Plus))
  M <- c(M, rep(M[length(M)],Plus))

  age <- c(age, max(age)+1:Plus)

  A <- length(M)

  L <- c(1,exp(-cumsum(M[-A])))

  G <- sum(age*L*maa)/sum(L*maa)

  return(G)
}
