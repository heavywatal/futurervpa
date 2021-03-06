% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est-msy.R
\name{est.MSY}
\alias{est.MSY}
\title{Estimate MSY}
\usage{
est.MSY(vpares, farg, seed = 1, n.imputation = 1, nyear = 50,
  eyear = 0, FUN = mean, N = 1000, onlylower.pgy = FALSE,
  is.small = FALSE, is.Kobe = FALSE, is.5perlower = FALSE,
  optim.method = "optimize", max.target = "catch.mean",
  calc.yieldcurve = TRUE, Blimit = 0, trace.multi = c(seq(from = 0,
  to = 0.9, by = 0.1), 1, seq(from = 1.1, to = 2, by = 0.1), 3:5, 7, 20,
  100), is.plot = TRUE, PGY = NULL, B0percent = NULL)
}
\arguments{
\item{eyear}{将来予測の最後のeyear+1年分を平衡状態とする}

\item{FUN}{漁獲量の何を最大化するか？}

\item{N}{stochastic計算するときの繰り返し回数}

\item{onlylower.pgy}{PGY計算するとき下限のみ計算する（計算時間省略のため）}

\item{is.small}{将来予測の結果を返さない。}

\item{is.Kobe}{Kobeの計算をするかどうか。順番に、HS, BH, RIの順}

\item{is.5perlower}{HSの折れ点を5\%の確率で下回るときの親魚資源量}

\item{max.target}{method="optimize"以外を使うとき、どの指標を最大化するか。
他のオプションとしては"catch.median" (漁獲量のmedianの最大化)}

\item{calc.yieldcurve}{yield curveを正確に計算するかどうか。
TRUEだと計算時間が余計にかかる。FALSEだと、yield curveは正確ではない}

\item{trace.multi}{Fmsyを探索したり、Yield curveを書くときにグリッドサーチをするときのFの刻み。
Fcurrentに対する乗数。Fが異常に大きい場合、親魚=0になって加入＝NA}

\item{PGY}{PGY管理基準値を計算するかどうか。
計算しない場合はNULLを、計算する場合はc(0.8,0.9,0.95)のように割合を入れる}

\item{B0percent}{B0_XX\%の管理基準値を計算するかどうか}
}
\description{
Estimate MSY
}
