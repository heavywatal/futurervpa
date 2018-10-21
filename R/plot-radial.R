#' Thin wrapper of plotrix::radial.plot
#'
#' @importFrom plotrix radial.plot
#' @importFrom RColorBrewer brewer.pal
#' @rdname plot-radial
#' @export
plotRadial <- function(index,base=1,col.tmp=NULL,lwd=2,...){
    old.par <- par()
    layout(matrix(c(1,2),2,1),heights=c(2,1))

    index2 <- sweep(matrix(unlist(index),nrow(index),ncol(index)),2,as.numeric(unlist(index[base,])),FUN="/")

    if(is.null(col.tmp)) col.tmp <- brewer.pal(nrow(index2-1),"Set1")

    radial.plot(index2,rp.type="p",lwd=lwd,show.grid.labels=FALSE,
                labels=colnames(index),
                radial.lim=c(0,1.5),clockwise=TRUE,start=1,
                line.col=c(NA,col.tmp),
                poly.col=c(rgb(40/255,96/255,163/255,0.2),rep(NA,nrow(index2)-1)), # MSYだけ色で塗る
                ...
                )
    refname <- rownames(index)
    par(mar=c(1,0,1,0))
    plot(0,10,type="n",axes=FALSE,ylab="")
    legend("topleft",legend=refname,
           col=c(rgb(40/255,96/255,163/255,0.2),col.tmp),
           ncol=2,lwd=c(10,rep(lwd,length(refname)-1)))
    layout(matrix(c(1),1,1),heights=c(1))
    par(old.par)
    invisible(index2)
}
