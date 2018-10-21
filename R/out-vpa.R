#' 結果の出力
#'
#' @rdname out-vpa
#' @export
out.vpa <- function(res=NULL, # VPA result
                    rres=NULL, # reference point
                    fres=NULL, # future projection result (not nessesarily)
                    ABC=NULL,
                    filename="vpa" # filename without extension
                    ){
  old.par <- par()
  exit.func <- function(){
#    par(old.par)
    dev.off()
    options(warn=0)
  }
  on.exit(exit.func())

  csvname <- paste(filename,".csv",sep="")
  pdfname <- paste(filename,".pdf",sep="")
  pdf(pdfname)
  par(mfrow=c(3,2),mar=c(3,3,2,1))
  options(warn=-1)

  write.table2 <- function(x,title.tmp="",is.plot=TRUE,...){
    if(is.plot){
	    if(!is.null(dim(x))){
    	  matplot(colnames(x),t(x),type="b",ylim=c(0,max(x,na.rm=T)),pch=substr(rownames(x),1,1))
	    }
    	else{
	      barplot(x)
    	}
    title(title.tmp)
    }
    if(!is.null(dim(x))){
      tmp <- matrix("",nrow(x)+1,ncol(x)+1)
      tmp[-1,-1] <- as.character(unlist(x))
      tmp[-1,1] <- rownames(x)
      tmp[1,-1] <- colnames(x)
    }
    else{
      tmp <- x
    }
    write.table(tmp,append=T,sep=",",quote=FALSE,file=csvname,col.names=F,row.names=F,...)
  }

  write(paste("# RVPA outputs at ",date()," & ",getwd()),file=csvname)

  if(!is.null(res)){
    write("# VPA results",file=csvname, append=T)

    write("\n# catch at age",file=csvname,append=T)
    write.table2(res$input$dat$caa,title.tmp="Catch at age")

    write("\n# maturity at age",file=csvname,append=T)
    write.table2(res$input$dat$maa,title.tmp="Maturity at age")

    write("\n# weight at age for biomass calculation",file=csvname,append=T)
    write.table2(res$input$dat$waa,title.tmp="Weight at age (for biomass)")

    if(!is.null(res$input$dat$waa.catch)){
    write("\n# weight at age for catch calculation",file=csvname,append=T)
    write.table2(res$input$dat$waa.catch,title.tmp="Weight at age (for catch)")
    }


    write("\n# M at age",file=csvname,append=T)
    write.table2(res$input$dat$M,title.tmp="M at age")

    write("\n# fishing mortality at age",file=csvname,append=T)
    write.table2(res$faa,title.tmp="F at age")

    write("\n# Current F",file=csvname,append=T)
    write.table2(res$Fc.at.age,title.tmp="Current F")

    write("\n# numbers at age",file=csvname,append=T)
    write.table2(res$naa,title.tmp="Numbers at age")

    write("\n# total and spawning biomass ",file=csvname,append=T)
    x <- rbind(colSums(res$ssb),colSums(res$baa),colSums(res$wcaa))
    rownames(x) <- c("Spawning biomass","Total biomass","Catch biomass")
    write.table2(x,title.tmp="Total and spawning biomass")
  }

  if(!is.null(rres)){
    write("\n# Reference points",file=csvname,append=T)
    write.table2(rres$summary,title.tmp="Future F at age",is.plot=F)
  }

  if(!is.null(fres)){
    write("\n# future projection results",file=csvname,append=T)
    write("\n# future F at age",file=csvname,append=T)
    write.table2(fres$faa[,,1],title.tmp="Future F at age")

    write("\n# future numbers at age",file=csvname,append=T)
    write.table2(fres$naa[,,1],title.tmp="Future numbers at age")

    write("\n# future total and spawning biomass",file=csvname,append=T)
    x <- rbind(fres$vssb[,1],fres$vbiom[,1],fres$vwcaa[,1])
    rownames(x) <- c("Spawning biomass","Total biomass","Catch biomass")
    write.table2(x,title.tmp="Future total, spawning and catch biomass")
  }

  if(!is.null(ABC)){
    write("\n# ABC summary",file=csvname,append=T)
    write.table2(ABC$ABC,title.tmp="Future F at age",is.plot=F)
    write("\n# Kobe matrix",file=csvname,append=T)
    for(i in 1:dim(ABC$kobe.matrix)[[3]]){
        write(paste("\n# ",dimnames(ABC$kobe.matrix)[[3]][i]),
              file=csvname,append=T)
        write.table2(ABC$kobe.matrix[,,i],
                     title.tmp=dimnames(ABC$kobe.matrix)[[3]][i],is.plot=T)
    }
  }
}
