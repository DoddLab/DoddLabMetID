% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FunMetaboliteId.R
\name{readMSP}
\alias{readMSP}
\title{readMSP}
\usage{
readMSP(
  file,
  mode = c("all", "standard"),
  source = c("MetAnalyzer", "MSDIAL", "Other")
)
}
\arguments{
\item{file}{the file name}

\item{mode}{standard: extract name, mz and RT from the file, which fit for MSP data exported from QI software; all: extract all information in the file}
}
\value{
A list object. info: MS1 information, including ms1, rt etc.; spec: the MS/MS spectrum
}
\description{
read MSP spectra files
}
\details{
#' @title readMSP
#' @description  Read a MSP file and return a list of spectra for all feature with feature information
#' @param file path of the msp file
#' @export

setGeneric('readMSP', function(file) {
  msp.data.list <- ListDB(file)
  nr.num.pk <- grep('Num Peaks', msp.data.list[[1]])
  info.spec <- lapply(msp.data.list, function(msp.data) {
    info.list <- strsplit(msp.data[1:nr.num.pk], split = ': ')
    info <- do.call(cbind, info.list)
    colnames(info) <- info[1, ]
    mz <- round(as.numeric(info[2,3]), digits = 4)
    rt <- round(as.numeric(strsplit(info[2,"Comment"], split = "_")[[1]][1])*60, digits = 0)
    name <- info[2,"Comment"]

    # info <- data.frame(mz=mz, rt=rt)
    rownames(info) <- name

    spec.list <- strsplit(msp.data[(nr.num.pk+1):length(msp.data)], ' |\\t')
    spec.list <- lapply(spec.list, as.numeric)
    spec <- do.call(rbind, spec.list)[, c(1, 2), drop = FALSE]
    colnames(spec) <- c('mz', 'intensity')

    list('info' = info,
         'spec' = spec)

  })

  info.spec
})

setGeneric('ListDB', function(file) {
  msp.data <- readLines(file)
  nl.db.new <- 1
  idx.db <- 1
  db.list <- list()
  len.data <- length(msp.data)
  for(nl in 1:len.data)
  {
    if(msp.data[nl]=="") {
      db.list[idx.db] <- list(Compound=msp.data[nl.db.new:(nl-1)])
      nl.db.new <- nl + 1
      idx.db <- idx.db + 1
    } else if (nl == len.data){
      db.list[idx.db] <- list(Compound=msp.data[nl.db.new:nl])
    }
  }
  db.list
})
}
\examples{
test <- readMSP(file = 'F:/01 MetIMMS/00 data processing/190515 external validation msms data extraction/zhumetlib_validation_pos_20v_190520.msp', mode = 'all')
}
\author{
Zhiwei Zhou
\email{zhouzw@sioc.ac.cn}
}
