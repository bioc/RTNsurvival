% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllMethods.R
\docType{methods}
\name{tnsCox,TNS-method}
\alias{tnsCox,TNS-method}
\alias{tnsCox}
\title{Cox regression analysis for TNS class objects}
\usage{
\S4method{tnsCox}{TNS}(tns, regs = NULL, qqkeycovar = FALSE, verbose = TRUE)
}
\arguments{
\item{tns}{A \linkS4class{TNS} object, which must have passed GSEA2 analysis.}

\item{regs}{An optional string vector listing regulons to be tested.}

\item{qqkeycovar}{A logical value. If TRUE, only the samples in the 2nd and 
3rd quartils of 'keycovar' are used in the analysis. If FALSE, all samples
are used (see \code{\link{tni2tnsPreprocess}}).}

\item{verbose}{A logical value specifying to display detailed messages 
(when verbose=TRUE) or not (when verbose=FALSE).}
}
\value{
Cox hazard models and statistics.
}
\description{
Run Cox multivariate regression for regulons and other covariates.
}
\examples{
# load survival data
data(survival.data)

# load TNI-object
data(stni, package = "RTN")

stns <- tni2tnsPreprocess(stni, survivalData = survival.data, 
keycovar = c('Age','Grade'), time = 1, event = 2)
stns <- tnsGSEA2(stns)
stns <- tnsCox(stns, regs = c('PTTG1','E2F2','FOXM1'))
tnsGet(stns, "coxTable")

}
