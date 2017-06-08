\name{IlluminaHumanMethylationEPICmanifest}
\alias{IlluminaHumanMethylationEPICmanifest}
\docType{data}
\title{
  Manifest object for the IlluminaHumanMethylationEPIC array
}
\description{
  This object contains the array design for Illumina's Human Methylation
  EPIC microarray. Based on the v1.0b2 version of the manifest file.
}
\usage{data(IlluminaHumanMethylationEPICmanifest)}
\details{
  This package is based on the file
  \sQuote{MethylationEPIC_v-1-0_B4.csv} from the zip archive
  \sQuote{infinium-methylationepic-v1-0-b4-manifest-file-csv.zip} from
  \url{www.illumina.com}.

  A script for creating the data object in this package is contained in
  \code{scripts/manifest.R}.
}
\format{
  An object of class \code{IlluminaMethylationManifest}
}
\references{
  None.
}
\seealso{
  The description of the class in the minfi package,
  \code{\link[minfi]{IlluminaMethylationManifest-class}}.
}
\examples{
data(IlluminaHumanMethylationEPICmanifest)
}
\keyword{datasets}
