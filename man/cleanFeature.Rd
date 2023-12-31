% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cleanFeature.R
\name{cleanFeature}
\alias{cleanFeature}
\title{Title}
\usage{
cleanFeature(
  peaks,
  chlorine_tb,
  rawlcms,
  rawfile_dir,
  lcmsfile,
  adducts = TRUE,
  isotopes = TRUE,
  inSourceFrag = TRUE,
  ISFtable,
  samNum = 1,
  rt_tol = 10,
  mz_tol = 0.01
)
}
\arguments{
\item{peaks}{data frame, all peaks}

\item{chlorine_tb}{data frame, all Cl-containing peaks}

\item{rawlcms}{lcms object extracted by xcmsRaw}

\item{rawfile_dir}{string, file folder storing the raw mzML/mzXML files}

\item{lcmsfile}{string, mzML/mzXML file}

\item{adducts}{logic, TRUE or FALSE}

\item{isotopes}{logic, TRUE or FALSE}

\item{inSourceFrag}{logic, TRUE or FALSE}

\item{ISFtable}{data frame, the primary in source fragment table}

\item{samNum}{numeric, sample number}

\item{rt_tol}{numeric, retention time in seconds}

\item{mz_tol}{numeric, retention time in seconds}
}
\value{
data frame with all cl features
}
\description{
Title
}
