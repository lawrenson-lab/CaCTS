% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/JSD.R
\name{run_CaCTS_score}
\alias{run_CaCTS_score}
\title{Calculate The Jensen-Shannon Divergence considering expression matrix of representative samples}
\usage{
run_CaCTS_score(rep.matrix, query.name)
}
\arguments{
\item{rep.matrix}{is a expression matrix of the representative samples from prepare_representaive_samples}

\item{query.index}{numeric index of the cancer type to be evaluated as query}
}
\value{
Return the list of CaCTS scores assigned to each TF.
}
\description{
This is a implementation of the Jensen-Shannon Divergence for Transcription Factors expression
    across diffrent cancer types. The expectation is that you have cancer cell-type expression specificity.
    That are considered as CaCTS scores.
}
