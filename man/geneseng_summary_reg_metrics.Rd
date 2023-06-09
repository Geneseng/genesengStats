% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geneseng_summary_reg_metrics.R
\name{geneseng_summary_reg_metrics}
\alias{geneseng_summary_reg_metrics}
\title{Summary of regression metrics}
\usage{
geneseng_summary_reg_metrics(data, group)
}
\arguments{
\item{data}{provide a data.frame}

\item{group}{group by a continuous variable}
}
\description{
Evaluate biomarker individual performances using best regression metrics
}
\details{
`Linear regression` was used for regression problem.
}
\examples{
library(genesengStats)
geneseng_summary_reg_metrics(data = mtcars, group = "mpg")
geneseng_summary_reg_metrics(data = iris, group = "Sepal.Length")

}
