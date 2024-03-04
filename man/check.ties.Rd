\name{check.ties}
\alias{check.ties}
\title{Report score ties}
\usage{
  check.ties(first,last)
}
\arguments{
  \item{first}{named vector including predictions when setting ties.method = "first"}

  \item{last}{named vector including predictions when setting ties.method = "last" }

}
\value{
  Message indicating sample with score tie and the highest scored and reported subtypes
}
\description{
  Checks score ties in the prediction scores and prints a message the sample
  where the tie occurred, the two subtypes with the tied scores and the subtype
  that is reported in the output object
}
