\name{wlrm_untransformed}
\alias{wlrm_untransformed}
\title{ function for species richness estimation }
\description{This function implements the untransformed version of the species richness estimation procedure outlined in Rocchetti, Bunge and Bohning (2011).}
\usage{wlrm_untransformed(data, print=TRUE, plot=FALSE, answers=FALSE, cutoff=NA)}
\arguments{
  \item{data}{
    The sample frequency count table for the population of interest. See \samp{breakaway} for details.
  }
  \item{print}{
    Logical: whether the results should be printed to screen.
  }
  \item{plot}{
    Logical: whether the data and model fit should be plotted.
  }
  \item{answers}{
    Logical: whether the function should return an argument.
  }
  \item{cutoff}{
    The maximum frequency count to use for prediction. Defaults to maximal.
  }}
\value{
  \item{name}{
    The ``name'' of the model.
  }
  \item{para}{
    Estimated model parameters and standard errors.
  }
  \item{est}{
    The estimate of total (observed plus unobserved) diversity.
  }
  \item{seest}{
    The standard error in the diversity estimate.
  }
  \item{full}{
    The chosen linear model for frequency ratios.
  }
  \item{ci}{
    An asymmetric 95\% confidence interval for diversity.
  }
}
\references{
  Rocchetti, I., Bunge, J. and Bohning, D. (2011). Population size estimation based upon ratios of recapture probabilities. \emph{Annals of Applied Statistics}, \bold{5}.
}
\author{Amy Willis}
\note{ This estimator is based on the negative binomial model and for that reason generally produces poor fits to microbial data. The result is usually artificially low standard errors. Caution is advised.
}
\seealso{ \code{\link{breakaway}}; \code{\link{apples}}; \code{\link{wlrm_transformed}} }
\examples{wlrm_untransformed(apples)}
\keyword{ diversity }
\keyword{ models }
