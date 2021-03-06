\name{wlrm_transformed}
\alias{wlrm_transformed}
\title{ function for species richness estimation }
\description{This function implements the transformed version of the species richness estimation procedure outlined in Rocchetti, Bunge and Bohning (2011).}
\usage{wlrm_transformed(data, print=TRUE, plot=FALSE, answers=FALSE, cutoff=NA)}
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
    \note{ While robust to many different structures, model is almost always misspecified. The result is usually implausible diversity estimates with artificially low standard errors. Extreme caution is advised.
    }
  \seealso{ \code{\link{breakaway}}; \code{\link{apples}}; \code{\link{wlrm_untransformed}} }
  \examples{wlrm_transformed(apples)
wlrm_transformed(apples,plot=FALSE,print=FALSE,answers=TRUE)
  }
  \keyword{ diversity }
  \keyword{ models }
