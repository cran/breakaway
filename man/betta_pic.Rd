\name{betta_pic}
\alias{betta_pic}
\title{ function for plotting total diversity  }
\description{A simple plotting interface for comparing total diversity across samples or a covariate gradient.}
\usage{betta_pic(y, se, x, ylimu, myy = NA, mymain = NA, 
          mycol = rep("black", length(y)), labs = NA, 
          mypch = rep(16, length(y)),
          myxlim = c(0, 1.1*max(x,na.rm=T)))}
\arguments{
  \item{y}{
  	A vector of estimates of total diversity. Other parameter estimates are accessible; this method may be used for plotting any parameter estimates.. 
	 }
  \item{se}{
	The standard errors in \samp{y}, the diversity (or other parameter's) estimates. 
	}
  \item{x}{
	A vector of covariates to form the x-coordinates of the intervals.
	}
	\item{ylimu}{
	The upper endpoint of the y-axis.
	}
	\item{myy}{
	A label for the y-axis. Default is none.
	}
	\item{mymain}{
	A main title for the plot. Default is none.
	}
	\item{mycol}{
	Colors for the plotting points. Default is black.
	}
	\item{labs}{
	x-axis labels. Default is none (non x-axis plotted).
	}
	\item{mypch}{
	Plotting characters for the estimates. Defaults to circles.
	}
	\item{myxlim}{
	A vector of x-axis limits. Default is 0 to 10\% greater than the maximum x value.
	}
}
\references{
	Willis, A., Bunge, J., and Whitman, T. (2015). Inference for changes in biodiversity. arXiv preprint. 
}
\author{Amy Willis}
\seealso{ \code{\link{breakaway}}; \code{\link{breakaway_nof1}}; \code{\link{apples}}; \code{\link{hawaii}} }
\examples{## test for heterogeneity of butterfly diversity estimates vs apples estimates
betta_pic(c(1552,1500,884),c(305,675,205),c(1,2,3),3000)
	}
\keyword{ diversity }