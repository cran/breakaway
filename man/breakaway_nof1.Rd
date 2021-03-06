\name{breakaway_nof1}
\alias{breakaway_nof1}
\title{ species richness estimation without singletons }
\description{This function permits estimation of total diversity based on a sample frequency count table. Unlike \code{\link{breakaway}}, it does not require an input for the number of species observed once, making it an excellent exploratory tool for microbial ecologists who believe that their sample may contain spurious singletons. The underlying estimation procedure is similar to that of \code{\link{breakaway}} and is outlined in Willis & Bunge (2014). The diversity estimate, standard error, estimated model coefficients and plot of the fitted model are returned.}
\usage{breakaway_nof1(data, print = TRUE, plot = TRUE, answers = FALSE, force = FALSE)}
\arguments{
  \item{data}{
  	The sample frequency count table for the population of interest. The first row must correspond to the doubletons. Acceptable formats include a matrix, data frame, or file path (csv or txt). The standard frequency count table format is used: two columns, the first of which contains the frequency of interest (eg. 1 for singletons, species observed once; 2 for doubletons, species observed twice, etc.) and the second of which contains the number of species observed this many times. Frequencies (first column) should be ordered least to greatest. At least 6 contiguous frequencies are necessary. Do not concatenate large frequencies. See dataset \code{\link{apples}} for sample formatting.
 }
  \item{print}{
Logical: whether the results should be printed to screen. If \samp{FALSE}, \samp{answers} should be set to \samp{TRUE} so that results will be returned.
}
  \item{plot}{
Logical: whether the data and model fit should be plotted.
}
  \item{answers}{
Logical: whether the function should return an argument. If \samp{FALSE}, \samp{print} should be set to \samp{TRUE}.
}
  \item{force}{
Logical: force the procedure to run in the presence of frequency count concatenation. A basic check procedure confirms that the user has not appeared to concatenate multiple upper frequencies. \samp{force=TRUE} will force the procedure to fit models in the presence of this. \samp{breakaway_nof1}'s diversity estimates cannot be considered reliable in this case.
	}
}
\value{
	\item{code}{
		A category representing algorithm behaviour. \samp{code=1} indicates no nonlinear models converged and the transformed WLRM diversity estimate of Rocchetti et. al. (2011) is returned. \samp{code=2} indicates that the iteratively reweighted model converged and was returned. \samp{code=3} indicates that iterative reweighting did not converge but a model based on a simplified variance structure was returned (in this case, the variance of the frequency ratios is assumed to be proportional to the denominator frequency index). Please peruse your fitted model before using your diversity estimate. 
		}
	\item{name}{
		The ``name'' of the selected model. The first integer represents the numerator polynomial degree and the second integer represents the denominator polynomial degree. See Willis & Bunge (2014) for details.
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
		The chosen nonlinear model for frequency ratios.
		}
}
\references{
	  Willis, A. (2015). Species richness estimation with high diversity but spurious singletons. \emph{Under review.}

	Willis, A. and Bunge, J. (2015). Estimating diversity via frequency ratios. \emph{Biometrics.}
}
\author{Amy Willis}
\note{ It is common for microbial ecologists to believe that their dataset contains false diversity. This often arises because sequencing errors result in classifying already observed organisms as new organisms. \samp{breakaway_nof1} was developed as an exploratory tool in this case. Practitioners can run \samp{breakaway} on their dataset including the singletons, and \samp{breakaway_nof1} on their dataset excluding the singletons, and assess if the estimated levels of diversity are very different. Great disparity may provide evidence of an inflated singleton count, or at the very least, that \samp{breakaway} is especially sensitive to the number of rare species observed. Note that \samp{breakaway_nof1} may be less stable than \samp{breakaway} due to predicting based on a reduced dataset, and have greater standard errors.
	}
\seealso{ \code{\link{breakaway}}; \code{\link{apples}} }
\examples{breakaway_nof1(apples[-1,])
breakaway_nof1(apples[-1,],plot=FALSE,print=FALSE,answers=TRUE) 
}
\keyword{ diversity }
\keyword{ microbial }
\keyword{ error }
\keyword{ models }
\keyword{ nonlinear }