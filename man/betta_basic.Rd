\name{betta_basic}
\alias{betta_basic}
\title{ function for modelling total diversity  }
\description{This function tests for heterogeneity of total diversity (observed plus unobserved) across multiple sites. It can account or test for fixed effects that may explain diversity. It returns the significance of the covariates in explaining diversity and a hypothesis test for heterogeneity.}
\usage{ betta_basic(Chats,seChats,X=NA) }
\arguments{
  \item{Chats}{
  	A vector of estimates of total diversity at different sampling locations. \samp{breakaway} estimates are suggested but not enforced.
	 }
  \item{seChats}{
	The standard errors in \samp{Chats}, the diversity estimates. 
	}
  \item{X}{
	A numeric matrix of covariates. If not supplied, an intercept-only model will be fit.
	}
}
\value{
	\item{table}{
		A coefficient table for the model parameters. The columns give the parameter estimates, standard errors, and p-values, respectively. Note that the significance tests are Wald-type (since the procedure is unrestricted maximum-likelihood) and one-sided for the random effects variance and two-sided for all others. Note that since the procedure is unrestricted, small variances and small samples may results\ in negative estimates for the random effect variance. Proceed with caution. This model is only as effective as your diversity estimation procedure; for this reason please confirm that your estimates are appropriate and that your model is well-specified. 
		}
	\item{estimates}{
		The parameter estimates. The first estimate is the estimated variance of the random effects, followed by the intercept estimate, followed by the covariates (if any) in the order they were supplied.
		}
	\item{variances}{
		The estimated variances in the parameter estimates.
		}
	\item{blups}{
		The estimated best linear unbiased predictors (eBLUPs). These are the realizations of the random variables corresponding to the random effects, so are normally distributed with zero mean and variance estimated by \samp{sigsq_u}. 
		}
	\item{y}{
		The supplied diversity estimates (Chats).
		}
	\item{yvars}{
		The supplied estimated variances in the diversity estimates.
		}
	\item{condfits}{
		The conditional expected values of the diversity estimates (conditional on the random effects). Sum of the linear term and the eBLUPs. The authors propose that if the practitioner believes that information from one diversity estimator may inform the others, then using the \samp{condfits} as estimators of total diversity rather than \samp{Chats} may reduce variance of diversity estimates by ``sharing strength'' across the samples. 
		}
	\item{interval}{
		A Wald-type 95\% confidence interval for the random effects variance. An interval containing zero suggests homogeneity of the sample. 
		}
}
\references{
	Willis, A. and Bunge, J. (2014). A first look at beta diversity. Manuscript in preparation.
  
	Willis, A. and Bunge, J. (2014). A ratio-based method for estimating an unknown number of classes. Manuscript submitted for publication.
}
\author{Amy Willis}
\note{ Ecologists who are interested in the way species richness varies with covariate information often run a regression-type analysis on the observed diversity using their covariate information as predictors. However, in many settings (especially microbial), rare and unobserved taxa play a hugely important role in explaining the subtleties of the ecosystem, however, a regression analysis on the observed diversity level fails to account for these unobserved taxa. By predicting the total level of diversity (for example, via \code{\link{breakaway}}) and estimating the standard error in the estimate, one can take account of these unobserved, but important, taxa. In order to account for the estimated nature of the response, a mixed model approach is taken, whereby the varying levels of confidence in the estimates contributes to a diagonal but heteroscedastic covariance matrix. Given covariates constitute the fixed effects in the mixed model, and significance of the random effect term ``sigsq_u'' reflects heterogeneity in the sample, that is, variability that cannot be explained by only the covariates. The authors believe this to be the first attempt at modelling total diversity in a way that accounts for its estimated nature. 
		}
\seealso{ \code{\link{breakaway}};\code{\link{breakaway_nof1}};\code{\link{apples}} }
\examples{
	## test for heterogeneity of butterfly diversity estimates vs apples estimates
	betta_basic(c(1552,1500,884),c(305,675,205))
	}
\keyword{ diversity }