\name{breakaway-package}
\alias{breakaway-package}
\docType{package}
\title{ Species richness estimation and modelling in the high-diversity setting }
\description{ \samp{breakaway} estimates total (observed plus unobserved) species richness. Microbial diversity datasets are characterized by a large number of rare species and a small number of highly abundant species. The class of models implemented by \samp{breakaway} is flexible enough to model both these features. \samp{breakaway_nof1} implements a similar procedure however does not require a singleton count. \samp{betta_basic} provides a method for modelling total diversity with covariates in a way that accounts for its estimated nature and thus permits accounting for unobserved taxa. }
\details{
\tabular{ll}{
Package: \tab breakaway\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-08-03\cr
License: \tab GPL-2 \cr
}
The function \code{\link{breakaway}} estimates the total (observed plus unobserved) number of classes (usually, distinct species) based on a sample of the frequency counts. Standard errors and model fits are also given. The algorithm is based on theory of characterization of distributions by ratios of their probabilities. Parameter estimation is done via nonlinear regression. The class of models available is usually broad enough to account for the high-diversity case, which is often observed in microbial diversity datasets. Since many classical estimation procedures either fail to provide an estimate or provide poor fits in the microbial setting, \code{\link{breakaway}} addresses this data structure. Additionally, since sequencing errors may result in an inflated singleton count, \code{\link{breakaway_nof1}} performs a similar procedure but does not require a singleton count. It can be used as an exploratory tool for investigating the plausibility of the given singleton count. \code{\link{betta_basic}} runs a regression-type analysis of estimated total diversity, thus permitting accounting for unobserved taxa. It does not enforce use of \samp{breakaway} for diversity estimation. A mixed-model approach accounts for the differing levels of confidence in the diversity estimates, and covariates constitute the fixed effects. Support of this work from Cornell University's Department of Statistical Sciences is gratefully acknowledged.
}
\author{
Amy Willis & John Bunge

Maintainer: Amy Willis <adw96@cornell.edu>
}
\references{
	Willis, A. and Bunge, J. (2014). A ratio-based method for estimating an unknown number of classes. Manuscript submitted for publication.
  
	Rocchetti, I., Bunge, J. and Bohning, D. (2011). Population size estimation based upon ratios of recapture probabilities. \emph{Annals of Applied Statistics}, \bold{5}.
  
	Willis, A. and Bunge, J. (2014). A first look at beta diversity. Manuscript in preparation.
}
\examples{
breakaway(apples)
breakaway(apples,plot=FALSE,print=FALSE,answers=TRUE)
breakaway_nof1(apples[-1,])
betta_basic(c(1552,1500,884),c(305,675,205))
}
