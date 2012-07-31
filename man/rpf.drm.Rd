\name{rpf.drm}
\alias{rpf.drm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
dichotomous response model
}
\description{
Specify an item to be scored using a dichotomous response model (1PL, 2PL, and 3PL)
}
\usage{
rpf.drm(D = 1, numAlternatives = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{
Use 1.7 for the normal-ogive version
}
  \item{numAlternatives}{
%%     ~~Describe \code{numAlternatives} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (D = 1, numAlternatives = 5) 
{
    guess.weight <- 20
    guessing <- (1/numAlternatives)
    new("rpf.drm", numChoices = 2, D = D, a.prior.meanlog = 0, 
        a.prior.sdlog = 0.5, c.prior.alpha = guess.weight * guessing + 
            1, c.prior.beta = guess.weight * (1 - guessing) + 
            1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
