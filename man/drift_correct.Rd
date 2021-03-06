% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/itrackR.R
\name{drift_correct}
\alias{drift_correct}
\title{Perform drift correction on fixations for itrackR object}
\usage{
drift_correct(obj, vars = c("ID"), threshold = 10)
}
\arguments{
\item{obj}{an itrackR object}

\item{vars}{variables to specify subsets of data for performing drift correction. default is \code{'ID'}. Can also name any column in \code{obj$beh} (e.g., 'Block')}

\item{threshold}{do not adjust fixations that deviate less than this threshold (in pixels). Default is 10 pixels}

\item{append}{set to TRUE if you've already added behavioral data and you're adding more columns (default = FALSE, overwrite existing data).}
}
\value{
All fixations in \code{obj$fixations} are adjusted, and adjustment amounts are stored in \code{obj$transform}.
To undo correction, use \code{\link{undrift}}.
}
\description{
Correct for drifts in fixation data. Algorithm computes the median x/y coordinates separately for each subject and behavioral variable specified,
then fixations are adjusted based on difference between this "true" center and the center of the screen based on \code{obj$resolution}. Fixations are only
adjusted if deviation is greater than some threshold.
}
\examples{
\dontrun{
# itrackr.data('edfs') returns full path to 2 edf files
z <- itrackr(edfs=itrackr.data('edfs'))

#Extracts messages like "BLOCK 1", and "TRIAL 7" and creates variables "Block" and "Trial", giving them values
# 1 and 7.
z <- set_index(z,c('Block', 'Trial'), patterns=c('BLOCK [0-9]*', 'TRIAL [0-9][0-9]*), numeric.only = TRUE)

beh <- itrackr.data('beh)

z <- add_behdata(z,beh)

#perform drift correction for each subject and Block.
z <- drift_correct(z, vars=c('Block'))

#if you're curious how adjustment was done, check it out:

View(z$transform)

}


}
\seealso{
\code{\link{undrift}}
}

