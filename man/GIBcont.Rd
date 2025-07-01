\name{GIBcont}
\alias{GIBcont}
\title{Cluster Continuous Data Using the Generalised Information Bottleneck Algorithm}
\description{
The \code{GIBcont} function implements the Generalised Information Bottleneck (GIB) algorithm
for fuzzy clustering of continuous data. This method optimizes an information-theoretic objective to
preserve relevant information while forming concise and interpretable cluster representations
\insertCite{strouse_ib_2019}{IBclust}.
}
\usage{
GIBcont(X, ncl, beta, alpha, randinit = NULL, s = -1, scale = TRUE,
        maxiter = 100, nstart = 100,
        verbose = FALSE)
}
\arguments{
  \item{X}{
    A numeric matrix or data frame containing the continuous data to be clustered. All variables should be of type \code{numeric}.
  }
  \item{ncl}{
    An integer specifying the number of clusters to form.
  }
  \item{beta}{
    Regularisation strength.
  }
  \item{alpha}{
    Strength of relative entropy term.
  }
  \item{randinit}{
    Optional. A vector specifying initial cluster assignments. If \code{NULL}, cluster assignments are initialized randomly.
  }
\item{s}{
A numeric value or vector specifying the bandwidth parameter(s) for continuous variables. The values must be greater than \eqn{0}. The default value is \eqn{-1}, which enables the automatic selection of optimal bandwidth(s).
}
  \item{scale}{
    A logical value indicating whether the continuous variables should be scaled to have unit variance before clustering. Defaults to \code{TRUE}.
  }
  \item{maxiter}{
    The maximum number of iterations allowed for the clustering algorithm. Defaults to \eqn{100}.
  }
  \item{nstart}{
    The number of random initializations to run. The best clustering result (based on the information-theoretic criterion) is returned. Defaults to \code{100}.
  }
  \item{verbose}{
    Logical. Default to \code{FALSE} to suppress progress messages. Change to \code{TRUE} to print.
  }
}

\value{
A list containing the following elements:
 \item{Cluster}{A cluster membership matrix.}
  \item{Entropy}{A numeric value representing the entropy of the cluster assignment, \eqn{H(T)}.}
  \item{RelEntropy}{A numeric value representing the relative entropy of cluster assignment, given the observation weights \eqn{H(X \mid T)}.}
  \item{MutualInfo}{A numeric value representing the mutual information, \eqn{I(Y;T)}, between the original labels (\eqn{Y}) and the cluster assignments (\eqn{T}).}
  \item{beta}{A numeric value of the regularisation strength beta used.}
  \item{alpha}{A numeric value of the strength of relative entropy used.}
  \item{s}{A numeric vector of bandwidth parameters used for the continuous variables.}
  \item{ht}{A numeric vector tracking the entropy value of the cluster assignments across iterations.}
  \item{hy_t}{A numeric vector tracking the relative entropy values between the cluster assignments and observations weights across iterations.}
  \item{iyt}{A numeric vector tracking the mutual information values between original labels and cluster assignments across iterations.}
  \item{losses}{A numeric vector tracking the final loss values across iterations.}
}

\details{
The \code{GIBcont} function applies the Generalised Information Bottleneck algorithm to do fuzzy clustering of datasets comprising only continuous variables. This method leverages an information-theoretic objective to optimize the trade-off between data compression and the preservation of relevant information about the underlying data distribution.
Set \eqn{\alpha = 1} and \eqn{\alpha = 0} to recover the Information Bottleneck and its Deterministic variant, respectively. If \eqn{\alpha = 0}, the algorithm ignores
the value of the regularisation parameter \eqn{\beta}.

The function utilizes the Gaussian kernel \insertCite{silverman_density_1998}{IBclust} for estimating probability densities of continuous features. The kernel is defined as:

\deqn{K_c\left(\frac{x - x'}{s}\right) = \frac{1}{\sqrt{2\pi}} \exp\left\{-\frac{\left(x - x'\right)^2}{2s^2}\right\}, \quad s > 0.}

The bandwidth parameter \eqn{s}, which controls the smoothness of the density estimate, is automatically determined by the algorithm if not provided by the user.
}


\examples{
# Generate simulated continuous data
set.seed(123)
X <- matrix(rnorm(200), ncol = 5)  # 200 observations, 5 features

# Run GIBcont with automatic bandwidth selection and multiple initializations
result <- GIBcont(X = X, ncl = 2, beta = 50, alpha = 0.75, s = -1, nstart = 20)

# Print clustering results
print(result$Cluster)       # Cluster membership matrix
print(result$Entropy)       # Entropy of final clustering
print(result$RelEntropy)    # Relative entropy of final clustering
print(result$MutualInfo)    # Mutual information between Y and T
}
\seealso{
\code{\link{GIBmix}}, \code{\link{GIBcat}}
}
\author{
Efthymios Costa, Ioanna Papatsouma, Angelos Markos
}
\references{
{
\insertRef{strouse_dib_2017}{IBclust}

\insertRef{silverman_density_1998}{IBclust}
}
}
\keyword{clustering}
