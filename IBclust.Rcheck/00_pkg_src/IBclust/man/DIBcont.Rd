\name{DIBcont}
\alias{DIBcont}
\title{Cluster Continuous Data Using the Deterministic Information Bottleneck Algorithm}
\description{
The \code{DIBcont} function implements the Deterministic Information Bottleneck (DIB) algorithm
for clustering continuous data. This method optimizes an information-theoretic objective to
preserve relevant information while forming concise and interpretable cluster representations
\insertCite{costa_dib_2025}{IBclust}.
}
\usage{
DIBcont(X, ncl, randinit = NULL, s = -1, scale = TRUE,
        maxiter = 100, nstart = 100, verbose = FALSE)
}
\arguments{
  \item{X}{
    A numeric matrix or data frame containing the continuous data to be clustered. All variables should be of type \code{numeric}.
  }
  \item{ncl}{
    An integer specifying the number of clusters to form.
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
  \item{Cluster}{An integer vector indicating the cluster assignment for each observation.}
  \item{Entropy}{A numeric value representing the entropy of the cluster assignments at convergence.}
  \item{MutualInfo}{A numeric value representing the mutual information, \eqn{I(Y;T)}, between the underlying data distribution and the cluster assignments.}
  \item{beta}{A numeric vector of the final beta values used during the iterative optimization.}
  \item{s}{A numeric value or vector of bandwidth parameters used for the continuous variables. Typically, this will be a single value if all continuous variables share the same bandwidth.}
  \item{ents}{A numeric vector tracking the entropy values over the iterations, providing insight into the convergence process.}
  \item{mis}{A numeric vector tracking the mutual information values over the iterations.}
}

\details{
The \code{DIBcont} function applies the Deterministic Information Bottleneck algorithm to cluster datasets comprising only continuous variables. This method leverages an information-theoretic objective to optimize the trade-off between data compression and the preservation of relevant information about the underlying data distribution.

The function utilizes the Gaussian kernel \insertCite{silverman_density_1998}{IBclust} for estimating probability densities of continuous features. The kernel is defined as:

\deqn{K_c\left(\frac{x - x'}{s}\right) = \frac{1}{\sqrt{2\pi}} \exp\left\{-\frac{\left(x - x'\right)^2}{2s^2}\right\}, \quad s > 0.}

The bandwidth parameter \eqn{s}, which controls the smoothness of the density estimate, is automatically determined by the algorithm if not provided by the user.
}


\examples{
# Generate simulated continuous data
set.seed(123)
X <- matrix(rnorm(1000), ncol = 5)  # 200 observations, 5 features

# Run DIBcont with automatic bandwidth selection and multiple initializations
result <- DIBcont(X = X, ncl = 3, s = -1, nstart = 50)

# Print clustering results
print(result$Cluster)       # Cluster assignments
print(result$Entropy)       # Final entropy
print(result$MutualInfo)    # Mutual information
}
\seealso{
\code{\link{DIBmix}}, \code{\link{DIBcat}}
}
\author{
Efthymios Costa, Ioanna Papatsouma, Angelos Markos
}
\references{
{
\insertRef{costa_dib_2025}{IBclust}

\insertRef{silverman_density_1998}{IBclust}
}
}
\keyword{clustering}
