\name{IBcont}
\alias{IBcont}
\title{Cluster Continuous Data Using the Information Bottleneck Algorithm}
\description{
The \code{IBcont} function implements the Information Bottleneck (IB) algorithm
for fuzzy clustering of continuous data. This method optimizes an information-theoretic objective to
preserve relevant information while forming concise and interpretable cluster representations
\insertCite{strouse_ib_2019}{IBclust}.
}
\usage{
IBcont(X, ncl, beta, randinit = NULL, s = -1, scale = TRUE,
       maxiter = 100, nstart = 100, contkernel = "gaussian",
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
  \item{contkernel}{
    Kernel used for continuous variables. Can be one of \code{gaussian} (default) or \code{epanechnikov}.
  }
  \item{verbose}{
    Logical. Default to \code{FALSE} to suppress progress messages. Change to \code{TRUE} to print.
  }
}

\value{
A list containing the following elements:
  \item{Cluster}{A cluster membership matrix.}
  \item{Entropy}{A numeric value representing the entropy of the cluster assignments at the end of the iterative procedure.}
  \item{CondEntropy}{A numeric value representing the conditional entropy of cluster assignment, given the observation weights \eqn{H(T \mid X)}.}
  \item{MutualInfo}{A numeric value representing the mutual information, \eqn{I(Y;T)}, between the data distribution and the cluster assignments.}
  \item{InfoXT}{A numeric value representing the mutual information, \eqn{I(X;T)}, between the original observations weights (\eqn{X}) and the cluster assignments (\eqn{T}).}
  \item{beta}{A numeric value of the regularisation strength beta used.}
  \item{s}{A numeric vector of bandwidth parameters used for the continuous variables.}
  \item{ixt}{A numeric vector tracking the mutual information values between original observation weights and cluster assignments across iterations.}
  \item{iyt}{A numeric vector tracking the mutual information values between original labels and cluster assignments across iterations.}
  \item{losses}{A numeric vector tracking the final loss values across iterations.}
}

\details{
The \code{IBcont} function applies the Information Bottleneck algorithm to do fuzzy clustering of datasets comprising only continuous variables. This method leverages an information-theoretic objective to optimize the trade-off between data compression and the preservation of relevant information about the underlying data distribution.

To estimate the distributions of continuous features, the function utilizes specialized kernel functions:

\itemize{
  \item \emph{Gaussian (RBF) kernel \insertCite{silverman_density_1998}{IBclust}:}
  \deqn{K_c\left(\frac{x - x'}{s}\right) = \frac{1}{\sqrt{2\pi}} \exp\left\{-\frac{\left(x - x'\right)^2}{2s^2}\right\}, \quad s > 0.}
  \item \emph{Epanechnikov kernel \insertCite{epanechnikov1969non}{IBclust}:}
  \deqn{K_c(x - x'; s) = \begin{cases}
    \frac{3}{4\sqrt{5}}\left(1 - \frac{(x-x')^2}{5s^2} \right), & \text{if } \frac{(x - x')^2}{s^2} < 5 \\
    0, & \text{otherwise}
\end{cases}, \quad s > 0.}
}
}


\examples{
# Generate simulated continuous data
set.seed(123)
X <- matrix(rnorm(200), ncol = 5)  # 200 observations, 5 features

# Run IBcont with automatic bandwidth selection and multiple initializations
result <- IBcont(X = X, ncl = 3, beta = 50, s = -1, nstart = 20)

# Print clustering results
print(result$Cluster)       # Cluster membership matrix
print(result$InfoXT)       # Mutual information between X and T
print(result$MutualInfo)    # Mutual information between Y and T
}
\seealso{
\code{\link{IBmix}}, \code{\link{IBcat}}
}
\author{
Efthymios Costa, Ioanna Papatsouma, Angelos Markos
}
\references{
{
\insertRef{strouse_ib_2019}{IBclust}

\insertRef{silverman_density_1998}{IBclust}

\insertRef{epanechnikov1969non}{IBclust}
}
}
\keyword{clustering}
