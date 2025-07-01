\name{GIBcat}
\alias{GIBcat}
\title{Cluster Categorical Data Using the Generalised Information Bottleneck Algorithm}
\description{
The \code{GIBcat} function implements the Generalised Information Bottleneck (GIB) algorithm
for fuzzy clustering of datasets containing categorical variables. This method balances information retention
and data compression to create meaningful clusters, leveraging bandwidth parameters to handle
different categorical data types (nominal and ordinal) effectively \insertCite{strouse_ib_2019}{IBclust}.
}
\usage{
GIBcat(X, ncl, beta, alpha, randinit = NULL, lambda = -1,
       maxiter = 100, nstart = 100, verbose = FALSE)
}
\arguments{
  \item{X}{
    A data frame containing the categorical data to be clustered. All variables should be categorical,
    either \code{factor} (for nominal variables) or \code{ordered} (for ordinal variables).
  }
  \item{beta}{
    Regularisation strength.
  }
  \item{alpha}{
    Strength of relative entropy term.
  }
  \item{ncl}{
    An integer specifying the number of clusters to form.
  }
  \item{randinit}{
    Optional. A vector specifying initial cluster assignments. If \code{NULL}, cluster assignments are initialized randomly.
  }
\item{lambda}{
A numeric value or vector specifying the bandwidth parameter for categorical variables. The default value is \eqn{-1}, which enables automatic determination of the optimal bandwidth. For nominal variables, the maximum allowable value of \code{lambda} is \eqn{(l - 1)/l}, where \eqn{l} represents the number of categories. For ordinal variables, the maximum allowable value of \code{lambda} is \eqn{1}.
}

  \item{maxiter}{
    The maximum number of iterations for the clustering algorithm. Defaults to \eqn{100}.
  }
  \item{nstart}{
    The number of random initializations to run. The best clustering result (based on the information-theoretic criterion) is returned. Defaults to \eqn{100}.
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
  \item{lambda}{A numeric vector of bandwidth parameters for categorical variables, controlling how categories are weighted in the clustering.}
  \item{ht}{A numeric vector tracking the entropy value of the cluster assignments across iterations.}
  \item{hy_t}{A numeric vector tracking the relative entropy values between the cluster assignments and observations weights across iterations.}
  \item{iyt}{A numeric vector tracking the mutual information values between original labels and cluster assignments across iterations.}
  \item{losses}{A numeric vector tracking the final loss values across iterations.}
}

\details{
The \code{GIBcat} function applies the Generalised Information Bottleneck algorithm to do fuzzy clustering of datasets containing only categorical variables, both nominal and ordinal. The algorithm optimizes an information-theoretic objective to balance the trade-off between data compression and the retention of information about the original distribution.
Set \eqn{\alpha = 1} and \eqn{\alpha = 0} to recover the Information Bottleneck and its Deterministic variant, respectively. If \eqn{\alpha = 0}, the algorithm ignores
the value of the regularisation parameter \eqn{\beta}.

To estimate the distributions of categorical features, the function utilizes specialized kernel functions, as follows:

\deqn{K_u(x = x'; \lambda) = \begin{cases}
    1 - \lambda, & \text{if } x = x' \\
    \frac{\lambda}{\ell - 1}, & \text{otherwise}
\end{cases}, \quad 0 \leq \lambda \leq \frac{\ell - 1}{\ell},}
where \eqn{\ell} is the number of categories, and \eqn{\lambda} controls the smoothness of the Aitchison & Aitken kernel for nominal variables \insertCite{aitchison_kernel_1976}{IBclust}.

\deqn{K_o(x = x'; \nu) = \begin{cases}
    1, & \text{if } x = x' \\
    \nu^{|x - x'|}, & \text{otherwise}
\end{cases}, \quad 0 \leq \nu \leq 1,}
where \eqn{\nu} is the bandwidth parameter for ordinal variables, accounting for the ordinal relationship between categories \insertCite{li_nonparametric_2003}{IBclust}.

Here, \eqn{\lambda}, and \eqn{\nu} are bandwidth or smoothing parameters, while \eqn{\ell} is the number of levels of the categorical variable. The lambda parameter is automatically determined by the algorithm if not provided by the user. For ordinal variables, the lambda parameter of the function is used to define \eqn{\nu}.
}


\examples{
# Simulated categorical data
set.seed(123)
X <- data.frame(
  Var1 = as.factor(sample(letters[1:3], 200, replace = TRUE)),  # Nominal variable
  Var2 = as.factor(sample(letters[4:6], 200, replace = TRUE)),  # Nominal variable
  Var3 = factor(sample(c("low", "medium", "high"), 200, replace = TRUE),
                levels = c("low", "medium", "high"), ordered = TRUE)  # Ordinal variable
)

# Run GIBcat with automatic lambda selection and multiple initializations
result <- GIBcat(X = X, ncl = 2, beta = 25, alpha = 0.75, lambda = -1, nstart = 20)

# Print clustering results
print(result$Cluster)       # Cluster membership matrix
print(result$Entropy)       # Entropy of final clustering
print(result$RelEntropy)    # Relative entropy of final clustering
print(result$MutualInfo)    # Mutual information between Y and T
}
\seealso{
\code{\link{GIBmix}}, \code{\link{GIBcont}}
}
\author{
Efthymios Costa, Ioanna Papatsouma, Angelos Markos
}
\references{
{
\insertRef{strouse_dib_2017}{IBclust}

\insertRef{aitchison_kernel_1976}{IBclust}

\insertRef{li_nonparametric_2003}{IBclust}
}
}
\keyword{clustering}
