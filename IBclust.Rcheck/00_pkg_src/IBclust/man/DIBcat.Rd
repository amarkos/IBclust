\name{DIBcat}
\alias{DIBcat}
\title{Cluster Categorical Data Using the Deterministic Information Bottleneck Algorithm}
\description{
The \code{DIBcat} function implements the Deterministic Information Bottleneck (DIB) algorithm
for clustering datasets containing categorical variables. This method balances information retention
and data compression to create meaningful clusters, leveraging bandwidth parameters to handle
different categorical data types (nominal and ordinal) effectively \insertCite{costa_dib_2025}{IBclust}.
}
\usage{
DIBcat(X, ncl, randinit = NULL, lambda = -1,
       maxiter = 100, nstart = 100,
       verbose = FALSE)
}
\arguments{
  \item{X}{
    A data frame containing the categorical data to be clustered. All variables should be categorical,
    either \code{factor} (for nominal variables) or \code{ordered} (for ordinal variables).
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
  \item{Cluster}{An integer vector indicating the cluster assignment for each data point at convergence.}
  \item{Entropy}{A numeric value representing the entropy of the cluster assignments at the end of the iterative procedure.}
  \item{MutualInfo}{A numeric value representing the mutual information, \eqn{I(Y;T)}, between the data distribution and the cluster assignments.}
  \item{lambda}{A numeric vector of bandwidth parameters for categorical variables, controlling how categories are weighted in the clustering.}
  \item{beta}{A numeric vector of the final beta values used during the iterative optimization.}
  \item{ents}{A numeric vector tracking the entropy values across iterations, providing insights into the convergence pattern.}
  \item{mis}{A numeric vector tracking the mutual information values across iterations.}
}

\details{
The \code{DIBcat} function applies the Deterministic Information Bottleneck algorithm to cluster datasets containing only categorical variables, both nominal and ordinal. The algorithm optimizes an information-theoretic objective to balance the trade-off between data compression and the retention of information about the original distribution.

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

# Run DIBcat with automatic lambda selection and multiple initializations
result <- DIBcat(X = X, ncl = 3, lambda = -1, nstart = 50)

# Print clustering results
print(result$Cluster)       # Cluster assignments
print(result$Entropy)       # Final entropy
print(result$MutualInfo)    # Mutual information
}
\seealso{
\code{\link{DIBmix}}, \code{\link{DIBcont}}
}
\author{
Efthymios Costa, Ioanna Papatsouma, Angelos Markos
}
\references{
{
\insertRef{costa_dib_2025}{IBclust}

\insertRef{aitchison_kernel_1976}{IBclust}

\insertRef{li_nonparametric_2003}{IBclust}
}
}
\keyword{clustering}
