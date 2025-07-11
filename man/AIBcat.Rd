\name{AIBcat}
\alias{AIBcat}
\title{Cluster Categorical Data Using the Agglomerative Information Bottleneck Algorithm}
\description{
The \code{AIBcat} function implements the Agglomerative Information Bottleneck (AIB) algorithm
for hierarchical clustering of datasets containing categorical variables. This method merges clusters
so that information retention is maximised at each step to create meaningful clusters,
leveraging bandwidth parameters to handle
different categorical data types (nominal and ordinal) effectively \insertCite{slonim_aib_1999}{IBclust}.
}
\usage{
AIBcat(X, lambda = -1, nomkernel = "aitchisonaitken", ordkernel = "liracine")
}
\arguments{
  \item{X}{
    A data frame containing the categorical data to be clustered. All variables should be categorical,
    either \code{factor} (for nominal variables) or \code{ordered} (for ordinal variables).
  }
  \item{lambda}{
    A numeric value or vector specifying the bandwidth parameter for categorical variables. The default value is \eqn{-1}, which enables automatic determination of the optimal bandwidth. For nominal variables, the maximum allowable value of \code{lambda} is \eqn{(l - 1)/l}, where \eqn{l} represents the number of categories. For ordinal variables, the maximum allowable value of \code{lambda} is \eqn{1}.
  }
  \item{nomkernel}{
    Kernel used for nominal (unordered categorical) variables. Can be one of \code{aitchisonaitken} (default) or \code{liracine}.
  }
  \item{ordkernel}{
    Kernel used for ordinal (ordered categorical) variables. Can be one of \code{liracine} (default) or \code{wangvanryzin}.
  }
}

\value{
A list containing the following elements:
  \item{merges}{A data frame with 2 columns and \eqn{n} rows, showing which observations are merged at each step.}
  \item{merge_costs}{A numeric vector tracking the cost incurred by each merge \eqn{I(Z_{m} ; Y) - I(Z_{m-1} ; Y)}.}
  \item{partitions}{A list containing \eqn{n} sub-lists. Each sub-list includes the cluster partition at each step.}
  \item{I_Z_Y}{A numeric vector including the mutual information \eqn{I(Z_{m}; Y)} as the number of clusters \eqn{m} increases.}
  \item{I_X_Y}{A numeric value of the mutual information \eqn{I(X; Y)} between observation indices and location.}
  \item{info_ret}{A numeric vector of length \eqn{n} including the fraction of the original information retained after each merge.}
  \item{dendrogram}{A dendrogram visualising the cluster hierarchy. The height is determined by the cost of cluster merges.}
}

\details{
The \code{AIBcat} function applies the Agglomerative Information Bottleneck algorithm to do hierarchical clustering of datasets containing only categorical variables, both nominal and ordinal. The algorithm uses an information-theoretic criterion to merge clusters so that information retention is maximised at each step to create meaningful clusters with maximal information about the original distribution.

To estimate the distributions of categorical features, the function utilizes specialized kernel functions.

For nominal (unordered categorical) variables, the kernel functions implemented are:

\itemize{
\item \emph{Aitchison & Aitken kernel \insertCite{aitchison_kernel_1976}{IBclust}:}
\deqn{K_u(x = x'; \lambda) = \begin{cases}
    1 - \lambda, & \text{if } x = x' \\
    \frac{\lambda}{\ell - 1}, & \text{otherwise}
\end{cases}, \quad 0 \leq \lambda \leq \frac{\ell - 1}{\ell}.}
\item \emph{Li & Racine kernel \insertCite{ouyang2006cross}{IBclust}:}
\deqn{K_u(x = x'; \lambda) = \begin{cases}
    1, & \text{if } x = x' \\
    \lambda, & \text{otherwise}
\end{cases}, \quad 0 \leq \lambda \leq 1.}
}

For ordinal (ordered categorical) variables, the kernel functions implemented are:

\itemize{
\item \emph{Li & Racine kernel \insertCite{li_nonparametric_2003}{IBclust}:}
\deqn{K_o(x = x'; \nu) = \begin{cases}
    1, & \text{if } x = x' \\
    \nu^{|x - x'|}, & \text{otherwise}
\end{cases}, \quad 0 \leq \nu \leq 1.}
\item \emph{Wang & van Ryzin kernel \insertCite{wang1981class}{IBclust}:}
\deqn{K_o(x = x'; \nu) = \begin{cases}
    1 - \nu, & \text{if } x = x' \\
    \frac{1-\nu}{2}\nu^{|x - x'|}, & \text{otherwise}
\end{cases}, \quad 0 \leq \nu \leq 1.}
}
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

# Run AIBcat with automatic lambda selection 
result <- AIBcat(X = X, lambda = -1)

# Print clustering results
plot(result$dendrogram, xlab = "", sub = "")  # Plot dendrogram
}
\seealso{
\code{\link{AIBmix}}, \code{\link{AIBcont}}
}
\author{
Efthymios Costa, Ioanna Papatsouma, Angelos Markos
}
\references{
{
\insertRef{slonim_aib_1999}{IBclust}

\insertRef{aitchison_kernel_1976}{IBclust}

\insertRef{li_nonparametric_2003}{IBclust}

\insertRef{ouyang2006cross}{IBclust}

\insertRef{wang1981class}{IBclust}
}
}
\keyword{clustering}
