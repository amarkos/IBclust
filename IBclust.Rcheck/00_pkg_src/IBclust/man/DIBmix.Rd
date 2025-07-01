\name{DIBmix}
\alias{DIBmix}
\title{Deterministic Information Bottleneck Clustering for Mixed-Type Data}
\description{
The \code{DIBmix} function implements the Deterministic Information Bottleneck (DIB) algorithm
for clustering datasets containing mixed-type variables, including categorical (nominal and ordinal)
and continuous variables. This method optimizes an information-theoretic objective to preserve
relevant information in the cluster assignments while achieving effective data compression
\insertCite{costa_dib_2025}{IBclust}.
}
\usage{
DIBmix(X, ncl, catcols, contcols, randinit = NULL,
       lambda = -1, s = -1, scale = TRUE,
       maxiter = 100, nstart = 100,
       verbose = FALSE)
}
\arguments{
\item{X}{A data frame containing the input data to be clustered. It should include categorical variables
(\code{factor} for nominal and \code{Ord.factor} for ordinal) and continuous variables (\code{numeric}).}

\item{ncl}{An integer specifying the number of clusters.}

\item{catcols}{A vector indicating the indices of the categorical variables in \code{X}.}

\item{contcols}{A vector indicating the indices of the continuous variables in \code{X}.}

\item{randinit}{An optional vector specifying the initial cluster assignments. If \code{NULL}, cluster assignments are initialized randomly.}

\item{lambda}{
A numeric value or vector specifying the bandwidth parameter for categorical variables. The default value is \eqn{-1}, which enables automatic determination of the optimal bandwidth. For nominal variables, the maximum allowable value of \code{lambda} is \eqn{(l - 1)/l}, where \eqn{l} represents the number of categories. For ordinal variables, the maximum allowable value of \code{lambda} is 1.
}

\item{s}{
A numeric value or vector specifying the bandwidth parameter(s) for continuous variables. The values must be greater than \eqn{0}. The default value is \eqn{-1}, which enables the automatic selection of optimal bandwidth(s).
}

\item{scale}{A logical value indicating whether the continuous variables should be scaled to have unit variance before clustering. Defaults to \code{TRUE}.}

\item{maxiter}{The maximum number of iterations allowed for the clustering algorithm. Defaults to \eqn{100}.}

\item{nstart}{The number of random initializations to run. The best clustering solution is returned. Defaults to \eqn{100}.}

\item{verbose}{Logical. Default to \code{FALSE} to suppress progress messages. Change to \code{TRUE} to print.}
}
\value{
A list containing the following elements:
  \item{Cluster}{An integer vector giving the cluster assignments for each data point.}
  \item{Entropy}{A numeric value representing the entropy of the cluster assignments at convergence.}
  \item{MutualInfo}{A numeric value representing the mutual information, \eqn{I(Y;T)}, between the original labels (\eqn{Y}) and the cluster assignments (\eqn{T}).}
  \item{beta}{A numeric vector of the final beta values used in the iterative procedure.}
  \item{s}{A numeric vector of bandwidth parameters used for the continuous variables.}
  \item{lambda}{A numeric vector of bandwidth parameters used for the categorical variables.}
  \item{ents}{A numeric vector tracking the entropy values across iterations.}
  \item{mis}{A numeric vector tracking the mutual information values across iterations.}
}


\details{
The \code{DIBmix} function clusters data while retaining maximal information about the original variable
distributions. The Deterministic Information Bottleneck algorithm optimizes an information-theoretic
objective that balances information preservation and compression. Bandwidth parameters for categorical
(nominal, ordinal) and continuous variables are adaptively determined if not provided. This iterative
process identifies stable and interpretable cluster assignments by maximizing mutual information while
controlling complexity. The method is well-suited for datasets with mixed-type variables and integrates
information from all variable types effectively.

The following kernel functions are used to estimate densities for the clustering procedure:
\itemize{
    \item \emph{Continuous variables: Gaussian kernel}
    \deqn{K_c\left(\frac{x-x'}{s}\right) = \frac{1}{\sqrt{2\pi}} \exp\left\{ - \frac{\left(x-x'\right)^2}{2s^2} \right\}, \quad s > 0.}
    \item \emph{Nominal categorical variables: Aitchison & Aitken kernel}
    \deqn{K_u\left(x = x' ; \lambda\right) = \begin{cases}
          1-\lambda & \text{if } x = x' \\
          \frac{\lambda}{\ell-1} & \text{otherwise}
        \end{cases}, \quad 0 \leq \lambda \leq \frac{\ell-1}{\ell}.}
    \item \emph{Ordinal categorical variables: Li & Racine kernel}
    \deqn{K_o\left(x = x' ; \nu\right) = \begin{cases}
          1 & \text{if } x = x' \\
          \nu^{|x - x'|} & \text{otherwise}
        \end{cases}, \quad 0 \leq \nu \leq 1.}
}

Here, \eqn{s}, \eqn{\lambda}, and \eqn{\nu} are bandwidth or smoothing parameters, while \eqn{\ell} is the number of levels of the categorical variable. \eqn{s} and \eqn{\lambda} are automatically determined by the algorithm if not provided by the user. For ordinal variables, the lambda parameter of the function is used to define \eqn{\nu}.
}


\examples{
# Example dataset with categorical, ordinal, and continuous variables
set.seed(123)
data <- data.frame(
  cat_var = factor(sample(letters[1:3], 100, replace = TRUE)),      # Nominal categorical variable
  ord_var = factor(sample(c("low", "medium", "high"), 100, replace = TRUE),
                   levels = c("low", "medium", "high"),
                   ordered = TRUE),                                # Ordinal variable
  cont_var1 = rnorm(100),                                          # Continuous variable 1
  cont_var2 = runif(100)                                           # Continuous variable 2
)

# Perform Mixed-Type Clustering
result <- DIBmix(X = data, ncl = 3, catcols = 1:2, contcols = 3:4)

# Print clustering results
print(result$Cluster)       # Cluster assignments
print(result$Entropy)       # Final entropy
print(result$MutualInfo)    # Mutual information
}
\seealso{
\code{\link{DIBcont}}, \code{\link{DIBcat}}
}
\author{
Efthymios Costa, Ioanna Papatsouma, Angelos Markos
}
\references{
{
\insertRef{costa_dib_2025}{IBclust}

\insertRef{aitchison_kernel_1976}{IBclust}

\insertRef{li_nonparametric_2003}{IBclust}

\insertRef{silverman_density_1998}{IBclust}
}
}
\keyword{clustering}
