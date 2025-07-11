\name{IBmix}
\alias{IBmix}
\title{Information Bottleneck Clustering for Mixed-Type Data}
\description{
The \code{IBmix} function implements the Information Bottleneck (IB) algorithm
for clustering datasets containing mixed-type variables, including categorical (nominal and ordinal)
and continuous variables. This method optimizes an information-theoretic objective to preserve
relevant information in the cluster assignments while achieving effective data compression
\insertCite{strouse_ib_2019}{IBclust}.
}
\usage{
IBmix(X, ncl, beta, catcols, contcols, randinit = NULL,
       lambda = -1, s = -1, scale = TRUE,
       maxiter = 100, nstart = 100,
       contkernel = "gaussian",
       nomkernel = "aitchisonaitken", ordkernel = "liracine",
       verbose = FALSE)
}
\arguments{
\item{X}{A data frame containing the input data to be clustered. It should include categorical variables
(\code{factor} for nominal and \code{Ord.factor} for ordinal) and continuous variables (\code{numeric}).}

\item{ncl}{An integer specifying the number of clusters.}

\item{beta}{Regularisation strength.}

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

\item{contkernel}{Kernel used for continuous variables. Can be one of \code{gaussian} (default) or \code{epanechnikov}.}

\item{nomkernel}{Kernel used for nominal (unordered categorical) variables. Can be one of \code{aitchisonaitken} (default) or \code{liracine}.}

\item{ordkernel}{Kernel used for ordinal (ordered categorical) variables. Can be one of \code{liracine} (default) or \code{wangvanryzin}.}

\item{verbose}{Logical. Defaults to \code{FALSE} to suppress progress messages. Change to \code{TRUE} to print.}
}

\value{
A list containing the following elements:
  \item{Cluster}{A cluster membership matrix.}
  \item{InfoXT}{A numeric value representing the mutual information, \eqn{I(X;T)}, between the original observations weights (\eqn{X}) and the cluster assignments (\eqn{T}).}
  \item{InfoYT}{A numeric value representing the mutual information, \eqn{I(Y;T)}, between the original labels (\eqn{Y}) and the cluster assignments (\eqn{T}).}
  \item{beta}{A numeric value of the regularisation strength beta used.}
  \item{s}{A numeric vector of bandwidth parameters used for the continuous variables.}
  \item{lambda}{A numeric vector of bandwidth parameters used for the categorical variables.}
  \item{ixt}{A numeric vector tracking the mutual information values between original observation weights and cluster assignments across iterations.}
  \item{iyt}{A numeric vector tracking the mutual information values between original labels and cluster assignments across iterations.}
  \item{losses}{A numeric vector tracking the final loss values across iterations.}
}


\details{
The \code{IBmix} function produces a fuzzy clustering of the data while retaining maximal information about the original variable
distributions. The Information Bottleneck algorithm optimizes an information-theoretic
objective that balances information preservation and compression. Bandwidth parameters for categorical
(nominal, ordinal) and continuous variables are adaptively determined if not provided. This iterative
process identifies stable and interpretable cluster assignments by maximizing mutual information while
controlling complexity. The method is well-suited for datasets with mixed-type variables and integrates
information from all variable types effectively.

The following kernel functions can be used to estimate densities for the clustering procedure. For continuous variables:

\itemize{
  \item \emph{Gaussian (RBF) kernel \insertCite{silverman_density_1998}{IBclust}:}
  \deqn{K_c\left(\frac{x - x'}{s}\right) = \frac{1}{\sqrt{2\pi}} \exp\left\{-\frac{\left(x - x'\right)^2}{2s^2}\right\}, \quad s > 0.}
  \item \emph{Epanechnikov kernel \insertCite{epanechnikov1969non}{IBclust}:}
  \deqn{K_c(x - x'; s) = \begin{cases}
    \frac{3}{4\sqrt{5}}\left(1 - \frac{(x-x')^2}{5s^2} \right), & \text{if } \frac{(x - x')^2}{s^2} < 5 \\
    0, & \text{otherwise}
\end{cases}, \quad s > 0.}
}

For nominal (unordered categorical variables):

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

For ordinal (ordered categorical) variables:

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

The bandwidth parameters \eqn{s}, \eqn{\lambda}, and \eqn{\nu} control the smoothness of the density estimate and are automatically determined by the algorithm if not provided by the user. \eqn{\ell} is the number of levels of the categorical variable. For ordinal variables, the lambda parameter of the function is used to define \eqn{\nu}.
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

# Perform Mixed-Type Fuzzy Clustering
result <- IBmix(X = data, ncl = 3, beta = 2, catcols = 1:2, contcols = 3:4, nstart = 20)

# Print clustering results
print(result$Cluster)       # Cluster membership matrix
print(result$InfoXT)       # Mutual information between X and T
print(result$InfoYT)    # Mutual information between Y and T
}
\seealso{
\code{\link{DIBcont}}, \code{\link{DIBcat}}
}
\author{
Efthymios Costa, Ioanna Papatsouma, Angelos Markos
}
\references{
{
\insertRef{strouse_ib_2019}{IBclust}

\insertRef{aitchison_kernel_1976}{IBclust}

\insertRef{li_nonparametric_2003}{IBclust}

\insertRef{silverman_density_1998}{IBclust}

\insertRef{ouyang2006cross}{IBclust}

\insertRef{wang1981class}{IBclust}

\insertRef{epanechnikov1969non}{IBclust}
}
}
\keyword{clustering}
