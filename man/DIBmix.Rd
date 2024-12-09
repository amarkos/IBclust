\name{DIBmix}
\alias{DIBmix}
\title{Deterministic Information Bottleneck Clustering for Mixed-Type Data}
\description{
The \code{DIBmix} function implements the Deterministic Information Bottleneck (DIB) algorithm
for clustering datasets containing mixed-type variables, including categorical (nominal and ordinal)
and continuous variables. This method optimizes an information-theoretic objective to preserve
relevant information in the cluster assignments while achieving effective data compression
(Costa, Papatsouma & Markos, 2024).
}
\usage{
DIBmix(X, ncl, catcols, contcols, randinit = NULL,
       lambda = -1, s = -1, scale = TRUE,
       maxiter = 100, nstart = 100,
       select_features = FALSE)
}
\arguments{
\item{X}{A data frame containing the input data to be clustered. It should include categorical variables
(\code{factor} for nominal and \code{ordered} for ordinal) and continuous variables (\code{numeric}).}

\item{ncl}{An integer specifying the number of clusters.}

\item{catcols}{A vector indicating the indices of the categorical variables in \code{X}.}

\item{contcols}{A vector indicating the indices of the continuous variables in \code{X}.}

\item{randinit}{An optional vector specifying the initial cluster assignments. If \code{NULL}, cluster assignments are initialized randomly.}

\item{lambda}{A numeric value or vector specifying the bandwidth parameter for categorical variables. Defaults to \code{-1}, where the optimal bandwidth is automatically determined.}

\item{s}{A numeric value or vector specifying the bandwidth parameter(s) for continuous variables. Defaults to \code{-1}, where the optimal bandwidth(s) are selected automatically.}

\item{scale}{A logical value indicating whether the continuous variables should be scaled to have unit variance before clustering. Defaults to \code{TRUE}.}

\item{maxiter}{The maximum number of iterations allowed for the clustering algorithm. Defaults to \code{100}.}

\item{nstart}{The number of random initializations to run. The best clustering solution is returned. Defaults to \code{100}.}

\item{select_features}{A logical value indicating whether to perform feature selection based on the eigengap heuristic. Defaults to \code{FALSE}.}
}
\value{
A list with the following elements:
\item{Cluster}{An integer vector giving the cluster assignments for each data point.}
\item{Entropy}{A numeric value representing the entropy of the cluster assignments at convergence.}
\item{MutualInfo}{A numeric value representing the mutual information (\code{I(Y;T)}) between the original labels (\code{Y}) and the cluster assignments (\code{T}).}
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
}
\examples{
# Example dataset with categorical, ordinal, and continuous variables
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
Costa, E., Papatsouma, I., & Markos, A. (2024). A Deterministic Information Bottleneck Method for Clustering Mixed-Type Data. \emph{arXiv:2407.03389 [stat.ME]}. Retrieved from https://arxiv.org/abs/2407.03389
}
\keyword{clustering}
