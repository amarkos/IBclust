\name{DIBcat}
\alias{DIBcat}
\title{Cluster Categorical Data Using the Deterministic Information Bottleneck Algorithm}
\description{
The \code{DIBcat} function implements the Deterministic Information Bottleneck (DIB) algorithm
for clustering datasets containing categorical variables. This method balances information retention
and data compression to create meaningful clusters, leveraging bandwidth parameters to handle
different categorical data types (nominal and ordinal) effectively (Costa, Papatsouma & Markos, 2024).
}
\usage{
DIBcat(X, ncl, randinit = NULL, lambda = -1,
       maxiter = 100, nstart = 100, select_features = FALSE)
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
    A numeric value, vector, or \code{-1} specifying the bandwidth parameters for categorical variables.
    If \code{-1}, the function automatically determines suitable values based on the number of levels in each variable.
  }
  \item{maxiter}{
    The maximum number of iterations for the clustering algorithm. Defaults to \code{100}.
  }
  \item{nstart}{
    The number of random initializations to run. The best clustering result (based on the information-theoretic criterion) is returned. Defaults to \code{100}.
  }
  \item{select_features}{
    Logical. If \code{TRUE}, uses an eigengap heuristic for feature selection, potentially improving clustering quality by reducing dimensionality. Defaults to \code{FALSE}.
  }
}
\details{
The \code{DIBcat} function applies the Deterministic Information Bottleneck algorithm to cluster datasets
containing categorical variables, including both nominal and ordinal data. The algorithm iteratively refines
cluster assignments to maximize mutual information between the original data distribution and the cluster
structure, while adaptively selecting or using provided bandwidth parameters (\code{lambda}) for categorical
variables. This method is particularly useful for handling complex categorical datasets and supports automatic
parameter selection to simplify its application.
}
\value{
A list containing the following elements:
  \itemize{
    \item{\code{Cluster}: An integer vector indicating the cluster assignment for each data point at convergence.}
    \item{\code{Entropy}: A numeric value representing the entropy of the cluster assignments at the end of the iterative procedure.}
    \item{\code{MutualInfo}: A numeric value representing the mutual information (\code{I(Y;T)}) between the data distribution and the cluster assignments.}
    \item{\code{lambda}: A numeric vector of bandwidth parameters for categorical variables, controlling how categories are weighted in the clustering.}
    \item{\code{ents}: A numeric vector tracking the entropy values across iterations, providing insights into the convergence pattern.}
    \item{\code{mis}: A numeric vector tracking the mutual information values across iterations.}
  }
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
Costa, E., Papatsouma, I., & Markos, A. (2024). A Deterministic Information Bottleneck Method for Clustering Mixed-Type Data. \emph{arXiv:2407.03389 [stat.ME]}. Retrieved from https://arxiv.org/abs/2407.03389
}
\keyword{clustering}
