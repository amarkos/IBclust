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
       contkernel = "gaussian",
       nomkernel = "aitchisonaitken", ordkernel = "liracine",
       cat_first = FALSE, verbose = FALSE)
}
\arguments{
\item{X}{A data frame containing the input data to be clustered. It should include categorical variables
(\code{factor} for nominal and \code{Ord.factor} for ordinal) and continuous variables (\code{numeric}).}

\item{ncl}{An integer specifying the number of clusters.}

\item{catcols}{A vector indicating the indices of the categorical variables in \code{X}.}

\item{contcols}{A vector indicating the indices of the continuous variables in \code{X}.}

\item{randinit}{An optional vector specifying the initial cluster assignments. If \code{NULL}, cluster assignments are initialized randomly.}

\item{lambda}{A numeric value or vector specifying the bandwidth parameter for categorical variables. The default value is \eqn{-1}, which enables automatic determination of the optimal bandwidth. For nominal variables and \code{nomkernel = 'aitchisonaitken'}, the maximum allowable value of \code{lambda} is \eqn{(l - 1)/l}, where \eqn{l} represents the number of categories, whereas for \code{nomkernel = 'liracine'} the maximum allowable value is \eqn{1}. For ordinal variables, the maximum allowable value of \code{lambda} is \eqn{1}, regardless of what \code{ordkernel} is being used.}

\item{s}{
A numeric value or vector specifying the bandwidth parameter(s) for continuous variables. The values must be greater than \eqn{0}. The default value is \eqn{-1}, which enables the automatic selection of optimal bandwidth(s).
}

\item{scale}{A logical value indicating whether the continuous variables should be scaled to have unit variance before clustering. Defaults to \code{TRUE}.}

\item{maxiter}{The maximum number of iterations allowed for the clustering algorithm. Defaults to \eqn{100}.}

\item{nstart}{The number of random initializations to run. The best clustering solution is returned. Defaults to \eqn{100}.}

\item{contkernel}{Kernel used for continuous variables. Can be one of \code{gaussian} (default) or \code{epanechnikov}.}

\item{nomkernel}{Kernel used for nominal (unordered categorical) variables. Can be one of \code{aitchisonaitken} (default) or \code{liracine}.}

\item{ordkernel}{Kernel used for ordinal (ordered categorical) variables. Can be one of \code{liracine} (default) or \code{wangvanryzin}.}

\item{cat_first}{A logical value indicating whether bandwidth selection is prioritised for the categorical variables, instead of the continuous. Defaults to \code{FALSE}. Set to \code{TRUE} if you suspect that the continuous variables are not informative of the cluster structure. Can only be \code{TRUE} when all bandwidths are selected automatically (i.e. \code{s = -1}, \code{lambda = -1}).}

\item{verbose}{Logical. Defaults to \code{FALSE} to suppress progress messages. Change to \code{TRUE} to print.}
}
\value{
A list containing the following elements:
  \item{Cluster}{An integer vector giving the cluster assignments for each data point.}
  \item{Entropy}{A numeric value representing the entropy of the cluster assignments at convergence.}
  \item{CondEntropy}{A numeric value representing the conditional entropy of cluster assignment, given the observation weights \eqn{H(T \mid X)}.}
  \item{MutualInfo}{A numeric value representing the mutual information, \eqn{I(Y;T)}, between the original labels (\eqn{Y}) and the cluster assignments (\eqn{T}).}
  \item{InfoXT}{A numeric value representing the mutual information, \eqn{I(X;T)}, between the original observations weights (\eqn{X}) and the cluster assignments (\eqn{T}).}
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
# Example 1: Basic Mixed-Type Clustering
set.seed(123)

# Create a more realistic dataset with mixed variable types
data <- data.frame(
  # Categorical variables
  education = factor(sample(c("High School", "Bachelor", "Master", "PhD"), 150, 
                           replace = TRUE, prob = c(0.4, 0.3, 0.2, 0.1))),
  employment = factor(sample(c("Full-time", "Part-time", "Unemployed"), 150, 
                            replace = TRUE, prob = c(0.6, 0.25, 0.15))),
  
  # Ordinal variable
  satisfaction = factor(sample(c("Low", "Medium", "High"), 150, replace = TRUE),
                       levels = c("Low", "Medium", "High"), ordered = TRUE),
  
  # Continuous variables  
  income = rlnorm(150, meanlog = 10, sdlog = 0.5),  # Log-normal income
  age = rnorm(150, mean = 35, sd = 10),             # Normally distributed age
  experience = rpois(150, lambda = 8)               # Years of experience
)

# Perform Mixed-Type Clustering
result <- DIBmix(X = data, ncl = 3, 
                catcols = c(1, 2, 3),    # education, employment, satisfaction
                contcols = c(4, 5, 6),   # income, age, experience
                nstart = 50)

# View results
print(paste("Number of clusters found:", length(unique(result$Cluster))))
print(paste("Mutual Information:", round(result$MutualInfo, 3)))
table(result$Cluster)

# Example 2: Comparing cat_first parameter
# When categorical variables are more informative
result_cat_first <- DIBmix(X = data, ncl = 3,
                          catcols = c(1, 2, 3), 
                          contcols = c(4, 5, 6),
                          cat_first = TRUE,  # Prioritize categorical variables
                          nstart = 50)

# When continuous variables are more informative (default)
result_cont_first <- DIBmix(X = data, ncl = 3,
                           catcols = c(1, 2, 3), 
                           contcols = c(4, 5, 6),
                           cat_first = FALSE,
                           nstart = 50)

# Compare clustering performance
library(mclust)  # For adjustedRandIndex
print(paste("Agreement between approaches:", 
           round(adjustedRandIndex(result_cat_first$Cluster, 
                                 result_cont_first$Cluster), 3)))                        
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

\insertRef{ouyang2006cross}{IBclust}

\insertRef{wang1981class}{IBclust}

\insertRef{epanechnikov1969non}{IBclust}
}
}
\keyword{clustering}