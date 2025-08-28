\name{AIBmix}
\alias{AIBmix}
\title{Agglomerative Information Bottleneck Clustering for Mixed-Type Data}
\description{
The \code{AIBmix} function implements the Agglomerative Information Bottleneck (AIB) algorithm
for hierarchical clustering of datasets containing mixed-type variables, including categorical (nominal and ordinal)
and continuous variables. This method merges clusters so that information retention is maximised at each step to create meaningful clusters,
leveraging bandwidth parameters to handle different categorical data types (nominal and ordinal) effectively \insertCite{slonim_aib_1999}{IBclust}.
}
\usage{
AIBmix(X, s = -1, lambda = -1,
       scale = TRUE, contkernel = "gaussian",
       nomkernel = "aitchisonaitken", ordkernel = "liracine",
       cat_first = FALSE)
}
\arguments{
  \item{X}{
    A data frame containing the data to be clustered. Variables should be of type \code{numeric} (for continuous variables),
    \code{factor} (for nominal variables) or \code{ordered} (for ordinal variables).
  }
  \item{s}{
  A numeric value or vector specifying the bandwidth parameter(s) for continuous variables. The values must be greater than \eqn{0}. The default value is \eqn{-1}, which enables the automatic selection of optimal bandwidth(s). Argument is ignored when no variables are continuous.
  }
  \item{lambda}{
    A numeric value or vector specifying the bandwidth parameter for categorical variables. The default value is \eqn{-1}, which enables automatic determination of the optimal bandwidth. For nominal variables and \code{nomkernel = 'aitchisonaitken'}, the maximum allowable value of \code{lambda} is \eqn{(l - 1)/l}, where \eqn{l} represents the number of categories, whereas for \code{nomkernel = 'liracine'} the maximum allowable value is \eqn{1}. For ordinal variables, the maximum allowable value of \code{lambda} is \eqn{1}, regardless of what \code{ordkernel} is being used. Argument is ignored when all variables are continuous.
  }
  \item{scale}{
    A logical value indicating whether the continuous variables should be scaled to have unit variance before clustering. Defaults to \code{TRUE}. Argument is ignored when all variables are categorical.
  }
  \item{contkernel}{
    Kernel used for continuous variables. Can be one of \code{gaussian} (default) or \code{epanechnikov}. Argument is ignored when no variables are continuous.
  }
  \item{nomkernel}{
    Kernel used for nominal (unordered categorical) variables. Can be one of \code{aitchisonaitken} (default) or \code{liracine}. Argument is ignored when no variables are nominal.
  }
  \item{ordkernel}{
    Kernel used for ordinal (ordered categorical) variables. Can be one of \code{liracine} (default) or \code{wangvanryzin}. Argument is ignored when no variables are ordinal.
  }
  \item{cat_first}{
    A logical value indicating whether bandwidth selection is prioritised for the categorical variables, instead of the continuous. Defaults to \code{FALSE}. Set to \code{TRUE} if you suspect that the continuous variables are not informative of the cluster structure. Can only be \code{TRUE} when all bandwidths are selected automatically (i.e. \code{s = -1}, \code{lambda = -1}). 
  }
}

\value{
An object of class \code{"aibclust"} representing the final clustering result. The returned object is a list with the following components:
  \item{merges}{A data frame with 2 columns and \eqn{n} rows, showing which observations are merged at each step.}
  \item{merge_costs}{A numeric vector tracking the cost incurred by each merge \eqn{I(T_{m} ; Y) - I(T_{m-1} ; Y)}.}
  \item{partitions}{A list containing \eqn{n} sub-lists. Each sub-list includes the cluster partition at each step.}
  \item{I_T_Y}{A numeric vector including the mutual information \eqn{I(T_{m}; Y)} as the number of clusters \eqn{m} increases.}
  \item{I_X_Y}{A numeric value of the mutual information \eqn{I(X; Y)} between observation indices and location.}
  \item{info_ret}{A numeric vector of length \eqn{n} including the fraction of the original information retained after each merge.}
  \item{s}{A numeric vector of bandwidth parameters used for the continuous variables. A value of \eqn{-1} is returned if all variables are categorical.}
  \item{lambda}{A numeric vector of bandwidth parameters used for the categorical variables. A value of \eqn{-1} is returned if all variables are continuous.}
  \item{call}{The matched call.}
  \item{n}{Number of observations.}
  \item{contcols}{Indices of continuous columns in \code{X}.}
  \item{catcols}{Indices of categorical columns in \code{X}.}
  \item{kernels}{List with names of kernels used for continuous, nominal, and ordinal features.}
  \item{obs_names}{Names of rows in \code{X}; used for plotting the cluster hierarchy using a dendrogram.}

Objects of class \code{"aibclust"} support the following methods:
  \itemize{
    \item \code{print.aibclust}: Display a concise description of the cluster hierarchy.
    \item \code{summary.aibclust}: Show detailed information including cluster sizes for 2, 3, and 5 clusters, 
          information-theoretic metrics, and hyperparameters.
    \item \code{plot.aibclust}: Produce diagnostic plots:
      \itemize{
        \item \code{type = "dendrogram"}: dendrogram visualising the hierarchy of partitions obtained.
        \item \code{type = "info"}: information retention curve; the proportion of information preserved \eqn{I(T_m;Y)/I(X;Y)} by the clustering \eqn{T_m} is plotted against the number of clusters \eqn{m}.
    }
  }

}

\details{
The \code{AIBmix} function produces a hierarchical agglomerative clustering of the data while retaining maximal information about the original variable
distributions. The Agglomerative Information Bottleneck algorithm uses an information-theoretic criterion to merge clusters so that information retention is maximised at each step,
hence creating meaningful clusters with maximal information about the original distribution. Bandwidth parameters for categorical
(nominal, ordinal) and continuous variables are adaptively determined if not provided. This process identifies stable and interpretable cluster assignments by maximizing mutual information while
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

The bandwidth parameters \eqn{s}, \eqn{\lambda}, and \eqn{\nu} control the smoothness of the density estimate and are automatically determined by the algorithm if not provided by the user using the approach in \insertCite{costa_dib_2025;textual}{IBclust}. \eqn{\ell} is the number of levels of the categorical variable. For ordinal variables, the lambda parameter of the function is used to define \eqn{\nu}.
}



\examples{
# Example dataset with categorical, ordinal, and continuous variables
set.seed(123)
data_mix <- data.frame(
  cat_var = factor(sample(letters[1:3], 100, replace = TRUE)),      # Nominal categorical variable
  ord_var = factor(sample(c("low", "medium", "high"), 100, replace = TRUE),
                   levels = c("low", "medium", "high"),
                   ordered = TRUE),                                # Ordinal variable
  cont_var1 = rnorm(100),                                          # Continuous variable 1
  cont_var2 = runif(100)                                           # Continuous variable 2
)

# Perform Mixed-Type Hierarchical Clustering with Agglomerative IB
result_mix <- AIBmix(X = data_mix, lambda = -1, s = -1, scale = TRUE)

# Print clustering results
plot(result_mix, type = "dendrogram", xlab = "", sub = "", cex = 0.5)  # Plot dendrogram
plot(result_mix, type = "info", col = "black", pch = 16)  # Plot dendrogram

# Simulated categorical data example
set.seed(123)
data_cat <- data.frame(
  Var1 = as.factor(sample(letters[1:3], 200, replace = TRUE)),  # Nominal variable
  Var2 = as.factor(sample(letters[4:6], 200, replace = TRUE)),  # Nominal variable
  Var3 = factor(sample(c("low", "medium", "high"), 200, replace = TRUE),
                levels = c("low", "medium", "high"), ordered = TRUE)  # Ordinal variable
)

# Run AIBmix with automatic lambda selection 
result_cat <- AIBmix(X = data_cat, lambda = -1)

# Print clustering results
plot(result_cat, type = "dendrogram", xlab = "", sub = "", cex = 0.5)  # Plot dendrogram

# Results summary
summary(result_cat)

# Simulated continuous data example
set.seed(123)
# Continuous data with 200 observations, 5 features
data_cont <- as.data.frame(matrix(rnorm(1000), ncol = 5))

# Run AIBmix with automatic bandwidth selection 
result_cont <- AIBmix(X = data_cont, s = -1, scale = TRUE)

# Print concise summary ofoutput
print(result_cont)

# Print clustering results
plot(result_cont, type = "dendrogram", xlab = "", sub = "", cex = 0.5)  # Plot dendrogram
}

\author{
Efthymios Costa, Ioanna Papatsouma, Angelos Markos
}
\references{
{
\insertRef{slonim_aib_1999}{IBclust}

\insertRef{aitchison_kernel_1976}{IBclust}

\insertRef{li_nonparametric_2003}{IBclust}

\insertRef{silverman_density_1998}{IBclust}

\insertRef{ouyang2006cross}{IBclust}

\insertRef{wang1981class}{IBclust}

\insertRef{epanechnikov1969non}{IBclust}

\insertRef{costa_dib_2025}{IBclust}
}
}
\keyword{clustering}
