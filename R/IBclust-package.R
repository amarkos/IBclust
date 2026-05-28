#' IBclust: Information Bottleneck Clustering for Mixed-Type Data
#'
#' The \pkg{IBclust} package provides clustering methods based on the
#' Information Bottleneck (IB) principle, supporting continuous, categorical
#' (nominal and ordinal), and mixed-type data. Four algorithmic variants are
#' implemented: \code{\link{AIBmix}} for hierarchical (agglomerative)
#' clustering, \code{\link{DIBmix}} for hard partitional clustering,
#' \code{\link{IBmix}} for soft (fuzzy) clustering, and \code{\link{GIBmix}}
#' for the generalised case interpolating between the two.
#' 
#' @section Choosing a method:
#' All four methods share the same kernel similarity machinery; they
#' differ in the form of the partition produced and the objective optimised:
#'
#' \itemize{
#'   \item \code{\link{AIBmix}}: hierarchical, hard clustering via
#'         greedy bottom-up merging that minimises Jensen-Shannon
#'         divergence between cluster distributions. Returns a full
#'         hierarchy from \eqn{n} singletons to one cluster.
#'   \item \code{\link{DIBmix}}: hard partitional clustering for a
#'         fixed number of clusters; uses dynamic regularisation to
#'         guarantee non-empty clusters.
#'   \item \code{\link{IBmix}}: soft (fuzzy) partitional clustering;
#'         each observation has a probability vector over clusters.
#'   \item \code{\link{GIBmix}}: generalised case with a tunable
#'         fuzziness parameter \eqn{\alpha \in [0, 1]} that interpolates
#'         between \code{DIBmix} (\eqn{\alpha = 0}) and \code{IBmix}
#'         (\eqn{\alpha = 1}).
#' }
#'
#' For mathematical details of each algorithm, see the references on
#' the individual function help pages, or
#' \insertCite{costa_dib_2025;textual}{IBclust} for a unified treatment.
#'
#' @section Kernel functions:
#' The functions in \pkg{IBclust} use kernel density estimation with
#' generalised product kernels to construct the matrix \eqn{P_{Y|X}} of
#' kernel weights. Bandwidth parameters control the smoothness of the
#' estimate and can be either user-supplied or automatically selected.
#'
#' For continuous variables:
#' \itemize{
#'   \item \emph{Gaussian (RBF) kernel \insertCite{silverman_density_1998}{IBclust}:}
#'   \deqn{K_c\left(\frac{x - x'}{s}\right) = \frac{1}{\sqrt{2\pi}} \exp\left\{-\frac{\left(x - x'\right)^2}{2s^2}\right\}, \quad s > 0.}
#'   \item \emph{Epanechnikov kernel \insertCite{epanechnikov1969non}{IBclust}:}
#'   \deqn{K_c(x - x'; s) = \begin{cases}
#'     \frac{3}{4\sqrt{5}}\left(1 - \frac{(x-x')^2}{5s^2} \right), & \text{if } \frac{(x - x')^2}{s^2} < 5 \\
#'     0, & \text{otherwise}
#' \end{cases}, \quad s > 0.}
#' }
#'
#' For nominal (unordered categorical) variables:
#' \itemize{
#' \item \emph{Aitchison & Aitken kernel \insertCite{aitchison_kernel_1976}{IBclust}:}
#' \deqn{K_u(x = x'; \lambda) = \begin{cases}
#'     1 - \lambda, & \text{if } x = x' \\
#'     \frac{\lambda}{\ell - 1}, & \text{otherwise}
#' \end{cases}, \quad 0 \leq \lambda \leq \frac{\ell - 1}{\ell}.}
#' \item \emph{Li & Racine kernel \insertCite{ouyang2006cross}{IBclust}:}
#' \deqn{K_u(x = x'; \lambda) = \begin{cases}
#'     1, & \text{if } x = x' \\
#'     \lambda, & \text{otherwise}
#' \end{cases}, \quad 0 \leq \lambda \leq 1.}
#' }
#'
#' For ordinal (ordered categorical) variables:
#' \itemize{
#' \item \emph{Li & Racine kernel \insertCite{li_nonparametric_2003}{IBclust}:}
#' \deqn{K_o(x = x'; \nu) = \begin{cases}
#'     1, & \text{if } x = x' \\
#'     \nu^{|x - x'|}, & \text{otherwise}
#' \end{cases}, \quad 0 \leq \nu \leq 1.}
#' \item \emph{Wang & van Ryzin kernel \insertCite{wang1981class}{IBclust}:}
#' \deqn{K_o(x = x'; \nu) = \begin{cases}
#'     1 - \nu, & \text{if } x = x' \\
#'     \frac{1-\nu}{2}\nu^{|x - x'|}, & \text{otherwise}
#' \end{cases}, \quad 0 \leq \nu \leq 1.}
#' }
#'
#' \eqn{\ell} is the number of levels of the categorical variable. For
#' ordinal variables, the \code{lambda} argument is used to define \eqn{\nu}.
#' Bandwidths are automatically selected using the approach in
#' \insertCite{costa_dib_2025;textual}{IBclust} when set to \eqn{-1}.
#' 
#' @section Nystr\enc{ö}{o}m approximation:
#' Constructing the kernel similarity matrix \eqn{P_{Y|X}} has cost
#' \eqn{O(n^2)}, which becomes prohibitive for large datasets. The
#' partitional methods (\code{\link{DIBmix}}, \code{\link{IBmix}},
#' \code{\link{GIBmix}}) support a Nystr\enc{ö}{o}m approximation
#' \insertCite{williams2000using}{IBclust} of the form
#' \eqn{P_{Y|X} \approx C W C^\top}, where \eqn{C} is an \eqn{n \times m}
#' matrix of kernel weights between all observations and \eqn{m}
#' landmark points, and \eqn{W} is the \eqn{m \times m} kernel matrix
#' over the landmarks alone. This reduces the cost of constructing
#' \eqn{P_{Y|X}} from \eqn{O(n^2)} to \eqn{O(nm)}, with \eqn{m \ll n}.
#' The heuristic \eqn{m \approx \sqrt{n}} is used as the default.
#' The approximation is enabled by setting \code{nystrom = TRUE} in the fit
#' functions; the number of landmark points is controlled via \code{n_landmarks},
#' and specific landmark observations can be requested via \code{landmark_indices}.
#' The Nystr\enc{ö}{o}m approximation is not available for
#' \code{\link{AIBmix}}, since its greedy merging procedure has cost
#' \eqn{O(n^3)} regardless of how \eqn{P_{Y|X}} is constructed.
#'
#'
#' @references
#' \insertRef{aitchison_kernel_1976}{IBclust}
#'
#' \insertRef{li_nonparametric_2003}{IBclust}
#'
#' \insertRef{silverman_density_1998}{IBclust}
#'
#' \insertRef{ouyang2006cross}{IBclust}
#'
#' \insertRef{wang1981class}{IBclust}
#'
#' \insertRef{epanechnikov1969non}{IBclust}
#'
#' \insertRef{costa_dib_2025}{IBclust}
#' 
#' \insertRef{williams2000using}{IBclust}
#'
#' @useDynLib IBclust, .registration = TRUE
#' @import stats
#' @import np
#' @import Rdpack
#' @importFrom utils head flush.console txtProgressBar setTxtProgressBar tail
#' @importFrom Rcpp evalCpp
#' @importFrom rje powerSet
#' @importFrom graphics barplot points box image rect axis lines
#' @importFrom grDevices colorRampPalette hcl.colors adjustcolor
"_PACKAGE"