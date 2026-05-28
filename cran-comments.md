## Submission 1.4
* Added `fitted`, `info_metrics`, and `coef` methods for objects of both `aibclust` and `gibclust` class.
* Added `predict` method for `gibclust` class objects.
* Objects of class `aibclust` can now be converted to `hclust` objects through `as.hclust`.
* Added optional `keep_data` argument to all clustering functions; when `TRUE`, the input data is stored in the returned object.
* Added new information-theoretic metrics (H(T), H(T|X), I(T;X)) to `aibclust` objects, accessible via `info_metrics`.
* Included a new type of plot for both `aibclust` and `gibclust` objects that draws the similarity matrix.
* Included a parallel coordinates membership plot for `gibclust` objects to visualise fuzzy cluster memberships.
* Added the function `find_elbow` that detects the knee/elbow of a curve (useful for information retention curves).
* Updated documentation to avoid repetitions.

## Submission 1.3
* Added an implementation of the Nyström approximation for large kernel Gram matrices.
* Optimised continuous bandwidth search based on average nearest and furthest neighbour heuristics.
* Added variable importance plot in plotting method.
* Added progress bar for AIBmix.
* Minor bug fixes.

## Submission 1.2.1
* Added logo.
* Added vignette.
* More kernel functions supported.
* Input argument `cat_first` added for `AIBmix`, `DIBmix`, `GIBmix` & `IBmix` functions.
* Changed license; MIT to GPL >= 3.
* Added input checking & preprocessing function `input_checks.R`.
* Only exporting `AIBmix`, `DIBmix`, `GIBmix`, `IBmix` now for all cases of continuous/categorical/mixed-type data; rest of the functions deprecated.
* Outputs changed to be of `aibclust` or `gibclust` class supporting `print`, `summary`, and `plot` methods.

## Submission 1.2
* Added Agglomerative IB.
* Added Generalised IB.
* Added standard IB.