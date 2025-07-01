make_dendrogram <- function(merges, costs, labels = NULL) {
  # merges: (n-1)×2 matrix of cluster‐IDs merged at each step
  # costs : length-(n-1) vector of merge "heights"
  # labels : optional vector of length n for leaf names
  n <- nrow(merges) + 1
  # mapping from static cluster‐IDs to hclust convention
  # unmerged leaves: -i
  # a new cluster at step s: +s
  cur_label <- -seq_len(n)
  # build the hclust merging matrix
  merge_mat <- matrix(NA_integer_, nrow = n-1, ncol = 2)
  for(s in seq_len(n-1)) {
    i_idx <- merges[s, 1]
    j_idx <- merges[s, 2]
    # look up their current hclust‐labels
    a <- cur_label[i_idx]
    b <- cur_label[j_idx]
    merge_mat[s, ] <- c(a, b)
    # assign the newly formed cluster its label +s
    cur_label[i_idx] <- s
    # j_idx is now dead but we leave its label alone
  }
  # heights come straight from merge costs
  heights <- costs
  # compute a default plotting order by walking the tree
  get_order <- function(node){
    if (node < 0) return(-node)
    children <- merge_mat[node, ]
    c(get_order(children[1]), get_order(children[2]))
  }
  ord <- get_order(n-1)
  
  # assemble the hclust object
  hc <- list(
    merge       = merge_mat,
    height      = heights,
    order       = ord,
    labels      = labels %||% seq_len(n),
    method      = "IB",
    call        = match.call(),
    dist.method = "information_bottleneck"
  )
  class(hc) <- "hclust"
  hc
}
