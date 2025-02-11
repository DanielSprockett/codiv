
#' Collapses monophyletic tips of a phylogenetic tree prior to `codiv()`
#'
#' @param symbiont_subtree A sub-node of the complete `Symbiont_tree`
#' @param Host_to_Symbiont_df A `data.frame()` that pairs all of the tips in the `Symbiont_tree` with the tip in the `Host_tree` it was isolated from.
#'
#' @returns
#' A phylogenetic tree with monophyletic nodes (i.e. all tips from the same host) are collapsed into 1 tip.
#'
#' @export
#' Host_tree <- ape::rtree(10)
#' Symbiont_tree <- ape::rtree(100)
#' Host_to_Symbiont_df <- data.frame("Host" = rep(Host_tree$tip.label, 10), 
#'                                     "Symbiont" = Symbiont_tree$tip.label)
#' @examples
#' collapse_monophyletic(Symbiont_tree, Host_to_Symbiont_df)
#'
collapse_monophyletic <- function(symbiont_subtree, Host_to_Symbiont_df) {
  tips_to_drop <- c()
  
  # traverse tree
  symbiont_subtree_list <- lapply(1:ape::Nnode(symbiont_subtree), function(x)
    castor::get_subtree_at_node(symbiont_subtree, x)$subtree)
  names(symbiont_subtree_list) <- symbiont_subtree$node.label
  
  # list all of the monophyletic tips from the same host, then select one to keep
  for (j in 1:length(symbiont_subtree_list)) {
    j_tree <- symbiont_subtree_list[[j]]
    j_hosts <- Host_to_Symbiont_df$Host[match(
      j_tree$tip.label,
      Host_to_Symbiont_df$Symbiont)]
    if (length(unique(j_hosts)) == 1) {
      tip_distances <- castor::get_all_distances_to_root(j_tree)[1:length(j_tree$tip.label)]
      tip_to_drop <- j_tree$tip.label[which(tip_distances != max(tip_distances))]
      tips_to_drop <- unique(c(tips_to_drop, tip_to_drop))
    }}
  
  collapsed_symbiont_subtree <- ape::drop.tip(symbiont_subtree, tips_to_drop)
  return(collapsed_symbiont_subtree)
}

