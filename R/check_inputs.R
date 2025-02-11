#' Checks inputs to `codiv()`
#'
#' @param Host_tree A binary phylogenetic tree of hosts
#' @param Symbiont_tree A binary phylogenetic tree of symbionts
#' @param Host_to_Symbiont_df A `data.frame()` that pairs all of the tips in the `Symbiont_tree` with the tip in the `Host_tree` it was isolated from.
#' @param min_hosts The minimum number hosts for a node to be evaluated (must be greater than 3)
#' @param min_symbiont_tips The minimum number symbiont tips for a node to be evaluated (recommend 7 or more)
#' @param span_fraction The fraction of the `Symbiont_tree` included in the codiv scan
#' @param permutations The number of permutations (recommend at least 99)
#'
#' @returns Tests to make sure the inputs are acceptable
#' @export
#'
#' @examples
#' Host_tree <- ape::rtree(10)
#' Symbiont_tree <- ape::rtree(100)
#' Host_to_Symbiont_df <- data.frame("Host" = rep(Host_tree$tip.label, 10), 
#'                                     "Symbiont" = Symbiont_tree$tip.label)
#' check_inputs(Host_tree, Symbiont_tree, Host_to_Symbiont_df, 3, 7, 0.25, 99)
#'
check_inputs <- function(Host_tree, Symbiont_tree, Host_to_Symbiont_df, 
                         min_hosts, min_symbiont_tips, span_fraction, 
                         permutations) {
  if (methods::is(Host_tree, "phylo") | methods::is(Symbiont_tree, "phylo")) {
    message("`Host_tree` and `Symbiont_tree` Inputs must be binary phylogentic trees")
  }
  
  if (min_hosts < 3) {message("`min_hosts` must be greater than 3")}
  if (span_fraction > 1) {message("`span_fraction` must be between 0 and 1")}
  if (min_symbiont_tips < 7) {message("We recommend setting `min_symbiont_tips` to 7 or more")}
  if (permutations < 99) {message("We recommend setting `permutations` to 99 or more")}
  
  if (sum(c("Host", "Symbiont") %in% colnames(Host_to_Symbiont_df)) > 2) {
    message("`Host_to_Symbiont_df` must contain columns labeled `Host` and `Symbiont`")
  }
  
  if (sum(Symbiont_tree$tip.label %in% Host_to_Symbiont_df$Host) > 0) {
    message("`Symbiont_tree` contains tips that are not in`Host_to_Symbiont_df`")
  }
  
}
