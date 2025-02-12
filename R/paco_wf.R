
#' Runs PACo
#'
#' @param i_host_subtree 
#' @param i_symbiont_subtree 
#' @param i_host_to_symbiont_mat 
#' @param permutations 
#' @param seed 
#'
#' @returns
#' @export
#'
#' @examples
paco_wf <- function(i_host_subtree, i_symbiont_subtree, i_host_to_symbiont_mat, permutations, seed) {
  gdist <- cophenetic(i_host_subtree)
  ldist <- cophenetic(i_symbiont_subtree)
  D <- prepare_paco_data(gdist, ldist, i_host_to_symbiont_mat)  
  D <- add_pcoord(D)
  D <- suppressWarnings(PACo(D, nperm = permutations, seed = seed, method = "swap"))
  return(D)
}

