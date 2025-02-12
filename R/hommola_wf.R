
#' Calculates the Hommola Correlation with permutations
#'
#' @param i_host_subtree_dist 
#' @param i_symbiont_subtree_dist 
#' @param i_host_to_symbiont_df 
#' @param permutations The number of permutations (recommend at least 99)
#' @param seed A numeric, used to `set.seed()`
#'
#' @returns
#' @export
#'
#' @examples
hommola_wf <- function(i_host_subtree_dist, i_symbiont_subtree_dist, 
                       i_host_to_symbiont_df, permutations, seed) {
  
  hommola_res <- hommola(i_host_subtree_dist, i_symbiont_subtree_dist, i_host_to_symbiont_df)
  
  # PERMUTATAION
  # generate a matrix of randomize indicies that effectively permute the 
  # link between the symbionts and their hosts, while maintaining
  # the distribution of genetic distances
  
  # This generates a matrix n permutations wide and n MAGs long.
  # Used as indexes for shuffling the Host_to_Symbiont_df relationships
  if (!is.na(seed)) set.seed(seed)
  permutations_mat <- replicate(permutations, sample(1:nrow(i_host_to_symbiont_df)))
  
  # This generates a list of length (n permutations)
  # each entry in the list is a data.frame similar to the Host_to_Symbiont_df, 
  # but the column of symbionts has been randomly shuffled,
  # according to the previous indexes in the permutations_mat
  i_host_to_symbiont_permutations_df <- lapply(1:permutations, function(x) 
    data.frame("Host" = i_host_to_symbiont_df[, "Host"],
               "Symbiont" = i_host_to_symbiont_df[permutations_mat[,x], "Symbiont"]))
  
  # aside: does this work? N<-unique(N); N[sample(nrow(N), 10000), ] to get 10000 unique permutations.
  # see also permutation:Permutation In Rfast: A Collection of Efficient and Extremely Fast R Functions
  
  hommola_res_perms <- vector("list", permutations)
  for (i in 1:permutations) {
    hommola_res_perms[[i]] <- hommola(i_host_subtree_dist, i_symbiont_subtree_dist,
                                      i_host_to_symbiont_permutations_df[[i]])
  }
  p_value = (sum(unlist(hommola_res_perms) >= hommola_res) + 1) / (permutations + 1)
  hommola_out <- vector("list", 2)
  hommola_out[[1]] <- hommola_res
  hommola_out[[2]] <- p_value
  names(hommola_out) <- c("Hommola_r", "Hommola_pvalue")
  return(hommola_out)
  
}

