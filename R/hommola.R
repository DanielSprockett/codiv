
#' Calculates the correlation between host and symbiont trees
#'
#' @param i_host_subtree_dist 
#' @param i_symbiont_subtree_dist 
#' @param i_host_to_symbiont_df 
#'
#' @returns
#' @export
#'
#' @examples
hommola <- function(i_host_subtree_dist, i_symbiont_subtree_dist, i_host_to_symbiont_df) {
  
  dist_df <- as_tibble(as.matrix(i_symbiont_subtree_dist), rownames = "Symbiont_1") %>% 
    gather(Symbiont_2, Symbiont_Distance, -Symbiont_1)
  # remove self comparisons
  dist_df <- dist_df[!dist_df$Symbiont_1 == dist_df$Symbiont_2, ]
  # look up hosts using the matrix
  dist_df$Host_1 <- i_host_to_symbiont_df$Host[match(dist_df$Symbiont_1, i_host_to_symbiont_df$Symbiont)]
  dist_df$Host_2 <- i_host_to_symbiont_df$Host[match(dist_df$Symbiont_2, i_host_to_symbiont_df$Symbiont)]
  
  Host_dist <- as_tibble(as.matrix(i_host_subtree_dist), rownames = "Host_1") %>% 
    gather(Host_2, Host_Distance, -Host_1)
  dist_df <- suppressMessages(left_join(dist_df, Host_dist))
  
  # calcualte correlation between Symbiont Distances and Host Distances
  # address a bug where the permutation randomly assigned all the same hosts, making all Host_Distances = 0
  if (length(unique(dist_df$Host_Distance)) == 1) {
    i_correlation <- 0
  } else {
    i_correlation <- cor(dist_df$Symbiont_Distance, dist_df$Host_Distance, method = 'pearson')
  }
  return(i_correlation)
}

