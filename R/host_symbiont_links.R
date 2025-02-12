
#' Changes the host-to-symbiont `data.frame()` to wide-format
#'
#' @param Host_to_Symbiont_df A `data.frame()` that includes columns named 
#' "Host" and "Symbiont"
#'
#' @returns
#' A wide-format host-to-symbiont `data.frame()`
#' 
#' @export
#'
#' @examples
#' host_symbiont_links(Host_to_Symbiont_df)
#' 
host_symbiont_links <- function(Host_to_Symbiont_df) {
  if (!sum(colnames(Host_to_Symbiont_df) %in% c("Host", "Symbiont")) == 2) {
    message("colnames(Host_to_Symbiont_df) must include both 'Host' and 'Symbiont'")
  }
  
  if (length(unique(Host_to_Symbiont_df$Symbiont)) != length(Host_to_Symbiont_df$Symbiont)) {
    message("One or more 'Symbionts' are duplicated in your Host_to_Symbiont_df file.")
  }
  # make a dataframe to hold results
  host_symb_df <- data.frame(matrix(0, 
                                    nrow = length(unique(Host_to_Symbiont_df$Host)),
                                    ncol = length(unique(Host_to_Symbiont_df$Symbiont))))
  rownames(host_symb_df) <- unique(Host_to_Symbiont_df$Host)
  colnames(host_symb_df) <- unique(Host_to_Symbiont_df$Symbiont)
  for (i in unique(Host_to_Symbiont_df$Symbiont)) {
    i_host <- Host_to_Symbiont_df$Host[which(Host_to_Symbiont_df$Symbiont == i)]
    host_symb_df[
      which(rownames(host_symb_df) == i_host),
      which(colnames(host_symb_df) == i)
    ] <- 1
  }
  return(host_symb_df)
}

