% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collapse_monophyletic.R
\name{collapse_monophyletic}
\alias{collapse_monophyletic}
\title{Collapses monophyletic tips of a phylogenetic tree prior to \code{codiv()}}
\usage{
collapse_monophyletic(symbiont_subtree, Host_to_Symbiont_df)
}
\arguments{
\item{symbiont_subtree}{A sub-node of the complete \code{Symbiont_tree}}

\item{Host_to_Symbiont_df}{A \code{data.frame()} that pairs all of the tips in the \code{Symbiont_tree} with the tip in the \code{Host_tree} it was isolated from.}
}
\value{
A phylogenetic tree with monophyletic nodes (i.e. all tips from the same host) are collapsed into 1 tip.
}
\description{
Collapses monophyletic tips of a phylogenetic tree prior to \code{codiv()}
}
\examples{
collapse_monophyletic(Symbiont_tree, Host_to_Symbiont_df)

}
