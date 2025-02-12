#' Co-Diversification Scan Between Hosts and Symbionts
#'
#' @param Host_tree 
#' @param Symbiont_tree 
#' @param Host_to_Symbiont_df 
#' @param min_hosts 
#' @param min_symbiont_tips 
#' @param span_fraction 
#' @param permutations 
#' @param seed 
#' @param verbose 
#' @param Save_fp 
#' @param methods 
#' @param focus_hosts 
#' @param subtree_features 
#' @param cores 
#' @param continue 
#'
#' @returns
#' @export
#'
#' @examples
codiv <- function(Host_tree, Symbiont_tree, Host_to_Symbiont_df,
                  min_hosts = 3, min_symbiont_tips = 7, span_fraction = 0.1, 
                  permutations = 99, seed = 8675309,
                  verbose = TRUE, Save_fp, methods = c("hommola", "paco", "parafit"), 
                  focus_hosts, subtree_features = TRUE, cores = 4, continue = TRUE) {
  
  ## check inputs ##
  check_inputs(Host_tree, Symbiont_tree, Host_to_Symbiont_df, min_hosts, min_symbiont_tips, span_fraction, permutations)
  
  ## SYMBIONT TREE ##
  
  # Fix possible issues with input Symbiont_tree
  
  # make node names unique if they are not
  if (sum(duplicated(Symbiont_tree$node.label)) > 0) {
    if (verbose) {print("Prepending Node Labels with unique numbers")}
    Symbiont_tree$node.label <- paste0(1:Symbiont_tree$Nnode, "_",
                                       Symbiont_tree$node.label)
  }
  
  # add node labels if tree does not have them
  if (is.null(Symbiont_tree$node.label)) {
    if (verbose) {message("Creating unique Node Labels")}
    Symbiont_tree$node.label <- paste0("Node_", 1:Symbiont_tree$Nnode)
  }
  
  # Identify which nodes are going to be included in scan
  
  # get list of all Symbiont_tree subtrees
  subtree_list <- lapply(1:Nnode(Symbiont_tree), function(x)
    castor:::get_subtree_at_node(Symbiont_tree, x)$subtree)
  names(subtree_list) <- Symbiont_tree$node.label
  
  # calculate max tree span for all Symbiont_tree subtrees
  span_list <- unlist(lapply(1:Nnode(Symbiont_tree), function(x)
    castor:::get_tree_span(subtree_list[[x]], 
                           as_edge_count = FALSE)$max_distance))
  
  symbiont_tip_count_list <- castor:::count_tips_per_node(Symbiont_tree)
  
  host_tip_count_list <- unlist(lapply(1:Nnode(Symbiont_tree), function(x)
    length(unique(Host_to_Symbiont_df$Host[
      match(subtree_list[[x]]$tip.label, 
            Host_to_Symbiont_df$Symbiont)]))))
  
  # generate a data.frame of all symbiont nodes
  Symbiont_df <- data.frame(
    "Node_Label" = Symbiont_tree$node.label,
    "node" = 1:length(Symbiont_tree$node.label),
    "Subtree_Span" = span_list,
    "Subtree_Symbiont_Tips" = symbiont_tip_count_list,
    "Subtree_Host_Tips" = host_tip_count_list)
  
  # identify nodes that have sufficient tree span to include in the codiv scan
  Symbiont_Max_Span <- max(Symbiont_df$Subtree_Span)
  Span_CutOff <- span_fraction * Symbiont_Max_Span
  # span_3rd_quant <- quantile(Symbiont_df$Subtree_Span) 
  
  n1 <- dim(Symbiont_df)[1]
  Symbiont_df <- Symbiont_df[Symbiont_df$Subtree_Span <= Span_CutOff, ]
  n2 <- dim(Symbiont_df)[1]
  Symbiont_df <- Symbiont_df[Symbiont_df$Subtree_Symbiont_Tips >= min_symbiont_tips, ]
  n3 <- dim(Symbiont_df)[1]
  Symbiont_df <- Symbiont_df[Symbiont_df$Subtree_Host_Tips >= min_hosts, ]
  n4 <- dim(Symbiont_df)[1]
  
  if (verbose) {
    print(paste0("Of the ", n1, 
                 " internal nodes in the symbiont tree,"))
    print(paste0("   ", n2,
                 " (or ", label_percent(big.mark = ",", suffix = "%")(n2/n1), 
                 ") of them have a span at least " , 
                 label_percent(big.mark = ",", suffix = "%")(span_fraction),
                 " of that of the entire tree."))
    print(paste0("   Of those, ", n3,
                 " (or ", label_percent(big.mark = ",", suffix = "%")(n3/n1), 
                 " of total) of them have at least " , 
                 min_symbiont_tips,
                 " decendant tips."))
    print(paste0("   Of those, ", n4,
                 " (or ", label_percent(big.mark = ",", suffix = "%")(n4/n1), 
                 " of total) of them have tips from at least " , 
                 min_hosts,
                 " different hosts."))
    print(paste0("These ", n4, " Nodes in Symbiont_Tree will be scanned for evidence of codiversification."))
  }
  
  Symbiont_df <- Symbiont_df[order(Symbiont_df$Subtree_Symbiont_Tips, Symbiont_df$Subtree_Host_Tips),]
  nodes_to_scan <- Symbiont_df$Node_Label
  
  # create Results_df data.frame which will be output
  
  # set up which columns should be in the output
  basic_vars <- c("Node_ID", 
                  "Symbiont_Tree", "N_Symbionts", "Symbiont_Colless", "Symbiont_Sackin", 
                  "Host_Tree", "N_Hosts", "Host_Colless", "Host_Sackin")
  Hommola_vars <- c()
  PACo_vars <- c()
  ParaFit_vars <- c()
  subtree_vars <- c()
  if ("hommola" %in% methods) 
  {Hommola_vars <- c("Hommola_r", "Hommola_pvalue", 
                     "Collapsed_Hommola_r", "Collapsed_Hommola_pvalue")}
  if ("paco" %in% methods) 
  {PACo_vars <- c("Collapsed_PACo_ss", "Collapsed_PACo_pvalue")}
  if ("parafit" %in% methods) 
  {ParaFit_vars <- c("Collapsed_ParaFitGlobal", "Collapsed_ParaFit_pvalue")}
  if (subtree_features)
  {subtree_vars <- c("Collapsed_TreeDistance", "Collapsed_SharedPhylogeneticInfo", 
                     "Collapsed_DifferentPhylogeneticInfo", "Collapsed_NyeSimilarity", 
                     "Collapsed_JaccardRobinsonFoulds", "Collapsed_MatchingSplitDistanc", 
                     "Collapsed_MatchingSplitInfoDistance", "Collapsed_MutualClusteringInfo")}
  if (length(focus_hosts) > 0)
  {focus_hosts_vars <- paste0(focus_hosts, c("_PRESENT"))}
  
  col_vars <- c(basic_vars, Hommola_vars, PACo_vars, ParaFit_vars, subtree_vars, focus_hosts_vars)
  
  # create the output data.frame to hold the results of the codiv scan
  Results_df <- data.frame(matrix(nrow = length(nodes_to_scan), ncol = length(col_vars)))
  colnames(Results_df) <- col_vars
  rownames(Results_df) <- nodes_to_scan
  Results_df$Node_ID <- nodes_to_scan
  
  if (continue && file.exists(Save_fp)) {
    temp_results <- suppressWarnings(read.table(Save_fp, sep = "\t", header = TRUE, fill = TRUE))
    completed <- temp_results$Node_ID[which(!is.na(temp_results$Hommola_pvalue))]
    nodes_to_scan <- nodes_to_scan[!(nodes_to_scan %in% completed)]
    temp_results <- temp_results[temp_results$Node_ID %in% completed, ]
    temp_res <- data.frame(matrix(nrow = length(nodes_to_scan), ncol = length(col_vars)))
    colnames(temp_res) <- col_vars
    rownames(temp_res) <- nodes_to_scan
    temp_res$Node_ID <- nodes_to_scan
    Results_df <- rbind(temp_results, temp_res)
  }
  
  write.table(Results_df, Save_fp, sep = "\t", quote = FALSE, row.names = FALSE)
  n_iter <- length(nodes_to_scan)
  
  # progress bar
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = n_iter,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar
  
  # remove after testing 
  #i <- nodes_to_scan[500]
  for (i in nodes_to_scan) { 
    
    # Updates the current state for progress bar
    pb$tick()
    
    # isolate each node as a subtree
    i_symbiont_subtree <- subtree_list[[i]]
    # reduce Host_to_Symbiont_df and Host_tr to only hosts relevant for this subtree
    i_host_to_symbiont_df <- Host_to_Symbiont_df[Host_to_Symbiont_df$Symbiont %in%  i_symbiont_subtree$tip.label,]
    i_host_subtree <- keep.tip(Host_tr, unique(i_host_to_symbiont_df$Host))
    
    # convert data.frame to matrix of 0's or 1's denoting which symbionts are in which hosts
    i_host_to_symbiont_mat <- host_symbiont_links(i_host_to_symbiont_df)
    
    # calculate distance matrices for i_symbiont_subtree and i_host_subtree
    # inputs to functions
    i_symbiont_subtree_dist <- adephylo:::distTips(i_symbiont_subtree, tips = "all", 
                                                   method = c("patristic"), 
                                                   useC = TRUE)
    i_host_subtree_dist <- adephylo:::distTips(i_host_subtree, tips = "all", 
                                               method = c("patristic"), 
                                               useC = TRUE)
    
    # collapse all nodes in the i_symbiont_subtree that are symbionts from the same host
    i_collapsed_symbiont_subtree <- collapse_monophyletic(i_symbiont_subtree, i_host_to_symbiont_df)
    i_collapsed_symbiont_subtree_dist <- adephylo:::distTips(i_collapsed_symbiont_subtree, tips = "all", 
                                                             method = c("patristic"), 
                                                             useC = TRUE)
    i_collapsed_host_to_symbiont_df <- Host_to_Symbiont_df[Host_to_Symbiont_df$Symbiont %in% 
                                                             i_collapsed_symbiont_subtree$tip.label,]
    i_collapsed_host_to_symbiont_mat <- host_symbiont_links(i_collapsed_host_to_symbiont_df)
    
    if ("hommola" %in% methods) {
      hommola_results <- hommola_wf(i_host_subtree_dist, i_symbiont_subtree_dist, 
                                    i_host_to_symbiont_df, permutations, seed)
      hommola_r <- hommola_results$Hommola_r
      hommola_pvalue <- hommola_results$Hommola_pvalue
      
      Collapsed_hommola_results <- hommola_wf(i_host_subtree_dist, i_collapsed_symbiont_subtree_dist,
                                              i_host_to_symbiont_df, permutations, seed)
      Collapsed_hommola_r <- Collapsed_hommola_results$Hommola_r
      Collapsed_hommola_pvalue <- Collapsed_hommola_results$Hommola_pvalue
    }
    
    if ("paco" %in% methods) {
      paco_results <- paco_wf(i_host_subtree, i_collapsed_symbiont_subtree, 
                              i_collapsed_host_to_symbiont_mat,
                              permutations, seed)
      Collapsed_PACo_ss <- paco_results$gof$ss
      Collapsed_PACo_pvalue <- paco_results$gof$p
    }
    
    if ("parafit" %in% methods) {
      parafit_results <- parafit(i_host_subtree_dist, i_collapsed_symbiont_subtree_dist, 
                                 i_collapsed_host_to_symbiont_mat, 
                                 nperm = permutations, test.links = FALSE,
                                 seed = seed, correction = "cailliez", silent = TRUE)
      Collapsed_parafit_ParaFitGlobal <- parafit_results$ParaFitGlobal
      Collapsed_parafit_pvalue <- parafit_results$p.global
    }
    
    
    # calculate features of i_symbiont_subtree and i_host_subtree 
    i_symbiont_subtree_newick <- write.tree(i_symbiont_subtree)
    i_symbiont_subtree_n_tips <- length(i_symbiont_subtree$tip.label)
    i_symbiont_subtree_colless <- colless(as.treeshape(i_symbiont_subtree))
    i_symbiont_subtree_sackin <- sackin(as.treeshape(i_symbiont_subtree))
    
    i_host_subtree_newick <- write.tree(i_host_subtree)
    i_host_subtree_n_tips <- length(i_host_subtree$tip.label)
    i_host_subtree_colless <- colless(as.treeshape(i_host_subtree))
    i_host_subtree_sackin <- sackin(as.treeshape(i_host_subtree))
    
    # Identify if the focus_hosts are in the host_subtree
    if (length(focus_hosts_vars) > 0) {
      for (j in 1:length(focus_hosts_vars)) {
        focus_hosts_present <- ifelse(focus_hosts[j] %in% i_host_to_symbiont_df$Host, TRUE, FALSE)
      }}
    
    if (subtree_features) {
      ##  tree distance metrics
      # see https://github.com/ms609/TreeDist for details:
      # must first rename the symbiont subtree tips with the hosts
      f_symbiont_subtree <- i_collapsed_symbiont_subtree 
      f_symbiont_subtree <- reorder(f_symbiont_subtree, order = "cladewise")
      f_symbiont_subtree$tip.label <- i_collapsed_host_to_symbiont_df$Host[
        match(f_symbiont_subtree$tip.label,
              i_collapsed_host_to_symbiont_df$Symbiont)]
      
      i_TreeDistance <- TreeDistance(i_host_subtree, f_symbiont_subtree)
      i_SharedPhylogeneticInfo <- SharedPhylogeneticInfo(i_host_subtree, f_symbiont_subtree)
      i_DifferentPhylogeneticInfo <- DifferentPhylogeneticInfo(i_host_subtree, f_symbiont_subtree)
      i_NyeSimilarity <- NyeSimilarity(i_host_subtree, f_symbiont_subtree)
      i_JaccardRobinsonFoulds <- JaccardRobinsonFoulds(i_host_subtree, f_symbiont_subtree)
      i_MatchingSplitDistance <- MatchingSplitDistance(i_host_subtree, f_symbiont_subtree)
      i_MatchingSplitInfoDistance <- MatchingSplitInfoDistance(i_host_subtree, f_symbiont_subtree) 
      i_MutualClusteringInfo <- MutualClusteringInfo(i_host_subtree, f_symbiont_subtree)
    }
    
    ## WRITE TO FILE
    # Tree Features
    Results_df[which(Results_df$Node_ID == i), "Symbiont_Tree"] <- i_symbiont_subtree_newick
    Results_df[which(Results_df$Node_ID == i), "Symbiont_Colless"] <- i_symbiont_subtree_colless
    Results_df[which(Results_df$Node_ID == i), "Symbiont_Sackin"] <- i_symbiont_subtree_colless
    Results_df[which(Results_df$Node_ID == i), "N_Symbionts"] <- i_symbiont_subtree_n_tips
    
    Results_df[which(Results_df$Node_ID == i), "Host_Tree"] <- i_host_subtree_newick
    Results_df[which(Results_df$Node_ID == i), "Host_Colless"] <- i_host_subtree_colless
    Results_df[which(Results_df$Node_ID == i), "Host_Sackin"] <- i_host_subtree_sackin
    Results_df[which(Results_df$Node_ID == i), "N_Hosts"] <- i_host_subtree_n_tips
    
    # Hommola Results
    if ("hommola" %in% methods) {
      Results_df[which(Results_df$Node_ID == i), "Hommola_r"] <- hommola_r
      Results_df[which(Results_df$Node_ID == i), "Hommola_pvalue"] <- hommola_pvalue
      Results_df[which(Results_df$Node_ID == i), "Collapsed_Hommola_r"] <- Collapsed_hommola_r
      Results_df[which(Results_df$Node_ID == i), "Collapsed_Hommola_pvalue"] <- Collapsed_hommola_pvalue
    }
    
    # PACo Results
    if ("paco" %in% methods) {
      Results_df[which(Results_df$Node_ID == i), "Collapsed_PACo_ss"] <- Collapsed_PACo_ss
      Results_df[which(Results_df$Node_ID == i), "Collapsed_PACo_pvalue"] <- Collapsed_PACo_pvalue
    }
    
    # parafit Results
    if ("parafit" %in% methods) {
      Results_df[which(Results_df$Node_ID == i), "Collapsed_ParaFitGlobal"] <- Collapsed_parafit_ParaFitGlobal
      Results_df[which(Results_df$Node_ID == i), "Collapsed_ParaFit_pvalue"] <- Collapsed_parafit_pvalue
    }
    
    if (subtree_features) {
      # Alternative Tree Distance Metrics
      Results_df[which(Results_df$Node_ID == i), "Collapsed_TreeDistance"] <- i_TreeDistance
      Results_df[which(Results_df$Node_ID == i), "Collapsed_SharedPhylogeneticInfo"] <- i_SharedPhylogeneticInfo
      Results_df[which(Results_df$Node_ID == i), "Collapsed_DifferentPhylogeneticInfo"] <- i_DifferentPhylogeneticInfo
      Results_df[which(Results_df$Node_ID == i), "Collapsed_NyeSimilarity"] <- i_NyeSimilarity
      Results_df[which(Results_df$Node_ID == i), "Collapsed_JaccardRobinsonFoulds"] <- i_JaccardRobinsonFoulds
      Results_df[which(Results_df$Node_ID == i), "Collapsed_MatchingSplitDistanc"] <- i_MatchingSplitDistance
      Results_df[which(Results_df$Node_ID == i), "Collapsed_MatchingSplitInfoDistance"] <- i_MatchingSplitInfoDistance
      Results_df[which(Results_df$Node_ID == i), "Collapsed_MutualClusteringInfo"] <- i_MutualClusteringInfo
    }
    
    # focus_hosts
    Results_df[which(Results_df$Node_ID == i), focus_hosts_vars] <- focus_hosts_present
    
    if (file.exists(Save_fp)) {
      write.table(Results_df, Save_fp, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
  return(Results_df)
}
