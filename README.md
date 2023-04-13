# codiv
Performs Host-Microbe Codiversification Scans 

Description: The codiv package performs codiversification scans between pairs of 
    host and symbiont phylogenetic trees. The main functions accepts a host tree,
    a symbiont tree, and a data.frame that links which symbionts were isolated
    from which hosts. It outputs a correlation coefficient for each node of the 
    symbiont tree denoted how well it's topology is correlated with the host tree. 
    Nodes with a high degree of correlation are consistent with codiversification
    between hosts and those symbionts. Other helper functions are also included. 
