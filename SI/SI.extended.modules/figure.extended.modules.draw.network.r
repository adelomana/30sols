###
### script to build interaction networks with RPG corems
###

draw.subnetwork <- function(listgenes) {
  library(RCy3)
  g <- new('graphNEL')
  ##Draw nodes
  ##gene nodes
  genes <- unlist(listgenes)
  for (i in genes) {
    g <- graph::addNode(node = i, object = g)
  }
  #cluster nodes
  clusters <- unique(unlist(lapply(listgenes, corem.finder)))
  for (j in clusters) {
    g <- graph::addNode(node = j, object = g)
  }
  cw <- CytoscapeWindow(title = 'RStudio', graph =  g, overwriteWindow = TRUE)
  displayGraph(cw)
  g <- cw@graph
  ##Node Attributes
  g <- initNodeAttribute (graph=g, attribute.name='label', attribute.type='char', default.value='undefined')
  g <- initNodeAttribute (graph=g, attribute.name='nodeType', attribute.type='char', default.value='undefined')
  for (i in genes) {
    nodeData(g, i, 'label') <- i
    nodeData(g, i, 'nodeType') <- 'gene'
  }
  for (j in clusters) {
    nodeData(g, j, 'label') <- j
    nodeData(g, j, 'nodeType') <- 'module'
  }
  ## gene names
  for (k in genes) {
    if (k %in% geneNames[,1]) { nodeData(g, k, 'label') <- geneNames[k, 2] }
    if (k %in% old2newNames[1:35,2]) { nodeData(g, k, 'label') <- geneNames[old2newNames[old2newNames[,2] %in% k,1], 2] }
  }
  ##Exp corems
  for (i in RPcoremExp) {
    nodeData(g, i, 'nodeType') <- 'exp.module'
  }
  ##Sat corems
  for (i in RPcoremSat) {
    nodeData(g, i, 'nodeType') <- 'sat.module'
  }
  ## aaRS 
  for (i in aarsVNG) {
    nodeData(g, i, 'nodeType') <- 'aars.gene'
  }
  ## processivity factors
  for (i in EIRfVNG) {
    nodeData(g, i, 'nodeType') <- 'eirf.gene'
  }
  ## RNase
  for (i in RNaseVNG) {
    nodeData(g, i, 'nodeType') <- 'rnase.gene'
  }
  ## RNA bp
  for (i in RBPVNG) {
    nodeData(g, i, 'nodeType') <- 'rbp.gene'
  }
  ## protease
  for (i in proteaseVNG) {
    nodeData(g, i, 'nodeType') <- 'protease.gene'
  }
  ## chap.skel
  for (i in c(secVNG,chap.skelVNG)) {
    nodeData(g, i, 'nodeType') <- 'chap.skel.sec.gene'
  }
  ## DNA subunits
  #for (i in dnaVNG) {
  #  nodeData(g, i, 'nodeType') <- 'dna.gene'
  #}
  ## RNAPol subunits
  for (i in RNApolVNG) {
    nodeData(g, i, 'nodeType') <- 'RNApol.gene'
  }
  ## ribosome subunit numbers
  for (i in riboVNG) {
    nodeData (g, i, 'label') <- riboSubunits[i,3]
    nodeData(g, i, 'nodeType') <- 'rp.gene'
  }
  cw <- setGraph(cw, g)
  displayGraph(cw)
  #print(noa(g, 'label'))
  #print(noa(g, 'nodeType'))
  #print(noa(g, 'edges'))
  ### Node attribute viz rules
  attribute.values <- c('gene', 'module', 'exp.module', 'sat.module','rp.gene', 'RNApol.gene', 'aars.gene', 'eirf.gene', 'rnase.gene', 'rbp.gene', 'protease.gene', 'chap.skel.sec.gene')
  node.shapes <- c('ROUND_RECTANGLE', 'DIAMOND', 'DIAMOND', 'DIAMOND', 'HEXAGON', 'HEXAGON', 'HEXAGON', 'HEXAGON', 'HEXAGON', 'HEXAGON', 'HEXAGON', 'HEXAGON')
  node.colors <- c('#85C1E9', '#F7DC6F', '#ff00ff', '#ff00bf', '#1D8348', '#ff4000', '#40ff00', '#0040ff', '#ff00ff', '#ff00bf', '#1a75ff', '	#003d99')
  setNodeShapeRule(cw, 'nodeType', attribute.values, node.shapes)
  setNodeColorRule(cw, 'nodeType', attribute.values, node.colors, mode = 'lookup')
  ##Draw edges
  g <- cw@graph
  g <- initEdgeAttribute(graph = g, attribute.name = 'edgeType', attribute.type = 'char', default.value = 'unspecified')
  for (i in genes) {
    for(j in clusters) {
      if((i %in% o$corem_list$genes[[j]]) == TRUE) {
        g <- graph::addEdge(from = j, to = i, graph = g)}
    } 
  }
  ##output graph to cytoscape
  cw@graph <- g
  print(noa.names(cw@graph))
  displayGraph(cw)
  #layoutNetwork(cw, layout.name = getLayoutNames(cw)[1])
  #redraw(cw)
}
