getDAG <- function(listChildren=GO.db::GOBPCHILDREN){
	ch <- lapply(AnnotationDbi::as.list(listChildren), function(x) x[names(x)=="part_of" | names(x)=="is_a"])
	l <- sapply(ch, length)
	ch[l>0]
}

getGeneSet <- function(listGene, min.genes=1, filter="IEA", available.genes){
    go <- AnnotationDbi::as.list(listGene)
	# filter unwanted annotation
	if(length(filter)==1) go <- lapply(go, function(x) x[names(x) != filter])
	if(length(filter)>1) go <- lapply(go, function(x) x[!(names(x) %in% filter)])
	# keep only available.genes or unique features, will lose all annotation names
    if(!missing(available.genes)){
        go <- lapply(go, function(x) base::intersect(x, available.genes))
    }else{
    	go <- lapply(go, function(x) unique(x))
    }
	l <- sapply(go, length)

	# minimum size
	go <- go[l>=min.genes]
}

buildTree <- function(root="GO:0008150", DAG, geneSet){
	tree <- tree_transformation(DAG, geneSet, root)
	#45min
	
	nDAG = length(tree$mapping)
	nTree = length(tree$parent_index)
	parentIndex = tree$parent_index
	sum(1:nTree > parentIndex) == nTree
	# zero based
	depth <- rep(0, nTree)
	for (i in 2:nTree) {
	    depth[i] <- depth[parentIndex[i]] + 1
	}
	childrenIndexMat = tree$children_index
	childrenCount = colSums(childrenIndexMat != 0)
	componMappingMat = matrix(0, nrow = nDAG, ncol = nTree)
	###
	componMappingMat[1, 1] = 1
	for (i in 2:nDAG) {
	    componMappingMat[i, tree$mapping[[i]]] = 1
	}
	geneSet = lapply(1:nTree, function(x) tree$all_genes[tree$genes_index[, x]])
	GO.ID <- sapply(tree$clones, function(x) x[1])
	listClone <- tree$clones
	gene_count <- sapply(geneSet, length)
	propInParent = rep(1, nTree)
	propInParent[-1] = (gene_count[-1])/(gene_count[tree$parent_index[-1]])
	leaf = which(childrenCount == 0)
	#nonleaf = unlist(lapply((max(depth) - 1):0, function(x) (1:nTree)[depth == x & childrenCount > 0]))
	nonleaf = unlist(lapply(max(depth):1, function(x) unique(parentIndex[depth == x])))

	treedata <- list(nDAG = nDAG, nTree = nTree, depth = depth, 
	                 parentIndex = parentIndex, childrenCount = childrenCount, 
	                 childrenIndexMat = childrenIndexMat, componMappingMat = componMappingMat, 
	                 geneSet = geneSet, GO.ID = GO.ID, listClone = listClone, 
	                 propInParent = propInParent, leaf = leaf, nonleaf = nonleaf)
	
	#15.5min				 
	#compute the mutual ancestor with the largest index of all the components
	MA <- MAC(treedata)
	formula <- vector(mode = "list", length = nDAG)
	sepLocation <- vector(mode = "list", length = nDAG)
	for (i in 2:nDAG) {
	    componList <- tree$mapping[[i]]
	    while (length(componList) > 1) {
	        componMatrix <- MA[componList, componList]
	        temp <- unique(setdiff(as.vector(componMatrix), -1))
	        maxValue <- max(temp[depth[temp] == max(depth[temp])])
	        maxIndexes <- union(which(componMatrix == maxValue, 
	                                  arr.ind = TRUE), c())
	        formula[[i]] <- c(formula[[i]], c(maxValue, sort(componList[maxIndexes])))
	        sepLocation[[i]] <- c(sepLocation[[i]], length(formula[[i]]))
	        componList <- union(componList[-maxIndexes], maxValue)
	    }
	}
	treedata$formula <- formula
	treedata$sepLocation <- sepLocation
	
	return(treedata)
}
