# methods for working with a GeneMatrix

setMethod("dim", "GeneMatrix", function(x) dim(x@matrix))
setMethod("nrow", "GeneMatrix", function(x) nrow(x@matrix))
setMethod("ncol", "GeneMatrix", function(x) ncol(x@matrix))
setMethod("rownames", "GeneMatrix", function(x) rownames(x@matrix))
setMethod("colnames", "GeneMatrix", function(x) colnames(x@matrix))
setMethod("as.matrix", "GeneMatrix", function(x, ...) as.matrix(x@matrix))


setMethod("print", signature(x = "GeneMatrix"), function(x) {
  cat(paste("Gene Matrix with", nrow(x), "genes and", ncol(x), "columns"))
})


setMethod("show", signature(object = "GeneMatrix"), function(object) {
  print(object)
})


setMethod("[", c("GeneMatrix", "ANY", "ANY", "ANY"),
          function(x, i, j, ..., drop = TRUE)
          {
            if (!is.null(x@fit)) {
              warning("Subsetting a model that has already been tested")        
            }
            
            # subset the matrix
            if (missing(i)) {
              i <- seq_len(nrow(x@matrix))
            }
            if (missing(j)) {
              j <- seq_len(ncol(x@matrix))
            }
            
            mat <- x@matrix[i, j]
            if (length(i) == 1) {
              mat <- Matrix(mat, nrow = 1)
            } else if (length(j) == 1) {
              mat <- Matrix(mat, ncol = 1)
            }
            x@matrix <- mat
            
            # subset gene and column data
            if (!is.character(i)) {
              x@geneData = x@geneData[i, ]
            } else {
              x@geneData = x@geneData[match(i, x@geneData$ID), ]
            }
            
            if (!is.character(j)) {
              x@colData = x@colData[j, ]
            } else {
              x@colData = x@colData[match(j, x@colData$ID)]
            }    
            
            x@geneData$Size = rowSums(x@matrix != 0)
            x@colData$Size = colSums(x@matrix != 0)
            stopifnot(all(rownames(x@matrix) == x@geneData$ID))
            stopifnot(all(colnames(x@matrix) == x@colData$ID))
            x
          })
