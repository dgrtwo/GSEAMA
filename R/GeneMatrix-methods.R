# methods for working with a GeneMatrix

setMethod("dim", "GeneMatrix", function(x) dim(x@matrix))
setMethod("nrow", "GeneMatrix", function(x) nrow(x@matrix))
setMethod("ncol", "GeneMatrix", function(x) ncol(x@matrix))

setMethod("rownames", "GeneMatrix", function(x) rownames(x@matrix))
setMethod("colnames", "GeneMatrix", function(x) colnames(x@matrix))

setMethod("as.matrix", "GeneMatrix", function(x, ...) as.matrix(x@matrix))

setMethod("print", signature(x="GeneMatrix"), function(x) {
    cat(paste("Gene Matrix with", nrow(x), "genes and", ncol(x), "columns"))
})

setMethod("[", c("GeneMatrix", "ANY", "ANY", "ANY"),
          function(x, i, j, ..., drop=TRUE)
          {
              if (!is.null(x@fit)) {
                  warning("Subsetting a model that has already been tested")        
              }
              
              if (missing(i)) {
                  i = 1:nrow(x@matrix)
              }
              if (missing(j)) {
                  j = 1:ncol(x@matrix)
              }
              
              mat = x@matrix[i, j]
              if (length(i) == 1) {
                  mat = Matrix(mat, nrow=1)
              }
              else if (length(j) == 1) {
                  mat = Matrix(mat, ncol=1)
              }
              x@matrix = mat
              
              if (!is.character(i)) {
                  x@geneData = x@geneData[i, ]
              } else {
                  x@geneData = x@geneData[match(i, ID), ]
              }
              if (!is.character(j)) {
                  x@colData = x@colData[j, ]
              } else {
                  x@colData = x@colData[match(j, ID), ]
              }    
              
              x@geneData$Count = rowSums(x@matrix != 0)
              x@colData$Count = colSums(x@matrix != 0)
              stopifnot(all(rownames(x@matrix) == x@geneData$ID))
              stopifnot(all(colnames(x@matrix) == x@colData$ID))
              x
          })
