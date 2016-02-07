##################################
######## AllClasses.R
########
########    all classes in GSEAMA
##################################
####################################################################

setClass("GeneMatrix",
    representation(
        matrix = "Matrix",
        colData = "ANY",
        geneData = "ANY",
        fit = "ANY",
        rankingMetric = "character",
        effectMetric = "character"
    )
)
