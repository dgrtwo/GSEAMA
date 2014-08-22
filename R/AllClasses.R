##################################
######## AllClasses.R
########
########    all classes in GSEAMA
##################################
####################################################################

setClass("GeneMatrix",
    representation(
        matrix = "Matrix",
        colData = "data.table",
        geneData = "data.table",
        fit = "ANY",
        rankingMetric = "character"
    )
)
