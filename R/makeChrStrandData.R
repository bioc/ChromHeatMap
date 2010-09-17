## $Id: makeChrStrandData.R 2068 2009-10-05 11:01:44Z tfrayner $

retrieveAnnot <- function ( genes, envir ) {

    require('AnnotationDbi')

    ## We break the query into chunks so that it doesn't overwhelm SQLite.
    chunk   <- 100000
    results <- list()
    iters   <- floor((length(genes)-1)/chunk)
    for ( i in 0:iters ) {
        start <- ( i * chunk ) + 1
        end   <- min( ( i + 1 ) * chunk, length(genes) )
        results <- append(results, AnnotationDbi::mget(genes[start:end], envir, ifnotfound=NA))
    }

    return(results)
}

setGeneric('makeChrStrandData', def=function(expr, lib) standardGeneric('makeChrStrandData'))

require(Biobase)
setMethod('makeChrStrandData', signature(expr='ExpressionSet'), function(expr , lib ) makeChrStrandData(exprs(expr), lib))

## The below was cut and pasted from the geneplotter Makesense code in
## svn, and then modified to remove the call to lowess and add support for CHRLOCEND.
setMethod('makeChrStrandData', signature(expr='matrix'), function(expr, lib){

    require('AnnotationDbi')
    require('annotate')

    if ( missing(expr) )
      stop("Error: expr argument is required")

    if ( ! inherits(expr, 'matrix') )
      stop("Error: expr must be a matrix")

    if ( missing(lib) )
      stop("Error: lib argument is required")

    ## Load the specified library (add the .db extension).
    liblen <- nchar(lib)
    if ( substr(lib, liblen-2, liblen) == '.db' ) {
        require( lib, character.only=TRUE )
    }
    else {
        require( paste(lib, '.db', sep=''), character.only=TRUE )
    }

    if (length(lib) != 1 || nchar(lib) < 1)
        stop("'lib' argument must be length one")

    genes <- rownames(expr)
    if ( length(genes) != nrow(expr) )
        stop("Error: The data matrix must have row names corresponding to the annotation package probe or gene IDs.")

    libCHR <- getAnnMap("CHR", lib)
    libCHRLOC <- getAnnMap("CHRLOC", lib)
    libCHRLOCEND <- getAnnMap("CHRLOCEND", lib)

    ## Select genes that are annotated at exactly _one_ chromosome.
    chr <- retrieveAnnot(genes, libCHR)

    if( all(is.na(chr)) )
        stop("No annotation returned from the library. Are you using the right one?")

    oneC <- sapply(chr, function(x)
                   if (length(x) == 1 && !is.na(x)) TRUE else FALSE)

    ## Select genes that are annotated at exactly _one_ chrom location
    chrL <- retrieveAnnot(genes, libCHRLOC)
    chrLE <- retrieveAnnot(genes, libCHRLOCEND)

    if( all(is.na(chrL)) )
        stop("No annotation returned from the library. Are you using the right one?")

    ## FIXME we should probably check that length(na.rm(chrL)) and
    ## length(na.rm(chrLE)) are the same etc. Or add a "oneLE" requirement to wanted, below.

    oneL <- sapply(chrL, function(x)
                   if (length(x) == 1 && !is.na(x)) TRUE else FALSE)
    wanted <- oneC & oneL
    chrName <- unlist(chr[wanted])
    chrPos <- unlist(chrL[wanted])
    chrPosEnd <- unlist(chrLE[wanted])

    cP <- split(chrPos, chrName)
    cPE <- split(chrPosEnd, chrName)

    gE <- expr[wanted, ]
    ans2 <- vector("list", length=ncol(gE))

    for( j in 1:ncol(gE) ) {
        s1 <- split(gE[,j], chrName)
        ans <- NULL
        for (i in names(cP)) {
            d1 <- s1[[i]]
            cL <- cP[[i]]
            cLE <- cPE[[i]]

            dp <- d1[cL>0]
            lp <- cL[cL>0]
            lpe <- cLE[cL>0]  # deliberate. ensures length(lp) == length(lpe)

            dn <- d1[cL<0]
            ln <- -cL[cL<0]
            lne <- -cLE[cL<0]  # deliberate. ensures length(ln) == length(lne)

            if (length(lp)) {

                lw1 <- list()
                lw1$x <- sort(lp)
                lw1$xe <- lpe[order(lp)]
                names(lw1$x) <- names(dp)[order(lp)]
                names(lw1$xe) <- names(dp)[order(lp)]
                lw1$y <- dp[order(lp)]
            }
            else {
                lw1 <- list(x=numeric(0), y=numeric(0), xe=numeric(0))
            }
            if (length(ln)) {
                lw2 <- list()
                lw2$x <- sort(ln)
                lw2$xe <- lne[order(ln)]
                names(lw2$x) <- names(dn)[order(ln)]
                names(lw2$xe) <- names(dn)[order(ln)]
                lw2$y <- dn[order(ln)]
            }
            else {
                lw2 <- list(x=numeric(0), y=numeric(0), xe=numeric(0))
            }
            ans[[i]] <- list(posS = lw1, negS =lw2)
        }
        ans2[[j]] <- ans
    }

    ## Copy over the data matrix row names.
    names(ans2) <- colnames( expr )

    ## Generate a list of chromosomes for later inspection.
    chrs <- Reduce(union, lapply(ans2, names))

    data <- new('ChrStrandData', data=ans2, lib=lib, chrs=chrs)

    return(data)
})

