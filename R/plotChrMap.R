### $Id: plotChrMap.R 2075 2009-10-06 09:55:10Z tfrayner $

calculateCytoband <- function(species,
                              chr,
                              cytoband) {

    cytobands <- NULL # avoiding the "no visible binding" warning
    rm(cytobands)
    data('cytobands')

    chrom.bands <- cytobands[[species]]
    if ( is.null(chrom.bands) ) {
        stop(paste("Error: Species", species, "unknown. Cannot calculate chromosome region to plot."))
    }    
    chr <- paste('chr', chr, sep='')
  
    ## Parse arm, band from cytoband
    arm <- tolower(substr(cytoband, 1, 1))
    stopifnot( arm %in% c('p','q') )
    band <- substr(cytoband, 2, 1000000)

    arm.data  <- chrom.bands[chrom.bands$chr == chr & chrom.bands$arm == arm,]
    pattern <- paste( '^', band, '\\b', sep='' )
    band.data <- arm.data[grep( pattern, arm.data$band, perl=TRUE),]

    if ( nrow(band.data) == 0 )
        stop("Error: no coordinates found for the specified band.")

    start <- min(band.data$start)
    end   <- max(band.data$end)

    return(c(start, end))
}

calculateLimits <- function( data, chr, start=1, end, cytoband ) {

    if ( ! inherits(data, 'ChrStrandData') )
        stop("Error: data must be a ChrStrandData object, e.g. the output of makeChrStrandData()")

    require(annotate)    ## imports getAnnMap

    ## Load the specified library (add the .db extension).
    lib <- data@lib
    if (is.null(lib))
        stop("Error: data has no lib attribute.")
    liblen <- nchar(lib)
    if ( substr(lib, liblen-2, liblen) == '.db' ) {
        require( lib, character.only=TRUE )
    }
    else {
        require( paste(lib, '.db', sep=''), character.only=TRUE )
    }

    chrlength.env <- getAnnMap('CHRLENGTHS', lib)
    species <- getAnnMap('ORGANISM', lib)
    chr <- as.character(chr)
    maxchr <- chrlength.env[[chr]]
    if (missing(end))
        end <- maxchr

    range <- c(start, end)
    if ( ! missing(cytoband) ) {
        range <- calculateCytoband(species = species,
                                   chr = chr,
                                   cytoband = cytoband)
    }

    return(list(range, species))
}

plotChrMap <- function(data,
                       chr,
                       start = 1,
                       end,
                       subset = NULL,
                       cytoband,
                       interval = ceiling((end-start)/500),
                       strands = c('forward', 'reverse'),
                       ...) {

    if ( ! inherits(data, 'ChrStrandData') )
        stop("Error: data must be a ChrStrandData object, e.g. the output of makeChrStrandData()")

    if ( ! (length(strands) <= 2 && length(strands) >= 1) )
        stop("Error: Only one or two strands can currently be plotted.")

    lims <- calculateLimits( data, chr, start, end, cytoband )

    start <- lims[[1]][1]
    end   <- lims[[1]][2]
    species <- lims[[2]]

    message(sprintf("Using interval argument set to %d", interval))

    ## Generate the matrices we need to plot.
    matrices <- lapply(strands,
                       function(x) {
                           createChrMatrix(data = data,
                                           chr = chr,
                                           start = start,
                                           end = end,
                                           interval = interval,
                                           subset = subset,
                                           strand = x )
                       })

    ## Create a closure around our metadata to be used later to plot
    ## the idiogram.
    cytopaint <- function (boxwidth, cexCyto, srt) {
        paintCytobands(species,
                       chr,
                       pos=c(0,9),
                       length.out=boxwidth,
                       width=7,
                       cex.leg=cexCyto,
                       srt=srt,
                       start=start,
                       end=end)
    }

    chrHeatMap(strand.data    = matrices,
               cytopaint.func = cytopaint,
               start = start,
               end   = end,
               ...)

    ## Return values mapping plotted coordinate to probe IDs.
    ids <- unlist(lapply(matrices, function(m) { m@probeID }))
    ids <- ids[ order(as.numeric(names(ids))) ]

    plotmap <- new('ChrMapPlot', labels = ids, start = 1, end = (end-start)/interval)
    
    return(plotmap)
}

