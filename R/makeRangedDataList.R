## $Id: makeRangedDataList.R 2068 2009-10-05 11:01:44Z tfrayner $

makeRangedDataList <- function (data, chr, start=1, end, genome, subset=NULL,
                                cytoband, plot=FALSE, session) {

    .makeGRanges <- function( sample, mat, chr, targetRanges, genome ) {
        score <- mat@data[sample,]
        as(GenomicData(targetRanges,
                       score,
                       chrom = chr,
                       genome = genome), 'GRanges')
    }
                    
    if ( ! inherits(data, 'ChrStrandData') )
        stop("Error: data must be a ChrStrandData object, e.g. the output of makeChrStrandData()")

    lims <- calculateLimits( data, chr, start, end, cytoband )

    start <- lims[[1]][1]
    end   <- lims[[1]][2]
    species <- lims[[2]]

    ## N.B. the interval arg is probably not needed here, given the
    ## width of the average UCSC browser plot.
    mat <- createChrMatrix(data = data,
                           chr = chr,
                           start = start,
                           end = end,
                           subset = subset,
                           strand = 'both' )

    starts <- mat@start
    startord <- order(starts)
    overlap <- starts[startord][-1] - mat@end[startord][-length(starts)] < 1
    if ( any(overlap) ) {
        val.st <- mat@start[startord[-1][overlap]]
        val.en <- mat@end[startord[-length(starts)][overlap]]
        val.mid <- val.st+((val.st - val.en) / 2)
        
        mat@start[startord[-1][overlap]] <- val.mid
        mat@end[startord[-length(starts)][overlap]] <- val.mid-1
    }

    targetRanges <- IRanges(mat@start, mat@end)

    chr <- paste('chr', chr, sep='')
    tracks <- lapply( names(data@data), .makeGRanges, mat, chr, targetRanges, genome )
    names( tracks ) <- names(data@data)
    tracks <- do.call(GRangesList, tracks)

    min <- min(mat@data) - (0.1 * abs(min(mat@data)))
    max <- max(mat@data) + (0.1 * abs(max(mat@data)))

    if ( plot == TRUE ) {
        if ( missing(session) )
            session<-browserSession()

        for ( sample in names(tracks) ) {
            track(session=session, trackName=sample, maxHeightPixels=12,
                  autoScale=FALSE, viewLimits=c(min, max)) <- tracks[[sample]]
        }

        view <- browserView(session, range(tracks[[1]]), dense=names(tracks))
    }

    return(tracks)
}
