# $Id: chrHeatMap.R 2075 2009-10-06 09:55:10Z tfrayner $

chrHeatMap <- function (strand.data,
                        cytopaint.func = NULL,
                        col = "heat.colors",
                        start,
                        end,
                        breaks,
                        RowSideColors,
                        title=TRUE,
                        margins = c(6, 6),
                        cexCyto = 0.8,
                        srtCyto = 90,
                        lmat = NULL,
                        lhei = NULL,
                        lwid = NULL,
                        ...) {

    ## Allow the user to pass just a single matrix devoid of a list wrapper.
    if ( inherits(strand.data, 'ChrStrandMatrix' ) )
        strand.data <- list(strand.data)

    if ( length(strand.data) > 2 || length(strand.data) < 1 )
        stop("Error: Attempting to plot an unsupported number of strands.")

    if ( ! all(unlist(lapply( strand.data, class )) == 'ChrStrandMatrix') )
        stop("Error: data must be formatted as output from createChrMatrix.")

    if ( missing( start ) )
        start <- min( unlist(lapply(strand.data, function(x) { min(x@start) } ) ) )

    if ( missing( end ) )
        end   <- max( unlist(lapply(strand.data, function(x) { max(x@end) } ) ) )

    if ( ! is.null(title) & is.character(title) )
        if ( length(title) != length(strand.data) )
            stop("User-supplied title list must be of the same length as the strand list.")

    ## We want to plot the first specified strand at the top of the
    ## figure, but that heatmap is plotted last.
    strand.data <- rev(strand.data)
    title <- rev(title)

    ## Layout panel sizes. These will probably end up being hard-coded.
    if ( length(strand.data) == 2 ) {

        ## Plotting both strands.
        if (missing(lhei) || is.null(lhei)) 
            lhei <- c(1, 0.1, 1)
        if (missing(lwid) || is.null(lwid)) 
            lwid <- c(0.1, 1)
        if (missing(lmat) || is.null(lmat)) {
            lmat <- rbind(4:5, c(NA, 1), 2:3)
            if (!missing(RowSideColors)) {
                lmat <- rbind(c(5:7), c(NA, NA, 1), c(2:4) )
                lwid <- c(lwid[1], 0.02, lwid[2])
            }
            lmat[is.na(lmat)] <- 0
        }
    }
    else {

        ## Plotting just one strand.
        if (missing(lhei) || is.null(lhei)) 
            lhei <- c(0.1, 1)
        if (missing(lwid) || is.null(lwid)) 
            lwid <- c(0.1, 1)
        if (missing(lmat) || is.null(lmat)) {
            lmat <- rbind(c(NA, 1), 2:3)
            if (!missing(RowSideColors)) {
                lmat <- rbind(c(NA, NA, 1), c(2:4) )
                lwid <- c(lwid[1], 0.02, lwid[2])
            }
            lmat[is.na(lmat)] <- 0
        }
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) = ", ncol(lmat))

#    on.exit(layout(pmat, widths = lwid, heights = sum(lhei), respect = FALSE))
#    on.exit(par( c(margins[1], 0, 1, margins[2]) ), add=TRUE)

    ## Create a gridded layout for output.
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

    ## Draw the idiogram.
    par( mar=c(0,0,0,margins[2]))

    if ( ! is.null( cytopaint.func ) ) {

        boxwidth<-ncol(strand.data[[1]]@data)
        plot(c(1, boxwidth), c(1, 10), type='n', axes=FALSE, xaxs='i')

        cytopaint.func( boxwidth, cexCyto, srtCyto )
    }
    else {
        plot.new()
    }
    
    ## Calculate colors, in case they've not been supplied.
    col <- calculateColors( col, breaks )
    
    ## Calculate default color breakpoints. We space the color breaks according to the data quantiles.
    if ( missing(breaks) ) {
        i <- seq(1:(length(col)+1))
        i <- i/i[length(i)]
        alldata <- unlist(lapply(strand.data, function(d) { d@data } ))
        if ( any ( !is.na(alldata) ) )
            breaks <- quantile(alldata, probs=i, na.rm=TRUE)
        else
            breaks <- i
    }

    ## Do the actual drawing.
    for ( n in 1:length(strand.data) ) {
        strand <- strand.data[[n]]

        ## Figure out what we want for heatmap titles.
        main <- NULL
        if ( isTRUE(title) )
            main <- sprintf('chromosome %s, %s %s',
                            chrNames(strand),
                            strandName(strand),
                            ifelse(strandName(strand) == 'both', 'strands','strand'))
        else
            if ( ! is.null(title) & is.character(title) )
                main <- title[n]

        drawMapDendro(strand,
                      RowSideColors = RowSideColors,
                      margins = margins,
                      col = col,
                      breaks = breaks,
                      start = start,
                      end = end,
                      main = main,
                      ... );
    }
    
    return()
}

