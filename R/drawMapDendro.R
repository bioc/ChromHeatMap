# $Id: drawMapDendro.R 2075 2009-10-06 09:55:10Z tfrayner $

calculateColors <- function( col, breaks ) {

    if ( ! missing(breaks) )
        number <- length(breaks) - 1
    else
        number <- 100
    
    if (class(col) == "function") 
        col <- col(number)
    else if (is.character(col) && length(col) == 1) 
        col <- do.call(col, list(number))
    return(col)
}

drawMapDendro <- function (x,
                           start,
                           end,
                           col = "heat.colors", 
                           dendrogram = TRUE,
                           Rowv = TRUE,
                           margins = c(6, 6), 
                           na.rm = TRUE,
                           hclustfun = hclust,
                           distfun = dist,
                           breaks,
                           RowSideColors,
                           cexRow,
                           cexCol,
                           xlab, 
                           ylab,
                           labRow,
                           labCol,
                           na.color = 'gray',
                           ...) {

    ## FUNCTIONS.
    invalid <- function (x) {
        if (missing(x) || is.null(x) || length(x) == 0) 
            return(TRUE)
        if (is.list(x)) 
            return(all(sapply(x, invalid)))
        else if (is.vector(x)) 
            return(all(is.na(x)))
        else return(FALSE)
    }

    calculateDendro <- function( Rowv, x, na.rm = TRUE, hclustfun = hclust, distfun = dist ) {
        if (inherits(Rowv, "dendrogram")) {
            ddr <- Rowv
        }
        else if (is.integer(Rowv)) {
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            ddr <- reorder(ddr, Rowv)
        }
        else if (isTRUE(Rowv)) {
            Rowv <- rowMeans(x, na.rm = na.rm)
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            ddr <- reorder(ddr, Rowv)
        }
        else {
            ddr <- NULL
        }

        return(ddr);
    }

    stopifnot( inherits(x, 'ChrStrandMatrix') )

    if ( missing(start) )
        start <- min(x@start)
    if ( missing(end) )
        end <- max(x@end)
    
    ## Check and sanitize our dendrogram arguments. Only row dendrograms supported. 
    if (!inherits(Rowv, "dendrogram")) {
        if ((!isTRUE(Rowv) || is.null(Rowv)) && isTRUE(dendrogram)) {
            dendrogram <- FALSE
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                    dendrogram, "'. Omitting row dendogram.")
        }
    }

    data <- x@data

    ## Get nr, nc; we'll be using them below.
    if (length(di <- dim(data)) != 2 || !is.numeric(data)) 
        stop("The data slot of `x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]

    ## Basic argument checking on our matrix.
    if (nr < 1 || nc < 1) 
        stop("The data slot of `x' must have at least 1 row and 1 column")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    
    ## Calculate the dendrogram if needed.
    if ( any( ! is.na(data) ) ) {
        ddr <- calculateDendro(Rowv,
                               data,
                               na.rm = na.rm,
                               hclustfun = hclustfun,
                               distfun = distfun  )
    } else {
        ddr <- NULL
    }
    
    ## rowInd is derived from the dendrogram if possible.
    if (is.null(ddr))
        rowInd<-nr:1
    else rowInd<-order.dendrogram(ddr)
    if (nr != length(rowInd)) 
        stop("row dendrogram ordering gave index of wrong length")

    ## There will never be a column dendrogram.
    colInd <- 1:nc

    ## Reorder the rows of our matrix appropriately.
    data <- data[rowInd, , drop=FALSE]

    ## Row and column labels, where provided.
    if ( missing(labRow) ) 
        labRow <- if (is.null(rownames(data))) 
            (1:nr)[rowInd]
        else rownames(data)
    else labRow <- labRow[rowInd]
    if ( missing(labCol) )
        labCol <- as.integer(seq(start, end, by=(end-start)/20))
    else
        stopifnot( is.numeric(labCol) )

    ## Calculate colors, in case they've not been supplied.
    col <- calculateColors( col, breaks )
    
    ## Sort out our default color breaks.
    if ( missing(breaks) ) {
        i <- seq(1:(length(col)+1))
        i <- i/i[length(i)]
        if ( any( ! is.na(data) ) ) {
            breaks <- quantile(data, probs=i, na.rm=TRUE)
        } else {

            ## Dummy for all-NA data.
            breaks <- 1:(length(col)+1)
        }
    }

    ## Anything falling outwith breaks is moved to the boundary of
    ## min.breaks or max.breaks.
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    data[] <- ifelse(data < min.breaks, min.breaks, data)
    data[] <- ifelse(data > max.breaks, max.breaks, data)

    ## Top margin width.
    top.margin = 2

    ## Removed any chance of a column dendrogram; here we handle rows.
    if (isTRUE(dendrogram) && ! is.null(ddr)) {
        par(mar = c(margins[1], 0, top.margin, 0))
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()    ## Skip the dendrogram plot.

    ## Set our RowSideColors.
    if (!missing(RowSideColors)) {
        if (!is.character(RowSideColors) || length(RowSideColors) != nr) 
            stop("'RowSideColors' must be a character vector of length nrow(data)")
        par(mar = c(margins[1], 0, top.margin, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }

    ## Set our margins for the main heatmap.
    par(mar = c(margins[1], 0, top.margin, margins[2]))

    ## Get our X coordinates.
    xc <- c(start,
            unlist(lapply(1:length(x@start),
                       function(n) { c(x@start[n], x@end[n]) })),
            end)

    ## Fudge to avoid crashes on overlapping genes. Not elegant, but
    ## it should help make this more robust.
    for ( n in 2:length(xc) ) {
        if ( xc[n] <= xc[n-1] ) {
            xc[n] <- xc[n-1] <- xc[n] + (abs(xc[n]-xc[n-1])/2)
            xc[n] <- xc[n]+2  # must be *increasing*. Why +1 doesn't work here is totally beyond me.
        }
    }

    stopifnot( ! is.unsorted(xc) )

    ## Rewrite the data matrix such that gaps between plotted regions
    ## are covered with an NA.
    plotdata <- matrix(rep(NA,nr), byrow=FALSE, nrow=nr)
    for ( n in 1:ncol(data) ) {
        nextrow <- data[,n]
        nextrow <- append(nextrow, rep(NA,nr))
        plotdata <- matrix(append(plotdata, nextrow), byrow=FALSE, nrow=nr)
    }
    plotdata <- t(plotdata)

    ## Plot the actual heatmap.
    image(xc, 1:nr, plotdata, xlim = 0.5 + c(start, end), ylim = 0.5 + 
          c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
          breaks = breaks, ...)

    ## Overlay na.color on the image for any(is.na(x)).
    if (!invalid(na.color) & any(is.na(plotdata))) {
        mmat <- ifelse(is.na(plotdata), 1, NA)
        image(xc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
              col = na.color, add = TRUE)
    }

    ## Draw the x axis.
    if ( missing(cexCol) )
        cexCol <- max(c(0.4 + 1/log10(length(labCol)+1), 0.8))   # nc+1 in case of all-NA
    axis(1, labCol, labels = labCol, las = 2, line = -0.5, tick = 0, 
         cex.axis = cexCol)
    if ( ! missing(xlab) ) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)

    ## Draw the y axis.
    if ( missing(cexRow) )
        cexRow <- max(c(0.2 + 1/log10(length(labRow)+1), 0.8))
    axis(4, 1:nr, labels = labRow, las = 2, line = -0.5, tick = 0, 
         cex.axis = cexRow)
    if ( ! missing(ylab) ) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)

    return();
}

