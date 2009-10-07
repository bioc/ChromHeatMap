# $Id: createChrMatrix.R 2073 2009-10-05 17:54:40Z tfrayner $

createChrMatrix <- function(data,
                            chr,
                            strand = c('forward','reverse','both'),
                            subset = NULL,
                            start = 1,
                            end,
                            interval = ceiling((end-start)/500)) {

    getData <- function(sample, chr, allx, strands) {
        z <- lapply(strands,
                    function(x) {

                        ## N.B. if allx[[x]] is NA tapply will ignore it.
                        tapply(sample[[chr]][[x]]$y, allx[[x]], mean)
                    })
        z <- unlist(z)
        z <- z[order(as.numeric(names(z)))]
        return(z)
    }

    if ( ! inherits(data, 'ChrStrandData') )
        stop("Error: data must be a ChrStrandData object, e.g. the output of makeChrStrandData()")

    strand <- match.arg(strand)
    if ( ! (length(strand) <= 2 && length(strand) >= 1) )
        stop("Error: Required strand argument can only be one or two elements long.")

    if ( start > end )
        stop("Error: end position must be greater than start position.")
    if ( end - start <= interval )
        stop("Error: interval must be smaller than the distance between start and end.")

    sampledata <- data@data
    if ( is.null( subset ) )
        subset <- 1:length(sampledata)
    chr <- as.character(chr)

    ## Generate the matrix to be plotted.
    query.strands <- c()
    if ( strand == 'both' ) {
        query.strands <- c('posS','negS')
    }
    else {
        query.strands <- unlist(lapply(strand, switch, forward='posS', reverse='negS'))
    }

    ## Sample 1 is used as the standard for x coordinates
    allx <- lapply(query.strands, function(x) { sampledata[[1]][[chr]][[x]]$x } )
    names(allx) <- query.strands
    if ( is.null(allx) )
        stop("Error: no data for this chromosome")

    ## Here we surreptitiously add 1 to any reverse strand coords which
    ## share a coordinate with the forward strand. We have to do this
    ## because later we assume there are no collisions. Possible bug:
    ## if the colliding coordinate corresponds to the first of two
    ## adjacent forward strand coordinates (this seems unlikely though).
    if ( strand == 'both' ) {
        shared <- allx[[2]] %in% allx[[1]]
        allx[[2]][ shared ] <- allx[[2]][ shared ] + 1
    }

    ## CHRLOCENDS
    allxends <- lapply(query.strands, function(x) { sampledata[[1]][[chr]][[x]]$xe } )
    names(allxends) <- query.strands

    ## The core data matrix ("data.mat") is built here.
    sampledata <- sampledata[subset]
    data.mat <- getData(sampledata[[1]], chr, allx, query.strands)
    for (sample in sampledata[-1])
        data.mat <- rbind(data.mat, getData(sample, chr, allx, query.strands))

    ## Combine all probe IDs for each x coordinate. 
    probeIds <- lapply(query.strands,
                       function(x) {

                           ## N.B. if allx[[x]] is NA tapply will ignore it.
                           tapply(names(allx[[x]]), allx[[x]], paste, collapse=';')
                       })
    probeIds <- unlist(probeIds)
    probeIds <- probeIds[order(as.numeric(names(probeIds)))]

    ## Sort the X coordinates and remove any NAs (N.B. data.mat and
    ## probeIds have both had values linked to allx==NA removed by
    ## tapply(), so this is valid).
    allx <- sort(as.numeric(unlist(allx)))
    allxends <- sort(as.numeric(unlist(allxends)))

    ## End coordinates are divided up according to their corresponding
    ## start coords.
    endcoords <- unlist(lapply(split(allxends, allx), mean))
    coords <- as.numeric(unique(na.omit(allx)))


    if ( ! length(coords) == length(endcoords) )
        stop("Error: mismatch between CHRLOC and CHRLOCEND data.")

    if ( any( endcoords < coords ) )
        stop("Error: some CHRLOCENDs are smaller than the CHRLOCs.")

    if ( ! length(coords) == ncol(data.mat) )
        stop("Error: mismatch between CHRLOC and the data matrix dimensions.")

    ## Map the data.mat, via coords, onto a new matrix representing the
    ## actual chromosome coordinates. Data is binned according to the
    ## "interval" attribute.

    ## Define our breaks. This is designed to give an irregular
    ## distribution where genes bigger than the interval are kept
    ## intact.
    breaks <- coords[1]
    running <- 0
    for ( n in 2:length(coords) ) {

        seg.length <- coords[n] - coords[n-1]
        if ( seg.length + running < interval ) {

            ## Region not big enough yet; skip to the next gene.
            running <- running + seg.length
        } else {

            ## Region is big enough; record the following start coord.
            breaks <- append(breaks, coords[n])
            running <- 0
        }
    }
    breaks <- append(breaks, endcoords[length(endcoords)])

    ## Cut the data.mat and probeIds objects according to those breaks.
    bins <- cut(coords, breaks, right=FALSE)

    ## Filter out bins on either side of the region of interest
    ## (tapply will ignore NA bins). Only remove regions which are
    ## entirely outwith the start-end interval.
    s <- tapply(coords, bins, min)
    e <- tapply(endcoords, bins, max)
    bins[bins %in% levels(bins)[s>=end | e<=start]] <- NA
    bins <- factor(bins)

    if( all( is.na(bins) ) ) {
        new.mat <- matrix(rep(as.numeric(NA), nrow(data.mat)),
                          nrow=nrow(data.mat))
        coords <- start
        endcoords <- end
        probeIds <- NA
    } else {

        ## Calculate means for the bins (and also probeIds).
        new.mat <- matrix(t(apply(data.mat, 1, tapply, bins, mean)),
                          nrow=nrow(data.mat))

        ## Annotate with x, xe, probeIds.
        coords <- tapply(coords, bins, min)
        endcoords <- tapply(endcoords, bins, max)
        probeIds <- tapply(probeIds, bins, paste, collapse=';')
    }

    rownames(new.mat) <- names(sampledata)
    names(probeIds) <- coords

    stopifnot(! is.unsorted(coords) )
    stopifnot(! is.unsorted(endcoords) )

    ## Final run through to make sure that all bands are at least
    ## "interval" in width (otherwise very small bands drop out of the
    ## display due to pixel width limits). Note that overlapping
    ## sections will be dealt with in drawMapDendro.
    endcoords <- ifelse( endcoords - coords < interval, coords + interval, endcoords )

    ## Sanity check
    if ( ncol(new.mat) == 0 )
        stop("Error: at least one probe needed on each strand in the specified chromosome region.")
    stopifnot( length(probeIds) == ncol(new.mat) )

    ## Set the probeID, start and end attributes.
    chr.mat <- new('ChrStrandMatrix',
                   data=new.mat,
                   probeID=probeIds,
                   chr=chr,
                   strand=strand,
                   start=coords,
                   end=endcoords)
    
    return(chr.mat)
}

