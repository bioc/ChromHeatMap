# $Id: paintCytoBands.R 2075 2009-10-06 09:55:10Z tfrayner $

## Taken from idiogram and quantsmooth packages, and rejigged to allow
## subsetting. Also added a mapping from organism name (as specified
## in the AnnotationDbi packages) to cytoband information stored
## within the package. The cytoband data is derived from the files
## available at UCSC.

paintCytobands <- function (species,
                            chrom,
                            pos = c(0, 0), 
                            width = 0.4,
                            length.out,
                            bands = "major",
                            legend = TRUE,
                            cex.leg = 0.8,
                            srt = 90,
                            bleach = 0,
                            start = NULL,
                            end = NULL) {
    
    bleacher <- function(x) {
        (x * (1 - bleach)) + bleach
    }

    cytobands <- NULL # avoiding the "no visible binding" warning
    stains    <- NULL
    rm(cytobands, stains)
    data('cytobands')

    chrom.bands <- cytobands[[species]]
    if ( is.null(chrom.bands) ) {
        warning(paste("Warning: Species", species, "unknown. Cannot plot idiogram."))
        return()
    }    

    chrom <- switch(as.character(chrom), "98" = "X", "99" = "Y", 
                    as.character(chrom))
    chrom <- paste('chr', chrom, sep='')

    if (length(pos) == 1) 
        pos <- c(0, pos)

    chromdata <- subset(chrom.bands, chrom.bands$chr == chrom)

    if (nrow(chromdata) > 0) {
        lc <- nchar(chromdata$band)
        sel <- !(substr(chromdata$band, lc, lc) %in% letters)
        if (bands != "major") 
            sel <- !sel
        chromdata <- chromdata[sel, ]
        rm(lc, sel)
        bandpos <- chromdata[, c("start", "end")]
        type.b <- match(chromdata$stain, stains$type)
        if ( any(isNA(type.b)) )
            warning("Warning: Some cytoband stains not recognised.")
        bandcol <- gray(bleacher(stains$bandcol))[type.b]
        textcol <- gray(bleacher(stains$textcol))[type.b]
        banddens <- stains$banddens[type.b]
        bandbord <- gray(bleacher(stains$bandbord))[type.b]
        n <- nrow(chromdata)
        centromere <- which(chromdata$arm[-n] != chromdata$arm[-1])
        if ( length(centromere) > 0 )
            idx <- c(2:(centromere - 1), (centromere + 2):(n-1))
        else
            idx <- c(2:(n-1))

        ## Calculate the idiogram segment to draw
        plotmax <- max(bandpos)
        if ( !is.null(start) && !is.null(end) ) {
            if ( end < start )
                stop("Error: end position must be greater than start position.")
            bandpos <- bandpos - start
            bandpos <- (bandpos/(end-start)) * plotmax
        }

        if (!missing(length.out)) {
            bandpos <- (bandpos/plotmax) * length.out
        }

        ## Start the drawing. Here we draw the main bands for each arm:
        rect(pos[1] + bandpos[idx, 1],
             pos[2],
             pos[1] + bandpos[idx, 2],
             pos[2] - width,
             col = bandcol[idx],
             density = banddens[idx], 
             border = bandbord[idx])

        ## Chromosome ends.
        qs.semicircle(pos[1] + bandpos[1, 2],
                      pos[2] - width, 
                      width,
                      bandpos[1, 2] - bandpos[1, 1],
                      2,
                      col = bandcol[1], 
                      density = banddens[1],
                      border = bandbord[1])
        qs.semicircle(pos[1] + bandpos[n, 1],
                      pos[2] - width, 
                      width,
                      bandpos[n, 2] - bandpos[n, 1],
                      4,
                      col = bandcol[n], 
                      density = banddens[n],
                      border = bandbord[n])

        ## If we want a centromere, draw it here.
        if ( length(centromere) > 0 ) {
            if ( centromere > idx[1] & centromere < idx[length(idx)] ) {
                qs.semicircle(pos[1] + bandpos[centromere, 1],
                              pos[2] - width,
                              width,
                              bandpos[centromere, 2] - bandpos[centromere, 1],
                              4,
                              col = bandcol[centromere],
                              density = banddens[centromere], 
                              border = bandbord[centromere])
                qs.semicircle(pos[1] + bandpos[centromere + 1, 2], 
                              pos[2] - width,
                              width,
                              bandpos[centromere + 1, 2] - bandpos[centromere + 1, 1],
                              2,
                              col = bandcol[centromere + 1],
                              density = banddens[centromere + 1],
                              border = bandbord[centromere + 1])
                centromere.size = 0.6 * 0.5 * width/yinch(1)
                symbols(pos[1] + bandpos[centromere, 2],
                        pos[2] - 0.5 * width,
                        circles = 1,
                        inches = centromere.size, 
                        add = TRUE,
                        fg = gray(bleacher(0)),
                        bg = "white")
            }
        }
        if (legend) {

            ## A little fudging to get the text position right.
            xpos <- pos[2]
            if ( srt < 60 )
                xpos <- xpos * 0.3
            else
                if ( srt < 120 )
                    xpos <- xpos * 0.4
                else
                    xpos <- xpos * 0.5

            ## Label the cytobands.
            text(pos[1] + (bandpos[, 1] + bandpos[, 2])/2,
                 xpos,
                 paste(chromdata[, "arm"], 
                       chromdata[, "band"],
                       sep = ""),
                 pos = 3,
                 srt = srt,
                 col = textcol,
                 cex = cex.leg)
        }
    }
    else {
        warning(paste("Chromosome", chrom, "is not plotted because cytoband data is not available"))
    }
}

## Taken verbatim from quantsmooth package.
qs.semicircle <- function(base.x, base.y, base.length, height=base.length, side=1, orientation=NULL,plottype="poly",...) {

    ## based on lodplot package
    ## - col is now propagated through ..., other plotting parameters can now also be given
    ## - different types poly/line 
    radius<-base.length/2
    x<-radius*seq(-1,1,length=40)
    y<-height/radius*sqrt(radius^2-x^2)
    if (is.null(orientation)) {
        co<-as.integer(cos(pi*(3-side)/2))
        so<-as.integer(sin(pi*(3-side)/2))
    }
    else {
        co<-cos(orientation)
        so<-sin(orientation)
    }
    tx<-co*x - so*y
    ty<-so*x + co*y
    if (is.null(orientation)) {
        if (side==1 || side==3) {
            base.x<-base.x+radius
        }
        else if (side==2 || side==4) {
            base.y<-base.y+radius
        }
    }
    x<-base.x+tx
    y<-base.y+ty
    switch(plottype,
           poly=polygon(x,y,...),
           line=lines(x,y,...))
}

