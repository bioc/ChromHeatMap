# $Id: grabChrMapProbes.R 2063 2009-10-02 14:37:48Z tfrayner $

## Basic idea ripped off from the idiogram::idiograb function. Removed some
## things that aren't really necessary on a one-dimensional chromosome
## plot (i.e. all the y values are ignored). This means that pos/neg
## strand plots don't support grabbing within a panel (all hits from
## both panels are returned), this may be a problem at some point
## FIXME.
grabChrMapProbes <- function (plotmap) {

    stopifnot ( inherits(plotmap, 'ChrMapPlot') )

    ## Save our par options, set them to be restored on
    ## exit. N.B. this screws with the layout such that repeated
    ## grabChrMapProbes on the same heatmap will fail. So this is
    ## deactivated for now.
    
#    op <- par(no.readonly = TRUE)
#    on.exit(par(op))

    ## Note that this function assumes that the layout() and par()
    ## status has been left in a state where the plotmap values still
    ## make sense.
    coords <- as.numeric(names(plotmap@labels))

    cat("Please click on two points to define a chromosome region.\n")

    ## Set this to type='o' once we figure out how to reset the
    ## layout() after the heatmap plot. At the moment this works fine,
    ## but plotting the selected points/lines is subject to clipping
    ## and other artifacts.
    pick <- locator(2)

    pick <- sort( pick$x )

    between <- coords > pick[1] & coords < pick[2]

    a <- plotmap@labels[between]
    return(as.character(a[!is.na(a)]))
}
