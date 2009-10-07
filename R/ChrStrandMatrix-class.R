## $Id: ChrStrandMatrix-class.R 2073 2009-10-05 17:54:40Z tfrayner $

setClass('ChrStrandMatrix',
         representation=list(
           data='matrix',
           probeID='array',
           chr='character',
           strand='character',
           start='array',
           end='array')
         )

setMethod('show',
          signature=signature(object='ChrStrandMatrix'),
          function(object) {
              cat(class( object ))
              cat('\n')
              cat('Dimensions:         ', nrow(object@data), 'samples,',
                                          ncol(object@data), 'features\n')
              cat('Chromosome:         ', chrNames(object),  '\n')
              cat('Strand:             ', strandName(object),'\n')
              cat('Starting coordinate:', min(object@start), '\n')
              cat('Ending coordinate:  ', max(object@end),   '\n')
          })

setMethod('summary',
          signature=signature(object='ChrStrandMatrix'),
          function(object) {
              cat(class( object ), 'with', (nrow(object@data)), 'samples:\n\n')
              print(summary(t(object@data)))
          })


setMethod('featureNames', signature=signature(object='ChrStrandMatrix'),
          function(object) object@probeID)

setMethod('chrNames', signature=signature(object='ChrStrandMatrix'),
          function(object) object@chr)

setGeneric('strandName', def=function(object) standardGeneric('strandName'))

setMethod('strandName', signature=signature(object='ChrStrandMatrix'),
          function(object) object@strand)

setMethod('sampleNames', signature=signature(object='ChrStrandMatrix'),
          function(object) rownames(object@data))

setMethod('exprs', signature=signature(object='ChrStrandMatrix'),
          function(object) {
            x <- t(object@data)
            rownames(x) <- object@probeID
            return(x)
          })
