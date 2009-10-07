## $Id: ChrStrandData-class.R 2072 2009-10-05 16:00:38Z tfrayner $

setClass('ChrStrandData', representation=list(data='list', lib='character', chrs='character'))

setMethod('show',
          signature=signature(object='ChrStrandData'),
          function(object) {
              cat(class( object ))
              cat('\n')
              cat('Annotation: ', annotation(object),                       '\n')
              cat('Samples:    ', Biobase::selectSome(sampleNames(object)), '\n')
              cat('Chromosomes:', Biobase::selectSome(chrNames(object)),    '\n')
          })

setMethod('summary',
          signature=signature(object='ChrStrandData'),
          function(object) {
              cat(class( object ))
              cat('\n')

              chrs <- list()
              for ( n in 1:length(object@data) ) {
                  sampchrs <- object@data[[n]]
                  for ( m in 1:length(sampchrs) ) {
                      chr <- names(sampchrs)[m]
                      chrs[[chr]] <- union(chrs[[chr]],
                                           unlist(lapply(sampchrs[[m]],
                                                         function(p) names(p$x)) ))
                  }
              }
              chrs <- data.frame(total=unlist(lapply(chrs, length)))

              cur.warn <- options('warn'=-1)   ## warnings off briefly
              chrs <- chrs[order(as.numeric(rownames(chrs))),, drop=FALSE]
              options('warn'=cur.warn$warn)

              cat('Features per chromosome:\n')
              print(chrs)
          })

setMethod('annotation', signature=signature(object='ChrStrandData'),
          function(object) object@lib)

setGeneric('chrNames', def=function(object) standardGeneric('chrNames'))

setMethod('chrNames', signature=signature(object='ChrStrandData'),
          function(object) object@chrs)

setMethod('sampleNames', signature=signature(object='ChrStrandData'),
          function(object) names(object@data))
