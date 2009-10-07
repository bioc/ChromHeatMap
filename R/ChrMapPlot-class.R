## $Id: ChrMapPlot-class.R 2072 2009-10-05 16:00:38Z tfrayner $

setClass('ChrMapPlot',
         representation=list(
           labels='character',
           start='numeric',
           end='numeric')
         )

setMethod('show',
          signature=signature(object='ChrMapPlot'),
          function(object) {
              cat(class( object ))
              cat('\n')
              cat('Number of features plotted:', length(object@labels), '\n')
          })


