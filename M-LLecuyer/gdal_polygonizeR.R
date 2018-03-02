



gdal_polygonizeR <- function(x, 
                             outshape   = NULL, 
                             gdalformat = 'ESRI Shapefile',
                             pypath     = NULL, 
                             readpoly   = TRUE, 
                             quiet      = TRUE, 
                             overwrite  = FALSE) {
  
  if (is.null(pypath)) {
    pypath <- Sys.which('gdal_polygonize.py')
  }
  
  if (!file.exists(pypath)) stop("Can't find gdal_polygonize.py on your system.")
  owd <- getwd()
  on.exit(setwd(owd))
  setwd(dirname(pypath))
  if (!is.null(outshape)) {
    outshape <- sub('\\.shp$', '', outshape)
    f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep = '.'))
    if (any(f.exists)) {
      if (overwrite == FALSE) {
        stop(sprintf('File already exists: %s',
                     toString(paste(outshape, c('shp', 'shx', 'dbf'),
                                    sep = '.')[f.exists])), call.=FALSE)
      } else (
        unlink(paste(outshape, c('shp', 'shx', 'dbf')))
      )
    }
    if (methods::is(x, 'Raster')) {
      
      raster::writeRaster(x, {f <- tempfile(fileext = '.tif')})
      rastpath <- normalizePath(f)
    } else if (is.character(x)) {
      rastpath <- normalizePath(x)
    } else stop('x must be a file path (character string), or a Raster object.')
    
    system2('python', args = (paste0(sprintf('"%1$s" "%2$s" -f "%3$s" "%4$s.shp"',
                                             pypath, rastpath, gdalformat, outshape), " -fieldname id")))
    if (isTRUE(readpoly)) {
      shp <- readshape(outshape)
      return(shp)
    }
    return(NULL)
  }
}