#' Produce a model of irradiance for an entire study area or selected control points
#'
#' @param dem A [sp::SpatialGridDataFrame] or [raster::RasterLayer] object giving the digital
#' 		elevation model of the study area.
#' @param output_points A set out output points; either a [sp::SpatialPoints] object or a raster mask.
#' @param times A vector of [POSIXct] objects giving the dates and times at which to compute irradiance. These should be given in local (to the DEM) time
#' @param timezone An integer giving the offset in hours from UTC that applies to the DEM; see the manual for r.sun for details
#' @param gisBase `gisBase` parameter for call to [rgrass7::initGRASS()].
#' @param home Location to write temporary grass scratch files; defaults to tempdir().
#' @param horizon_pars List of arguments for GRASS74's r.horizon; see the GRASS GIS help files if you wish to change the defaults.

#' @details Requires an existing installation of GRASS7, as well as the rgrass7 package.
#' 
#' Note that running the light model can be extremely slow for large study areas and/or high-resolution DEMs and/or long time series. 

#' It is strongly recommended to use the `output_points` argument to set a list of locations at which to compute the output, otherwise memory usage can be extreme. For rasters, output will be returned for all non-NA cells.
#' @return A matrix, one row per output point, one column per date/time for computation. Column names are the integer values of the date/time, rownames are the coordinates of the point
irradiance <- function(dem, output_points, times, timezone, gisBase, home=tempdir(), 
	horizon_pars = list(step=30, bufferzone=200, maxdistance=5000, horBaseName = "horizonAngle")) 
{
	if(!requireNamespace('rgrass7')) 
		stop("The rgrass7 package must be installed to use this function")

	if("RasterLayer" %in% class(dem))
		demSPDF <- rasterToSPDF(dem)

	## init grass and provide data
	if(nchar(Sys.getenv("GISRC") == 0)) {
		rgrass7::initGRASS(gisBase, home=home, mapset="PERMANENT")
		rgrass7::execGRASS("g.proj", flags = "c", proj4 = proj4string(dem))
		ext <- as.character(as.vector(raster::extent(dem)))
		rasres <- as.character(res(dem))
		rgrass7::execGRASS("g.region", n = ext[4], s = ext[3], e = ext[2], w = ext[1], 
			rows=nrow(dem), cols=ncol(dem), nsres = rasres[2], ewres = rasres[1])
	}
	rgrass7::writeRAST(demSPDF, "dem")

	## build slope and aspect maps
	err <- rgrass7::execGRASS("r.slope.aspect", flags=c("overwrite"), elevation="dem", 
		aspect="aspect", slope="slope")

	## produce a longitute raster
	err <- rgrass7::execGRASS("r.latlong", flags=c('overwrite', 'l'), input='dem', 
		output='longitude')

	## pre-compute horizon angles, slow but saves significant time
	err <- rgrass7::execGRASS("r.horizon", flags=c("overwrite"), elevation='dem', 
		step=horizon_pars$step, bufferzone=horizon_pars$bufferzone, 
		output = horizon_pars$horBaseName, maxdistance=horizon_pars$maxdistance)

	## preallocate a large matrix to store results
	irrad_mat <- matrix(NA, ncol=length(times), nrow=length(output_points))
	colnames(irrad_mat) <- as.integer(times)
	rownames(irrad_mat) <- apply(coordinates(output_points), 1, paste, collapse=',')
	## for each time/date combo, run the light model
	for(i in 1:length(times))
	{
		dt <- times[i]
		dy <- lubridate::yday(dt)
		tm <- lubridate::hour(dt) + lubridate::minute(dt)/60
		err <- rgrass7::execGRASS("r.sun", flags=c("overwrite"), elevation="dem", 
			horizon_basename = horizon_pars$horBaseName, horizon_step = horizon_pars$step, 
			day=dy, time=tm, glob_rad="irradiance_out", aspect="aspect", slope="slope", 
			long="longitude",  civil_time=timezone)

		# bring the raster back into R
		irr <- raster::raster(rgrass7::readRAST("irradiance_out"))
		irrad_mat[,i] <- extract(irr, output_points)
	}
	irrad_mat
}


#' Produce a [sp::SpatialPixelsDataFrame] from a [raster::RasterLayer]
#'
#' @param x A [raster::RasterLayer] object
#' @return A [sp::SpatialPixelsDataFrame]
#' @keywords internal
rasterToSPDF <- function(x)
{
	coords <- sp::coordinates(x)
	gr <- data.frame(x=coords[,1], y=coords[,2], val=raster::values(x))
	sp::coordinates(gr) <- c(1,2)
	sp::proj4string(gr) <- sp::proj4string(x)
	sp::gridded(gr) <- TRUE
	return(gr)
}