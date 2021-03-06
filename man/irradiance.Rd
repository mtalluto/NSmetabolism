% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/light.r
\name{irradiance}
\alias{irradiance}
\title{Produce a model of irradiance for an entire study area or selected control points}
\usage{
irradiance(
  dem,
  output_points,
  times,
  timezone,
  use_existing = FALSE,
  dem_name = "dem",
  slope_name = "slope",
  aspect_name = "aspect",
  lon_name = "longitude",
  horizon_pars = list(step = 30, bufferzone = 200, maxdistance = 5000, horBaseName =
    "horizonAngle"),
  newSession = FALSE,
  gisBase
)
}
\arguments{
\item{dem}{A \link[sp:SpatialGridDataFrame]{sp::SpatialGridDataFrame} or \link[raster:RasterLayer]{raster::RasterLayer} object giving the digital
elevation model of the study area.}

\item{output_points}{A set out output points; either a \link[sp:SpatialPoints]{sp::SpatialPoints} object or a raster mask.}

\item{times}{A vector of \link{POSIXct} objects giving the dates and times at which to compute irradiance. These should be given in local (to the DEM) time}

\item{timezone}{An integer giving the offset in hours from UTC that applies to the DEM; see the manual for r.sun for details}

\item{use_existing}{logical, if TRUE then existing slope, aspect, longitute, and horizon
rasters will be used, if FALSE then they will be recomputed}

\item{dem_name}{character; this parameter (and similar name parameters) tell GRASS what names to
use for the rasters needed for the analysis; they must exist in the GRASS location/mapset if
\code{use_existing == TRUE}.}

\item{horizon_pars}{List of arguments for GRASS74's r.horizon; see the GRASS GIS help files ]
if you wish to change the defaults.}

\item{newSession}{Boolean, if TRUE a new grass session will be started even if one already
exists}

\item{gisBase}{character; the location of the GRASS installation
(see WatershedTools::GrassSession)}
}
\value{
A matrix, one row per output point, one column per date/time for computation. Column names are the integer values of the date/time, rownames are the coordinates of the point
}
\description{
Produce a model of irradiance for an entire study area or selected control points
}
\details{
Requires an existing installation of GRASS7, as well as the rgrass7 package. It is
first necessary to run \link[rgrass7:initGRASS]{rgrass7::initGRASS}. Recommended settings for initGRASS include
\code{home = tempdir()} and \code{SG = SpatialGrid(dem_gridded)}. To produce the gridded DEM
from a raster, use:
\code{dem_gridded <- SpatialPoints(dem)}
\code{gridded(dem_gridded) <- TRUE}

Note that running the light model can be extremely slow for large study areas and/or high-resolution DEMs and/or long time series.

It is strongly recommended to use the \code{output_points} argument to set a list of locations at which to compute the output, otherwise memory usage can be extreme. For rasters, output will be returned for all non-NA cells.
}
