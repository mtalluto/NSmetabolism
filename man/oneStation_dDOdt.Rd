% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/onestation.r
\name{oneStation_dDOdt}
\alias{oneStation_dDOdt}
\title{Time derivative of DO concentration}
\usage{
oneStation_dDOdt(t, y, parms, data, return_list = FALSE)
}
\arguments{
\item{t}{Time step for the ode integration}

\item{y}{Current state (dissolved oxygen concentration)}

\item{parms}{a named vector of parameters; see details}

\item{data}{a list of data for the ode; see details}

\item{return_list}{logical, should a list (suitable for deSolve) be returned?}
}
\value{
Time derivative of dissolved oxygen
}
\description{
Time derivative of DO concentration
}
\details{
\code{params} must be a vector including the following named items:
\itemize{
\item \code{P1} $(W min g^{-1} O_2)$ inverse of the slope of a photosynthesis–irradiance curve at low
light intensity
\item \code{P2} $(m^2 min g^{-1} O_2)$ inverse maximum photosynthesis rate; can be zero to assume
GPP is linear with light intensity instead of saturating
\item \code{k600} coefficient of gas exchange for a gas with a Schmidt number of 600 (see \code{\link[=kT]{kT()}}).
\item \code{ER24_20} daily ecosystem respiration rate, standardized at 20 degrees C
}

\code{data} is a named list of (constant) data items, including the following:
\itemize{
\item \code{PAR}  a function that takes time as a parameter and returns light intensity $(W/m^2)$
\item \code{temp} a function that takes time as a parameter and returns water temperature (degrees C)
\item \code{P}    pressure, in atmospheres
\item \code{z}    Depth, in meters
}
}
\references{
Fuß, T. et al. (2017). Land use controls stream ecosystem metabolism by shifting
dissolved organic matter and nutrient regimes. \emph{Freshw Biol} \strong{62}:582–599.\preformatted{ Van de Bogert MC et al. 2007. Assessing pelagic and benthic metabolism using free
	 water measurements. *Limnol. Oceanogr. Methods* **5**:145–155. 
}
}
