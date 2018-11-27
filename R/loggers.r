#' Get data from a directory of minidot files in txt (csv) format
#'
#' @param fname Directory name containing minidot files to read
#' @param skip How many header (non-data) lines will be at the beginning of the files? Normally
#'		the default of 3 is fine.
#' @param ... Additional named arguments to be passed to [data.table::fread()] or [read.csv()]
#'
#' @details This function searches a directory for minidot data files, reads them in and returns
#'		a data.table or data.frame; all files found in the directory will be collected into a 
#'		single table.
#' @return A [data.table::data.table] or `data.frame`
read_minidot <- function(fname, skip = 3, ...)
{
	files <- list.files(fname, full.names=TRUE)
	if(requireNamespace('data.table', quietly=TRUE))
	{
		dat <- data.table::rbindlist(lapply(files, data.table::fread, skip = skip, ...))
	} else {
		dat <- do.call(rbind, lapply(files, read.csv, skip = skip, ...))
	}
	
	# drop battery column
	if(ncol(dat) == 5) 
		dat <- dat[,-2]

	colnames(dat) <- c('time_sec', 'temperature', 'DO', 'q')
	dat$timestamp <- as.POSIXct(dat[[1]], origin="1970-01-01")
	dat
}


#' Read hobo light/temperature logger data
#' 
#' @param fname file or directory name, see 'details'
#' @param pattern grep-compatible search pattern for the filename if fname is a directory
#'
#' @details If \code{fname} is a directory, then \code{pattern} is also required; this allows
#'   choosing the hobo by name
#' @return A [data.table::data.table] or `data.frame`
read_hobo <- function(fname, pattern)
{
	if(dir.exists(fname))
	{
		if(missing(pattern))
			stop("pattern must be specified because ", fname, " is a directory")
		files <- list.files(fname)
		file <- file.path(fname, files[grep(pattern, files)])
		if(length(file) == 0)
			stop("No files matching the pattern '", pattern, "'")
		if(length(file) > 1)
			stop("The pattern '", pattern, "' is ambiguous, matches:\n  ", paste(file, collapse='\n  '))
	} else if(file.exists(fname))
	{
		file <- fname
	} else {
		stop("file not found: ", fname)
	}
	if(requireNamespace('data.table', quietly=TRUE))
	{
		dat <- data.table::fread(file)
	} else {
		dat <- read.csv(file)
	}
	

	# we drop the row numnber and other columns as needed
	dropCols <- c(1, grep("Host Connected|Coupler|Stopped|End of File", colnames(dat), ignore.case=TRUE))
	dat[,dropCols] <- NULL

	# process dates into a date format, fix time zone
	datecol <- grep("Date", names(dat))
	tz_add <- -1 * as.numeric(sub("^.+GMT([+-]\\d?\\d):.+", "\\1", names(dat)[datecol]))
	dat$timestamp <- lubridate::parse_date_time(dat[[datecol]], orders="mdyIMSp") + lubridate::hours(tz_add)
	dat[[datecol]] <- NULL

	# capture the serial number, rename temperature column
	tempcol <- grep("Temp", names(dat))
	sn <- sub(".+LGR S/N: (\\d+).+", "\\1", names(dat)[tempcol])
	names(dat)[tempcol] <- "temperature"
	dat$serial_number <- sn

	# light
	lightcol <- grep("Lux", names(dat))
	names(dat)[lightcol] <- "light"

	dat
}