
dat <- readRDS("~/work/projects/vjosa_biodiversity_ef_ms/dat/one_station_data.rds")
dat <- dat[c("26-VjosaSpring2018", "65-VjosaSpring2018")]


light <- do.call(rbind, lapply(dat, function(x) x$light))
water_temp <- do.call(rbind, lapply(dat, function(x) x$temp))
pressure <- do.call(rbind, lapply(dat, function(x) x$pressure))
depth <- sapply(dat, function(x) x$z)

# these parameters were snatched from model runs on these two sites
DO_i <- c(8.57150491, 8.76083385)
lP1 <- c(6.75929342, 7.82380862)
lP2 <- c(-3.50763173, -2.15861567)
er24_20 <- c(-7.95076799, -25.17492411)
k600 <- c(74.73116152, 12.87476450)

saveRDS(list(light = light, water_temp = water_temp, pressure = pressure, DO_i = DO_i, 
	depth = depth, lP1 = lP1, lP2 = lP2, er24_20 = er24_20, k600 = k600), 
	"inst/testdata/site_metab_data/site_metab_data.rds")
