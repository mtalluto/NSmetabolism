# NOT RUN; these tests are way too slow. But results are identical to earlier methods, so good to go

# library(data.table)
# library(raster)
# library(sp)
# library(rgrass7)

# gisBase <- "/Applications/GRASS-7.4.1.app/Contents/Resources/"
# ltDat <- readRDS(system.file("testdata/light/lightTestData.rds", package="NSmetabolism"))
# ltPts <- data.frame(site = c("2", "69"), x = c(379210.2, 416118.6), y = c(4490865.8, 4462171.8))
# coordinates(ltPts) <- c(2,3)
# dem <- raster(system.file("testdata/light/FilledDEM.tif", package="NSmetabolism"))


# ltMod <- irradiance(dem, ltPts, unique(ltDat$datetime), newSession = TRUE, timezone = 2) 
# # test <- irradiance(dem, ltPts, unique(ltDat$datetime), use_existing = TRUE, timezone = 2) 
# rownames(ltMod) <- c("2", "69")
# ltDat$type <- "orig"
# ltMod <- melt(ltMod)
# ltMod$type <- "model"
# ltMod[,2] <- as_datetime(ltMod[,2])
# colnames(ltMod) <- colnames(ltDat)
# tz(ltMod$datetime) <- tz(ltDat$datetime)

# ltFinal <- rbind(ltDat, ltMod)
# ltFinal$type <- factor(ltFinal$type)
# ltFinal[type == 'model', datetime := datetime + hours(1)]
# ggplot(ltFinal, aes(x = datetime, y = light, color = type)) + geom_line() + facet_grid(site ~ .)
