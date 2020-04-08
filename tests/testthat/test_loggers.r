context("Logger input")
library("NSmetabolism")

test_that("New minidots are parsed correctly", {
	expect_silent(dat <- read_minidot(system.file("testdata/minidot_new", 
		package="NSmetabolism")))
	expect_equal(ncol(dat), 6)
	expect_equal(nrow(dat), 145*3) # 145 rows in 3 files
})

test_that("Old minidots are parsed correctly", {
	expect_silent(dat <- read_minidot(system.file("testdata/minidot_old", 
		package="NSmetabolism")))
	expect_equal(ncol(dat), 6)
	expect_equal(nrow(dat), 36*5) # 36 rows in 5 files
})

test_that("Hobo light/temperature pendants are parsed correctly", {
	# ambiguous name
	expect_error(read_hobo(system.file("testdata/hobo_lt", package="NSmetabolism"), "hobo")) 
	# not found
	expect_error(read_hobo(system.file("testdata/hobo_lt", package="NSmetabolism"), "herbo")) 
	expect_silent(dat <- 
		read_hobo(system.file("testdata/hobo_lt", package="NSmetabolism"), "hobo1"))
	expect_equal(ncol(dat), 4)
	expect_equal(nrow(dat), 1308)
})