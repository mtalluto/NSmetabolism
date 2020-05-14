context("IDW")
library("NSmetabolism")

test_that("IDW-river works on single and multiple input values", {
	q <- c(1, 5, 6, 10, 16)
	dist <- matrix(c(
		 0, NA, -1, NA, -2,
		NA,  0, -1, NA, -2,
		 1,  1,  0, NA, -1,
		NA, NA, NA,  0, -1,
		 2,  2,  1,  1,  0), byrow=TRUE, nrow=length(q))
	vals <- c(1, 2, NA, NA, 3)

	# set up the lists
	idwVals <- nbQ <- idwDist <- list()
	for(i in 1:length(q)) {
		j <- !is.na(vals) & !is.na(dist[i,])
		idwVals[[i]] <- vals[j]
		nbQ[[i]] <- q[j]
		idwDist[[i]] <- dist[i,j]
	}

	expect_error(newVals <- idw_river(idwVals, idwDist, nbQ, q), regex = NA)
	j <- is.na(vals)
	expect_identical(newVals[!j], vals[!j])
	i <- c(1,2,5)
	j <- 3
	qr <- q[i] / q[j]
	qr[3] <- 1/qr[3]
	wt <- qr * 1/(dist[j,i]^2)
	wt <- wt/sum(wt)
	expect_equal(newVals[3], sum(vals[i] * wt))
	expect_equal(newVals[4], vals[5]) # 5 is the only neighbor, so they are equal
})

test_that("Regular IDW", {
	vals = matrix(c(5,10, 1, 3, 2, 4), nrow=2) 
	dmat = c(1, 10)
	expect_error(idw_matrix(vals, dmat[1]), regex="length")
	expect_error(idw_matrix(vals, dmat), regex=NA)
})
