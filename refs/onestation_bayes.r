library(rstan)
library(data.table)
dat <- fread("dat/1In_Welse_May.txt")
standat <- list(nt = nrow(dat), dt = 10, PAR = dat$Light, temp = dat$Temp, do_initial = dat$DO[1], z = 0.47, BP = 0.998815, DO = dat$DO)

fit <- stan(file = 'R/onestation.stan', data = standat, iter = 5000, chains = 3)

samples <- extract(fit, permuted = FALSE)
mcmc_trace(samples, pars = c('P1', 'P2', 'ER24', 'K600', 'sigma', 'gpp24', 'insituER24'))

dims <- grep("do_pr", dimnames(fit)$parameters)
dopr <- apply(samples[,,dims], 3, mean)
dopr_qu <- apply(samples[,,dims], 3, quantile, c(0.025, 0.975))

quartz(w=7, h=7, file="do.png", dpi=72, type='png')
plot(1:145, dat$DO, pch=16, ylim=c(7.2,9.2) , xlab="Time Step", ylab = "DO", cex=0.7)
polygon(x=c(1:145, rev(1:145)), y = c(dopr_qu[1,], rev(dopr_qu[2,])), border=NA, col='#3333ff55')
lines(1:145, dopr, col='#3333ff')
dev.off()