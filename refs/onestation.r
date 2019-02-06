#############  
## matt notes
# these are kept for historical significance
# no longer using the nll-based calibration, using a bayesian method instead
# GPP24 and insituER24 may be useful later on, but they should be made more flexible and they
# should be able to integrate properly across any measurement period


# oxyfile has the empirical data: original times and measured oxygen values (determine the fitting quality of the model) and observed values for forcing functions light and temp
# parameters P1, P2 and ER24 [MET] are estimated by minimizing the negative log likelihood. 
# optionally K600 is fixed (and is an exchange velocity, not a rate!).
DO_nll<-function(MET, K600, start, end, timestep, oxyfile, parfun, tempfun, z, BP, prFunction = DO_predict_discrete)
{
	# discrete approach
	if(!missing(K600))
		MET[4] <- K600
	names(MET) <- c('P1', 'P2', 'ER24', 'K600')
	times <- seq(start, end, timestep)
	oxy.m <- oxyfile$DO # "measured" oxygen
	oxy.p <- prFunction(params=MET, DO_initial=oxy.m[1], intTimes=times, parfun=parfun, tempfun=tempfun, z, BP)
	# print(oxy.p)
	oxy.p2 <- oxy.p[times %in% oxyfile$TimeMin]
	
	sqdiff <- (oxy.m - oxy.p2)^2
	return(length(oxy.m)*(log(((sum(sqdiff)/length(oxy.m))^0.5)) + 0.5 * log(6.28)) + ((2*sum(sqdiff)/length(oxy.m))^-1)*sum(sqdiff) )  # Negative Loglikelhood (Hilborn and Mangel 1997, page 137)
}



#GPP calculation - TAKE CARE: integrates over whole measurement period, not just 24 hours, therefore only correct when measurement period starts anytime at night and ends anytime during following night (and both nights are dark)
# timestep... must be timestep of datafile!
GPP24calc<-function(P1, P2, oxyfile, timestep){
	PAR<-(oxyfile$Light) 
	GPP24<-sum(  PAR/(P1+P2*PAR)*timestep ); names(GPP24)<-"GPP24"
	return(GPP24)
}


# ER24 is modelled standardized to 20°C, for later NEP and P/R computation ER at in-situ conditions may be useful
# TAKE CARE: integrates over whole measurement period, not just 24 hours
# timestep... must be timestep of datafile!
insituER24calc<-function(ER24,oxyfile,timestep) {
	sum( ER24/(24*60)*(1.045^(oxyfile$Temp-20))*timestep )  
}