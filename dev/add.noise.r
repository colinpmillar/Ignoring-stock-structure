###############################################################################
# EJ,CM,IM,CO(20120430)
# Functions to introduce variability in FLR objects and simulate abundance
# indices.
###############################################################################


genStochN <- function(object, cv=0.1, niters){
	require(MASS)
	mu <- log(stock.n(object))
	Rho <- cor(t(mu[drop=TRUE]))
	flq <- mvrnorm(niters*dim(mu)[2], rep(0,nrow(Rho)), cv^2*Rho)
	mu <- propagate(mu, niters)
	flq <- FLQuant(c(t(flq)), dimnames=dimnames(mu))
	exp(mu+flq)
}

setMethod("z", "FLStock", function(object, ...){
	f <- harvest(object)
	if(units(f)=='harvest'){
		stop("Your exploitation is defined in harvest rates, can not be added to natural mortality")
	} else { 
		m <- m(object)
		m+f
	}
})

genStochF <- function(object, cv=0.1, niters){
	require(MASS)
	mu <- log(harvest(object))
	Rho <- cor(t(mu[drop=TRUE]))
	flq <- mvrnorm(niters*dim(mu)[2], rep(0,nrow(Rho)), cv^2*Rho)
	mu <- propagate(mu, niters)
	flq <- FLQuant(c(t(flq)), dimnames=dimnames(mu))
	flq <- exp(mu+flq)
	units(flq) <- 'f'
	flq
}

genStochStk <- function(object, n.cv, f.cv=n.cv, niters){

	nyrs <- dim(catch.n(object))[2]
	nages <- dim(catch.n(object))[1]

	# initiate the object with iters
	sstk <- propagate(object, niters)

	# generate stochastic [N] with lognormal and cv
	Ns <- genStochN(object, cv=n.cv, niters)

	# update R
	stock.n(sstk)[1] <- Ns[1]

	# generate stochastic [F] with lognormal and cv
	Fs <- genStochF(object, cv=f.cv, niters)

	# compute cumulative Z
	Z <- FLCohort(Fs+m(object))
	Z[] <- apply(Z, c(2:6), cumsum)

	# expand variability into [N] by R*[F] 
	Ns <- FLCohort(Ns)[rep(1,nages)]
	Ns <- Ns*exp(-Z)
	Ns <- flc2flq(Ns)

	# Update object
	# [N]
	stock.n(sstk)[-1,-1] <- Ns[-nages,-nyrs] 
	# plus group
	stock.n(sstk)[nages,-1] <- Ns[nages,-1] + stock.n(sstk)[nages,-1]
	# [F]
	harvest(sstk) <- Fs
	# [C]
	catch.n(sstk) <- harvest(sstk)/z(sstk)*(1-exp(-z(sstk)))*stock.n(sstk)
	catch(sstk) <- computeCatch(sstk)
	# [L]
	landings.n(sstk) <- catch.n(sstk)
	landings(sstk) <- catch(sstk)
	# out
	sstk
}

