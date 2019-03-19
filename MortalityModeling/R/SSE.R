## function for building up a block diagonal
## from a list of matrices
bdiag <- function(x){ 
	if(!is.list(x)) stop("x not a list") 
	n 			<- length(x) 
	if(n==0) return(NULL) 
	x      <- lapply(x, function(y) 
				if(length(y)) as.matrix(y) else 
					stop("Zero-length component in x")) 
	d      		<- array(unlist(lapply(x, dim)), c(2, n)) 
	rr 			<- d[1,] 
	cc 			<- d[2,] 
	rsum 		<- sum(rr) 
	csum 		<- sum(cc) 
	out 		<- array(0, c(rsum, csum)) 
	ind 		<- array(0, c(4, n)) 
	rcum 		<- cumsum(rr) 
	ccum 		<- cumsum(cc) 
	ind[1,-1] 	<- rcum[-n] 
	ind[2,] 	<- rcum 
	ind[3,-1] 	<- ccum[-n] 
	ind[4,] 	<- ccum 
	imat 		<- array(1:(rsum * csum), c(rsum, csum)) 
	iuse <- apply(ind, 2, function(y, imat)
				imat[(y[1]+1):y[2],
						(y[3]+1):y[4]], imat=imat)
	iuse 		<- as.vector(unlist(iuse)) 
	out[iuse] 	<- unlist(x) 
	return(out) 
}

# main SSE function. Could still be modularized quite a bit
SSEdecomposition <- function(Z,E,ages,years,sex="Males",print=F){
	## packages
	require(MortalitySmooth)
	flag.stop <- FALSE
	
	# -----------------------------------
	# TR: FYI you can make this not need a matrix. Say Z and E given as 1d vectors,
	# then just force to matrix() with ncol = 1,
	# and be sure to not inadvertently lose this dimensionality
	# in later steps by using the arg drop = FALSE
	# -----------------------------------
	# OR make it run on a vector, for the 1d case, since years are independant,
	# and then write a wrapper for the 2d case that uses apply() or a for loop
	# over years...
	# -----------------------------------
	
	## age axis
	m 		<- length(ages)
	n 		<- length(years)
	x0 		<- ages
	m0 		<- length(x0)
	x 		<- x0[-1]  ## no infant
	m 		<- length(x)
	
	## log mortality in matrix
	LHAZact <- log(Z/E)
	
	## select over which ages each component should be restricted
	## infant component
	staINF 	<- 1
	endINF 	<- 50
	## aging/senescence component
	staAGI 	<- 1
	endAGI 	<- m
	## accident-hump (middle-mortality) component
	staACC 	<- 1
	endACC 	<- 80
	
	## select data and create weights for each age-window
	## infant mortality
	p1 		<- c(staINF, endINF)
	x1 		<- x[x>=p1[1] & x<=p1[2]]
	m1 		<- length(x1)
	w1 		<- as.numeric(x%in%x1)
	x1A 	<- 0:x1[m1]
	## aging mortality
	p2 		<- c(staAGI, endAGI)
	x2 		<- x[x>=p2[1] & x<=p2[2]]
	m2 		<- length(x2)
	w2 		<- as.numeric(x%in%x2)
	## accident mortality
	p3 		<- c(staACC, endACC)
	x3 		<- x[x>=p3[1] & x<=p3[2]]
	m3 		<- length(x3)
	w3 		<- as.numeric(x%in%x3)
	
	## constructing the bases
	
	## B-splines for the infant part
	deg1 	<- floor(length(x1)/3)
	xl1 	<- min(x1)
	xr1 	<- max(x1)
	xmax1 	<- xr1 + 0.01 * (xr1 - xl1)
	xmin1 	<- xl1 - 0.01 * (xr1 - xl1)
	B1 		<- MortSmooth_bbase(x1, xmin1, xmax1, deg1, 3)
	nb1 	<- ncol(B1)
	## difference matrices
	D1 		<- diff(diag(nb1), diff=2)
	tDD1 	<- t(D1) %*% D1
	
	## B-splines for the aging part
	deg2 	<- floor(length(x2)/3)
	xl2 	<- min(x2)
	xr2 	<- max(x2)
	xmax2 	<- xr2 + 0.01 * (xr2 - xl2)
	xmin2 	<- xl2 - 0.01 * (xr2 - xl2)
	B2 		<- MortSmooth_bbase(x2, xmin2, xmax2, deg2, 3)
	nb2 	<- ncol(B2)
	## difference matrices
	D2 		<- diff(diag(nb2), diff=2)
	tDD2 	<- t(D2) %*% D2
	
	## B-splines for the accident-hump
	deg3 	<- floor(length(x3)/3)
	xl3 	<- min(x3)
	xr3 	<- max(x3)
	xmax3 	<- xr3 + 0.01 * (xr3 - xl3)
	xmin3 	<- xl3 - 0.01 * (xr3 - xl3)
	B3 		<- MortSmooth_bbase(x3, xmin3, xmax3, deg3, 3)
	nb3 	<- ncol(B3)
	## difference matrices
	D3 		<- diff(diag(nb3), diff=3)
	tDD3 	<- t(D3) %*% D3
	
	## shape penalties
	## including monotonicy for infant part
	D1mon 	<- diff(diag(nb1), diff=1)
	w1mon 	<- rep(0, nrow(D1mon))
	W1mon 	<- diag(w1mon)
	## including monotonicy for aging part
	D2mon 	<- diff(diag(nb2), diff=1)
	w2mon 	<- rep(0, nrow(D2mon))
	W2mon 	<- diag(w2mon)
	## including log-concaveness for accident-hump
	D3con 	<- diff(diag(nb3), diff=2)
	w3con 	<- rep(0, nrow(D3con))
	W3con 	<- diag(w3con)
	kappa 	<- 10^6
	
	## smoothing parameters
	nl 			<- 7
	lambdas1 	<- 10^seq(2.5,5,length.out=nl)
	nl1 		<- length(lambdas1)
	lambdas2 	<- 10^seq(4,8,length.out=nl)
	nl2 		<- length(lambdas2)
	lambdas3 	<- 10^seq(3,6.5,length.out=nl)
	nl3 		<- length(lambdas3)
	
	## axis and basis for smooth components
	xx1A 	<- seq(0, max(x1), by=0.1)
	xx1 	<- seq(min(x1), max(x1), by=0.1)
	ms1 	<- length(xx1)
	BB1 	<- MortSmooth_bbase(xx1, xmin1, xmax1, deg1, 3)
	BB1A 	<- MortSmooth_bbase(xx1A, xmin1, xmax1, deg1, 3)
	xx2 	<- seq(min(x2), max(x2), by=0.1)
	ms2 	<- length(xx2)
	BB2 	<- MortSmooth_bbase(xx2, xmin2, xmax2, deg2, 3)
	xx3 	<- seq(min(x3), max(x3), by=0.1)
	ms3 	<- length(xx3)
	BB3 	<- MortSmooth_bbase(xx3, xmin3, xmax3, deg3, 3)
	
	## complete model matrix as a list
	## in their order and with dimensions
	XX 		<- list(X1=B1, X2=B2, X3=B3)
	nx 		<- length(XX)
	nc 		<- unlist(lapply(XX, ncol))
	## indicators for the coefficients of each basis
	ind2 	<- cumsum(nc)
	ind1 	<- c(1, ind2[1:(nx-1)]+1)
	ind 	<- NULL
	for(i in 1:nx){
		ind[[i]] <- ind1[i]:ind2[i]
	}
	## indicators for the fitted values
	indF 	<- cbind(w1, w2, w3)
	
	## RESULTS NEEDED
	E1 		<- LHAZ1hat <- Z1hat <- matrix(NA,(m1+1),n)
	E2 		<- LHAZ2hat <- Z2hat <- matrix(NA,m2,n)
	E3 		<- LHAZ3hat <- Z3hat <- matrix(NA,m3,n)
	LHAZhat <- matrix(NA,m0,n)
	GAMMA1hatSMO 	<- matrix(NA,ms1,n)
	GAMMA2hatSMO 	<- matrix(NA,ms2,n)
	GAMMA3hatSMO 	<- matrix(NA,ms3,n)
	GAMMA1hatSMO_0 	<- matrix(NA,ms1+10,n)
	SMOOTHpar 		<- POSpar <- matrix(NA,3,n)
	ED1 <- ED2 <- ED3 <- numeric(n)
	
	## LOOP OVER YEARS
	PLOT <- FALSE
	jj   <- 1
	for (jj in 1:n){
		#cat("Fitting year",years[jj],"\n")
		## subset years, exposures and weights
		y 			<- Z[-1,jj]
		e 			<- E[-1,jj]
		y[e==0] 	<- 0  ## to avoid issues with NA
		lhaz.act 	<- log(y/e)
		# plot(x,lhaz.act)
		haz.act 	<- y/e
		inf.age <- x[which.min(haz.act[haz.act>0])]
		wei <- rep(1, m)
		wei[e==0] <- 0
		
		## update deaths and exposures
		## infant mortality
		y1 <- y[which(w1==1)]
		e1 <- e[which(w1==1)]
		## aging mortality
		y2 <- y[which(w2==1)]
		e2 <- e[which(w2==1)]
		## accident mortality
		y3 <- y[which(w3==1)]
		e3 <- e[which(w3==1)]
		
		if (jj==1){
			## starting values
			
			## fitting P-splines for infant mortality
			## taking only the first 8 ages and extrapolating
			## for the successive ages
			ww1 <- rep(0,length(y1))
			ww1[x1<=inf.age] <- 1
			y1a <- c(y)
			y1a[which(y1a<=0)] <- 0
			y1 <- y1a[which(w1==1)]
			options(warn = -1)
			fit1 <- Mort1Dsmooth(x=x1, y=y1, offset=log(e1),
					w=ww1,
					method=3, lambda=10^1,
					ndx = deg1, deg = 3, pord = 2)
			options(warn = 0)
			y1.st0 <- exp(XX[[1]] %*% fit1$coef)*e1
			y1.st <- rep(0,m)
			y1.st[which(w1==1)] <- y1.st0
			lhaz1.st <- log(y1.st/e)
			# plot(x, lhaz.act, ylim=c(-12,2))
			# lines(x, lhaz1.st, col=3, lwd=2, lty=3)
			
			## fitting P-splines for the aging
			## taking only ages 50+ and extrapolating
			## for the previous ages
			y2 <- y[which(w2==1)]
			ww2 <- rep(0,length(y2))
			ww2[x2>=40] <- 1
			options(warn = -1)
			fit2 <- Mort1Dsmooth(x=x2, y=y2, offset=log(e2),
					w=ww2*wei[w2==1],
					method=3, lambda=10^3,
					ndx = deg2, deg = 3, pord = 2)
			options(warn = 0)
			y2.st0 <- exp(XX[[2]] %*% fit2$coef)*e2
			y2.st <- rep(0,m)
			y2.st[which(w2==1)] <- y2.st0
			
			## fitting P-splines for the accident component
			## by fitting the differences between the
			## overall mortality and the infant+aging components
			## and over ages 8:45 and extrapolating
			## previous and successive ages
			y3a <- c(y-y1.st-y2.st)
			y3a[which(y3a<=0)] <- 0
			y3 <- y3a[which(w3==1)]
			ww3 <- rep(0, length(x3))
			ww3[x3>inf.age & x3<40] <- 1
			wei3 <- wei[which(w3==1)]
			options(warn = -1)
			fit3 <- Mort1Dsmooth(x=x3, y=y3, offset=log(e3),
					w=ww3*wei3,
					method=3, lambda=10^4,
					ndx = deg3, deg = 3, pord = 3)
			options(warn = 0)
			y3.st0 <- exp(XX[[3]] %*% fit3$coef)*e3
			y3.st <- rep(0,m)
			y3.st[which(w3==1)] <- y3.st0
			
			## starting values at the log-mortality scale
			lhaz1.st <- log(y1.st/e)
			lhaz2.st <- log(y2.st/e)
			lhaz3.st <- log(y3.st/e)
			y.st <- y1.st + y2.st + y3.st
			lhaz.st <- log(y.st/e)
			
			## plotting starting values
			if(PLOT){
				plot(x, lhaz.act, ylim=c(-12,2))
				lines(x, lhaz.st, col=2, lwd=3)
				lines(x, lhaz1.st, col=3, lwd=2, lty=3)
				lines(x, lhaz2.st, col=4, lwd=2, lty=3)
				lines(x, lhaz3.st, col=5, lwd=2, lty=3)
			}
			## concatenating starting coefficients
			coef.st <- as.vector(c(fit1$coef,
							fit2$coef,
							fit3$coef))
			
		}else{
			## take previous year optimal value
			coef.st <- coef.hat
		}
		
		## objects for saving outcomes
		## for each lambdas combination
		COEFS <- array(NA, dim=c(sum(nc), nl1, nl2, nl3),
				dimnames=list(NULL, lambdas1, lambdas2, lambdas3))
		COEFS[,nl1,nl2,nl3] <- coef.st
		BICs <- array(NA, dim=c(nl1,nl2,nl3),
				dimnames=list(lambdas1, lambdas2, lambdas3))
		EDs1 <- array(NA, dim=c(nl1,nl2,nl3),
				dimnames=list(lambdas1, lambdas2, lambdas3))
		EDs2 <- array(NA, dim=c(nl1,nl2,nl3),
				dimnames=list(lambdas1, lambdas2, lambdas3))
		EDs3 <- array(NA, dim=c(nl1,nl2,nl3),
				dimnames=list(lambdas1, lambdas2, lambdas3))
		ITs <- array(NA, dim=c(nl1,nl2,nl3),
				dimnames=list(lambdas1, lambdas2, lambdas3))
		
		
		## matrix which create the position
		## in a 3-dimensional search
		## for getting "warm" starting values,
		## i.e. previously estimated coefficients
		POS <- array(0, dim=c(nl1,nl2,nl3))
		for(i in 1:nl3){
			POS[,,i] <- matrix(seq(prod(dim(POS)),nl3,-nl3)-i+1,
					nl1,nl2,byrow=TRUE)
		}
		pos0 <- rep(  c(rep(nl2*nl3,nl3), rep(nl3,(nl2-1)*nl3) ), nl1)
		pos <- 1:prod(dim(POS))- pos0
		pos[pos<0] <- c(1,1:(nl3-1))
		
		
		## solving step
		d <- 0.3
		## maximum number of iteration
		max.it <- 200
		
		## should we plot current estimations?
		PLOT <- F
		## starting objects for
		## getting properly "warm" starting values
		coef <- 1
		it <- 1
		
		## iterations for each smoothing parameters
		## combinations
		## this takes about 22 seconds on a 
		## a 2.6 GHz portable PC with
		## an Intel?? Core??? i5-3320M processor
		## and 4 GB of RAM
		l1 <- 1
		l2 <- 1
		l3 <- 1
		
		for(l1 in 1:nl1){
			## smoothing parameter for the infant
			lambda1 <- lambdas1[l1]      
			for(l2 in 1:nl2){
				## smoothing parameter for the aging
				lambda2 <- lambdas2[l2]    
				for(l3 in 1:nl3){
					## smoothing parameter for the accident-hump
					lambda3 <- lambdas3[l3]
					
					## component specific penalty terms
					P1 <- lambda1 * tDD1
					P2 <- lambda2 * tDD2
					P3 <- lambda3 * tDD3
					## overall penalty term
					P <- bdiag(list(P1, P2, P3))
					
					## "warm" starting values:
					## previously estimated coefficients
					## additionally whenever either
					## max.it is reached or one encountered
					## issue in the estimation (NA coefficients)
					## it takes the original starting values
					cur.pos <- POS[l1,l2,l3]
					whi.coe <- which(POS==pos[cur.pos],
							arr.ind=TRUE)
					if(it==max.it | any(is.na(coef))){
						whi.coe <- c(nl1,nl2,nl3)
					}
					coef.st <- COEFS[,whi.coe[1],
							whi.coe[2],
							whi.coe[3]]
					coef <- coef.st
					
					## if the first iteration does not work, then problem
					if(any(is.na(coef))&class(fit3)=="Mort1Dsmooth"){
						coef <- as.vector(c(fit1$coef,fit2$coef,fit3$coef))
					}else if(any(is.na(coef))&class(fit3)!="Mort1Dsmooth"){
						coef <- coef.hat
					}
					## plotting actual log-mortality
					if(PLOT) plot(x, lhaz.act, ylim=c(-12,2))
					
					## convergence flag is true
					conv <- TRUE
					
					## iterations for a given lambdas-combination
					for(it in 1:max.it){
						
						## penalty for the shape constraints
						P1mon <- kappa * t(D1mon) %*%W1mon%*% D1mon
						P2mon <- kappa * t(D2mon) %*%W2mon%*% D2mon
						P3con <- kappa * t(D3con) %*%W3con%*% D3con
						Psha <- bdiag(list(P1mon, P2mon, P3con))
						
						## linear predictor
						eta <- numeric(nx*m)
						for(i in 1:nx){
							eta0 <- rep(0, m)
							eta0[which(indF[,i]==1)] <- XX[[i]] %*% coef[ind[[i]]]
							eta[1:m+(i-1)*m] <- eta0
							if(PLOT){ ## plotting each component 
								lines(x[which(indF[,i]==1)],
										eta0[which(indF[,i]==1)], col=i+2)
							}
						}
						## components
						gamma <- exp(eta)*c(indF)
						## expected values
						mu <- numeric(m)
						for(i in 1:nx){
							mu <- (e * gamma[1:m+(i-1)*m]) + mu
						}
						## plotting overall log-mortality
						if(PLOT) lines(x, log(mu/e), col=2, lwd=4)
						## weights for the IWLS
						w <- mu
						## modified model matrix for a CLM
						U <- matrix(NA, m, sum(nc))
						for(i in 1:nx){
							u <- gamma[1:m+(i-1)*m]/mu * e
							XXi <- matrix(0, nrow=m, ncol=nc[i])
							XXi[which(indF[,i]==1), ] <- XX[[i]]
							U0 <- u * XXi
							U[,ind[[i]]] <- U0
						}
						U[is.nan(U)] <- 0
						## regression parts for the P-CLM
						tUWU <- t(U) %*% (w*wei * U)
						tUWUpP <- tUWU + P + Psha
						r <- y - mu
						tUr <- t(U) %*% r
						## check eventual problem with convergence?
						tryInv <- try(solve(tUWUpP,
										tUr + tUWU %*% coef),
								silent=TRUE)
						## if so, set "conv" equal to FALSE
						## and break the loop w/o interrupting
						## the overall lambdas optimization        
						if(!is.matrix(tryInv) |
								any(is.na(tryInv)) |
								any(is.nan(tryInv))){
							conv <- FALSE
							break
						}
						## updating coefficients with a d-step 
						coef.old <- coef
						coef <- tryInv
						coef <- d*coef.old + (1-d)*coef
						## update weights for shape constraints
						## infant, monotonicity
						W1mon.old <- W1mon
						coef1 <- coef[ind[[1]]]
						W1mon <- diag(diff(coef1) >= 0)
						## aging, monotonicity
						W2mon.old <- W2mon
						coef2 <- coef[ind[[2]]]
						W2mon <- diag(diff(coef2) <= 0)
						## accident-hump, log-concaveness
						W3con.old <- W3con
						coef3 <- coef[ind[[3]]]
						W3con <- diag(diff(coef3, diff=2) >= 0)
						## convergence criterion for coefficients
						dif.coef <- max(abs(coef.old-coef))/max(abs(coef))
						## stopping loop at convergence
						if(dif.coef < 1e-04 & it > 4) break
					}
					
					## if we reached convergence
					if(conv){
						## compute devaince
						yy <- y
						yy[y==0] <- 10^-8
						mumu <- mu
						mumu[mu==0] <- 10^-8
						dev <- 2*sum(y * log(yy/mumu))
						## effective dimensions
						H <- solve(tUWUpP, tUWU)
						diagH <- diag(H)
						ed1 <- sum(diagH[ind[[1]]])
						ed2 <- sum(diagH[ind[[2]]])
						ed3 <- sum(diagH[ind[[3]]])
						## BIC
						bic <- dev + log(m)*sum(diagH)
					}else{ ## if we did NOT reached convergence
						coef <- rep(NA, ncol(P))        
						ed1 <- NA
						ed2 <- NA
						ed3 <- NA
						bic <- NA
					}
					## saving estimated values
					COEFS[,l1,l2,l3] <- as.vector(coef)
					BICs[l1,l2,l3] <- bic
					EDs1[l1,l2,l3] <- ed1
					EDs2[l1,l2,l3] <- ed2
					EDs3[l1,l2,l3] <- ed3
					ITs[l1,l2,l3] <- it
					# ## monitoring part
					# cat(l1, "of", nl1, "\n",
					#     l2, "of", nl2, "\n",
					#     l3, "of", nl3, "\n",
					#     it, dif.coef, "\n")
					if(conv & it != max.it){
						if(print) cat(it, dif.coef, bic, "\n")
						break
					}
				}
				if(conv & it != max.it){
					break
				}
			}
			if(conv & it != max.it){
				break
			}  
		}  
		
		if (!conv){
			cat("No convergence, change starting Z","\n")
			flag.stop <- T
			break
		} 
		
		## take BICs values only when
		## convergence is reached 
		BICs0 <- BICs
		BICs[ITs==max.it] <- NA
		## selecting optimal smoothing parameters
		pmin <- which(BICs==min(BICs, na.rm=TRUE),
				arr.ind=TRUE)
		lambda1.hat <- lambdas1[pmin[1]]
		lambda2.hat <- lambdas2[pmin[2]]
		lambda3.hat <- lambdas3[pmin[3]]
		
		## fitted coefficients for the optimal lambdas
		coef.hat <- COEFS[,pmin[1],pmin[2],pmin[3]]
		
		## effective dimensions at the optimal lambdas
		ed1.hat <- EDs1[pmin[1],pmin[2],pmin[3]] 
		ed2.hat <- EDs2[pmin[1],pmin[2],pmin[3]]
		ed3.hat <- EDs3[pmin[1],pmin[2],pmin[3]]
		
		## linear predictor for each component
		etas <- NULL
		for(i in 1:nx){
			etas[[i]] <- XX[[i]] %*% coef.hat[ind[[i]]]
		}
		eta1.hat <- etas[[1]]
		eta2.hat <- etas[[2]]
		eta3.hat <- etas[[3]]
		
		## linear predictor in a matrix over the whole x
		ETA.hat <- matrix(NA, m, nx)
		for(i in 1:nx){
			ETA.hat[which(indF[,i]==1),i] <- etas[[i]]
		}
		## linear predictor for overall mortality
		eta.hat <- log(apply(exp(ETA.hat), 1, sum, na.rm=TRUE))
		
		## fitted values for each component
		gamma1.hat <- exp(eta1.hat)
		gamma2.hat <- exp(eta2.hat)
		gamma3.hat <- exp(eta3.hat)
		
		## expected values
		mu.hat <- exp(eta.hat)*e
		
		
		## plotting actual and fitted log-mortality
		testofit <- c(expression(paste("ln(",y, "/e)")),
				expression(paste("ln(", mu, "/e)")),
				expression(paste("ln(",gamma[1],")")),
				expression(paste("ln(",gamma[2],")")),
				expression(paste("ln(",gamma[3],")")))
		if (PLOT){
			plot(x, lhaz.act, col=8, pch=16, ylim=c(-12, 2),
					xlab="age",main=years[jj], 
					ylab="log-mortality")
			lines(x, eta.hat, col=2, lwd=2)
			lines(x1, eta1.hat, col=3)
			lines(x, eta2.hat, col=4)
			lines(x3, eta3.hat, col=5)
			legend("topleft", legend=testofit, col=c(8,2:5),
					lty=c(-1,1,1,1,1), pch=c(16,-1,-1,-1,-1))
		}
		
		## evaluate components over a finer grid
		gamma1.hatAUG <- exp(BB1 %*% coef.hat[ind[[1]]])
		gamma2.hatAUG <- exp(BB2 %*% coef.hat[ind[[2]]])
		gamma3.hatAUG <- exp(BB3 %*% coef.hat[ind[[3]]])
		gamma1.hatAUG0 <-  exp(c(seq(log(Z[1,jj]/E[1,jj]),log(gamma1.hatAUG[1]),length.out = 11)[1:10],
						log(gamma1.hatAUG)))
		
		## save results
		E1[,jj] <- c(E[1,jj],e1)
		E2[,jj] <- e2
		E3[,jj] <- e3
		LHAZhat[,jj] <- c(log(Z[1,jj]/E[1,jj]),eta.hat)
		LHAZ1hat[,jj] <- c(log(Z[1,jj]/E[1,jj]),eta1.hat)
		LHAZ2hat[,jj] <- eta2.hat
		LHAZ3hat[,jj] <- eta3.hat
		GAMMA1hatSMO[,jj] <- gamma1.hatAUG
		GAMMA1hatSMO_0[,jj] <- gamma1.hatAUG0
		GAMMA2hatSMO[,jj] <- gamma2.hatAUG
		GAMMA3hatSMO[,jj] <- gamma3.hatAUG
		Z1hat[,jj] <- c(Z[1,jj],gamma1.hat*e1)
		Z2hat[,jj] <- gamma2.hat*e2
		Z3hat[,jj] <- gamma3.hat*e3
		SMOOTHpar[,jj] <- c(lambda1.hat,lambda2.hat,lambda3.hat)
		POSpar[,jj] <- pmin
		ED1[jj] <- ed1.hat
		ED2[jj] <- ed2.hat
		ED3[jj] <- ed3.hat
	}
	
	# coly <- rainbow_hcl(n)
	# i <- n
	# for (i in 1:n){
	#   par(mar=c(4, 4, 1.5, 1.5) +0.1)
	#   plot(ages,log(Z[,jj]/E[,jj]),t="p",col=coly[i],pch=1,main=paste("SSE -",years[i]),ylim=c(-12,0),
	#        xlab="Ages",ylab="Log(m(x))")
	#   lines(ages,LHAZhat[,i],col=1,lwd=2)
	#   lines(0:50,LHAZ1hat[,i],col=2,lty=2,lwd=2)
	#   lines(1:ages[m0],LHAZ2hat[,i],col=3,lty=2,lwd=2)
	#   lines(1:80,LHAZ3hat[,i],col=4,lty=2,lwd=2)
	#   Sys.sleep(0.1)
	# }
	
	if (!flag.stop){
		
		xINF <- x1A
		xSEN <- x2
		xACC <- x3
		xINFs <- xx1A
		xSENs <- xx2
		xACCs <- xx3
		
		return.list <- list(E=E,Z=Z,LHAZact=LHAZact,ages=ages,years=years,
				LHAZhat=LHAZhat,LHAZ1hat=LHAZ1hat,LHAZ2hat=LHAZ2hat,LHAZ3hat=LHAZ3hat,
				xINF=xINF,xSEN=xSEN,xACC=xACC,xINFs=xINFs,xSENs=xSENs,xACCs=xACCs)
	}else{
		return.list <- NULL
	}
	# "x","xs",
	# "m","m1","m2","m3","ms","m1s","m2s","m3s",
	# "n","coly","cou","sex","countries",
	# "GAMMA1hatSMO","GAMMA2hatSMO","GAMMA3hatSMO")
	
	return(return.list)
}
