
## Given the vectors of sire and dam with/without marker information, 
## return a gametic relationship matrix 'G'.

## Arguments
## s:  a vector of sire
## d:  a vector of dam 
## sM: a vector of marker inheritance for sire. 0=unknown, 1=father, 2=mother
## dM: a vector of marker inheritance for dam. 0=unknown, 1=father, 2=mother
## r:  recombination rate between marker locus and marked QTL

## Note
## Unknown parents should be coded as zero.
## You can omit last three arguments if you don't have any marker information.

## Literature 1: Fernando, R. L. and Grossman, M. 1999. Marker-Assisted 
##                Selection Using Best Linear Unbiased Prediction. Genetic 
##                Selection Evolution 21, 467-477.
## Literature 2: Schaeffer, L. R. 2009. Animal Models Course Note.
## Literature 3: Mrode, R.A. 2005. Linear Models for the Prediction of 
##                Animal Breeding Values. CAB International, Oxon, UK.

## Author: Gota Morota <morota at wisc dot edu>
## Create: 3-Apr-2010
## Last-Modified: 4-Apr-2010
## License: GPLv3 or later

`createG` <-
function(s, d, sM, dM, r){
	
	if (nargs() < 2){
		stop("at least sire vector and dam vector are required")
	}
	
	if (length(s) != length(d)){
		stop("size of s and d is different!")
	}
	
	if (nargs()==2){
		sM <- rep(0, length(s))
		dM <- rep(0, length(s))
	}

	n <- length(s)
	s <- (s == 0)*n + s
	d <- (d == 0)*n + d
	G <- matrix(0, ncol=2*(n+1), nrow=2*(n+1))
	diag(G) <- 1

	for (i in 1:n){
		# off-diagonal elementss of a block
		
		# Paternal
		if (sM[i] == 0){
			R <- 0.5	
		}
		if (sM[i] == 1){
			R <- r	
		}
		if (sM[i] == 2){
			R <- (1- r)
		}
		G[2*i, 2*i-1] <-(1-R)*G[2*i,2*s[i]-1] + R*G[2*i,2*s[i] ] 
		
		# Maternal
		if (dM[i] == 0){
			R <- 0.5	
		}
		if (dM[i] == 1){
			R <- r	
		}
		if (dM[i] == 2) {
			R <- 1-r	
		}
		G[2*i-1, 2*i] <- (1-R)*G[2*i-1,2*d[i]-1] + R*G[2*i-1,2*d[i]]
	
		#############j starts#################
		for (j in (i+1):n){
			
			if (j > n) break
			
			# Paternal
			if (sM[j] == 0){
				R <- 0.5	
			}
			if (sM[j] == 1){
				R <- r	
			}
			if (sM[j] == 2) {
				R <- 1- r	
			}
			G[2*i-1,2*j-1] <- (1-R)*G[2*i-1,2*s[j]-1] + R*G[2*i-1,2*s[j]]
			G[2*j-1,2*i-1] <- G[2*i-1,2*j-1]
		
			
			# Maternal
			if (dM[j] == 0){
				R <- 0.5	
			}
			 if (dM[j] == 1){
				R <- r	
			}
			if (dM[j] == 2){
				R <- 1-r	
			}
			G[2*i-1,2*j ] <- (1-R)*G[2*i-1,2*d[j]-1] + R*G[2*i-1,2*d[j]]
			G[2*j,2*i-1] <- G[2*i-1  ,2*j]
		
			# Paternal
			if (sM[j] == 0){
				R <- 0.5	
			}
			if (sM[j] == 1){
				R <- r	
			}
			if (sM[j] == 2){
				R <- 1- r	
			}
			G[2*i  ,2*j-1] <- (1-R)*G[2*i, 2*s[j]-1] + R*G[2*i  ,2*s[j]]
			G[2*j-1  ,2*i] <- G[2*i,2*j-1]
		
			# Maternal
			if (dM[j] == 0){
				R <- 0.5	
			}
			if (dM[j] == 1){
				R <- r	
			}
			if (dM[j] == 2){
				R <- 1-r	
			}
			G[2*i  ,2*j] <- (1-R)*G[2*i  ,2*d[j]-1] + R*G[2*i  ,2*d[j]]
			G[2*j  ,2*i] <- G[2*i  ,2*j]

		}
	}
	
	return(G[1:(2*n),1:(2*n)])
}






