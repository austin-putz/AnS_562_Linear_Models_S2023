
# Given the vectors of sire and dam, return a lower triangular matrix 'L' 
# such that LL' = 'A'.

# Arguments
# s: a vector of sire
# d: a vector of dam 

# Note: Unknown parents should be coded as zero.

# Literature: Henderson, C. R. 1976. Simple Method for Computing the 
# Inverse of a Numerator Relationship Matrix Used in Prediction of 
# Breeding Values. Biometrics 32:69-83.

# Author: Gota Morota
# Create: 16-Apr-2009
# Last-Modified: 1-Apr-2010
# License: GPLv3 or later

# create function
`createL` <-
function(s, d){
	
	if (nargs()==1){
		stop("sire vector and dam vector are required")
	}
	
	if (length(s) != length(d)){
		stop("size of s and d are different!")
	}
	
	n <- length(s)
	N <- n
	L <- matrix(0.0, ncol=N, nrow=N)
		
	s <- (s == 0)*(N) + s
	d <- (d == 0)*N + d
				
	for(t in 1:n){	
		
		if (s[t] == N && d[t] == N){
			L[t,t] <- 1.0
			if(t!=1) {
				for(j in 1:(t-1)){
					L[t,j] <- 0				
			    }
			}
		}
					
		if (s[t] != N && d[t] == N){
			for (j in 1:s[t]){
				L[t,j] <- 0.5*L[s[t], j] 
			}
			tmp <- 0.0
			for(j in 1:s[t]){
				tmp <- tmp + L[t,j]^2
			}
			L[t,t] <- sqrt(1- tmp)				
		}
			
		if (s[t] == N && d[t] != N){
			for (j in 1:d[t]){
				L[t,j] <- 0.5*L[d[t], j] 
			}
			tmp <- 0.0
			for(j in 1:d[t]){
				tmp <- tmp + L[t,j]^2
			}
			L[t,t] <- sqrt(1- tmp)			
		}

		if (s[t] != N && d[t] != N ){
			if (s[t] < d[t]){
				p <- s[t]
				q <- d[t]
			}
			else{
				p <- d[t]
				q <- s[t]
			}				
			for(j in 1:t-1){
				L[t,j] <- 0.5*(L[p,j] + L[q,j])
			}
			tmp <- 0.0
			for (j in 1:p){
				tmp <- tmp + L[p,j]*L[q,j]
			}			
			tmp2 <- 0.0
			for (j in 1:q){
				tmp2 <- tmp2 + L[t,j]^2
			}
			
			L[t,t] <- sqrt(1+0.5*tmp - tmp2) 	
		}
	
 	}
		
		return(L)
}


