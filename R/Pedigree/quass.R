
# Given the vectors of sire and dam, 
# directly return inverse of additive relationship matrix 'A' 
# without creating the 'A' itself. This is a modification of Henderson's 
# method and unlike createAinv.r, this can be used in inbred populations. 

# Arguments
# s: a vector of sire
# d: a vector of dam 

# Note 
# Unknown parents should be coded as zero.

# Literature: Quass, R. L. 1976. Computing the Diagonal Elements and 
# Inverse of a Large Numerator Relationship Matrix. Biometrics 32:949-953.

# Author: Gota Morota
# Create: 17-Apr-2009
# Last-Modified: 2-Apr-2010
# License: GPLv3 or later

`quass` <-
function(s, d){
	
	if (nargs() == 1){
		stop("sire vector and dam vector are required")
	}
	
	if (length(s) != length(d)){
		stop("size of s and d is different!")
	}
	else 
		n <- length(s)
		N <- n
		L <- matrix(0.0, ncol=N, nrow=N)
		
		s <- (s == 0)*(N) + s
		d <- (d == 0)*N + d
 				
  # construct L				
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
			
			L[t,t] <- sqrt(1 + 0.5*tmp - tmp2) 	
		}
	
	}	
  
  #========================================#
  # calculate A inverse based on L
  #========================================#
  
	# allocate the A matrix
	A <- matrix(0.0, ncol=N, nrow=N)
		
	for(i in 1:n){		
			
		tmp <- 1/L[i,i]^2	
		
		if (s[i] != N && d[i] != N ){
		  
			A[i,i]     <- A[i,i] + tmp
			A[i, s[i]] <- A[i, s[i]] - tmp/2.0
			A[s[i], i] <- A[s[i], i] - tmp/2.0
			A[i,d[i]]  <- A[i, d[i]] - tmp/2.0
			A[d[i], i] <- A[d[i], i] - tmp/2.0
			
			A[s[i], s[i]] <- A[s[i], s[i]] + tmp/4.0
			A[s[i], d[i]] <- A[s[i], d[i]] + tmp/4.0
			A[d[i], s[i]] <- A[d[i], s[i]] + tmp/4.0
			A[d[i], d[i]] <- A[d[i], d[i]] + tmp/4.0
			
		}
			
		if (s[i] != N && d[i] == N){
		  
			A[i,i]        <- A[i,i] + tmp
			A[s[i], i]    <- A[s[i], i] - tmp/2.0
			A[i, s[i]]    <- A[i, s[i]] - tmp/2.0
			A[s[i], s[i]] <- A[s[i], s[i]] + tmp/4.0
			
		}
			
		if (s[i] == N && d[i] != N){
		  
			A[i,i]        <- A[i,i] + tmp
			A[d[i], i]    <- A[d[i], i] - tmp/2.0
			A[i, d[i]]    <- A[i, d[i]] - tmp/2.0
			A[d[i], d[i]] <- A[d[i], d[i]] + tmp/4.0
			
		}
		
		if (s[i] == N && d[i] == N){
			A[i,i] <- A[i,i] + tmp
		}								
							
	}
		
		# return A
		return(list(L=L, Ainv=A))

}


