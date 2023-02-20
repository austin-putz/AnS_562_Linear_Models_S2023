
# Given the vectors of sire and dam, directly return inverse of additive 
# relationship matrix 'A' without creating the 'A' and the lower traingular 
# matrix 'L' themselves. This can only be applied to noninbred populations.

# Arguments
#   1) Pedigree (animal, sire, dam)

# Note:
# Unknown parents should be coded as zero.
# Does not take inbreeding into account so this can only be applied to 
# noninbred population or ignoring inbreeding on purpose

# Literature: Henderson, C. R. 1976. Simple Method for Computing the Inverse 
# of a Numerator Relationship Matrix Used in Prediction of Breeding Values. 
# Biometrics 32:69-83.

# Author: Gota Morota
# Create: 17-Apr-2009
# Last-Modified: 2-Apr-2010
# License: GPLv3 or later

# Modified by: Austin Putz

# create function
createAinv <- function(ped){
	
	if (nargs()>1){
		stop("Too many arguments")
	}
  
  # subset sire and dam from pedigree
  s = ped[, 2]
  d = ped[, 3]
	
  # check lengths
	if (length(s) != length(d)){
		stop("size of s and d is different!")
	}
	else 
	  
		n <- length(s)
		N <- n + 1
		
		# need to intially set A inverse to 0
		A <- matrix(0.0, ncol=N, nrow=N)
		
		s <- (s == 0)*(N) + s
		d <- (d == 0)*N + d			
				
	for(i in 1:n){
		
		if (s[i] != N && d[i] != N ){
			tmp        <- 2.0
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
			tmp           <- 4.0/3.0
			A[i,i]        <- A[i,i] + tmp
			A[s[i], i]    <- A[s[i], i] - tmp/2.0
			A[i, s[i]]    <- A[i, s[i]] - tmp/2.0
			A[s[i], s[i]] <- A[s[i], s[i]] + tmp/4.0
		}
			
		if (s[i] == N && d[i] != N){
			tmp           <- 4.0/3.0
			A[i,i]        <- A[i,i] + tmp
			A[d[i], i]    <- A[d[i], i] - tmp/2.0
			A[i, d[i]]    <- A[i, d[i]] - tmp/2.0
			A[d[i], d[i]] <- A[d[i], d[i]] + tmp/4.0
		}
		
		if (s[i] == N && d[i] == N){
			tmp    <- 1.0
			A[i,i] <- A[i,i] + tmp
		}									
							
	}
		
		# return the matrix
		return(A[1:n, 1:n])
		
}


