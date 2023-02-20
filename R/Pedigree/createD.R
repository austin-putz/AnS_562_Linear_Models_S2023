
## Given the vectors of sire and dam, return a dominance relationship matrix 'D'.

## Arguments
## s: a vector of sire
## d: a vector of dam 

## Note: Unknown parents should be coded as zero.

## Literature: Mrode, R.A. 2005. Linear Models for the Prediction of Animal Breeding Values. CAB International, Oxon, UK.

## Author: Gota Morota <morota at wisc dot edu>
## Create: 25-Mar-2010
## Last-Modified: 2-Apr-2010
## License: GPLv3 or later

`createD` <-
function(s, d){
	
	if (nargs()==1){
		stop("sire vector and dam vector are required")
	}
	
	if (length(s) != length(d)){
		stop("size of s and d is different!")
	}
	else 
		n <- length(s)	
		N <- n + 1
		A <- matrix(0, ncol=N, nrow=N)
		D <- matrix(0, ncol=N, nrow=N)
		
		s <- (s == 0)*(N) + s
		d <- (d == 0)*N + d
 		
	# compute A			
	for(i in 1:n){
				
		A[i,i] <- 1 + A[s[i], d[i]]/2
			
		for(j in (i+1):n){
			if (j > n) break
			A[i,j] <- ( A[i, s[j]] + A[i,d[j]] )/2
			A[j,i] <- A[i,j] 
		}
	}
	
	# compute D
	for(i in 1:n){
	  
		D[i,i] <- 1 
		
		for(j in (i+1):n){
			if (j > n) break
			D[i,j] <- ( A[s[i], d[j]]*A[d[i],s[j] ] + A[s[i],s[j]]*A[d[i], d[j]] )/4
	  	    D[j,i] <- D[i,j] 
		}
	}
		return(D[1:n, 1:n])
}


