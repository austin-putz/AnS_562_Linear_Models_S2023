## Given a gametic realationship matrix, return an additive relationship matrix 'A'.

## Argument
## G: gametic relationship matrix.

## Literature: Schaeffer, L. R. 2009. Animal Models Course Note.

## Author: Gota Morota <morota at wisc dot edu>
## Create: 3-Apr-2010
## Last-Modified: 4-Apr-2010
## License: GPLv3 or later

`createAfromG` <-
function(G){
		 
	n <- dim(G)[1]/2 
	A <- matrix(0,ncol=n, nrow=n)
	
	for (i in 1:n){
		
        A[i,i] <- (G[2*i-1,2*i-1] + G[2*i-1,2*i] + G[2*i,2*i-1] + G[2*i,2*i])/2
        
        for (j in (i+1):n){
        		
			if (j > n) break
			A[i,j] <- (G[2*i-1,2*j-1] + G[2*i-1,2*j] + G[2*i,2*j-1] + G[2*i,2*j])/2
			A[j,i] <- A[i,j]
        
		}
	}

	return(A)	
}


