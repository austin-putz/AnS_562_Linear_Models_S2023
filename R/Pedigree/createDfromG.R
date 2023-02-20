
## Given a gametic realationship matrix 'G', return a dominance relationship matrix 'D'.

## Argument
## G: gametic relationship matrix

## Literature: Schaeffer, L. R. 2009. Animal Models Course Note.

## Author: Gota Morota <morota at wisc dot edu>
## Create: 3-Apr-2010
## Last-Modified: 4-Apr-2010
## License: GPLv3 or later

`createDfromG` <-
function(G){
		 
	n <- dim(G)[1]/2 
	D <- matrix(0,ncol=n, nrow=n)
	
	for (i in 1:n){
		
		D[i,i] <- (G[2*i-1,2*i-1]*G[2*i,2*i] + G[2*i-1,2*i]*G[2*i,2*i-1])
        
		for (j in (i+1):n){
        		
			if (j > n) break
			D[i,j] <- (G[2*i-1,2*j-1]*G[2*i,2*j] + G[2*i-1,2*j]*G[2*i,2*j-1])
			D[j,i] <- D[i,j]
        
    	  	}
	}
	
	return(D)	
}


