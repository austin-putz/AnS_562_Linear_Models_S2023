################################################################################
# Functions
################################################################################

# Functions ---------------------------------------------------------------

# Given the pedigree, return additive relationship matrix 'A'.

# Arguments
#		1) ped
#			ped needs to be in the format
#				column 1 = animal
# 			column 2 = sire
#				column 3 = dam

####### NOTE: Unknown parents should be coded as ZERO (NOT 'NA')

# Literature: Henderson, C. R. 1976. 
# Simple Method for Computing the Inverse of a Numerator Relationship 
# Matrix Used in Prediction of Breeding Values. Biometrics 32:69-83.

# Author: Gota Morota <morota at wisc dot edu>
# Create: 16-Apr-2009
# Last-Modified: 1-Apr-2010
# License: GPLv3 or later

	`createA` <-
	  
		function(ped){
		  
		  # load compiler package
		  #library(compiler)
		  #enableJIT(3)
		  
			if (nargs() > 1 ) {
				stop("Only the pedigree is required (Animal, Sire, Dam)")
			}
			
			# This is changed from Gota's function
			s = ped[, 2]
			d = ped[, 3]
			
			# check for NA
			
			# Stop if they are different lengths
			if (length(s) != length(d)){
				stop("\nSize of the sire vector and dam vector are different!")
			}
			
			# set number of animals and empty vector
			n <- length(s)
			N <- n + 1
			A <- matrix(0, ncol=N, nrow=N)
			
			# set sires and dams
			s <- (s == 0)*(N) + s
			d <- (d == 0)*N + d
			
			# start time
			start_time <- Sys.time()
			
			# Begin for loop
			for(i in 1:n){
				
			  # diagonal elements of A
				A[i,i] <- 1 + A[s[i], d[i]]/2
				
				# loop over j
				for(j in (i+1):n){
					if (j > n) break    # only do half of the matrix (symmetric)
					A[i,j] <- ( A[i, s[j]] + A[i, d[j]] ) / 2  # half relationship of parents
					A[j,i] <- A[i,j] 	# symettric matrix
				}
			}
		
    cat("\t", sprintf("%-30s:%f", "Time it took (sec)", as.numeric(Sys.time() - start_time)), "\n")
			
			return(A[1:n, 1:n])
      
		  enableJIT(0)
		  
		}

