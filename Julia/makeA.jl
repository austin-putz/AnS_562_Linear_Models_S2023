#==============================================================================#
# makeA4.jl
#==============================================================================#

# Austin:  Austin Putz
# Created: 2017-09-27
# Updated: 2023-02-18
# License: MIT
# Version: 0.0.3

#==============================================================================#
# Description
#==============================================================================#

# Create A with the tabular method

# Make Sure:
#   - Animals are ordered 1,2,...,n
#   - Missing parents are coded as 0

# Example: 
#ped = DataFrame( 
#	animal = [1, 2, 3, 4, 5, 6], 
#	sire   = [0, 0, 1, 1, 4, 5],
#	dam    = [0, 0, 2, 0, 3, 2]
#)

#==============================================================================#
# Function
#==============================================================================#

# Begin Function
function makeA(ped::DataFrame)

    # check eltypes (should be Int64, not a string)
    if eltype(ped[:,1]) == Int64
        @info "Gave the correct type (Int64)"
    else 
        @info "Gave wrong column type for this algorithm"
        error("Error: Need to provide Int64 numbered 1 to n")
    end

    # pull out animal, sire, dam
    a = ped[:,1]
 	s = ped[:,2]
	d = ped[:,3]

    # add 1 to the sire and dam if ordered 1,2,3,...,n
    sp1 = s .+ 1
    dp1 = d .+ 1

    # I do this because I will add an extra row and column so that
    # when parents are referenced it refers to row 1 or column 1 and
    # I won't need to write a bunch of if/then statements

    # set number of animals
	n = size(ped, 1)

    # add 1 for use in loop (will add a padding of 0's in row 1 and column 1)
	N = n + 1

    # check to make sure the pedigree is ordered 1 to n
    if a[1] != 1
        # throw an error if the first animal isn't listed as animal #1
        error("Pedigree must be sequential 1 to n")
    elseif a[n] != n
        error("Pedigree must be sequential 1 to n")
    end

    # allocate the whole A Matrix
    # This matrix will be padded with a row and column of 0's
	A = zeros(N, N)

    # Begin FOR loop
	for i in 2:N

        # calculate the diagonal element of A
        # diag = 1 + half the relationship of the sire and dam
        A[i, i] = 1 + (A[dp1[i-1], sp1[i-1]] / 2)

        # loop to calculate off-diagonal elements of A
        # Loop down a column for memory efficiency
        for j in (i+1):N
            
            # bottom left triangle of A
            # average the relationship between animal and parents of other individual
            A[j, i] = 0.5 * A[sp1[j-1], i] + 
                      0.5 * A[dp1[j-1], i]
            
            # Symmetric
            A[i, j] = A[j, i]

	    end
    end

	# return the A matrix (remove first row and column)
	return A[2:end, 2:end]

end # end function




