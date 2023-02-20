#--------------------------------------------------------------------------------#
# 05-Relationship_Matrices.jl
#--------------------------------------------------------------------------------#

# Relationship Matrices:
# Lab 05 for AnS 562 Linear Models Class

# Author: Juan Steibel and Austin Putz
# Created: Feb 18, 2023
# Modified: Feb 18, 2023
# Version: 0.0.3
# License: MIT

# Original code by Austin Putz. 

#------------------------------------------------------------------------------#
# Description:
#------------------------------------------------------------------------------#

# In our LHS of the mixed model equations:
# 
#         | X'X    X'Z                 |
#         | Z'X    Z'Z + inv(A)*lambda | 
# 
# We need this very important matrix _A_ and namely it's inverse only to solve
# the MME. 
# 
# We can either set up A directly and invert it at the cost of O(n^3)
# or we can setup inv(A) directly. We'll do some of both below. 
# 

#--------------------------------------------------------------------------------#
# Packages
#--------------------------------------------------------------------------------#

# If you don't have the following packages downloaded:
#      using Pkg
#      Pkg.add("PackageName")
# OR
#    go into the package manager by clicking ']' (right bracket) in REPL (interactive mode)
#    then:
#      pkg> add PackageName

# load Pkg (package manager for Julia)
using Pkg

# list packages to install if not installed
my_packages = ["Pipe", "LinearAlgebra", 
               "CSV", "DataFrames", "DataFramesMeta",
               "StatsBase", "FreqTables", 
               "Plots", "UnicodePlots",
			   "SparseArrays"];

# list current installed packages
cur_packages_installed = [ p.name for p in values(Pkg.dependencies()) ];

# loop over to install packages or load them
for package in my_packages
	if !(package in cur_packages_installed)
		@info "Installing: $package"
		Pkg.add(package)
	else
		@info "Package $package is installed"
		@eval using $(Symbol(package))
	end
end

#--------------------------------------------------------------------------------#
# Pedigree
#--------------------------------------------------------------------------------#

# create small pedigree
ped = DataFrame( 
	animal = [1, 2, 3, 4, 5, 6], 
	sire   = [0, 0, 1, 1, 4, 5], 
	dam    = [0, 0, 2, 0, 3, 2]
)

#--------------------------------------------------------------------------------#
# Create A with the tabular method
#--------------------------------------------------------------------------------#

# load makeA.jl
include("/Users/austinputz/Documents/Programming/Julia/Animal_Breeding/Pedigrees/makeA.jl")

# make A
A = makeA(ped)

# inbreeding
F = diag(A) .- 1

# inverse of A
Ainv = inv(A)

# make sparse (SparseArrays package)
sparse(Ainv)






