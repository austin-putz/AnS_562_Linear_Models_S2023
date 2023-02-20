#------------------------------------------------------------------------------#
# File: linear_algebra_intro.jl
#------------------------------------------------------------------------------#

# Linear Algebra:
# Lab 02 for AnS 562 Linear Models Class

# Author: Austin Putz
# Created: Jan 23, 2023
# Modified: Feb 18, 2023
# Version: 0.0.3
# License: MIT

#------------------------------------------------------------------------------#
# Description:
#------------------------------------------------------------------------------#

#  - This file introduces students to linear algebra in Julia.
#  - Topics:
#     - addition
#     - subtraction
#     - multiplication
#     - scalers
#     - vectors
#     - parallel computing
#     - inverses
#     - determinant
#     - inverses
#     - much more

#------------------------------------------------------------------------------#
# Background
#------------------------------------------------------------------------------#

# Linear algebra became common place in animal breeding in the 1960's and 70's
# when Shayle Searle started integrating linear (matrix) algebra into
# animal breeding with Charles (Chuck) Henderson at Cornell.

# Linear models can be used for solving animal models for EBVs and many other
# topics in animal breeding such as optimal contribution selection, 
# economic selection index, genomics, and many other topics.

# Matrices are 2 dimensional Arrays.
# They are unique from data frames in that all elements within the Matrix is
# of the same type, while data frames can have heterogeneous types within. 

# Much of this can be found here:
# https://docs.julialang.org/en/v1/manual/arrays/

#------------------------------------------------------------------------------#
# Packages
#------------------------------------------------------------------------------#

# load Pkg (package manager for Julia)
using Pkg

# list packages to install if not installed
my_packages = ["LinearAlgebra", "SparseArrays", "Distributions"]

# list current installed packages
cur_packages_installed = [ p.name for p in values(Pkg.dependencies()) ]

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

# load packages
using LinearAlgebra
using SparseArrays
using Distributions

#------------------------------------------------------------------------------#
# Setup
#------------------------------------------------------------------------------#

# print basic julia and computer
versioninfo()

# print info on variables - equivalent to ls() in R with more info
varinfo()

#------------------------------------------------------------------------------#
# Matrix functions
#------------------------------------------------------------------------------#

# define matrix A
A = [
 1 5
 2 10
 ]

# print A
A

# ask for type of A
typeof(A)
#STOUT: Matrix{Int64} (alias for Array{Int64, 2})

# need eltype() to determine type of elements within an array
eltype(A)
#STOUT: Int64

# number of total elements in A
length(A)

# number of dimensions
ndims(A)

# dimension of A
size(A)

# number of rows
size(A, 1)

# number of columns
size(A, 2)

# tuple containing the valid indices of A
axes(A)

# get index of each - use for loops
eachindex(A)

# if you want to understand types better, you can view
# the string of bits with bitstring()
bitstring(Int8.(1))
bitstring(Int32.(1))
bitstring(Int64.(1))

#------------------------------------------------------------------------------#
# Initializing Arrays
#------------------------------------------------------------------------------#

# initialize array - integers
Array{Int8}(undef, 2, 2)
Array{Int16}(undef, 2, 2)
Array{Int32}(undef, 2, 2)
Array{Int64}(undef, 2, 2)
Array{Int128}(undef, 2, 2)

# UInt is unsigned integers
Array{UInt64}(undef, 2, 2)

# initialize array - floats
Array{Float16}(undef, 2, 2)
Array{Float32}(undef, 2, 2)
Array{Float64}(undef, 2, 2)

# initalize with zeros() function
zeros(Int64, 5, 5)
zeros((2,3))
zeros(Float32, (2,3))

# ones() to initialize matrix
ones(Int64, 5, 5)

# trues
trues(2, 2)

# reshape A to a [1 x 4] Matrix
reshape(A, 1, 4)

# initialize a random matrix of floats
rand(5, 5)

# sample random normal variables in 5x5 matrix
rand(Normal(0,1), 5, 5)

# initialize Matrix
W = Matrix{Int64}(undef, 2, 4)

# range - seq from 1 to 5
range(1, 5)

# fill W with 2
fill!(W, 2)

# fill a new matrix with 8's [2, 5]
fill(8, 2, 5)

# initialize column vector
[1, 2, 3]

# initialize row vector
[1 2 3]

# promote 1st integer to Float64
promote(1, 2.3)

#------------------------------------------------------------------------------#
# Basic Operations
#------------------------------------------------------------------------------#

# print A
A

# define matrix B
B = [
 3 5
 4 6
 ]

# define matrix C
C = [
  2 4 6
  3 5 7
  ]

# sum matrices element-by-element
A + B

# subtract matrices element-by-element
A - B

# concatenate A and B horizontally
[A B]
hcat(A, B)

# concatenate vertically
[A; B]
vcat(A, B)

# concatenate horizontally with vectors
[1:2 4:5 6:7]

#------------------------------------------------------------------------------#
# Comprehensions
#------------------------------------------------------------------------------#

# initialize x
x = rand(8)

# weighted average of last element, current, next element in vector
[ 0.25*x[i-1] + 0.5*x[i] + 0.25*x[i+1] for i=2:length(x)-1 ]

#------------------------------------------------------------------------------#
# Generator Expressions
#------------------------------------------------------------------------------#

# sum all from 1 to 100
sum(n for n=1:100)

# using map
map(tuple, (1/(i+j) for i=1:2, j=1:2), [1 3; 2 4])

# example with allele frequencies
p = [0.25, 0.5, 0.75, .9]
q = 1 .- p

# sum(2pq) like in G matrix calculation
sum(2 * p[i] * q[i] for i=1:length(p))

# double loop
[(i,j) for i=1:3 for j=1:i]

# insert if statement
[(i,j) for i=1:3 for j=1:i if i+j == 4]


#------------------------------------------------------------------------------#
# Indexing (slicing)
#------------------------------------------------------------------------------#

# you can slice with bracket [ ] notation

# extract 1st row
A_row_1 = A[1,:]

# extract 1st column
A[:,1]

# one can also use the special 'end' variable in slicing

# extract last 2 elements in 1st row
C[1, 2:end]

# you can also return a matrix with this notation (not a vector)
A[[1],:]

# notice other operations before returned a vector

# convert seq to vector with collect (Vector{Int64})
collect(1:16)

# reshape vector into matrix
reshape(collect(1:16), (4, 4))

# 3 dimensional array
reshape(collect(1:27), (3, 3, 3))

# set U with uniform scaling
u = UniformScaling(2)

# set A
M = [1 2; 3 4]

# add u to diagonal elements
M + u

# multiply all elements by u
M * u

#------------------------------------------------------------------------------#
# Advanced Operations
#------------------------------------------------------------------------------#

# https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/

# random matrix
X = rand(5, 5)

# trace
tr(X)

# determinant of X
det(X)

# log of determinant (here is negative so we get an error)
log(det(X))

# or use function logdet (negative so we get an error)
logdet(X)

# inverse of X
inv(X)

# pseudoinverse (generalize inverse)
pinv(X)

# extract only eigenvalues
eigenX.values

# eigenvectors only
eigenX.vectors

# xtract into variables directly
Xvals, Xvecs = eigen(X)

# eigmax() and eigmin() but not working

# eigen values and vectors
eigvals(X)
eigvecs(X)

# SVD = singular value decomposition of X
u, s, v = svd(X)

# factorize X
factorize(X)

# create A
A = [4. 12. -16.; 12. 37. -43.; -16. -43. 98.]

# cholesky factorization
C = cholesky(A)

# extract parts of cholesky
C.U
C.L

# re-create A
C.L * C.U == A

# extract into variables directly
L, U = C;

# LU factorization
lu(X)

# QR factorization
qr(X)

# create S
S = SymTridiagonal([3., 4., 5.], [1., 2.])

# LDLt
ldlt(S)

# rank of S
rank(S)

# ranks
rank(Matrix(I, 3, 3))

# diagm will create a diagonal matrix
rank(diagm(0 => [1, 0, 2]))
rank(diagm(0 => [1, 0.001, 2]), rtol=0.1)
rank(diagm(0 => [1, 0.001, 2]), rtol=0.00001)
rank(diagm(0 => [1, 0.001, 2]), atol=1.5)

# create v
v = [3, -2, 6]

# norm of v (p = 2)
norm(v)

# norm of vector
norm([1 2 3 4 5 6 7 8 9])

# make a
a = [1,2,4];

# normalize a
normalize(a)

# create a
a = [1 2 ; 5 9]

# condition number of matrix
cond(a)

# create A
A = [1 2; 3 4]

# create B
B = [im 1; 1 -im]

# kronecker product
kron(A, B)

#------------------------------------------------------------------------------#
# Other Matrices
#------------------------------------------------------------------------------#

# create A
A = [1 0 2 0 3; 0 4 0 5 0; 6 0 7 0 8; 0 9 0 1 0; 2 0 3 0 4]

# make Symmetric
Aupper = Symmetric(A, :U)
Alower = Symmetric(A, :L)

# create A
A = [1 0 2+2im 0 3-3im; 0 4 0 5 0; 6-6im 0 7 0 8+8im; 0 9 0 1 0; 2+2im 0 3-3im 0 4];

# hermitian 
Hupper = Hermitian(A)

# hermitian
Hlower = Hermitian(A, :L)

# need to use Hermitian to get real eigenvalues and vectors with eigen!!

# create A
A = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]

# lower triangular
LowerTriangular(A)

# upper triangular
UpperTriangular(A)

#------------------------------------------------------------------------------#
# Sparse Matrices
#------------------------------------------------------------------------------#

# identity matrix
I5 = I(5)

# print
I5

# multiply identity matrix by 5
5 * I5

# create new diagonal matrix
(0.7*I)(3)

# create diagonal matrix
Diagonal([5, 6, 7, 8])

# diagonal matrix dense
diagm([7, 13])

# create A
A = permutedims(reshape(1:15, 5, 3))

# create diagonal matrix from A
Diagonal(A)

# extract diagonals
diag(A)

# extract 1 above diagonal
diag(A, 1)

# create A
A = [2 1 1; 1 2 0; 1 0 2]

# make sparse
sparse(A)

# create A
A = sparse([1,2,3,4], [1,1,2,2], [1.0,1.0,1.0,1.0])

#------------------------------------------------------------------------------#
# Understanding memory management in Julia
#------------------------------------------------------------------------------#

# Julia is very memory efficient, part of this is not making
# copies of matrices and subsets of matrices in memory (RAM).
# This is good for memory management but bad if you forget
# this happens in Julia. See examples below:

# define matrix D
D = [
	2 5 9
	1 9 3
	3 7 2
	]

# assign D to E
E = D

# change one element of E
E[1,1] = 100

# print
E
D

# Note:
#   Now both D and E have element [1,1] changed to 100
#   You have to be very careful altering elements in a matrix or vector
#   Make sure you understand when something is copied and when
#   it is not in memory (RAM).

# to make a copy in memory of another matrix
F = copy(D)

# alter element [1,1]
F[1,1] = 50

# print F
F
D






















