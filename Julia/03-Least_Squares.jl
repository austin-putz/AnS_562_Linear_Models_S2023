#--------------------------------------------------------------------------------#
# 03-Linear_Models.jl
#--------------------------------------------------------------------------------#

# Linear Models:
# Lab 03 for AnS 562 Linear Models Class

# Author: Juan Steibel and Austin Putz
# Created: Feb 07, 2023
# Modified: Feb 18, 2023
# Version: 0.0.3
# License: MIT

# Original R code by Juan Steibel, re-written in Julia by Austin Putz. 

#------------------------------------------------------------------------------#
# Description:
#------------------------------------------------------------------------------#

# In this file we will cover basic linear models in terms of least squares
# and show you how to solve the equations by hand. 

#------------------------------------------------------------------------------#
# Packages
#------------------------------------------------------------------------------#

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
my_packages = ["RDatasets", "LinearAlgebra", "StatsBase", 
               "CSV", "DataFrames", "IterativeSolvers"];

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

# load packages
using RDatasets
using LinearAlgebra
using StatsBase
using CSV
using DataFrames
using IterativeSolvers

#--------------------------------------------------------------------------------#
# Dataset
#--------------------------------------------------------------------------------#

# load langley dataset from RDatasets package
longley = RDatasets.dataset("datasets", "longley")

#--------------------------------------------------------------------------------#
# Solve Linear Equations
#--------------------------------------------------------------------------------#

# print regressor variable (:GNP is a symbol we can use for column names)
longley[:, :GNP]

# generate y vector (Employed)
y = longley[:, :Employed]

# generate X matrix
X = [ones(size(longley, 1), 1) longley[:,:GNP]]

# rank of X
rank(X)

# X'X Matrix
XpX = X'X

# (X'X)^1 (inverse of X'X)
XpXi = inv(XpX)

# X'y
Xpy = X'y

# solve for b-hat
bhat = XpXi * Xpy

# y-hat = X * bhat
yhat = X * bhat

# e-hat (y - y_hat) because y = y_hat + e
ehat = y .- yhat

# model residual variance (2nd part is the 'degrees of freedom')
sigma2_e_hat = sum(ehat.^2) / (size(X, 1) - size(X, 2))

# sqrt that variance
sqrt(sigma2_e_hat)

# (co)variance matrix of coefficients (variance in estimates)
var_hat_beta_hat = XpXi * sigma2_e_hat

# equivalent to vcov() of the model later...

# calculate SE of coef estimates
sqrt.(diag(var_hat_beta_hat))
# or
diag(var_hat_beta_hat).^0.5

#--------------------------------------------------------------------------------#
# Use GLM Package
#--------------------------------------------------------------------------------#

# We now want to simply use the lm

# load GLM package
using GLM

# Fit linear model:
#   y = Xb + e
#   (Employed on GNP)
lm_Emp_GNP = lm(@formula(Employed ~ GNP), longley)

# equivalent to this:
lm_Emp_GNP = fit(LinearModel, 
                    @formula(Employed ~ GNP), 
                    longley)

# the latter is more 'Julia'

# print variables within
propertynames(lm_Emp_GNP)

# print table with output (coef, std error, t, p-val, confid interval)
lm_Emp_GNP.model

# extract X (model) matrix (within mm object)
X = lm_Emp_GNP.mm.m

# calc Rank of X
rank(X)

# vectors of y and x (column)
lm_Emp_GNP.mf.data.Employed
lm_Emp_GNP.mf.data.GNP

# print model (variables fit)
lm_Emp_GNP.mf.f

#----------------------------------------#
# Extract from model fit
#----------------------------------------#

# coefficients
coef(lm_Emp_GNP)

# predict y-hat values
predict(lm_Emp_GNP)

# residuals into a vector
residuals(lm_Emp_GNP)

# deviance
deviance(lm_Emp_GNP)  # 6.036

# Residual DF (convert to Int64 after calculation)
dfresid = Int64.(dof_residual(lm_Emp_GNP)) # 14

# r2
r2(lm_Emp_GNP)

# std error of coefficients
stderror(lm_Emp_GNP)

# vcov - (co)variance matrix of model estimates for coefficients
vcov(lm_Emp_GNP)

# get Standard Errors the hard way
diag(vcov(lm_Emp_GNP)).^0.5

# rank of X
rank(lm_Emp_GNP.mm.m)

# use predict and calculate yourself the y-hat values
[predict(lm_Emp_GNP)   X*bhat]

# use residuals and calculate e-hat yourself
[residuals(lm_Emp_GNP)   y .- yhat]

# test if X'X is positive def
minimum(eigvals(Hermitian(X'X))) > 0
# Or:
isposdef(X'X)
# Or:
all(eigen(X'X).values .> 0)

#--------------------------------------------------------------------------------#
# QR Factorization to Solve
#--------------------------------------------------------------------------------#

# qr factorization of X
qrX = qr(X)

# access Q and R
qrX.Q
qrX.R

# QR factorization of X
Q, R1 = qr(X)
Q1 = Q[:, 1:2]
# OR use:
Q1 = Matrix(qrX.Q)
Q2 = Q[:, 3:end]
Rbottom = zeros(14, 2)
R = vcat(R1, Rbottom)

# re-create X
X
Q * R1

# re-create X'X
X'X
R1' * R1

# re-create (X'X)
inv(X'X)
inv(R1) * inv(R1')

# re-create b-hat
inv(X'X) * X'y
inv(R1) * Q1' * y

# load BenchmarkTools
using BenchmarkTools

# calc b-hat with inverse
@benchmark inv(X'X) * X'y

# calc b-hat with QR decomposition
# we must account for the QR decomposition first as well to compare equally
@benchmark Q, R1 = qr(X)
@benchmark Q1 = Q[:, 1:2]
@benchmark inv(R1) * Q1' * y

# test that Q has orthogonal columns (inner products are = 0)
round.(Q' * Q)

#--------------------------------------------------------------------------------#
# Standardizing
#--------------------------------------------------------------------------------#

# From StatsBase package

# standardize X
sX = fit(ZScoreTransform, X, dims=1)

# transform back
StatsBase.transform(sX, X)

# scale predictor
sy = fit(ZScoreTransform, longley[:, :GNP], dims=1)
StatsBase.transform(sy, longley[:, :GNP])

#--------------------------------------------------------------------------------#
# Simulate
#--------------------------------------------------------------------------------#

# load Random package
using Random

function sample_data(N, p)
    # sample X
    X = randn(N, p)
    # sample y
    y = randn(N, 1)
    # return values
    return X, y
end

# return X and y from sample_data() function
X, y = sample_data(100, 10)

# you can also sample random normals this way
# rand(Normal(0,1), 100, 10)

# Benchmark different solvers

#--------------------------------------------------------------------------------#
# npk Dataset
#--------------------------------------------------------------------------------#

# load NPK dataset
npk = RDatasets.dataset("MASS", "npk")

# fit model
mm1 = lm(@formula(Yield ~ N + P + K), npk)

# print model
print(mm1)

#--------------------------------------------------------------------------------#
# load read dataset
#--------------------------------------------------------------------------------#

# load swine data
dts = CSV.read("/Users/austinputz/Documents/ISU/Classes/AnS_562/2023/Part_A/Julia/swine_data.csv", DataFrame,
        header=true, delim=',', missingstring="NA")

# replace missing
#dts.car_wt = map(x -> x == "NA" ? missing : x, dts.car_wt)
#dts.num_ribs = map(x -> x == "NA" ? missing : x, dts.num_ribs)
#dts.car_bf10 = map(x -> x == "NA" ? missing : x, dts.car_bf10)

#--------------------------------------------------------------------------------#
# Iterative Solvers 
#--------------------------------------------------------------------------------#

# Iterative solvers may be used when the solution is too large and we cannot
# do direct inversion of X'X or do the QR solve as above. We can simply setup
# X'X and X'y and solve the linear equation with these solvers. 

# load package
#using IterativeSolvers

# Must give these X'X and X'y to work

# solve with conjugate gradient solver
cgX = cg(X'X, X'y)

# or with GMRES
gmresX = gmres(X'X, X'y)











