#--------------------------------------------------------------------------------#
# Lab_02
#--------------------------------------------------------------------------------#

# If you don't have the following packages downloaded:
#  using Pkg
#  Pkg.add("PackageName")
# OR
#  go into the package manager by clicking ']' (right closed bracket) in REPL (interactive mode)
#  then:
#    pkg> add PackageName

# load packages
using RDatasets
using LinearAlgebra

#--------------------------------------------------------------------------------#
# Dataset
#--------------------------------------------------------------------------------#

# read in langley dataset from RDatasets package
longley = RDatasets.dataset("datasets", "longley")

#--------------------------------------------------------------------------------#
# Solve Linear Equations
#--------------------------------------------------------------------------------#

# print regressor variable (:GNP is a symbol we can use for column names)
longley[:, :GNP]

# generate y vector
y = longley[:, :Employed]

# generate X matrix
X = [ones(size(longley, 1), 1) longley[:,:GNP]]

# X'X Matrix
XpX = X'X

# (X'X)^1 (inverse of X'X)
XpXi = inv(XpX)

# X'y
Xpy = X'y

# solve for b-hat
bhat = XpXi * Xpy

# y-hat
yhat = X * bhat

# e-hat
ehat = y .- yhat

# model residual variance
sigma2_e_hat = sum(ehat.^2) / (size(X, 1) - size(X, 2))

# sqrt that variance
sqrt(sigma2_e_hat)

# (co)variance matrix of coefficients (variance in estimates)
var_hat_beta_hat = XpXi * sigma2_e_hat

# calculate SE of coef estimates
diag(var_hat_beta_hat).^0.5

#--------------------------------------------------------------------------------#
# Use GLM Package
#--------------------------------------------------------------------------------#

# load GLM package
using GLM

# Fit linear model:
#   y = Xb + e
#   (Employed on GNP)
lm_Emp_GNP = lm(@formula(Employed ~ GNP), longley)

# print variables within
propertynames(lm_Emp_GNP)

# print table with output (coef, std error, t, p-val, confid interval)
lm_Emp_GNP.model

# extract X (model) matrix (within mm object)
lm_Emp_GNP.mm.m

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

# load StatsBase
using StatsBase

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


#--------------------------------------------------------------------------------#
# Iterative Solvers 
#--------------------------------------------------------------------------------#

# load package
using IterativeSolvers

# solve with conjugate gradient solver
cgX = cg(X'X, X'y)

# or with GMRES
gmresX = gmres(X'X, X'y)











