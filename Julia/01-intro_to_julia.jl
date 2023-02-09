#------------------------------------------------------------------------------#
# intro_to_modeling.jl
#------------------------------------------------------------------------------#

# Author: Austin Putz
# Created: Jan 19, 2023
# Modified: Jan 20, 2023
# License: MIT

#------------------------------------------------------------------------------#
# Description
#------------------------------------------------------------------------------#

# This script will show you how to do a simple model and how to extract 
# output from such model using GLM. 

#------------------------------------------------------------------------------#
# Packages
#------------------------------------------------------------------------------#

# load Pkg package (package manager)
using Pkg

# download package if not installed
if !("DataFrames" in keys(Pkg.dependencies()))
	println("----- DataFrames does not exist and needs to be installed. Please wait...")
	Pkg.add("DataFrames")
else 
	println("----- DataFrames is downloaded already!")
end
if !("CategoricalArrays" in keys(Pkg.dependencies()))
	println("----- CategoricalArrays does not exist and needs to be installed. Please wait...")
	Pkg.add("CategoricalArrays")
else 
	println("----- CategoricalArrays is downloaded already!")
end
if !("GLM" in keys(Pkg.dependencies()))
	println("----- GLM does not exist and needs to be installed. Please wait...")
	Pkg.add("GLM")
else 
	println("----- GLM is downloaded already!")
end
if !("MixedModels" in keys(Pkg.dependencies()))
	println("----- MixedModels does not exist and needs to be installed. Please wait...")
	Pkg.add("MixedModels")
else 
	println("----- MixedModels is downloaded already!")
end
if !("UnicodePlots" in keys(Pkg.dependencies()))
	println("----- UnicodePlots does not exist and needs to be installed. Please wait...")
	Pkg.add("UnicodePlots")
else 
	println("----- UnicodePlots is downloaded already!")
end
if !("RDatasets" in keys(Pkg.dependencies()))
	println("----- RDatasets does not exist and needs to be installed. Please wait...")
	Pkg.add("RDatasets")
else 
	println("----- RDatasets is downloaded already!")
end
if !("Statistics" in keys(Pkg.dependencies()))
	println("----- Statistics does not exist and needs to be installed. Please wait...")
	Pkg.add("Statistics")
else 
	println("----- Statistics is downloaded already!")
end
if !("LinearAlgebra" in keys(Pkg.dependencies()))
	println("----- LinearAlgebra does not exist and needs to be installed. Please wait...")
	Pkg.add("LinearAlgebra")
else 
	println("----- LinearAlgebra is downloaded already!")
end

# load packages
using Random
using Distributions
using UnicodePlots
using CategoricalArrays
using MixedModels
using GLM
using DataFrames
using RDatasets
using Plots
using StatsPlots
using Statistics
using LinearAlgebra

#------------------------------------------------------------------------------#
# Load longley from R through Arrow package
#------------------------------------------------------------------------------#

# read in langley dataset
longley = RDatasets.dataset("datasets", "longley")

# print
longley

# dimensions (dim in R)
size(longley)

# print first 5 lines (head in R)
first(longley, 5)

#------------------------------------------------------------------------------#
# Model
#------------------------------------------------------------------------------#

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

# predict y-hat values
predict(lm_Emp_GNP)

# residuals into a vector
residuals(lm_Emp_GNP)

#----------------------------------------#
# F-Test - 2 Models
#----------------------------------------#

# fit:
#   y = int + x*GNP + x^2*GNP + e
lm_Emp_GNP_2 = lm(@formula(Employed ~ GNP + GNP^2), longley)

# F-test of both models
ftest(lm_Emp_GNP.model, lm_Emp_GNP_2.model)

#----------------------------------------#
# Test Categorical Model (later)
#----------------------------------------#

# convert to categorical
#categorical(longley[:, :Year])

#------------------------------------------------------------------------------#
# Simulate
#------------------------------------------------------------------------------#

#-------------------------------------------------------------#
# Simple
#-------------------------------------------------------------#

# Here we want to generate new y response vectors with current estimates of the 
# model coefficients and the residual variance. 

# extract coefficients of model
coefs = coef(lm_Emp_GNP)

# extract residual variance from model
sigma_e = stderror(lm_Emp_GNP)[1]

# calculate y_hat from coefficients extracted
y_hat = coefs[1] .+ (coefs[2] .* lm_Emp_GNP.mf.data.GNP)

# plot histogram of y-hat values
UnicodePlots.histogram(y_hat, nbins=6)

# reps
reps = 100

# sample random normal (Random and Distributions packages)
rnorm = Normal(0, sigma_e)
e_m = rand(rnorm, size(longley, 1), reps)

# sum predicted y with randomly sampled residuals to get new y vector's
# note we need .+ notation as julia does not recycle so you need to tell it to sum each 
# y_hat over each column of e_m. 
p_m = y_hat .+ e_m

# combine y variable with resampled y variable into full Matrix/DataFrame ([A B] will concat)
matrix_samples = [lm_Emp_GNP.mf.data.GNP    p_m]
data_samples   = DataFrame([lm_Emp_GNP.mf.data.GNP    p_m], :auto)

# plot new sampled y on GNP (loop and add plots)
for i in 1:reps
	if i == 1
		p = plot(matrix_samples[:,1], matrix_samples[:, 2], 
				seriestype=:scatter, mc=:gray, legend=false)
		p = Plots.title!("Samples on GNP")
		p = Plots.xlabel!("GNP")
		p = Plots.ylabel!("Sampled y")
		display(p)
	end
	j = i + 1
	p = plot!(matrix_samples[:,1], matrix_samples[:, j], 
				seriestype=:scatter, mc=:gray, legend=false)
	display(p)
end

#-------------------------------------------------------------#
# Advanced
#-------------------------------------------------------------#

# vcv WARNING THIS IS DIFFERENT THAN JUAN'S R CODE
vcv = vcov(lm_Emp_GNP)

# set nreps
nrep = 100

# choleski decomposition of vcv
C = LinearAlgebra.cholesky(vcv)
L = C.U

# simulate new coefs
n01 = rand(Normal(0, 1), nrep, 2)

# multiply
beta_sim = n01 * L

# add coefs to each column
beta_sim[:, 1] = beta_sim[:,1] .+ coefs[1]
beta_sim[:, 2] = beta_sim[:,2] .+ coefs[2]

# head beta_sim
first(beta_sim, 5)

# get covariance of beta_sim
cov(beta_sim)

# column means of beta_sim
Statistics.mean(beta_sim, dims=1)

# plot with regression lines simulated
for i in 1:nrep
	if i == 1
		global p = scatter(matrix_samples[:,1], matrix_samples[:, 2], 
				mc=:red, legend=false)
		p = Plots.title!("Samples on GNP")
		p = Plots.xlabel!("GNP")
		p = Plots.ylabel!("Employed")
		display(p)
	end
	# set up intercept and slope that was sampled
	int   = beta_sim[i, 1]
	slope = beta_sim[i, 2]
	Plots.abline!(slope, int, lc=:gray)
	display(p)
end













