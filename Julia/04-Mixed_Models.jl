#--------------------------------------------------------------------------------#
# 04-Animal_Models.jl
#--------------------------------------------------------------------------------#

# Animal Models:
# Lab 04 for AnS 562 Linear Models Class

# Author: Juan Steibel and Austin Putz
# Created: Feb 07, 2023
# Modified: Feb 18, 2023
# Version: 0.0.3
# License: MIT

# Original R code by Juan Steibel, re-written in Julia by Austin Putz. 

#------------------------------------------------------------------------------#
# Description:
#------------------------------------------------------------------------------#

# In this file we will cover basic animal models y = Xb + Zu + e
# We'll do some testing of models as well. 

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
my_packages = ["Pipe", "LinearAlgebra", "RDatasets", 
               "CSV", "DataFrames", "DataFramesMeta",
               "StatsBase", "FreqTables", 
               "GLM", "MixedModels", "Effects", 
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
# Read Data
#--------------------------------------------------------------------------------#

# directory and data file name
working_dir = "/Users/austinputz/Documents/ISU/Classes/AnS_562/2023/Part_A/Julia/"
data_file   = "swine_data.csv"

# load swine data
dts = CSV.read(working_dir * data_file, 
                DataFrame,
                header=true, 
                delim=',', 
                missingstring="NA")

#----------------------------------------#
# Change columns
#----------------------------------------#

# change scale of birth weight trait
dts.wt_birth_h = dts.wt_birth * 10

# convert to string (will fit as categorical later)
dts.ID       = string.(dts.ID)
dts.sex      = string.(dts.sex)
dts.sire     = string.(dts.sire)
dts.dam      = string.(dts.dam)
dts.litter   = string.(dts.litter)
dts.slgdt_cd = string.(dts.slgdt_cd)

# print describe of dataset to get details on all columns
describe(dts)

#----------------------------------------#
# Show how to get tables
#----------------------------------------#

# count sires
countmap(dts.sire)

# count dam
countmap(dts.dam)

# or count sires with DataFrames:
dts_group_sire = groupby(dts, :sire);
combine(dts_group_sire, nrow)

# or count sires with FreqTables:
freqtable(dts.sire)

# 2 way table - Sire & Sex
freqtable(dts.sire, dts.sex)

# unique sires
unique(dts.sire)

# unique dams
unique(dts.dam)

# histogram of birth wt
UnicodePlots.histogram(dts.wt_birth_h, nbins=20)

#--------------------------------------------------------------------------------#
# Fit Model 1: Full Model - Sire + Dam as Random
#--------------------------------------------------------------------------------#

# Note on terminology:
#   We often call the 'full model' the model with all terms. Basically we assume
#   everything that may impact the trait could be statistically significant. 
#   We then do 'model selection' to reduce the model. So we call these models
#   'reduced'. 

# formula for full model to test first
# Fitting
#   Fixed: Sex and % Duroc
#   Random: Sire and Dam
mod_full = @formula(wt_birth_h ~ sex + perc_duroc + (1|sire) + (1|dam))

# (1|sire) means we are just fitting an intercept (mean) for each sire
# later we can fit random regression models where we could fit a linear trend
# per sire or something (1+Age|sire) would be a random regression later... 

# fit mixed model
mme = fit(MixedModel,     # tell it to fit a Mixed Model
              mod_full,   # this is the formula from above (you don't need it separate)
              dts)        # name of the dataset    

# Adjusted means:
# Effects in R or 'LSMeans' from SAS
emmeans(mme)

# y-hat values
predict(mme)

#--------------------------------------------------------------------------------#
# Fit Model 2: Reduced - only dam as random
#--------------------------------------------------------------------------------#

# formula for reduced model
mod_red_1 = @formula(wt_birth_h ~ sex + perc_duroc + (1|dam))

# fit mixed model
mmer = fit(MixedModel, mod_red_1, dts)

# LRT on full vs reduced model
lrtest(mmer, mme)

# decided they are not significantly different (p > 0.05)
# so we can likely drop the sire effect from the model as the dam is absorbind
# the variance

#--------------------------------------------------------------------------------#
# Fit Model 3: Reduced - only dam as random + sire fixed
#--------------------------------------------------------------------------------#

# formula for 1st full model
mod_red_2 = @formula(wt_birth_h ~ sex + perc_duroc + sire + (1|dam))

# fit mixed model
mmer2 = fit(MixedModel, mod_red_2, dts)

# LRT on full vs reduced model
lrtest(mmer2, mme) # does not work on non-nested models like in R

# calculate AIC
aic(mmer2)
aic(mme)

# check BIC
bic(mmer2)
bic(mme)

# mme has a lower BIC so therefore a better model adjusting for complexity
# BIC adds a penalty that scales with the number of rows unlike AIC. 
# Often I suggest BIC for any large models as AIC will tend to select very 
# complex models that don't crossvalidate as well. 

#--------------------------------------------------------------------------------#
# Final Model: Dam as random
#--------------------------------------------------------------------------------#

# fill in missing with 0.5 for % Duroc
dts.perc_duroc = coalesce.(dts.perc_duroc, 0.5)

# formula for full model
mod_final = @formula(wt_birth_h ~ sex + perc_duroc + (1|dam))

# fit mixed model
mmer = fit(MixedModel, mod_final, dts)

#--------------------------------------------------------------------------------#
# MME
#--------------------------------------------------------------------------------#

# APPROXIMATION here:
varu = 1.6
vare = 8
lambda = vare / varu
lambda  

# fill y
y = dts.wt_birth_h

# extract X matrix
X = mmer.X

# formula for dam as fixed (to get Z matrix)
formula_dam_fixed = @formula(wt_birth_h ~ -1 + dam)

# fit mixed model
mm_dam_fixed = fit(LinearModel, formula_dam_fixed, dts)

# Z matrix (fit dam as fixed)
Z = mm_dam_fixed.mm.m

# set up MME components
XpX = X'X
XpZ = X'Z
ZpX = Z'X
ZpZ = Z'Z
XpY = X'y
ZpY = Z'y

# A matrix here will just be diagonal (no relationships and no inbreeding)
A   = I(size(Z, 2))

# build Left-hand side
LHS = [XpX XpZ
       ZpX ZpZ .+ inv(A)*lambda]

# build Right-hand side
RHS = [XpY
       ZpY]

# solve MME (inverse LHS * RHS)
solutions = inv(LHS) * RHS

# fixed effects (here first 3 rows)
fixed_effects = solutions[1:size(X, 2)]

# random effects (dams)
random_effects = solutions[size(X,2)+1:end]

#--------------------------------------------------------------------------------#
# PEVs
#--------------------------------------------------------------------------------#

# inverse LHS
iLHS = inv(LHS)

# set size of X
n = size(X, 2)

# extract diagonals of iLHS
EBV_PEV = diag(iLHS[n+1:end, n+1:end] .* vare)
EBV_SEP = EBV_PEV .^0.5

# hist of PEVs
#histogram(EBV_PEV)

# number of progeny
nprog = freqtable(dts.dam)

# PEV
DFpev = DataFrame(nprog = nprog, PEV = EBV_PEV, SEP = EBV_SEP)

# plot PEV on number of prog
Plots.scatter(DFpev.nprog, DFpev.PEV)
Plots.title!("PEV on Number of Progeny")
Plots.xlabel!("Number of Progeny")
Plots.ylabel!("PEV")

# extract PEV num 5
EBV_PEV[5]

# plot PEV on number of prog
Plots.scatter(DFpev.SEP, DFpev.PEV)
Plots.title!("PEV on SEP")
Plots.xlabel!("SEP")
Plots.ylabel!("PEV")

# calc reliability based on PEV
reliability = 1 .- (EBV_PEV ./ varu)

# calc acc
accuracy = sqrt.(reliability)

# add both to DataFrame
DFpev.Rel = reliability
DFpev.Acc = accuracy

# plot
Plots.scatter(DFpev.Acc, DFpev.Rel)
Plots.title!("Accuracy on Reliability")
Plots.xlabel!("Accuracy")
Plots.ylabel!("Reliability")




