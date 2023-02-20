
#--------------------------------------------------------------------------------#
# Lab 03.R
#--------------------------------------------------------------------------------#

# Let's obtain the OLS of the longley dataset using matrix algebra operations 
# and LM Work with a dataset from Gualdr et al. (2014). 
# (doi: 10.1186/1471-2105-15-246)

#--------------------------------------------------------------------------------#
# Load Libraries
#--------------------------------------------------------------------------------#

# load libraries
library(lme4)
library(lmerTest)
library(tidyverse)

# file
data_file = "pheno_ready.txt"

# platform
my_platform = .Platform$OS.type
print(my_platform)

# get current working directory
getwd()

# check if you can find the current data file in your working directory
if (data_file %in% list.files()){
  # print message
  cat("Found file in your current working directory!")
} else {
  # find file within OS
  if (my_platform %in% c("unix", "linux")) {
    file_location = system(command = paste0("find ~/Documents/ -name \"", data_file, "\" | head -n 1"), intern=TRUE)
    file_location = gsub("pheno_ready.txt", "", x = file_location)
    setwd(file_location)
  } else {
    system(comamnd = paste0("Get-ChildItem -Path C:\ -Filter ", data_file, " -File"))
  }
}

# set working directory
getwd()

#--------------------------------------------------------------------------------#
# Read Data
#--------------------------------------------------------------------------------#

# read dataset
dts <- read.table("pheno_ready.txt")

# head
head(dts)

# size of dataset
dim(dts)

# table of Dam and Sire
table(dts$dam)
table(dts$sire)

#------------------------------------------------------------------------------#
# Manage Data
#------------------------------------------------------------------------------#

# manage data
dts <- dts %>%
  mutate( 
    wt_birth_h = wt_birth * 10, # change to hectograms
    dam        = as.factor(dam),
    litter     = as.factor(litter),
    sire       = as.factor(sire),
    slgdt_cd   = as.factor(slgdt_cd)
  ) 

#------------------------------------------------------------------------------#
# Full Model - Both sire and dam as random
#------------------------------------------------------------------------------#

# First: fit a simple linear model
# lmer is from the lme4 package (lme4::lmer())
mme <- lmer(wt_birth_h ~ sex + perc_duroc + (1|sire) + (1|dam), 
            data=dts)

# print model
mme

# summary of model
summary(mme)

# calculate & print variance components estimates
vr <- VarCorr(mme)
vr
print(vr, comp="Variance")

# notice sire variance! is it zero?

#------------------------------------------------------------------------------#
# Fit Reduced Model 1 - Dam Only as Random
#------------------------------------------------------------------------------#

# fit model with only Dam as random
mmer <- lmer(wt_birth_h ~ sex + perc_duroc + (1|dam),
           data=dts)

# print model
mmer

# calculate & print variance components estimates
vrr <- VarCorr(mmer)
print(vrr, comp="Variance")

# Use a likelihood ratio test (LRT) to check for significance of sire variance
# we will learn about this in the 2nd part of this course
# but for now, all you it's necessary to follow this lab is
# we are testing for significance of that one variance component 

# print test for full model (sire + dam as random) and reduced (dam only)
anova(mme, mmer)

#------------------------------------------------------------------------------#
# Fit Reduced Model 2 - Dam Only as Random + sire fixed
#------------------------------------------------------------------------------#

# model with sire as fixed effect...
mmer2 <- lmer(wt_birth_h ~ sex + perc_duroc + sire + (1|dam),
              data=dts)

# print summary of model
summary(mmer2)

# print test
anova(mmer, mmer2) 

# not the best way to test this!!!

# Does this mean we can discount sire as an important effect?

#------------------------------------------------------------------------------#
# Back to reduced model 1 - Dam Only as Random
#------------------------------------------------------------------------------#

# Let's settle on a model:
mmer <- lmer(wt_birth_h ~ sex + perc_duroc + (1|dam),
             data=dts)

# print model
mmer

# calculate & print variance components estimates
vrr <- VarCorr(mmer)
print(vrr, comp="Variance")

#------------------------------------------------------------------------------#
# Check missing
#------------------------------------------------------------------------------#

# check for missing values BEFORE proceeding -Class-
sum(is.na(dts$wt_birth_h))
sum(is.na(dts$sex))
sum(is.na(dts$sire)) 
sum(is.na(dts$dam)) 
sum(is.na(dts$perc_duroc)) 

# missing! what to do?

# solution: naive imputation for % Duroc (50%)
dts$perc_duroc[is.na(dts$perc_duroc)] <- 0.5

# check how many are missing now
sum(is.na(dts$perc_duroc))

#------------------------------------------------------------------------------#
# Solve MME
#------------------------------------------------------------------------------#

# APPROXIMATION here:
varu <- 1.6
vare <- 8

# ratio we need for MME
lambda <- vare / varu
lambda

# build X Matrix
X <- model.matrix(~ sex + perc_duroc, data=dts)
head(X)

# build Z Matrix
Z <- model.matrix(~dam-1,data=dts)
print(Z)

# response (y) variable
y <- dts$wt_birth_h

# blocks
XtX <- t(X) %*% X
XtZ <- t(X) %*% Z
Zty <- t(Z) %*% y
Xty <- t(X) %*% y
ZtZ <- t(Z) %*% Z
A   <- diag(ncol(Z))

# LHS
LHS <- rbind(
          cbind(  XtX,  XtZ), 
          cbind(t(XtZ), ZtZ + solve(A)*lambda)
        )

# dim
dim(LHS)

# row/column names
colnames(LHS)
rownames(LHS)

# set up RHS
RHS <- rbind(Xty,Zty)
RHS

# solve solutions
slt <- solve(LHS) %*% RHS

# head
head(slt)

#------------------------------------------------------------------------------#
# PEVs
#------------------------------------------------------------------------------#

# fit model
mmer <- lmer(wt_birth_h ~ sex + perc_duroc + (1|dam), 
             data=dts)
mmer

# extract variance components
vrr <- VarCorr(mmer)
print(vrr, comp="Variance")

# extract random effects from model
uhat_lmer <- ranef(mmer)

# random effects solutions we solved for by hand
uhat_mme <- slt[-(1:3)]

# breeding values
compa_u <- cbind(uhat_lmer$dam$`(Intercept)`, uhat_mme)  

# plot
plot(compa_u, 
     pch=19, 
     xlab = "lmer u_hat", 
     ylab = "MME solve u_hat")
abline(0,1, lty="dashed")

# 3) Obtain PEV of u_hat and obtain variance of beta_hat (demo + hands on)

# Hint: BEWARE of MME where error variance was factored out!!!

# do not show this:
dim(LHS)

# inverse LHS
Cm <- solve(LHS)
Cm

# dim
dim(Cm)

# solve for PEV
PEV <- diag(Cm[-(1:3), -(1:3)] * vare)

# hist of PEVs
hist(PEV)

# get SEP (standard error)
SEP <- sqrt(PEV)

# check: represent PEV as function of number of progeny
nprog <- as.numeric(table(dts$dam))
nprog

# plot PEV on number of progeny
plot(PEV ~ nprog)
plot(SEP ~ nprog)

# ??
#identify(nprog, PEV)

# print PEV and nprog
PEV[5]
nprog[5]

# dam 47
PEV[47]
nprog[47]

# 
dts[dts$dam==313,]
dts[dts$dam==491,]

# SD
sd(dts[dts$dam==313, "wt_birth_h"])
sd(dts[dts$dam==491, "wt_birth_h"])

# I can't find a reason for this!

# reliability
reliability <- 1 - (PEV/varu)

# accuracy
accuracy <- sqrt(reliability)

# hist rel
hist(reliability, 
     col="dodgerblue3", 
     main = "Reliability of EBVs")

# hist acc
hist(accuracy, 
     col="dodgerblue3", 
     main = "Accuracy of EBVs")

# plot rel vs acc
plot(x = accuracy, 
     y = reliability, 
     type = "p", 
     col = "orange",
     xlim = c(0, 1), 
     ylim = c(0, 1), 
     main = "Reliability on Accuracy")


#------------------------------------------------------------------------------#
# Homework
#------------------------------------------------------------------------------#

# 4) fit an equivalent model using GLS: show beta_hat, show u_hat

# 5) repeat the prediction of random effects using other variance ratios.
lambda <- c(0.2, 1.0, 5.0, 20)

# compare these estimates to each other

# explain the results
  
# 6) another equivalent model: Z* = I and A*=ZZt sigma2_u 

# rebuild MME and refit model, show equivalence

# which MME is larger? explain dimensions










