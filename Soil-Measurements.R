##############################################################################################
##############################################################################################
##############################################################################################

#This script evaluates whether management intensive grazing (MIG) changes soil carbon and 
#nitrogen stocks in three organic dairy farms located across New England.

#Code developed by A. Contosta.
#Most recent version 9/13/2019

##############################################################################################
##############################################################################################
##############################################################################################
#Initial Set Up

#call libraries
library(data.table)
library(splitstackshape)
library(zoo)
library(matrixStats)
library(stringr)
#library(dunn.test)
library(multcomp)
library(plyr)
library(multcomp)

#run script to introduce multiple comparison methods for 'gls' objects.
model.matrix.gls <- function(object, ...) {
  model.matrix(terms(object), data = getData(object), ...)
}
model.frame.gls <- function(object, ...) {
  model.frame(formula(object), data = getData(object), ...)
}
terms.gls <- function(object, ...) {
  terms(model.frame(object), ...)
}

#set directories and import files
setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\USDA_ORG\\R Projects\\All-Farm-Tradeoffs\\Organic-Dairy\\Data")
dat <- read.table("soil_dat.csv", head = TRUE, sep = ",")
mig.bd = read.table("MIG_bd.csv", head = T, sep = ",")

##############################################################################################
##############################################################################################
##############################################################################################

#data preprocessing

##############
#zero values #
##############

#correct BD and perC, and perN for zero values

mig.bd$BD = ifelse(mig.bd$BD == 0, NA, mig.bd$BD)

dat$perC = ifelse(dat$perC == 0, NA, dat$perC)
dat$perN = ifelse(dat$perN == 0, NA, dat$perN)


##################################################################################
#caluclate C, N, and labile C stocks from bulk density and C  and Nconcentrations#
##################################################################################

#make ID fields for merging bulk density and soil data

mig.bd$id = paste(mig.bd$Farm, mig.bd$MF, mig.bd$Plot, mig.bd$Depth.Fraction)
dat$id = paste(dat$Farm, dat$Field, dat$Plot, dat$Elevation)

all.dat = merge(mig.bd, dat, by.x = "id", by.y = "id", all.x = T, all.y = T)

#select columns of interest for calculating stocks
all.d = all.dat[ , c("id", "Farm.x", "Year", "Mgmt", "Field.x", "MF", "Plot.x", "Subsample", "Depth.Fraction", "Core_Length",
        "BD.x", "BD.y", "perN", "perC", "POXc", "cstock", "nstock")]

names(all.d) = c("id", "Farm", "Year", "Mgmt", "Field", "MF", "Plot", "Subsample", "Depth.Fraction", "Core_Length",
        "BD.1", "BD.2", "perN", "perC", "POXc", "cstock.2", "nstock.2")

#calculate C and N stocks based on actual BD, not averaged across cores (which is how it was 
#done previously)

all.d$cstock.1 = (all.d$perC / 100) * all.d$BD.1 * all.d$Core_Length * (1 / 1000) * (1e8 / 1)
all.d$nstock.1 = (all.d$perN / 100) * all.d$BD.1 * all.d$Core_Length * (1 / 1000) * (1e8 / 1)
all.d$labc.1 = (all.d$POXc / 1e6) * all.d$BD.1 * all.d$Core_Length * (1 / 1000) * (1e8 / 1)

#calculate average stocks

all.d1 = data.table(all.d)

all.d2 = all.d1[ , list(cstock.1 = mean(cstock.1, na.rm = T), nstock.1 = mean(nstock.1, na.rm = T), labc.1 = mean(labc.1, na.rm = T),
                        cstock.2 = mean(cstock.2, na.rm = T) * 10, nstock.2 = mean(nstock.2, na.rm = T) * 10), by = id]

#merge with dat to have all soils variables in one place

dat.2 = merge(dat, all.d2, by.x = "id", by.y = "id", all.x = T, all.y = T)

#select columns
dat.3 = dat.2[ , c("id", "Farm", "Field", "Treatment", "Field.ID", "Plot", "Elevation", "Core.Segment.Length",
        "perN", "perC", "CN", "pH", "BD", "POXc", "moist", "sand", "ag2000", "ag250", "ag53", "ag0",
        "mwd", "cstock.1", "nstock.1", "labc.1", "cstock.2", "nstock.2")]

#####################################
#Calculate total stocks from 0-50 cm#
#####################################

#add POXc, total C, and total N stocks across the measured profile

dat.4 = data.table(dat.3)
  
dat.4$id2 = paste(dat.4$Farm, dat.4$Field, dat.4$Plot, sep = " ")

dat.5 = dat.4[, list(Farm = unique(Farm), Field = unique(Field), Treatment = unique(Treatment), Plot = unique(Plot),
                    POXc = sum(POXc, na.rm = T), cstock = sum(cstock.1, na.rm = T), 
                    nstock = sum(nstock.1, na.rm = T)), by = id2]

dat.5$CN = dat.5$cstock / dat.5$nstock

##############################################################################################
##############################################################################################
##############################################################################################

#Data characterization (normality, homogeneity of variance, correlation)

#####################
#Normal distribution#
#####################

par(mfrow = c(2,2))

hist(dat.3$labc.1)
hist(dat.3$cstock.1)
hist(dat.3$nstock.1)
hist(dat.3$CN)

hist(dat.5$POXc)
hist(dat.5$cstock)
hist(dat.5$nstock)
hist(dat.5$CN)

################
#Equal variance#
################

boxplot(labc.1 ~ Treatment * Elevation * Farm, data = dat.3)
boxplot(cstock.1 ~ Treatment * Elevation * Farm, data = dat.3)
boxplot(nstock.1 ~ Treatment * Elevation * Farm, data = dat.3)
boxplot(CN ~ Treatment * Elevation * Farm, data = dat.3)

boxplot(POXc ~ Treatment * Farm, data = dat.5)
boxplot(cstock ~ Treatment * Farm, data = dat.5)
boxplot(nstock ~ Treatment * Farm, data = dat.5)
boxplot(CN ~ Treatment * Farm, data = dat.5)

##############################################################################################
##############################################################################################
##############################################################################################

#Statistical analysis

############################
#Total C 0 to 15 cm        #
############################

#select random effects
gls.ctot15 = gls(cstock.1 ~ Farm + Treatment + Farm:Treatment, data = dat.3, na.action = na.omit, subset = dat.3$Elevation == 15)
lme.ctot15 = lme(fixed = cstock.1 ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, 
                 random = ~1|id, subset = dat.3$Elevation == 15, data = dat.3)
anova(gls.ctot15, lme.ctot15)
#random effects do not improve overall model fit

#select variance structure
var.ctot15.1 = update(gls.ctot15, weights = varIdent(form = ~1|Farm))
var.ctot15.2 = update(gls.ctot15, weights = varIdent(form = ~1|Treatment))
var.ctot15.3 = update(gls.ctot15, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.ctot15, var.ctot15.1, var.ctot15.2, var.ctot15.3)
#variance structures do not improve model fit

#examine fit of fixed effects
anova(gls.ctot15, type = "marginal")
#only Farm is significant

#refit with REML to select fixed effects
ctot15.1 = update(gls.ctot15, method = "ML")

#remove interaction between Farm * Treatment
ctot15.2 = update(ctot15.1, cstock.1 ~ Farm + Treatment)

#compare model fits
anova(ctot15.1, ctot15.2)
#no difference in model fit

#examine fit of remaining single effects
anova(ctot15.2, type = "marginal")
#only Farm is significant. ctot.2 is now full model

#remove each single fixed effect in turn
#remove the effect of Treatment
ctot15.3 = update(ctot15.2, cstock.1 ~ Farm)

#remove the effect of Treatment
ctot15.4 = update(ctot15.2, cstock.1 ~ Treatment)

#compare model fits
anova(ctot15.2, ctot15.3)
anova(ctot15.2, ctot15.4)
#removal of Treatment has no effect on model fit; removal of Farm significantly worsens model fit
#ctot15.3 final model

#refit with REML
ctot15.fin = update(ctot15.3, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(ctot15.fin)

#make qqplot (verify normality)
qqnorm(ctot15.fin)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(ctot15.fin, linfct = mcp(Farm = "Tukey"))))

############################
#Total C 15 to 30 cm        #
############################

#select random effects
gls.ctot30 = gls(cstock.1 ~ Farm + Treatment + Farm:Treatment, data = dat.3, na.action = na.omit, subset = dat.3$Elevation == 30)
lme.ctot30 = lme(fixed = cstock.1 ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, 
                 random = ~1|id, subset = dat.3$Elevation == 30, data = dat.3)
anova(gls.ctot30, lme.ctot30)
#random effects do not improve overall model fit

#select variance structure
var.ctot30.1 = update(gls.ctot30, weights = varIdent(form = ~1|Farm))
var.ctot30.2 = update(gls.ctot30, weights = varIdent(form = ~1|Treatment))
var.ctot30.3 = update(gls.ctot30, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.ctot30, var.ctot30.1, var.ctot30.2, var.ctot30.3)
#variance structures do not improve model fit

#examine fit of fixed effects
anova(gls.ctot30, type = "marginal")
#only Farm is significant

#refit with REML to select fixed effects
ctot30.1 = update(gls.ctot30, method = "ML")

#remove interaction between Farm * Treatment
ctot30.2 = update(ctot30.1, cstock.1 ~ Farm + Treatment)

#compare model fits
anova(ctot30.1, ctot30.2)
#no difference in model fit

#examine fit of remaining single effects
anova(ctot30.2, type = "marginal")
#neither effect is significant

#remove each single fixed effect in turn
#remove the effect of Treatment
ctot30.3 = update(ctot30.2, cstock.1 ~ Farm)

#remove the effect of Treatment
ctot30.4 = update(ctot30.2, cstock.1 ~ Treatment)

#compare model fits
anova(ctot30.2, ctot30.3)
anova(ctot30.2, ctot30.4)
#no effects were significant

############################
#Total C 30 to 50 cm        #
############################

#select random effects
gls.ctot50 = gls(cstock.1 ~ Farm + Treatment + Farm:Treatment, data = dat.3, na.action = na.omit, subset = dat.3$Elevation == 50)
lme.ctot50 = lme(fixed = cstock.1 ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, 
                 random = ~1|id, subset = dat.3$Elevation == 50, data = dat.3)
anova(gls.ctot50, lme.ctot50)
#random effects do not improve overall model fit

#select variance structure
var.ctot50.1 = update(gls.ctot50, weights = varIdent(form = ~1|Farm))
var.ctot50.2 = update(gls.ctot50, weights = varIdent(form = ~1|Treatment))
var.ctot50.3 = update(gls.ctot50, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.ctot50, var.ctot50.1, var.ctot50.2, var.ctot50.3)
#Farm as a variance structures do improves model fit

#examine fit of fixed effects
anova(var.ctot50.1, type = "marginal")
#only Farm is significant

#refit with REML to select fixed effects
ctot50.1 = update(var.ctot50.1, method = "ML")

#remove interaction between Farm * Treatment
ctot50.2 = update(ctot50.1, cstock.1 ~ Farm + Treatment)

#compare model fits
anova(ctot50.1, ctot50.2)
#no difference in model fit

#examine fit of remaining single effects
anova(ctot50.2, type = "marginal")
#only Farm is significant. ctot50.2 is now full model

#remove each single fixed effect in turn
#remove the effect of Treatment
ctot50.3 = update(ctot50.2, cstock.1 ~ Farm)

#remove the effect of Treatment
ctot50.4 = update(ctot50.2, cstock.1 ~ Treatment)

#compare model fits
anova(ctot50.2, ctot50.3)
anova(ctot50.2, ctot50.4)
#removal of Treatment has no effect on model fit; removal of Farm significantly worsens model fit
#ctot50.3 final model

#refit with REML
ctot50.fin = update(ctot50.3, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(ctot50.fin)

#make qqplot (verify normality)
qqnorm(ctot50.fin)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(ctot50.fin, linfct = mcp(Farm = "Tukey"))))

#########################
#Total C from 0 to 50 cm#
#########################

#select random effects
gls.ctot = gls(cstock ~ Farm + Treatment + Farm:Treatment, data = dat.5)
lme.ctot = lme(fixed = cstock ~ Farm + Treatment + Farm:Treatment, random = ~1|id2, data = dat.5)
anova(gls.ctot, lme.ctot)
#random effects do not improve overall model fit

#select variance structure
var.ctot.1 = update(gls.ctot, weights = varIdent(form = ~1|Farm))
var.ctot.2 = update(gls.ctot, weights = varIdent(form = ~1|Treatment))
var.ctot.3 = update(gls.ctot, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.ctot, var.ctot.1, var.ctot.2, var.ctot.3)
#variance structures do not improve model fit

#examine fit of fixed effects
anova(gls.ctot, type = "marginal")
#only Farm is significant

#refit with REML to select fixed effects
ctot.1 = update(gls.ctot, method = "ML")

#remove interaction between Farm * Treatment
ctot.2 = update(ctot.1, cstock ~ Farm + Treatment)

#compare model fits
anova(ctot.1, ctot.2)
#no difference in model fit

#examine fit of remaining single effects
anova(ctot.2, type = "marginal")
#only Farm is significant. ctot.2 is now full model

#remove each single fixed effect in turn
#remove the effect of Treatment
ctot.3 = update(ctot.2, cstock ~ Farm)

#remove the effect of Treatment
ctot.4 = update(ctot.2, cstock ~ Treatment)

#compare model fits
anova(ctot.2, ctot.3)
anova(ctot.2, ctot.4)
#removal of Treatment has no effect on model fit; removal of Farm significantly worsens model fit
#ctot.3 final model

#refit with REML
ctot.fin = update(ctot.3, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(ctot.fin)

#make qqplot (verify normality)
qqnorm(ctot.fin)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(ctot.fin, linfct = mcp(Farm = "Tukey"))))

############################
#Labile C 0 to 15 cm      #
############################

#select random effects
gls.labc15 = gls(labc.1 ~ Farm + Treatment + Farm:Treatment, data = dat.3, na.action = na.omit, subset = dat.3$Elevation == 15)
lme.labc15 = lme(fixed = labc.1 ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, 
                 random = ~1|id, subset = dat.3$Elevation == 15, data = dat.3)
anova(gls.labc15, lme.labc15)
#random effects do not improve overall model fit

#select variance structure
var.labc15.1 = update(gls.labc15, weights = varIdent(form = ~1|Farm))
var.labc15.2 = update(gls.labc15, weights = varIdent(form = ~1|Treatment))
var.labc15.3 = update(gls.labc15, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.labc15, var.labc15.1, var.labc15.2, var.labc15.3)
#variance structures do not improve model fit

#examine fit of fixed effects
anova(gls.labc15, type = "marginal")
#only Farm is significant

#refit with REML to select fixed effects
labc15.1 = update(gls.labc15, method = "ML")

#remove interaction between Farm * Treatment
labc15.2 = update(labc15.1, labc.1 ~ Farm + Treatment)

#compare model fits
anova(labc15.1, labc15.2)
#no difference in model fit

#examine fit of remaining single effects
anova(labc15.2, type = "marginal")
#only Farm is significant. labc.2 is now full model

#remove each single fixed effect in turn
#remove the effect of Treatment
labc15.3 = update(labc15.2, labc.1 ~ Farm)

#remove the effect of Treatment
labc15.4 = update(labc15.2, labc.1 ~ Treatment)

#compare model fits
anova(labc15.2, labc15.3)
anova(labc15.2, labc15.4)
#removal of Treatment has no effect on model fit; removal of Farm significantly worsens model fit
#labc15.3 final model

#refit with REML
labc15.fin = update(labc15.3, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(labc15.fin)

#make qqplot (verify normality)
qqnorm(labc15.fin)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(labc15.fin, linfct = mcp(Farm = "Tukey"))))

############################
#Labile C 15 to 30 cm      #
############################

#select random effects
gls.labc30 = gls(labc.1 ~ Farm + Treatment + Farm:Treatment, data = dat.3, na.action = na.omit, subset = dat.3$Elevation == 30)
lme.labc30 = lme(fixed = labc.1 ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, 
                 random = ~1|id, subset = dat.3$Elevation == 30, data = dat.3)
anova(gls.labc30, lme.labc30)
#random effects do not improve overall model fit

#select variance structure
var.labc30.1 = update(gls.labc30, weights = varIdent(form = ~1|Farm))
var.labc30.2 = update(gls.labc30, weights = varIdent(form = ~1|Treatment))
var.labc30.3 = update(gls.labc30, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.labc30, var.labc30.1, var.labc30.2, var.labc30.3)
#variance structure of farm improves model fit

#examine fit of fixed effects
anova(var.labc30.1, type = "marginal")
#both Farm and the Farm:Treatment interaction are significant

#refit with REML to select fixed effects
labc30.1 = update(gls.labc30, method = "ML")

#remove interaction between Farm * Treatment
labc30.2 = update(labc30.1, labc.1 ~ Farm + Treatment)

#compare model fits
anova(labc30.1, labc30.2)
#model fit significant worse with interaction removed. var.labc30.1 is final model

labc30.fin = var.labc30.1

#plot model residuals (verify homogeneity of variance)
plot(labc30.fin)

#make qqplot (verify normality)
qqnorm(labc30.fin)

#make new column that combines Farm and Treatment
dat.3$FarmTrt = as.factor(paste(dat.3$Farm, dat.3$Treatment, sep = " "))

#fit models with Farm, Treatment, Farm + Treatment, and Farm:Treatment as fixed effects
MC.labc.30.Farm = update(labc30.fin, labc.1 ~ Farm)
MC.labc.30.Treatment = update(labc30.fin, labc.1 ~ Treatment)
MC.labc.30.FarmTrt = update(labc30.fin, labc.1 ~ FarmTrt)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(MC.labc.30.Farm, linfct = mcp(Farm = "Tukey"))))
summary((glht(MC.labc.30.Treatment, linfct = mcp(Treatment = "Tukey"))))
summary((glht(MC.labc.30.FarmTrt, linfct = mcp(FarmTrt = "Tukey"))))

############################
#Labile C 30 to 50 cm      #
############################

#select random effects
gls.labc50 = gls(labc.1 ~ Farm + Treatment + Farm:Treatment, data = dat.3, na.action = na.omit, subset = dat.3$Elevation == 50)
lme.labc50 = lme(fixed = labc.1 ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, 
                 random = ~1|id, subset = dat.3$Elevation == 50, data = dat.3)
anova(gls.labc50, lme.labc50)
#random effects do not improve overall model fit

#select variance structure
var.labc50.1 = update(gls.labc50, weights = varIdent(form = ~1|Farm))
var.labc50.2 = update(gls.labc50, weights = varIdent(form = ~1|Treatment))
var.labc50.3 = update(gls.labc50, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.labc50, var.labc50.1, var.labc50.2, var.labc50.3)
#Farm as a constant variance structures improves model fit

#examine fit of fixed effects
anova(var.labc50.1, type = "marginal")
#only Farm is significant

#refit with REML to select fixed effects
labc50.1 = update(var.labc50.1, method = "ML")

#remove interaction between Farm * Treatment
labc50.2 = update(labc50.1, labc.1 ~ Farm + Treatment)

#compare model fits
anova(labc50.1, labc50.2)
#no difference in model fit

#examine fit of remaining single effects
anova(labc50.2, type = "marginal")
#only Farm is significant. labc50.2 is now full model

#remove each single fixed effect in turn
#remove the effect of Treatment
labc50.3 = update(labc50.2, labc.1 ~ Farm)

#remove the effect of Treatment
labc50.4 = update(labc50.2, labc.1 ~ Treatment)

#compare model fits
anova(labc50.2, labc50.3)
anova(labc50.2, labc50.4)
#removal of Treatment has no effect on model fit; removal of Farm significantly worsens model fit
#labc50.3 final model

#refit with REML
labc50.fin = update(labc50.3, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(labc50.fin)

#make qqplot (verify normality)
qqnorm(labc50.fin)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(labc50.fin, linfct = mcp(Farm = "Tukey"))))

#########################
#Total C from 0 to 50 cm#
#########################

#select random effects
gls.POXc = gls(POXc ~ Farm + Treatment + Farm:Treatment, data = dat.5)
lme.POXc = lme(fixed = POXc ~ Farm + Treatment + Farm:Treatment, random = ~1|id2, data = dat.5)
anova(gls.POXc, lme.POXc)
#random effects do not improve overall model fit

#select variance structure
var.POXc.1 = update(gls.POXc, weights = varIdent(form = ~1|Farm))
var.POXc.2 = update(gls.POXc, weights = varIdent(form = ~1|Treatment))
var.POXc.3 = update(gls.POXc, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.POXc, var.POXc.1, var.POXc.2, var.POXc.3)
#variance structures do not improve model fit

#examine fit of fixed effects
anova(gls.POXc, type = "marginal")
#only Farm is significant

#refit with REML to select fixed effects
POXc.1 = update(gls.POXc, method = "ML")

#remove interaction between Farm * Treatment
POXc.2 = update(POXc.1, POXc ~ Farm + Treatment)

#compare model fits
anova(POXc.1, POXc.2)
#no difference in model fit

#examine fit of remaining single effects
anova(POXc.2, type = "marginal")
#only Farm is significant. POXc.2 is now full model

#remove each single fixed effect in turn
#remove the effect of Treatment
POXc.3 = update(POXc.2, POXc ~ Farm)

#remove the effect of Treatment
POXc.4 = update(POXc.2, POXc ~ Treatment)

#compare model fits
anova(POXc.2, POXc.3)
anova(POXc.2, POXc.4)
#removal of Treatment has no effect on model fit; removal of Farm significantly worsens model fit
#POXc.3 final model

#refit with REML
POXc.fin = update(POXc.3, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(POXc.fin)

#make qqplot (verify normality)
qqnorm(POXc.fin)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(POXc.fin, linfct = mcp(Farm = "Tukey"))))

############################
#Total N 0 to 15 cm        #
############################

#select random effects
gls.ntot15 = gls(nstock.1 ~ Farm + Treatment + Farm:Treatment, data = dat.3, na.action = na.omit, subset = dat.3$Elevation == 15)
lme.ntot15 = lme(fixed = nstock.1 ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, 
                 random = ~1|id, subset = dat.3$Elevation == 15, data = dat.3)
anova(gls.ntot15, lme.ntot15)
#random effects do not improve overall model fit

#select variance structure
var.ntot15.1 = update(gls.ntot15, weights = varIdent(form = ~1|Farm))
var.ntot15.2 = update(gls.ntot15, weights = varIdent(form = ~1|Treatment))
var.ntot15.3 = update(gls.ntot15, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.ntot15, var.ntot15.1, var.ntot15.2, var.ntot15.3)
#variance structures do not improve model fit

#examine fit of fixed effects
anova(gls.ntot15, type = "marginal")
#only Farm is significant

#refit with REML to select fixed effects
ntot15.1 = update(gls.ntot15, method = "ML")

#remove interaction between Farm * Treatment
ntot15.2 = update(ntot15.1, nstock.1 ~ Farm + Treatment)

#compare model fits
anova(ntot15.1, ntot15.2)
#no difference in model fit

#examine fit of remaining single effects
anova(ntot15.2, type = "marginal")
#only Farm is significant. ntot.2 is now full model

#remove each single fixed effect in turn
#remove the effect of Treatment
ntot15.3 = update(ntot15.2, nstock.1 ~ Farm)

#remove the effect of Treatment
ntot15.4 = update(ntot15.2, nstock.1 ~ Treatment)

#compare model fits
anova(ntot15.2, ntot15.3)
anova(ntot15.2, ntot15.4)
#removal of either Farm or Treatment significantly worsens model fit
#ntot15.2 final model

#refit with REML
ntot15.fin = update(ntot15.2, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(ntot15.fin)

#make qqplot (verify normality)
qqnorm(ntot15.fin)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(ntot15.fin, linfct = mcp(Farm = "Tukey"))))
summary((glht(ntot15.fin, linfct = mcp(Treatment = "Tukey"))))

############################
#Total N 15 to 30 cm        #
############################

#select random effects
gls.ntot30 = gls(nstock.1 ~ Farm + Treatment + Farm:Treatment, data = dat.3, na.action = na.omit, subset = dat.3$Elevation == 30)
lme.ntot30 = lme(fixed = nstock.1 ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, 
                 random = ~1|id, subset = dat.3$Elevation == 30, data = dat.3)
anova(gls.ntot30, lme.ntot30)
#random effects do not improve overall model fit

#select variance structure
var.ntot30.1 = update(gls.ntot30, weights = varIdent(form = ~1|Farm))
var.ntot30.2 = update(gls.ntot30, weights = varIdent(form = ~1|Treatment))
var.ntot30.3 = update(gls.ntot30, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.ntot30, var.ntot30.1, var.ntot30.2, var.ntot30.3)
#variance structures do not improve model fit

#examine fit of fixed effects
anova(gls.ntot30, type = "marginal")
#only Farm is significant

#refit with REML to select fixed effects
ntot30.1 = update(gls.ntot30, method = "ML")

#remove interaction between Farm * Treatment
ntot30.2 = update(ntot30.1, nstock.1 ~ Farm + Treatment)

#compare model fits
anova(ntot30.1, ntot30.2)
#no difference in model fit

#examine fit of remaining single effects
anova(ntot30.2, type = "marginal")
#neither effect is significant

#remove each single fixed effect in turn
#remove the effect of Treatment
ntot30.3 = update(ntot30.2, nstock.1 ~ Farm)

#remove the effect of Treatment
ntot30.4 = update(ntot30.2, nstock.1 ~ Treatment)

#compare model fits
anova(ntot30.2, ntot30.3)
anova(ntot30.2, ntot30.4)
#no effects were significant

############################
#Total N 30 to 50 cm        #
############################

#select random effects
gls.ntot50 = gls(nstock.1 ~ Farm + Treatment + Farm:Treatment, data = dat.3, na.action = na.omit, subset = dat.3$Elevation == 50)
lme.ntot50 = lme(fixed = nstock.1 ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, 
                 random = ~1|id, subset = dat.3$Elevation == 50, data = dat.3)
anova(gls.ntot50, lme.ntot50)
#random effects do not improve overall model fit

#select variance structure
var.ntot50.1 = update(gls.ntot50, weights = varIdent(form = ~1|Farm))
var.ntot50.2 = update(gls.ntot50, weights = varIdent(form = ~1|Treatment))
var.ntot50.3 = update(gls.ntot50, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.ntot50, var.ntot50.1, var.ntot50.2, var.ntot50.3)
#Treatment as a variance structures do improves model fit

#examine fit of fixed effects
anova(var.ntot50.2, type = "marginal")
#both Farm and Treatment are significant

#refit with REML to select fixed effects
ntot50.1 = update(var.ntot50.1, method = "ML")

#remove interaction between Farm * Treatment
ntot50.2 = update(ntot50.1, nstock.1 ~ Farm + Treatment)

#compare model fits
anova(ntot50.1, ntot50.2)
#no difference in model fit

#examine fit of remaining single effects
anova(ntot50.2, type = "marginal")
#only Farm is significant. ntot50.2 is now full model

#remove each single fixed effect in turn
#remove the effect of Treatment
ntot50.3 = update(ntot50.2, nstock.1 ~ Farm)

#remove the effect of Treatment
ntot50.4 = update(ntot50.2, nstock.1 ~ Treatment)

#compare model fits
anova(ntot50.2, ntot50.3)
anova(ntot50.2, ntot50.4)
#removal of Treatment has no effect on model fit; removal of Farm significantly worsens model fit
#ntot50.3 final model

#refit with REML
ntot50.fin = update(ntot50.3, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(ntot50.fin)

#make qqplot (verify normality)
qqnorm(ntot50.fin)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(ntot50.fin, linfct = mcp(Farm = "Tukey"))))

#########################
#Total N from 0 to 50 cm#
#########################

#select random effects
gls.ntot = gls(nstock ~ Farm + Treatment + Farm:Treatment, data = dat.5)
lme.ntot = lme(fixed = nstock ~ Farm + Treatment + Farm:Treatment, random = ~1|id2, data = dat.5)
anova(gls.ntot, lme.ntot)
#random effects do not improve overall model fit

#select variance structure
var.ntot.1 = update(gls.ntot, weights = varIdent(form = ~1|Farm))
var.ntot.2 = update(gls.ntot, weights = varIdent(form = ~1|Treatment))
var.ntot.3 = update(gls.ntot, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.ntot, var.ntot.1, var.ntot.2, var.ntot.3)
#variance structures do not improve model fit

#examine fit of fixed effects
anova(gls.ntot, type = "marginal")
#only Farm is significant

#refit with REML to select fixed effects
ntot.1 = update(gls.ntot, method = "ML")

#remove interaction between Farm * Treatment
ntot.2 = update(ntot.1, nstock ~ Farm + Treatment)

#compare model fits
anova(ntot.1, ntot.2)
#no difference in model fit

#examine fit of remaining single effects
anova(ntot.2, type = "marginal")
#only Farm is significant. ntot.2 is now full model

#remove each single fixed effect in turn
#remove the effect of Treatment
ntot.3 = update(ntot.2, nstock ~ Farm)

#remove the effect of Treatment
ntot.4 = update(ntot.2, nstock ~ Treatment)

#compare model fits
anova(ntot.2, ntot.3)
anova(ntot.2, ntot.4)
#removal of either Treatment or Farm significantly worsens model fit
#ntot.2 final model

#refit with REML
ntot.fin = update(ntot.2, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(ntot.fin)

#make qqplot (verify normality)
qqnorm(ntot.fin)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(ntot.fin, linfct = mcp(Farm = "Tukey"))))
summary((glht(ntot.fin, linfct = mcp(Treatment = "Tukey"))))

############################
#CN Ratio 0 to 15 cm        #
############################

#select random effects
gls.cn15 = gls(CN ~ Farm + Treatment + Farm:Treatment, data = dat.3, na.action = na.omit, subset = dat.3$Elevation == 15)
lme.cn15 = lme(fixed = CN ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, 
                 random = ~1|id, subset = dat.3$Elevation == 15, data = dat.3)
anova(gls.cn15, lme.cn15)
#random effects do not improve overall model fit

#select variance structure
var.cn15.1 = update(gls.cn15, weights = varIdent(form = ~1|Farm))
var.cn15.2 = update(gls.cn15, weights = varIdent(form = ~1|Treatment))
var.cn15.3 = update(gls.cn15, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.cn15, var.cn15.1, var.cn15.2, var.cn15.3)
#variance structures do not improve model fit

#examine fit of fixed effects
anova(gls.cn15, type = "marginal")
#only Farm is significant

#refit with REML to select fixed effects
cn15.1 = update(gls.cn15, method = "ML")

#remove interaction between Farm * Treatment
cn15.2 = update(cn15.1, CN ~ Farm + Treatment)

#compare model fits
anova(cn15.1, cn15.2)
#no difference in model fit

#examine fit of remaining single effects
anova(cn15.2, type = "marginal")
#both Farm and Treatment are significant. cn.2 is now full model

#remove each single fixed effect in turn
#remove the effect of Treatment
cn15.3 = update(cn15.2, CN ~ Farm)

#remove the effect of Treatment
cn15.4 = update(cn15.2, CN ~ Treatment)

#compare model fits
anova(cn15.2, cn15.3)
anova(cn15.2, cn15.4)
#removal of either Farm or Treatment significantly worsens model fit
#cn15.2 final model

#refit with REML
cn15.fin = update(cn15.2, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(cn15.fin)

#make qqplot (verify normality)
qqnorm(cn15.fin)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(cn15.fin, linfct = mcp(Farm = "Tukey"))))
summary((glht(cn15.fin, linfct = mcp(Treatment = "Tukey"))))

############################
#CN Ratio 15 to 30 cm        #
############################

#select random effects
gls.cn30 = gls(CN ~ Farm + Treatment + Farm:Treatment, data = dat.3, na.action = na.omit, subset = dat.3$Elevation == 30)
lme.cn30 = lme(fixed = CN ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, 
                 random = ~1|id, subset = dat.3$Elevation == 30, data = dat.3)
anova(gls.cn30, lme.cn30)
#random effects do not improve overall model fit

#select variance structure
var.cn30.1 = update(gls.cn30, weights = varIdent(form = ~1|Farm))
var.cn30.2 = update(gls.cn30, weights = varIdent(form = ~1|Treatment))
var.cn30.3 = update(gls.cn30, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.cn30, var.cn30.1, var.cn30.2, var.cn30.3)
#treatment as a variance structure improves model fit

#examine fit of fixed effects
anova(var.cn30.2, type = "marginal")
#only Farm is significant

#refit with REML to select fixed effects
cn30.1 = update(var.cn30.2, method = "ML")

#remove interaction between Farm * Treatment
cn30.2 = update(cn30.1, CN ~ Farm + Treatment)

#compare model fits
anova(cn30.1, cn30.2)
#no difference in model fit

#examine fit of remaining single effects
anova(cn30.2, type = "marginal")
#effect of farm is not significant

#remove each single fixed effect in turn
#remove the effect of Treatment
cn30.3 = update(cn30.2, CN ~ Farm)

#remove the effect of Treatment
cn30.4 = update(cn30.2, CN ~ Treatment)

#compare model fits
anova(cn30.2, cn30.3)
anova(cn30.2, cn30.4)
#only the effect of Treatment is significant. cn30.4 is final model

#refit with REML
cn30.fin = update(cn30.4, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(cn30.fin)

#make qqplot (verify normality)
qqnorm(cn30.fin)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(cn30.fin, linfct = mcp(Treatment = "Tukey"))))

############################
#CN Ratio 30 to 50 cm      #
############################

#select random effects
gls.cn50 = gls(CN ~ Farm + Treatment + Farm:Treatment, data = dat.3, na.action = na.omit, subset = dat.3$Elevation == 50)
lme.cn50 = lme(fixed = CN ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, 
                 random = ~1|id, subset = dat.3$Elevation == 50, data = dat.3)
anova(gls.cn50, lme.cn50)
#random effects do not improve overall model fit

#select variance structure
var.cn50.1 = update(gls.cn50, weights = varIdent(form = ~1|Farm))
var.cn50.2 = update(gls.cn50, weights = varIdent(form = ~1|Treatment))
var.cn50.3 = update(gls.cn50, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.cn50, var.cn50.1, var.cn50.2, var.cn50.3)
#Farm as a variance structures improves model fit

#examine fit of fixed effects
anova(var.cn50.1, type = "marginal")
#both Farm and Treatment are significant

#refit with REML to select fixed effects
cn50.1 = update(var.cn50.1, method = "ML")

#remove interaction between Farm * Treatment
cn50.2 = update(cn50.1, CN ~ Farm + Treatment)

#compare model fits
anova(cn50.1, cn50.2)
#model fit significantly worse. cn50.1 final model

#refit with REML
cn50.fin = update(cn50.1, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(cn50.fin)

#make qqplot (verify normality)
qqnorm(cn50.fin)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
#fit models with Farm, Treatment, Farm + Treatment, and Farm:Treatment as fixed effects
MC.cn.50.Farm = update(cn50.fin, CN ~ Farm)
MC.cn.50.Treatment = update(cn50.fin, CN ~ Treatment)
MC.cn.50.FarmTrt = update(cn50.fin, CN ~ FarmTrt)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(MC.cn.50.Farm, linfct = mcp(Farm = "Tukey"))))
summary((glht(MC.cn.50.Treatment, linfct = mcp(Treatment = "Tukey"))))
summary((glht(MC.cn.50.FarmTrt, linfct = mcp(FarmTrt = "Tukey"))))

#########################
#CN Ratio from 0 to 50 cm#
#########################

#select random effects
gls.cn = gls(CN ~ Farm + Treatment + Farm:Treatment, data = dat.5)
lme.cn = lme(fixed = CN ~ Farm + Treatment + Farm:Treatment, random = ~1|id2, data = dat.5)
anova(gls.cn, lme.cn)
#random effects do not improve overall model fit

#select variance structure
var.cn.1 = update(gls.cn, weights = varIdent(form = ~1|Farm))
var.cn.2 = update(gls.cn, weights = varIdent(form = ~1|Treatment))
var.cn.3 = update(gls.cn, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.cn, var.cn.1, var.cn.2, var.cn.3)
#variance structures do not improve model fit

#examine fit of fixed effects
anova(gls.cn, type = "marginal")
#no effects are significant

#refit with REML to select fixed effects
cn.1 = update(gls.cn, method = "ML")

#remove interaction between Farm * Treatment
cn.2 = update(cn.1, CN ~ Farm + Treatment)

#compare model fits
anova(cn.1, cn.2)
#no difference in model fit

#examine fit of remaining single effects
anova(cn.2, type = "marginal")
#both Farm and Treatment are significant. cn.2 is now full model

#remove each single fixed effect in turn
#remove the effect of Treatment
cn.3 = update(cn.2, CN ~ Farm)

#remove the effect of Treatment
cn.4 = update(cn.2, CN ~ Treatment)

#compare model fits
anova(cn.2, cn.3)
anova(cn.2, cn.4)
#removal of either Treatment or Farm significantly worsens model fit
#cn.2 final model

#refit with REML
cn.fin = update(cn.2, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(cn.fin)

#make qqplot (verify normality)
qqnorm(cn.fin)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(cn.fin, linfct = mcp(Farm = "Tukey"))))
summary((glht(cn.fin, linfct = mcp(Treatment = "Tukey"))))

#############plot differences in C and N 

#setwd
setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\USDA_ORG\\R Projects\\All-Farm-Tradeoffs\\Organic-Dairy\\Figures")

#call pdf.options to define graphical parameters in pdf export file
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="Figure_1_revised.pdf")

par(mfrow = c (2,2), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(5,5,5,5))

boxplot(cstock ~ Treatment * Farm, data = dat.5, axes = F, col = c("coral", "gray80"), ylim = c(50000,240000))
box(lty = 1)
#axis(1, at = c(1:6), lab = c("Graze", "Hay", "Graze", "Hay", "Graze", "Hay"))
axis(2, at = c(5e4, 10e4, 15e4, 20e4), lab = c(50, 100, 150, 200))
abline(v = 2.5, lty = 2, lwd = 0.5)
abline(v = 4.5, lty = 2, lwd = 0.5)
text(0.65, 23.8e4, "(A)")

braces(xfrom = 0.65, xto = 2.35, yfrom = 21e4, yto = 23e4)
braces(xfrom = 2.65, xto = 4.35, yfrom = 17.5e4, yto = 19.4e4)
braces(xfrom = 4.65, xto = 6.35, yfrom = 17.5e4, yto = 19.4e4)

text(1.5, 23.5e4, "a")
text(3.5, 19.9e4, "b")
text(5.5, 19.9e4, "ab")

mtext(side = 2, expression(paste("Total C (t "*ha^{-1}*")")), line = 3, outer = F)

boxplot(POXc ~ Treatment * Farm, data = dat.5, axes = F, col = c("coral", "gray80"), ylim = c(500,3000))
box(lty = 1)
axis(4, at = c(500, 1000, 1500, 2000, 2500, 3000), lab = c("0.5", "1.0", "1.5", "2.0", "2.5", "3.0"))
abline(v = 2.5, lty = 2, lwd = 0.5)
abline(v = 4.5, lty = 2, lwd = 0.5)
text(0.65, 2900, "(B)")

braces(xfrom = 0.65, xto = 2.35, yfrom = 2600, yto = 2800)
braces(xfrom = 2.65, xto = 4.35, yfrom = 2000, yto = 2200)
braces(xfrom = 4.65, xto = 6.35, yfrom = 2000, yto = 2200)

text(1.5, 2900, "a")
text(3.5, 2300, "b")
text(5.5, 2300, "ab")

mtext(side = 4, expression(paste("Active C (t "*ha^{-1}*")")), line = 3, outer = F)
legend(5, 2900, fill = c("coral", "gray80"), c("Grazed", "Hayed"), bty = "n",  cex = 0.85)

boxplot(nstock ~ Treatment * Farm, data = dat.5, axes = F, col = c("coral", "gray80"), ylim = c(3000,18000))
box(lty = 1)
axis(1, at = c(1.5, 3.5, 5.5), lab = c("VT", "NH", "ME"))
axis(2, at = c(3e3, 6e3, 9e3, 12e3, 15e3, 18e3), lab = c(3, 6, 9, 12, 15, 18))
text(0.55, 17.6e3, "(C)")
abline(v = 2.5, lty = 2, lwd = 0.5)
abline(v = 4.5, lty = 2, lwd = 0.5)

braces(xfrom = 0.65, xto = 2.35, yfrom = 16e3, yto = 17e3)
braces(xfrom = 2.65, xto = 4.35, yfrom = 14e3, yto = 15e3)
braces(xfrom = 4.65, xto = 6.35, yfrom = 14e3, yto = 15e3)

text(1.5, 17.5e3, "a")
text(3.5, 15.5e3, "b")
text(5.5, 15.5e3, "ab")

braces(xfrom = 0.65, xto = 5.35, yfrom = 6e3, yto = 5e3)
braces(xfrom = 1.65, xto = 6.35, yfrom = 4.5e3, yto = 3.5e3)

text(3, 4.5e3, "x")
text(4, 3e3, "y")

mtext(side = 2, expression(paste("Total N (t "*ha^{-1}*")")), line = 3, outer = F)

boxplot(CN ~ Treatment * Farm, data = dat.5, axes = F, col = c("coral", "gray80"), ylim = c(9, 16))
box(lty = 1)
axis(1, at = c(1.5, 3.5, 5.5), lab = c("VT", "NH", "ME"))
axis(4, at = c(10, 12, 14, 16))
abline(v = 2.5, lty = 2, lwd = 0.5)
abline(v = 4.5, lty = 2, lwd = 0.5)
text(0.55, 16, "(D)")

braces(xfrom = 0.65, xto = 2.35, yfrom = 14, yto = 14.5)
braces(xfrom = 2.65, xto = 4.35, yfrom = 15.25, yto = 15.75)
braces(xfrom = 4.65, xto = 6.35, yfrom = 14, yto = 14.5)

text(1.5, 14.75, "a")
text(3.5, 16, "b")
text(5.5, 14.75, "ab")

braces(xfrom = 0.65, xto = 5.35, yfrom = 10.5, yto = 10)
braces(xfrom = 1.65, xto = 6.35, yfrom = 9.75, yto = 9.25)

text(3, 9.75, "x")
text(4, 9, "y")

mtext(side = 4, "C:N Ratio", line = 3, outer = F)
mtext(side = 1, "Farm", outer = T, line = 2.5)

dev.off()

