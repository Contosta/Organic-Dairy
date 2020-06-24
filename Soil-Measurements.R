##############################################################################################
##############################################################################################
##############################################################################################

#This script evaluates whether management intensive grazing (MIG) changes soil carbon and 
#nitrogen stocks in three organic dairy farms located across New England.

#Code developed by A. Contosta.
#Most recent version 6/19/20

##############################################################################################
##############################################################################################
##############################################################################################
#Initial Set Up

#call libraries
library(data.table)
library(Hmisc)
library(splitstackshape)
library(zoo)
library(nlme)
library(matrixStats)
library(stringr)
library(multcomp)
library(plyr)
library(multcomp)
library(swfscMisc)
library(lattice)

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
cncyc = read.table("cncyc.csv", head = T, sep = ",")

##############################################################################################
##############################################################################################
##############################################################################################

#Data preprocessing

##############
#zero values #
##############

#correct BD and perC, and perN for zero values

mig.bd$BD = ifelse(mig.bd$BD == 0, NA, mig.bd$BD)

dat$perC = ifelse(dat$perC == 0, NA, dat$perC)
dat$perN = ifelse(dat$perN == 0, NA, dat$perN)


##################################################################################
#caluclate C, N, and labile C stocks from bulk density and C and N concentrations#
##################################################################################

#make ID fields for merging bulk density and soil data and for specifying random effects

mig.bd$id = paste(mig.bd$Farm, mig.bd$MF, mig.bd$Plot, mig.bd$Depth.Fraction)
dat$id = paste(dat$Farm, dat$Field, dat$Plot, dat$Elevation)
cncyc$id = paste(cncyc$Farm, cncyc$Field, cncyc$Plot)

all.dat = merge(mig.bd, dat, by.x = "id", by.y = "id", all.x = T, all.y = T)

#select columns of interest for calculating stocks
all.d = all.dat[ , c("id", "Farm.x", "Year", "Mgmt", "Field.x", "MF", "Plot.x", "Subsample", "Depth.Fraction", "Core_Length",
        "BD.x", "BD.y", "perN", "perC", "POXc", "cstock", "nstock")]

names(all.d) = c("id", "Farm", "Year", "Mgmt", "Field", "MF", "Plot", "Subsample", "Depth.Fraction", "Core_Length",
        "BD.1", "BD.2", "perN", "perC", "POXc", "cstock.2", "nstock.2")

#calculate C and N stocks based on actual BD, not averaged across cores (which is how it was 
#done previously). Units are kg / ha

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
                    POXc = sum(POXc, na.rm = T), cstock = sum(cstock.2, na.rm = T), 
                    nstock = sum(nstock.2, na.rm = T)), by = id2]

dat.5$CN = dat.5$cstock / dat.5$nstock

#merge cncyc and soils data from 0 - 15 cm to look at patterns in total stocks and 
#potential turnover

#subset to dat.4 to only include 0-15 cm
dat.4sub = dat.4[dat.4$Elevation == "15", ]

soilcn = merge(cncyc, dat.4sub, by.x = "id", by.y = "id2")

##############################################################################################
##############################################################################################
##############################################################################################

#Data characterization (normality, homogeneity of variance, correlation)

#####################
#Normal distribution#
#####################

#soil enzymes
par(mfrow = c(2,2))

hist(cncyc$CBH)
hist(cncyc$BG)
hist(cncyc$LAP)
hist(cncyc$NAG)

hist(log(cncyc$CBH))
hist(log(cncyc$BG))
hist(log(cncyc$LAP))
hist(log(cncyc$NAG))
#data more normal with log transformation

#net C turnover
par(mfrow = c(1,2))
hist(cncyc$C.flux)
hist(cncyc$Net.C.Enz)

hist(log(cncyc$C.flux))
hist(log(cncyc$Net.C.Enz))
#Net.C.Enz more normal with log transformation

par(mfrow = c(1,3))
hist(cncyc$N.flux)
hist(cncyc$Net.NO3)
hist(cncyc$Net.N.Enz)

hist(log(cncyc$N.flux))
hist(log(cncyc$Net.NO3))
hist(log(cncyc$Net.N.Enz))
#N.Flux and Net.N.Enz more normal with log transformation

#soil BD and agreggates
par(mfrow = c(2,3))

hist(dat.3$BD)
hist(dat.3$ag2000)
hist(dat.3$ag250)
hist(dat.3$ag0)
hist(dat.3$ag0)
hist(dat.3$mwd)

hist(log(dat.3$BD))
hist(log(dat.3$ag2000))
hist(log(dat.3$ag250))
hist(log(dat.3$ag53))
hist(log(dat.3$ag0))
hist(log(dat.3$mwd))
#data more normal with log transformation

#soil C and N
par(mfrow = c(2,2))

hist(dat.3$labc.1)
hist(dat.3$cstock.1)
hist(dat.3$nstock.1)
hist(dat.3$CN)

hist(log(dat.3$labc.1))
hist(log(dat.3$cstock.1))
hist(log(dat.3$nstock.1))
hist(log(dat.3$CN))
#log transformation may not be needed

hist(dat.5$POXc)
hist(dat.5$cstock)
hist(dat.5$nstock)
hist(dat.5$CN)

hist(log(dat.5$POXc))
hist(log(dat.5$cstock))
hist(log(dat.5$nstock))
hist(log(dat.5$CN))
#log transformation may not be needed

################
#Equal variance#
################

#soil enzymes
par(mfrow = c(2,2))

boxplot(log(BG) ~ Treatment * Farm, data = cncyc)
boxplot(log(CBH) ~ Treatment * Farm, data = cncyc)
boxplot(log(NAG) ~ Treatment * Farm, data = cncyc)
boxplot(log(LAP) ~ Treatment * Farm, data = cncyc)
#all enzymes show unequal variace among groups

#soil C cycling
par(mfrow = c(1,2))
boxplot(C.flux ~ Treatment * Farm, data = cncyc)
boxplot(Net.C.Enz ~ Treatment * Farm, data = cncyc)

#soil N cycling
par(mfrow = c(1,3))
boxplot(log(N.flux) ~ Treatment * Farm, data = cncyc)
boxplot(Net.N.Enz ~ Treatment * Farm, data = cncyc)
boxplot(Net.NO3 ~ Treatment * Farm, data = cncyc)

#soil BD and agreggates
par(mfrow = c(2,3))

boxplot(BD ~ Treatment * Elevation * Farm, data = dat.3)
boxplot(ag2000 ~ Treatment * Farm, subset = dat.3$Elevation == 15, data = dat.3)
boxplot(ag250 ~ Treatment * Farm, subset = dat.3$Elevation == 15, data = dat.3)
boxplot(ag53 ~ Treatment * Farm, subset = dat.3$Elevation == 15, data = dat.3)
boxplot(ag0 ~ Treatment * Farm, subset = dat.3$Elevation == 15, data = dat.3)
boxplot(mwd ~ Treatment * Farm, subset = dat.3$Elevation == 15, data = dat.3)
#some differences in variation among groups

#soil C and N
par(mfrow = c(2,2))

boxplot(labc.1 ~ Treatment * Elevation * Farm, data = dat.3)
boxplot(cstock.1 ~ Treatment * Elevation * Farm, data = dat.3)
boxplot(nstock.1 ~ Treatment * Elevation * Farm, data = dat.3)
boxplot(CN ~ Treatment * Elevation * Farm, data = dat.3)

boxplot(POXc ~ Treatment * Farm, data = dat.5)
boxplot(cstock ~ Treatment * Farm, data = dat.5)
boxplot(nstock ~ Treatment * Farm, data = dat.5)
boxplot(CN ~ Treatment * Farm, data = dat.5)
#some differences in variation among groups

################
#Correlations  #
################

#plot correlations between variables

MyNames <- c("C.flux", "N.flux", "Net.N.Enz", "Net.NO3", "Net.C.Enz",
             "perC", "perN", "CN", "POXc", "mwd", 
             "cstock.1", "nstock.1")

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs(soilcn[,c("C.flux", "N.flux", "Net.N.Enz", "Net.NO3", "Net.C.Enz",
              "perC", "perN", "CN", "POXc", "mwd", 
              "cstock.1", "nstock.1")],
      upper.panel = panel.smooth,
      lower.panel = panel.cor,
      cex.labels=1.3,
      labels=MyNames)

cor(soilcn[,c("C.flux", "N.flux", "Net.N.Enz", "Net.NO3", "Net.C.Enz",
              "perC", "perN", "CN", "POXc", "mwd", 
              "cstock.1", "nstock.1")])


##############################################################################################
##############################################################################################
##############################################################################################

#Statistical analysis

#########
#Enzymes#
#########

#############
#C degrading#
#############

#select random effects
gls.ccncyc = gls(Net.C.Enz ~ Farm + Treatment + Farm:Treatment, data = cncyc, na.action = na.omit)
lme.ccncyc = lme(fixed = Net.C.Enz ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, 
              random = ~1|id, data = cncyc)
anova(gls.ccncyc, lme.ccncyc)
#random effects do not improve overall model fit

#select variance structure
var.ccncyc.1 = update(gls.ccncyc, weights = varIdent(form = ~1|Farm))
var.ccncyc.2 = update(gls.ccncyc, weights = varIdent(form = ~1|Treatment))
var.ccncyc.3 = update(gls.ccncyc, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.ccncyc, var.ccncyc.1, var.ccncyc.2, var.ccncyc.3)
#variance structure does not improve model fit

#examine fit of fixed effects
anova(gls.ccncyc, type = "marginal")
#all three are significant

#refit with REML to select fixed effects
ccncyc.1 = update(gls.ccncyc, method = "ML")

#remove interaction between Farm * Treatment
ccncyc.2 = update(ccncyc.1, Net.C.Enz ~ Farm + Treatment)

#compare model fits
anova(ccncyc.1, ccncyc.2)
#model fit significantly worse. gls.ccnyc is final model

#refit with REML
ccncyc.fin = update(ccncyc.1, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(ccncyc.fin)

#make qqplot (verify normality)
qqnorm(ccncyc.fin)

#make new column for Farm * Treatment for pairwise comparisons
cncyc$farm_trt = as.factor(paste(cncyc$Farm, cncyc$Treatment, sep = " "))

#create models for all single factor and interaction terms
MCccncyc.farm = gls(Net.C.Enz ~ Farm, na.action = na.omit, data = cncyc)
MCccncyc.trt = gls(Net.C.Enz ~ Treatment, na.action = na.omit, data = cncyc)
MCccncyc.farm_trt = gls(Net.C.Enz ~ farm_trt, na.action = na.omit, data = cncyc)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(MCccncyc.farm, linfct = mcp(Farm = "Tukey"))))
summary((glht(MCccncyc.trt, linfct = mcp(Treatment = "Tukey"))))
summary((glht(MCccncyc.farm_trt, linfct = mcp(farm_trt = "Tukey"))))

#############
#N degrading#
#############

#select random effects
gls.ncncyc = gls(log(Net.N.Enz) ~ Farm + Treatment + Farm:Treatment, data = cncyc, na.action = na.omit)
lme.ncncyc = lme(fixed = Net.N.Enz ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, 
                 random = ~1|id, data = cncyc)
anova(gls.ncncyc, lme.ncncyc)
#random effects do not improve overall model fit

#select variance structure
var.ncncyc.1 = update(gls.ncncyc, weights = varIdent(form = ~1|Farm))
var.ncncyc.2 = update(gls.ncncyc, weights = varIdent(form = ~1|Treatment))
var.ncncyc.3 = update(gls.ncncyc, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.ncncyc, var.ncncyc.1, var.ncncyc.2, var.ncncyc.3)
#variance structure of Farm*Treatment improves model fit

#examine fit of fixed effects
anova(var.ncncyc.3, type = "marginal")
#only Farm is significant

#refit with REML to select fixed effects
ncncyc.1 = update(var.ncncyc.3, method = "ML")

#remove interaction between Farm * Treatment
ncncyc.2 = update(ncncyc.1, Net.N.Enz ~ Farm + Treatment)

#compare model fits
anova(ncncyc.1, ncncyc.2)
#no difference in model fit

#examine fit of remaining single effects
anova(ncncyc.2, type = "marginal")
#both Farm and Treatment are significant. ncncyc.2 is now full model

#remove each single fixed effect in turn
#remove the effect of Treatment
ncncyc.3 = update(ncncyc.2, Net.N.Enz ~ Farm)

#remove the effect of Farm
ncncyc.4 = update(ncncyc.2, Net.N.Enz ~ Treatment)

#compare model fits
anova(ncncyc.2, ncncyc.3)
anova(ncncyc.2, ncncyc.4)
#removal of Farm significantly worsens model fit.
#ncncyc.3 final model

#refit with REML
ncncyc.fin = update(ncncyc.3, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(ncncyc.fin)

#make qqplot (verify normality)
qqnorm(ncncyc.fin)

#create models for all single factor and interaction terms
Mnccncyc.farm = gls(Net.N.Enz ~ Farm, na.action = na.omit, data = cncyc)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(Mnccncyc.farm, linfct = mcp(Farm = "Tukey"))))

##########################
#Potential Mineralization#
##########################

#######
#C Min#
#######

#select random effects
gls.Cmin = gls(C.flux ~ Farm + Treatment + Farm:Treatment, data = cncyc, na.action = na.omit)
lme.Cmin = lme(fixed = C.flux ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, 
                 random = ~1|id, data = cncyc)
anova(gls.Cmin, lme.Cmin)
#random effects do not improve overall model fit

#select variance structure
var.Cmin.1 = update(gls.Cmin, weights = varIdent(form = ~1|Farm))
var.Cmin.2 = update(gls.Cmin, weights = varIdent(form = ~1|Treatment))
var.Cmin.3 = update(gls.Cmin, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.Cmin, var.Cmin.1, var.Cmin.2, var.Cmin.3)
#variance structure of Farm improves model fit

#examine fit of fixed effects
anova(var.Cmin.1, type = "marginal")
#only Farm is significant

#refit with REML to select fixed effects
Cmin.1 = update(var.Cmin.1, method = "ML")

#remove interaction between Farm * Treatment
Cmin.2 = update(Cmin.1, C.flux ~ Farm + Treatment)

#compare model fits
anova(Cmin.1, Cmin.2)
#no difference in model fit

#examine fit of remaining single effects
anova(Cmin.2, type = "marginal")
#both Farm and Treatment are significant. Cmin.2 is now full model

#remove each single fixed effect in turn
#remove the effect of Treatment
Cmin.3 = update(Cmin.2, C.flux ~ Farm)

#remove the effect of Farm
Cmin.4 = update(Cmin.2, C.flux ~ Treatment)

#compare model fits
anova(Cmin.2, Cmin.3)
anova(Cmin.2, Cmin.4)
#removal of Farm significantly worsens model fit.
#Cmin.3 final model

#refit with REML
Cmin.fin = update(Cmin.3, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(Cmin.fin)

#make qqplot (verify normality)
qqnorm(Cmin.fin)

#create models for all single factor and interaction terms
Mnccncyc.farm = gls(C.flux ~ Farm, na.action = na.omit, weights = varIdent(form = ~1|Farm), data = cncyc)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(Mnccncyc.farm, linfct = mcp(Farm = "Tukey"))))

#######
#N Min#
#######

#select random effects
gls.Nmin = gls(N.flux ~ Farm + Treatment + Farm:Treatment, data = cncyc, na.action = na.omit)
lme.Nmin = lme(fixed = N.flux ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, 
               random = ~1|id, data = cncyc)
anova(gls.Nmin, lme.Nmin)
#random effects do not improve overall model fit

#select variance structure
var.Nmin.1 = update(gls.Nmin, weights = varIdent(form = ~1|Farm))
var.Nmin.2 = update(gls.Nmin, weights = varIdent(form = ~1|Treatment))
var.Nmin.3 = update(gls.Nmin, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.Nmin, var.Nmin.1, var.Nmin.2, var.Nmin.3)
#variance structure of Farm*Treatment improves model fit

#examine fit of fixed effects
anova(var.Nmin.3, type = "marginal")
#only Farm is significant

#refit with REML to select fixed effects
Nmin.1 = update(var.Nmin.3, method = "ML")

#remove interaction between Farm * Treatment
Nmin.2 = update(Nmin.1, N.flux ~ Farm + Treatment)

#compare model fits
anova(Nmin.1, Nmin.2)
#no difference in model fit

#examine fit of remaining single effects
anova(Nmin.2, type = "marginal")
#both Farm and Treatment are significant. Nmin.2 is now full model

#remove each single fixed effect in turn
#remove the effect of Treatment
Nmin.3 = update(Nmin.2, N.flux ~ Farm)

#remove the effect of Farm
Nmin.4 = update(Nmin.2, N.flux ~ Treatment)

#compare model fits
anova(Nmin.2, Nmin.3)
anova(Nmin.2, Nmin.4)
#removal of Farm significantly worsens model fit.
#Nmin.3 final model

#refit with REML
Nmin.fin = update(Nmin.3, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(Nmin.fin)

#make qqplot (verify normality)
qqnorm(Nmin.fin)

#create models for all single factor and interaction terms
Mnccncyc.farm = gls(N.flux ~ Farm, na.action = na.omit, weights = varIdent(form = ~1|Farm*Treatment), data = cncyc)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(Mnccncyc.farm, linfct = mcp(Farm = "Tukey"))))

########
#Nitrif#
########

#select random effects
gls.nitr = gls(Net.NO3 ~ Farm + Treatment + Farm:Treatment, data = cncyc, na.action = na.omit)
lme.nitr = lme(fixed = Net.NO3 ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, 
               random = ~1|id, data = cncyc)
anova(gls.nitr, lme.nitr)
#random effects do not improve overall model fit

#select variance structure
var.nitr.1 = update(gls.nitr, weights = varIdent(form = ~1|Farm))
var.nitr.2 = update(gls.nitr, weights = varIdent(form = ~1|Treatment))
var.nitr.3 = update(gls.nitr, weights = varIdent(form = ~1|Farm*Treatment))

BIC(gls.nitr, var.nitr.1, var.nitr.2, var.nitr.3)
#variance structure does not improve model fit

#examine fit of fixed effects
anova(gls.nitr, type = "marginal")
#only Farm is significant

#refit with REML to select fixed effects
nitr.1 = update(gls.nitr, method = "ML")

#remove interaction between Farm * Treatment
nitr.2 = update(nitr.1, Net.NO3 ~ Farm + Treatment)

#compare model fits
anova(nitr.1, nitr.2)
#no difference in model fit

#examine fit of remaining single effects
anova(nitr.2, type = "marginal")
#both Farm and Treatment are significant. nitr.2 is now full model

#remove each single fixed effect in turn
#remove the effect of Treatment
nitr.3 = update(nitr.2, Net.NO3 ~ Farm)

#remove the effect of Farm
nitr.4 = update(nitr.2, Net.NO3 ~ Treatment)

#compare model fits
anova(nitr.2, nitr.3)
anova(nitr.2, nitr.4)
#removal of Farm significantly worsens model fit.
#nitr.3 final model

#refit with REML
nitr.fin = update(nitr.3, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(nitr.fin)

#make qqplot (verify normality)
qqnorm(nitr.fin)

#create models for all single factor and interaction terms
Mnnitr.farm = gls(Net.NO3 ~ Farm, na.action = na.omit, data = cncyc)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(Mnccncyc.farm, linfct = mcp(Farm = "Tukey"))))

############
#Aggregates#
############

########
#ag2000#
########
gls.ag2000 = gls(ag2000 ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, subset = dat.3$Elevation == 15, data = dat.3)
lme.ag2000 = lme(fixed = ag2000 ~ Farm + Treatment + Farm:Treatment, random = ~1|id, na.action = na.omit, subset = dat.3$Elevation == 15, data = dat.3)
anova(gls.ag2000, lme.ag2000)
#random effect does not improve overall model fit

#select variance structure
var.ag2000.1 = update(gls.ag2000, weights = varIdent(form = ~ 1 | Farm))
var.ag2000.2 = update(gls.ag2000, weights = varIdent(form = ~ 1 | Treatment))
var.ag2000.3 = update(gls.ag2000, weights = varIdent(form = ~ 1 | Farm * Treatment))

BIC(gls.ag2000, var.ag2000.1, var.ag2000.2, var.ag2000.3)
#model fit not improved with a constant variance structure

#examine fit of fixed effects
anova(gls.ag2000, type = "marginal")
#all effects are significat. gls.ag2000 is final model

#plot model residuals (verify homogeneity of variance)
plot(gls.ag2000)

#make qqplot (verify normality)
qqnorm(gls.ag2000)

#make new column for Farm * Treatment for pairwise comparisons
dat.3$farm_trt = as.factor(paste(dat.3$Farm, dat.3$Treatment, sep = " "))

#create models for all single factor and interaction terms
MCag2000.farm = gls(ag2000 ~ Farm, na.action = na.omit, subset = dat.3$Elevation == 15, data = dat.3)
MCag2000.trt = gls(ag2000 ~ Treatment, na.action = na.omit, subset = dat.3$Elevation == 15, data = dat.3)
MCag2000.farm_trt = gls(ag2000 ~ farm_trt, na.action = na.omit, subset = dat.3$Elevation == 15, data = dat.3)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(MCag2000.farm, linfct = mcp(Farm = "Tukey"))))
summary((glht(MCag2000.trt, linfct = mcp(Treatment = "Tukey"))))
summary((glht(MCag2000.farm_trt, linfct = mcp(farm_trt = "Tukey"))))

########
#ag250#
########
gls.ag250 = gls(ag250 ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, subset = dat.3$Elevation == 15, data = dat.3)
lme.ag250 = lme(fixed = ag250 ~ Farm + Treatment + Farm:Treatment, random = ~1|id, na.action = na.omit, subset = dat.3$Elevation == 15, data = dat.3)
anova(gls.ag250, lme.ag250)
#random effect does not improve overall model fit

#select variance structure
var.ag250.1 = update(gls.ag250, weights = varIdent(form = ~ 1 | Farm))
var.ag250.2 = update(gls.ag250, weights = varIdent(form = ~ 1 | Treatment))
var.ag250.3 = update(gls.ag250, weights = varIdent(form = ~ 1 | Farm * Treatment))

BIC(gls.ag250, var.ag250.1, var.ag250.2, var.ag250.3)
#model fit not improved with a constant variance structure

#examine fit of fixed effects
anova(gls.ag250, type = "marginal")
#all effects are significat. gls.ag250 is final model

#plot model residuals (verify homogeneity of variance)
plot(gls.ag250)

#make qqplot (verify normality)
qqnorm(gls.ag250)

#create models for all single factor and interaction terms
MCag250.farm = gls(ag250 ~ Farm, na.action = na.omit, subset = dat.3$Elevation == 15, data = dat.3)
MCag250.trt = gls(ag250 ~ Treatment, na.action = na.omit, subset = dat.3$Elevation == 15, data = dat.3)
MCag250.farm_trt = gls(ag250 ~ farm_trt, na.action = na.omit, subset = dat.3$Elevation == 15, data = dat.3)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(MCag250.farm, linfct = mcp(Farm = "Tukey"))))
summary((glht(MCag250.trt, linfct = mcp(Treatment = "Tukey"))))
summary((glht(MCag250.farm_trt, linfct = mcp(farm_trt = "Tukey"))))

########
#ag53#
########
gls.ag53 = gls(ag53 ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, subset = dat.3$Elevation == 15, data = dat.3)
lme.ag53 = lme(fixed = ag53 ~ Farm + Treatment + Farm:Treatment, random = ~1|id, na.action = na.omit, subset = dat.3$Elevation == 15, data = dat.3)
anova(gls.ag53, lme.ag53)
#random effect does not improve overall model fit

#select variance structure
var.ag53.1 = update(gls.ag53, weights = varIdent(form = ~ 1 | Farm))
var.ag53.2 = update(gls.ag53, weights = varIdent(form = ~ 1 | Treatment))
var.ag53.3 = update(gls.ag53, weights = varIdent(form = ~ 1 | Farm * Treatment))

BIC(gls.ag53, var.ag53.1, var.ag53.2, var.ag53.3)
#model fit improved with Farm as a constant variance structure

#examine fit of fixed effects
anova(var.ag53.1, type = "marginal")
#only Farm is significant

#refit with REML to select fixed effects
ag53.1 = update(var.ag53.1, method = "ML")

#remove interaction between Farm * Treatment
ag53.2 = update(ag53.1, ag53 ~ Farm + Treatment)

#compare model fits
anova(ag53.1, ag53.2)
#no difference in model fit

#examine fit of remaining single effects
anova(ag53.2, type = "marginal")
#Both Farm and Treatment are significant. ag53.2 is final model

#refit with REML
ag53.fin = update(ag53.2, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(ag53.fin)

#make qqplot (verify normality)
qqnorm(ag53.fin)

#create models for all single factor and interaction terms
MCag53.farm = gls(ag53 ~ Farm, na.action = na.omit, weights = varIdent(form = ~ 1 | Farm), subset = dat.3$Elevation == 15, data = dat.3)
MCag53.trt = gls(ag53 ~ Treatment, na.action = na.omit, weights = varIdent(form = ~ 1 | Farm), subset = dat.3$Elevation == 15, data = dat.3)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(MCag53.farm, linfct = mcp(Farm = "Tukey"))))
summary((glht(MCag53.trt, linfct = mcp(Treatment = "Tukey"))))

########
#ag0#
########
gls.ag0 = gls(ag0 ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, subset = dat.3$Elevation == 15, data = dat.3)
lme.ag0 = lme(fixed = ag0 ~ Farm + Treatment + Farm:Treatment, random = ~1|id, na.action = na.omit, subset = dat.3$Elevation == 15, data = dat.3)
anova(gls.ag0, lme.ag0)
#random effect does not improve overall model fit

#select variance structure
var.ag0.1 = update(gls.ag0, weights = varIdent(form = ~ 1 | Farm))
var.ag0.2 = update(gls.ag0, weights = varIdent(form = ~ 1 | Treatment))
var.ag0.3 = update(gls.ag0, weights = varIdent(form = ~ 1 | Farm * Treatment))

BIC(gls.ag0, var.ag0.1, var.ag0.2, var.ag0.3)
#model fit improved with Farm as a constant variance structure

#examine fit of fixed effects
anova(var.ag0.1, type = "marginal")
#only Farm is significant

#refit with REML to select fixed effects
ag0.1 = update(var.ag0.1, method = "ML")

#remove interaction between Farm * Treatment
ag0.2 = update(ag0.1, ag0 ~ Farm + Treatment)

#compare model fits
anova(ag0.1, ag0.2)
#no difference in model fit

#examine fit of remaining single effects
anova(ag0.2, type = "marginal")
#Only Farm is significant

#remove each single fixed effect in turn
#remove the effect of Treatment
ag0.3 = update(ag0.2, ag0 ~ Farm)

#remove the effect of Treatment
ag0.4 = update(ag0.2, ag0 ~ Treatment)

#compare model fits
anova(ag0.2, ag0.3)
anova(ag0.2, ag0.4)
#removal of Farm worsens model fit
#ag0.3 final model

#refit with REML
ag0.fin = update(ag0.3, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(ag0.fin)

#make qqplot (verify normality)
qqnorm(ag0.fin)

#create models for all single factor and interaction terms
MCag0.farm = gls(ag0 ~ Farm, na.action = na.omit, subset = dat.3$Elevation == 15, data = dat.3)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(MCag0.farm, linfct = mcp(Farm = "Tukey"))))

########
#mwd#
########
gls.mwd = gls(mwd ~ Farm + Treatment + Farm:Treatment, na.action = na.omit, subset = dat.3$Elevation == 15, data = dat.3)
lme.mwd = lme(fixed = mwd ~ Farm + Treatment + Farm:Treatment, random = ~1|id, na.action = na.omit, subset = dat.3$Elevation == 15, data = dat.3)
anova(gls.mwd, lme.mwd)
#random effect does not improve overall model fit

#select variance structure
var.mwd.1 = update(gls.mwd, weights = varIdent(form = ~ 1 | Farm))
var.mwd.2 = update(gls.mwd, weights = varIdent(form = ~ 1 | Treatment))
var.mwd.3 = update(gls.mwd, weights = varIdent(form = ~ 1 | Farm * Treatment))

BIC(gls.mwd, var.mwd.1, var.mwd.2, var.mwd.3)
#model fit not improved with constant variance structure

#examine fit of fixed effects
anova(gls.mwd, type = "marginal")
#only Farm is significant

#refit with REML to select fixed effects
mwd.1 = update(gls.mwd, method = "ML")

#remove interaction between Farm * Treatment
mwd.2 = update(mwd.1, mwd ~ Farm + Treatment)

#compare model fits
anova(mwd.1, mwd.2)
#model fit significantly worse. gls.mwd final model

#examine fit of remaining single effects
anova(mwd.2, type = "marginal")
#Only Farm is significant

#remove each single fixed effect in turn
#remove the effect of Treatment
mwd.3 = update(mwd.2, mwd ~ Farm)

#remove the effect of Treatment
mwd.4 = update(mwd.2, mwd ~ Treatment)

#compare model fits
anova(mwd.2, mwd.3)
anova(mwd.2, mwd.4)
#removal of Farm worsens model fit
#mwd.3 final model

#refit with REML
mwd.fin = update(mwd.3, method = "REML")

#plot model residuals (verify homogeneity of variance)
plot(mwd.fin)

#make qqplot (verify normality)
qqnorm(mwd.fin)

#create models for all single factor and interaction terms
MCmwd.farm = gls(mwd ~ Farm, na.action = na.omit, subset = dat.3$Elevation == 15, data = dat.3)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(MCmwd.farm, linfct = mcp(Farm = "Tukey"))))


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

##############################################################################################
##############################################################################################
##############################################################################################

#Plot Results

##########
#Figure 1#
##########

#calculate average C, N, and CN values over entire profile to add to boxplots
#make new ID column that groups management within farms

dat.5$id3 = paste(dat.5$Farm, dat.5$Treatment, sep = " ")

avgCN = dat.5[ , list(totc = mean(cstock), labc = mean(POXc),
                      totn = mean(nstock), cn = mean(CN)), by = id3]


#setwd
setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\USDA_ORG\\R Projects\\All-Farm-Tradeoffs\\Organic-Dairy\\Figures")

#call pdf.options to define graphical parameters in pdf export file
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="Figure_1_revised_v3.pdf")

par(mfrow = c (2,2), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(5,5,5,5))

boxplot(cstock ~ Treatment * Farm, data = dat.5, axes = F, col = c("coral", "gray80"), ylim = c(50000,240000))
box(lty = 1)
#axis(1, at = c(1:6), lab = c("Graze", "Hay", "Graze", "Hay", "Graze", "Hay"))
points(avgCN$totc, pch = 23, cex = 1.5, bg = "black")
axis(2, at = c(5e4, 10e4, 15e4, 20e4), lab = c(50, 100, 150, 200))
abline(v = 2.5, lty = 2, lwd = 0.5)
abline(v = 4.5, lty = 2, lwd = 0.5)
text(0.65, 23.7e4, "(A)")

braces(xfrom = 0.65, xto = 2.35, yfrom = 21e4, yto = 23e4)
braces(xfrom = 2.65, xto = 4.35, yfrom = 17.5e4, yto = 19.4e4)
braces(xfrom = 4.65, xto = 6.35, yfrom = 17.5e4, yto = 19.4e4)

text(1.5, 23.5e4, "a")
text(3.5, 19.9e4, "b")
text(5.5, 19.9e4, "ab")

mtext(side = 2, expression(paste("Total C (t "*ha^{-1}*")")), line = 3, outer = F)

boxplot(POXc ~ Treatment * Farm, data = dat.5, axes = F, col = c("coral", "gray80"), ylim = c(500,3000))
box(lty = 1)
points(avgCN$labc, pch = 23, cex = 1.5, bg = "black")
axis(4, at = c(500, 1000, 1500, 2000, 2500, 3000), lab = c("0.5", "1.0", "1.5", "2.0", "2.5", "3.0"))
abline(v = 2.5, lty = 2, lwd = 0.5)
abline(v = 4.5, lty = 2, lwd = 0.5)
text(0.65, 2950, "(B)")

braces(xfrom = 0.65, xto = 2.35, yfrom = 2600, yto = 2800)
braces(xfrom = 2.65, xto = 4.35, yfrom = 2000, yto = 2200)
braces(xfrom = 4.65, xto = 6.35, yfrom = 2000, yto = 2200)

text(1.5, 2900, "a")
text(3.5, 2300, "b")
text(5.5, 2300, "c")

mtext(side = 4, expression(paste("Active C (t "*ha^{-1}*")")), line = 3, outer = F)
legend(0.5, 1000, fill = c("coral", "gray80"), c("Grazed", "Hayed"), bty = "n",  cex = 0.85)

boxplot(nstock ~ Treatment * Farm, data = dat.5, axes = F, col = c("coral", "gray80"), ylim = c(3000,19000))
box(lty = 1)
points(avgCN$totn, pch = 23, cex = 1.5, bg = "black")
axis(1, at = c(1.5, 3.5, 5.5), lab = c("VT", "NH", "ME"))
axis(2, at = c(3e3, 6e3, 9e3, 12e3, 15e3, 18e3), lab = c(3, 6, 9, 12, 15, 18))
text(0.65, 18.7e3, "(C)")
abline(v = 2.5, lty = 2, lwd = 0.5)
abline(v = 4.5, lty = 2, lwd = 0.5)

braces(xfrom = 0.65, xto = 2.35, yfrom = 17.e3, yto = 18e3)
braces(xfrom = 2.65, xto = 4.35, yfrom = 14e3, yto = 15e3)
braces(xfrom = 4.65, xto = 6.35, yfrom = 14e3, yto = 15e3)

text(1.5, 18.5e3, "a")
text(3.5, 15.5e3, "b")
text(5.5, 15.5e3, "ab")

text(1, 16.75e3, "x")
text(3, 12e3, "x")
text(5, 14e3, "x")

text(2, 12e3, "y")
text(4, 13e3, "y")
text(6, 10.6e3, "y")

mtext(side = 2, expression(paste("Total N (t "*ha^{-1}*")")), line = 3, outer = F)

boxplot(CN ~ Treatment * Farm, data = dat.5, axes = F, col = c("coral", "gray80"), ylim = c(9, 16))
box(lty = 1)
points(avgCN$cn, pch = 23, cex = 1.5, bg = "black")
axis(1, at = c(1.5, 3.5, 5.5), lab = c("VT", "NH", "ME"))
axis(4, at = c(10, 12, 14, 16))
abline(v = 2.5, lty = 2, lwd = 0.5)
abline(v = 4.5, lty = 2, lwd = 0.5)
text(0.65, 15.8, "(D)")

braces(xfrom = 0.65, xto = 2.35, yfrom = 14, yto = 14.5)
braces(xfrom = 2.65, xto = 4.35, yfrom = 15.25, yto = 15.75)
braces(xfrom = 4.65, xto = 6.35, yfrom = 14, yto = 14.5)

text(1.5, 14.75, "a")
text(3.5, 16, "a")
text(5.5, 14.75, "b")

text(1, 13, "x")
text(3, 13, "x")
text(5, 13, "x")

text(2, 13.5, "y")
text(3.75, 14, "y")
text(6, 13, "y")

mtext(side = 4, "C:N Ratio", line = 3, outer = F)
mtext(side = 1, "Farm", outer = T, line = 2.5)

dev.off()

##########
#Figure 2#
##########

#calculate average values to add to boxplots
cncyc.1 = data.table(cncyc)

avgturn = cncyc.1[ , list(cenz = mean(Net.C.Enz), nenz = mean(Net.N.Enz),
                      C.flux = mean(C.flux), netn = mean(Net.NO3)), by = farm_trt]


#setwd
setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\USDA_ORG\\R Projects\\All-Farm-Tradeoffs\\Organic-Dairy\\Figures")

#call pdf.options to define graphical parameters in pdf export file
pdf.options(width= 6.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="CN_Turnover.pdf")

par(mfrow = c (2,2), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(5,5,5,5))

boxplot(Net.C.Enz ~ Treatment * Farm, data = cncyc, axes = F, col = c("coral", "gray80"), ylim = c(-1.8, 4.2))
box(lty = 1)
points(avgturn$cenz, pch = 23, cex = 1.5, bg = "black")
axis(2, at = c(-1.8, 0, 1.8, 3.6), lab = c("-1.8", "0.0", "1.8", "3.6"))
abline(v = 2.5, lty = 2, lwd = 0.5)
abline(v = 4.5, lty = 2, lwd = 0.5)
text(0.65, 4.175, "(A)")

braces(xfrom = 0.65, xto = 2.35, yfrom = 2.6, yto = 3.0)
braces(xfrom = 2.65, xto = 4.35, yfrom = -0.6, yto = -0.2)
braces(xfrom = 4.65, xto = 6.35, yfrom = 3.6, yto = 4)

text(1.5, 3.2, "a")
text(3.5, 0, "b")
text(5.5, 4.2, "c")

text(1, 2.1, "x")
text(3, -0.7, "x")
text(5, 3.525, "x")

text(2, 2.4, "y")
text(4, -1.0, "x")
text(6, 3.425, "x")

mtext(side = 2, "Net C-Degrading Enzyme Activity", line = 3.5, outer = F)
mtext(side = 2, expression(paste("(pmol "*~g^{-1}*~min^{-1}*")")), line = 2, outer = F)

#N-Degrading Enzymes

boxplot(Net.N.Enz ~ Treatment * Farm, data = cncyc, axes = F, col = c("coral", "gray80"), ylim = c(-1.8, 4.2))
box(lty = 1)
points(avgturn$nenz, pch = 23, cex = 1.5, bg = "black")
axis(4, at = c(-1.8, 0, 1.8, 3.6), lab = c("-1.8", "0.0", "1.8", "3.6"))
abline(v = 2.5, lty = 2, lwd = 0.5)
abline(v = 4.5, lty = 2, lwd = 0.5)
text(0.65, 4.175, "(B)")

braces(xfrom = 0.65, xto = 2.35, yfrom = 3.7, yto = 4.1)
braces(xfrom = 2.65, xto = 4.35, yfrom = -0.3, yto = 0.1)
braces(xfrom = 4.65, xto = 6.35, yfrom = 2.0, yto = 2.4)

text(1.5, 4.3, "a")
text(3.5, 0.3, "b")
text(5.5, 2.6, "c")

mtext(side = 4, "Net N-Degrading Enzyme Activity", line = 2, outer = F)
mtext(side = 4, expression(paste("(pmol "*~g^{-1}*~min^{-1}*")")), line = 3.5, outer = F)
legend(0.5, -1, fill = c("coral", "gray80"), c("Grazed", "Hayed"), bty = "n",  cex = 0.85)

#C Mineralization

boxplot(C.flux ~ Treatment * Farm, data = cncyc, axes = F, col = c("coral", "gray80"), ylim = c(1e4, 5e4))
box(lty = 1)
points(avgturn$C.flux, pch = 23, cex = 1.5, bg = "black")
axis(1, at = c(1.5, 3.5, 5.5), lab = c("VT", "NH", "ME"))
axis(2, at = c(1e4, 2e4, 3e4, 4e4,5e4), lab = c(10, 20, 30, 40, 50))
abline(v = 2.5, lty = 2, lwd = 0.5)
abline(v = 4.5, lty = 2, lwd = 0.5)
text(0.65, 5e4, "(C)")
abline(v = 2.5, lty = 2, lwd = 0.5)
abline(v = 4.5, lty = 2, lwd = 0.5)

braces(xfrom = 0.65, xto = 2.35, yfrom = 2.9e4, yto = 3.1e4)
braces(xfrom = 2.65, xto = 4.35, yfrom = 2e4, yto = 2.2e4)
braces(xfrom = 4.65, xto = 6.35, yfrom = 4e4, yto = 4.2e4)

text(1.5, 3.2e4, "a")
text(3.5, 2.3e4, "b")
text(5.5, 4.3e4, "c")

mtext(side = 2, "Cumulative Potential C Mineralization", line = 3.5, outer = F)
mtext(side = 2, expression(paste("(g "~~CO[2]*~g~soil^{-1}*")")), line = 2, outer = F)

#Nitrification

boxplot(Net.NO3 ~ Treatment * Farm, data = cncyc, axes = F, col = c("coral", "gray80"), ylim = c(-10, 15))
box(lty = 1)
points(avgturn$netn, pch = 23, cex = 1.5, bg = "black")
axis(1, at = c(1.5, 3.5, 5.5), lab = c("VT", "NH", "ME"))
axis(4, at = c(-10, -5, 0, 5, 10, 15))
abline(v = 2.5, lty = 2, lwd = 0.5)
abline(v = 4.5, lty = 2, lwd = 0.5)
text(0.65, 15, "(D)")

braces(xfrom = 0.65, xto = 2.35, yfrom = 2.5, yto = 3.75)
braces(xfrom = 2.65, xto = 4.35, yfrom = 9.5, yto = 10.75)
braces(xfrom = 4.65, xto = 6.35, yfrom = 12.5, yto = 13.75)

text(1.5, 4.75, "a")
text(3.5, 11.75, "b")
text(5.5, 14.75, "b")

mtext(side = 4, "Net Nitrification", line = 2, outer = F)
mtext(side = 4, expression(paste("("*mu*g~NO[3]*~g~soil^{-1}*~day^{-1}*")")), line = 3.5, outer = F)
mtext(side = 1, "Farm", outer = T, line = 2.5)

dev.off()


