##############################################################################################
##############################################################################################
##############################################################################################

#This script evaluates whether management intensive grazing (MIG) changes aboveground biomass
#production in three organic dairy farms located across New England.

#Code developed by A. Contosta.

##############################################################################################
##############################################################################################
##############################################################################################
#Initial Set Up

#call libraries

library(data.table)
library(splitstackshape)
library(zoo)
library(matrixStats)
library(MASS)
library(stringr)
library(dplyr)
library(nlme)
library(multcomp)
library(pBrackets)
library(reshape2)  # to draw pattern in boxplot


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

#import file
dat <- read.table("MIG_GHG_Flux_v2.csv", head = TRUE, sep = ",", na.strings=c("NA", "NAN"))

##############################################################################################
##############################################################################################
##############################################################################################

#Data preprocessing

#check that variable names and classes are correct

str(dat)

#convert date to year and doy
dat$year = as.numeric(strftime(as.POSIXct(strptime(dat$Date.Sampled, "%m/%d/%Y", tz = "EST")), "%Y"))
dat$doy = as.numeric(strftime(as.POSIXct(strptime(dat$Date.Sampled, "%m/%d/%Y", tz = "EST")), "%j"))

##################
#Calculate fluxes#
##################

#open containers for variables to be included post-regression
farm = c()
year = c()
dates = c()
site = c()
reps = c()
subp = c()

t_air = c()
t_5cm = c()
t_10cm = c()
moist = c()

cflmod = c()
cflslope = c()
cflpval = c()
cflr2 = c()

cflux = c()

nflmod = c()
nflslope = c()
nflpval = c()
nflr2 = c()

nflux = c()

#############################################
#run loop for calculating CO2 and N2O fluxes#
#############################################

for (i in seq(1,nrow(dat))) {

    #write information about sampling
    farm[i] = dat$Farm[i]
    year[i] = dat$year[i]
    dates[i] = dat$doy[i]
    site[i] = dat$Site[i]
    reps[i] = dat$Rep[i]
    subp[i] = dat$Plot[i]

    #write environmental data
    t_air[i] = (dat$T_Air[i])
    t_5cm[i] = (dat$T_5cm[i])
    t_10cm[i] = (dat$T_10cm[i])
    moist[i] = (dat$AVG_Moist[i])

    #run CO2 regressions with all possbile combinations of 3, 4, and 5 sampling points
    lmc1 = lm(t(dat[i, c(12:16)]) ~ c(0, 5, 10, 15, 20))
    lmc2 = lm(t(dat[i, c(12:15)]) ~ c(0, 5, 10, 15))
    lmc3 = lm(t(dat[i, c(12, 14:16)]) ~ c(0, 10, 15, 20))
    lmc4 = lm(t(dat[i, c(12, 13, 15, 16)]) ~ c(0, 5, 15, 20))
    lmc5 = lm(t(dat[i, c(12, 13, 14, 16)]) ~ c(0, 5, 10, 20))
    lmc5 = lm(t(dat[i, c(13, 14, 15, 16)]) ~ c(5, 10, 15, 20))
    lmc6 = lm(t(dat[i, c(12:14)]) ~ c(0, 5, 10))
    lmc7 = lm(t(dat[i, c(12, 13, 15)]) ~ c(0, 5, 15))
    lmc8 = lm(t(dat[i, c(12, 13, 16)]) ~ c(0, 5, 20))
    lmc9 = lm(t(dat[i, c(12, 14, 15)]) ~ c(0, 10, 15))
    lmc10 = lm(t(dat[i, c(12, 15, 16)]) ~ c(0, 15, 20))
    lmc11 = lm(t(dat[i, c(13, 14, 15)]) ~ c(5, 10, 15))
    lmc12 = lm(t(dat[i, c(13, 14, 16)]) ~ c(5, 10, 20))
    lmc13 = lm(t(dat[i, c(13, 15, 16)]) ~ c(5, 15, 20))
    lmc14 = lm(t(dat[i, c(14, 15, 16)]) ~ c(10, 15, 20))

    #extract CO2 regression results
    cflmods = c("lmc1", "lmc2", "lmc3", "lmc4", "lmc5", "lmc6", "lmc7", "lmc8", "lmc9", "lmc10", "lmc11", "lmc12", "lmc13", "lmc14")
    cflslopes = c(summary(lmc1)[[4]][2], summary(lmc2)[[4]][2],summary(lmc3)[[4]][2], summary(lmc4)[[4]][2], summary(lmc5)[[4]][2],
              summary(lmc6)[[4]][2], summary(lmc7)[[4]][2], summary(lmc8)[[4]][2], summary(lmc9)[[4]][2], summary(lmc10)[[4]][2],
              summary(lmc11)[[4]][2], summary(lmc12)[[4]][2], summary(lmc13)[[4]][2], summary(lmc14)[[4]][2])
    cflpvals = c(summary(lmc1)[[4]][8], summary(lmc2)[[4]][8],summary(lmc3)[[4]][8], summary(lmc4)[[4]][8],  summary(lmc5)[[4]][8],
             summary(lmc6)[[4]][8], summary(lmc7)[[4]][8], summary(lmc8)[[4]][8], summary(lmc9)[[4]][8], summary(lmc10)[[4]][8],
             summary(lmc11)[[4]][8], summary(lmc12)[[4]][8], summary(lmc13)[[4]][8], summary(lmc14)[[4]][8])

    cflr2s = c(summary(lmc1)[[9]], summary(lmc2)[[9]],summary(lmc3)[[9]], summary(lmc4)[[9]],  summary(lmc5)[[9]],
           summary(lmc6)[[9]], summary(lmc7)[[9]], summary(lmc8)[[9]], summary(lmc9)[[9]], summary(lmc10)[[9]],
           summary(lmc11)[[9]], summary(lmc12)[[9]], summary(lmc13)[[9]], summary(lmc14)[[9]])

    #choose the CO2 model with the highest r2
    cflmod[i] = (cflmods[cflr2s == max(cflr2s, na.rm = TRUE)][1])
    cflslope[i] = (cflslopes[cflr2s == max(cflr2s, na.rm = TRUE)][1])
    cflpval[i] = (cflpvals[cflr2s == max(cflr2s, na.rm = TRUE)][1])
    cflr2[i] = (max(cflr2s, na.rm = TRUE)[1])

    #calculate CO2 flux using the slope of the regression with the highest r2
    #units are in mg CO2-C m2 / h

    cflux[i] =  ((((cflslope[i] * 60) * 12.012 ) / 1000) / 0.0729)

    #run CO2 regressions with all possbile combinations of 3, 4, and 5 sampling points
    lmn1 = lm(t(dat[i, c(17:21)]) ~ c(0, 5, 10, 15, 20))
    lmn2 = lm(t(dat[i, c(17:20)]) ~ c(0, 5, 10, 15))
    lmn3 = lm(t(dat[i, c(17, 19:21)]) ~ c(0, 10, 15, 20))
    lmn4 = lm(t(dat[i, c(17, 18, 19, 20)]) ~ c(0, 5, 15, 20))
    lmn5 = lm(t(dat[i, c(17, 18, 19, 20)]) ~ c(0, 5, 10, 20))
    lmn5 = lm(t(dat[i, c(18, 19, 20, 21)]) ~ c(5, 10, 15, 20))
    lmn6 = lm(t(dat[i, c(17:19)]) ~ c(0, 5, 10))
    lmn7 = lm(t(dat[i, c(17, 18, 19)]) ~ c(0, 5, 15))
    lmn8 = lm(t(dat[i, c(17, 18, 20)]) ~ c(0, 5, 20))
    lmn9 = lm(t(dat[i, c(17, 19, 20)]) ~ c(0, 10, 15))
    lmn10 = lm(t(dat[i, c(17, 20, 21)]) ~ c(0, 15, 20))
    lmn11 = lm(t(dat[i, c(18, 19, 20)]) ~ c(5, 10, 15))
    lmn12 = lm(t(dat[i, c(18, 19, 21)]) ~ c(5, 10, 20))
    lmn13 = lm(t(dat[i, c(18, 20, 21)]) ~ c(5, 15, 20))
    lmn14 = lm(t(dat[i, c(19, 20, 21)]) ~ c(10, 15, 20))

    #extract CO2 regression results
    nflmods = c("lmn1", "lmn2", "lmn3", "lmn4", "lmn5", "lmn6", "lmn7", "lmn8", "lmn9", "lmn10", "lmn11", "lmn12", "lmn13", "lmn14")
    nflslopes = c(summary(lmn1)[[4]][2], summary(lmn2)[[4]][2],summary(lmn3)[[4]][2], summary(lmn4)[[4]][2], summary(lmn5)[[4]][2],
              summary(lmn6)[[4]][2], summary(lmn7)[[4]][2], summary(lmn8)[[4]][2], summary(lmn9)[[4]][2], summary(lmn10)[[4]][2],
              summary(lmn11)[[4]][2], summary(lmn12)[[4]][2], summary(lmn13)[[4]][2], summary(lmn14)[[4]][2])
    nflpvals = c(summary(lmn1)[[4]][8], summary(lmn2)[[4]][8],summary(lmn3)[[4]][8], summary(lmn4)[[4]][8],  summary(lmn5)[[4]][8],
             summary(lmn6)[[4]][8], summary(lmn7)[[4]][8], summary(lmn8)[[4]][8], summary(lmn9)[[4]][8], summary(lmn10)[[4]][8],
             summary(lmn11)[[4]][8], summary(lmn12)[[4]][8], summary(lmn13)[[4]][8], summary(lmn14)[[4]][8])

    nflr2s = c(summary(lmn1)[[9]], summary(lmn2)[[9]],summary(lmn3)[[9]], summary(lmn4)[[9]],  summary(lmn5)[[9]],
           summary(lmn6)[[9]], summary(lmn7)[[9]], summary(lmn8)[[9]], summary(lmn9)[[9]], summary(lmn10)[[9]],
           summary(lmn11)[[9]], summary(lmn12)[[9]], summary(lmn13)[[9]], summary(lmn14)[[9]])

    #choose the N2O model with the highest r2
    nflmod[i] = (nflmods[nflr2s == max(nflr2s, na.rm = TRUE)][1])
    nflslope[i] = (nflslopes[nflr2s == max(nflr2s, na.rm = TRUE)][1])
    nflpval[i] = (nflpvals[nflr2s == max(nflr2s, na.rm = TRUE)][1])
    nflr2[i] = max(nflr2s, na.rm = TRUE)[1]

    #calculate N2O flux using the slope of the regression with the highest r2
    #units are in mg N2O-N m2 / h
    nflux[i] =  nflslope[i] * ((dat$AVG_Depth[i] * (15.42^2) * pi) / 0.0729) * (1 / (273.15 + dat$T_Air[i])) * (1/0.08206) *
                (1 / 1000000) * ((14.0067 * 2) / 1) * (1 / 1000) * (60 / 1)


    #write data table
    flux.calc = cbind.data.frame(farm, year, dates, site, reps, subp, t_air, t_5cm, t_10cm, moist,
                cflmod, cflslope, cflpval, cflr2, cflux, nflmod, nflslope, nflpval, nflr2, nflux)


    }

######################
#QA calculated fluxes#
######################

#remove values that are below the threshold of r2 >= 0.90 / CO2 fluxes that are below 0
flux.calc$finco = ifelse(flux.calc$cflr2 < 0.9 | flux.calc$cflux < 0, NA, flux.calc$cflux)
flux.calc$finno = ifelse(flux.calc$nflr2 < 0.9, NA, flux.calc$nflux)

flux.calc$posno = ifelse(flux.calc$nflr2 < 0.9 | flux.calc$nflux < 0, NA, flux.calc$nflux)

#write flux table
write.table(flux.calc, file = paste("flux_calc.csv"), sep = ",", na="NA", append=FALSE, col.names=TRUE, row.names = FALSE)

#calculate mean and median values for examining trends among farms, mgmt, and years
avg.flux = data.table(flux.calc)
avg.flux$id = paste(avg.flux$farm, avg.flux$year, avg.flux$reps, avg.flux$subp)

avg.flux.1 = avg.flux[ , list(avgc = mean(finco, na.rm = T), medc = median(finco, na.rm = T),
             avgn = mean(posno, na.rm = T), medn = median(posno, na.rm = T)), by = id]

#separate id column
sps <- data.frame(do.call(rbind, str_split(avg.flux.1$id, " ")))
names(sps) <- c("farm", "year", "reps", "subp")

avg.flux.2 = cbind(sps, avg.flux.1)

#recode management system
avg.flux.2$mgmt = ifelse(avg.flux.2$reps == 1 | avg.flux.2$reps == 2, "G", "H")

########################################################
#interpolate values for calculating total seasonal flux#
########################################################

#make matrix with all possible dates for each farm, management system, field, and plot combination 
farm = c(1:3)
year = c(2016, 2017)
dates = c(139:304)
site = c(1:2)
reps = c(1:4)
subp = c(1:3)

all.flux = data.frame(matrix(vector(), nrow = length(farm) * length(year) * length(dates) * length(reps) * length(subp),
           ncol = 6, dimnames = list(c(), c("farm", "year", "dates", "site", "reps", "subp"))))

all.flux$farm = farm
all.flux = all.flux[order(all.flux$farm), ]

all.flux$subp = c(1:3)
all.flux = all.flux[order(all.flux$farm, all.flux$subp), ]

all.flux$year =  rep(year, each = 664)
all.flux = all.flux[order(all.flux$farm, all.flux$year, all.flux$subp), ]

all.flux$dates = rep(dates, each = 4)
all.flux = all.flux[order(all.flux$farm, all.flux$year, all.flux$dates, all.flux$subp), ]

all.flux$reps = reps
all.flux = all.flux[order(all.flux$farm, all.flux$year, all.flux$reps, all.flux$subp), ]

all.flux$site = ifelse(all.flux$reps < 3, 1, 2)

#create unique IDs for merging all.flux and flux.calc
flux.calc$id = paste(flux.calc$farm, flux.calc$year, flux.calc$dates, flux.calc$site, flux.calc$reps, flux.calc$subp)
all.flux$id = paste(all.flux$farm, all.flux$year, all.flux$dates, all.flux$site, all.flux$reps, all.flux$subp)

#merge tables together
all.calc = merge(all.flux, flux.calc, by.x = "id", by.y = "id", all.x = T, all.y = T)

#gap fill missing data with linear interpolation

#make sure data are in plot by date order
all.calc = all.calc[order(all.calc$farm.x, all.calc$reps.x, all.calc$subp.x, all.calc$year.x, all.calc$dates.x), ]

#make ID field for grouping interpolations
all.calc$int = paste(all.calc$farm.x, all.calc$reps.x, all.calc$subp.x, all.calc$year.x)

#use na.approx inside mutate to interpolate within plots

#make all.calc a data table
all.calc.1 = data.table(all.calc)

#use na.approx within mutate to do the interpolation
all.calc.2 = data.table(all.calc %>%
   group_by(int) %>%
   mutate(cflux_int =  na.approx(finco, na.rm=FALSE), nflux_int = na.approx(posno, na.rm = F)))

#scale to 24 h and convert units to get kg ha d-1
all.calc.2$cint_d = (all.calc.2$cflux_int * 24) * (1 / 1e6) * (1e5 / 1)
all.calc.2$nint_d = (all.calc.2$nflux_int * 24) * (1 / 1e6) * (1e5 / 1)

#add values within sites and determine the number of days within each sum to scale total seasonal fluxes
all.calc.3 = all.calc.2[ , list(totc = sum(cint_d, na.rm = T), lenc = length(cint_d[!is.na(cint_d)]),
             totn = sum(nint_d, na.rm = T), lenn = length(nint_d[!is.na(nint_d)])), by = int]


#median number of fluxes for N2O is 128 days. So all fluxes will be scaled to that
all.calc.3$totc_cor = all.calc.3$totc * 128 / all.calc.3$lenc
all.calc.3$totn_cor = all.calc.3$totn * 128 / all.calc.3$lenn

#separate id column
sps <- data.frame(do.call(rbind, str_split(all.calc.3$int, " ")))
names(sps) <- c("farm", "reps", "subp", "year")

#add id column to data frame
all.calc.4 = cbind(sps, all.calc.3)

#add column for management
all.calc.4$mgmt = ifelse(all.calc.4$reps == 1 | all.calc.4$reps == 2, "G", "H")

#convert to t / ha / y
all.calc.4$totc_g = all.calc.4$totc_cor / 1000
all.calc.4$totn_g = all.calc.4$totn_cor / 1000

##############################################################################################
##############################################################################################
##############################################################################################

#statistical analysis

##############################################################################################
##############################################################################################
##############################################################################################

#Data characterization (normality, homogeneity of variance, correlation)

#####################
#Normal distribution#
#####################

par(mfrow = c(2,2))

hist(flux.calc$finco)
hist(flux.calc$posno)

hist(log(flux.calc$finco))
hist(log(flux.calc$posno))

hist(all.calc.4$totc_g)
hist(all.calc.4$totn_g)

hist(log(all.calc.4$totc_g))
hist(log(all.calc.4$totn_g))

################
#Equal variance#
################

boxplot(log(finco) ~ year * farm *site, data = flux.calc)
boxplot(log(posno) ~ year * farm *site, data = flux.calc)

boxplot(log(totc_g) ~ year * farm *mgmt, data = all.calc.4)
boxplot(log(totn_g) ~ year * farm *mgmt, data = all.calc.4)


##############################################################################################
##############################################################################################
##############################################################################################

#Statistical analysis

#make new column in flux.calc for year_doy (for modeling possible autocorrelation structure in data
flux.calc$ydoy = as.integer(paste(flux.calc$year, flux.calc$dates, sep = ""))

#make new ID column for specifying random intercept effects
flux.calc$ID2 = paste(flux.calc$farm, flux.calc$site, flux.calc$reps, flux.calc$subp, sep = " ")
all.calc.4$ID = paste(all.calc.4$farm, all.calc.4$reps, all.calc.4$subp, sep = " ")

#factorize other variables to run ANOVA-type analysis
flux.calc$year = factor(flux.calc$year)
flux.calc$farm = factor(flux.calc$farm)
flux.calc$site = factor(flux.calc$site)
flux.calc$reps = factor(flux.calc$reps)
flux.calc$subp = factor(flux.calc$subp)
flux.calc$resu = as.factor(paste(flux.calc$reps, flux.calc$subp))

##############################
#Repeated Measures CO2 Fluxes#
##############################

#repeated measures
gls.cf = gls(log(finco) ~ year + farm + site + year:site + year:farm + farm:site + year:farm:site, na.action = na.omit, data = flux.calc)
lme.cf = lme(fixed = log(finco) ~ year + farm + site + year:site + year:farm + farm:site + year:farm:site, random = ~1|ID2, na.action = na.omit, data = flux.calc)
anova(gls.cf, lme.cf)

#random effect does not improve overall model fit

var.cf.1 = update(gls.cf, weights = varIdent(form = ~ 1 | year))
var.cf.2 = update(gls.cf, weights = varIdent(form = ~ 1 | farm))
var.cf.3 = update(gls.cf, weights = varIdent(form = ~ 1 | site))
var.cf.4 = update(gls.cf, weights = varIdent(form = ~ 1 | year*farm))
var.cf.5 = update(gls.cf, weights = varIdent(form = ~ 1 | year*site))

BIC(gls.cf, var.cf.1, var.cf.2, var.cf.3, var.cf.4, var.cf.5)
#variance structure of farm improves overall model fit

#examine the effect of an autocorrelation structure
ac.cf = update(gls.cf, correlation = corAR1(form = ~ 1|ydoy))
BIC(gls.cf, ac.cf)
#autocorrelation improves model fit

#combine variance and autocorrelation structure
var.ac.cf = update(gls.cf, weights = varIdent(form = ~ 1 | farm), correlation = corAR1(form = ~ 1|ydoy))
BIC(gls.cf, var.cf.2, var.ac.cf)
#including both variance and autocorrelation structure improves overall model fit

#examine fit of fixed effects
anova(var.ac.cf, type = "marginal")

#three way interaction is not significant.
#refit with REML to select fixed effects
cf.1 = update(var.ac.cf, method = "ML")

#refit without three-way interaction
cf.2 = update(cf.1, log(finco) ~ year + farm + site + year:site + year:farm + farm:site)
anova(cf.1, cf.2)
#no difference in model fit

#remove two-way interaction to test significance

#remove the year:site interaction
cf.3 = update(cf.1, log(finco) ~ year + farm + site + year:farm + farm:site)

#remove the year:farm interaction
cf.4 = update(cf.1, log(finco) ~ year + farm + site + year:site + farm:site)

#remove the farm:site interaction
cf.5 = update(cf.1, log(finco) ~ year + farm + site + year:site + year:farm)

anova(cf.2, cf.3)
anova(cf.2, cf.4)
anova(cf.2, cf.5)

#significant differences when year:site interaction is removed
#cf.6 now full model
cf.6 = update(cf.1, log(finco) ~ year + farm + site + year:farm + farm:site)
#all effects significant.

#refit with REML for obtaining summary statistics
cf.fin = update(cf.6, method = 'REML')


#refit with REML for obtaining summary statistics
cf.fin = var.ac.cf

anova(cf.fin, type = "marginal")
summary(cf.fin)

plot(cf.fin)
qqnorm(cf.fin)

#pairwise comparisons

#create new columns to create interaction terms
flux.calc$farm_site = as.factor(paste(flux.calc$farm, flux.calc$site, sep = " "))
flux.calc$farm_year = as.factor(paste(flux.calc$farm, flux.calc$year, sep = " "))
flux.calc$year_site = as.factor(paste(flux.calc$year, flux.calc$site, sep = " "))
flux.calc$farm_site_year = as.factor(paste(flux.calc$farm, flux.calc$site, flux.calc$year, sep = " "))

#create models for all single factor and interaction terms
MCtscarb.farm <- update(cf.fin, log(finco) ~ farm)
MCtscarb.year <- update(cf.fin, log(finco) ~ year)
MCtscarb.site <- update(cf.fin, log(finco) ~ site)

MCtscarb.farm_year <- update(cf.fin, log(finco) ~ farm_year)
MCtscarb.farm_site <- update(cf.fin, log(finco) ~ farm_site)
MCtscarb.year_site <- update(cf.fin, log(finco) ~ year_site)

MCtscarb.farm_site_year <- update(cf.fin, log(finco) ~ farm_site_year)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(MCtscarb.farm, linfct = mcp(farm = "Tukey"))))
summary((glht(MCtscarb.year, linfct = mcp(year = "Tukey"))))
summary((glht(MCtscarb.site, linfct = mcp(site = "Tukey"))))

summary((glht(MCtscarb.farm_year, linfct = mcp(farm_year = "Tukey"))))
summary((glht(MCtscarb.farm_site, linfct = mcp(farm_site = "Tukey"))))
summary((glht(MCtscarb.year_site, linfct = mcp(year_site = "Tukey"))))

summary((glht(MCtscarb.farm_site_year, linfct = mcp(farm_site_year = "Tukey"))))

###########################
#Total Seasonal CO2 Fluxes#
###########################

gls.scf = gls(log(totc_g) ~ year + farm + mgmt + year:mgmt + year:farm + farm:mgmt + year:farm:mgmt, na.action = na.omit, data = all.calc.4)
lme.scf = lme(fixed = log(totc_g) ~ year + farm + mgmt + year:mgmt + year:farm + farm:mgmt + year:farm:mgmt, random = ~1|ID, na.action = na.omit, data = all.calc.4)
anova(gls.scf, lme.scf)
#random effect improves overall model fit

#select variance structure
var.scf.1 = update(lme.scf, weights = varIdent(form = ~ 1 | year))
var.scf.2 = update(lme.scf, weights = varIdent(form = ~ 1 | farm))
var.scf.3 = update(lme.scf, weights = varIdent(form = ~ 1 | mgmt))
var.scf.4 = update(lme.scf, weights = varIdent(form = ~ 1 | year*farm))
var.scf.5 = update(lme.scf, weights = varIdent(form = ~ 1 | year*mgmt))
var.scf.6 = update(lme.scf, weights = varIdent(form = ~ 1 | year*farm*mgmt))

BIC(lme.scf, var.scf.1, var.scf.2, var.scf.3, var.scf.4, var.scf.5, var.scf.6)
#variance structure of farm improves overall model fit

#examine fit of fixed effects
anova(var.scf.2, type = "marginal")
#three way interaction is not significant. 

#refit with ML for selecting fixed effects
scf.1 = update(var.scf.2, method = "ML")

#remove three-way interaction
scf.2 = update(scf.1, log(totc_g) ~ year + farm + mgmt + year:mgmt + year:farm + farm:mgmt)
anova(scf.1, scf.2)
#no difference in model fit


#remove two-way interaction to test significance

#remove the year:mgmt interaction
scf.3 = update(scf.1, log(totc_g) ~ year + farm + mgmt + year:farm + farm:mgmt)

#remove the year:farm interaction
scf.4 = update(scf.1, log(totc_g) ~ year + farm + mgmt + year:mgmt + farm:mgmt)

#remove the farm:mgmt interaction
scf.5 = update(scf.1, log(totc_g) ~ year + farm + mgmt + year:mgmt + year:farm)

anova(scf.2, scf.3)
anova(scf.2, scf.4)
anova(scf.2, scf.5)

#difference when year:mgmt is removed, but there are no differences when year:farm and farm:site removed
#all effects now significant. 
scf.6 = update(scf.1, log(totc_g) ~ year + farm + mgmt + year:farm + farm:mgmt)

#refit with REML for obtaining summary statistics
scf.fin = update(scf.6, method = 'REML')

anova(scf.fin, type = "marginal")
summary(scf.fin)

plot(scf.fin)
qqnorm(scf.fin)

#pairwise comparisons

#create new columns to create interaction terms
all.calc.4$farm_mgmt = as.factor(paste(all.calc.4$farm, all.calc.4$mgmt, sep = " "))
all.calc.4$farm_year = as.factor(paste(all.calc.4$farm, all.calc.4$year, sep = " "))
all.calc.4$year_mgmt = as.factor(paste(all.calc.4$year, all.calc.4$mgmt, sep = " "))
all.calc.4$farm_mgmt_year = as.factor(paste(all.calc.4$farm, all.calc.4$mgmt, all.calc.4$year, sep = " "))

#create models for all single factor and interaction terms
MCtsseascarb.farm <- update(scf.fin, log(totc_g) ~ farm)
MCtsseascarb.year <- update(scf.fin, log(totc_g) ~ year)
MCtsseascarb.mgmt <- update(scf.fin, log(totc_g) ~ mgmt)

MCtsseascarb.farm_year <- update(scf.fin, log(totc_g) ~ farm_year)
MCtsseascarb.farm_mgmt <- update(scf.fin, log(totc_g) ~ farm_mgmt)
MCtsseascarb.year_mgmt <- update(scf.fin, log(totc_g) ~ year_mgmt)

MCtsseascarb.farm_mgmt_year <- update(scf.fin, log(totc_g) ~ farm_mgmt_year)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(MCtsseascarb.farm, linfct = mcp(farm = "Tukey"))))
summary((glht(MCtsseascarb.year, linfct = mcp(year = "Tukey"))))
summary((glht(MCtsseascarb.mgmt, linfct = mcp(mgmt = "Tukey"))))

summary((glht(MCtsseascarb.farm_year, linfct = mcp(farm_year = "Tukey"))))
summary((glht(MCtsseascarb.farm_mgmt, linfct = mcp(farm_mgmt = "Tukey"))))
summary((glht(MCtsseascarb.year_mgmt, linfct = mcp(year_mgmt = "Tukey"))))

summary((glht(MCtsseascarb.farm_mgmt_year, linfct = mcp(farm_mgmt_year = "Tukey"))))

##############################
#Repeated measures N2O fluxes#
##############################

gls.nf = gls(log(posno) ~ year + farm + site + year:site + year:farm + farm:site + year:farm:site, na.action = na.omit, data = flux.calc)
lme.nf = lme(fixed = log(posno) ~ year + farm + site + year:site + year:farm + farm:site + year:farm:site, random = ~1|ID2, na.action = na.omit, data = flux.calc)
anova(gls.nf, lme.nf)
#random effect does not improve overall model fit

var.nf.1 = update(gls.nf, weights = varIdent(form = ~ 1 | year))
var.nf.2 = update(gls.nf, weights = varIdent(form = ~ 1 | farm))
var.nf.3 = update(gls.nf, weights = varIdent(form = ~ 1 | site))
var.nf.4 = update(gls.nf, weights = varIdent(form = ~ 1 | year*farm))
var.nf.5 = update(gls.nf, weights = varIdent(form = ~ 1 | year*site))
var.nf.6 = update(gls.nf, weights = varIdent(form = ~ 1 | year*farm*site))

BIC(gls.nf, var.nf.1, var.nf.2, var.nf.3, var.nf.4, var.nf.5, var.nf.6)
#variance structure of farm improves overall model fit

#examine the effect of an autocorrelation structure
ac.nf = update(gls.nf, correlation = corAR1(form = ~ 1|ydoy))
BIC(gls.nf, ac.nf)
#autocorrelation improves model fit

#combine variance and autocorrelation structures
var.ac.nf = update(gls.nf, weights = varIdent(form = ~ 1 | farm), correlation = corAR1(form = ~ 1|ydoy))
BIC(gls.nf, var.nf.2, ac.nf, var.ac.nf)
#model fit improved with both variance and autocorrelation structures

#examine fit of fixed effects
anova(var.ac.nf, type = "marginal")

#refit with ML to select fixed effects
nf.1 = update(var.nf.2, method = "ML")

#remove three-way interaction to test significance
nf.2 = update(nf.1, log(posno) ~ year + farm + site + year:site + year:farm + farm:site)
anova(nf.1, nf.2)
#no significant difference in model fit

#remove two-way interaction to test significance

#remove the year:site interaction
nf.3 = update(nf.1, log(posno) ~ year + farm + site + year:farm + farm:site)

#remove the year:farm interaction
nf.4 = update(nf.1, log(posno) ~ year + farm + site + year:site + farm:site)

#remove the farm:site interaction
nf.5 = update(nf.1, log(posno) ~ year + farm + site + year:site + year:farm)

anova(nf.2, nf.3)
anova(nf.2, nf.4)
anova(nf.2, nf.5)

#significant differences when year:farm interaction is removed
#nf.6 now full model
nf.6 = update(nf.1, log(posno) ~ year + farm + site + year:farm)
#all effects significant.

#refit with REML for obtaining summary statistics
nf.fin = update(nf.6, method = 'REML')

anova(nf.fin, type = "marginal")
summary(nf.fin)

plot(nf.fin)
qqnorm(nf.fin)

#pairwise comparisons

#create models for all single factor and interaction terms

#create models for all single factor and interaction terms
MCtsnitr.farm <- update(nf.fin, log(posno) ~ farm)
MCtsnitr.year <- update(nf.fin, log(posno) ~ year)
MCtsnitr.site <- update(nf.fin, log(posno) ~ site)

MCtsnitr.farm_year <- update(nf.fin, log(posno) ~ farm_year)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.

summary((glht(MCtsnitr.farm, linfct = mcp(farm = "Tukey"))))
summary((glht(MCtsnitr.site, linfct = mcp(site = "Tukey"))))
summary((glht(MCtsnitr.year, linfct = mcp(year = "Tukey"))))

summary((glht(MCtsnitr.farm_year, linfct = mcp(farm_year = "Tukey"))))

#####################
#Seasonal N2O fluxes#
#####################

gls.snf = gls(log(totn_g) ~ year + farm + mgmt + year:mgmt + year:farm + farm:mgmt + year:farm:mgmt, na.action = na.omit, data = all.calc.4)
lme.snf = lme(fixed = log(totn_g) ~ year + farm + mgmt + year:mgmt + year:farm + farm:mgmt + year:farm:mgmt, random = ~1|ID, na.action = na.omit, data = all.calc.4)
anova(gls.snf, lme.snf)
#random effect does not improve overall model fit

var.snf.1 = update(gls.snf, weights = varIdent(form = ~ 1 | year))
var.snf.2 = update(gls.snf, weights = varIdent(form = ~ 1 | farm))
var.snf.3 = update(gls.snf, weights = varIdent(form = ~ 1 | mgmt))
var.snf.4 = update(gls.snf, weights = varIdent(form = ~ 1 | year*farm))
var.snf.5 = update(gls.snf, weights = varIdent(form = ~ 1 | year*mgmt))
var.snf.6 = update(gls.snf, weights = varIdent(form = ~ 1 | year*farm*mgmt))

BIC(gls.snf, var.snf.1, var.snf.2, var.snf.3, var.snf.4, var.snf.5, var.snf.6)
#variance structure of year*mgmt improves overall model fit

#examine fit of fixed effects
anova(var.snf.5, type = "marginal")

#refit with ML to select fixed effects
snf.1 = update(var.snf.5, method = "ML")

#remove three-way interaction to test significance
snf.2 = update(snf.1, log(totn_g) ~ year + farm + mgmt + year:mgmt + year:farm + farm:mgmt)
anova(snf.1, snf.2)
#no significant difference in model fit


#remove two-way interaction to test significance

#remove the year:mgmt interaction
snf.3 = update(snf.1, log(totn_g) ~ year + farm + mgmt + year:farm + farm:mgmt)

#remove the year:farm interaction
snf.4 = update(snf.1, log(totn_g) ~ year + farm + mgmt + year:mgmt + farm:mgmt)

#remove the farm:mgmt interaction
snf.5 = update(snf.1, log(totn_g) ~ year + farm + mgmt + year:mgmt + year:farm)

anova(snf.2, snf.3)
anova(snf.2, snf.4)
anova(snf.2, snf.5)

#difference when year:mgmt is removed, but there are no differences when year:farm and farm:site removed
#all effects now significant. 
snf.6 = update(snf.1, log(totn_g) ~ year + farm + mgmt + year:mgmt)

#refit with REML for obtaining summary statistics
snf.fin = update(snf.6, method = 'REML')

anova(snf.fin, type = "marginal")
summary(snf.fin)

plot(snf.fin)
qqnorm(snf.fin)

#pairwise comparisons

#create models for all single factor and interaction terms

MCseanitr.farm <- update(snf.fin, log(totn_g) ~ farm)
MCseanitr.mgmt <- update(snf.fin, log(totn_g) ~ mgmt)
MCseanitr.year <- update(snf.fin, log(totn_g) ~ year)
MCseanitr.mgmt_year <- update(snf.fin, log(totn_g) ~ year_mgmt)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.

summary((glht(MCseanitr.farm, linfct = mcp(farm = "Tukey"))))
summary((glht(MCseanitr.mgmt, linfct = mcp(mgmt = "Tukey"))))
summary((glht(MCseanitr.year, linfct = mcp(year = "Tukey"))))

summary((glht(MCseanitr.mgmt_year, linfct = mcp(year_mgmt = "Tukey"))))

##############################
#Repeated measures moisture  #
##############################

lme.moist = gls(moist ~ year + farm + site + year:site + year:farm + farm:site + year:farm:site, na.action = na.omit, data = flux.calc)
lme.moist = lme(fixed = moist ~ year + farm + site + year:site + year:farm + farm:site + year:farm:site, random = ~1|ID2, na.action = na.omit, data = flux.calc)
anova(gls.moist, lme.moist)
#random effect improves overall model fit

var.moist.1 = update(lme.moist, weights = varIdent(form = ~ 1 | year))
var.moist.2 = update(lme.moist, weights = varIdent(form = ~ 1 | farm))
var.moist.3 = update(lme.moist, weights = varIdent(form = ~ 1 | site))
var.moist.4 = update(lme.moist, weights = varIdent(form = ~ 1 | year*farm))
var.moist.5 = update(lme.moist, weights = varIdent(form = ~ 1 | year*site))
var.moist.6 = update(lme.moist, weights = varIdent(form = ~ 1 | year*farm*site))

BIC(lme.moist, var.moist.1, var.moist.2, var.moist.3, var.moist.4, var.moist.5, var.moist.6)
#variance structure of farm improves overall model fit

#examine the effect of an autocorrelation structure
ac.moist = update(lme.moist, correlation = corAR1(form = ~ ydoy|ID2))
BIC(lme.moist, ac.moist)
#autocorrelation improves model fit

#combine variance and autocorrelation structures
var.ac.moist = update(lme.moist, weights = varIdent(form = ~ 1| farm), correlation = corAR1(form = ~ ID2|ydoy))
BIC(gls.moist, var.moist.2, ac.moist, var.ac.moist)
#model fit improved with both variance and autocorrelation structures

#examine fit of fixed effects
anova(var.moist.1, type = "marginal")

#refit with ML to select fixed effects
moist.1 = update(var.moist.1, method = "ML")

#remove three-way interaction to test significance
moist.2 = update(moist.1, moist ~ year + farm + site + year:site + year:farm + farm:site)
anova(moist.1, moist.2)
#no significant difference in model fit

#remove two-way interaction to test significance

#remove the year:site interaction
moist.3 = update(moist.1, moist ~ year + farm + site + year:farm + farm:site)

#remove the year:farm interaction
moist.4 = update(moist.1, moist ~ year + farm + site + year:site + farm:site)

#remove the farm:site interaction
moist.5 = update(moist.1, moist ~ year + farm + site + year:site + year:farm)

anova(moist.2, moist.3)
anova(moist.2, moist.4)
anova(moist.2, moist.5)

#significant differences when year:farm interaction is removed
#moist.6 now full model
moist.6 = update(moist.1, moist ~ year + farm + site + year:farm)
#all effects significant.

#refit with REML for obtaining summary statistics
moist.fin = update(moist.6, method = 'REML')

anova(moist.fin, type = "marginal")
summary(moist.fin)

plot(moist.fin)
qqnorm(moist.fin)

#pairwise comparisons

#create models for all single factor and interaction terms

#create models for all single factor and interaction terms
MCtsnitr.farm <- update(moist.fin, moist ~ farm)
MCtsnitr.year <- update(moist.fin, moist ~ year)
MCtsnitr.site <- update(moist.fin, moist ~ site)

MCtsnitr.farm_year <- update(moist.fin, moist ~ farm_year)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.

summary((glht(MCtsnitr.farm, limoistct = mcp(farm = "Tukey"))))
summary((glht(MCtsnitr.site, limoistct = mcp(site = "Tukey"))))
summary((glht(MCtsnitr.year, limoistct = mcp(year = "Tukey"))))

summary((glht(MCtsnitr.farm_year, limoistct = mcp(farm_year = "Tukey"))))



########################################
#SI Table 2 (Instaneous Fluxes)
########################################

flux.calcdt = data.table(flux.calc)
flux.calcdt$id3 = paste(flux.calcdt$farm, flux.calcdt$year, flux.calcdt$site, sep = " ")

sitab2 = flux.calcdt[ , list(mfinco = mean(finco, na.rm = T), lofinco = quantile(finco, 0.025, na.rm = T),
                     hifinco = quantile(finco, 0.975, na.rm = T),
                     mposno = mean(posno, na.rm = T), loposno = quantile(posno, 0.025, na.rm = T),
                     hiposno = quantile(posno, 0.975, na.rm = T),
                     mt_5cm = mean(t_5cm, na.rm = T), lot_5cm = quantile(t_5cm, 0.025, na.rm = T),
                     hit_5cm = quantile(t_5cm, 0.975, na.rm = T),
                     mt_10cm = mean(t_10cm, na.rm = T), lot_10cm = quantile(t_10cm, 0.025, na.rm = T),
                     hit_10cm = quantile(t_10cm, 0.975, na.rm = T),
                     mmoist = mean(moist, na.rm = T), lomoist = quantile(moist, 0.025, na.rm = T),
                     himoist = quantile(moist, 0.975, na.rm = T)), by = id3]

#separate out just the mean values and rename
sitab2.mn = sitab2[ , c("id3", "mfinco", "mposno", "mt_5cm", "mt_10cm", "mmoist")]
names(sitab2.mn) = c("id3", "finco", "posno", "t_5cm", "t_10cm", "moist")

#separate out lower quartile
sitab2.lo = sitab2[ , c("id3", "lofinco", "loposno", "lot_5cm", "lot_10cm", "lomoist")]
names(sitab2.lo) = c("id3", "finco", "posno", "t_5cm", "t_10cm", "moist")

#separate out upper quartile
sitab2.hi = sitab2[ , c("id3", "hifinco", "hiposno", "hit_5cm", "hit_10cm", "himoist")]
names(sitab2.hi) = c("id3", "finco", "posno", "t_5cm", "t_10cm", "moist")

#melt tables and add column names
sitab2.mnmel = melt(sitab2.mn, id.vars = "id3")
names(sitab2.mnmel) = c("id3", "var", "avg")

sitab2.lomel = melt(sitab2.lo, id.vars = "id3")
names(sitab2.lomel) = c("id3", "var", "lo")

sitab2.himel = melt(sitab2.hi, id.vars = "id3")
names(sitab2.himel) = c("id3", "var", "hi")

#round avg, lo, and hi values to two sig. digs.
sitab2.mnmel$avg = round(sitab2.mnmel$avg, 2)
sitab2.lomel$lo = round(sitab2.lomel$lo, 2)
sitab2.himel$hi = round(sitab2.himel$hi, 2)

#make new id column that contains Farm, Management, and var
sitab2.mnmel$id = paste(sitab2.mnmel$id3, sitab2.mnmel$var, sep = " ")
sitab2.lomel$id = paste(sitab2.lomel$id3, sitab2.lomel$var, sep = " ")
sitab2.himel$id = paste(sitab2.himel$id3, sitab2.himel$var, sep = " ")

#remove id3 and year columns
sitab2.mnmel = sitab2.mnmel[ , -c(1:2)]
sitab2.lomel = sitab2.lomel[ , -c(1:2)]
sitab2.himel = sitab2.himel[ , -c(1:2)]

#merge tables together
sitab2.1 = merge(sitab2.mnmel, sitab2.lomel, by.x = "id", by.y = "id")
sitab2.2 = merge(sitab2.1, sitab2.himel, by.x = "id", by.y = "id")

#make new column that contains the mean +- the upper and lower quartiles
sitab2.2$all = paste(sitab2.2$avg, " (", sitab2.2$lo, ", ", sitab2.2$hi, ")", sep = "")

#split the id string
sps <- data.frame(do.call(rbind, str_split(sitab2.2$id, " ")))
names(sps) <- c("farm", "year", "mgmt", "var")

#add to dataframe
sitab2.3 = cbind(sps, sitab2.2)

#remove columns for id, avg, lo, and hi
sitab2.3 = sitab2.3[ , -c(5:8)]

#select rows that are just FRF Graze
sitab2.4 = sitab2.3[sitab2.3$farm == "1" & sitab2.3$mgmt == "1", ]
names(sitab2.4) = c("farm", "year", "mgmt", "var", "FRF Graze")

#add columns for 1 2, etc.
sitab2.5 = cbind(sitab2.4, "FRF Hay" = sitab2.3[sitab2.3$farm == "1" & sitab2.3$mgmt == "2",]$all, 
               "ORG Graze" = sitab2.3[sitab2.3$farm == "2" & sitab2.3$mgmt == "1",]$all,
               "ORG Hay" = sitab2.3[sitab2.3$farm == "2" & sitab2.3$mgmt == "2",]$all,
               "WNF Graze" = sitab2.3[sitab2.3$farm == "3" & sitab2.3$mgmt == "1",]$all,
               "WNF Hay" = sitab2.3[sitab2.3$farm == "3" & sitab2.3$mgmt == "2",]$all)

#export table

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\USDA_ORG\\R Projects\\All-Farm-Tradeoffs\\Organic-Dairy\\Data")
write.table(sitab2.5, "SI_Table_2.csv", col.names = T, row.names = F, sep = ",")


########################################
#SI Table 3 (Seasonal Fluxes)
########################################

all.calc.4dt = data.table(all.calc.4)
all.calc.4dt$id3 = paste(all.calc.4dt$farm, all.calc.4dt$year, all.calc.4dt$mgmt, sep = " ")

sitab3 = all.calc.4dt[ , list(mtotc_g = mean(totc_g, na.rm = T), lototc_g = quantile(totc_g, 0.025, na.rm = T),
                             hitotc_g = quantile(totc_g, 0.975, na.rm = T),
                             mtotn_g = mean(totn_g, na.rm = T), lototn_g = quantile(totn_g, 0.025, na.rm = T),
                             hitotn_g = quantile(totn_g, 0.975, na.rm = T)), by = id3]

#separate out just the mean values and rename
sitab3.mn = sitab3[ , c("id3", "mtotc_g", "mtotn_g")]
names(sitab3.mn) = c("id3", "totc_g", "totn_g")

#separate out lower quartile
sitab3.lo = sitab3[ , c("id3", "lototc_g", "lototn_g")]
names(sitab3.lo) = c("id3", "totc_g", "totn_g")

#separate out upper quartile
sitab3.hi = sitab3[ , c("id3", "hitotc_g", "hitotn_g")]
names(sitab3.hi) = c("id3", "totc_g", "totn_g")

#melt tables and add column names
sitab3.mnmel = melt(sitab3.mn, id.vars = "id3")
names(sitab3.mnmel) = c("id3", "var", "avg")

sitab3.lomel = melt(sitab3.lo, id.vars = "id3")
names(sitab3.lomel) = c("id3", "var", "lo")

sitab3.himel = melt(sitab3.hi, id.vars = "id3")
names(sitab3.himel) = c("id3", "var", "hi")

#round avg, lo, and hi values to two sig. digs.
sitab3.mnmel$avg = round(sitab3.mnmel$avg, 2)
sitab3.lomel$lo = round(sitab3.lomel$lo, 2)
sitab3.himel$hi = round(sitab3.himel$hi, 2)

#make new id column that contains Farm, Management, and var
sitab3.mnmel$id = paste(sitab3.mnmel$id3, sitab3.mnmel$var, sep = " ")
sitab3.lomel$id = paste(sitab3.lomel$id3, sitab3.lomel$var, sep = " ")
sitab3.himel$id = paste(sitab3.himel$id3, sitab3.himel$var, sep = " ")

#remove id3 and year columns
sitab3.mnmel = sitab3.mnmel[ , -c(1:2)]
sitab3.lomel = sitab3.lomel[ , -c(1:2)]
sitab3.himel = sitab3.himel[ , -c(1:2)]

#merge tables together
sitab3.1 = merge(sitab3.mnmel, sitab3.lomel, by.x = "id", by.y = "id")
sitab3.2 = merge(sitab3.1, sitab3.himel, by.x = "id", by.y = "id")

#make new column that contains the mean +- the upper and lower quartiles
sitab3.2$all = paste(sitab3.2$avg, " (", sitab3.2$lo, ", ", sitab3.2$hi, ")", sep = "")

#split the id string
sps <- data.frame(do.call(rbind, str_split(sitab3.2$id, " ")))
names(sps) <- c("farm", "year", "mgmt", "var")

#add to dataframe
sitab3.3 = cbind(sps, sitab3.2)

#remove columns for id, avg, lo, and hi
sitab3.3 = sitab3.3[ , -c(5:8)]

#select rows that are just FRF Graze
sitab3.4 = sitab3.3[sitab3.3$farm == "1" & sitab3.3$mgmt == "G", ]
names(sitab3.4) = c("farm", "year", "mgmt", "var", "FRF Graze")

#add columns for 1 H, etc.
sitab3.5 = cbind(sitab3.4, "FRF Hay" = sitab3.3[sitab3.3$farm == "1" & sitab3.3$mgmt == "H",]$all, 
                 "ORG Graze" = sitab3.3[sitab3.3$farm == "2" & sitab3.3$mgmt == "G",]$all,
                 "ORG Hay" = sitab3.3[sitab3.3$farm == "2" & sitab3.3$mgmt == "H",]$all,
                 "WNF Graze" = sitab3.3[sitab3.3$farm == "3" & sitab3.3$mgmt == "G",]$all,
                 "WNF Hay" = sitab3.3[sitab3.3$farm == "3" & sitab3.3$mgmt == "H",]$all)

#export table

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\USDA_ORG\\R Projects\\All-Farm-Tradeoffs\\Organic-Dairy\\Data")
write.table(sitab3.5, "SI_Table_3.csv", col.names = T, row.names = F, sep = ",")


################################################################################
################################################################################
#Make Figures


#calculate average seasonal fluxes to add to boxplots
ac.5 = data.table(all.calc.4)

ac.5$id = paste(ac.5$farm, ac.5$year, ac.5$mgmt, sep = " ")

mseas = ac.5[ , list(farm = unique(farm), year = unique(year), mgmt = unique(mgmt), 
                     carb = mean(totc_g), nitr = mean(totn_g)), by = id]

mseas = mseas[order(mseas$farm, mseas$year, mseas$mgmt), ]
#Seasonal CO2 fluxes 

#setwd
setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\USDA_ORG\\R Projects\\All-Farm-Tradeoffs\\Organic-Dairy\\Figures")

#name pdf export file
#call pdf.options to define graphical parameters in pdf export file
pdf.options(width= 4.5, height= 6.5, paper="letter", pointsize=10)

#name pdf export file
pdf(file="Seasonal_Fluxes.pdf")

par(mfrow = c(2,1), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(5,5,5,5))

(bx = boxplot(log(totc_g) ~ mgmt * year * farm, data = all.calc.4, axes = F, ylab = " ", ylim = log(c(25, 110)),# log = "y",
        xlab = " ", col = c("coral", "gray80")))
box(lty = 1)

rect(c(3-0.4, 4-0.4, 7-0.4, 8-0.4, 11-0.4, 12-0.4), # left margin of box, position at x-axis 
     bx$stats[2,                   # lower limit of box, may not need to change
              c(3,4,7,8,11,12)],      # position of box at x-axis, subject to change
     c(3+0.4, 4+0.4, 7+0.4, 8+0.4, 11+0.4, 12+0.4),  # right margin of box, position at x-axis 
     bx$stats[4,  # upper limit of box, may not need to change
              c(3,4,7,8,11,12)],      # position of box at x-axis
     density=12, angle=45)
points(log(mseas$carb), pch = 23, bg = "black")

axis(2, at = log(c(25, 35, 50, 70, 100)), lab = c(25, 35, 50, 70, 100))
text(0.55, log(109.5), "(A)")
abline(v = 4.5, lty = 2, lwd = 0.5)
abline(v = 8.5, lty = 2, lwd = 0.5)

braces(xfrom = 0.65, xto = 4.35, yfrom = 4.55, yto = 4.65)
braces(xfrom = 4.65, xto = 8.35, yfrom = 4.55, yto = 4.65)
braces(xfrom = 8.65, xto = 12.35, yfrom = 4.55, yto = 4.65)

text(2.5, 4.7, "a", cex = 0.75)
text(6.5, 4.7, "b", cex = 0.75)
text(10.5, 4.7, "a", cex = 0.75)

braces(xfrom = 0.65, xto = 2.35, yfrom = 4, yto = 4.1)
braces(xfrom = 4.65, xto = 6.35, yfrom = 4, yto = 4.1)
braces(xfrom = 8.65, xto = 10.35, yfrom = 4.3, yto = 4.4)

braces(xfrom = 2.65, xto = 4.35, yfrom = 4.35, yto = 4.45)
braces(xfrom = 6.65, xto = 8.35, yfrom = 4.35, yto = 4.45)
braces(xfrom = 10.65, xto = 12.35, yfrom = 4.3, yto = 4.4)

text(1.5, log(65), "d", cex = 0.75)
text(5.5, log(65), "d", cex = 0.75)
text(9.5, log(85), "d", cex = 0.75)

text(3.5, log(90), "e", cex = 0.75)
text(7.5, log(90), "e", cex = 0.75)
text(11.5, log(90), "e", cex = 0.75)


#mtext(side = 2, expression(Total~CO[2]* Flux ~ (t ~ CO[2]*-C~ha^{-1}~season^{-1})), cex = 1,
#line = 2.5, outer = F)

legend("bottomleft", legend=c("Grazed 2016", "Hayed 2016", "Grazed 2017", "Hayed 2017"), 
       fill=c("coral", "gray80", "coral", "gray80"),
       density=c(NA, NA, NA, NA), bty="n", border = "black", cex = 0.65) 
legend("bottomleft", legend=c("Grazed 2016", "Hayed 2016", "Grazed 2017", "Hayed 2017"), 
       fill=c("coral", "gray80", "black", "black"),
       density=c(NA, NA, 20, 20), bty="n",border=c("black"), cex = 0.65) 

(bx = boxplot(log(totn_g) ~ mgmt * year * farm, data = all.calc.4, axes = F, ylab = " ", ylim = log(c(0.004, 0.21)),# log = "y",
              xlab = " ", col = c("coral", "gray80")))
box(lty = 1)

rect(c(3-0.4, 4-0.4, 7-0.4, 8-0.4, 11-0.4, 12-0.4), # left margin of box, position at x-axis 
     bx$stats[2,                   # lower limit of box, may not need to change
              c(3,4,7,8,11,12)],      # position of box at x-axis, subject to change
     c(3+0.4, 4+0.4, 7+0.4, 8+0.4, 11+0.4, 12+0.4),  # right margin of box, position at x-axis 
     bx$stats[4,  # upper limit of box, may not need to change
              c(3,4,7,8,11,12)],      # position of box at x-axis
     density=12, angle=45)
points(log(mseas$nitr), pch = 23, bg = "black")

text(0.55, log(0.20), "(B)")
abline(v = 4.5, lty = 2, lwd = 0.5)
abline(v = 8.5, lty = 2, lwd = 0.5)

axis(2, at = log(c(0.005, 0.025, 0.100)), lab = c("0.005", "0.025", "0.100"))
axis(1, at = c(2.5, 6.5, 10.5), lab = c("VT", "NH", "ME"))

braces(xfrom = 0.65, xto = 4.35, yfrom = log(0.13), yto =log(0.18))
braces(xfrom = 4.65, xto = 8.35, yfrom = log(0.13), yto =log(0.18))
braces(xfrom = 8.65, xto = 12.35, yfrom = log(0.13), yto =log(0.18))

text(2.5, log(0.205), "a", cex = 0.75)
text(6.5, log(0.205), "b", cex = 0.75)
text(10.5, log(0.205), "b", cex = 0.75)

text(1, log(0.022), "x", cex = 0.75)
text(3, log(0.050), "x", cex = 0.75)
text(5, log(0.040), "x", cex = 0.75)
text(7, log(0.060), "x", cex = 0.75)
text(9, log(0.055), "x", cex = 0.75)
text(11, log(0.135), "x", cex = 0.75)

text(2, log(0.05), "y", cex = 0.75)
text(4, log(0.018), "y", cex = 0.75)
text(5.65, log(0.08), "y", cex = 0.75)
text(8, log(0.022), "y", cex = 0.75)
text(10, log(0.035), "y", cex = 0.75)
text(12, log(0.035), "y", cex = 0.75)

mtext(side = 2, expression(Total~CO[2]*~or~N[2]*O ~ Flux ~ (t ~CO[2]*-C~or~N[2]*O-N ~ ha^{-1} ~season^{-1})), cex = 1,
      line = 2.5, outer = T) 

dev.off()
