##############################################################################################
##############################################################################################
##############################################################################################

#This script evaluates whether management intensive grazing (MIG) changes aboveground biomass
#production in three organic dairy farms located across New England.

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
library(multcomp)
library(plyr)
library(multcomp)
library(swfscMisc)

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

bio = read.table("bio.csv", head = TRUE, sep = ",")

##############################################################################################
##############################################################################################
##############################################################################################

#Data preprocessing

#fill in gaps with linear interpolation
bio$totdry_2 = na.approx(bio$totdry, na.rm = F)

#convert dry matter to kg C and kg C / ha
bio$kgha = bio$totdry_2 * (10 / 1) * (10000 / 1000)
bio$kgcha = bio$totdry_2 * (10 / 1) * (10000 / 1000) * 0.45

#export table with interpolated values
write.table(bio, file = paste("bio_int.csv"), sep = ",", na="NA", append=FALSE, col.names=TRUE, row.names = FALSE)

#calculate the cumulative sum by Farm / Field / Rep
bio$ID = paste(bio$Farm, bio$Field, bio$Trt, bio$Rep, sep = " ")
bio$csum = ave(bio$kgcha, bio$ID, FUN = cumsum)

#determine the average field level value per sampling date
#and the max value per farm-field-rep for running stats
bio$ID = paste(bio$Farm, bio$Trt, bio$Rep, bio$Date.Sampled, sep = " ")
bio$ID.2 = paste(bio$Farm, bio$Field, bio$Trt, bio$Rep, sep = " ")

bio.1 = data.table(bio)
bio.2 = bio.1[ , list(Farm = unique(Farm), Field = unique(Field), Trt = unique(Trt), Rep = unique(Rep),
                      Date.Sampled = unique(Date.Sampled), 
                      mcsum = median(csum), maxcsum = quantile(csum, 0.75), mincsum = quantile(csum, 0.25)), by = ID]

#convert date into posix object
bio.2$DATETIME = as.POSIXct(strptime(bio.2$Date.Sampled, "%m/%d/%Y", tz = "EST"))

#reduce data further to summarize by farm
bio.2$ID2 = paste(bio.2$Farm, bio.2$Trt, bio.2$Date.Sampled, sep = " ")
bio.2 = bio.2[ , list(Farm = unique(Farm), Trt = unique(Trt), Date.Sampled = unique(Date.Sampled), mcsum = median(mcsum), 
                     maxcsum = quantile(mcsum, 0.75), mincsum = quantile(mcsum, 0.25)), by = ID2]

#convert date into posix object
bio.2$DATETIME = as.POSIXct(strptime(bio.2$Date.Sampled, "%m/%d/%Y", tz = "EST"))

#extract the last value of the cumsum for analyzing total seasonal growth
bio.3 = bio.1[ , list(Farm = unique(Farm), Field = unique(Field), Trt = unique(Trt), 
                      Rep = unique(Rep), csum = max(csum)), by = ID.2]

##############################################################################################
##############################################################################################
##############################################################################################

#Data characterization (normality, homogeneity of variance, correlation)


#####################
#Normal distribution#
#####################

par(mfrow = c(2,2))

hist(bio$kgcha)
hist(bio.3$csum)

hist(log(bio$kgcha))
hist(log(bio.3$csum))

################
#Equal variance#
################

boxplot(log(kgcha) ~ Trt * Farm, data = bio)
boxplot(log(csum) ~ Trt * Farm, data = bio.3)

##############################################################################################
##############################################################################################
##############################################################################################

#Statistical analysis

#make Field a factor variable and combine Farm and Field for pairwise comparisons
bio$Field = factor(bio$Field)
bio$FarmField = paste(bio$Farm, bio$Field, sep = " ")

bio.3$Field = factor(bio.3$Field)
bio.3$FarmField = paste(bio.3$Farm, bio.3$Field, sep = " ")

###################
#repeated measures#
###################

#model random effects
gls.bio.ts = gls(log(kgcha) ~ Trt + Farm + Trt:Farm, data = bio, na.action = na.omit)
lme.bio.ts.1 = lme(fixed = log(kgcha) ~ Trt + Farm + Trt:Farm, random = ~1 | ID.2, data = bio, na.action = na.omit)

anova(gls.bio.ts, lme.bio.ts.1)
#random intercept does not improve model fit

#model variance structure
var.bio.ts.1 = update(gls.bio.ts, weights = varIdent(form = ~ 1 | Trt))
var.bio.ts.2 = update(gls.bio.ts, weights = varIdent(form = ~ 1 | Farm))
var.bio.ts.3 = update(gls.bio.ts, weights = varIdent(form = ~ 1 | Trt*Farm))

BIC(gls.bio.ts, var.bio.ts.1, var.bio.ts.2, var.bio.ts.3)
#Farm as a constant variance structure improves model fit

#model autocorrelation structure
ac.bio.ts.1 = update(gls.bio.ts, correlation = corAR1(form = ~ 1))
BIC(gls.bio.ts, ac.bio.ts.1)
#autocorrelation structure improves model fit

#combine variance and autocorrelation structures
var.ac.bio.ts = update(gls.bio.ts, weights = varIdent(form = ~ 1| Farm), 
                       correlation = corAR1(form = ~ 1))
BIC(gls.bio.ts, var.bio.ts.2, ac.bio.ts.1, var.ac.bio.ts)
#model with combined variance and autocorrelation structures has best model fit

#examine fit of fixed effects
anova(var.ac.bio.ts, type = "marginal")
#only the effect of Farm is significant

#refit with ML select fixed effects
bio.ts.1 = update(var.ac.bio.ts, method = "ML")

#remove two-way interaction
bio.ts.2 = update(bio.ts.1, log(kgcha) ~ Trt + Farm)
anova(bio.ts.1, bio.ts.2)
#no difference in model fit

#remove each fixed effect in turn
bio.ts.3 = update(bio.ts.1, log(kgcha) ~ Trt)
bio.ts.4 = update(bio.ts.1, log(kgcha) ~ Farm)

anova(bio.ts.2, bio.ts.3)
anova(bio.ts.2, bio.ts.4)

#model fit significantly worse with the effects of Farm and Trt removed. 
#bio.ts.2 final model

bio.ts.fin = update(bio.ts.2, method = "REML")

summary(bio.ts.fin)
anova(bio.ts.fin)

plot(bio.ts.fin)
qqnorm(bio.ts.fin)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.

summary((glht(bio.ts.fin, linfct = mcp(Farm = "Tukey"))))
summary((glht(bio.ts.fin, linfct = mcp(Trt = "Tukey"))))

#########################
#seasonal values

#model random effects
gls.bio = gls(log(csum) ~ Trt + Farm + Trt:Farm, data = bio.3, na.action = na.omit)
lme.bio.1 = lme(fixed = log(csum) ~ Trt + Farm + Trt:Farm, random = ~1 | ID.2, data = bio.3, na.action = na.omit)

anova(gls.bio, lme.bio.1)
#random intercept does not improve model fit

#model variance structure
var.bio.1 = update(gls.bio, weights = varIdent(form = ~ 1 | Trt))
var.bio.2 = update(gls.bio, weights = varIdent(form = ~ 1 | Farm))
var.bio.3 = update(gls.bio, weights = varIdent(form = ~ 1 | Trt * Farm))
AIC(gls.bio, var.bio.1, var.bio.2, var.bio.3)
#Trt as a constant variance structure improves model fit

#examine fit of fixed effects
anova(var.bio.1, type = "marginal")
#only the effect of farm is significant

#Refit with ML to select fixed effects
bio.seas.1 = update(var.bio.1, method = "ML")

#remove two-way interaction
bio.seas.2 = update(bio.seas.1, log(csum) ~ Trt + Farm)
anova(bio.seas.1, bio.seas.2)

#no difference in model fit
#remove each fixed effect in turn
bio.seas.3 = update(bio.seas.1, log(csum) ~ Trt)
bio.seas.4 = update(bio.seas.1, log(csum) ~ Farm)

anova(bio.seas.2, bio.seas.3)
anova(bio.seas.2, bio.seas.4)

#model fit significantly worse with the effect of Farm removed. 
#no significant difference with the effect of Trt removed
#bio.seas.4 final model

bio.fin = update(bio.seas.4, method = "REML")

summary(bio.fin)
anova(bio.fin)

plot(bio.fin)
qqnorm(bio.fin)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.
summary((glht(bio.fin, linfct = mcp(Farm = "Tukey"))))

##############################################################################################
##############################################################################################
##############################################################################################

#Plot Results

#order data by date and make new ID variable for plotting
bio.2 = bio.2[order(bio.2$Farm, bio.2$Trt, bio.2$DATETIME), ]
bio.2$ID = paste(bio.2$Farm, bio.2$Trt, sep = " ")

#setwd
setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\USDA_ORG\\R Projects\\All-Farm-Tradeoffs\\Organic-Dairy\\Figures")

#call pdf.options to define graphical parameters in pdf export file
pdf.options(width= 7.5, height= 3.5, paper="letter", pointsize=12)

#name pdf export file
pdf(file="Figure_2_Revised.pdf")

par(mfrow = c (1,3), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(5,5,5,5))

#FRF
with(bio.2[bio.2$ID == "FRF Graze",], plot(DATETIME, mcsum, type = "n", lwd = 2, lty = 1, col = "coral", yaxt = "n",
                    ylim = c(0,6000), xlab = " ", ylab = " ", main = " ", cex.axis = 1.25, cex.main = 1.5))

with(bio.2[bio.2$ID == "FRF Graze",], lines(DATETIME, mincsum, col = 'coral', lwd = 1, lty = 2))
with(bio.2[bio.2$ID == "FRF Graze",], lines(DATETIME, maxcsum, col = 'coral', lwd = 1, lty = 2))
with(bio.2[bio.2$ID == "FRF Graze",], polygon(c(DATETIME, rev(DATETIME)), c(mincsum, rev(maxcsum)),
                    col = adjustcolor("coral", alpha = 0.15), border = NA))
with(bio.2[bio.2$ID == "FRF Graze",], lines(DATETIME, mcsum, type = "l", lwd = 2, lty = 1, col = "coral"))


with(bio.2[bio.2$ID == "FRF Hay",], lines(DATETIME, mincsum, col = 'gray30', lwd = 1, lty = 2))
with(bio.2[bio.2$ID == "FRF Hay",], lines(DATETIME, maxcsum, col = 'gray30', lwd = 1, lty = 2))
with(bio.2[bio.2$ID == "FRF Hay",], polygon(c(DATETIME, rev(DATETIME)), c(mincsum, rev(maxcsum)),
                    col = adjustcolor("gray30", alpha = 0.15), border = NA))
with(bio.2[bio.2$ID == "FRF Hay",], lines(DATETIME, mcsum, lwd = 2, lty = 1, col = "gray30"))

axis(2, at = c(1000, 3000, 5000), cex.axis = 1.5)
mtext(side = 1, "Sampling Date", cex = 1, line = 3, outer = T)
mtext(side = 2, expression("Biomass Carbon (kg C" *~ha^{-1}*")"), cex = 1, line = 3)

legend( x="topleft", 
        legend=c("Grazed", "Hayed", NA, NA),
        col=c("coral","gray30"), lwd = 2, lty = c(1, 1, NA, NA), bty = "n",
        pch=c(NA,NA, NA, NA), pt.cex = 3, cex = 1.25)

legend( x="topleft", 
        legend=c(NA, NA, NA, NA), 
        col =c(adjustcolor("coral", alpha = 0.15), 
                 adjustcolor("gray30", alpha = 0.15)), lwd=2, lty=c(0, 0, 0, 0), bty = "n",
        pch=c(15,15, NA, NA), pt.cex = 3, cex = 1.25)

#ORG
with(bio.2[bio.2$ID == "ORG Graze",], plot(DATETIME, mcsum, type = "n", lwd = 2, lty = 1, col = "coral", yaxt = "n",
                                           ylim = c(0,6000), xlab = " ", ylab = " ", main = " ", cex.axis = 1.25, cex.main = 1.5))

with(bio.2[bio.2$ID == "ORG Graze",], lines(DATETIME, mincsum, col = 'coral', lwd = 1, lty = 2))
with(bio.2[bio.2$ID == "ORG Graze",], lines(DATETIME, maxcsum, col = 'coral', lwd = 1, lty = 2))
with(bio.2[bio.2$ID == "ORG Graze",], polygon(c(DATETIME, rev(DATETIME)), c(mincsum, rev(maxcsum)),
                                              col = adjustcolor("coral", alpha = 0.15), border = NA))
with(bio.2[bio.2$ID == "ORG Graze",], lines(DATETIME, mcsum, type = "l", lwd = 2, lty = 1, col = "coral"))


with(bio.2[bio.2$ID == "ORG Hay",], lines(DATETIME, mincsum, col = 'gray30', lwd = 1, lty = 2))
with(bio.2[bio.2$ID == "ORG Hay",], lines(DATETIME, maxcsum, col = 'gray30', lwd = 1, lty = 2))
with(bio.2[bio.2$ID == "ORG Hay",], polygon(c(DATETIME, rev(DATETIME)), c(mincsum, rev(maxcsum)),
                                            col = adjustcolor("gray30", alpha = 0.15), border = NA))
with(bio.2[bio.2$ID == "ORG Hay",], lines(DATETIME, mcsum, lwd = 2, lty = 1, col = "gray30"))

#WNF
with(bio.2[bio.2$ID == "WNF Graze",], plot(DATETIME, mcsum, type = "n", lwd = 2, lty = 1, col = "coral", yaxt = "n",
                                           ylim = c(0,6000), xlab = " ", ylab = " ", main = " ", cex.axis = 1.25, cex.main = 1.5))

with(bio.2[bio.2$ID == "WNF Graze",], lines(DATETIME, mincsum, col = 'coral', lwd = 1, lty = 2))
with(bio.2[bio.2$ID == "WNF Graze",], lines(DATETIME, maxcsum, col = 'coral', lwd = 1, lty = 2))
with(bio.2[bio.2$ID == "WNF Graze",], polygon(c(DATETIME, rev(DATETIME)), c(mincsum, rev(maxcsum)),
                                              col = adjustcolor("coral", alpha = 0.15), border = NA))
with(bio.2[bio.2$ID == "WNF Graze",], lines(DATETIME, mcsum, type = "l", lwd = 2, lty = 1, col = "coral"))


with(bio.2[bio.2$ID == "WNF Hay",], lines(DATETIME, mincsum, col = 'gray30', lwd = 1, lty = 2))
with(bio.2[bio.2$ID == "WNF Hay",], lines(DATETIME, maxcsum, col = 'gray30', lwd = 1, lty = 2))
with(bio.2[bio.2$ID == "WNF Hay",], polygon(c(DATETIME, rev(DATETIME)), c(mincsum, rev(maxcsum)),
                                            col = adjustcolor("gray30", alpha = 0.15), border = NA))
with(bio.2[bio.2$ID == "WNF Hay",], lines(DATETIME, mcsum, lwd = 2, lty = 1, col = "gray30"))

dev.off()

