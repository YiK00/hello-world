require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)


setwd("C:/Users/KUA9TE/Desktop")
dat <- read.dta("https://stats.idre.ucla.edu/stat/data/ologit.dta")
head(dat)

# t value is simply the ratio of the coefficient to its standard error.



# EXAMPLE 2 ---------------------------------------------------------------

# ref 
# https://data.library.virginia.edu/fitting-and-interpreting-a-proportional-odds-model/

party <- factor(rep(c("Rep","Dem"), c(407, 428)), 
                levels=c("Rep","Dem"))  
rpi <- c(30, 46, 148, 84, 99) # cell counts
dpi <- c(80, 81, 171, 41, 55) # cell counts
ideology <- c("Very Liberal","Slightly Liberal","Moderate","Slightly Conservative","Very Conservative")
pol.ideology <- factor(c(rep(ideology, rpi), 
                         rep(ideology, dpi)), levels = ideology)
dat <- data.frame(party,pol.ideology)

library(MASS)
pom <- polr(pol.ideology ~ party, data = dat, Hess = TRUE)
summary(pom)

# The type="p" argument says we want probabilities.
predict(pom, newdata = data.frame(party = "Dem"), type = "p")
predict(pom, newdata = data.frame(party = "Rep"), type = "p")

# Democrat cumulative probabilities
dcp <- exp(pom$zeta - pom$coefficients)/(1 + exp(pom$zeta - pom$coefficients))
# Republican cumulative probabilities
rcp <- exp(pom$zeta)/(1 + exp(pom$zeta))

# Democrat odds
dodds <- (dcp/(1-dcp))
dodds
# Repuclican odds
rodds <- (rcp/(1-rcp))
rodds

# “Proportional” means that two ratios are equal.
odds_ratio <- dodds/rodds
odds_ratio

# check the assumption of protional odds model by comparing POM to multinomial logit model.
library(nnet)
mlm <- multinom(pol.ideology ~ party, data = dat)

M1 <- logLik(pom)
M1
M2 <- logLik(mlm)
M2
G <- -2*(M1[1] - M2[1])
G
pchisq(G, 3, lower.tail = FALSE)
