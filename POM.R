setwd("C:/Users/KUA9TE/Desktop")

# brislte data
df1 <- read.csv("Workbook1.csv", header = TRUE)
is.factor(df1$Phenotype)
df1$Phenotype <- factor(df1$Phenotype, ordered = TRUE, 
                       levels = c("Normal", "Extra1", "Extra2", "Extra3"))


library(MASS)
POM1 <- polr(Phenotype ~ Genotype, 
             weights = Count, data = df1, Hess = TRUE)
summary(POM1)


ctable <- coef(summary(POM1))
p_value <- pnorm(abs(ctable[, "t value"]), 
                 lower.tail = FALSE) *2
ctable <- cbind(ctable, "p value" = p_value)
ctable

# wing nick
df2 <- read.csv("Workbook2.csv", header = T)
is.factor(df2$Phenotype)
df2$Phenotype <- factor(df2$Phenotype, ordered = T,
                        levels = c("Normal", "Nick1", "Nick2", "Nick3"))

POM2 <- polr(Phenotype ~ Genotype, weights = Count, df2, Hess = TRUE)
summary(POM2)

ctable2 <- coef(summary(POM2))
p_value2 <- pnorm(abs(ctable2[, "t value"]), 
                 lower.tail = FALSE) *2
ctable2 <- cbind(ctable2, "p value" = p_value2)
ctable2

predict(pom, newdata = data.frame(party = "Dem"), type = "p")
predict(pom, newdata = data.frame(party = "Rep"), type = "p")

# Democrat cumulative probabilities
Cdk8_cp <- exp(POM2$zeta - POM2$coefficients)/(1 + exp(POM2$zeta - POM2$coefficients))
# Republican cumulative probabilities
N_cp <- exp(POM2$zeta)/(1 + exp(POM2$zeta))
N_cp

# Democrat odds
Cdk_odds <- (dcp/(1-dcp))
dodds
# Repuclican odds
N_odds <- (N_cp/(1-N_cp))
N_odds




N <- c(111,	36,	13,	2)
NxCdk8 <- c(80,	18,	2,	0)
sum(N)
sum(NxCdk8)
pheno <- c("Normal", "nick1", "nick2", "nick3")
data <- factor(c(rep(pheno, N), 
                 rep(pheno, NxCdk8)), levels = pheno)
Geno <- factor(rep(c("N","NxCdk8"), c(162, 100)), 
                levels=c("N","NxCdk8"))  
dataf <- data.frame(Geno, data)

library(MASS)
pom_cdk8 <- polr(data ~ Geno, data = dataf)
summary(pom_cdk8)

coef <- coef(summary(pom_cdk8))
p_cdk8 <- pnorm(abs(coef[, "t value"]), 
                  lower.tail = FALSE) *2
p_cdk8

# Democrat cumulative probabilities
Cdk8_cp <- exp(pom_cdk8$zeta - pom_cdk8$coefficients)/(1 + exp(pom_cdk8$zeta - pom_cdk8$coefficients))
# Republican cumulative probabilities
N_cp <- exp(pom_cdk8$zeta)/(1 + exp(pom_cdk8$zeta))
N_cp

# Democrat odds
Cdk8_odds <- (Cdk8_cp/(1-Cdk8_cp))
Cdk8_odds
# Repuclican odds
N_odds <- (N_cp/(1-N_cp))
N_odds

# odds ratio
OR <- N_odds/Cdk8_odds
OR
