setwd("C:/Users/KUA9TE/Desktop")


# bristle test N x med/ago -------------------------------------------

# load bristle data
df_bristle <- read.csv("Workbook1.csv", 
                header = TRUE)
is.factor(df_brislte$Phenotype)
df_bristle$Phenotype <- factor(df_bristle$Phenotype, 
                        ordered = TRUE, 
                        levels = c("Normal", "Extra1", "Extra2", "Extra3"))


library(MASS)
pom_bristle <- polr(Phenotype ~ Genotype, 
             weights = Count, 
             data = df_bristle, 
             Hess = TRUE)
summary(pom_bristle)

# Wald–statistic test
se <- summary(pom_bristle)[1][[1]][1:8,2] # get standard error
wald <- -coef(pom_bristle)/se 
1-pchisq(wald^2, df=1) # not sure if df=1


# fit multinomial logistic model
mlm <- multinom(Phenotype ~ Genotype, data = df_bristle)
M1 <- logLik(pom_bristle)
M2 <- logLik(mlm)
G <- -2*(M1[1] - M2[1])
pchisq(G, 3,  lower.tail = FALSE) # not sure if df=3

# odds ratio
or_brislte <- exp(-coef(pom_bristle))
or_brislte

# obtain raw p value
ctable_bristle <- coef(summary(pom_bristle))
p_value_bristle <- pnorm(abs(ctable_bristle[, "t value"]), 
                 lower.tail = FALSE) *2
ctable_bristle <- cbind(ctable_bristle, 
                "p value" = p_value_bristle)
ctable_bristle

# p value adjustment
library(FSA)
geno <- row.names(ctable_bristle)
p_value_bristle <- data.frame(geno, p_value_bristle)
headtail(p_value_bristle)

# perform p value adjustment and add to dataframe
p_value_bristle$Bonferroni <- p.adjust(p_value_bristle$p_value_bristle, 
                                method = "bonferroni")
p_value_bristle$BH <- p.adjust(p_value_bristle$p_value_bristle, 
                        method = "BH")
p_value_bristle$Holm <- p.adjust(p_value_bristle$p_value_bristle, 
                          method = "holm")
p_value_bristle$Hommel <- p.adjust(p_value_bristle$p_value_bristle, 
                            method = "hommel")
p_value_bristle$BY <- p.adjust(p_value_bristle$p_value_bristle, 
                        method = "BY")
p_value_bristle


# plot
X = p_value_bristle$p_value_bristle
Y = cbind(p_value_bristle$Bonferroni,
          p_value_bristle$BH,
          p_value_bristle$Holm,
          p_value_bristle$Hochberg,
          p_value_bristle$Hommel,
          p_value_bristle$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2)

legend('bottomright', 
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"), 
       col = 1:6, 
       cex = 1,    
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)



# wing nicking test N x med/ago --------------------------------------


# wing nick
df_nick <- read.csv("Workbook2.csv", 
                header = T)
is.factor(df_nick$Phenotype)
df_nick$Phenotype <- factor(df_nick$Phenotype, 
                        ordered = T,
                        levels = c("Normal", "Nick1", "Nick2", "Nick3"))

library(MASS)
pom_nick <- polr(Phenotype ~ Genotype, 
             weights = Count, 
             df_nick, 
             Hess = TRUE)
summary(pom_nick)

ctable_nick <- coef(summary(pom_nick))
p_value_nick <- pnorm(abs(ctable_nick[, "t value"]), 
                  lower.tail = FALSE) *2
ctable_nick <- cbind(ctable_nick, 
                 "p value" = p_value_nick)
ctable_nick



# p value adjustment
library(FSA)
geno <- row.names(ctable_nick)
p_value_nick <- data.frame(geno, p_value_nick)
headtail(p_value_nick)

# perform p value adjustment and add to dataframe
p_value_nick$Bonferroni <- p.adjust(p_value_nick$p_value_nick, 
                                    method = "bonferroni")
p_value_nick$BH <- p.adjust(p_value_nick$p_value_nick, 
                            method = "BH")
p_value_nick$Holm <- p.adjust(p_value_nick$p_value_nick, 
                              method = "holm")
p_value_nick$Hommel <- p.adjust(p_value_nick$p_value_nick, 
                                method = "hommel")
p_value_nick$BY <- p.adjust(p_value_nick$p_value_nick, 
                            method = "BY")
p_value_nick


# plot
X = p_value_nick$p_value_nick
Y = cbind(p_value_nick$Bonferroni,
          p_value_nick$BH,
          p_value_nick$Holm,
          p_value_nick$Hochberg,
          p_value_nick$Hommel,
          p_value_nick$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2)

legend('bottomright', 
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"), 
       col = 1:6, 
       cex = 1,    
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

# compare N vs Nxcdk8hets only
df_cdk8nick <- df_nick[1:8,]
is.factor(df_cdk8nick$Phenotype)
df_cdk8nick$Phenotype <- factor(df_cdk8nick$Phenotype, 
                            ordered = T,
                            levels = c("Normal", "Nick1", "Nick2", "Nick3"))

library(MASS)
pom_cdk8nick <- polr(Phenotype ~ Genotype, 
                 weights = Count, 
                 df_cdk8nick, 
                 Hess = TRUE)
summary(pom_cdk8nick)

ctable_cdk8 <- coef(summary(pom_cdk8nick))
p_value_cdk8 <- pnorm(abs(ctable_cdk8[, "t value"]), 
                      lower.tail = FALSE) *2
ctable_cdk8 <- cbind(ctable_cdk8, 
                     "p value" = p_value_cdk8)
ctable_cdk8


# 1 variable test (N vs NxCdk8) ---------------------------------------------------------


N <- c(111,	36,	13,	2)
NxCdk8 <- c(80,	18,	2,	0)
sum(N)
sum(NxCdk8)
pheno <- c("Normal", "nick1", "nick2", "nick3")
data <- factor(c(rep(pheno, N), 
                 rep(pheno, NxCdk8)), 
               levels = pheno)
Geno <- factor(rep(c("N","NxCdk8"), c(162, 100)), 
               levels=c("N","NxCdk8"))  
dataf <- data.frame(Geno, data)

library(MASS)
pom_cdk8 <- polr(data ~ Geno, 
                 data = dataf)
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




# Chi square test ---------------------------------------------------------


df3 <- read.csv("Book3.csv")


# perform chi-square test for wing nicking data

chisq.test(df3$Cdk8, p = df3$N, rescale.p = TRUE)
chisq.test(df3$cycC, p = df3$N, rescale.p = TRUE)
chisq.test(df3$kto241, p = df3$N, rescale.p = TRUE)
chisq.test(df3$kto631, p = df3$N, rescale.p = TRUE)
chisq.test(df3$skd13, p = df3$N, rescale.p = TRUE)
chisq.test(df3$skd413, p = df3$N, rescale.p = TRUE)
chisq.test(df3$ago1, p = df3$N, rescale.p = TRUE)
chisq.test(df3$ago3, p = df3$N, rescale.p = TRUE)


df4 <- read.csv("Book4.csv")
chisq.test(df4$Cdk8, p = df4$N, rescale.p = TRUE)
chisq.test(df4$cycC, p = df4$N, rescale.p = TRUE)
chisq.test(df4$kto241, p = df4$N, rescale.p = TRUE)
chisq.test(df4$kto631, p = df4$N, rescale.p = TRUE)
chisq.test(df4$skd13, p = df4$N, rescale.p = TRUE)
chisq.test(df4$skd413, p = df4$N, rescale.p = TRUE)
chisq.test(df4$ago1, p = df4$N, rescale.p = TRUE)
chisq.test(df4$ago3, p = df4$N, rescale.p = TRUE)
