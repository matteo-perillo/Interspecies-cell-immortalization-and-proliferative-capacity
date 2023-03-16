#BEFORE STARTING, RUN THE SCRIPT "01_IMPORT DATA AND TREES", OR LOAD THE
#WORKING ENVIRONMENT "IMPORTED DATA AND TREES"

library(ape)
library(phytools)
library(phangorn)
library(ppcor)
#library(stargazer)
#stargazer(reg_cp_peso, type = "text")
library(ape)
library(car)
library(phytools)
library(geiger)
library(nlme)
library(flexmix)
library(AICcmodavg)

#Scatterplot
pairs(cbind(log(data_pd[,1:3]), data_pd[,4]),
      labels = c("log(mass)", "log(long)", "log(PD)", "BMRr"), pch=16)

#Raw and partial correlations considering also BMRr
rawcor <- cor(cbind(log(data_pd[,1:3]), data_pd[,4])[!is.na(data_pd$BMRr),])
rawcor
partcor = pcor(cbind(log(data_pd[,1:3]), data_pd[,4])[!is.na(data_pd$BMRr),])
partcor$estimate[-3, 3]
partcor$p.value[-3, 3] #statistically significant correlation mass-PD and BMRr-PD

#Linear regression with only intercept ("null model")
reg_pd0 <- lm(log(PD)~1, data = data_pd)
summary(reg_pd0)
plot(log(data_pd[,c(1,3)]), xlab = "ln(mass)", ylab = "ln(PD)", pch=16)
text(x = log(data_pd$mass)+1, y = log(data_pd$PD), labels = rownames(data_pd), cex=0.5)
abline(reg_pd0, col=2, lwd=2)
phylosig(tree, reg_pd0$residuals, method = "lambda", test = T)
#p-val 0.3: there is phylo sig in the residuals of the null model

#Simple linear regressions
#Mass vs PD
reg_pdm <- lm(log(PD)~log(mass), data = data_pd)
summary(reg_pdm) #R2 high (0.68), strongly significant regressor (p=0.0001)
plot(log(data_pd[,c(1,3)]), xlab = "ln(mass)", ylab = "ln(PD)", pch=16)
abline(reg_pdm, col=2, lwd=2)
phylosig(tree, reg_pdm$residuals, method = "lambda", test = T) #pval = 0.56, no phylosig

#Mass vs PD - bis
reg_pdm_2 <- lm(log(PD)~log(mass), data = na.omit(data_pd))
summary(reg_pdm_2)

#Longevity vs PD
reg_pdl <- lm(log(PD)~log(long), data = data_pd)
summary(reg_pdl) #R2 moderate (0.46), regressor correlated (p=0.006)
plot(log(data_pd[,c(2,3)]), xlab = "ln(long)", ylab = "ln(PD)", pch=16)
abline(reg_pdl, col=2, lwd=2)
phylosig(tree, reg_pdl$residuals, method = "lambda", test = T) #pval = 1, no phylosig

#BMRr vs PD
reg_pdb <- lm(log(PD)~BMRr, data = data_pd)
summary(reg_pdb) #R2 moderate/low (0.34), not significant regressor (p=0.08)
plot(data_pd$BMRr, log(data_pd$PD), xlab = "BMRr", ylab = "ln(PD)", pch=16)
abline(reg_pdb, col=2, lwd=2)
phylosig(tree, reg_pdb$residuals, method = "lambda", test = T) #pval = 0.29, no phylosig

#For residual-residual plot
reg_ml <- lm(log(mass)~log(long), data = data_pd)
reg_lm <- lm(log(long)~log(mass), data = data_pd)
reg_mb <- lm(log(mass)~BMRr, data = na.omit(data_pd))
reg_bm <- lm(BMRr~log(mass), data = na.omit(data_pd))

#Partial correlation mass-PD cont long
partreg_pdm <- lm(reg_pdl$residuals~reg_ml$residuals)
summary(partreg_pdm) #significant relation (pval=0.001)

plot(reg_ml$residuals, reg_pdl$residuals, pch = 16,
     xlab = "Residuals of ln(long) vs ln(mass)",
     ylab = "Residuals of ln(long) vs ln(PD)")
#abline(partreg_pdm, col=2, lwd=2)

#Partial correlation long-PD cont mass
partreg_pdl <- lm(reg_pdm$residuals~reg_lm$residuals)
summary(partreg_pdl) #not significant relation (pval=0.064)

plot(reg_lm$residuals, reg_pdm$residuals, pch = 16,
     xlab = "Residuals of ln(mass) vs ln(long)",
     ylab = "Residuals of ln(mass) vs ln(PD)")
#abline(partreg_pdl, col=2, lwd=2)

#Partial correlation mass-PD cont BMRr
partreg_pdm2 <- lm(reg_pdb$residuals~reg_mb$residuals)
summary(partreg_pdm2) #significant relation (pval=0.0002)

plot(reg_mb$residuals, reg_pdb$residuals, pch = 16,
     xlab = "Residuals of BMRr vs ln(mass)",
     ylab = "Residuals of BMRr vs ln(PD)")
#abline(partreg_pdm2, col=2, lwd=2)

#Partial correlation BMRr-PD cont mass
partreg_pdb <- lm(reg_pdm_2$residuals~reg_bm$residuals)
summary(partreg_pdb) #not significant relation (pval=0.005)

plot(reg_bm$residuals, reg_pdm_2$residuals, pch = 16,
     xlab = "Residuals of ln(mass) vs BMRr",
     ylab = "Residuals of ln(mass) vs ln(PD)")
#abline(partreg_pdl, col=2, lwd=2)

#Multiple linear regression
#PD vs Body mass and Longevity
reg_pdmalo <- lm(log(PD)~log(long)+log(mass), data = data_pd)
summary(reg_pdmalo) #pretty high adj-R2 (0.72), mass significant (0.002), long not significant (0.075)

#PD vs Body mass and Longevity - bis
reg_pdmalo_2 <- lm(log(PD)~log(long)+log(mass), data = data_pd[!is.na(data_pd[4]),])
summary(reg_pdmalo_2)

#PD vs Body mass and BMRr
reg_pdmabr <- lm(log(PD)~log(mass)+BMRr, data = data_pd)
summary(reg_pdmabr) #high adj-R2 (0.77), mass (0.003) and BMRr (0.029) significant
#plot(reg_pdmabr)
phylosig(tree, reg_pdmabr$residuals, method = "lambda", test = T) #pval = 1, no phylosig

#PD vs Body mass, Longevity and BMRr
reg_pdfull <- lm(log(PD)~log(long)+log(mass)+BMRr, data = data_pd)
summary(reg_pdfull) #pretty high adj-R2 (0.73), mass significant (0.009), BMRr (0.079) and long (0.989) not sig

#Comparison of different regression models
BIC(reg_pdm)
BIC(reg_pdmalo)

BIC(reg_pdmabr)
BIC(reg_pdfull)
BIC(reg_pdmalo_2)
BIC(reg_pdm_2)

aictab(list(reg_pdm, reg_pdmalo), modnames = c("Mass", "Mass+Longevity"), second.ord = F)
aictab(list(reg_pdmabr, reg_pdfull, reg_pdm_2, reg_pdmalo_2), second.ord = F,
       modnames = c("Mass+BMRr", "Mass+Longevity+BMRr", "Mass sub", "Mass+Longevity sub"))

#Analysis taking phylogenesis into account
sp_pd <- row.names(data_pd)[-(13:14)]

pdm_form = log(PD)~log(mass)
pdform = log(PD)~log(mass)+BMRr

#Pagel lambda phylogenetic signal evolution model
reg_pdpagl <- gls(pdform, data=data_pd[-c(13,14),],
                  correlation=corPagel(value=1, phy=tree, form = ~sp_pd))
intervals(reg_pdpagl, which = "var-cov") #optimal value=???

#Restrict the range to 0-1
lambda <- seq(0, 1, length.out = 500)
lik <- sapply(lambda,
              function(lambda) logLik(gls(pdform, data=data_pd[-c(13,14),], method = "ML",
                                          correlation=corPagel(value = lambda, phy = tree,
                                                               fixed = TRUE, form = ~sp_pd))))
plot(lik~lambda, type = "l", ylab = "Log-Likelihood", xlab = expression(lambda))
abline(v = lambda[which.max(lik)], col = 3, lwd = 2, lty = 2)
points(lambda[which.max(lik)], max(lik), col = 3, pch = 16)

#Refit the best model
reg_pdpagl_best <- gls(pdform, data=data_pd[-c(13,14),], method = "ML",
                   correlation=corPagel(value=lambda[which.max(lik)], phy=tree, form = ~sp_pd, fixed=T))
summary(reg_pdpagl_best)
reg_pdpagl_best$logLik

reg_pdpagl0 <- gls(pdform, data=data_pd[-c(13,14),], method = "ML",
                       correlation=corPagel(value=0, phy=tree, form = ~sp_pd, fixed=T))
summary(reg_pdpagl0)
reg_pdpagl0$logLik

#Refit a model with lambda=0.1, lambda=1
reg_pdpagl_01 <- gls(pdform, data=data_pd[-c(13,14),], method = "ML",
                       correlation=corPagel(value=0.1, phy=tree, form = ~sp_pd, fixed=T))
summary(reg_pdpagl_01)
reg_pdpagl_01$logLik

reg_pdpagl_1 <- gls(pdform, data=data_pd[-c(13,14),], method = "ML",
                     correlation=corPagel(value=1, phy=tree, form = ~sp_pd, fixed=T))
summary(reg_pdpagl_1)
reg_pdpagl_1$logLik

BIC(reg_pdpagl_01)
BIC(reg_pdpagl_1)
aictab(list(reg_pdpagl_best, reg_pdpagl_01, reg_pdpagl_1),
       modnames = c("Mass+BMRr, lambda=0.0", "Mass+BMRr, lambda=0.1", "Mass+BMRr, lambda=1"),
       second.ord = F)
