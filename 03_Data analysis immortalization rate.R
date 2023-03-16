#BEFORE STARTING, RUN THE SCRIPT "01_IMPORT DATA AND TREES", OR LOAD THE
#WORKING ENVIRONMENT "IMPORTED DATA AND TREES"

library(ape)
library(phytools)
library(phangorn)
library(ppcor)
library(stargazer)
#stargazer(reg_cp_peso, type = "text")
library(ape)
library(car)
library(phytools)
library(geiger)
library(nlme)
library(flexmix)
library(AICcmodavg)

fitted_values = function(model, xvar){
  b0 = model$coefficients[1]
  b1 = model$coefficients[2]
  x = seq(min(xvar)-0.5, max(xvar)+0.5, 0.02)
  y_hat = exp(b0+b1*x)/(1+exp(b0+b1*x))
  return(data.frame(x, y_hat))
}
#Well done!

#Scatterplot
pairs(cbind(log(data_if[,1:2]), car::logit(data_if$IF), data_if$BMRr),
      labels = c("log(mass)", "log(long)", "logit(IF)", "BMRr"), pch = 16)

#Raw and partial correlations considering also BMRr
rawcor <- cor(cbind(log(data_if[,1:2]), car::logit(data_if[,3]), data_if[,4])[!is.na(data_if$BMRr),])
rawcor
partcor = pcor(cbind(log(data_if[,1:2]), car::logit(data_if[,3]), data_if[,4])[!is.na(data_if$BMRr),])
partcor$estimate
partcor$p.value #the only statistically significant correlation is mass-IF

#Null model
reg_if0_glm <- glm(IF~1, data = data_if, family = quasibinomial)
reg_if0_glm_2 <- glm(IF~1, data = na.omit(data_if), family = quasibinomial)

#Simple linear regressions
#Mass vs IF
reg_ifm_glm <- glm(IF~log(mass), data = data_if, family = quasibinomial)
summary(reg_ifm_glm) #mass significant (pval=0.004)
fit_ifm_glm = fitted_values(model = reg_ifm_glm, xvar = log(data_if$mass))
plot(cbind(log(data_if[,1]), data_if[,3]),
     xlab = "ln(mass)", ylab = "Immortalization Frequency", pch=16)
lines(fit_ifm_glm, col = 2, lwd = 2)
#plot(reg_ifm_glm, pch=16) #mus musculus outlier

reg_ifm_glm_2 <- glm(IF~log(mass), data = na.omit(data_if), family = quasibinomial)
summary(reg_ifm_glm_2)

#Longevity vs IF
reg_ifl_glm <- glm(IF~log(long), data = data_if, family = quasibinomial)
summary(reg_ifl_glm) #long not significant (pval=0.068)
fit_ifl_glm = fitted_values(model = reg_ifl_glm, xvar = log(data_if$long))
plot(cbind(log(data_if[,2]), data_if[,3]),
     xlab = "ln(long)", ylab = "IF", pch=16)
lines(fit_ifl_glm, col = 2, lwd = 2)
#plot(reg_ifl_glm)

#BMRr vs IF
reg_ifb_glm <- glm(IF~BMRr, data = data_if, family = quasibinomial)
summary(reg_ifb_glm) #BMRr not significant (pval=0.565)
fit_ifb_glm = fitted_values(model = reg_ifb_glm, xvar = data_if[!is.na(data_if$BMRr),4])
plot(cbind(data_if$BMRr, data_if[,3]),
     xlab = "BMRr", ylab = "IF", pch=16)
lines(fit_ifb_glm, col = 2, lwd = 2)
#plot(reg_ifb_glm)

#Multiple linear regression
#Body mass and Longevity vs IF
reg_ifmalo_glm <- glm(IF~log(long)+log(mass), data = data_if, family=quasibinomial)
summary(reg_ifmalo_glm) #mass significant (0.021), long not significant (0.725)

reg_ifmalo_glm_2 <- glm(IF~log(long)+log(mass), data = na.omit(data_if), family=quasibinomial)
summary(reg_ifmalo_glm_2)

#Body mass and BMRr vs IF
reg_ifmabr_glm <- glm(IF~log(mass)+BMRr, data = data_if, family=quasibinomial)
summary(reg_ifmabr_glm) #mass significant, BMRr not

#Body mass, Longevity and BMRr vs IF
reg_iffull_glm <- glm(IF~log(mass)+log(long)+BMRr, data = data_if, family=quasibinomial)
summary(reg_iffull_glm) #nothing significant (mass lowest pval)

#F-tests
anova(reg_ifm_glm, reg_if0_glm, test = "F")
anova(reg_ifmalo_glm, reg_if0_glm, test = "F")
anova(reg_ifm_glm_2, reg_if0_glm_2, test = "F")
anova(reg_ifmabr_glm, reg_if0_glm_2, test = "F")
anova(reg_ifmalo_glm_2, reg_if0_glm_2, test = "F")
anova(reg_iffull_glm, reg_if0_glm_2, test = "F")

#Analysis taking phylogenesis into account
ifm_form = car::logit(IF)~log(mass)

#Pagel lambda phylogenetic signal evolution model
sp_if <- row.names(data_if)
lambda <- seq(0, 1, length.out = 500)
lik <- sapply(lambda,
              function(lambda) logLik(gls(ifm_form, data = data_if, method = "ML",
                                          correlation = corPagel(value = lambda, phy = tree,
                                                                 fixed = TRUE, form = ~sp_if))))
plot(lik~lambda, type = "l", ylab = "Log Likelihood", xlab = expression(lambda))
abline(v = lambda[which.max(lik)], col = 3, lwd = 2)
points(c(lambda[which.max(lik)]), max(lik), col = 3, pch = 16)

#Best model (Pagel lambda = 0)
reg_ifpagl0 <- gls(ifm_form, data=data_if, method = "ML",
                   correlation=corPagel(value=0, phy=tree, form = ~sp_if, fixed=T))
summary(reg_ifpagl0)
#THE END
