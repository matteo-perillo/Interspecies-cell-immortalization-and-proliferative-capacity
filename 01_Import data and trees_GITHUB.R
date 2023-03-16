#BEFORE STARTING, MAKE SURE THAT THE CURRENT WORKING DIRECTORY
#CONTAINS THE FILE "TABLE FOR ANALYSIS.XLSX" AND THE FOLDER "TREES"

library(readxl)
library(ppcor)
library(ape)
library(phytools)
library(phangorn)

f = function(in_s){
  out_s = gsub(" ", "_", in_s)
  return(out_s)
}

#Import data and reset column names
data <- read_excel("Table for analysis.xlsx", sheet = "Table for analysis")
Species <- sapply(data$Species, f)
data <- data[-1]
row.names(data) <- Species

#Create the variable "Residual BMR"
cor(cbind(log(data$mass), log(data$BMR))[!is.na(log(data$BMR)),])[1,2] #0.99 - almost perfectly correlated
reg_massbmr = lm(log(mass)~log(BMR), data = data)
summary(reg_massbmr) #linear relation
plot(log(data[,c(3,1)]), xlab = "ln(BMR)", ylab = "ln(mass)", pch=16)
text(x = log(data$BMR)+1, y = log(data$mass), labels = rownames(data)[!is.na(data$BMR)], cex=0.5)
abline(reg_massbmr, col=2, lwd=2)
BMRr = rep(NA, length(data$BMR))
j = 1
for(i in 1:length(data$BMR)){
  if(!is.na(data$BMR[i])){
    BMRr[i] <- reg_massbmr$residuals[j]
    j = j+1
  }
}

data$BMRr <- BMRr
row.names(data) <- Species

plot(log(data$mass), data$BMRr, pch=16)
cor(na.omit(cbind(log(data$mass), data$BMRr)))

#Scatterplot
pairs(cbind(log(data[,c(1,2,4)]), car::logit(data$IF), data$BMRr),
      labels = c("log(mass)", "log(long)", "log(PD)", "logit(IP)", "BMRr"), pch = 16)

#Pairwise raw correlation mass-long
cor(log(data$mass), log(data$long)) #mass-long=0.54
#For Table 2
reg_ml <- lm(log(mass)~log(long), data = data)
summary(reg_ml)

#Create sub-dataset for PD analysis
data_pd = data[!is.na(data$PD),-c(3,5)]
row.names(data_pd) = row.names(data)[!is.na(data$PD)]

cor(log(data_pd[-4]))[3,-3] #PD-mass=0.83; PD-long=0.68
cor(cbind(log(na.omit(data_pd)[,1:3]), na.omit(data_pd)[,4]))[4,-4]
#BMRr vs: mass=0.15; long=0.60, PD=0.57
#For Table 2
reg_pl <- lm(log(PD)~log(long), data = data_pd)
summary(reg_pl)
reg_pm <- lm(log(PD)~log(mass), data = data_pd)
summary(reg_pm)
reg_mb <- lm(log(mass)~BMRr, data = data_pd)
summary(reg_mb)
reg_lb <- lm(log(long)~BMRr, data = data_pd)
summary(reg_lb)
reg_pb <- lm(log(PD)~BMRr, data = data_pd)
summary(reg_pb)

#Mass-longevity-PD (log) partial correlations
partcor_pd = pcor(log(data_pd[,-4]))
partcor_pd$estimate
partcor_pd$p.value #only mass-PD correlated when controlling for the third variable

#Create sub-dataset for IF analysis
data_if = data[!is.na(data$IF),-c(3,4)]
row.names(data_if) = row.names(data)[!is.na(data$IF)]

#Mass-longevity-IF (log and logit) partial correlations
partcor_if = pcor(cbind(log(data_if[1:2]), car::logit(data_if$IF)))
partcor_if$estimate
partcor_if$p.value #only mass-IF correlated when controlling for the third variable

#Import trees
filename <- "trees.nex"
trees <- ape::read.nexus(filename) #trees from VertLife

#Average tree
rf.tree <- averageTree(trees,method="symmetric.difference")
#Consensus tree
c.tree <- consensus(trees,p=0.95)
#Median tree
ss.d <- colSums(as.matrix(RF.dist(trees))^2)
ii <- which(ss.d==min(ss.d))
set.seed(10)
i = sample(ii,1)

#Plots
plotTree(rf.tree)
plotTree(root(rf.tree,outgroup=c("Heterocephalus_glaber","Sciurus_carolinensis"), resolve.root=TRUE))
plotTree(c.tree)
plotTree(root(c.tree,outgroup=c("Heterocephalus_glaber","Sciurus_carolinensis"), resolve.root=TRUE))
plotTree(trees[[i]])
plotTree(compute.brlen(root(trees[[i]],outgroup=c("Heterocephalus_glaber","Sciurus_carolinensis"), resolve.root=TRUE)))

tree=trees[[100]]
plot(tree, direction="right", srt=0, label.offset = 2, cex = 1)

wd = getwd()
save.image(paste(wd, "Imported data and trees.RData", sep = "/"))
