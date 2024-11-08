rm(list=ls())
graphics.off()

library(vegan)
library(gam)
library(scam)
library(gdm)
library(car)
library(betapart)
library(tidyverse)
library(zetadiv)

quartz()

load("islands_Soc_Haw_null_models_practical.RData")



################
##Permutations##
################

##Society
beta.soc.obs <- c(mean(beta.pair(dat.soc.pa)$beta.sor),mean(beta.pair(dat.soc.pa)$beta.sim))

Sys.time()
dat.soc.pa.perm.none <- permatfull(m=dat.soc.pa,fixedmar = "none",times=99,mtype="prab")
Sys.time()
dat.soc.pa.perm.row <- permatfull(m=dat.soc.pa,fixedmar = "rows",times=99,mtype="prab")
Sys.time()
dat.soc.pa.perm.col <- permatfull(m=dat.soc.pa,fixedmar = "columns",times=99,mtype="prab")
Sys.time()
dat.soc.pa.perm.both <- permatfull(m=dat.soc.pa,fixedmar = "both",times=99,mtype="prab")
Sys.time()


beta.mean <- function(dat){
  beta.dat <- beta.pair(dat)
  return(c(mean(beta.dat$beta.sor),mean(beta.dat$beta.sim)))
}

beta.rand.soc.none <- data.frame(matrix(unlist(lapply(dat.soc.pa.perm.none$perm,beta.mean)),99,2,byrow = TRUE))
names(beta.rand.soc.none) <- c("Sorensen","Simpson")

beta.rand.soc.row <- data.frame(matrix(unlist(lapply(dat.soc.pa.perm.row$perm,beta.mean)),99,2,byrow = TRUE))
names(beta.rand.soc.row) <- c("Sorensen","Simpson")

beta.rand.soc.col <- data.frame(matrix(unlist(lapply(dat.soc.pa.perm.col$perm,beta.mean)),99,2,byrow = TRUE))
names(beta.rand.soc.col) <- c("Sorensen","Simpson")

beta.rand.soc.both <- data.frame(matrix(unlist(lapply(dat.soc.pa.perm.both$perm,beta.mean)),99,2,byrow = TRUE))
names(beta.rand.soc.both) <- c("Sorensen","Simpson")


par(mfrow=c(2,4))
hist(beta.rand.soc.none$Sorensen,breaks=seq(0.5,1,0.01),main="None fixed",xlab="Sorensen")
abline(v=quantile(beta.rand.soc.none$Sorensen,0.025),col="blue")
abline(v=quantile(beta.rand.soc.none$Sorensen,0.0975),col="blue")
abline(v=beta.soc.obs[1],col="red")
hist(beta.rand.soc.row$Sorensen,breaks=seq(0.5,1,0.01),main="Rows fixed",xlab="Sorensen")
abline(v=quantile(beta.rand.soc.row$Sorensen,0.025),col="blue")
abline(v=quantile(beta.rand.soc.row$Sorensen,0.0975),col="blue")
abline(v=beta.soc.obs[1],col="red")
hist(beta.rand.soc.col$Sorensen,breaks=seq(0.5,1,0.01),main="Columns fixed",xlab="Sorensen")
abline(v=quantile(beta.rand.soc.col$Sorensen,0.025),col="blue")
abline(v=quantile(beta.rand.soc.col$Sorensen,0.0975),col="blue")
abline(v=beta.soc.obs[1],col="red")
hist(beta.rand.soc.both$Sorensen,breaks=seq(0.5,1,0.01),main="Both fixed",xlab="Sorensen")
abline(v=quantile(beta.rand.soc.both$Sorensen,0.025),col="blue")
abline(v=quantile(beta.rand.soc.both$Sorensen,0.0975),col="blue")
abline(v=beta.soc.obs[1],col="red")

hist(beta.rand.soc.none$Simpson,breaks=seq(0.2,1,0.01),main="",xlab="Simpson")
abline(v=quantile(beta.rand.soc.none$Simpson,0.025),col="blue")
abline(v=quantile(beta.rand.soc.none$Simpson,0.0975),col="blue")
abline(v=beta.soc.obs[2],col="red")
hist(beta.rand.soc.row$Simpson,breaks=seq(0.2,1,0.01),main="",xlab="Simpson")
abline(v=quantile(beta.rand.soc.row$Simpson,0.025),col="blue")
abline(v=quantile(beta.rand.soc.row$Simpson,0.0975),col="blue")
abline(v=beta.soc.obs[2],col="red")
hist(beta.rand.soc.col$Simpson,breaks=seq(0.2,1,0.01),main="",xlab="Simpson")
abline(v=quantile(beta.rand.soc.col$Simpson,0.025),col="blue")
abline(v=quantile(beta.rand.soc.col$Simpson,0.0975),col="blue")
abline(v=beta.soc.obs[2],col="red")
hist(beta.rand.soc.both$Simpson,breaks=seq(0.2,1,0.01),main="",xlab="Simpson")
abline(v=quantile(beta.rand.soc.both$Simpson,0.025),col="blue")
abline(v=quantile(beta.rand.soc.both$Simpson,0.0975),col="blue")
abline(v=beta.soc.obs[2],col="red")



##Hawai'i
beta.haw.obs <- c(mean(beta.pair(dat.haw.pa)$beta.sor),mean(beta.pair(dat.haw.pa)$beta.sim))

# Sys.time()
# dat.haw.pa.perm.none <- permatfull(m=dat.haw.pa,fixedmar = "none",times=99,mtype="prab")
# Sys.time()
# dat.haw.pa.perm.row <- permatfull(m=dat.haw.pa,fixedmar = "rows",times=99,mtype="prab")
# Sys.time()
# dat.haw.pa.perm.col <- permatfull(m=dat.haw.pa,fixedmar = "columns",times=99,mtype="prab")
# Sys.time()
# dat.haw.pa.perm.both <- permatfull(m=dat.haw.pa,fixedmar = "both",times=99,mtype="prab")
# Sys.time()
# 
# save(dat.haw.pa.perm.none,dat.haw.pa.perm.row,dat.haw.pa.perm.col,dat.haw.pa.perm.both,file="permutations_hawaii.RData")
load("permutations_hawaii.RData")

beta.rand.haw.none <- data.frame(matrix(unlist(lapply(dat.haw.pa.perm.none$perm,beta.mean)),99,2,byrow = TRUE))
names(beta.rand.haw.none) <- c("Sorensen","Simpson")

beta.rand.haw.row <- data.frame(matrix(unlist(lapply(dat.haw.pa.perm.row$perm,beta.mean)),99,2,byrow = TRUE))
names(beta.rand.haw.row) <- c("Sorensen","Simpson")

beta.rand.haw.col <- data.frame(matrix(unlist(lapply(dat.haw.pa.perm.col$perm,beta.mean)),99,2,byrow = TRUE))
names(beta.rand.haw.col) <- c("Sorensen","Simpson")

beta.rand.haw.both <- data.frame(matrix(unlist(lapply(dat.haw.pa.perm.both$perm,beta.mean)),99,2,byrow = TRUE))
names(beta.rand.haw.both) <- c("Sorensen","Simpson")



par(mfrow=c(2,4))
hist(beta.rand.haw.none$Sorensen,breaks=seq(0,1,0.025),main="None fixed",xlab="Sorensen")
abline(v=quantile(beta.rand.haw.none$Sorensen,0.025),col="blue")
abline(v=quantile(beta.rand.haw.none$Sorensen,0.0975),col="blue")
abline(v=beta.haw.obs[1],col="red")
hist(beta.rand.haw.row$Sorensen,breaks=seq(0,1,0.025),main="Rows fixed",xlab="Sorensen")
abline(v=quantile(beta.rand.haw.row$Sorensen,0.025),col="blue")
abline(v=quantile(beta.rand.haw.row$Sorensen,0.0975),col="blue")
abline(v=beta.haw.obs[1],col="red")
hist(beta.rand.haw.col$Sorensen,breaks=seq(0,1,0.025),main="Columns fixed",xlab="Sorensen")
abline(v=quantile(beta.rand.haw.col$Sorensen,0.025),col="blue")
abline(v=quantile(beta.rand.haw.col$Sorensen,0.0975),col="blue")
abline(v=beta.haw.obs[1],col="red")
hist(beta.rand.haw.both$Sorensen,breaks=seq(0,1,0.025),main="Both fixed",xlab="Sorensen")
abline(v=quantile(beta.rand.haw.both$Sorensen,0.025),col="blue")
abline(v=quantile(beta.rand.haw.both$Sorensen,0.0975),col="blue")
abline(v=beta.haw.obs[1],col="red")

hist(beta.rand.haw.none$Simpson,breaks=seq(0,1,0.025),main="",xlab="Simpson")
abline(v=quantile(beta.rand.haw.none$Simpson,0.025),col="blue")
abline(v=quantile(beta.rand.haw.none$Simpson,0.0975),col="blue")
abline(v=beta.haw.obs[2],col="red")
hist(beta.rand.haw.row$Simpson,breaks=seq(0,1,0.025),main="",xlab="Simpson")
abline(v=quantile(beta.rand.haw.row$Simpson,0.025),col="blue")
abline(v=quantile(beta.rand.haw.row$Simpson,0.0975),col="blue")
abline(v=beta.haw.obs[2],col="red")
hist(beta.rand.haw.col$Simpson,breaks=seq(0,1,0.025),main="",xlab="Simpson")
abline(v=quantile(beta.rand.haw.col$Simpson,0.025),col="blue")
abline(v=quantile(beta.rand.haw.col$Simpson,0.0975),col="blue")
abline(v=beta.haw.obs[2],col="red")
hist(beta.rand.haw.both$Simpson,breaks=seq(0,1,0.025),main="",xlab="Simpson")
abline(v=quantile(beta.rand.haw.both$Simpson,0.025),col="blue")
abline(v=quantile(beta.rand.haw.both$Simpson,0.0975),col="blue")
abline(v=beta.haw.obs[2],col="red")




#####################
##Richness analyses##
#####################

load("islands_FP_regression_models_practical.RData")

##alpha diversity
islands.pred.FP$alpha <- rowSums(dat.FP.pa)

cor(islands.pred.FP[,2:7])

##GLMs
mod.glm.FP <- glm(alpha~elev_max+IslandArea+temp.mean+prec.an+temp.seas+prec.seas,data=islands.pred.FP,family = poisson())
mod.glm.FPb <- glm(alpha~elev_max+IslandArea+temp.seas+prec.an+prec.seas,data=islands.pred.FP,family = poisson())
mod.glm.FPc <- glm(alpha~elev_max+IslandArea+temp.seas+prec.an,data=islands.pred.FP,family = poisson())
mod.glm.FPd <- glm(alpha~elev_max+IslandArea,data=islands.pred.FP,family = poisson())

vif(mod.glm.FP)
vif(mod.glm.FPb)

AIC(mod.glm.FPb)
AIC(mod.glm.FPc)
AIC(mod.glm.FPd)


##GAMs
mod.gam.FP <- gam(alpha~s(elev_max)+s(IslandArea)+s(temp.seas)+s(prec.an)+s(prec.an)+s(prec.seas),data=islands.pred.FP,family = poisson(),method = "REML")

##Smoother GAMs
mod.gam.FP2 <- gam(alpha~s(elev_max,k=4)+s(IslandArea,k=4)+s(temp.seas,k=4)+s(prec.an,k=4)+s(prec.an,k=4)+s(prec.seas,k=4),data=islands.pred.FP,family = poisson(),method = "REML")

##Model summary
summary(mod.glm.FP)
summary(mod.glm.FPb)
summary(mod.glm.FPc)
summary(mod.gam.FP)
summary(mod.gam.FP2)

##model performance
cor(islands.pred.FP$alpha,predict(mod.glm.FP,type="response"))
cor(islands.pred.FP$alpha,predict(mod.glm.FPb,type="response"))
cor(islands.pred.FP$alpha,predict(mod.glm.FPc,type="response"))
cor(islands.pred.FP$alpha,predict(mod.gam.FP,type="response"))
cor(islands.pred.FP$alpha,predict(mod.gam.FP2,type="response"))

par(mfrow=c(1,5))
plot(islands.pred.FP$alpha,predict(mod.glm.FP,type="response"),type="p")
plot(islands.pred.FP$alpha,predict(mod.glm.FPb,type="response"),type="p")
plot(islands.pred.FP$alpha,predict(mod.glm.FPc,type="response"),type="p")
plot(islands.pred.FP$alpha,predict(mod.gam.FP,type="response"),type="p")
plot(islands.pred.FP$alpha,predict(mod.gam.FP2,type="response"),type="p")

##Plot GAM splines
par(mfrow=c(2,5))
plot(mod.gam.FP,ylim=c(-20,20),residuals = T,pch=1)
plot(mod.gam.FP2,ylim=c(-20,20),residuals = T,pch=1)



######
###Same as above, log transforming elevation and area
#######

##GLMs
mod.glm.FP <- glm(alpha~log(elev_max)+log(IslandArea)+temp.seas+prec.an+prec.seas,data=islands.pred.FP,family = poisson())

##GAMs
mod.gam.FP <- gam(alpha~s(log(elev_max))+s(log(IslandArea))+s(temp.seas)+s(prec.an)+s(prec.seas),data=islands.pred.FP,family = poisson(),method = "REML")

##Smoother GAMs
mod.gam.FP2 <- gam(alpha~s(log(elev_max),k=4)+s(log(IslandArea),k=4)+s(temp.seas,k=4)+s(prec.an,k=4)+s(prec.seas,k=4),data=islands.pred.FP,family = poisson(),method = "REML")

##Model summary
summary(mod.glm.FP)
summary(mod.gam.FP)
summary(mod.gam.FP2)

##model performance
cor(islands.pred.FP$alpha,predict(mod.glm.FP,type="response"))
cor(islands.pred.FP$alpha,predict(mod.gam.FP,type="response"))
cor(islands.pred.FP$alpha,predict(mod.gam.FP2,type="response"))

par(mfrow=c(1,3))
plot(islands.pred.FP$alpha,predict(mod.glm.FP,type="response"),type="p")
plot(islands.pred.FP$alpha,predict(mod.gam.FP,type="response"),type="p")
plot(islands.pred.FP$alpha,predict(mod.gam.FP2,type="response"),type="p")

##Plot GAM splines
par(mfrow=c(2,5))
plot(mod.gam.FP,ylim=c(-20,20),residuals = T,pch=1)
plot(mod.gam.FP2,ylim=c(-20,20),residuals = T,pch=1)




#####################
##Turnover analyses##
#####################

load("island_coordinates_FP.RData") ##load coordinates

##GDMs
GDM.FP.Sor <- Zeta.msgdm(data.spec = dat.FP.pa, data.env = islands.pred.FP[c(2,3,5,6,7)], xy = FP.coords[2:3],order=2,reg.type = "ispline",method.glm = "glm.fit.cons",cons = -1, cons.inter = -1,distance.type="ortho",normalize="Sorensen",family = binomial("log"))
GDM.FP.Sim <- Zeta.msgdm(data.spec = dat.FP.pa, data.env = islands.pred.FP[c(2,3,5,6,7)], xy = FP.coords[2:3],order=2,reg.type = "ispline",method.glm = "glm.fit.cons",cons = -1, cons.inter = -1,distance.type="ortho",normalize="Simpson",family = binomial("log"))
# GDM.FP.Sor <- Zeta.msgdm(data.spec = dat.FP.pa, data.env = islands.pred.FP[c(2,3,5,6,7)],order=2,reg.type = "ispline",method.glm = "glm.fit.cons",cons = -1, cons.inter = -1,distance.type="ortho",normalize="Sorensen",family = binomial("log"))
# GDM.FP.Sim <- Zeta.msgdm(data.spec = dat.FP.pa, data.env = islands.pred.FP[c(2,3,5,6,7)],order=2,reg.type = "ispline",method.glm = "glm.fit.cons",cons = -1, cons.inter = -1,distance.type="ortho",normalize="Simpson",family = binomial("log"))

##Model performance
Dat.FP.Sor.pred <- Predict.msgdm(model.msgdm = GDM.FP.Sor$model, reg.type = "ispline", newdata = GDM.FP.Sor$predictors)
cor(GDM.FP.Sor$val,Dat.FP.Sor.pred)^2
Dat.FP.Sim.pred <- Predict.msgdm(model.msgdm = GDM.FP.Sor$model, reg.type = "ispline", newdata = GDM.FP.Sor$predictors)
cor(GDM.FP.Sim$val,Dat.FP.Sim.pred)^2

##Plot outputs
par(mfrow=c(2,2))
plot(GDM.FP.Sor$val,Dat.FP.Sor.pred,pch=20,col="blue")
lines(sort(GDM.FP.Sor$val),predict(lm(Dat.FP.Sor.pred~GDM.FP.Sor$val))[order(GDM.FP.Sor$val)])
plot(GDM.FP.Sim$val,Dat.FP.Sim.pred,pch=20,col="blue")
lines(sort(GDM.FP.Sim$val),predict(lm(Dat.FP.Sim.pred~GDM.FP.Sim$val))[order(GDM.FP.Sim$val)])
Plot.ispline(msgdm=GDM.FP.Sor,data.env = islands.pred.FP[c(2,3,5,6,7)],distance=TRUE)
Plot.ispline(msgdm=GDM.FP.Sim,data.env = islands.pred.FP[c(2,3,5,6,7)],distance=TRUE)
# Plot.ispline(msgdm=GDM.FP.Sor,data.env = islands.pred.FP[c(2,3,5,6,7)],distance=FALSE)
# Plot.ispline(msgdm=GDM.FP.Sim,data.env = islands.pred.FP[c(2,3,5,6,7)],distance=FALSE)


##log transform elevation and area

islands.pred.FP2 <- islands.pred.FP
islands.pred.FP2$elev_max <- log(islands.pred.FP2$elev_max)
islands.pred.FP2$IslandArea <- log(islands.pred.FP2$IslandArea)

GDM.FP.Sor2 <- Zeta.msgdm(data.spec = dat.FP.pa, data.env = islands.pred.FP2[c(2,3,5,6,7)], xy = FP.coords[2:3],order=2,reg.type = "ispline",method.glm = "glm.fit.cons",cons = -1, cons.inter = -1,distance.type="ortho",normalize="Sorensen",family = binomial("log"))
GDM.FP.Sim2 <- Zeta.msgdm(data.spec = dat.FP.pa, data.env = islands.pred.FP2[c(2,3,5,6,7)], xy = FP.coords[2:3],order=2,reg.type = "ispline",method.glm = "glm.fit.cons",cons = -1, cons.inter = 91,distance.type="ortho",normalize="Simpson",family = binomial("log"))
# GDM.FP.Sor2 <- Zeta.msgdm(data.spec = dat.FP.pa, data.env = islands.pred.FP2[c(2,3,5,6,7)], order=2,reg.type = "ispline",method.glm = "glm.fit.cons",cons = -1, cons.inter = -1,distance.type="ortho",normalize="Sorensen",family = binomial("log"))
# GDM.FP.Sim2 <- Zeta.msgdm(data.spec = dat.FP.pa, data.env = islands.pred.FP2[c(2,3,5,6,7)], order=2,reg.type = "ispline",method.glm = "glm.fit.cons",cons = -1, cons.inter = 91,distance.type="ortho",normalize="Simpson",family = binomial("log"))

Dat.FP.Sor.pred2 <- Predict.msgdm(model.msgdm = GDM.FP.Sor2$model, reg.type = "ispline", newdata = GDM.FP.Sor2$predictors)
cor(GDM.FP.Sor2$val,Dat.FP.Sor.pred2)^2
Dat.FP.Sim.pred2 <- Predict.msgdm(model.msgdm = GDM.FP.Sor2$model, reg.type = "ispline", newdata = GDM.FP.Sor2$predictors)
cor(GDM.FP.Sim2$val,Dat.FP.Sim.pred2)^2


par(mfrow=c(2,2))
plot(GDM.FP.Sor2$val,Dat.FP.Sor.pred2,pch=20,col="blue")
lines(sort(GDM.FP.Sor2$val),predict(lm(Dat.FP.Sor.pred2~GDM.FP.Sor2$val))[order(GDM.FP.Sor2$val)])
plot(GDM.FP.Sim2$val,Dat.FP.Sim.pred2,pch=20,col="blue")
lines(sort(GDM.FP.Sim2$val),predict(lm(Dat.FP.Sim.pred2~GDM.FP.Sim2$val))[order(GDM.FP.Sim2$val)])
Plot.ispline(msgdm=GDM.FP.Sor2,data.env = islands.pred.FP2[c(2,3,5,6,7)],distance=TRUE)
Plot.ispline(msgdm=GDM.FP.Sim2,data.env = islands.pred.FP2[c(2,3,5,6,7)],distance=TRUE)
# Plot.ispline(msgdm=GDM.FP.Sor2,data.env = islands.pred.FP2[c(2,3,5,6,7)],distance=FALSE)
# Plot.ispline(msgdm=GDM.FP.Sim2,data.env = islands.pred.FP2[c(2,3,5,6,7)],distance=FALSE)





####
##using the gdm package
####

dat.FP.pa2 <- cbind(row.names(dat.FP.pa),dat.FP.pa)
names(dat.FP.pa2)[1] <- "island"
dat.FP.pa.long <- dat.FP.pa2 %>% 
  pivot_longer(cols=!island,
               names_to = "Species",
               values_to = "Presence")

dat.FP.pa.long <- dat.FP.pa.long[-which(dat.FP.pa.long$Presence==0),1:2]

dat.FP.pa.long$IslandArea <- numeric(nrow(dat.FP.pa.long))
dat.FP.pa.long$elev_max <- numeric(nrow(dat.FP.pa.long))
dat.FP.pa.long$temp.mean <- numeric(nrow(dat.FP.pa.long))
dat.FP.pa.long$temp.seas <- numeric(nrow(dat.FP.pa.long))
dat.FP.pa.long$precan <- numeric(nrow(dat.FP.pa.long))
dat.FP.pa.long$prec.seas <- numeric(nrow(dat.FP.pa.long))
dat.FP.pa.long$long <- numeric(nrow(dat.FP.pa.long))
dat.FP.pa.long$lat <- numeric(nrow(dat.FP.pa.long))

for(i in row.names(dat.FP.pa)){
  toto <- islands.pred.FP[which(islands.pred.FP$Name_USGSO==i),2:7]
  dat.FP.pa.long[which(dat.FP.pa.long$island==i),3:8] <- toto[1:6]
  toto <- FP.coords[which(FP.coords$island==i),2:3]
  dat.FP.pa.long[which(dat.FP.pa.long$island==i),9:10] <- toto
}

dat.FP.pa.long <- as.data.frame(dat.FP.pa.long)

head(dat.FP.pa.long)

sppData.FP <- dat.FP.pa.long[c(1,2,9,10)] 
envTab.FP <- dat.FP.pa.long[c(1,3:10)] 
#envTab.FP <- dat.FP.pa.long[c(1,4)] 
sitePairTab.FP <- formatsitepair(sppData.FP,2,XColumn="long",YColumn="lat",sppColumn="Species", siteColumn="island",predData=envTab.FP) 
gdmTabMod.FP<-gdm(sitePairTab.FP,geo=TRUE)
summary(gdmTabMod.FP)
plot(gdmTabMod.FP)






