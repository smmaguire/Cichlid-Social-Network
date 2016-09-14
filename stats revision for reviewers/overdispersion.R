#checking for overdispersion
library(MASS)
library(vcd)
data(quine)
rando<-rpois(100000,25)
fit<-goodfit(quine$Days,type="poisson")
fit<-goodfit(rando)
summary(fit)
rootogram(fit)
glmer(outdegree_from_F~day*manip+(1|unique)+(1|community),family=poisson,data=ind.data)

data("HorseKicks")
HK.fit <- goodfit(HorseKicks)
summary(HK.fit)
plot(HK.fit)

Ord_plot(ind.data$outdegree_from_F)
Ord_plot(quine$Days)
Ord_plot(rando,type="poisson")

distplot(quine$Days,type="poisson")
distplot(ind.data$outdegree_from_F,type="poisson")
distplot(rando,type="poisson")

#install.packages('AER')
library(AER)
mod1<-glm(Days~Age+Sex, data=quine, family="poisson")
deviance(mod1)/mod1$df.residual
mod.r<-glmer(outdegree_from_F~day*manip+(1|unique)+(1|community)+(1|overd),family=poisson,data=ind.data)
mod.r2<-glmer(outdegree_from_F~day*manip+(1|unique)+(1|community),family=poisson,data=ind.data)
#mod.r<-glm(outdegree_from_F~day*manip,family=poisson,data=ind.data)

mod.glm<-glm(outdegree_from_M~day*manip,family="poisson",data=ind.data)
mod.zfl<-zeroinfl(outdegree_from_M~day*manip,dist="poisson",data=ind.data)
summary(mod.zfl)
AIC(mod.glm,mod.zfl)

summary(mod.r)
deviance(mod.r)/
  mod.rando$df.residual


767.8/161
mod.rando<-glm(rando~1,family="poisson")
deviance(mod.rando)/mod.rando$df.residual

res<-residuals(mod.r,type='deviance')
plot(log(predict(mod.r)+1),res)
abline(h=0,lty=2)
qqnorm(res)
qqline(res)

res<-residuals(mod.r2,type='deviance')
plot(log(predict(mod.r2)+1),res)
abline(h=0,lty=2)
qqnorm(res)
qqline(res)

res<-residuals(mod.glm,type='deviance')
plot(log(predict(mod.glm)+1),res)
abline(h=0,lty=2)
qqnorm(res)
qqline(res)

res<-residuals(mod.zfl,type='pearson')
plot(log(predict(mod.zfl)+1),res)
abline(h=0,lty=2)
qqnorm(res)
qqline(res)

install.packages('pscl')
qp.mod1<-glmmPQL(indegree_from_F~day+manip+day:manip,random =list(~ 1|unique, ~ 1|community), family =quasipoisson(link = "log"),data=ind.data)
qp.mod2<-glmmPQL(indegree_from_F~day+manip+day:manip,random =list(~ 1|community,~ 1|unique), family =negative.binomial(theta=5,link = "log"),data=ind.data)

res<-residuals(qp.mod,type='normalized')
plot(log(predict(qp.mod)),res)
abline(h=0,lty=2)
qqnorm(res)
qqline(res)
