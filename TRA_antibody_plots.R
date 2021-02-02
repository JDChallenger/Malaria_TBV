#This script generates Figure 4 (all panels), Figure 3D & 3E, Suppl. Figures 3, 5A & 5B

library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(gridExtra)
library(dplyr)
library(grid)

themeJDC <- theme(panel.background = element_blank(),panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black" ,fill = NA), legend.background = element_blank(),
                  legend.key = element_blank(), axis.text = element_text(size = 11.2),
                  axis.title = element_text(size=12), legend.text = element_text(size = 10.5),
                  legend.title = element_text(size = 10.8))
source("vaccine_params.R") #Parameters from the vaccine model

TRA <- function(titre){

  ret <- ((titre / mu25)^hill1)/(((titre / mu25)^hill1) + hill2)

  return(ret)
}

TBA <- function(titre){
  #TBA model (via SMFA, Bompard [2017])
  ret <- 0.917*(((titre/12.628)^1.09)/(((titre/12.628)^1.09) + 0.72))

  return(ret)
}

traTitre <- seq(1.5,40.5,0.5)
tbaTitre <- seq(1.5,40.5,0.5)

for(i in 1:length(traTitre)){

  traTitre[i] <- TRA((i/2))
  tbaTitre[i] <- TBA((i/2))

}
tbaTitre
ttr <- seq(1.5,40.5,0.5)
dftr = data.frame("titre" = ttr, "TRA" = traTitre, "TBA" = tbaTitre)
dftrm <- melt(dftr, id.vars = "titre")
pl1 <- ggplot(dftrm, aes(x = titre, y = value, linetype = variable)) + geom_line(size=0.9,color='darkblue') +
  themeJDC + ylim(c(0,1)) + theme(legend.position = c(0.75,0.3), legend.key.size = unit(0.95,"cm")) +
  xlim(c(0,31.5)) + xlab(expression(paste("Antibody Titre [",mu,"g / ml]"))) + ylab("Vaccine Activity") +
  scale_linetype_manual(values=c("solid", "dotted")) + labs(linetype="Lab-estimated activity")
#pl1

antibody <- function(time,ds,dl){

  ret <- tau25*(rho25*(exp(-time*log(2)/ds)) +
                  (1-rho25)*(exp(-time*log(2)/dl)))
}

AB <- seq(0,364,1)
ABb <- seq(0,364,1)
ABw <- seq(0,364,1)
tra <- seq(0,364,1)
traB <- seq(0,364,1)
traW <- seq(0,364,1)

for(i in seq(1,365,1)){
  AB[i] <- antibody(i-1,ds25,dl25)
  ABb[i] <- antibody(i-1,ds25better,dl25better)
  ABw[i] <- antibody(i-1,ds25worse,dl25worse)

  tra[i] <- TRA(AB[i])
  traB[i] <- TRA(ABb[i])
  traW[i] <- TRA(ABw[i])
}

tt <- seq(0,364,1)
abod <- data.frame(time = tt, "Intermediate titre" = AB, "Higher titre" = ABb, "Lower titre" = ABw)
#names(abod) <- c("time", "AB")
abodm <- melt(abod, id.vars = "time")
pl2 <- ggplot(abodm, aes(x=time, y=value, color = variable)) + geom_line(size=1) + themeJDC +
  ylab(expression(paste("Antibody Titre [ ",mu,"g / ml]"))) + scale_color_manual("Duration of antibody response",
    values=c("#1b9e77","#d95f02","#7570b3"),
    labels=c("RTS,S","Longer","Shorter")) + ylim(c(0.0,22.7)) + xlab("Time Post Vaccination (Days)") +
  theme(legend.position = c(0.7,0.8), plot.margin = unit(c(5.5,5.5,5.5,13), "pt"))
#pl2

traDF <- data.frame(time = tt, "Intermediate titre" = tra, "Higher titre" = traB, "Lower titre" = traW)
traDFm <- melt(traDF, id.vars = "time")
pl3 <- ggplot(traDFm, aes(x=time, y=value, color = variable)) + geom_line(size=1) + themeJDC +
  xlab("Time Post Vaccination (Days)") +
  ylab("Transmission Reducing Activity") +
  scale_color_manual("Duration of vaccine activity", values=c("#1b9e77","#d95f02","#7570b3"),
      labels=c("RTS,S","Longer","Shorter")) + ylim(c(0.0,1)) +
  theme(legend.position = c(0.28,0.18))

#Figure 2 in article
grid.arrange(pl1,pl2,pl3,nrow=1)

print(c(TRA(AB[183]), TRA(ABb[183]), TRA(ABw[183])))

TBAz <- function(r,m,TRAz){

  ret <- (1/(1-(r/(r+m))^r)) * ( (r/(r+m*(1-TRAz)))^r - (r/(r+m))^r)

  return(ret)

}
mT <- 30.3 #Negative binomial parameters
rrT <- 0.91 #Negative binomial parameters

mD <- 46.7 #Negative binomial parameters
rrD <- 0.91 #Negative binomial parameters

mA <- 3.55  #Negative binomial parameters
rrA <- 0.91 #Negative binomial parameters

mU <- 0.84 #Negative binomial parameters
rrU <- 0.91 #Negative binomial parameters

mALL <- 0.000157#0.4 #Negative binomial parameters (guess)
rrALL <- 0.00000495#0.01 #Negative binomial parameters (guess)

tbaD <- seq(0.01,1,0.01)
tbaT <- seq(0.01,1,0.01)
tbaA <- seq(0.01,1,0.01)
tbaU <- seq(0.01,1,0.01)
tbaALL <- seq(0.01,1,0.01) #new

for(i in seq(1,100,1)){
  tbaD[i] <- TBAz(rrD, mD, i*0.01)
  tbaT[i] <- TBAz(rrT, mT, i*0.01)
  tbaA[i] <- TBAz(rrA, mA, i*0.01)
  tbaU[i] <- TBAz(rrU, mU, i*0.01)
  tbaALL[i] <- TBAz(rrALL, mALL, i*0.01)
}

#eff_v <- 0.16*df4simp3[,1] + 0.23*df4simp3[,2] + 0.65*df4simp3[,3] + 0.82*df4simp3[,4]
#mean(TBAz(rrD, mD, 0.9)*df4simp3[,1] + TBAz(rrT, mT, 0.9)*df4simp3[,2] +
#  TBAz(rrA, 4.0, 0.9)*df4simp3[,3] + TBAz(rrU, mU, 0.9)*df4simp3[,4])

# klee <- rep(0,99) #Klesso
# for(i in 1:99){
#   klee[i] <- mean(TBAz(rrD, 50, 0.01*i)*df4simp3[,1] + TBAz(rrT, mT, 0.01*i)*df4simp3[,2] +
#                     TBAz(rrA, 4.3, 0.01*i)*df4simp3[,3] + TBAz(rrU, mU, 0.01*i)*df4simp3[,4])
# }
# klee
#
# loon <- rep(0,99) #Klesso
# for(i in 1:99){
#   loon[i] <- mean(TBAz(rrD, 50, 0.01*i)*df4simp3[,1] + TBAz(rrT, mT, 0.01*i)*df4simp3[,2] +
#                     TBAz(rrA, 4.3, 0.01*i)*df4simp3[,3] + TBAz(rrU, mU, 0.01*i)*df4simp3[,4])
# }
# loon
# dxx <- data.frame('klee'=klee, 'loon' = loon)
# write.csv(dxx,"stuff.csv")

traX <- seq(0.01,1,0.01)
dfTBA <- data.frame(TRA = traX, D = tbaD, T = tbaT, A = tbaA, U = tbaU)
dfTBAm <- melt(dfTBA, id.vars = "TRA")
dfTBA2 <- data.frame(TRA = traX, D = tbaD, T = tbaT, A = tbaA, U = tbaU, ALL = tbaALL)
dfTBA2m <- melt(dfTBA2, id.vars = "TRA")

cb <- brewer.pal(n = 7, name = "BuPu") # Colour scheme
##Fig 1E in article
pl4b <- ggplot() + geom_line(data=dfTBAm, aes(x=TRA, y=value, color=variable), size=1) + themeJDC +
  xlab("Transmission Reducing Activity") + theme(legend.position = 'none') + #labs(color="State") +
  ylab("Transmission Blocking Activity") + #theme(legend.position = c(0.27,0.73), legend.title = element_blank()) +
  #guides(colour = guide_legend(title = "Infection State", override.aes = list(linetype=c(1,1,1,1),shape=c(NA,NA,NA,NA)))) +
  #geom_point(data = shapezM, aes(x = tra, y = value, color = variable, pch = Q), size = 2.9) +
  scale_color_manual(values = c("#c51b8a", cb[7], cb[5], cb[3]),
      labels = c("Symptomatic (Untreated)","Symptomatic (Treated)","Asymptomatic","Subpatent")) +
  annotate('segment', x = 0.9, xend = 0.9, y = 1, yend = 0, linetype = 'dashed')
pl4b

#Supplementary Figure 5A
eqpl4b <- ggplot() + geom_line(data=dfTBA2m, aes(x=TRA, y=value, color=variable), size=1) + themeJDC +
  xlab("Transmission Reducing Activity") + #theme(legend.position = 'none') + #labs(color="State") +
  ylab("Transmission Blocking Activity") + theme(legend.position = c(0.27,0.73), legend.title = element_blank()) +
  #guides(colour = guide_legend(title = "Infection State", override.aes = list(linetype=c(1,1,1,1),shape=c(NA,NA,NA,NA)))) +
  #geom_point(data = shapezM, aes(x = tra, y = value, color = variable, pch = Q), size = 2.9) +
  scale_color_manual(values = c("#c51b8a", cb[7], cb[5], cb[3], 'black'),
        labels = c("Symptomatic (Untreated)","Symptomatic (Treated)","Asymptomatic","Subpatent","ALL"))
eqpl4b

#Make probability mass functions
pmf <- function(i,m,k){

  result <- (1/((1-(((1+m/k))**(-k)))))*((factorial(k+i-1))/((factorial(i))*(factorial((k-1)))))*((1+m/k)**(-i-k))*(m/k)**i

  return <- result
}
#Test function
sq <- 160
t <- pmf(1,mA,rrA)
lu <- seq(1,sq,1)
la <- seq(1,sq,1)
lt <- seq(1,sq,1)
ld <- seq(1,sq,1)
lALL <- seq(1,sq,1)
for(i in 1:sq){
  lu[i] <- pmf(i,mU,rrU)
  la[i] <- pmf(i,mA,rrA)
  lt[i] <- pmf(i,mT,rrT)
  ld[i] <- pmf(i,mD,rrD)
  lALL[i] <- pmf(i,mALL,rrALL)
}
#cumulative
clu <- lu
cla <- la
clt <- lt
cld <- ld
clALL <- lALL
for(i in 2:sq){
  clu[i] <- clu[i-1] + lu[i]
  cla[i] <- cla[i-1] + la[i]
  clt[i] <- clt[i-1] + lt[i]
  cld[i] <- cld[i-1] + ld[i]
  clALL[i] <- clALL[i-1] + lALL[i]
}
dfpmf <- data.frame("i"=seq(1,sq,1), "U"=lu, "A"=la,"T"=lt,"D"=ld)
dfpmfm <- melt(dfpmf,id.vars = "i")
#ggplot(dfpmfm, aes(x=i,y=value, line = variable, color = variable)) + geom_line() + themeJDC +
#  scale_color_manual(values = c("#c51b8a", cb[7], cb[5], cb[3]))

#Same but cumulative
breaks2 <- c(1,10,20,30,40)
dfc <- data.frame("i"=seq(1,sq,1), "U"=clu, "A"=cla,"T"=clt,"D"=cld)
dfcm <- melt(dfc,id.vars = "i")
dfc2 <- data.frame("i"=seq(1,sq,1), "U"=clu, "A"=cla,"T"=clt,"D"=cld,"ALL"=clALL)
dfc2m <- melt(dfc2,id.vars = "i")
#Fig 1D in article
pl4c <- ggplot(dfcm, aes(x=i,y=value, line = variable, color = variable) ) + geom_line(size=1.1) + themeJDC +
  scale_color_manual(values = c(cb[3], cb[5], cb[7],"#c51b8a")) + xlab("Oocyst Count") + ylab("Cumulative Distribution") +
  theme(legend.position = "none", legend.title = element_blank()) +
  scale_x_continuous(breaks = breaks2, limits = c(1,40))
pl4c

#Supplementary Figure 5B
eqpl4c <- ggplot(dfc2m, aes(x=i,y=value, line = variable, color = variable) ) + geom_line(size=1.1) + themeJDC +
  scale_color_manual(values = c(cb[3], cb[5], cb[7],"#c51b8a","black")) + xlab("Oocyst Count") +
  ylab("Cumulative Distribution") +
  theme(legend.position = "none", legend.title = element_blank()) +
  scale_x_continuous(breaks = breaks2, limits = c(1,40))
eqpl4c

 #################################### This section for Figure S3 (line 321) ############################################

#dfpmf <- data.frame("i"=seq(1,sq,1), "U"=lu, "A"=la,"T"=lt,"D"=ld)
#ggplot() + geom_point(data=dfpmf, aes(x=i,y=U),color = cb[3]) + themeJDC + xlim(c(1,30)) +
#  xlab("Oocyst Count") + ylab("Probability Mass Function")
  #scale_color_manual(values = c(cb[3]))

trz <- function(z,n){
  z**n
}
#size of text annotate
sz = 3.6

#U
scU <- 0.55
bray <- trz(0.9,seq(1,sq,1))
bray2 <- bray*lu
bray3 <- data.frame(cent = seq(1,sq,1), prz = scU * bray2) # scaled to fit
bray4 <- data.frame(cent = seq(1,sq,1), prz = scU * bray) # scaled to fit
sm <- sum(bray2)
probU2 <- ggplot() + #dKd, aes(x=num, y=value, fill=variable)
  labs(x = "Oocyst counts", y = "Probability mass function") +
  geom_line(data=dfpmf, aes(x=i,y=U), color = cb[3], size=1.25) + #xlim(c(0.5,15)) +
  themeJDC + geom_line(data = bray4, aes(x=cent, y=prz), color = 'black', size=1.25) +
  scale_y_continuous(sec.axis = sec_axis(~.*(1/scU), name = "Proportion of Transmission Blocked")) +
  scale_x_continuous(breaks = c(1,5,10,15), limits = c(0.8,15)) +
  annotate("text",x=11, y=0.47,label = paste0("For TRA = 0.9, TBA = ",round(sm,2)), size = sz) + ggtitle("U")


#A
scA <- 0.23
bray <- trz(0.9,seq(1,sq,1))
bray2 <- bray*la
bray3 <- data.frame(cent = seq(1,sq,1), prz = scA * bray2) # scaled to fit
bray4 <- data.frame(cent = seq(1,sq,1), prz = scA * bray) # scaled to fit
sm <- sum(bray2)
probA2 <- ggplot() +
  labs(x = "Oocyst counts", y = "Probability mass function") +
  geom_line(data=dfpmf, aes(x=i,y=A), color = cb[5], size=1.25) + #xlim(c(0.5,15)) +
  themeJDC + geom_line(data = bray4, aes(x=cent, y=prz), color = 'black', size=1.25) +
  scale_y_continuous(sec.axis = sec_axis(~.*(1/scA), name = "Proportion of Transmission Blocked")) +
  scale_x_continuous(breaks = c(1,10,20), limits = c(0.8,22)) +
  annotate("text",x=16, y=0.2,label = paste0("For TRA = 0.9, TBA = ",round(sm,2)), size = sz) + ggtitle("A")

#T
scT <- 0.038
bray <- trz(0.9,seq(1,sq,1))
bray2 <- bray*lt
bray3 <- data.frame(cent = seq(1,sq,1), prz = scT * bray2) # scaled to fit
bray4 <- data.frame(cent = seq(1,sq,1), prz = scT * bray) # scaled to fit
sm <- sum(bray2)
probT2 <- ggplot() +
  labs(x = "Oocyst counts", y = "Probability mass function") +
  geom_line(data=dfpmf, aes(x=i,y=T), color = cb[7], size=1.25) + #xlim(c(0.5,120)) +
  themeJDC + geom_line(data = bray4, aes(x=cent, y=prz), color = 'black', size=1.25) +
  scale_y_continuous(sec.axis = sec_axis(~.*(1/scT), name = "Proportion of Transmission Blocked")) +
  scale_x_continuous(breaks = c(1,25,50,75,100), limits = c(0.8,111)) +
  annotate("text",x=75, y=0.03,label = paste0("For TRA = 0.9, TBA = ",round(sm,2)), size = sz) + ggtitle("T")


#D
scD <- 0.026
bray <- trz(0.9,seq(1,sq,1))
bray2 <- bray*ld
bray3 <- data.frame(cent = seq(1,sq,1), prz = scD * bray2) # scaled to fit
bray4 <- data.frame(cent = seq(1,sq,1), prz = scD * bray) # scaled to fit
sm <- sum(bray2)
probD2 <- ggplot() +
  labs(x = "Oocyst counts", y = "Probability mass function") +
  geom_line(data=dfpmf, aes(x=i,y=D), color = "#c51b8a", size=1.25) + #xlim(c(0.5,120)) +
  themeJDC + geom_line(data = bray4, aes(x=cent, y=prz), color = 'black', size=1.25) +
  scale_y_continuous(sec.axis = sec_axis(~.*(1/scD), name = "Proportion of Transmission Blocked")) +
  scale_x_continuous(breaks = c(1,50,100,150), limits = c(0.8,160)) +
  annotate("text",x=100, y=0.022,label = paste0("For TRA = 0.9, TBA = ",round(sm,2)), size = sz) + ggtitle("D")

grid.arrange(probU2, probA2, probT2, probD2, nrow = 2)
