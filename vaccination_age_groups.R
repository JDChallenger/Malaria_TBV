#This script generates the lower panels of fig. 4.
#This might be slow/problematic to run on a laptop
#Set the EIR in line 36

library(ICDMM)
library(ggplot2)
library(gganimate)
library(RColorBrewer)
library(reshape2)
library(gridExtra)
library(stringr)
library(wesanderson)

themeJDC <- theme(panel.background = element_blank(),panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black" ,fill = NA), legend.background = element_blank(),
                  legend.key = element_blank(), axis.text = element_text(size = 11),
                  axis.title = element_text(size=12), legend.text = element_text(size = 10),
                  legend.title = element_text(size = 10))

#Age vector: 45 bins with the same number of people in
#In each vaccine scenario, vaccinate 3 bins (one 15th of the population)
vvv <- c(0, 0.45, 0.92, 1.40, 1.89, 2.39, 2.90, 3.43, 3.97, 4.53, 5.10, 5.69, 6.29, 6.92, 7.56,
  8.23, 8.92, 9.63, 10.4, 11.13, 11.93, 12.76, 13.62, 14.52, 15.46, 16.46, 17.5, 18.6, 19.76,
  20.99, 22.3, 23.70, 25.2, 26.82, 28.59, 30.52, 32.66, 35.05, 37.76, 40.89, 44.59, 49.12, 54.96,
  63.19, 77.25)
vl <- length(vvv)
init_age <- vvv

# provide the length of time (in days) that you want to run the model for
time_period <- 365*6

# provide a value for the proportion of cases that are treated
prop_treated <- 0.4

# provide a value of the annual EIR for this model run
init_EIR <- 15
vacc_lag <- 0*365
ITN_IRS_on <- 4 * 365 + 80
vac_cov <- 0.8 #What should this be

age_min_rts <- 3
age_max_rts <- which(init_age > 60)[1]-1
age_min_tbv <- 3
age_max_tbv = which(init_age > 60)[1]-1

source("scripts/vaccine_params.R") # Vaccine params that won't change with intervention details

country_str <- "Burkina Faso"	#Note: if seasonality is off in the odin file, doesn't matter what you put here
admin_str <- "Houet" # Hauts-Bassins, bioassay = 0.3.

### 1. No interventions
out <- run_model(              model = 'odin_model_JDC',
                                time = time_period,
                                het_brackets = 5,
                                age = init_age,
                                num_int = 4,
                                init_EIR = init_EIR,
                                init_ft = prop_treated,
                                country = country_str,
                                admin2 = admin_str,
                                switch_TBV = 0,
                                switch_TRA_to_TBA = 0,
                                irs_cov = 0, # USE irs_cov for vaccine(s) now
                                itn_cov = 0,
                                ITN_IRS_on = ITN_IRS_on,
                                vacc_lag = vacc_lag,
                                hill1 = hill1,
                                hill2 = hill2,
                                mu25 = mu25,
                                tau25 = tau25,
                                rho25 = rho25,
                                ds25 = ds25,
                                dl25 = dl25,
                                v_interval = v_interval,
                                age_min_tbv = age_min_tbv,
                                age_max_tbv = age_max_tbv,
                                RTS_switch = 0,
                                age_min_rts = age_min_rts,
                                age_max_rts = age_max_rts,
                                t_boost_rts = t_boost_rts,
                                CS_peak = CS_peak,
                                CS_boost = CS_boost,
                                p_peak = p_peak,
                                p_boost = p_boost,
                                ds = ds,
                                dl = dl,
                                beta_RTS = beta_RTS,
                                alpha_RTS = alpha_RTS,
                                V_max_RTS = V_max_RTS)

#cumulative incidence, stratified by age
cumul <- out$inc_age
#incidence in whole population
cumulJ <- out$inc
cumulJ05 <- out$inc05
for(k in 2:length(out$inc) ){
  cumul[k,] <-  cumul[k-1,] + cumul[k,]
  cumulJ[k] <- cumulJ[k] + cumulJ[k-1]
  cumulJ05[k] <- cumulJ05[k] + cumulJ05[k-1]
}
#plot(out$t,cumulJ)
#plot(out$t,cumulJ05)

x1 <- (cumul[ITN_IRS_on + 365,] - cumul[ITN_IRS_on,])
dfxx <- data.frame("age"=init_age, "x1"=x1,"x1mod"=x1/out$den[11,])
dfxxm <- melt(dfxx, id.vars = "age")

dfC <- data.frame(cumul)
dfCJ <- data.frame(t=out$t, inc=cumulJ, inc2=cumulJ2)
dfC$t <- out$t
dfCm <- melt(dfC, id.vars = "t")
dfCJm <- melt(dfCJ, id.vars = "t")

#cases in a one-year period after vaccination
cases <- seq(1,length(vvv), 1)
cases <- cumul[ITN_IRS_on + 365,] - cumul[ITN_IRS_on,]
cumulJ[ITN_IRS_on + 365] - cumulJ[ITN_IRS_on]

#
vll <- 15
init_age3 <- init_age[seq(1,vl,3)]
init_age3
CAmatrix <- matrix(0, nrow = vll-1, ncol = vll)
CAmatrixPC <- matrix(0, nrow = vll-1, ncol = vll) # Percentage of cases averted in each age group
CAvecPC <- matrix(0, nrow = vll-1, ncol = 1) # Percentage of cases averted in each age group
CAPC <- seq(1,vll-1,1) # % cases averted full stop
for(i in 2:vll){
  age_min_tbv <- 3*(i-2) + 1
  age_max_tbv <- 3*(i-2) + 3
  print(paste0(age_min_tbv," ",age_max_tbv))

  out2 <- run_model(              model = 'odin_model_JDC',
                                  time = time_period,
                                   het_brackets = 5,
                                   age = init_age,
                                   num_int = 4,
                                   init_EIR = init_EIR,
                                   init_ft = prop_treated,
                                   country = country_str,
                                   admin2 = admin_str,
                                   switch_TBV = 1,
                                   switch_TRA_to_TBA = 1,
                                   irs_cov = 0.999, # USE irs_cov for vaccine(s) now
                                   itn_cov = 0,
                                   ITN_IRS_on = ITN_IRS_on,
                                   vacc_lag = vacc_lag,
                                   hill1 = hill1,
                                   hill2 = hill2,
                                   mu25 = mu25,
                                   tau25 = tau25,
                                   rho25 = rho25,
                                   ds25 = ds25,
                                   dl25 = dl25,
                                   v_interval = v_interval,
                                   age_min_tbv = age_min_tbv,
                                   age_max_tbv = age_max_tbv,
                                   RTS_switch = 0,
                                   age_min_rts = age_min_rts,
                                   age_max_rts = age_max_rts,
                                   t_boost_rts = t_boost_rts,
                                   CS_peak = CS_peak,
                                   CS_boost = CS_boost,
                                   p_peak = p_peak,
                                   p_boost = p_boost,
                                   ds = ds,
                                   dl = dl,
                                   beta_RTS = beta_RTS,
                                   alpha_RTS = alpha_RTS,
                                   V_max_RTS = V_max_RTS)

  #cumulative incidence, stratified by age
  cumul2 <- out2$inc_age
  #Usual incidence
  cumulJ2 <- out2$inc
  for(k in 2:length(out2$inc) ){
    cumul2[k,] <-  cumul2[k-1,] + cumul2[k,]
    cumulJ2[k] <- cumulJ2[k] + cumulJ2[k-1]
  }
  CAPC[i-1] <- 100*(1-((cumulJ2[ITN_IRS_on + 365] - cumulJ2[ITN_IRS_on])/(cumulJ[ITN_IRS_on + 365] - cumulJ[ITN_IRS_on])))
  cases <- cumul[ITN_IRS_on + 365,] - cumul2[ITN_IRS_on + 365,]
  #Group in threes
  cases3 <- seq(1,vll,1)
  cumul_3 <- matrix(0,vll,nrow=(length(out2$inc)))
  cumul2_3 <- matrix(0,vll,nrow=(length(out2$inc)))
  for(j in 1:vll){
    cases3[j] <- cases[3*(j-1) + 1] +cases[3*(j-1) + 2] +cases[3*(j-1) + 3]
    for(k in 1:length(out2$inc)){
      cumul_3[k,j] <- cumul[k,j] + cumul[k,j] + cumul[k,j]
      cumul2_3[k,j] <- cumul2[k,j] + cumul2[k,j] + cumul2[k,j]
    }
  }
  CAmatrix[i-1,] <- cases3
  CAmatrixPC[i-1,]<- 100*(1-(cumul2_3[ITN_IRS_on + 365,] - cumul2_3[ITN_IRS_on,])/(
    cumul_3[ITN_IRS_on + 365,] - cumul_3[ITN_IRS_on,]))
  CAvecPC[i-1] <- 100*(1-(sum(cumul2_3[ITN_IRS_on + 365,]) - sum(cumul2_3[ITN_IRS_on,]))/(
    sum(cumul_3[ITN_IRS_on + 365,]) - sum(cumul_3[ITN_IRS_on,])))
  print(i)
}

CAmatrix <- 1000*CAmatrix # rescale
max(CAmatrix)

rowsmz <- rowSums(CAmatrix)
#now, we'll need to match each of these to the midpoint of each age group
midz <- seq(1,length(init_age3)-1,1)
for(i in 1:(length(init_age3)-1)){
  midz[i] <- 0.5*(init_age3[i] + init_age3[i+1])
}
midz
rmin <- min(rowsmz)
rmax <- max(rowsmz)

#cumulative version of age vector
vchar <- as.character(init_age3)

colnames(CAmatrix) <- vchar[1:vll]
rownames(CAmatrix) <- vchar[1:vll-1]
colnames(CAmatrixPC) <- vchar[1:vll]
rownames(CAmatrixPC) <- vchar[1:vll-1]
CAmatrixM <- melt(CAmatrix)
CAmatrixPCM <- melt(CAmatrixPC)

CAmatrixM$Var1 <- as.numeric(CAmatrixM$Var1)
CAmatrixM$Var2 <- as.numeric(CAmatrixM$Var2)
CAmatrixM$xmin <- NA
CAmatrixM$xmax <- NA
CAmatrixM$ymin <- NA
CAmatrixM$ymax <- NA

CAmatrixPCM$Var1 <- as.numeric(CAmatrixPCM$Var1)
CAmatrixPCM$Var2 <- as.numeric(CAmatrixPCM$Var2)
CAmatrixPCM$xmin <- NA
CAmatrixPCM$xmax <- NA
CAmatrixPCM$ymin <- NA
CAmatrixPCM$ymax <- NA

for(i in 1:(vll*(vll-1))){
  CAmatrixM[i,4] <- CAmatrixM[i,1]
  CAmatrixPCM[i,4] <- CAmatrixPCM[i,1]
  CAmatrixM[i,6] <- CAmatrixM[i,2]
  CAmatrixPCM[i,6] <- CAmatrixPCM[i,2]
}

for(i in 1:(vll*(vll-1))){

  CAmatrixM[i,7] <- CAmatrixM[(i+vll-1),2]
  CAmatrixPCM[i,7] <- CAmatrixPCM[(i+vll-1),2]

}

xx <- init_age3[2:vll]
xxx <- rep(xx,vll)
CAmatrixM$xmax <- xxx
CAmatrixPCM$xmax <- xxx

CAmatrixM2 <- CAmatrixM[-seq(((vll-1)*(vll-1))+1,vll*(vll-1),1),]
CAmatrixPCM2 <- CAmatrixPCM[-seq(((vll-1)*(vll-1)),vll*(vll-1),1),]

#If you want to wrap percentage figures in brackets and add a '%' you'll need a function
wrap <- function(x){
  a <- paste0("(",as.character(round(x)),"%)")
}
wrap2 <- function(x){
  a <- paste0(as.character(round(x)),"%")
}
b <- wrap(5)
pcv <- lapply(CAPC,wrap2)
pcva <- unlist(pcv, use.names = FALSE)

# Gradient color
pal <- wes_palette("Zissou1", 100, type = "continuous")

dfmz <- data.frame(mz = midz, ca = round(CAPC[1:14]), pcvX = pcva, rs = round(rowsmz))

#One Lower panel of Figure 4 (set EIR to desired value at the beginning of this script)
ggplot() + labs(fill = "Cases \naverted per\n1000 people") +
  geom_rect(data = CAmatrixM2, aes(fill = value, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), colour = "white")+
  xlim(c(0,27)) + ylim(c(0,24.2)) +
  xlab("Age of vaccinees (years)") + ylab("Age group in which cases are averted (years)") + themeJDC +
  scale_fill_distiller(palette = "Spectral")+#, limits = c(0.0, 8.8)) + #ggtitle("TBV Only (EIR=20)") +
  theme(plot.margin = unit(c(.2,0.01,.3,.1), "cm")) +
  geom_text(data=data.frame(mz = midz, rs = round(rowsmz)), aes(label = rs, x = mz, y = 23.9,alpha = rs)) +
  geom_text(data=dfmz,
            aes(label = pcvX, x = mz, y = 22.9,alpha = rs), size = 3.75) +
  scale_alpha_continuous(range = c(0.5,1.0), guide = 'none')


