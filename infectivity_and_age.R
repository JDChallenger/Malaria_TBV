#This script will recreate Figs 4A & B: set the EIR (line 28) to those values

library(ICDMM)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(gtable)
library(cowplot)

themeJDC <- theme(panel.background = element_blank(),panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black" ,fill = NA), legend.background = element_blank(),
                  legend.key = element_blank(), legend.title = element_text(), axis.text = element_text(size = 11),
                  axis.title = element_text(size=12), legend.text = element_text(size = 10))

init_age <- seq(0,55,1)#c(seq(0,40,1), seq(50,80,10)) #This plot needs an even spacing in the age vector
#Note: this will be slow, especially on a laptop

# provide the length of time (in days) that you want to run the model for
time_period <- 365*6

# provide a value for the proportion of cases that are treated (referred to as ft in the paper)
prop_treated <- 0.4 # was 0.45

# provide a value of the annual EIR for this model run
init_EIR <- 15
vacc_lag <- 0
ITN_IRS_on <- 4 * 365 + 80
vac_cov <- 0.8 #What should this be

#Note: vaccine turned off for this figure
age_min_rts <- 3
age_max_rts <- which(init_age > 30)[1] - 1 # index of age vector before age is 40 years
age_min_tbv = 3
age_max_tbv = which(init_age > 30)[1] #- 1 # index of age vector before age is 40 years

source("scripts/vaccine_params.R") # Vaccine params that won't change with intervention details

country_str <- "Burkina Faso"	#Note: if seasonality is off in the odin file, doesn't matter what you put here
admin_str <- "Houet" # Hauts-Bassins, bioassay = 0.3.

# creates the odin model
out0 <- run_model(              model = 'odin_model_JDC',
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
                                irs_cov = 0,
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

#Average infectivity over the transmission season (over 6 months peak transmission)
meanT <- init_age
meanD <- init_age
meanU <- init_age
meanAv <- init_age
meanAnv <- init_age
foivSv <- init_age
foiv <- init_age
span <- 180 # time period over which to average
for(i in 1:length(init_age)){
  meanT[i] <- mean(out$tinf_Age[ITN_IRS_on:(ITN_IRS_on+span),i])
  meanD[i] <- mean(out$dinf_Age[ITN_IRS_on:(ITN_IRS_on+span),i])
  meanU[i] <- mean(out$uinf_Age[ITN_IRS_on:(ITN_IRS_on+span),i])
  meanAv[i] <- mean(out$avisinf_Age[ITN_IRS_on:(ITN_IRS_on+span),i])
  meanAnv[i] <- mean(out$ainf_Age[ITN_IRS_on:(ITN_IRS_on+span),i]) -
    mean(out$avisinf_Age[ITN_IRS_on:(ITN_IRS_on+span),i])
  foivSv[i] <- mean(out$FOIv_ages2[ITN_IRS_on:(ITN_IRS_on+span),i])
  foiv[i] <- mean(out$FOIv_ages[ITN_IRS_on:(ITN_IRS_on+span),i])
}

sm <- sum(meanD)+sum(meanT)+sum(meanAv)+sum(meanAnv)+sum(meanU)
#
MeanDpc <- 100*(meanD/sm)
MeanTpc <- 100*(meanT/sm)
MeanAvpc <- 100*(meanAv/sm)
MeanAnvpc <- 100*(meanAnv/sm)
MeanUpc <- 100*(meanU/sm)
#
#FOI per person... Just 'divide' by demography
foipp <- foivSv / out$den[888,]
norm <- sum(foipp)
foipp <- (1/norm)*foipp

cb <- brewer.pal(n = 7, name = "BuPu") # alpha("#c51b8a",1.0)
mx <- 0.0016
#Use geom_rect to separate each 'bin'
df2pc <- data.frame(age = init_age, "Asymptomatic (Subpatent)" = MeanTpc+MeanDpc+MeanAvpc+MeanAnvpc + MeanUpc,
                    "Asymtomatic" = MeanTpc+MeanDpc+MeanAvpc,"Symptomatic (Treated)" = MeanTpc+MeanDpc,
                    "Symptomatic (Untreated)" = MeanDpc)
df2pcm <- melt(df2pc, id.vars = "age")
main_plotXpc <- ggplot(df2pcm)  + xlim(c(0,60.1)) + themeJDC +ylim(c(0,6.1)) +
  geom_rect(aes(xmin = age+0, xmax = age+1, ymin=0, ymax = value, fill = variable)) +
  ylab("Percentage of total infectivity") + xlab("Age (Years)") + theme(legend.position = 'none') +
  scale_fill_manual(values = c(alpha(cb[3],1.0),alpha(cb[5],1.0),alpha(cb[7],1.0),alpha("#c51b8a",1.0)),
  labels = c("Asymptomatic\n(subpatent)","Asymptomatic","Symptomatic\n(Treated)","Symptomatic\n(Untreated)"))+
  geom_rect(data = data.frame(age = init_age, mdetect = (MeanDpc+MeanTpc+MeanAvpc)),
    aes(xmin = age+0.03, xmax = age+0.97, ymin = mdetect, ymax = mdetect), color='black', size =0.8)
main_plotXpc

dfpp <- data.frame(age = init_age, foipp = foipp)
insetpp2 <- ggplot(dfpp, aes(x=age,y=foipp)) + geom_line(color="aquamarine3") + themeJDC + xlim(c(0,60.1)) +
  xlab("Age (Years)") +
  ylab("Per-person contribution") + scale_y_continuous(breaks = c(0,2), limits = c(0,0.036)) +
  theme(legend.position = "none",axis.text = element_text(size = 10),axis.title = element_text(size=10))

#One upper panel of Fig 4
text4bXpc <- ggdraw() + draw_plot(main_plotXpc) + draw_plot(insetpp2, x=0.45, y=0.45, width = 0.5,height = 0.5)
text4bXpc

