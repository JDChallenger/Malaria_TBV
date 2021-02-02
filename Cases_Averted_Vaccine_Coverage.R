#script generates fig 5c

library(ICDMM)
library(ggplot2)
library(reshape2)
library(gridExtra)

v1 <- seq(0,1.75,0.25)
v2 <- seq(2,5,1)
v3 <- seq(6,18,2)
v4 <- seq(20,80,10)
v <- c(v1,v2,v3,v4)
length(v)
init_age <- v

# provide the length of time (in days) that you want to run the model for
time_period <- 365*9

# provide a value for the proportion of cases that are treated (referred to as ft in the paper)
prop_treated <- 0.4 # was 0.45

# provide a value of the annual EIR for this model run
init_EIR <- 50 # down from 150 (100)
vacc_lag <- 0
ITN_IRS_on <- 6 * 365 + 100
itn_cov <- 0

age_min_rts <- 12
age_max_rts <- which(init_age > 70)[1] - 1 # index of age vector before age is 40 years
age_min_tbv = 3
age_max_tbv = which(init_age > 79)[1] #- 1 # index of age vector before age is 40 years

source("scripts/vaccine_params.R") # Vaccine params that won't change with intervention details

country_str <- "Burkina Faso"
admin_str <- "Houet"

## No-intervention model
out0 <- run_model(   model = 'odin_model_TBV',
                                 het_brackets = 5,
                                 age = init_age,
                                 num_int = 4,
                                 time = time_period,
                                 init_EIR = init_EIR,
                                 init_ft = prop_treated,
                                 country = country_str,
                                 admin2 = admin_str,
                                 switch_TBV = 1,
                                 switch_TRA_to_TBA = 1,
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

#Cumulative incidence
cumul0 <- out0$inc
for(k in 2:length(out0$inc)){
  cumul0[k]<-  cumul0[k-1] + cumul0[k]
}
cases0 <- cumul0[ITN_IRS_on + 365] - cumul0[ITN_IRS_on]
cases0

vac_coverage <- seq(0.2,0.9,0.10)
CA1 <- seq(0.2,0.9,0.10)
CA2 <- seq(0.2,0.9,0.10)
CA3 <- seq(0.2,0.9,0.10)
CA4 <- seq(0.2,0.9,0.10)
CA5 <- seq(0.2,0.9,0.10)
CA6 <- seq(0.2,0.9,0.10)

### Loop over vaccine coverage for a number of modelling scenarios. list here:
# 1. TBA Model (normal titre)
# 2. Naive TBA Model (normal titre)
# 3. TBA Model (better titre)
# 4. Naive TBA Model (better titre)
# 5. TBA Model (worse titre)
# 6. Naive TBA Model (worse titre)

age_min_tbvB = 3
age_max_tbvB = which(init_age > 79)[1] - 1

for(i in 1:length(vac_coverage)){
  vac_cov <- vac_coverage[i]
  print(vac_cov)

  #1.
  out1 <- run_model(   model = 'odin_model_TBV',
                                  het_brackets = 5,
                                  age = init_age,
                                  num_int = 4,
                                  init_EIR = init_EIR,
                                  time = time_period,
                                  init_ft = prop_treated,
                                  country = country_str,
                                  admin2 = admin_str,
                                  switch_TBV = 1,
                                  switch_TRA_to_TBA = 1,
                                  irs_cov = vac_cov,
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
                                  age_min_tbv = age_min_tbvB,
                                  age_max_tbv = age_max_tbvB,
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

  #2.
  out2 <- run_model(   model = 'odin_model_TBV',
                                   het_brackets = 5,
                                   age = init_age,
                                   num_int = 4,
                                   init_EIR = init_EIR,
                                   time = time_period,
                                   init_ft = prop_treated,
                                   country = country_str,
                                   admin2 = admin_str,
                                   switch_TBV = 1,
                                   switch_TRA_to_TBA = 0,
                                   irs_cov = vac_cov,
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
                                   age_min_tbv = age_min_tbvB,
                                   age_max_tbv = age_max_tbvB,
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

  #3.
  out3 <- run_model(   model = 'odin_model_TBV',
                                   het_brackets = 5,
                                   age = init_age,
                                   num_int = 4,
                                   time = time_period,
                                   init_EIR = init_EIR,
                                   init_ft = prop_treated,
                                   country = country_str,
                                   admin2 = admin_str,
                                   switch_TBV = 1,
                                   switch_TRA_to_TBA = 1,
                                   irs_cov = vac_cov,
                                   itn_cov = 0,
                                   ITN_IRS_on = ITN_IRS_on,
                                   vacc_lag = vacc_lag,
                                   hill1 = hill1,
                                   hill2 = hill2,
                                   mu25 = mu25,
                                   tau25 = tau25,
                                   rho25 = rho25,
                                   ds25 = ds25better,
                                   dl25 = dl25better,
                                   v_interval = v_interval,
                                   age_min_tbv = age_min_tbvB,
                                   age_max_tbv = age_max_tbvB,
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

  #4.
  out4 <- run_model(   model = 'odin_model_TBV',
                                   het_brackets = 5,
                                   age = init_age,
                                   num_int = 4,
                                   time = time_period,
                                   init_EIR = init_EIR,
                                   init_ft = prop_treated,
                                   country = country_str,
                                   admin2 = admin_str,
                                   switch_TBV = 1,
                                   switch_TRA_to_TBA = 0,
                                   irs_cov = vac_cov,
                                   itn_cov = 0,
                                   ITN_IRS_on = ITN_IRS_on,
                                   vacc_lag = vacc_lag,
                                   hill1 = hill1,
                                   hill2 = hill2,
                                   mu25 = mu25,
                                   tau25 = tau25,
                                   rho25 = rho25,
                                   ds25 = ds25better,
                                   dl25 = dl25better,
                                   v_interval = v_interval,
                                   age_min_tbv = age_min_tbvB,
                                   age_max_tbv = age_max_tbvB,
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

  #5.
  out5 <- run_model(   model = 'odin_model_TBV',
                                   het_brackets = 5,
                                   age = init_age,
                                   num_int = 4,
                                   time = time_period,
                                   init_EIR = init_EIR,
                                   init_ft = prop_treated,
                                   country = country_str,
                                   admin2 = admin_str,
                                   switch_TBV = 1,
                                   switch_TRA_to_TBA = 1,
                                   irs_cov = vac_cov,
                                   itn_cov = 0,
                                   ITN_IRS_on = ITN_IRS_on,
                                   vacc_lag = vacc_lag,
                                   hill1 = hill1,
                                   hill2 = hill2,
                                   mu25 = mu25,
                                   tau25 = tau25,
                                   rho25 = rho25,
                                   ds25 = ds25worse,
                                   dl25 = dl25worse,
                                   v_interval = v_interval,
                                   age_min_tbv = age_min_tbvB,
                                   age_max_tbv = age_max_tbvB,
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

  #6.
  out6 <- run_model(   model = 'odin_model_TBV',
                                   het_brackets = 5,
                                   age = init_age,
                                   num_int = 4,
                                   time = time_period,
                                   init_EIR = init_EIR,
                                   init_ft = prop_treated,
                                   country = country_str,
                                   admin2 = admin_str,
                                   switch_TBV = 1,
                                   switch_TRA_to_TBA = 0,
                                   irs_cov = vac_cov,
                                   itn_cov = 0,
                                   ITN_IRS_on = ITN_IRS_on,
                                   vacc_lag = vacc_lag,
                                   hill1 = hill1,
                                   hill2 = hill2,
                                   mu25 = mu25,
                                   tau25 = tau25,
                                   rho25 = rho25,
                                   ds25 = ds25worse,
                                   dl25 = dl25worse,
                                   v_interval = v_interval,
                                   age_min_tbv = age_min_tbvB,
                                   age_max_tbv = age_max_tbvB,
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


  #Cumulative incidence
  cumul1 <- out1$inc
  for(k in 2:length(out1$inc)){
    cumul1[k]<-  cumul1[k-1] + cumul1[k]
  }
  cases1 <- cumul1[ITN_IRS_on + 365] - cumul1[ITN_IRS_on]

  cumul2 <- out2$inc
  for(k in 2:length(out2$inc)){
    cumul2[k]<-  cumul2[k-1] + cumul2[k]
  }
  cases2 <- cumul2[ITN_IRS_on + 365] - cumul2[ITN_IRS_on]

  cumul3 <- out3$inc
  for(k in 2:length(out3$inc)){
    cumul3[k]<-  cumul3[k-1] + cumul3[k]
  }
  cases3 <- cumul3[ITN_IRS_on + 365] - cumul3[ITN_IRS_on]

  cumul4 <- out4$inc
  for(k in 2:length(out4$inc)){
    cumul4[k]<-  cumul4[k-1] + cumul4[k]
  }
  cases4 <- cumul4[ITN_IRS_on + 365] - cumul4[ITN_IRS_on]

  cumul5 <- out5$inc
  for(k in 2:length(out5$inc)){
    cumul5[k]<-  cumul5[k-1] + cumul5[k]
  }
  cases5 <- cumul5[ITN_IRS_on + 365] - cumul5[ITN_IRS_on]

  cumul6 <- out6$inc
  for(k in 2:length(out6$inc)){
    cumul6[k]<-  cumul6[k-1] + cumul6[k]
  }
  cases6 <- cumul6[ITN_IRS_on + 365] - cumul6[ITN_IRS_on]

  CA1[i] <- cases0 - cases1
  CA2[i] <- cases0 - cases2
  CA3[i] <- cases0 - cases3
  CA4[i] <- cases0 - cases4
  CA5[i] <- cases0 - cases5
  CA6[i] <- cases0 - cases6

  print(i)
}

dfc <- data.frame("Coverage" = vac_coverage, c1 = CA1, c2 = CA2, c3 = CA3, c4 = CA4, c5 = CA5, c6 = CA6)
dfc1000 <- data.frame("Coverage" = vac_coverage, c1 = 1000*CA1, c2 = 1000*CA2, c3 = 1000*CA3,
                      c4 = 1000*CA4, c5 = 1000*CA5, c6 = 1000*CA6)
namez2 <- c("Coverage","RTSS-like duration (Adjusted TBA)","RTSS-like duration (Lab TBA)",
           "Longer duration (Adjusted TBA)",
           "Longer duration (Lab TBA)","Shorter duration (Adjusted TBA)",
           "Shorter duration (Lab TBA)")
names(dfc) <- namez2
names(dfc1000) <- namez2
dfc_m <- melt(dfc, id.vars = "Coverage")
dfc1000_m <- melt(dfc1000, id.vars = "Coverage")

### Loop over vaccine coverage for a number of modelling scenarios. list here:
### S.A.C. throughout
# 1. TBA Model (normal titre)
# 2. Naive TBA Model (normal titre)
# 3. TBA Model (better titre)
# 4. Naive TBA Model (better titre)
# 5. TBA Model (worse titre)
# 6. Naive TBA Model (worse titre)

fig5c <- ggplot(data = dfc1000_m, aes(x = Coverage,y=value, color = variable, linetype = variable, size = variable)) +
  geom_line() + themeJDC +
  theme(legend.key.width = unit(3,"line"),legend.position = c(0.34,0.81), axis.text = element_text(size = 12),
        axis.title = element_text(size=12), plot.margin = unit(c(0.2,0.2,0.2,0.45),"cm"),
        legend.text = element_text(size = 10.5)) + ylim(c(0,445)) +
  ylab("Cases averted per 1000 people") + xlab("Vaccine coverage") +
  scale_color_manual(values = c("#1b9e77","#1b9e77","#d95f02","#d95f02","#7570b3","#7570b3"),
                     name = "TBV Model", labels = namez2[2:7]) +
  scale_linetype_manual(values = c("solid","dotted","solid","dotted","solid","dotted"),
                        labels = namez2[2:7], name = "TBV Model") +
  scale_size_manual(values = rep(0.9,9), labels = namez2[2:7], name = "TBV Model")
fig5c

#If you're combining with fig5a and fig5b (and you've already run that script-- It's called 'infectivity_after_vaccination.R')
#you can combine then together, as per the paper
#grid.arrange(fig3a, fig3b,fig3c, nrow=1)


