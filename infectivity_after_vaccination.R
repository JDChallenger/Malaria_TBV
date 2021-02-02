#Script generates Figure 5a and 5b

library(ICDMM)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

# create a vector of age categories
v1 <- seq(0,1.75,0.25)
v2 <- seq(2,5,1)
v3 <- seq(6,18,2)
v4 <- seq(20,80,10)
v <- c(v1,v2,v3,v4)
length(v)
init_age <- v

# provide the length of time (in days) that you want to run the model for
time_period <- 365*10

# provide a value for the proportion of cases that are treated
prop_treated <- 0.4

# provide a value of the annual EIR for this model run
init_EIR <- 50
vacc_lag <- 0
ITN_IRS_on <- 6 * 365 + 100
vac_cov <- 0.8

age_min_rts <- 12
age_max_rts <- which(init_age > 70)[1] - 1 # index of age vector before age is 40 years
age_min_tbv = 3
age_max_tbv = which(init_age > 79)[1] #- 1 # index of age vector before age is 40 years

source("scripts/vaccine_params.R") # Vaccine params that won't change with intervention details

country_str <- "Burkina Faso"	#Note: if seasonality is off in the odin file, doesn't matter what you put here
admin_str <- "Houet"

out0 <- run_model(model = 'odin_model_TBV',
                  het_brackets = 5,
                  age = init_age,
                  num_int = 4,
                  init_EIR = init_EIR,
                  time = time_period,
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

# 1. Model with intervention

out <- run_model(model = 'odin_model_TBV',
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

y1 <- 0
y2 <- 10
df4 <- data.frame(out$Dout[(y1*365+1):(y2*365)], out$Tout[(y1*365+1):(y2*365)], out$Aoutvis[(y1*365+1):(y2*365)],
                  out$Aout[(y1*365+1):(y2*365)] - out$Aoutvis[(y1*365+1):(y2*365)] + out$Uout[(y1*365+1):(y2*365)],
                  out$Sout[(y1*365+1):(y2*365)] + out$Pout[(y1*365+1):(y2*365)])

df_I5 <- data.frame(out$Dinf[(y1*365+1):(y2*365)] , out$Tinf[(y1*365+1):(y2*365)],
                    out$Avisinf[(y1*365+1):(y2*365)],
                    (out$Ainf[(y1*365+1):(y2*365)]-out$Avisinf[(y1*365+1):(y2*365)])+ out$Uinf[(y1*365+1):(y2*365)])

df4$day <- seq((1/365) - (ITN_IRS_on/365),(y2-y1) - (ITN_IRS_on/365) ,(1/365))

df_I5$day <- seq((1/365) - (ITN_IRS_on/365),(y2-y1) - (ITN_IRS_on/365) ,(1/365))
df_prev <- data.frame(day = (1/365)*(out$t-ITN_IRS_on), p0 = out0$prev_tot, p1 = out$prev_tot)
cb <- brewer.pal(n = 7, name = "BuPu")

namesI5 <- c("D","T","A+","A- + U","day")
colnames(df_I5) <- namesI5

dfm4 <- melt(df4, id.vars = "day")
dfmprev <- melt(df_prev, id.vars = "day")
dfmI5 <- melt(df_I5, id.vars = "day")

themeJDC <- theme(panel.background = element_blank(),panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), legend.title = element_blank(),
                  panel.border = element_rect(colour = "black" ,fill = NA), legend.background = element_blank(),
                  legend.key = element_blank(), axis.text = element_text(size = 12),
                  axis.title = element_text(size=11.2), legend.text = element_text(size = 10.2))#, legend.title = element_text())

#Figure 5b
sc <- 0.055
eirDs2 <- data.frame(day = df_I5$day, eir = sc * out$EIRout[(y1*365+1):(y2*365)],
                     micI = out$Dinf[(y1*365+1):(y2*365)] + out$Tinf[(y1*365+1):(y2*365)] + out$Avisinf[(y1*365+1):(y2*365)])
eirDs2m <- melt(eirDs2, id.vars = "day")
fig5b <- ggplot() +xlab("Time after vaccine introduced (Years)")+#ylim(c(0,0.035)) +
  geom_area(data=dfmI5, aes(x = day, y = value, fill = variable),position = position_stack(reverse=T)) +
  ylab("Infectivity of population") + themeJDC +# theme(aspect.ratio = 0.25) +
  theme(legend.title = element_blank(),legend.position = c(0.825,0.85), legend.direction = "vertical",
        legend.spacing.x = unit(0.3, 'cm'),  legend.key.size = unit(0.81, 'cm')) + #
  xlim(c(-2, 0.995)) + scale_fill_manual(values = c(alpha("#c51b8a",1.0),alpha(cb[7],1.0),alpha(cb[5],1.0),alpha(cb[3],1.0)),
                                         labels = c("Symptomatic\n(Untreated)","Symptomatic\n(Treated)","Asymptomatic","Asymptomatic\n(subpatent)")) +
  guides(colour = guide_legend(nrow=1)) + geom_line(data = eirDs2m, aes(x=day, y=value, linetype = variable),
                                                    size = 1.0, show.legend = F) +
  scale_y_continuous(sec.axis = sec_axis(~.*(1/sc), name = "Mean Daily EIR")) +
  scale_linetype_manual(values  = c("dashed","solid")) +
  annotate("segment", x = 0, xend = 0, y=0, yend = 0.029, colour = 'grey') + guides(fill = guide_legend(reverse = T))
fig5b

#Figure 5a
fig5a <- ggplot() + geom_area(data=dfm4, aes(x = day, y = value, fill = variable),position = position_stack(reverse=T)) +
  ylab("Proportion of population") + xlim(c(-2,1)) + xlab("Time after vaccine introduced (years)") + themeJDC +
  scale_fill_manual(values = c(alpha("#c51b8a",1.0),alpha(cb[7],1.0),alpha(cb[5],1.0),
                               alpha(cb[3],1.0),alpha(cb[1],1.0)), labels = c("Symptomatic \n(Untreated)","Symptomatic \n(Treated)","Asymptomatic",
                                                                              "Asymptomatic \n(subpatent)","Uninfected")) + geom_line(data = dfmprev, aes(x=day,y=value,linetype=variable),
                                                                                                                                      size = 1.0, show.legend = F) + guides(fill = guide_legend(reverse = T)) +
  scale_linetype_manual(values  = c("dotdash","solid")) + theme(legend.position = c(0.2,0.8), legend.background = element_rect()) +
  annotate("segment", x = 0, xend = 0, y=0, yend = 1, colour = "grey")#, linetype = "dashed")
fig5a
