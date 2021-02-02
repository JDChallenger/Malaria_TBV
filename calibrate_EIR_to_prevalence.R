library(ICDMM)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(cowplot)

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
prop_treated <- 0.4 #
#net coverage
net_cov <- 0.7

# PART 1: KLESSO

# provide a value of the annual EIR for this model run
init_EIR <- 490 #Fitted further down this script (lines 111-149)

country_str <- "Burkina Faso"
admin_str <- "Houet"

# creates the odin model
out <- run_model(   model = 'odin_model_JDC',
                       het_brackets = 5,
                       age = init_age,
                       time = time_period,
                       init_EIR = init_EIR,
                       country = country_str,
                       admin2 = admin_str,
                       itn_cov = net_cov,
                       irs_cov = 0,
                       num_int = 4,
                       bites_Bed = 0.79,
                       ITN_IRS_on = 1*365,
                       init_ft = prop_treated,
                       d_ITN0 = 0.450,
                       r_ITN0 = 0.359,
                       itn_half_life = 1.821*365, #half life in days not years
                       switch_TBV = 0,
                       switch_TRA_to_TBA = 0,
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

#Malaria prevalence (by microscopy)
dfx <- data.frame(tt=out$t[(8*365+1):(9*365)] - 8*365, prev=out$prev_tot[(8*365+1):(9*365)])

#Malaria status of human population over the course of one year
dfs <- data.frame(out$Dout[(8*365+1):(9*365)] , out$Tout[(8*365+1):(9*365)], out$Aoutvis[(8*365+1):(9*365)],
                      out$Aout[(8*365+1):(9*365)] - out$Aoutvis[(8*365+1):(9*365)] +
                        out$Uout[(8*365+1):(9*365)], out$Pout[(8*365+1):(9*365)] +
                        out$Sout[(8*365+1):(9*365)])

dfs$day <- seq(1,1*365,1)

cb <- brewer.pal(n = 7, name = "BuPu")

dfsm <- melt(dfs, id.vars = "day")

#Target data
target <- c(0.32,0.53,0.5,0.42)
CI <- 1.96 * sqrt(target*(1-target)*(1/300))
targetTimes <- c(120,210,265,330) #
tdf <- data.frame(target = target , targetTimes = targetTimes, ci = CI)

#Theme for plot
themeJDC <- theme(panel.background = element_blank(),panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black" ,fill = NA), legend.background = element_blank(),
                  legend.key = element_blank(), legend.title = element_text(), axis.text = element_text(size = 13),
                  axis.title = element_text(size=13), legend.text = element_text(size = 10))


#Fig 3A
KlessoSimple <- ggplot() + geom_area(data = dfsm, aes(x = day, y = value, fill = variable), #colour = "white",
                     position = position_stack(reverse = T)) +
  ylab("Proportion of population") + xlab("Month") +
  #ggtitle("Klesso [Houet]") +
  themeJDC + #guides(fill = guide_legend()) +
  geom_line(data = dfx, aes(x = tt, y = prev), size = 1.25) +
  geom_point(data = tdf,aes(x=targetTimes,y=target)) + theme(legend.direction = "vertical", legend.position = "right",
   legend.spacing.x = unit(0.25, 'cm'), legend.key.size = unit(0.99, 'cm'), legend.text = element_text(size = 11)) +
  geom_errorbar(data = tdf,width=.1, aes(x=targetTimes,ymin=target-ci,ymax=target+ci))+
  scale_fill_manual(values = c(alpha("#c51b8a",1.0),alpha(cb[7],1.0),alpha(cb[5],1.0),
                               alpha(cb[3],1.0),alpha(cb[1],1.0)),
                    labels = c("Symptomatic \n(Untreated)","Symptomatic \n(Treated)","Asymptomatic",
                               "Asymptomatic \n(subpatent)","Uninfected")) +
  guides(colour = guide_legend(nrow=1), fill=guide_legend(title="Infection State",reverse = T)) +
  scale_x_continuous(breaks = c(59,151,243,335), labels = c("Mar","Jun","Sept","Dec")) +
  annotate("segment", x = 61, xend = 120, y = 0.95, yend = 0.95, colour = "black", size=1) +
  annotate("text", x = 91, y = 0.98, colour = "black", size=5, label = "1") +
  annotate("segment", x = 180, xend = 210, y = 0.95, yend = 0.95, colour = "black", size=1) +
  annotate("text", x = 195, y = 0.98, colour = "black", size=5, label = "2") +
  annotate("segment", x = 235, xend = 265, y = 0.95, yend = 0.95, colour = "black", size=1) +
  annotate("text", x = 250, y = 0.98, colour = "black", size=5, label = "3") +
  annotate("segment", x = 305, xend = 335, y = 0.95, yend = 0.95, colour = "black", size=1) +
  annotate("text", x = 320, y = 0.98, colour = "black", size=5, label = "4")
KlessoSimple

#vary EIR, and store sum of squares to assess best EIR

vary_eir <- seq(370,590,10)
targetTimes <- c(120,210,265,330)
targetTimes <- targetTimes + (8*365)
store_ssq <- vary_eir

for(j in 1:length(vary_eir)){

    out <- run_model(   model = 'odin_model_JDC',
                      het_brackets = 5,
                      age = init_age,
                      time = time_period,
                      init_EIR = vary_eir[j],
                      country = country_str,
                      admin2 = admin_str,
                      itn_cov = net_cov,
                      irs_cov = 0,
                      num_int = 4,
                      bites_Bed = 0.79,
                      ITN_IRS_on = 1*365,
                      init_ft = prop_treated,
                      d_ITN0 = 0.450,
                      r_ITN0 = 0.359,
                      itn_half_life = 1.821*365, #half life in days not years
                      switch_TBV = 0,
                      switch_TRA_to_TBA = 0,
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

  store_ssq[j] <- sqrt((target[1] - out$prev_tot[targetTimes[1]])^2 + (target[2] - out$prev_tot[targetTimes[2]])^2 +
    (target[3] - out$prev_tot[targetTimes[3]])^2 + (target[4] - out$prev_tot[targetTimes[4]])^2)

  print(j)
}#End of for loop
plot(vary_eir, store_ssq)
ve <- which.min(store_ssq)
vary_eir[ve]


# PART 2: Longo

# Provide a value of the annual EIR for this model run
init_EIR <- 195 #Fitted further down (lines 232-271)
time_period <- 365*10

country_str <- "Burkina Faso"	#Note: if seasonality is off in the odin file, doesn't matter what you put here
admin_str <- "Koulpelogo" # Centre-Est bioassay=0.4

# creates the odin model
out <- run_model(   model = 'odin_model_JDC',
                                het_brackets = 5,
                                age = init_age,
                                init_EIR = init_EIR,
                                country = country_str,
                                admin2 = admin_str,
                                itn_cov = net_cov,
                                time = time_period,
                                num_int = 4,
                                ITN_IRS_on = 1*365,
                                bites_Bed = 0.79,
                                init_ft = prop_treated,
                                d_ITN0 = 0.427,
                                r_ITN0 = 0.372,
                                itn_half_life = 1.485*365,
                                switch_TBV = 0,
                                switch_TRA_to_TBA = 0,
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

#Malaria prevalence (by microscopy)
dfx <- data.frame(tt=out$t[(8*365+1):(9*365)] - 8*365, prev=out$prev_tot[(8*365+1):(9*365)])

#Malaria status of human population over the course of one year
dfs2 <- data.frame(out$Dout[(8*365+1):(9*365)] , out$Tout[(8*365+1):(9*365)], out$Aoutvis[(8*365+1):(9*365)],
                      out$Aout[(8*365+1):(9*365)] - out$Aoutvis[(8*365+1):(9*365)] +
                         out$Uout[(8*365+1):(9*365)], out$Pout[(8*365+1):(9*365)] +
                         out$Sout[(8*365+1):(9*365)])

dfs2$day <- seq(1,1*365,1)

dfs2m <- melt(dfs2, id.vars = "day")

target <- c(0.21,0.43,0.52,0.35)
CI <- 1.96 * sqrt(target*(1-target)*(1/300))
targetTimes <- c(120,210,265,330)
tdf <- data.frame(target = target , targetTimes = targetTimes, ci = CI)

#Suppl Figure 1
LongoSimp <- ggplot() + geom_area(data = dfs2m, aes(x = day, y = value, fill = variable), #colour = "white",
                     position = position_stack(reverse = T)) +
  ylab("Proportion of population") + xlab("Month") +
  #ggtitle("Longo") +
  themeJDC + guides(fill = guide_legend(reverse = T)) +
  geom_line(data = dfx, aes(x = tt, y = prev), size = 1.25) +
  geom_point(data = tdf,aes(x=targetTimes,y=target)) + theme(legend.title = element_blank(),
    legend.direction = "vertical", legend.position = "right", legend.spacing.x = unit(0.25, 'cm'),
    legend.key.size = unit(0.99, 'cm')) +
  geom_errorbar(data = tdf,width=.1, aes(x=targetTimes,ymin=target-ci,ymax=target+ci))+
  scale_fill_manual(values = c(alpha("#c51b8a",1.0),alpha(cb[7],1.0),alpha(cb[5],1.0),
                               alpha(cb[3],1.0),alpha(cb[1],1.0)),
                    labels = c("Symptomatic \n(Untreated)","Symptomatic \n(Treated)","Asymptomatic",
                               "Asymptomatic \n(subpatent)","Uninfected")) +
  guides(colour = guide_legend(nrow=1)) +
  scale_x_continuous(breaks = c(59,151,243,335), labels = c("Mar","Jun","Sept","Dec")) +
  annotate("segment", x = 61, xend = 120, y = 0.95, yend = 0.95, colour = "black", size=1) +
  annotate("text", x = 91, y = 0.98, colour = "black", size=5, label = "1") +
  annotate("segment", x = 180, xend = 210, y = 0.95, yend = 0.95, colour = "black", size=1) +
  annotate("text", x = 195, y = 0.98, colour = "black", size=5, label = "2") +
  annotate("segment", x = 235, xend = 265, y = 0.95, yend = 0.95, colour = "black", size=1) +
  annotate("text", x = 250, y = 0.98, colour = "black", size=5, label = "3") +
  annotate("segment", x = 305, xend = 335, y = 0.95, yend = 0.95, colour = "black", size=1) +
  annotate("text", x = 320, y = 0.98, colour = "black", size=5, label = "4")
LongoSimp

vary_eir <- seq(170,220,5)
targetTimes <- c(120,210,265,330)
targetTimes <- targetTimes + (8*365)

store_ssq <- vary_eir

for(j in 1:length(vary_eir)){
  out <- run_model(   model = 'odin_model_JDC',
                      het_brackets = 5,
                      age = init_age,
                      init_EIR = vary_eir[j],
                      country = country_str,
                      admin2 = admin_str,
                      itn_cov = net_cov,
                      time = time_period,
                      num_int = 4,
                      ITN_IRS_on = 1*365,
                      bites_Bed = 0.79,
                      init_ft = prop_treated,
                      d_ITN0 = 0.427,
                      r_ITN0 = 0.372,
                      itn_half_life = 1.485*365,
                      switch_TBV = 0,
                      switch_TRA_to_TBA = 0,
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

  store_ssq[j] <- sqrt((target[1] - out$prev_tot[targetTimes[1]])^2 + (target[2] - out$prev_tot[targetTimes[2]])^2 +
                         (target[3] - out$prev_tot[targetTimes[3]])^2 + (target[4] - out$prev_tot[targetTimes[4]])^2)

  print(j)
}#End of for loop

plot(vary_eir, store_ssq)
ve <- which.min(store_ssq)
vary_eir[ve]
