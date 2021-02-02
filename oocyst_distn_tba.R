#Script for figure 2

library(ggplot2)
library(reshape2)
library(gridExtra)
library(plyr)
library(cowplot)

themeJDC <- theme(panel.background = element_blank(),panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black" ,fill = NA),
  legend.background = element_blank(), legend.key = element_blank(),
  axis.text = element_text(size = 11),
  axis.title = element_text(size=12), legend.text = element_text(size = 10),
  legend.position = c(0.7,0.8))

#Data available through Dryad. You'll need to provide the folder than contains the data file.
dat <- read.delim("/Users/jchallen/Dropbox/NC_data_code/Oocyst_Data/act2datasetRpcr.txt", stringsAsFactors = FALSE)
dr <-dat[dat$Location == 'Klesso' & dat$day=="D7",]$Normal
head(dr)
table(dr)
dr1 <- dr[dr>0 ] #& dr$X0<40
hist(dr1)
str(dr1)
len <- length(dr1)

sum(dr1>99)
dr2 <- c(dr1,rep(100,8))

#histogram
kl <- ggplot(data = data.frame("x"=dr2 ), aes(x=x)) + geom_histogram(bins=100, fill = 'chartreuse3') + themeJDC +
  ylab("Mosquitoes Collected") + xlab("No. of oocysts counted") + ylim(0,95) +
  scale_x_continuous(breaks = c(1,25,50,75,100), labels = c("1","25","50","75",">100"), limits = c(0.5,101.8))
kl

rp <- 1000
tra <- seq(0.1,0.9,0.1)
block <- rep(0,rp*length(tra))
tra2 <- rep(0,rp*length(tra))
count2 <- 1
for(i in 1:length(tra)){
  for(k in 1:rp){
    count <- 0
    for(j in 1:len){
      aux <- rbinom(n=1,dr1[j],prob=1-tra[i])
      if(aux==0)
        count <- count + 1
    }
    #print(count/len)
    block[count2] <- count/len
    tra2[count2] <- tra[i]
    count2 <- count2 + 1
  }
}
dff <- data.frame('tba'=block,tra=tra2)
getmean <- function(trax){
  mean(dff$tba[dff$tra==trax])
}
v <- c(getmean(0.1), getmean(0.2), getmean(tra[3]), getmean(0.4), getmean(0.5), getmean(0.6),
       getmean(tra[7]), getmean(0.8), getmean(0.9))


######## Repeat for Longo ###########

dr <-dat[dat$Location == 'Longo' & dat$day=="D7",]$Normal
head(dr)
table(dr)
dr1 <- dr[dr>0] #& dr$X0<40
hist(dr1)
str(dr1)
len <- length(dr1)
mean(dr1)

sum(dr1>99)
dr2 <- c(dr1,rep(100,3))

#histogram
longo <- ggplot(data = data.frame("x"=dr2 ), aes(x=x)) + geom_histogram(bins=100, fill = 'orange2') + themeJDC +
  ylab("Mosquitoes Collected") + xlab("No. of oocysts counted") + ylim(0,95) +
  scale_x_continuous(breaks = c(1,25,50,75,100), labels = c("1","25","50","75",">100"), limits = c(0.5,101.8))
longo


rp <- 1000
tra <- seq(0.1,0.9,0.1)
block <- rep(0,rp*length(tra))
tra2 <- rep(0,rp*length(tra))
count2 <- 1
for(i in 1:length(tra)){
  for(k in 1:rp){
    count <- 0
    for(j in 1:len){
      aux <- rbinom(n=1,dr1[j],prob=1-tra[i])
      if(aux==0)
        count <- count + 1
    }
    #print(count/len)
    block[count2] <- count/len
    tra2[count2] <- tra[i]
    count2 <- count2 + 1
  }
}
dffL <- data.frame('tba'=block,tra=tra2)
getmean <- function(trax){
  mean(dffL$tba[dff$tra==trax])
}
vL <- c(getmean(0.1), getmean(0.2), getmean(tra[3]), getmean(0.4), getmean(0.5), getmean(0.6),
       getmean(tra[7]), getmean(0.8), getmean(0.9))

fd <- 1 #Change to 100 for a %
pan2 <- ggplot() + geom_point(data = dff,aes(x=1*(tra+0.0085),y=tba*fd,color = 'Klesso'),alpha=.3, size = 1.7, shape = 16) +
  geom_point(data = dffL,aes(x=1*(tra-0.0085),y=tba*fd,color = 'Longo'),alpha=.3, size = 1.7, shape = 16) +
  themeJDC + geom_point(data = data.frame(x=1*(tra+0.0085),y=v), aes(x=x,y=y*fd),color = 'black', size = 1.9, shape = 17) +
  geom_point(data = data.frame(x=1*(tra-0.0085),y=vL), aes(x=x,y=y*fd),color = 'black', size = 1.9, shape = 15) +
  xlab("Transmission Reducing Activity") + ylab("Transmission Blocking Activity") + ylim(0,1) +
  scale_color_manual("Site of\nMosquito\nCollection", values = c('chartreuse3','orange2')) +
  scale_x_continuous(breaks = c(0.0,.25,.50,.75,1), limits = c(0,1)) + theme(legend.position = 'none')
pan2
#Make a better legend.
fake <- ggplot() + geom_point(data = dff,aes(x=1*(tra+0.0076),y=tba*fd,color = 'Klesso'), size = 2.6, shape = 16) +
  geom_point(data = dffL,aes(x=1*(tra-0.0076),y=tba*fd,color = 'Longo'), size = 2.6, shape = 16) +
  themeJDC + geom_point(data = data.frame(x=1*(tra+0.0076),y=v), aes(x=x,y=y*fd),color = 'black', size = 3, shape = 17) +
  geom_point(data = data.frame(x=1*(tra-0.0076),y=vL), aes(x=x,y=y*fd),color = 'black', size = 3, shape = 15) +
  xlab("Transmission Reducing Activity") + ylab("Transmission Blocking Activity") + ylim(0,1) +
  scale_color_manual("Village", values = c('chartreuse3','orange2')) +
  scale_x_continuous(breaks = c(0.1,.30,.50,.70,.90), limits = c(0,1)) + theme(legend.position = c(0.15,0.8))

leg <- get_legend(fake)
text4a <- ggdraw() + draw_plot(pan2) + draw_plot(leg, x=0.23, y=0.6, width = 0.3,height = 0.3)

#Figure 2
grid.arrange(kl,longo,text4a, nrow = 1)

