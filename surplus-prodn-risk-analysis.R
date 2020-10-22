# Surplus production model and risk analysis 
#First fits surplus production model to CPUE series
# estimates probability of overfishing if there
# is a future productivity decline 
# Coral trout 
#CJ Brown 2020-09-26

#The model is from Meyer and Millar 2000 (Applied Statistics)
# and Winker et al. 2018 (Fisheries Research)
# But modified to use Pella thomlinson production.
# Fits a surplus production model to CPUE and catch timeseries
# Using the Stan program. 
# We allow q (catchability) to vary before and after the 
#2004 restructure and rezoning on the GBR. 
# m in the PT model set to 1.19 as per meta-analysis of 
# Thorson et al. 2012 CJFAS. This gives BMSY/B0 of 0.400
# (BMSY/B0 = m^(-1/(m-1)))
# see: Winker et al. 2018 JABBA

rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(rstan)
library(rethinking)
library(cowplot)

load("data-raw/qfish-dat.rda")

#calculates params for lognormal based on mean and SD
logNormalParams <- function(mB, sigma_B){
  a <-2*log(mB) - 0.5 * log(sigma_B^2 + mB^2) 
  b <- sqrt(-2*log(mB)+log(sigma_B^2 + mB^2))
  return(list(lmean = a, lsigma = b))
}

spkeep <- data.frame(COMMON_NAME = c(
  "coral trout", #plectropomus/variola, rank 1 by total catch 
  "Redthroat Emperor", #Lethrinus miniatus, rank 2
  "Saddletail Snapper", #Lutjanus malabaricus, rank 6
  "Red Emperor", #Lutjanus sabae, rank 7
  "Goldband snapper", #Pristipomoides multidens 
  "RSK prawn"
),
r_est = c(0.43, 0.42, 0.38, 0.28, 0.57, 1.4),
mod_rsq = NA,
mod_rsq_catch = NA)
nspp <- nrow(spkeep)
catchability_increase <- 1.0 #1% per year 
#also find and repleace to save with +0 etc.... 

spp <- unique(dat$COMMON_NAME)
spkeep$catchability_increase <- catchability_increase

#
# Prior params 
#
#Guess B1 and B0
bf_guess <- 0.8 #fraction biomass is of Bzero in 1989
b1_fract <- 50 # Multiple B1 is of catch in 1989
B0_CV <- 0.5 # CV for lognormal on B0 - prior 

sigma_shape <- 5 # shape parameter for process errors 


#Plot prior for mu as multiples of B[t-1]
mu <- (exp(rnorm(100000, mean = 0,
                 # sd = rexp(100000, 14))))
                 sd = rgamma(100000, 6,scale = 0.5))))
hist(mu,30)
quantile(mu)

#
# Prep data 
#

drsim <- NULL
modout <- NULL

for (ispp in 1:nspp){
  
  print(spkeep$COMMON_NAME[ispp])
  #Set up data 
  
  datct <- filter(dat, COMMON_NAME == spkeep$COMMON_NAME[ispp]) %>%
    filter(!is.na(catch)) %>%
    mutate(catch_std = catch/sd(catch), 
           #reducing days so lncpue is postiie helps with fitting 
           days_std = (days/(max(days)*5))*catchability_increase^(year-1989),
           ln_cpue = log(catch_std/days_std)) 
  sdcatch <- sd(datct$catch)
  sddays <- sd(datct$days)
  
  
  g1 <- ggplot(datct) + aes(x = year, y = ln_cpue) + geom_line()+ 
    theme_classic()
  g2 <- ggplot(datct) + aes(x = year, y = catch) + geom_line()+ 
    theme_classic()
  g3 <- ggplot(datct) + aes(x = year, y = days) + geom_line()+ 
    theme_classic()
  gall <- plot_grid(g1, g2, g3, nrow = 1)
  ggsave2(gall, file = 
            paste0("model-outputs/", spkeep$COMMON_NAME[ispp], "_dat+0.jpg"),
          width = 12)
  
  
  r_est <- spkeep$r_est[ispp]

  B1_guess <- datct$catch_std[1]*b1_fract
  B0_guess <- B1_guess*(1/bf_guess) #assume B1 is XX% of B0
  B0params <- logNormalParams(B0_guess, B0_guess * B0_CV)
  
  rparams <- logNormalParams(r_est, r_est*0.1)

  #Setup data 
  datstan <- list(N = nrow(datct),
                  lnC = log(datct$catch_std),
                  lnCPUE = datct$ln_cpue,
                  catches = datct$catch_std,
                  logK_mean = B0params$lmean,
                  logK_sd = B0params$lsigma,
                  init_fraction = bf_guess,
                  sigma_shape = sigma_shape,
                  logr_mean = rparams$lmean,
                  logr_sd = rparams$lsigma
  )
  
  options(mc.cores = 3) #parallel::detectCores())
  
  
  fitm1 <- stan(file = "surplus-prodn-model-cpue.stan", data = datstan,
                iter=5000, chains=3, 
                control = list(max_treedepth = 12), 
                init = list(
                  list(lnr = log(r_est), lnK = log(b1_fract), lnq1 = log(2),
                       lnq2 = log(2),
                       sigma = 0.5, tau = 0.5),
                  list(lnr = log(r_est), lnK = log(b1_fract/2), lnq1 = log(2),
                       lnq2 = log(2),
                       sigma = 0.05, tau = 0.15),
                  list(lnr = log(r_est), lnK = log(b1_fract*2), lnq1 = log(2),
                       lnq2 = log(2),
                       sigma = 0.3, tau = 0.2)
                ))
  
  modout <- c(modout, list(fitm1))
  
  #
  # Relative biomass, relative to BMSY
  #
  
  post <- extract.samples(fitm1) %>% data.frame()
  ib <- grep("\\bBrel[.]", names(post))
  post_expected_bio <- as.matrix(post[,ib])
  
  bio_est <- t(apply((post_expected_bio), 2, quantile, probs = c(0.05, 0.5, 0.975))) %>%
    data.frame() %>% bind_cols(datct)
  
  go1 <- ggplot(bio_est) + 
    geom_hline(yintercept = 0.5, color = "grey20")+ 
    geom_hline(yintercept = 1, color = "grey20")+ 
    aes(x = year, color = NULL) +
    geom_line(aes(y = X50.)) + 
    geom_ribbon(aes(ymin = X5., ymax = X97.5.), alpha = 0.5, color = NA) + 
    theme_classic()  + 
    ylim(0,1)+ 
    theme_classic()
  
  #
  #Plot CPUE over time
  #
  
  ib <- grep("lnCPUE_hat[.]", names(post))
  post_expected_bio <- as.matrix(post[,ib])
  iq <- grep("lnq2", names(post))
  
  cpue_post <- t(apply((post_expected_bio), 2, quantile, probs = c(0.05, 0.5, 0.975))) %>%
    data.frame() %>% bind_cols(datct)
  
  acf(cpue_post$X50. - datct$ln_cpue)
  
  go2 <- ggplot(cpue_post) + 
    aes(x = year, y = ln_cpue, color = NULL) +
    geom_point() + 
    geom_line(aes(y = X50.)) + 
    geom_ribbon(aes(ymin = X5., ymax = X97.5.), alpha = 0.5, color = NA) + 
    theme_classic()+ 
    theme_classic()
  
  go3 <- ggplot(cpue_post) + 
    aes(x = ln_cpue, y = X50., color = NULL) +
    geom_point() + 
    theme_classic()
  
  spkeep$mod_rsq[ispp] <-  cor(cpue_post$ln_cpue, cpue_post$X50.)^2

  gall <- plot_grid(go1, go2, go3, nrow = 1)
  ggsave2(gall, file = 
            paste0("model-outputs/", spkeep$COMMON_NAME[ispp], "_model+0.jpg"),
          width = 12)
  
  
  #
  # 'Risk analysis': Impact of reduction in r
  #
  
  catch_today <- mean(datct$catch_std[(nrow(datct)-3):nrow(datct)])
  effort_today <- mean(datct$days_std[(nrow(datct)-3):nrow(datct)])
  
  # ir <- grep("r", names(post))
  ir <- 38
  post_r <- as.matrix(post[,ir])
  iK <- grep("\\K", names(post))
  iK <- 37
  post_K <- as.matrix(post[,iK])
  iq <- grep("\\lnq2", names(post))
  post_q <- as.matrix(exp(post[iq]))
  
  
  #reduce r incrementally
  rper <- seq(0, 0.95, by = 0.01)
  nrper <- length(rper)
  prob_emsy <- prob_cmsy<- numeric(nrper)
  for (i in 1:nrper){
    
    post_r2 <- post_r *(1-rper[i])
    post_FMSY <-  (post_r2/0.19) * (1 - (1/1.19))
    post_BMSY <- post_K * (1.19^(-1/(0.19)))
    post_cmsy <- post_FMSY * post_BMSY
    post_emsy <-  post_FMSY*post_q
    
    #Calculate probability that catches and effort 
    # are over MSY and EMSY with a given level of productivity decline
    prob_emsy[i] <- sum(post_emsy<effort_today)/nrow(post_emsy)
    prob_cmsy[i] <- sum(post_cmsy<catch_today)/nrow(post_cmsy)
  }
  
  drsim <- c(drsim, list(data.frame(
    rper = rper,
    post_cmsy = prob_cmsy,
    spp = spkeep$COMMON_NAME[ispp]
  )))
  
}

save(spkeep, drsim, modout, file = "model-outputs/model-runs-finfish+0.rda")


drplot <- bind_rows(drsim) 

#Adjust position slightly to avoid overlapping in plots
i <- drplot$spp == "Goldband snapper"
drplot$post_cmsy[i] <- drplot$post_cmsy[i] + 0.02

i <- drplot$spp == "Saddletail Snapper"
drplot$post_cmsy[i] <- drplot$post_cmsy[i] + 0.02

i <- drplot$spp == "Redthroat Emperor"
drplot$post_cmsy[i] <- drplot$post_cmsy[i] - 0.02


gfinal <- ggplot(drplot) + 
  aes(x = rper*100, y = post_cmsy, color = spp) +
  geom_line(size = 1) + 
  theme_classic() + 
  scale_color_brewer(palette = "Dark2") + 
  xlab("Reduction in r (%)") + 
  ylab("Probability \n catch > MSY")

gfinal
ggsave(gfinal, file = "model-outputs/hab-loss-sensitivty_model+0.jpg",
        width = 12)


