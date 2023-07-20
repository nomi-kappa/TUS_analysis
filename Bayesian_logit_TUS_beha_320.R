#Bayesian logistic regression for TUS beha only 320 trials

library(readr)
tus <- read_csv("~/R/R/TUS_analysis/GNG_TUS_S1.csv")
names(tus)

library (tidyverse)
library(bridgesampling)
library(rstanarm)
#For a Bayesian, fit basically the same model using the rstanarm package and then use "bridge sampling" 
#to compare them to get a Bayes Factor 

#I should also actually use posterior predictions from these models to quantify the size of the effect
#(ie. how many more/less errors in sham/ai/acc cond do we expect)
#and to plot the effect, rather than use a Bayes Factor, 
#but I  prefer that as a simpler summary of whether the effect is "really there".

#correct with rstarm + bridging

# install htmltools
install.packages("htmltools", type = "source")
options(mc.cores = 4)

tus_dropped <-tus %>% drop_na();  

#model comparison 1
cor.m1.stan <- rstanarm::stan_glmer(correct ~ condition * OutValence * req_action + (1|ID), 
                                    weights=correct, data  = tus_dropped, family="binomial", chains=4, iter=2e4, 
                                    diagnostic_file = "__tus1.csv")
summary(cor.m1.stan)


cor.m1.0.stan <- rstanarm::stan_glmer(correct ~ OutValence * req_action + (1|ID),
                                      weights = correct, data = tus_dropped, family = "binomial", chains = 4, iter = 2e4,
                                      diagnostic_file = "__tus0.csv")

summary(cor.m1.0.stan)

bridge_1 <-bridgesampling::bridge_sampler(cor.m1.stan)
bridge_0 <- bridgesampling::bridge_sampler(cor.m1.0.stan)

bayes_factor(bridge_1, bridge_0, log = FALSE) # Estimated BF in favor of x1 over x2: 0.00532
  

#model comparison 2 - condition alone
cor.m2.stan <- rstanarm::stan_glmer(correct ~ condition + (1|ID),
                                      weights = correct, data = tus_dropped, family = "binomial", chains = 4, iter = 2e4,
                                      diagnostic_file = "__tus0.csv")

cor.m2.0.stan <- rstanarm::stan_glmer(correct ~   (1|ID),
                                    weights = correct, data = tus_dropped, family = "binomial", chains = 4, iter = 2e4,
                                    diagnostic_file = "__tus2.0.csv")

bridge_2 <-bridgesampling::bridge_sampler(cor.m2.stan) 
bridge_2.0 <- bridgesampling::bridge_sampler(cor.m2.0.stan) 

bayes_factor(bridge_2, bridge_2.0, log = FALSE)#Estimated Bayes factor in favor of x1 over x2: 0.00000


#model comparison 3 - condition req_action
cor.m3.stan <- rstanarm::stan_glmer(correct ~ condition * req_action +  (1|ID),
                                      weights = correct, data = tus_dropped, family = "binomial", chains = 4, iter = 2e4,
                                      diagnostic_file = "__tus3.csv")


cor.m3.0.stan <- rstanarm::stan_glmer(correct ~ condition + req_action  + (1|ID),
                                    weights = correct, data = tus_dropped, family = "binomial", chains = 4, iter = 2e4,
                                    diagnostic_file = "__tus3.0csv")


bridge_3 <-bridgesampling::bridge_sampler(cor.m3.stan) 
bridge_3.0 <- bridgesampling::bridge_sampler(cor.m3.0.stan) 

bayes_factor(bridge_3, bridge_3.0, log = FALSE)

#model comparison 4 - condition cue
cor.m4.stan <- rstanarm::stan_glmer(correct ~ condition * req_action * Cue  (1|ID),
                                    weights = correct, data = tus_dropped, family = "binomial", chains = 4, iter = 2e4,
                                    diagnostic_file = "__tus3.csv")


cor.m4.0.stan <- rstanarm::stan_glmer(correct ~ (condition + Cue) *req_action  + (1|ID),
                                      weights = correct, data = tus_dropped, family = "binomial", chains = 4, iter = 2e4,
                                      diagnostic_file = "__tus3.0csv")


bridge_4 <-bridgesampling::bridge_sampler(cor.m4.stan) 
bridge_4.0 <- bridgesampling::bridge_sampler(cor.m4.0.stan) 

bayes_factor(bridge_4, bridge_4.0, log = FALSE)