
---
  title: "GNG_TUS_Study_1"
author: "nomi"
date: "2023-04-05"
output: html_document
---
  
  ############################merge changes with github then delete
  #### LOAD LIBRARIES
  ```{r}
#library(Rccp)
#library(rlang)
library(tibbletime)
library(purrr)
library(rstanarm)
library(stringr)
library(report)
#library(ggplot)
library(ggplot2)
library(dplyr)
#library(lmer)
library(lmerTest)
library(BayesFactor)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(jtools) #for explanation of lmer
library(yarr)
library (yarrr)

#install.packages("devtools")
#devtools::install_github("mikabr/ggpirate")
library(ggpirate)
#p_load ( rms, htmltools, bridgesampling, rstanarm, magrittr, scales, RColorBrewer, 
# parsnip, zoo, wesanderson, yarr, BiocManager, afex, GGally, lmerTest, aod, splines)
```


```{r}
theme_APA <- theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),                                        
        text = element_text(
          size = 20,
          face="bold"), 
        axis.text = element_text(
          size = 20,
          face="plain",
          colour="black"),
        legend.title = element_blank(),
        legend.position = 'top',
        legend.direction = "horizontal",
        legend.text = element_text(
          face="plain",
          colour="black",
          size=20),
        strip.text.x = element_text(
          size = 20,
          face = "bold"),
        panel.background = element_rect(
          fill='white',
          colour='white'),
        strip.background = element_rect(
          fill='white',
          colour='white'),
        axis.title.x = element_text(
          margin = margin(t = 10)),
        plot.margin = grid::unit(c(5, 5, 5, 5), "mm"),
        plot.caption = element_text(size=15)
  )
```

####insert .txt files
```{r}
a.sham<- tibble(filename = list.files("1.SHAM", "*.txt", full.names=TRUE)) %>% 
  group_by(filename) %>% 
  do(read.table(.$filename))
c.ai<- tibble(filename = list.files("3.AI", "*.txt", full.names=TRUE)) %>% 
  group_by(filename) %>% 
  do(read.table(.$filename))
b.dacc<- tibble(filename = list.files("2.dACC", "*.txt", full.names=TRUE)) %>% 
  group_by(filename) %>% 
  do(read.table(.$filename))


```

####add condition column

```{r}
a.sham <- a.sham %>%
  mutate (condition="a.sham")
c.ai <- c.ai %>%
  mutate (condition="c.ai")
b.dacc <- b.dacc %>%
  mutate (condition="b.dacc")
```




#### Set variable names
```{r}
a.sham <- a.sham %>% set_names(c("ID", "stim_ID", "GW", "GAL", "NGW", "NGL", "jitter1", 
                                 "jitter2", "RT", "outcomefeed", "condition"))
b.dacc <- b.dacc %>% set_names(c("ID", "stim_ID", "GW", "GAL", "NGW", "NGL", "jitter1", 
                                 "jitter2", "RT", "outcomefeed", "condition"))
c.ai <- c.ai %>% set_names(c("ID", "stim_ID", "GW", "GAL", "NGW", "NGL", "jitter1", 
                             "jitter2", "RT", "outcomefeed", "condition"))
```


####combine datasets
```{r}
tus <-bind_rows(a.sham, c.ai, b.dacc)

```

### add session column
```{r}
# Extract numbers at the end of each level and create a new column
tus<-tus %>%
  mutate(session = as.factor(str_extract(ID,"(?<=\\D)\\d+(?=_log.txt)")))
View(tus)
```


#### Change name of ID rows to match the subj_ID (i.e.  "pilot data/DBGE0324_log.txt" --> "DBGE0324")
```{r}
tus <- tus %>%
  mutate_at("ID", str_replace, "1.SHAM/", "")# & "log.txt", "")
tus <- tus %>%
  mutate_at("ID", str_replace, "3.AI/", "")
tus <- tus %>%
  mutate_at("ID", str_replace, "2.dACC/", "")
tus <- tus %>%
  mutate_at("ID", str_replace, "_SHAM_2", "")
tus <- tus %>%
  mutate_at("ID", str_replace, "_SHAM_3", "")
tus <- tus %>%
  mutate_at("ID", str_replace, "_SHAM_1", "")
#tus <- tus %>%
#   mutate_at("ID", str_replace, "_SHAM__1", "_sham")
tus <- tus %>%
  mutate_at("ID", str_replace, "_dACC_2", "")
tus <- tus %>%
  mutate_at("ID", str_replace, "_dACC_3", "")
tus <- tus %>%
  mutate_at("ID", str_replace, "_dACC_1", "")
#tus <- tus %>%
#  mutate_at("ID", str_replace, "_ACC_2", "_dACC")
#tus <- tus %>%
#      mutate_at("ID", str_replace, "_ACC_1", "_dACC")
tus <- tus %>%
  mutate_at("ID", str_replace, "_AI_2", "")
tus <- tus %>%
  mutate_at("ID", str_replace, "_AI_3", "")
tus <- tus %>%
  mutate_at("ID", str_replace, "_AI_1", "")
tus <- tus %>%
  mutate_at("ID", str_replace, "_log.txt", "")
#tus <- tus %>%
#    mutate_at("ID", str_replace, "SINB_0180", "SINB0180")
#copied
#tus <- tus %>%
#   mutate_at("ID", str_replace, "_Copy", "")
```



#### Add response col (press - no_press)
```{r}
tus$response <- with(tus, ifelse(RT>0, "press", "no_press"))
tus$GoResponse <- with(tus, ifelse(RT>0, 1, 0)) 
```

### Add Cue column (if stim_ID=NGW then NGW, etc. for  NGL, GW, GAL)
```{r}
tus$Cue <- ifelse(tus$stim_ID == tus$NGL, 'NGL',
                  ifelse(tus$stim_ID == tus$NGW, 'NGW', 
                         ifelse(tus$stim_ID == tus$GW, 'GW', 
                                ifelse(tus$stim_ID == tus$GAL, 'GAL', NA))))
```

#### Add GNG col (Go-noGO)
```{r}
tus$req_action <- ifelse(tus$stim_ID == tus$NGL, 'noGo',
                         ifelse(tus$stim_ID == tus$NGW, 'noGo', 'Go'))   
#NB's way: 
#pilot$ReqAction <- ifelse(pilot$Cue == 'NGL', 'NoGo',
#       ifelse(pilot$Cue == 'NGW', 'NoGo',
#      ifelse(pilot$Cue == 'GW', 'Go',
#      felse(pilot$Cue == 'GAL', 'Go', NA))))

```

#### Add correct column (0-1)  
```{r}
tus <- tus %>%
  mutate(correct = ifelse(req_action == "noGo", ifelse(response == "no_press", 1, 0), 
                          
                          ifelse(req_action == "Go", ifelse(response == "press", 1, 0), 1)))
```

###Trial column 
```{r}
#1,2,3...100...320
tus <- tus %>%
  group_by(ID, condition) %>%
  mutate(trial_number = row_number())
```


#block colums 1, 2, 3, 4   =====sos with the trials that have 400 instead of 320
```{r}
tus$block <- ifelse(tus$trial_number == 1:80, '1',
                    ifelse(tus$trial_number == 81:160, '2', 
                           ifelse(tus$trial_number == 161:240, '3', 
                                  ifelse(tus$trial_number == 241:320, '4', NA))))
# counts trials of the same condition within each block
count_trials <- rollify(function(x) sum(last(x) == x), window=80)#turns a function 
#into a rolling version of itself for use inside of a call to dplyr::mutate() ,
# however it works equally as well when called from purrr::map()
```

#trial_count
```{r}
tus <- tus %>%
  group_by(ID, block, condition) %>%
  mutate(TrialCount = count_trials(Cue)) %>%
  group_by(ID, block,condition, Cue) %>%
  mutate(TrialCount = ifelse(is.na(TrialCount), row_number(ID), TrialCount))
```


### Add win_lose column (win - lose) 
```{r}
tus$feedback <- recode_factor(tus$outcomefeed, "-1" = "lose", 
                              "1" = "win", "0"= "neutral")

```

### Add Win/Avoid condition column (if stim_ID=NGW then NGW, etc. for  NGL, GW, GAL)
```{r}
tus$OutValence <- ifelse(tus$Cue == 'NGL', 'Avoid',
                         ifelse(tus$Cue == 'NGW', 'Win',
                                ifelse(tus$Cue == 'GW', 'Win',
                                       ifelse(tus$Cue == 'GAL', 'Avoid', NA))))
```

#### Add salient feedback column only
```{r}
tus <- tus %>%
  mutate(salient = ifelse(feedback == "win","win",
                          ifelse(feedback == "lose", "lose", NA)))
```


write.csv
```{r}
write.csv(tus, "GNG_TUS_S1.csv")
```


######Processing


#MEAN RT - PRESSES (distribution of RT)
```{r}
tus_presses <- subset(tus, RT != "0") # exlude 0s (no presses)
RT_presses <-tus_presses%>%  group_by(condition) %>% 
  summarise(Mean_RT = mean(RT, na.rm = TRUE))
RT_presses_many_variables <-tus_presses %>%  group_by(condition, req_action, correct, response) %>% 
  summarise(Mean_RT = mean(RT, na.rm = TRUE))
```

#MEAN CORRECT
```{r}
#tus_presses <- subset(tus, RT != "0") # exlude 0s (no presses)
correct_res <-tus%>%  group_by(condition, ID) %>% 
  summarise(correct = mean(correct))
RT_presses_many_variables <-tus_presses %>%  group_by(condition, req_action, correct, response) %>% 
  summarise(Mean_RT = mean(RT, na.rm = TRUE))

```

#MEAN RT - PRESSES (distribution of RT) 
```{r}
tus_presses <- subset(tus, RT != "0") # exlude 0s (no presses)
RT_presses <-tus_presses%>%  group_by(RT, ID) %>% 
  summarise(Mean_RT = mean(RT, na.rm = TRUE));RT_presses #%>%
ggplot(data = RT_presses, 
       aes(x = Mean_RT)) + 
  geom_histogram()+ ggtitle ("Distribution of RT for presses") +theme_APA

# get mean RT per participant
MeanRTpersub <- tus %>% filter(response == "press")%>% group_by(Cue, ID, req_action, OutValence,condition) %>% summarise(RT = mean(RT))

yarrr::pirateplot(formula = RT ~ req_action + OutValence + condition,    # DV = reaction time, IV1 = required action, IV2 = outcome valence
                  data = MeanRTpersub,           
                  theme = 2,
                  # main = "Condition",
                  ylab = "Reaction Time (ms)",
                  ylim = c(400, 1000),
                  bean.f.o = .4, # Bean fill
                  bean.b.o = .2, # Light bean border
                  point.o = .5, # Points (opacity)
                  inf.disp = "line",
                  inf.f.o = 0.8, # Inference fill
                  inf.b.o = 0.8, # Inference border
                  avg.line.o = 0, # Average line
                  point.pch = 21,
                  point.bg = "white",
                  point.cex = .7) + theme_APA
```

#proba of correct per participant#get mean RT per participant - probability of correct instead
```{r}
tus_correct<-tus %>%  group_by(ID, OutValence,condition, req_action) %>% summarise(correct = mean(correct))
yarrr::pirateplot(formula = correct ~ condition  +req_action +OutValence,    
                  data = tus_correct,           
                  theme = 2,
                  #main = "Probability of being correct per condition",
                  ylab = " (p) of correct",
                  ylim = c(0.1, 1.15),
                  bean.f.o = .4, # Bean fill
                  bean.b.o = .2, # Light bean border
                  point.o = .5, # Points (opacity)
                  inf.disp = "line",
                  inf.f.o = 0.8, # Inference fill
                  inf.b.o = 0.8, # Inference border
                  avg.line.o = 0, # Average line
                  point.pch = 21,
                  point.bg = "white",
                  point.cex = .7)
```

#proportion of GO response per condition 
```{r}
# compute proportion of go responses for each participant (first sum per ID, then feed that into general summary)
MeanGoResppersub <- tus %>% group_by(ID, Cue, req_action, OutValence, block, TrialCount) %>% summarise(GoResponse = mean(GoResponse))
MeanGoResppersub <- na.omit(MeanGoResppersub)
MeanGoResp <- MeanGoResppersub %>% group_by(Cue, req_action, block, OutValence, TrialCount) %>% summarise(GoResponse = mean(GoResponse));View(MeanGoResp)

#plot of p(Go) responses
ggplot(MeanGoResp) + 
  geom_smooth(aes(TrialCount, GoResponse,colour=OutValence, lty = req_action)) +
  scale_colour_brewer(palette = "Set1") +
  labs(y = "p(GO)", x = "Trial", title = "Trial-by-trial behaviour") 
+ theme_bw() 
#plot of p(Go) responses per block
ggplot(MeanGoResp) + 
  geom_smooth(aes(TrialCount,GoResponse,colour=OutValence, lty = req_action)) +
  scale_colour_brewer(palette = "Set1") +
  facet_grid(~block)+
  theme_bw() +
  labs(y = "p(GO)", x = "Trial", title = "Trial-by-trial behaviour") 
view(tus)
```


#proportion of GO response per type of Stimulation (sham, ai, dACC) - NADEIDE  and ELSA to check
```{r}

#AI - it does not summarize find why --> I excluded the block. Not sure it did not like it
ai <- subset(tus, condition == "c.ai")

# compute proportion of go responses for each participant
MeanGoResppersub_ai <- ai %>% group_by(ID, Cue, req_action, OutValence, TrialCount) %>% summarise(GoResponse = mean(GoResponse))

MeanGoResp.ai <- MeanGoResppersub_ai %>% group_by(Cue, req_action, OutValence, TrialCount) %>% summarise(GoResponse = mean(GoResponse))

#plot
ai_rolavg<- ggplot(MeanGoResp.ai) + 
  geom_smooth(aes(TrialCount, GoResponse,colour=OutValence, lty = req_action)) +
  scale_colour_brewer(palette = "Set1") +
  labs(y = "p(GO)", x = "Trial", title = "Trial-by-trial behaviour for AI") +
  xlim(NA, 20)+
  theme_bw() 


#dACC - 
dacc <- subset(tus, condition == "b.dacc")
#dacc <- na.omit(dacc)

# compute proportion of go responses for each participant
MeanGoResppersub_dacc <- dacc %>% group_by(ID, Cue, req_action, OutValence, TrialCount) %>% summarise(GoResponse = mean(GoResponse))
MeanGoResp.dacc <- MeanGoResppersub_dacc %>% group_by(Cue, req_action, OutValence, TrialCount ) %>% summarise(GoResponse = mean(GoResponse))

#plot
dacc_rolavg<- ggplot(MeanGoResp.dacc) + 
  geom_smooth(aes(TrialCount, GoResponse,colour=OutValence, lty = req_action)) +
  scale_colour_brewer(palette = "Set1") +
  labs(y = "p(GO)", x = "Trial", title = "Trial-by-trial behaviour for dACC")+
  xlim(NA, 20)+ theme_bw() 



#Sham
sham <- subset(tus, condition == "a.sham")

# compute proportion of go responses for each participant
MeanGoResppersub_sham <- sham %>% group_by(ID, Cue, req_action, OutValence, TrialCount) %>% summarise(GoResponse = mean(GoResponse))

MeanGoResp.sham <- MeanGoResppersub_sham %>% group_by(Cue, req_action, OutValence, TrialCount) %>% summarise(GoResponse = mean(GoResponse))

#plot
sham_rolavg<- ggplot(MeanGoResp.sham) + 
  geom_smooth(aes(TrialCount, GoResponse,colour=OutValence, lty = req_action)) +
  scale_colour_brewer(palette = "Set1") +
  labs(y = "p(GO)", x = "Trial", title = "Trial-by-trial behaviour for sham") +
  xlim(NA, 20)+
  theme_bw() 

#roll avg per CONDITION
rol_avg <- ggarrange(sham_rolavg, ai_rolavg, dacc_rolavg,
                     #labels = c("GW", "GAL", "NGW", "NGL"), 
                     ncol = 2, nrow = 2);rol_avg
```

#rolling average per CUE for all three conditions (plot with NGW for three conditions (sham, ai, dACC), plot GW..etc.)

***it looks like stimulation to ACC increases the propensity to press in the go condition. And stim to both areas slows down learning in the nogo condition***
  ```{r}
#Grey area is CI

#GW
GW <- subset(tus, Cue == "GW")

# compute proportion of go responses for each participant
MeanGoResppersub_gw <- GW %>% group_by(ID, req_action, OutValence, TrialCount, condition) %>% summarise(GoResponse = mean(GoResponse))
MeanGoResppersub_gw <- na.omit(MeanGoResppersub_gw)
MeanGoResp_gw <- MeanGoResppersub_gw %>% group_by(req_action, OutValence, TrialCount, condition) %>% summarise(GoResponse = mean(GoResponse))

#plot######################play with lines
GW<-ggplot(data=MeanGoResp_gw, mapping =aes(TrialCount, GoResponse, colour=condition)) + 
  geom_smooth(linetype = "longdash") +
  scale_colour_brewer(palette = "Set1") +
  #facet_grid(~condition)+
  xlim(NA, 20)+
  ylim(0.05, 0.9)+
  labs(y = "p(GO)", x = "Trial", title = "Trial-by-trial behaviour_GW") +theme_bw() 

#GAL
GAL <- subset(tus, Cue == "GAL")

# compute proportion of go responses for each participant
MeanGoResppersub_gal <- GAL %>% group_by(ID, TrialCount, condition) %>% summarise(GoResponse = mean(GoResponse))
MeanGoResppersub_gal <- na.omit(MeanGoResppersub_gal)
MeanGoResp_gal <- MeanGoResppersub_gal %>% group_by(TrialCount, condition) %>% summarise(GoResponse = mean(GoResponse))

#plot
GAL<-ggplot(data=MeanGoResp_gal, mapping =aes(TrialCount, GoResponse, colour=condition)) + 
  geom_smooth(linetype = "longdash") +
  scale_colour_brewer(palette = "Set1") +
  #facet_grid(~condition)+
  xlim(NA, 20)+
  ylim(0.05, 0.9)+
  labs(y = "p(GO)", x = "Trial", title = "Trial-by-trial behaviour_GAL") +theme_bw() 

#NGW
NGW <- subset(tus, Cue == "NGW")

# compute proportion of go responses for each participant
MeanGoResppersub_NGW <- NGW %>% group_by(ID, TrialCount, condition) %>% summarise(GoResponse = mean(GoResponse))
MeanGoResppersub_NGW <- na.omit(MeanGoResppersub_NGW)
MeanGoResp_NGW <- MeanGoResppersub_NGW %>% group_by(TrialCount, condition) %>% summarise(GoResponse = mean(GoResponse))

#plot - need to work on linetype
NGW<-ggplot(MeanGoResp_NGW) + 
  geom_smooth(aes(TrialCount, GoResponse, colour=condition )) +
  scale_colour_brewer(palette = "Set1") +
  # facet_grid(~condition)+
  xlim(NA, 20)+
  ylim(0.05, 0.9)+
  labs(y = "p(GO)", x = "Trial", title = "Trial-by-trial behaviour_NGW") +theme_bw() 



#NGW
NGL <- subset(tus, Cue == "NGL")

# compute proportion of go responses for each participant
MeanGoResppersub_NGL <- NGL %>% group_by(ID, TrialCount, condition) %>% summarise(GoResponse = mean(GoResponse))
MeanGoResppersub_NGL <- na.omit(MeanGoResppersub_NGL)
MeanGoResp_NGL <- MeanGoResppersub_NGL %>% group_by(TrialCount, condition) %>% summarise(GoResponse = mean(GoResponse))

#plot
NGL<-ggplot(MeanGoResp_NGL) + 
  geom_smooth(aes(TrialCount, GoResponse, colour=condition )) +
  scale_colour_brewer(palette = "Set1") +
  # facet_grid(~condition)+
  xlim(NA, 20)+
  ylim(0.05, 0.9)+
  labs(y = "p(GO)", x = "Trial", title = "Trial-by-trial behaviour_NGL") + theme_bw() 



#CUES PER CONDITION
cues <- ggarrange(GW, GAL, NGW, NGL,
                  #labels = c("GW", "GAL", "NGW", "NGL"), 
                  ncol = 2, nrow = 2);cues
#ylim(0.05, 0.9); 


```


#proportion of correct per condition ***IMPORTANT FOR THE GW FOR THE AI***
```{r}
# compute proportion of go responses for each participant
tus_correct_sub <- tus %>% group_by(ID, OutValence, condition, req_action) %>% summarise(correct = mean(correct))
# compute mean proportion of go responses 
tus_correct <- tus_correct_sub %>% group_by(OutValence, req_action, condition) %>% summarise(correct = mean(correct))


pirateplot(formula = correct ~ condition  + req_action +OutValence,
           data = tus_correct,
           theme = 1,
           bar.f.o = .5, 
           point.o = .5, 
           ylim = c(0.1, 1.0))

#breaks = seq(from = 0, to = 1, by = 0.1))

#THE SAME AS ABOVE BUT WITH VARIANCE/ERROR BARS
pirateplot(formula = correct ~ condition  + req_action +OutValence,
           data = tus_correct,
           theme = 0,
           #main = "Fully customized pirateplot",
           #pal = "southpark", # southpark color palette
           bean.f.o = .6, # Bean fill
           point.o = .3, # Points
           inf.f.o = .7, # Inference fill
           inf.b.o = .8, # Inference border
           avg.line.o = 1, # Average line
           bar.f.o = .5, # Bar
           inf.f.col = "white", # Inf fill col
           inf.b.col = "black", # Inf border col
           avg.line.col = "black", # avg line col
           #bar.f.col = gray(.8), # bar filling color
           point.pch = 21,
           point.bg = "white",
           point.col = "black",
           point.cex = .7)+theme_APA
# ylim = c(0.1, 1.0))

#with stat_summary
ggplot(data=tus_correct, aes(x=correct, y=condition, fill=correct))+ stat_summary(fun.y = mean, geom = "bar") + stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.3)



# Calculate summary statistics (mean and standard deviation)

n_ppt <-tus %>% group_by (ID)%>% count ()

MeanGoResppersub_cond <- tus %>% group_by(ID, Cue, req_action, OutValence, block, TrialCount, condition, correct) %>% summarise(GoResponse = mean(GoResponse))
MeanGoResppersub <- na.omit(MeanGoResppersub_conf)
MeanGoResp <- MeanGoResppersub_conf %>% group_by(Cue, req_action, block, OutValence, TrialCount) %>% summarise(GoResponse = mean(GoResponse))

summary_data <- MeanGoResppersub_cond %>%
  group_by(condition, req_action, OutValence) %>%
  summarise(
    mean_correct = mean(correct),
    sd_correct = sd(correct),
    count = n(),
    se_correct =sd_correct/sqrt(20),
  ) %>%
  na.omit()




# Create a ggplot with error bars, bars, and facets
ggplot(summary_data, aes(x = condition, y = mean_correct, ymin = mean_correct - se_correct, ymax = mean_correct + se_correct, fill = condition)) +
  geom_errorbar(position = position_dodge(width = 0.75), width = 0.2) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.5, alpha = 0.5) +
  geom_text(aes(label = round(mean_correct, 2)), vjust = -0.5, position = position_dodge(width = 0.75)) +
  facet_grid(req_action ~ OutValence) +
  theme_APA

# add points !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! NB


```




#presses for correct  and presses for incorrect (nogos) - ***SEE COMMENTS FOR FURTHER EXPLANATION*** ALSO
****IMPORTANT FOR THE GW FOR THE AI****
  ```{r}
#means - PRESSES FOR GOs (that is pressed when they supposed to be pressing or else pressed for the Go conditions)

presses_for_correct <- subset(tus, response != "no_press"& correct == "1") # that is for both Goes (GW and GAL)
press_mean_ID<-presses_for_correct %>% group_by(ID, OutValence, condition, req_action) %>% count(correct = "1")%>%summarize (n=mean(n))
press_mean_cr  <- press_mean_ID %>% group_by(OutValence, condition, req_action)%>%summarize (n=mean(n))

pirateplot(formula = n ~ condition  + req_action + OutValence,
           main = "presses for GAL(left) and GW (right)",
           data = press_mean_cr,
           theme = 2,
           bar.f.o = .5,
           ylim = c(0, 80),
           bean.f.o = .4, # Bean fill
           bean.b.o = .2, # Light bean border
           point.o = .5, # Points (opacity)
           inf.disp = "line",
           inf.f.o = 0.5, # Inference fill
           inf.b.o = 0.5, # Inference border
           avg.line.o = 0, # Average line
           point.pch = 21,
           point.bg = "black",
           point.cex = 1) #For GW, when AI was stimulated. participants pressed 


#presses for GW (as the one in previous chunk)
presses <- subset(tus,   Cue == "GW"& RT>0) 
press_mean_ID<-presses %>% group_by(ID, OutValence, condition, req_action) %>% count(req_action = "Go",  Cue = "GW" )%>%summarize (n=mean(n))
press_mean  <- press_mean_ID %>% group_by(OutValence, condition, req_action)%>%summarize (n=mean(n))

pirateplot(formula = n ~ condition  + req_action + OutValence,
           main = "presses",
           data = press_mean,
           theme = 2,
           bar.f.o = .5,
           ylim = c(0, 100),
           bean.f.o = .4, # Bean fill
           bean.b.o = .2, # Light bean border
           point.o = .5, # Points (opacity)
           inf.disp = "line",
           inf.f.o = 0.5, # Inference fill
           inf.b.o = 0.5, # Inference border
           avg.line.o = 0, # Average line
           point.pch = 21,
           point.bg = "black",
           point.cex = 1) #same as right graph above



#THE SAME AS ABOVE BUT WITH VARIANCE/ERROR BARS
pirateplot(formula = correct ~ condition  + req_action +OutValence,
           data = tus_correct,
           theme = 0,
           #main = "Fully customized pirateplot",
           #pal = "southpark", # southpark color palette
           bean.f.o = .6, # Bean fill
           point.o = .3, # Points
           inf.f.o = .7, # Inference fill
           inf.b.o = .8, # Inference border
           avg.line.o = 1, # Average line
           bar.f.o = .5, # Bar
           inf.f.col = "white", # Inf fill col
           inf.b.col = "black", # Inf border col
           avg.line.col = "black", # avg line col
           #bar.f.col = gray(.8), # bar filling color
           point.pch = 21,
           point.bg = "white",
           point.col = "black",
           point.cex = .7)+theme_APA
# ylim = c(0.1, 1.0))

``` 


#means - PRESSES FOR NOGOS - ERRORS
```{r}
presses_for_error <- subset(tus, response == "press" & correct == "0") 
press_mean_ID<-presses_for_error %>% group_by(ID, OutValence, condition, req_action) %>% count(correct )%>%summarize (n=mean(n))
press_mean_er  <- press_mean_ID %>% group_by(OutValence, condition, req_action)%>%summarize (n=mean(n))

pirateplot(formula = n ~ condition  + req_action +OutValence,
           main = "presses for noGos",
           data = press_mean_er,
           theme = 2,
           bar.f.o = .5,
           ylim = c(0, 80),
           bean.f.o = .4, # Bean fill
           bean.b.o = .2, # Light bean border
           point.o = .5, # Points (opacity)
           inf.disp = "line",
           inf.f.o = 0.5, # Inference fill
           inf.b.o = 0.5, # Inference border
           avg.line.o = 0, # Average line
           point.pch = 21,
           point.bg = "black",
           point.cex = 1 ) 

```

##performance per participant - Just numbers
```{r}
#summary - correct Cues
ggplot(tus,
       aes(x=Cue, y= correct, color = Cue, fill= Cue)) + 
  geom_bar(stat = "identity")

#correct counts per Cue
ccc<-tus %>% group_by (Cue)%>% count (correct);View (ccc)

#correct counts per Cue per ID
ccc_id<-tus %>% group_by (Cue, ID)%>% count (correct);View (ccc_id)

#correct per ID
c<-tus %>% filter(correct == "1")%>% group_by (ID)%>% count (correct);c

#response per ID
response_per_ppt <-tus %>% group_by (ID)%>% count (response);View (response_per_ppt)

mean(c$n) #avg of 294 correct responses--> 
#may be if a person has less than 230 or 210? correct responses --> CUT? or below 53? see below aggregate
```


#PROPORTIONs - correct PER PARTCIPANT
```{r}
#the proportion of 1s by participant, + more i.e. sound, and language.
#Because the proportion of 1s in a vector with only 0s and 1s is just the mean, this should work:

df <- aggregate(data=tus, correct ~ ID, FUN="mean")#fun = mean #https://r-coder.com/aggregate-r/
p <- ggplot(df)
p <- p + geom_bar(mapping=aes(x=correct, y=ID), stat='identity')+
  geom_vline(xintercept = 0.53, col = "red")+ 
  geom_vline(xintercept = 0.50, col = "grey")+
  geom_vline(xintercept = 0.80, col = "green")+
  geom_text(aes(x=0.54, label="53%", y=4), colour="black", angle=90) +
  geom_text(aes(x=0.49, label="50%", y=9), colour="black", angle=90)+
  geom_text(aes(x=0.79, label="80%", y=6), colour="black", angle=90) 
p + ylab('ID')

#same as above but with stat_summary
ggplot(tus) + 
  stat_summary(aes(x = ID, y = correct), 
               geom = "bar")

#performance per participant and condition
aggregate(data=tus, correct ~ ID+condition, FUN="mean") #fun = mean #https://r-coder.com/aggregate-r/

df_c <- aggregate(data=tus, correct ~ ID +condition, FUN="mean")

ggplot(data = df_c) + 
  geom_bar(mapping = aes(fill = condition, x = ID, y = correct), stat = "identity", position = "dodge")

#same but with lines and dots
df_c %>%ggplot(aes(condition, correct, color=ID)) +
  #geom_boxplot()+
  geom_point(aes(fill=ID),size=5) +
  # scale_x_log10()+
  geom_line(aes(group = ID),color="grey")
#ggsave("#performance per participant and condition.png")


#performance per  condition
aggregate(data=tus, correct ~ condition, FUN="mean") #fun = mean #https://r-coder.com/aggregate-r/

df_con <- aggregate(data=tus, correct ~ condition, FUN="mean")

ggplot(data = df_con) + 
  geom_bar(mapping = aes(fill = condition, x = condition, y = correct), stat = "identity", position = "dodge")
```

RTs for participants 
```{r}
#rt per participant per condition
aggregate(data=tus, RT ~ ID+condition, FUN="mean") #fun = mean #https://r-coder.com/aggregate-r/

df_rt <- aggregate(data=tus, RT ~ ID +condition, FUN="mean")

ggplot(data = df_rt) + 
  geom_bar(mapping = aes(fill = condition, x = ID, y = RT), stat = "identity", position = "dodge") 


#same but with lines and dots
df_rt %>%ggplot(aes(condition, RT, color=ID)) +
  #geom_boxplot()+
  geom_point(aes(fill=ID),size=5) +
  #scale_x_log10()+
  geom_line(aes(group = ID),color="grey")

#-------------------------------

#for NGW only (correct) - NO POINT
ngwin <- subset(tus, Cue == "NGW")

ngwin_cr <- aggregate(data=ngwin, correct ~ ID +condition, FUN="mean")

ngwin_cr %>%ggplot(aes(condition, correct, color=ID)) +
  #geom_boxplot()+
  geom_point(aes(fill=ID),size=5) +
  #scale_x_log10()+
  geom_line(aes(group = ID),color="grey")


#for NGW only (RT) - NO POINT

ngwin_rt <- aggregate(data=ngwin, RT ~ ID +condition, FUN="mean")

ngwin_rt %>%ggplot(aes(condition, RT, color=ID)) +
  #geom_boxplot()+
  geom_point(aes(fill=ID),size=5) +
  #scale_x_log10()+
  geom_line(aes(group = ID),color="grey")
```




Analysis

#Accuracy

#NB RUN A BASIC ANOVA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
```{r}
model_acc <- glmer(correct ~ req_action  + trial_number +OutValence + req_action*condition +(1|ID) , 
                   data  = tus, family="binomial", 
                   glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) #REML = FALSE
summary (model_acc) 
summ(model_acc, exp = T)# set "exp = T" to show esponentiated estimates; if you need standardised estimaets, set "scale = T"
summ(model_acc, scale = T)
report(model_acc)
report(model_acc) %>% summary()
report_table(model_acc)
report_performance(model_acc)
report_parameters(model_acc)


lm(correct ~ condition, data=tus)


#Cue impact
model_acc_cue <- glmer(correct ~ Cue*condition +(1|ID) , 
                       data  = tus, family="binomial", 
                       glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) #REML = FALSE
summary (model_acc_cue) 

#visualize how spread accuracy is 
ggplot(tus_correct_sub, aes(x = correct , y = req_action, color=req_action))+
  geom_jitter(stat = "identity")+
  facet_grid(~condition)
geom_smooth(method= "lm", se = TRUE)


#BF
tus$ID <- as.factor(tus$ID)
tus$req_action <- as.factor(tus$req_action)
tus$condition <- as.factor(tus$condition)
tus$trial_number <- as.factor(tus$trial_number)
tus$OutValence <- as.factor(tus$OutValence)

#drop NAs
df <- na.omit(tus)
df$ID <- as.factor(df$ID)
df$correct <- as.integer(df$correct)
df$req_action <- as.factor(df$req_action)
df$condition <- as.factor(df$condition)
df$trial_number <- as.factor(df$trial_number)
df$OutValence <- as.factor(df$OutValence)
#run anovaBF (not the correct analysis as the DV is categorical - just trying)
accu <-anovaBF(correct ~ req_action*condition  + ID, data = df, whichRandom = "ID",
               progress=FALSE)
accu[4]/accu[3]


#banova!!!!!!!! - not working atm
# Making sure JAGS is installed
# It can be downloaded here http://mcmc-jags.sourceforge.net (the size is only 2.3 MB)
library(BANOVA)
library (rstan)

tus_acc <-  BANOVA.run (correct ~ req_action, ~condition*req_action, data=df, model_name = 'Binomial', id = 'ID',
                        as.integer(16), num_trials = as.integer(16), iter = 100, thin = 5, chains = 2)

```

#RTs

```{r}
model_rt <- lmer(RT ~ req_action  + OutValence + req_action*condition +(1|ID), # + trial_number (excluded cause it gives me the outcome of all trials, it does not group them)
                 data  = tus, REML = FALSE) #REML = FALSE
summary (model_rt) 
#print(model_rt, correlation=FALSE)
report_table(model_rt)




#run anovaBF 
rt <-anovaBF(RT ~ req_action*condition  + ID, data = df, whichRandom = "ID",
             progress=FALSE)
accu[4]/accu[3]
summary(rt)

#Cue impact
model_rt_cue <- lmer(RT ~ Cue*condition +(1|ID) , 
                     data  = tus,REML = FALSE) #REML = FALSE
anova(model_rt)
summary (model_rt_cue) 

```




#SALIENT ONLY#
```{r}
#select salient only (win, lose) and keep variables ID through block and all columns between them
tus_salient <- subset(tus, salient=!"neutral",
                      select=ID:salient);  

salient_df <-tus_salient %>% drop_na();  View(salient_df) # WRONG IT DRPOPS ACC AS WELL

```  


#Accuracy
```{r}
model_acc_sal <- glmer(correct ~ req_action  + trial_number + OutValence*condition + req_action*condition +(1|ID) , 
                       data  = salient_df, family="binomial", 
                       glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) #REML = FALS
summary(model_acc_sal)
#Cue impact
sal_cue_accu <- glmer(correct ~ Cue*condition +(1|ID) , 
                      data  = salient_df, family="binomial", 
                      glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) #REML = FALSE
summary (sal_cue_accu) 

```

#RTs

```{r}
model_rt_sal <- lmer(RT ~ req_action  + trial_number + OutValence*condition + req_action*condition+(1|ID), 
                     data  = salient_df, REML = FALSE) #REML = FALSE
summary (model_rt_sal) 
report_table(model_rt_sal)
#Cue impact
sal_cue_rt <- lmer(RT ~ Cue*condition +(1|ID), 
                   data  = salient_df, REML = FALSE) #REML = FALSE
summary (sal_cue_rt) 
```