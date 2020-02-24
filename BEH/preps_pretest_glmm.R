library(ggplot2)
library(RColorBrewer)
library(data.table)
library(lme4)
library(dplyr)
library(sjPlot)
library(sjmisc)
library(plotly)
pastel_colors = brewer.pal(8, "Pastel1")
set_colors = brewer.pal(8, "Set1")

#load in previously preprocessed data
load('pretestdata_clean.RData')

#add type information to the full response matrix
df.responses_clean <- df.responses_clean %>% mutate(response = response_attachment=="Nomen")
df.responses_clean <- df.responses_clean %>% mutate(attach = attachment=="Noun")
mydata <- merge(df.responses_clean,df.stim_acc,by.x="items",by.y="fullID")
mydata <- merge(df.responses_clean,df.stim_acc,by.x="items",by.y="fullID")
mydata <- mydata[,c(1,2,3,4,5,7,8,9,10,11)]

#exclude missing data (not all items were seen by all subjects)
mydata <- mydata[complete.cases(mydata),]

#format data columns to be either factors or numeric
mydata$attach <- as.factor(as.numeric(mydata$attach))
mydata$response <- as.factor(as.numeric(mydata$response))
mydata$subjects <- as.factor(mydata$subjects)
mydata$rating_plausibility <- as.numeric(mydata$rating_plausibility)
mydata$hits <- as.factor(mydata$hits)
mydata$attachment <- as.factor(mydata$attachment)
mydata$response_attachment <- as.factor(mydata$response_attachment)
contrasts(mydata$response_attachment) = c(-1,1)
contrasts(mydata$attachment) <- c(-1,1) 
contrasts(mydata$hits) <- c(-1,1) 
mydata$rt_attachment <- as.numeric(mydata$rt_attachment)

mydata %>% group_by(hits) %>% summarise(n())

x <- df.responses_clean %>% group_by(items) %>% summarise(mean = mean(na.omit(rating_plausibility)))
hist(x$mean)
## How different are Type 1 and Type 2 items?
summary_acc <- mydata %>%
  group_by(items,Type,attachment,hits) %>%
  summarise(n=n()) %>% mutate(acc = n / sum(n))

p <- ggplot(data = summary_acc[summary_acc$hits==1,], aes(x = factor(Type), y = acc,item = items,fill=factor(Type))) +
  geom_boxplot() +
  facet_wrap(~attachment) +
  ggtitle("Accuracy comparison between manipulations")+
  labs(y="mean accuracy", x="Manipulation Type")+
  scale_x_discrete(labels = c("Type 1","Type 2","Type 1","Type 2")) +
  theme(legend.position = "none")
p <- ggplotly(p,tooltip = c("item","mean_acc"))
p

# same for RTs
summary_acc <- mydata %>%
  group_by(items,Type,attachment) %>%
  summarise(mean_rt = mean(rt_attachment)) 

p <- ggplot(data = summary_acc, aes(x = factor(Type), y = mean_rt,item = items,fill=factor(Type))) +
  geom_boxplot() +
  facet_wrap(~attachment) +
  ggtitle("Accuracy comparison between manipulations")+
  labs(y="mean reaction time", x="Manipulation Type")+
  scale_x_discrete(labels = c("Type 1","Type 2","Type 1","Type 2")) +
  theme(legend.position = "none")
p <- ggplotly(p)
p
# same for plausibility
summary_acc <- mydata %>%
  group_by(items,Type,attachment) %>%
  summarise(mean_plaus = mean(rating_plausibility)) 

p <- ggplot(data = summary_acc, aes(x = factor(Type), y = mean_plaus,item = items,fill=factor(Type))) +
  geom_boxplot() +
  facet_wrap(~attachment) +
  ggtitle("Accuracy comparison between manipulations")+
  labs(y="mean plausibility", x="Manipulation Type")+
  scale_x_discrete(labels = c("Type 1","Type 2","Type 1","Type 2")) +
  theme(legend.position = "none")
p <- ggplotly(p)
p

### glmm
# need random intercept for subject due to independence violation (1|subject)
# random intercept for items to account for by-item variation ( some items leading to larger differences between NA VA)
# not enough data to also include random slopes for rt_attachment, interaction and rating_plausibility.

m <- glmer(hits ~ 1 + attachment + rt_attachment + rating_plausibility + Type +
             rt_attachment*attachment + rating_plausibility*attachment + Type*attachment + 
             Type*rt_attachment + rt_attachment*rating_plausibility +  
             Type*rating_plausibility + 
             (1+attachment|subjects) + (1|items),
           data = mydata, family="binomial",control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(m)

#1) ppl had a slight Noun bias in this experiment (they were more often incorrect for VA items)
#2) acc decreased with longer rts (ppl were slower when incorrect)
#3) acc increased with higher plausibility (but mostly for verbs)
## plot interaction term
plot_model(m, type = "pred", terms = c("rating_plausibility", "attachment"))
# higher plausibility leads to higher accuracy. And lower plausibility leads to lower accuracy
# This effect is stronger in the verb condition
plot_model(m, type = "pred", terms = c("attachment", "Type"))
#For Type 1 items, NA sentences were more accurate than VA sentences
#But for Type 2 items, VA sentences were slightly more accurate than NA
plot_model(m, type = "pred", terms = c("rt_attachment [all]","rating_plausibility"))
#
mydata$rating_plausibility <- as.factor(mydata$rating_plausibility)
m2 <- glmer(rating_plausibility ~ rt_attachment + attachment + Type + hits + 
                                  rt_attachment*attachment + rt_attachment*Type + rt_attachment*hits +
                                  attachment*Type + attachment*hits +
                                  Type*hits +
            (1+attachment + rt_attachment|subjects) + (1|items),
            data = mydata,family="binomial",control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(m2)
#unfortunately cannot include random slope for both type and attachment due to singularity
# hits were rated more plausible compared to misses
# Verb attachments were rated as more plausible than noun attachments
# Sentences from Type 1 manipulations were rated as more plausible than those from type 2 
plot_model(m2,type="pred",terms=c("hits","Type","attachment"))

#check for singularity: "The definition of singularity is that some of the constrained parameters of the random effects theta parameters are on the boundary (equal to zero, or very very close to zero, say <10−6):"
tt <- getME(m,"theta")
ll <- getME(m,"lower")
min(tt[ll==0])


