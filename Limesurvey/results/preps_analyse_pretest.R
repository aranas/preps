
library(plotly)
library(reshape)
library(RColorBrewer)
library(data.table)
library(dplyr)
#Create color scheme
pastel_colors = brewer.pal(8, "Pastel1")
set_colors = brewer.pal(8, "Set1")
#load stimulus file
stimuli         <- read.csv(file="stimuliB.csv",head=TRUE,sep=",",na.strings = c("","NAN"), stringsAsFactors = FALSE)
stimuli$fullID  <- as.character(interaction(stimuli[,c(13,3)],sep = ""))
#load response file
file_names = list.files(pattern="^results",path = "Bresponses",full.names = TRUE)
temp <- lapply(file_names,read.csv,sep=",",na.strings = c("","NAN"), stringsAsFactors = TRUE)
df.responses <- Reduce(function(x,y) merge(x,y,all=TRUE,sort=TRUE),temp)


subjectnum      <- nrow(df.responses)
#transpose and fix column names and classes
df.responses                <- setDT(df.responses,keep.rownames = TRUE)
colnames(df.responses)[1]   <- "subjects"
df.responses                <- Filter(function(x) !(all(x=="")), df.responses) #delete blank columns

#Split all item specific data from subject specific data
ind_items   <- grep("[VN]A[0-9]",colnames(df.responses),value = TRUE)
df.subjectinfo <- df.responses
df.subjectinfo[,ind_items] <- NULL
df.responses  <- df.responses[,c("subjects",ind_items),with=FALSE]

#reshape items
ind             <- grep("^[VN]A[0-9]*$",colnames(df.responses),value = TRUE)
df.respAttach   <- df.responses[,c("subjects",ind),with=FALSE]
df.respAttach   <- melt(df.respAttach,id="subjects",value.name="response_attachment",variable.name='items')
ind             <- grep("^[VN]A[0-9]*Time$",colnames(df.responses),value = TRUE)
df.respTime     <- df.responses[,c("subjects",ind),with=FALSE]
df.respTime     <- melt(df.respTime,id="subjects",value.name="rt_attachment",variable.name='items')
ind             <- grep("\\.$",colnames(df.responses),value = TRUE)
df.respPAttach  <- df.responses[,c("subjects",ind),with=FALSE]
df.respPAttach  <- melt(df.respPAttach,id="subjects",value.name="rating_plausibility",variable.name='items')
ind             <- grep("P[VN]A[0-9]*Time$",colnames(df.responses),value = TRUE)
df.respPTime    <- df.responses[,c("subjects",ind),with=FALSE]
df.respPTime    <- melt(df.respPTime,id="subjects",value.name="rt_plausibility",variable.name='items')
df.responses    <- cbind(df.respAttach,df.respTime[,3],df.respPAttach[,3],df.respPTime[,3])
df.responses$hits <- as.numeric((grepl("N",df.responses$items) & df.responses$response_attachment == "Nomen") | 
                      (grepl("V",df.responses$items) & df.responses$response_attachment == "Verb"))
df.responses    <- df.responses %>%
                    mutate(attachment = ifelse(grepl("N",items),"Noun","Verb"))
                                                                       
#############ANALYSIS#################
summary(df.subjectinfo$age)
summary((df.subjectinfo$interviewtime)/60)

#find outlier (who did not pass the test?)
#number of wrong answers per participant x/4

df.subjectinfo %>%
  select(subjects,starts_with("correct")) %>%
    mutate(number_correct = rowSums(!is.na(.[,2:ncol(.)]))) %>%
      select(subjects,number_correct)

avg_rts <- df.responses %>%
                        group_by(subjects) %>%
                        summarise(mean_rt = mean(rt_attachment)) 
outlier <- boxplot(avg_arts[,2])$out
outlier_indx <- avg_rts$subjects[which(round(avg_rts$mean_rt,digits=2) %in% round(outlier,digits=2))]
#Remove Outlier
df.responses <- df.responses %>%
                filter(!(subjects %in% outlier_indx))
#average accuracy of unambiguous items
subset_unambiguous <- stimuli$fullID[(stimuli$Unambiguous.==1)]
sentences_unambiguous <- stimuli[(stimuli$Unambiguous.==1),4:12]

df.responses %>%
  filter(items %in% subset_unambiguous) %>%
    group_by(subjects) %>%
      summarise(accuracy = mean(hits))

#Biases: Reaction Times total and per attachment

summary_rts <- df.responses %>%
                 group_by(subjects,attachment,hits) %>%
                   summarise(mean_rt = mean(rt_attachment)) 
                      
p <- ggplot(summary_rts, aes(x = factor(hits), y = mean_rt,subject = subjects,fill=factor(hits))) +
    geom_boxplot() +
    facet_wrap(~attachment) +
    geom_jitter(size = 2) +
    ggtitle("RTs averaged over items")
p <- ggplotly(p,tooltip = c("subject","mean_rt"))
p
# Assess material
mean_accuracy <- df.responses %>%
                    group_by(items) %>%
                      summarise(percent_correct = round(sum(hits)/subjectnum,digits=2)) %>%
                        summarise(mean_correct = mean(percent_correct))
acc_per_item <- df.responses %>%
                  group_by(items) %>%
                    summarise(percent_correct = round(sum(hits)/subjectnum,digits=2))


items_reject  <- df.responses %>%
                  group_by(items) %>%
                    summarise(percent_correct = round(sum(hits)/subjectnum,digits=2)) %>%
                     filter(percent_correct <= 0.4)
stimuli_reject <- stimuli[(stimuli$fullID %in% items_reject$items),c(4:12,15)]
merge(stimuli_reject,items_reject,by.x = "fullID",by.y = "items")

#binomial distribution 
x <- seq(0,20,by = 1)
y <- dbinom(x,20,0.5)
plot(x,y)
#text mit uft8 ausgeben
##Plausibility analysieren
#first extract contingency table
plausibility_counts <- df.responses %>%
                        group_by(attachment) %>%
                        count(rating_plausibility)
plausibility_counts <- cbind(plausibility_counts[plausibility_counts$attachment=="Noun",3],
                             plausibility_counts[plausibility_counts$attachment=="Verb",3])
colnames(plausibility_counts) <- c("Noun","Verb")
#chi-square test
chisq.test(plausibility_counts)
#reduce dimensionality
plausibility_counts <- rbind(c(sum(plausibility_counts[c(1,2),1]),
                                sum(plausibility_counts[c(1,2),2])),
                             c(sum(plausibility_counts[c(4,5),1]),
                                sum(plausibility_counts[c(4,5),2])))

# Group RTs by plausibility ratings and attachment
summary_rts <- df.responses %>%
  group_by(items,attachment,rating_plausibility) %>%
  summarise(mean_rt = mean(rt_attachment)) 

p <- ggplot(summary_rts, aes(x = factor(rating_plausibility), y = mean_rt,items = items,fill=factor(rating_plausibility))) +
  geom_boxplot() +
  facet_wrap(~attachment) +
  ggtitle("RT distribution over items per plausibility bin")
p <- ggplotly(p,tooltip = c("items","mean_rt"))
p


summary_rts <- df.responses %>%
  group_by(items,rating_plausibility) %>%
  summarise(mean_hits = mean(hits)) 

p <- ggplot(summary_rts, aes(x = factor(rating_plausibility), y = mean_hits,items = items,fill=factor(rating_plausibility))) +
  geom_boxplot() +
  geom_jitter() +
  ggtitle("Hit distribution over items per plausibility bin")
p <- ggplotly(p,tooltip = c("items","y"))
p
# items that have low accuracy but partner has high accuracy might be repaired
# items where both have very low accuracy, attachment should be switched if makes sense
# Otherwise form completely new sentence

