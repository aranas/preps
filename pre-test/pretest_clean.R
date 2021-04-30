library(data.table)
library(DT)
library(dplyr)
library(plyr)

# Loading in data
#All following graphs and statistics are based on the results of the Limesurvey online pre-test. 
#40 participants filled out the online pre-test, 
#of which 20 saw the same sentences and each participant saw them in a different pseudo-randomized order.

stimuli         <- read.csv(file="stimuli_forpretest.csv",head=TRUE,sep=";",na.strings = c("","NAN"), stringsAsFactors = FALSE,encoding="UTF-8")
stimuli$fullID  <- as.character(interaction(stimuli[,c(13,3)],sep = ""))
stimuli$Condition <- as.factor(stimuli$Condition)
colnames(stimuli) <- c("Type","Pair","ID","Det.","N0","Verb","Det.2","N1","Prep.","Det.3","Adj.","N2","Attachment","Unambiguous","verb_number","fullID")
#load response file
file_names = list.files(pattern="^results",path = "allresponses",full.names = TRUE)
temp <- lapply(file_names,read.csv,sep=",",na.strings = c("","NAN"), stringsAsFactors = TRUE)
df.responses <- Reduce(function(x,y) merge(x,y,all=TRUE,sort=TRUE),temp)
#transpose and fix column names and classes
df.responses                <- setDT(df.responses,keep.rownames = TRUE)
colnames(df.responses)[1]   <- "subjects"
df.responses                <- Filter(function(x) !(all(x=="")), df.responses) #delete blank columns

#Data was subsequently split into subject info (demographics etc.) and responses to items

#Split all item specific data from subject specific data
ind_items   <- grep("[VN]A[0-9]",colnames(df.responses),value = TRUE)
df.subjectinfo <- df.responses
df.subjectinfo[,ind_items] <- NULL
df.responses  <- df.responses[,c("subjects",ind_items),with=FALSE]

#reshape items so that each row contains info for items per subject (items repeat over rows)
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


## Summary Statistics
#Age of participants:
summary(df.subjectinfo$age)

#Time taken to complete pre-test (in minutes):
summary((df.subjectinfo$interviewtime)/60)

## 1. Outlier
# Outlier subjects were identified according to three scores:  
# 1. Subjects who have more than one incorrect test item  
# 2. Subjects who have less than 60% correct for the unambiguous items  
# 3. Subjects who's average reaction time diverges extremely from the average

### 1a. Unambiguous sentences:
#average accuracy of unambiguous items
subset_unambiguous <- stimuli$fullID[(stimuli$Unambiguous==1)]
sentences_unambiguous <- stimuli[(stimuli$Unambiguous==1),4:12]
sentences_unambiguous

### 1b. Reaction Times
least_correct = 3
threshhold_acc = 0.60
avg_rts <- df.responses %>%
                        group_by(subjects) %>%
                        summarise(mean_rt = mean(rt_attachment,na.rm=TRUE)) 

outlier <- boxplot.stats(avg_rts$mean_rt)$out


test_correct   <- df.subjectinfo %>%
                   select(subjects,starts_with("correct")) %>%
                      mutate(number_correct = rowSums(!is.na(.[,2:ncol(.)]))) %>%
                        select(subjects,number_correct)

avg_unambiguous <-df.responses %>%
                    filter(items %in% subset_unambiguous) %>%
                      group_by(subjects) %>%
                        summarise(acc_unambiguous = mean(hits,na.rm=TRUE))

### 1c. Summary table per subject
summ <- join_all(list(avg_rts,test_correct,avg_unambiguous),by='subjects')
DT::datatable(summ) %>% formatRound(columns=colnames(summ),digits=2)


outlier_indx <- avg_rts$subjects[which(round(avg_rts$mean_rt,digits=2) %in% round(outlier,digits=2))]
outlier_indx <- c(outlier_indx,summ$subjects[which(summ$number_correct < least_correct | summ$acc_unambiguous < threshhold_acc)])
outlier_indx <- unique(outlier_indx)
#Remove Outlier
df.responses <- df.responses %>%
                filter(!(subjects %in% outlier_indx))

#Number of removed outliers: `r length(outlier_indx)`

## 2. Assess items

#Binomial distribution of probabilities if chance-performane is 50%
x <- seq(0,20,by = 1)
y <- dbinom(x,20,0.5)
plot(x,y,main="Probability distribution assuming per item probability of 50%",
      xlab="Number of correct answers per items (out of 20 answers)", ylab = "Probability")

### Reject stimuli with too few correct responses (less than 72% expected answers)
### 2a. Rejected items:
threshold = 0.72

acc_per_item <- df.responses %>%
                  group_by(items) %>%
                    summarise(acc = mean(hits,na.rm = TRUE), plaus = mean(rating_plausibility,na.rm = TRUE))
mean_accuracy <- acc_per_item %>%
                    summarise(mean_correct = mean(acc,na.rm = TRUE))

df.stim_acc <- merge(stimuli,acc_per_item,by.x="fullID",by.y="items")
df.stim_acc$acc <- round(df.stim_acc$acc,digits=2)
df.stim_acc$plaus <- round(df.stim_acc$plaus,digits=2)

items_reject <- df.stim_acc %>% arrange(Verb) %>% filter(acc < threshold)
                                                           
  #datatable(items_reject[,c(1,2,4,6:13,17)],
  #          filter = 'top',
  #          options = list("pageLength" = 10))
  datatable(items_reject[,c(1,2,4,6:13,17,18)])

###Sentences that have only around 50% expected responses
items_chance <- items_reject[items_reject$acc <= 0.6 & items_reject$acc >=0.4,]
pairs_chance <- df.stim_acc[df.stim_acc$ID %in% items_chance$ID,]  %>% arrange(ID)
#datatable(pairs_chance[,c(1,2,6:13,17)])
#It seems that for Type 2 items with low accuracy, that the semantics of N0 and N1 were both too related to N2, so that Verb attachment is preferred (high accuracy for verb attached sentences, low accuracy for noun attached sentences)

###Type 2 pairs that have only high accuracy (above 70%) across both Verb and Noun attached items
items_high <- df.stim_acc$ID[df.stim_acc$Type == "2" &  df.stim_acc$acc >= 0.7]
pairs_high <- items_high[duplicated(items_high)]
pairs_high <- df.stim_acc[df.stim_acc$ID %in% pairs_high,] %>% arrange(ID)
#datatable(pairs_high[,c(1,2,6:13,17)])

##Items with very low performance (less than 40%) for both items in a pair
items_low <- df.stim_acc$ID[df.stim_acc$acc <= 0.4]
pairs_low <- items_low[duplicated(items_low)]
pairs_low <- df.stim_acc[df.stim_acc$ID %in% pairs_low,] %>% arrange(ID)
#datatable(pairs_low[,c(1,2,6:13,17)])

##Pairs with very low performance (less than 40%) in one but high performance in paired item
pairs_onelow <- df.stim_acc[df.stim_acc$ID %in% items_low,] %>% arrange(ID)
#datatable(pairs_onelow[,c(1,2,6:13,17)])

### 2b. Remaining items:
#first reject items based on accuracy threshold with their paired items
pairs_reject <- df.stim_acc[df.stim_acc$ID %in% items_reject$ID,]
#for some 'repaired items there is a new sentence completing the pair
#for these items, we 1. only reject the 'bad' sentence and 2. relabel the others
#1.
pairs_reject <- pairs_reject[pairs_reject$fullID != "NA116",]
df.responses_clean <- df.responses %>%
                filter(!(items %in% pairs_reject$fullID))

#Then check if all verbs still occur twice (also for type 1 items) and reject remaining pair.
pairs_remain <- df.stim_acc[df.stim_acc$fullID %in% df.responses_clean$items,]
pairs_remain <-  unique(pairs_remain[duplicated(pairs_remain$verb_number),"verb_number"])
items_reject2 <- df.stim_acc[!(df.stim_acc$verb_number %in% pairs_remain),]
pairs_reject2 <- df.stim_acc[df.stim_acc$ID %in% items_reject2$ID,]
#1.
pairs_reject2 <- pairs_reject2[pairs_reject2$fullID != "NA116",]
pairs_reject2 <- pairs_reject2[pairs_reject2$fullID != "NA163",]
df.responses_clean <- df.responses_clean %>%
                filter(!(items %in% pairs_reject2$fullID))
#2.
levels(df.responses_clean$items) <- c(levels(df.responses_clean$items), "NA1110")
df.responses_clean$items[df.responses_clean$items == "NA116"] <- "NA1110"
df.responses_clean$items[df.responses_clean$items == "NA163"] <- "NA116"
df.responses_clean$items[df.responses_clean$items == "VA2105"] <- "VA116"
id <-"NA1110"
df.stim_acc$fullID[df.stim_acc$fullID == "NA116"] <- id
df.stim_acc$Pair[df.stim_acc$fullID == id] <- "110"
df.stim_acc$ID[df.stim_acc$fullID == id] <- "1110"
id <- "VA163"
df.stim_acc$fullID[df.stim_acc$fullID == "VA116"] <- id
df.stim_acc$Pair[df.stim_acc$fullID == id] <- "63"
df.stim_acc$ID[df.stim_acc$fullID == id] <- "163"
id <- "NA116"
df.stim_acc$fullID[df.stim_acc$fullID == "NA163"] <- id
df.stim_acc$Pair[df.stim_acc$fullID == id] <- "16"
df.stim_acc$ID[df.stim_acc$fullID == id] <- "116"
id <- "VA116"
df.stim_acc$fullID[df.stim_acc$fullID == "VA2105"] <- id
df.stim_acc$Pair[df.stim_acc$fullID == id] <- "16"
df.stim_acc$ID[df.stim_acc$fullID == id] <- "116"
df.stim_acc$Type[df.stim_acc$fullID == id] <- "1"

datatable(df.stim_acc[df.stim_acc$fullID %in% df.responses_clean$item,c(1,2,3,4,6:13,17,18)])


#manually reject superfluous items
items_reject3 = c("VA249","NA249","NA21","VA21","NA243","VA243","NA2104","VA2104","NA247","VA247","NA288","VA288","NA255","VA255","NA237","VA237","NA22","VA22","NA25","VA25","NA24","VA24","NA272","VA272","NA229","VA229","NA286","VA286","NA287","VA287","NA275","VA275","NA278","VA278","NA253","VA253","NA120","VA120","NA133","VA133","NA125","VA125","NA183","VA183","NA188","VA188","NA1101","VA1101","NA147","VA147","NA1117","VA1117","NA1118","VA1118","VA1109")
df.responses_clean <- df.responses_clean %>%
                filter(!(items %in% items_reject3))
datatable(df.stim_acc[df.stim_acc$fullID %in% df.responses_clean$item,c(1,2,4,6:13,17,18)])

#save the cleaned data
save(df.stim_acc,df.responses_clean,df.subjectinfo,file = 'pretestdata_clean.RData')
