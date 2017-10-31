library(plotly)
library(reshape)
library(RColorBrewer)
library(data.table)
library(dplyr)
#Create color scheme
pastel_colors = brewer.pal(8, "Pastel1")
set_colors = brewer.pal(8, "Set1")
#load stimulus file
stimuli         <- read.csv(file="stimuliA.csv",head=TRUE,sep=",",na.strings = c("","NAN"), stringsAsFactors = FALSE)
stimuli$fullID  <- as.character(interaction(stimuli[,c(13,3)],sep = ""))
#load response file
responses       <- read.csv(file="all_responses.csv",head=TRUE,sep=";", na.strings = c("","NA"),stringsAsFactors=FALSE)
#transpose and fix column names and classes
names                       <- responses$surveynum
df.responses                <- as.data.frame(t(responses[,-1]))
colnames(df.responses)      <- names
df.responses$age            <- as.numeric(df.responses$age)
subset_time                 <- grep("ime$",names(df.responses),value = TRUE)
df.responses[,subset_time]  <- lapply(df.responses[,subset_time], function(x) as.numeric(gsub(",",".",x)))
df.responses[sapply(df.responses, is.character)] <- lapply(df.responses[sapply(df.responses, is.character)], 
                                                           as.factor)
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
ind             <- grep("VA[0-9]*Time$",colnames(df.responses),value = TRUE)
df.respTime     <- df.responses[,c("subjects",ind),with=FALSE]
df.respTime     <- melt(df.respTime,id="subjects",value.name="rt_attachment",variable.name='items')
ind             <- grep("P[VN]A[0-9]*]$",colnames(df.responses),value = TRUE)
df.respPAttach  <- df.responses[,c("subjects",ind),with=FALSE]
df.respPAttach  <- melt(df.respPAttach,id="subjects",value.name="rating_plausibility",variable.name='items')
ind             <- grep("P[VN]A[0-9]*Time$",colnames(df.responses),value = TRUE)
df.respPTime    <- df.responses[,c("subjects",ind),with=FALSE]
df.respPTime    <- melt(df.respPTime,id="subjects",value.name="rt_plausibility",variable.name='items')
df.responses    <- cbind(df.respAttach,df.respTime[,3],df.respPAttach[,3],df.respPTime[,3])
df.responses$hits <- as.numeric((grepl("N",df.responses$items) & df.responses$response_attachment == "Nomen") | 
                      (grepl("V",df.responses$items) & df.responses$response_attachment == "Verb"))
#############ANALYSIS#################
summary(df.subjectinfo$age)
summary((df.subjectinfo$interviewtime)/60)

#find outlier (who did not pass the test?)
#number of wrong answers per participant x/4

df.subjectinfo %>%
  select(subjects,starts_with("correct")) %>%
    mutate(number_correct = rowSums(!is.na(.[,2:ncol(.)]))) %>%
      select(subjects,number_correct)

#average accuracy of unambiguous items
subset_unambiguous <- stimuli$fullID[(stimuli$Unambiguous.==1)]
sentences_unambiguous <- stimuli[(stimuli$Unambiguous.==1),4:12]

out_percentcorrect <- df.responses %>%
                      filter(items %in% subset_unambiguous) %>%
                      group_by(subjects) %>%
                        summarise(accuracy = mean(hits))

#Biases: Reaction Times total and per attachment

subset_timeV   <- grep("^VA[0-9]*Time$",names(df.responses),value = TRUE)
subset_timeN   <- grep("^NA[0-9]*Time$",names(df.responses),value = TRUE)
df.rtsV        <- df.responses[,subset_timeV]
df.rtsN        <- df.responses[,subset_timeN]
subset_attachV <- grep("^VA[0-9]*$",names(df.responses),value = TRUE)
subset_attachN <- grep("^NA[0-9]*$",names(df.responses),value = TRUE)
df.attachV     <- df.responses[,subset_attachV]
df.attachN     <- df.responses[,subset_attachN]
df.rtscorrectV <- df.
#first rearrange data.frame so that subjects, per attachment, per Rts, per correct (1/0)
#then compute summary as in average over trials
#then plot graph based on summary
y1            <- rowMeans(df.rtsV[,grep("^VA",names(df.rts))])
y2            <- rowMeans(df.rts[,grep("^NA",names(df.rts))])
y3            <- rowMeans(df.rts)
df.rts_means  <- data.frame(names(y1),y1,y2,y3)
colnames(df.rts_means) <- c("subject","Verb","Noun","All")
df.rts_means  <- melt(df.rts_means)
colnames(df.rts_means) <- c("Subject","Attachment","RTs")

p <- ggplot(df.rts_means, aes(x = Attachment, y = RTs, color = Attachment, subject = Subject)) +
    geom_boxplot() +
    geom_jitter(size = 2) +
    ggtitle("RTs averaged over items")

p <- ggplotly(p,tooltip = c("subject"))
p



