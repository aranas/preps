library(plotly)
library(reshape)
library("RColorBrewer")
pastel_colors = brewer.pal(8, "Pastel1")
set_colors = brewer.pal(8, "Set1")
#load stimulus file
stimuli         <- read.csv(file="stimuliA.csv",head=TRUE,sep=",",na.strings = "NAN", stringsAsFactors = FALSE)
stimuli$fullID  <- as.character(interaction(stimuli[,c(13,3)],sep = ""))
#load response file
responses       <- read.csv(file="all_responses.csv",head=TRUE,sep=";",stringsAsFactors=FALSE)
names           <- responses$surveynum


df.responses                <- as.data.frame(t(responses[,-1]))
colnames(df.responses)      <- names
df.responses$age            <- as.numeric(df.responses$age)
subset_time                 <- grep("ime$",names(df.responses),value = TRUE)
df.responses[,subset_time]  <- lapply(df.responses[,subset_time], function(x) as.numeric(gsub(",",".",x)))
df.responses[sapply(df.responses, is.character)] <- lapply(df.responses[sapply(df.responses, is.character)], 
                                                           as.factor)

summary(df.responses$age)
summary((df.responses$interviewtime)/60)

#find outlier (who did not pass the test?)
#number of wrong answers per participant x/4
subset_test   <- grep("correct[0-9]*Time$",names(df.responses),value=TRUE)
apply(df.responses[,subset_test], 1, function(x) sum(is.na(x)))
#average accuracy of unambiguous items
subset_unambiguous <- stimuli$fullID[(stimuli$Unambiguous.==1)]
sentences_unambiguous <- stimuli[(stimuli$Unambiguous.==1),4:12]
count_correct = apply(df.responses[,subset_unambiguous], 1, function(x) sum(x == "Nomen"))
count_correct/length(subset_unambiguous)
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



