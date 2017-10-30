library(plotly)
library(reshape)
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
subset_time   <- grep("^[VN]A[0-9]*Time$",names(df.responses),value = TRUE)
rts_persubj   <- melt(t(df.responses[,subset_time]))
colnames(rts_persubj) <- c("itemcode","subject","rt")
rts_persubj$attachment <- rep("Noun",nrow(rts_persubj))
rts_persubj[grep("^VA",rts_persubj$itemcode),4] <- rep("Verb",nrow(rts_persubj)/2)

p <- plot_ly(type = 'box') %>%
  add_boxplot(y = y1, jitter = 0.3, pointpos = -1.8, boxpoints = 'all',
              marker = list(color = 'rgb(7,40,89)'),
              line = list(color = 'rgb(7,40,89)'),
              name = "Verb attachment") %>%
  add_boxplot(y = y1, jitter = 0.3, pointpos = -1.8, boxpoints = 'all',
              marker = list(color = 'rgb(7,40,89)'),
              line = list(color = 'rgb(7,40,89)'),
              name = "Noun attachment") %>%
  add_boxplot(y = y1, jitter = 0.3, pointpos = -1.8, boxpoints = 'all',
              marker = list(color = 'rgb(7,40,89)'),
              line = list(color = 'rgb(7,40,89)'),
              name = "All trials") %>%
  layout(title = "Reaction times across items & subjects (data points are subject averages")


subset_attach <- grep("^NA[0-9]*$",names(df.responses),value = TRUE)
idx           <- match(subset_attach,names(df.responses))
df.Nattach    <- df.responses[,subset_attach]
subset_attach <- grep("^VA[0-9]*$",names(df.responses),value = TRUE)
idx           <- match(subset_attach,names(df.responses))
df.Vattach    <- df.responses[,subset_attach]