
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

responses <- read.csv(file="all_responses.csv",head=TRUE,sep=";")


names <- responses$surveynum

df.responses <- as.data.frame(t(responses[,-1]))
colnames(df.responses) <- names
df.responses$age <- as.numeric(df.responses$age)
#output mean age + sd of participants
#output mean time for survey taken



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
df.rts        <- df.responses[,subset_time]
y1            <- rowMeans(df.rts[,grep("^VA",names(df.rts))])
y2            <- rowMeans(df.rts[,grep("^NA",names(df.rts))])
y3            <- rowMeans(df.rts)

p <- plot_ly(type = 'box') %>%
  add_boxplot(y = y1, jitter = 0.3, pointpos = -1.8, boxpoints = 'all',
              marker = list(color = pastel_colors[1]),
              line = list(color = pastel_colors[1]),
              name = "Verb attachment",
              hoverinfo= row.names(df.rts)) %>%
  add_boxplot(y = y2, jitter = 0.3, pointpos = -1.8, boxpoints = 'all',
              marker = list(color = pastel_colors[2]),
              line = list(color = pastel_colors[2]),
              name = "Noun attachment") %>%
  add_boxplot(y = y3, name = "All trials", boxpoints = 'suspectedoutliers',
              marker = list(color = pastel_colors[3],
                            outliercolor = set_colors[1],
                            line = list(outliercolor = set_colors[1],
                                        outlierwidth = 2)),
              line = list(color = pastel_colors[3])) %>%
  layout(title = "Reaction times across items & subjects (data points are subject averages)",
         yaxis=list(title="RTs (seconds)"))




subset_attach <- grep("^NA[1-9]*$",names(df.responses),value = TRUE)
idx           <- match(subset_attach,names(df.responses))
df.Nattach    <- df.responses[,subset_attach]
subset_attach <- grep("^VA[1-9]*$",names(df.responses),value = TRUE)
idx           <- match(subset_attach,names(df.responses))
df.Vattach    <- df.responses[,subset_attach]

