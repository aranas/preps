library(data.table)
library(dplyr)
library(DT)
library(lme4)
#load log file
file_names = list.files(pattern = "*_log.txt",path = "/project/3011210.01/logfiles/",full.names = TRUE)
tmp <- ldply(file_names,read.table, skip = 6, header = FALSE , sep = "\t", fill = TRUE,
                  col.names = c("Subject","Trial","Word","Condition","PairNum",
                                "VerbNum","Attachment","Time","Duration","Answers",
                                "Response","ResponseTime"))
df <- as.data.frame(tmp)


# keep only questions with condition information
df$Condition[which(df$ResponseTime!=0)] = df$Condition[which(df$ResponseTime!=0)-1]
df <- df[which(df$ResponseTime!=0),]

x <- df %>% group_by(Subject, Response) %>% summarise(n = n())

# factorize and center
df$Subject <- as.factor(df$Subject)
df$Condition <- as.factor(df$Subject)
df$Response <- as.factor(df$Response)
df$ResponseTime <- as.numeric(df$ResponseTime)
contrasts(df$Response) <- c(-1,1) 



m <- glmer(Response ~ 1 + ResponseTime + Condition +
             ResponseTime*Condition +
             (1|Subjects),
           data = df, family="binomial",control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(m)