library(data.table)
library(plyr)
library(dplyr)
library(DT)
library(sjPlot)
library(ggplot2)
library(plotly)
library(lme4)
library(RColorBrewer)

pastel_colors = brewer.pal(8, "Pastel1")
set_colors = brewer.pal(10, "PiYG")
#load log file
file_names = list.files(pattern = "*_log.txt",path = "logfiles/",full.names = TRUE)
tmp <- ldply(file_names[5:14],read.table, skip = 6, header = FALSE , sep = "\t", fill = TRUE,
                  col.names = c("Subject","Trial","Word","Condition","PairNum",
                                "VerbNum","Attachment","Time","Duration","Answers",
                                "Response","ResponseTime"))
df <- as.data.frame(tmp)


# keep only questions with condition information
df$Condition[which(df$ResponseTime!="")] = df$Condition[which(df$ResponseTime!="")-1]
df <- df[which(df$ResponseTime!=""),]

x <- df %>% group_by(Subject,Response, Condition) %>% summarise(n = n())
y <- df %>% group_by(Subject,Condition) %>% summarise(n = n())
acc = merge(x,y,by=c("Condition","Subject")) %>% transform(new = (n.x/n.y*100))
acc[(acc$Response==1),]

acc$Condition <- as.factor(acc$Condition)
acc$Condition <- mapvalues(acc$Condition, from = c("1","2","3"), to = c("Verb condition", "Role condition","Filler"))

# plotting behavioral performance for all subjects
p <- ggplot(data = acc[acc$Response==1,], aes(x = Condition, y = new, fill = Condition)) +
  geom_boxplot(fatten = NULL) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax=..y.., ymin=..y..),width=0.75, size = 1, linetype = "solid") +
  ggtitle("Per subject accuracy on comprehension questions")+
  labs(y="Accuracy", x="Manipulation Condition")+
  scale_x_discrete(labels = c(" "," "," "," ")) +
  theme(plot.title = element_text(hjust=0.5))
p + geom_dotplot(binaxis='y',stackdir='center', dotsize=1)
p

x <- df %>% group_by(Response, Condition) %>% summarise(n = n())
y <- df %>% group_by(Condition) %>% summarise(n = n())
acc = merge(x,y,by=c("Condition")) %>% transform(new = (n.x/n.y*100))
acc[acc$Response==1,]

# factorize and center
df$Subject <- as.factor(df$Subject)
df$Condition <- as.factor(df$Condition)
df$Response <- as.factor(df$Response)
df$ResponseTime <- as.numeric(df$ResponseTime)
contrasts(df$Response) <- c(-1,1)
df$ResponseTime <- scale(df$ResponseTime)

m <- glmer(Response ~ ResponseTime + Condition +
             ResponseTime*Condition +
             (1|Subject),
           data = df, family="binomial",control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(m)

p <- plot_model(m, type = "pred", terms = c("Condition"))
p

## load txt file with post-test results:
tmp <- read.table("post_test_mat.txt", header = FALSE , sep = "\t", fill = TRUE,
             col.names = c("Condition","Subject","Response"))
df <- as.data.frame(tmp)

x <- df %>% group_by(Subject,Response, Condition) %>% summarise(n = n())
y <- df %>% group_by(Subject,Condition) %>% summarise(n = n())
acc = merge(x,y,by=c("Condition","Subject")) %>% transform(new = (n.x/n.y*100))
acc[(acc$Response==1),]

acc$Condition <- as.factor(acc$Condition)
acc$Condition <- mapvalues(acc$Condition, from = c("1","2"), to = c("Verb condition", "Role condition"))

# plotting behavioral performance for all subjects
set_colors = brewer.pal(10, "PiYG")
set_colors <- paste(c(set_colors,"#4292C6","#08519C"), sep=" ")
p <- ggplot(data = acc[acc$Response==1,], aes(x = Condition, y = new, fill = Condition)) +
  geom_boxplot(fatten = NULL) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax=..y.., ymin=..y..),width=0.75, size = 1, linetype = "solid") +
  ggtitle("Average per subject accuracy on post-test")+
  labs(y="Accuracy", x="Manipulation Condition")+
  theme(plot.title = element_text(hjust=0.5))
#p + geom_dotplot(aes(fill = as.factor(Subject)), binaxis='y', stackdir='center', dotsize=1, position = position_jitter(0.2)) +
  scale_fill_manual(values = set_colors)

p  +
geom_point(binaxis='y', stackdir='center', dotsize=1) +
geom_line(aes(group = Subject),color="black", size=1, alpha=0.5)
