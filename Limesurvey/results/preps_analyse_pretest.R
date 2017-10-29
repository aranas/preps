responses <- read.csv(file="all_responses.csv",head=TRUE,sep=";")

names <- responses$surveynum

df.responses <- as.data.frame(t(responses[,-1]))
colnames(df.responses) <- names
df.responses$age <- as.numeric(df.responses$age)
#output mean age + sd of participants
#output mean time for survey taken


subset_attach <- grep("^NA[1-9]*$",names(df.responses),value = TRUE)
idx           <- match(subset_attach,names(df.responses))
df.Nattach    <- df.responses[,subset_attach]
subset_attach <- grep("^VA[1-9]*$",names(df.responses),value = TRUE)
idx           <- match(subset_attach,names(df.responses))
df.Vattach    <- df.responses[,subset_attach]