#Written by Maheen Shermohammed

options(scipen=999)
library(xts)
require(stringr)
require(tools)
require(ggplot2)

upload_date <- '2016-10-31'

root <- sprintf("~/Documents/DeepPheno/430_smartphone/data_%s/",upload_date)
output_dir <- "~/Documents/DeepPheno/430_smartphone/social_analysis"
subj <- "mysubjID"
dir.create(sprintf("%s/%s/%s",output_dir,upload_date,subj),recursive = TRUE)

### ---- PREPPING SURVEY DATA ----

filenames <- list.files(path = paste(root,subj,"/surveyAnswers/5609db0797013e48c90597e6/",sep=""), pattern = "*\\.csv$")
timestamps <- file_path_sans_ext(filenames)
timestamps <- as.POSIXct(timestamps, format="%Y-%m-%d %H_%M_%S",tz="EST")-60*60*5

#import and format survey data
survey_data <- data.frame()
for (run in 1:length(filenames)) {
  data <- read.csv(paste(root,subj,"/surveyAnswers/5609db0797013e48c90597e6/",filenames[run],sep=""), stringsAsFactors = FALSE)
  
  ans_choices <- str_split_fixed(data[,4], "\\[|\\]", 3)[,2]
  
  #some answer choices have ";" in them but this is also a separator. Replace these with "," manually.
  ans_choices <- gsub("coffee; tea; or", "coffee, tea, or", ans_choices)
  ans_choices <- gsub("Typical; some moments", "Typical, some moments", ans_choices)
  ans_choices <- gsub("inimal; ", "inimal, ", ans_choices)
  ans_choices <- gsub("N/A; I", "N/A, I", ans_choices)
  data$answer <- gsub(";",",",data$answer)
  
  #convert answer choices to a scale
  answer_num <- NULL
  for (i in 1:dim(data)[1]) {
    ind<-which(strsplit(ans_choices[i],"; |;")[[1]] %in% data$answer[i])
    if ( !length(ind) ) {ind<- NA}
    answer_num[i] <- ind
  }
  survey_data <- rbind(survey_data,answer_num)
}
row.names(survey_data) <- timestamps
names(survey_data) <- c("energy","sleep","caffeine","alcohol","menstruation","cramping",
                        "anxiety_24hrs","anxiety_now","appetite","eating","exercise","walking",
                        "stomach","talkative","outgoing","reserved","social_person",
                        "social_digital","use_phone","phone_at_night","upset","hostile","alert",
                        "ashamed","inspired","nervous","determined","attentive","afraid","active",
                        "physically_active","stress","stress_managing","happy","lonely",
                        "extroverted","want_alone_time")
head(survey_data)

#remove data that was submitted between 3am-5pm
noDate <- strftime(timestamps, format="%H:%M:%S") #strip date from time
times <- as.POSIXct(noDate, format="%H:%M:%S") #reintroduce the same date across all times
limit_3am <- as.POSIXct("03:00:00", format="%H:%M:%S")
limit_5pm <- as.POSIXct("17:00:00", format="%H:%M:%S")
midnight <- as.POSIXct("00:00:01", format="%H:%M:%S")
after_midnight_ind <- (times<limit_3am & times>midnight)
survey_data <- survey_data[(times<limit_3am | times>limit_5pm),]
timestamps_cal <- timestamps[(times<limit_3am | times>limit_5pm)]
after_midnight_ind <- after_midnight_ind[(times<limit_3am | times>limit_5pm)]


### ---- PREPPING TEXT DATA ----

filenames <- list.files(path = paste(root,subj,"/textsLog/",sep=""), pattern = "*\\.csv$")
dataList <- lapply(filenames,function(x){read.csv(paste(root,subj,"/textsLog/",x,sep=""), stringsAsFactors = FALSE)}) 
text_data <- do.call(rbind, dataList) 
text_data <- text_data[-2] #remove UTC time column

#find duplicates
tbuffer <- 30000 #look within a time window of 30s
dupes <- NULL
for (i in 1:(dim(text_data)[1])) {
  if (i==1) {
    dupes[i] <- FALSE
  } else {
    if (text_data$timestamp[i]-text_data$timestamp[i-1] < tbuffer) {
      dupes[i] <- duplicated(text_data[(i-1):i,-c(1,5)])[2]
    } else {
      dupes[i] <- FALSE
    }
  }
}

text_data <- text_data[!dupes,] #only include non-duplicates in data
sent_texts <- text_data[text_data$sent.vs.received == "sent SMS",c(1,4)]
sent_times <- as.POSIXct(sent_texts$timestamp/1000, origin="1970-01-01")

sent_txts_xts <- xts(sent_texts$message.length,sent_times-(5*60*60)) #subtract 5 hours, so 5am is "midnight"
colnames(sent_txts_xts) <- "msglength"
sent_txts_xts$number <- rep(1,dim(sent_txts_xts)[1])
aggreg_sent_txts <- apply.daily(sent_txts_xts, colSums)

real_aggreg_times <- index(aggreg_sent_txts)+(5*60*60) #real time of last text for each day there was a text (end of day = 5am)
all_days_txt <- seq(as.Date(index(aggreg_sent_txts)[1]),as.Date(tail(index(aggreg_sent_txts),1)),"day")
day_of_study <- match(as.Date(index(aggreg_sent_txts), tz='EST'),all_days_txt)

sent_texts_number <- rep(0,length(all_days_txt))
sent_texts_number[day_of_study] <- aggreg_sent_txts$number
sent_texts_msglength <- rep(0,length(all_days_txt))
sent_texts_msglength[day_of_study] <- aggreg_sent_txts$msglength

png(sprintf("%s/%s/%s/NumberSent_sfn.png",output_dir,upload_date,subj))
ggplot(aes(x = 1:length(all_days_txt), y = sent_texts_number), data = data.frame(all_days_txt,sent_texts_number)) +  
  geom_point(size=3) + geom_line() + xlab("") + ylab("") + theme_classic() + ylim(0,50) + theme(axis.text.y = element_blank(),axis.text.x = element_text(size=15))
dev.off()

png(sprintf("%s/%s/%s/CharactersSent.png",output_dir,upload_date,subj))
plot(x=(1:length(all_days_txt)),y=sent_texts_msglength, pch = 16, type = "o", main = "Number of Characters Texted", xlab = "day", ylab="number of characters")
dev.off()


### ---- PREPPING CALL DATA ----

filenames <- list.files(path = paste(root,subj,"/callLog/",sep=""), pattern = "*\\.csv$")
dataList <- lapply(filenames,function(x){read.csv(paste(root,subj,"/callLog/","/",x,sep=""), stringsAsFactors = FALSE)}) 
call_data <- do.call(rbind, dataList) 

sum(duplicated(call_data)) #this should be 0, so no duplicates
answered_calls <- call_data[call_data$duration.in.seconds > 0,c(1,5,4)]
call_time <- as.POSIXct(answered_calls[,1]/1000, origin="1970-01-01")

calls_xts <- xts(answered_calls$duration.in.seconds,call_time-(5*60*60)) #subtract 5 hours, so 5am is "midnight"
colnames(calls_xts) <- "duration"
calls_xts$number <- rep(1,dim(calls_xts)[1])
aggreg_calls <- apply.daily(calls_xts, colSums)

real_aggreg_times <- index(aggreg_calls)+(5*60*60) #real time of last text for each day there was a text (end of day = 5am)
all_days_calls <- seq(as.Date(index(aggreg_calls)[1]),as.Date(tail(index(aggreg_calls),1)),"day")
day_of_study <- match(as.Date(index(aggreg_calls), tz='EST'),all_days_calls)

answered_calls_number <- rep(0,length(all_days_calls))
answered_calls_number[day_of_study] <- aggreg_calls$number
answered_calls_duration <- rep(0,length(all_days_calls))
answered_calls_duration[day_of_study] <- aggreg_calls$duration

png(sprintf("%s/%s/%s/NumberAnswered.png",output_dir,upload_date,subj))
plot(x=(1:length(all_days_calls)),y=answered_calls_number, pch = 16, type = "o", main = "Number of Answered Calls", xlab = "day", ylab="number of calls")
dev.off()

png(sprintf("%s/%s/%s/TimeSpent.png",output_dir,upload_date,subj))
plot(x=(1:length(all_days_calls)),y=answered_calls_duration, pch = 16, type = "o", main = "Time Spent on Phone", xlab = "day", ylab="time (s)")
dev.off()

### ---- COMBINED ANALYSES ----

calls_df <- data.frame(duration=answered_calls_duration,call_num=answered_calls_number,row.names=as.character(all_days_calls))
texts_df <- data.frame(msglength=sent_texts_msglength,text_num=sent_texts_number,row.names=as.character(all_days_txt))

timestamps_cal[after_midnight_ind] <- timestamps_cal[after_midnight_ind] - 3*3600 #subtract 3 hrs from submissions after midnight, so date is correct in next step
social_survey_df <- as.data.frame(survey_data[,c(14:19,35,37)],row.names=as.character(as.Date(timestamps_cal, tz='EST')))

passive_social <- merge(calls_df, texts_df, by="row.names", all=TRUE)
passive_social[is.na(passive_social)] <- 0
all_social <- merge(passive_social, social_survey_df, by.x="Row.names",by.y="row.names", all=TRUE)

#plot correlation matrix using first principal component to order it
library(corrplot)
social_cormat <- cor(all_social, use = "pairwise.complete.obs")
png(sprintf("%s/%s/%s/cormat_FPC.png",output_dir,upload_date,subj))
corrplot(social_cormat, method = "color", order = "FPC")
dev.off()

#make scree plot
require(psych)
scree(social_cormat)
png(sprintf("%s/%s/%s/scree_plot.png",output_dir,upload_date,subj))
scree(social_cormat,pc=FALSE)
dev.off()

#run factor analysis
FAfit <- factanal(all_social[complete.cases(all_social),], factors = 2)
print(FAfit, digits=2, cutoff=.3, sort=TRUE)

# plot factor 1 by factor 2 
load <- FAfit$loadings[,1:2] 
png(sprintf("%s/%s/%s/factor_plot.png",output_dir,upload_date,subj))
plot(load,type="n") # set up plot
text(load,labels=names(all_social[complete.cases(all_social),]),cex=.7) # add variable names
dev.off()

load2 <- load
row.names(load2)[1:4] <- c("call_length","call_number","text_length","text_number")
pdf(sprintf("%s/%s/%s/factor_plot_sfn.pdf",output_dir,upload_date,subj))
plot(load2,type="n",xlim=c(min(load[,1])-0.15,max(load[,1])+0.25),xlab="",ylab="") # set up plot 
text(load2,labels=row.names(load2),cex=1, col=c(rep("#A51C30",4),rep("black",8)),font=2) # add variable names
dev.off()

#make correlation matrix ordered by loadings from first factor of factor analysis
social_cormat2 <- cor(all_social[names(sort(load[,1]))], use = "pairwise.complete.obs")
png(sprintf("%s/%s/%s/cormat_FA1.png",output_dir,upload_date,subj))
corrplot(social_cormat2, method = "color")#, tl.col = "black")
dev.off()
