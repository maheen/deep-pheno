#!/usr/bin/env Rscript

##Written by Maheen Shermohammed

subjID2 = commandArgs(trailingOnly = TRUE)

if (length(subjID2)!=1) {
  stop("Exactly one subject ID must be supplied (beiwe ID)", call.=FALSE)
}

options(scipen=999)
library("NHPoisson", lib.loc="~/Rpackages")
library("effects")
library("mice",lib.loc = "~/Rpackages")
library(lattice)
library(caret,lib.loc = "~/Rpackages")



############# PREPARING DATA ############

root <- "~/beiwe/sleep_prediction"

dir.create(sprintf("%s/pred_outputs/%s/graphs",root,subjID2),recursive = TRUE)
dir.create(sprintf("%s/pred_outputs/%s/data",root,subjID2),recursive = TRUE)


if (file.exists(sprintf("%s/pred_outputs/%s/data/%s_acc_scrOn_sleep_bymin2.rds",root,subjID2,subjID2))) {
  acc_scrOn_sleep_bymin2 <- readRDS(sprintf("%s/pred_outputs/%s/data/%s_acc_scrOn_sleep_bymin2.rds",root,subjID2,subjID2))
} else {
  
  #Clean up power state data 
  pow_data <- read.csv(sprintf("%s/output_files/%s_pow_data.csv",root,subjID2), stringsAsFactors = FALSE)
  pow_data$utc <- as.POSIXct(pow_data$timestamp/1000, origin="1970-01-01", tz="EST")
  pow_data$utc_trunc <- trunc(pow_data$utc,"mins")
  pow_data$tdif_min <- difftime(pow_data$utc_trunc,pow_data[1,"utc_trunc"],units="mins")
  screen_on_data <- subset(pow_data,event=="Screen turned on")
  
  #Generate a vector of every minute in the study with how many times
  #the screen was on during that minute
  scrOnFreq<-table(screen_on_data$tdif_min)
  myvect <- numeric(as.numeric(names(tail(scrOnFreq,1)))+1)
  myvect[as.numeric(names(scrOnFreq))+1] <- scrOnFreq
  maxmin <- length(myvect) #how many minutes I want to analyze
  
  ####Use the NHPoisson package to generate a line graph that is high when there
  ####is a greater density of Screen On's and low when theres a low density
  scrOnMinsEv <- POTevents.fun(T = myvect[1:maxmin], thres = 0.5)
  my_lint <- 90
  emplambda_scrOnMins <- emplambda.fun(posE = scrOnMinsEv$Px, t = 1:maxmin,lint = my_lint, tit = "Screen On Mins") #takes ~1 min to run
  #points(myvect[1:maxmin]/max(myvect[1:maxmin])*max(emplambda_scrOnMins$emplambda,na.rm=T))
  scrOnDensity <- emplambda_scrOnMins$emplambda[(my_lint/2+1):(maxmin - my_lint/2)]
  #plot(scrOnDensity,type='l')
  
  ### Combine screen on info with accelerometer info
  acc_data_bymin <- read.csv(sprintf("%s/output_files/%s_accelerometer_ByMin_var.csv",root,subjID2), stringsAsFactors = FALSE)
  acc_data_bymin <- na.omit(acc_data_bymin)
  all_mins <- seq(pow_data$utc_trunc[1], tail(pow_data$utc_trunc,1), "min")
  acc_scrOn_bymin <- merge(acc_data_bymin,data.frame(scrOnDensity,timefact=as.character(all_mins[(my_lint/2+1):(maxmin - my_lint/2)])),by="timefact")
  acc_scrOn_bymin$timefact <- as.POSIXct(acc_scrOn_bymin$timefact, tz="EST")
  
  ### Add in sleep info
  sleep_data <- read.csv(sprintf("~/beiwe/sleep_prediction/subj_folders/%s/%s_sleep_w_dates.csv",subjID2,subjID2), stringsAsFactors = FALSE)
  
  #make a minute-by-minute boolean vector where 0=awake, 1=asleep
  ordered_end_dates <- as.Date(sleep_data$End.Date)[order(as.Date(sleep_data$End.Date))]
  all_mins_acc <- seq(as.POSIXct(paste(ordered_end_dates[1]-1,"18:00:00"),tz="EST"),
                      as.POSIXct(tail(paste(ordered_end_dates,"18:00:00"),1),tz="EST"), "min")
  start_inds <- which(substr(as.character(all_mins_acc),1,16) %in% substr(paste(sleep_data$Start.Date,sleep_data$sleep_onset),1,16))
  end_inds <- which(as.character(all_mins_acc) %in% paste(sleep_data$End.Date,sleep_data$sleep_offset))
  sleep_inds_list <- lapply(1:length(start_inds), function(u) position=start_inds[u]:end_inds[u])
  all_nighters <- unlist(lapply(sleep_inds_list, function(u) length(u)==1))
  sleep_inds <- unlist(sleep_inds_list[!all_nighters])
  sleep_bool_df <- data.frame(sleep_bool=numeric(length(all_mins_acc)),timefact=all_mins_acc)
  sleep_bool_df$sleep_bool[sleep_inds] <- 1
  
  #Remove days where we don't have sleep data
  all_sleep_days <- seq(ordered_end_dates[1],tail(ordered_end_dates,1),"day")
  diff(c(dim(sleep_data)[1],length(all_sleep_days))) #missing one day of sleep data
  missing_day <- all_sleep_days[which(!(all_sleep_days %in% as.Date(sleep_data$End.Date)))]
  if (length(missing_day)!=0) {
    missing_minutes <- as.vector(t(outer(X=which(sleep_bool_df$timefact %in% as.POSIXct(paste(missing_day-1,"18:00:00"),tz="EST")), 
                                         Y=c(0:1439), FUN="+")))
    sleep_bool_df <- sleep_bool_df[-c(missing_minutes),]
  }
  
  acc_scrOn_bymin_thr <- subset(acc_scrOn_bymin,total>=5)
  print(sprintf("%s minutes with less than 5 data points, removed",dim(acc_scrOn_bymin)[1]-dim(acc_scrOn_bymin_thr)[1]))
  
  acc_scrOn_sleep_bymin <- merge(acc_scrOn_bymin_thr,sleep_bool_df,by="timefact")
  
  acc_scrOn_sleep_bymin2 <- merge(acc_scrOn_bymin_thr,sleep_bool_df,by="timefact",all.y = T)
  acc_scrOn_sleep_bymin2$log_activ <- log(acc_scrOn_sleep_bymin2$activity+1)
  acc_scrOn_sleep_bymin2$log_scrOnDens <- log(acc_scrOn_sleep_bymin2$scrOnDensity+1)
  acc_scrOn_sleep_bymin2$Time <- as.numeric(difftime(acc_scrOn_sleep_bymin2$timefact,acc_scrOn_sleep_bymin2[1,1],units="mins"))
  acc_scrOn_sleep_bymin2$Time_circ <- as.numeric(substr(acc_scrOn_sleep_bymin2$timefact,12,13))*60+as.numeric(substr(acc_scrOn_sleep_bymin2$timefact,15,16))
  acc_scrOn_sleep_bymin2$Time_sin <- sin(2*pi/(60*24)*acc_scrOn_sleep_bymin2$Time)
  acc_scrOn_sleep_bymin2$Time_cos <- cos(2*pi/(60*24)*acc_scrOn_sleep_bymin2$Time)
  
  acc_scrOn_sleep_bymin2$day <- as.Date(acc_scrOn_sleep_bymin2$timefact+(60*60*6), tz="EST")
  all_days <- seq(acc_scrOn_sleep_bymin2$day[1], tail(acc_scrOn_sleep_bymin2$day,1), "day")
  acc_scrOn_sleep_bymin2$day <- match(acc_scrOn_sleep_bymin2$day,all_days)
  
  saveRDS(acc_scrOn_sleep_bymin2,sprintf("%s/pred_outputs/%s/data/%s_acc_scrOn_sleep_bymin2.rds",root,subjID2,subjID2))
}

summary(acc_scrOn_sleep_bymin2)

#what proportion of full data is missing?
sum(is.na(acc_scrOn_sleep_bymin2$log_activ))/(dim(acc_scrOn_sleep_bymin2)[1])


#Plot all data w/o imputed values
if (file.exists(sprintf("%s/pred_outputs/%s/data/%s_pmm_acc_scrOn_sleep_bymin_full.rds",root,subjID2,subjID2))) {
  acc_scrOn_sleep_bymin_full <- readRDS(sprintf("%s/pred_outputs/%s/data/%s_pmm_acc_scrOn_sleep_bymin_full.rds",root,subjID2,subjID2))
} else {
  acc_scrOn_sleep_bymin_full <- acc_scrOn_sleep_bymin2
  acc_scrOn_sleep_bymin_full$include <- 'y'
  pdf(sprintf("%s/pred_outputs/%s/graphs/%s_pmm_activity.pdf",root,subjID2,subjID2),width=20,height=10)
  #pdf(sprintf("/ncf/somerville_lab/studies/RegiStre/beiwe/%s_pmm_activity.pdf",subjID2),width=20,height=10)
  op <- par(mfrow=c(11,11), mar = c(2,2,2,0.5))
  for (i in unique(acc_scrOn_sleep_bymin_full$day)[1]:tail(unique(acc_scrOn_sleep_bymin_full$day),1)) {
    thisday_inds <- which(acc_scrOn_sleep_bymin_full$day==i)
    all_vals <- acc_scrOn_sleep_bymin_full[thisday_inds,c("log_activ","Time","Time_circ")]
    if (nrow(all_vals)!=0 & nrow(na.omit(all_vals))==0) {
      acc_scrOn_sleep_bymin_full$include[thisday_inds] <- 'n'
    }
    today_inds <- c(1080:1439,0:1079)
    today_times <- which(today_inds %in% all_vals$Time_circ)
    plot(x=today_times,y=all_vals$log_activ,ylim=c(0,1.75),xlim=c(1,1440),xaxt = 'n')
    axis(1, at=c(1,361,721,1081,1441), labels=c("6pm","M","6am","N","6pm"))
    title(sprintf("Day %s",i), line = 0.1)
  }
  par(op)
  dev.off()
  acc_scrOn_sleep_bymin_full <- subset(acc_scrOn_sleep_bymin_full,include=='y')
  saveRDS(acc_scrOn_sleep_bymin_full,sprintf("%s/pred_outputs/%s/data/%s_pmm_acc_scrOn_sleep_bymin_full.rds",root,subjID2,subjID2))
  
  pdf(sprintf("%s/pred_outputs/%s/graphs/%s_pmm_scrOn.pdf",root,subjID2,subjID2),width=20,height=10)
  #pdf(sprintf("/ncf/somerville_lab/studies/RegiStre/beiwe/%s_pmm_scrOn.pdf",subjID2),width=20,height=10)
  op <- par(mfrow=c(11,11), mar = c(2,2,2,0.5))
  for (i in unique(acc_scrOn_sleep_bymin_full$day)[1]:tail(unique(acc_scrOn_sleep_bymin_full$day),1)) {
    thisday_inds <- which(acc_scrOn_sleep_bymin_full$day==i)
    all_vals <- acc_scrOn_sleep_bymin_full[thisday_inds,c("log_scrOnDens","Time","Time_circ")]
    if (nrow(all_vals)!=0 & nrow(na.omit(all_vals))==0) {
      acc_scrOn_sleep_bymin_full$include[thisday_inds] <- 'n'
    }
    today_inds <- c(1080:1439,0:1079)
    today_times <- which(today_inds %in% all_vals$Time_circ)
    plot(x=today_times,y=all_vals$log_scrOnDens,ylim=c(0,max(acc_scrOn_sleep_bymin_full$log_scrOnDens,na.rm=T)),xlim=c(1,1440),xaxt = 'n')
    axis(1, at=c(1,361,721,1081,1441), labels=c("6pm","M","6am","N","6pm"))
    title(sprintf("Day %s",i), line = 0.1)
  }
  par(op)
  dev.off()
}

if (file.exists(sprintf("%s/pred_outputs/%s/data/%s_tempData_pmm.rds",root,subjID2,subjID2))) {
  tempData_pmm <- readRDS(sprintf("%s/pred_outputs/%s/data/%s_tempData_pmm.rds",root,subjID2,subjID2))
} else {
  tempData_pmm <- mice(acc_scrOn_sleep_bymin_full[c("log_activ","log_scrOnDens","Time_sin","Time_cos")],method=c("pmm","pmm","",""),maxit=50,seed=500)
  saveRDS(tempData_pmm,sprintf("%s/pred_outputs/%s/data/%s_tempData_pmm.rds",root,subjID2,subjID2))
  pdf(sprintf("%s/pred_outputs/%s/graphs/%s_pmm_densityplot.pdf",root,subjID2,subjID2))
  print(densityplot(tempData_pmm))
  dev.off()
  pdf(sprintf("%s/pred_outputs/%s/graphs/%s_pmm_chains.pdf",root,subjID2,subjID2))
  print(plot(tempData_pmm))
  dev.off()
}

#Plot all data with imputed values
#par(mar = c(5, 4, 4, 2) + 0.1)
if (!file.exists(sprintf("%s/pred_outputs/%s/graphs/%s_pmm_activity_imp.pdf",root,subjID2,subjID2))) {
  imputed_active <- tempData_pmm$imp$log_activ
  minutes_imputed <- NULL
  pdf(sprintf("%s/pred_outputs/%s/graphs/%s_pmm_activity_imp.pdf",root,subjID2,subjID2),width=20,height=10)
  #pdf(sprintf("/ncf/somerville_lab/studies/RegiStre/beiwe/%s_pmm_activity_imp.pdf",subjID2),width=20,height=10)
  op <- par(mfrow=c(11,11), mar = c(2,2,2,0.5))
  for (i in unique(acc_scrOn_sleep_bymin_full$day)[1]:tail(unique(acc_scrOn_sleep_bymin_full$day),1)) {
    thisday_inds <- which(acc_scrOn_sleep_bymin_full$day==i)
    all_vals <- acc_scrOn_sleep_bymin_full[thisday_inds,c("log_activ","Time","Time_circ")]
    imputed_vals <- rowMeans(imputed_active[rownames(imputed_active) %in% rownames(all_vals),])
    all_vals[names(imputed_vals),"log_activ"] <- imputed_vals
    today_inds <- c(1080:1439,0:1079)
    today_times <- which(today_inds %in% all_vals$Time_circ)
    plot(x=today_times,y=all_vals$log_activ,ylim=c(0,1.75),xlim=c(1,1440),type="l",xaxt = 'n')
    axis(1, at=c(1,361,721,1081,1441), labels=c("6pm","M","6am","N","6pm"))
    title(sprintf("Day %s",i), line = 0.1)
    today_times_imp <- which(today_inds %in% acc_scrOn_sleep_bymin_full$Time_circ[rownames(acc_scrOn_sleep_bymin_full) %in% names(imputed_vals)])
    points(x=today_times_imp,y=imputed_vals,col="red",pch=16,cex=.7)
    minutes_imputed <- append(minutes_imputed,length(imputed_vals))
  }
  par(op)
  dev.off()
  
  pdf(sprintf("%s/pred_outputs/%s/graphs/%s_pmm_minutes_imputed.pdf",root,subjID2,subjID2))
  hist(minutes_imputed)
  dev.off()
  
  imputed_active <- tempData_pmm$imp$log_scrOnDens
  pdf(sprintf("%s/pred_outputs/%s/graphs/%s_pmm_scrOn_imp.pdf",root,subjID2,subjID2),width=20,height=10)
  #pdf(sprintf("/ncf/somerville_lab/studies/RegiStre/beiwe/%s_pmm_scrOn_imp.pdf",subjID2),width=20,height=10)
  op <- par(mfrow=c(11,11), mar = c(2,2,2,0.5))
  for (i in unique(acc_scrOn_sleep_bymin_full$day)[1]:tail(unique(acc_scrOn_sleep_bymin_full$day),1)) {
    thisday_inds <- which(acc_scrOn_sleep_bymin_full$day==i)
    all_vals <- acc_scrOn_sleep_bymin_full[thisday_inds,c("log_scrOnDens","Time","Time_circ")]
    imputed_vals <- rowMeans(imputed_active[rownames(imputed_active) %in% rownames(all_vals),])
    all_vals[names(imputed_vals),"log_scrOnDens"] <- imputed_vals
    today_inds <- c(1080:1439,0:1079)
    today_times <- which(today_inds %in% all_vals$Time_circ)
    plot(x=today_times,y=all_vals$log_scrOnDens,ylim=c(0,max(acc_scrOn_sleep_bymin_full$log_scrOnDens,na.rm=T)),xlim=c(1,1440),type="l",xaxt = 'n')
    axis(1, at=c(1,361,721,1081,1441), labels=c("6pm","M","6am","N","6pm"))
    title(sprintf("Day %s",i), line = 0.1)
    today_times_imp <- which(today_inds %in% acc_scrOn_sleep_bymin_full$Time_circ[rownames(acc_scrOn_sleep_bymin_full) %in% names(imputed_vals)])
    points(x=today_times_imp,y=imputed_vals,col="red",pch=16,cex=.7)
  }
  par(op)
  dev.off()
}


#Function to do a cross-validation in parallel
cv_20fold <- function(myformula, dat, k, pred_cols, myseed) {
  
  fpr <- NULL # False positive rate
  fnr <- NULL # False negative rate
  acc <- NULL # Accuracy
  
  set.seed(myseed)
  
  for(i in 1:k)
  {
    # Train-test splitting
    # 95% of samples -> fitting
    # 5% of samples -> testing
    smp_size <- floor(0.95 * nrow(dat))
    index <- sample(seq_len(nrow(dat)),size=smp_size)
    train <- dat[index, ]
    test <- dat[-index, ]
    
    # Fitting
    model <- glm(myformula,family=binomial,data=train)
    
    # Predict results
    results_prob <- predict(model,subset(test,select=pred_cols),type='response')
    
    # If prob > 0.5 then 1, else 0
    results <- ifelse(results_prob > 0.5,1,0)
    
    # Actual answers
    answers <- test$sleep_bool
    
    # Accuracy calculation
    misClasificError <- mean(answers != results)
    
    # Collecting results
    acc[i] <- 1-misClasificError
    
    # Confusion matrix
    cm <- confusionMatrix(data=results, reference=answers)
    fpr[i] <- cm$table[2]/(nrow(dat)-smp_size)
    fnr[i] <- cm$table[3]/(nrow(dat)-smp_size)
  }
  return(list(acc=acc,fpr=fpr,fnr=fnr))
  
}

#Function to graphROC curve
library("ROCR", lib.loc="~/Rpackages")
library("ggplot2")
graph_roc <- function(mod,df,outcome) {
  prob <- predict(mod, newdata=df, type="response")
  pred <- prediction(prob, df[,outcome])
  perf <- performance(pred, measure = "tpr", x.measure = "fpr")
  tpr <- unlist(slot(perf, "y.values"))
  fpr <- unlist(slot(perf, "x.values"))
  roc <- data.frame(tpr, fpr)
  auc <- performance(pred, measure = "auc")
  auc <- auc@y.values[[1]]
  ggplot(roc) + geom_line(aes(x = fpr, y = tpr)) +
    geom_abline(intercept = 0, slope = 1, colour = "gray") + 
    ylab("Sensitivity") + 
    xlab("1 - Specificity") + ggtitle(paste("auc =",auc))
}


# Fit model with existing data, run cross-validation, determine accuracy for each imputated data set
my_seed <- 123
if (file.exists(sprintf("%s/pred_outputs/%s/data/%s_pmm_acc_scrOn_sleep_bymin_impAvg.rds",root,subjID2,subjID2))) {
  acc_scrOn_sleep_bymin_impAvg <- readRDS(sprintf("%s/pred_outputs/%s/data/%s_pmm_acc_scrOn_sleep_bymin_impAvg.rds",root,subjID2,subjID2))
} else {
  #Fit model with existing data
  acc_scrOn_sleep_bymin_real <- na.omit(acc_scrOn_sleep_bymin_full)
  fitLog_acc_scrOn_sin <- glm(sleep_bool ~ log_activ + log_scrOnDens + Time_sin + Time_cos,family=binomial,data=acc_scrOn_sleep_bymin_real)
  saveRDS(fitLog_acc_scrOn_sin,sprintf("%s/pred_outputs/%s/data/%s_pmm_fitLog_acc_scrOn_sin.rds",root,subjID2,subjID2))
  acc_scrOn_sleep_bymin_real$pred <- round(predict(fitLog_acc_scrOn_sin,acc_scrOn_sleep_bymin_real,type='response'))
  acc_real <- mean(acc_scrOn_sleep_bymin_real$sleep_bool == acc_scrOn_sleep_bymin_real$pred)
  print(paste("accuracy of model predicting only real data=",acc_real))
  
  #Plot ROC Curve
  pdf(sprintf("%s/pred_outputs/%s/graphs/%s_pmm_ROC_real.pdf",root,subjID2,subjID2,my_seed))
  print(graph_roc(fitLog_acc_scrOn_sin,acc_scrOn_sleep_bymin_real,"sleep_bool"))
  dev.off()
  
  #Run cross-validation on existing data
  cv_imp_output <- cv_20fold(as.formula("sleep_bool ~ log_activ + log_scrOnDens + Time_sin + Time_cos"),
                             acc_scrOn_sleep_bymin_real,200,c(6,7,10,11),my_seed)
  saveRDS(cv_imp_output,sprintf("%s/pred_outputs/%s/data/%s_pmm_cv_real_output_seed%s.rds",root,subjID2,subjID2,my_seed))
  acc_mean <- round(mean(cv_imp_output$acc),4)
  acc_sd <- round(sd(cv_imp_output$acc),4)
  fpr_mean <- round(mean(cv_imp_output$fpr),2)
  fnr_mean <- round(mean(cv_imp_output$fnr),2)
  pdf(sprintf("%s/pred_outputs/%s/graphs/%s_pmm_cv_real_acc_seed%s.pdf",root,subjID2,subjID2,my_seed))
  hist(cv_imp_output$acc,xlab='Accuracy',ylab='Freq',density=30,main=sprintf("acc=%s (fpr=%s,fnr=%s)\nacc sd =%s ",acc_mean,fpr_mean,fnr_mean,acc_sd))
  dev.off()
  
  #Determine accuracy for each imputed data set, as well as a data set that is the average of them
  all_log_activ <- NULL
  all_scrOnDens <- NULL
  acc_allImp <- NULL
  imputations  = 1:tempData_pmm$m
  for (m in imputations) {
    completedData <- complete(tempData_pmm,m)
    completedData$sleep_bool <- acc_scrOn_sleep_bymin_full$sleep_bool
    all_log_activ <- cbind(all_log_activ,completedData$log_activ)
    all_scrOnDens <- cbind(all_scrOnDens,completedData$log_scrOnDens)
    thisImp_pred <- round(predict(fitLog_acc_scrOn_sin,completedData,type='response'))
    acc_allImp[m] <- mean(completedData$sleep_bool == thisImp_pred)
    print(sprintf("accuracy of model predicting imputed dataset %s = %s",m,acc_allImp[m]))
  }
  
  acc_scrOn_sleep_bymin_impAvg <- acc_scrOn_sleep_bymin_full
  acc_scrOn_sleep_bymin_impAvg$log_activ <- rowMeans(all_log_activ)
  acc_scrOn_sleep_bymin_impAvg$log_scrOnDens <- rowMeans(all_scrOnDens)
  acc_scrOn_sleep_bymin_impAvg$pred <- round(predict(fitLog_acc_scrOn_sin,acc_scrOn_sleep_bymin_impAvg,type='response'))
  acc_impAvg <- mean(acc_scrOn_sleep_bymin_impAvg$sleep_bool == acc_scrOn_sleep_bymin_impAvg$pred)
  print(paste("accuracy of model predicting average imputed data=",acc_impAvg))
  
  thistitle <- sprintf("accuracy for avg imputed data = %s",round(acc_impAvg,4))
  pdf(sprintf("%s/pred_outputs/%s/graphs/%s_pmm_imp_acc_seed%s_full.pdf",root,subjID2,subjID2,my_seed))
  print(qplot(imputations,acc_allImp,ylim=c(0.5,1),main=thistitle,ylab="accuracy")+geom_errorbar(aes(x=imputations, ymin=acc_allImp-acc_sd, ymax=acc_allImp+acc_sd), width=0.25))
  dev.off()
  
  pdf(sprintf("%s/pred_outputs/%s/graphs/%s_pmm_ROC_impAvg.pdf",root,subjID2,subjID2,my_seed))
  print(graph_roc(fitLog_acc_scrOn_sin,acc_scrOn_sleep_bymin_impAvg,"sleep_bool"))
  dev.off()
  
  saveRDS(acc_scrOn_sleep_bymin_impAvg,sprintf("%s/pred_outputs/%s/data/%s_pmm_acc_scrOn_sleep_bymin_impAvg.rds",root,subjID2,subjID2))
}

if (!file.exists(sprintf("%s/pred_outputs/%s/graphs/%s_pmm_activity_pred.pdf",root,subjID2,subjID2))) {
  imputed_active <- tempData_pmm$imp$log_activ
  pdf(sprintf("%s/pred_outputs/%s/graphs/%s_pmm_activity_pred.pdf",root,subjID2,subjID2),width=20,height=10)
  #pdf(sprintf("/ncf/somerville_lab/studies/RegiStre/beiwe/%s_pmm_activity_pred.pdf",subjID2),width=20,height=10)
  op <- par(mfrow=c(11,11), mar = c(2,2,2,0.5))
  for (i in unique(acc_scrOn_sleep_bymin_impAvg$day)[1]:tail(unique(acc_scrOn_sleep_bymin_impAvg$day),1)) {
    thisday_inds <- which(acc_scrOn_sleep_bymin_impAvg$day==i)
    all_vals <- acc_scrOn_sleep_bymin_impAvg[thisday_inds,c("log_activ","Time","Time_circ","sleep_bool","pred")]
    imputed_vals <- rowMeans(imputed_active[rownames(imputed_active) %in% rownames(all_vals),])
    today_inds <- c(1080:1439,0:1079)
    today_times <- which(today_inds %in% all_vals$Time_circ)
    plot(x=today_times,y=all_vals$log_activ,ylim=c(0,1.75),xlim=c(1,1440),type="l",xaxt = 'n')
    axis(1, at=c(1,361,721,1081,1441), labels=c("6pm","M","6am","N","6pm"))
    title(sprintf("Day %s",i), line = 0.1)
    today_times_imp <- which(today_inds %in% acc_scrOn_sleep_bymin_impAvg$Time_circ[rownames(acc_scrOn_sleep_bymin_impAvg) %in% names(imputed_vals)])
    points(x=today_times_imp,y=imputed_vals,col="red",pch=16,cex=.7)
    points(x=today_times,y=(all_vals$sleep_bool*0.5),col="blue",type="l",lwd=3)
    points(x=today_times,y=(all_vals$pred*0.5)+1,col="green",type="l",lwd=3)
  }
  par(op)
  dev.off()
  
  imputed_active <- tempData_pmm$imp$log_scrOnDens
  pdf(sprintf("%s/pred_outputs/%s/graphs/%s_pmm_scrOn_pred.pdf",root,subjID2,subjID2),width=20,height=10)
  #pdf(sprintf("/ncf/somerville_lab/studies/RegiStre/beiwe/%s_pmm_scrOn_pred.pdf",subjID2),width=20,height=10)
  op <- par(mfrow=c(11,11), mar = c(2,2,2,0.5))
  for (i in unique(acc_scrOn_sleep_bymin_impAvg$day)[1]:tail(unique(acc_scrOn_sleep_bymin_impAvg$day),1)) {
    thisday_inds <- which(acc_scrOn_sleep_bymin_impAvg$day==i)
    all_vals <- acc_scrOn_sleep_bymin_impAvg[thisday_inds,c("log_scrOnDens","Time","Time_circ","sleep_bool","pred")]
    imputed_vals <- rowMeans(imputed_active[rownames(imputed_active) %in% rownames(all_vals),])
    today_inds <- c(1080:1439,0:1079)
    today_times <- which(today_inds %in% all_vals$Time_circ)
    scrOnMax <- max(acc_scrOn_sleep_bymin_full$log_scrOnDens,na.rm=T)
    plot(x=today_times,y=all_vals$log_scrOnDens,ylim=c(0,scrOnMax),xlim=c(1,1440),type="l",xaxt = 'n')
    axis(1, at=c(1,361,721,1081,1441), labels=c("6pm","M","6am","N","6pm"))
    title(sprintf("Day %s",i), line = 0.1)
    today_times_imp <- which(today_inds %in% acc_scrOn_sleep_bymin_impAvg$Time_circ[rownames(acc_scrOn_sleep_bymin_impAvg) %in% names(imputed_vals)])
    points(x=today_times_imp,y=imputed_vals,col="red",pch=16,cex=.7)
    points(x=today_times,y=(all_vals$sleep_bool*0.3*scrOnMax),col="blue",type="l",lwd=3)
    points(x=today_times,y=(all_vals$pred*0.3*scrOnMax)+.6*scrOnMax,col="green",type="l",lwd=3)
  }
  par(op)
  dev.off()
}

if (!file.exists(sprintf("%s/pred_outputs/%s/graphs/%s_pmm_duration_pred.pdf",root,subjID2,subjID2))) { 
  
  sleep_durations <- aggregate(cbind(sleep_bool,pred) ~ day,data=acc_scrOn_sleep_bymin_impAvg,FUN=sum)
  sleep_durations <- subset(sleep_durations,sleep_bool>0)
  
  cor_info <- cor.test(~ sleep_bool + pred,data=sleep_durations)
  
  pdf(sprintf("%s/pred_outputs/%s/graphs/%s_pmm_duration_pred.pdf",root,subjID2,subjID2))
  plot(sleep_bool ~ day,data=sleep_durations,type="l",xlab="Study Day",ylab="Sleep Duration",main=sprintf("r=%s,p=%s",round(cor_info$estimate,3),round(cor_info$p.value,3)),ylim=c(min(sleep_durations[2:3]),max(sleep_durations[2:3])))
  points(pred ~ day,data=sleep_durations,type="l",col="red")
  dev.off()
}

print(paste("done",subjID2))
