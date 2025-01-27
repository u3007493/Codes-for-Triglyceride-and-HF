# Codes-for-Triglyceride-and-HF
All Codes used in Triglyceride levels and its association with all-cause mortality and cardiovascular outcomes among patients with heart failure


##R codes#####
library(rms)
library(dplyr)
load("~/Desktop/Triglyceride.rdata")
indf <- Triglyceride  
#Set data environment
dd <- datadist(indf)
options(datadist="dd")
#library(survival)
S <- Surv(indf$survival.death,indf$All.death==1)
S <- Surv(indf$survival.death,indf$CV.Death==1)
S <- Surv(indf$HF.survival,indf$HF==1)
S <- Surv(indf$ASCVD.survival,indf$ASCVD==1)
#library(rms)
fit <- cph(S ~ rcs(Triglyceride, 4)
           +use.age+ use.Sex
           +use.Smoke+use.Alcohol+
             use.HT+use.DM+use.obesity+Stroke+
             use.CAD+ use.PVD+use.Dyslipidemia+ 
             use.AF+ use.Anemia+ 
             use.cirrhosis+use.CKF+other.CKD+cancer+
             use.parkinson+  use.rheumatism+
             use.Depression+use.sleep.apnea+
             use.Aspirin+use.ACEI.base+use.ARB.base+
             use.Bblocker.base+
             use.CCB.base+use.diuretics.base+
             statin_base+Insulin_base+other_antiDM,
           x=TRUE, y=TRUE,data=indf)

Pre_HR <-rms::Predict(fit,Triglyceride,
                      fun=exp,type="predictions",
                      ref.zero=T,conf.int = 0.95,digits=2)

ggplot(Pre_HR)

#  Sets the background color of the density curve 
violet <- "#437f9b"
#  Draw the left and right double axes base plot
par(mar = c(5, 4, 4, 4) + 0.3)
par(xpd=NA)

ylim.bot <- min(Pre_HR$lower)
ylim.top <- max(Pre_HR$upper)

par(new=TRUE) #  New canvas 

plot(Pre_HR[,colnames(Pre_HR)=="Triglyceride"],Pre_HR[,colnames(Pre_HR)=="yhat"],
      main = "All-cause mortality", 
     font.main = 1,
      xlab = "Triglyceride (mmol/L)",
     xlab = "TyG",
     ylab = paste0( "All-cause mortality (HR, 95% CI)"), 
     type = "l",
      ylim = c(ylim.bot,ylim.top),
     col="darkslateblue",lwd=2) 

lines(Pre_HR[,colnames(Pre_HR)=="Triglyceride"],Pre_HR[,colnames(Pre_HR)=="lower"],lty=2,lwd=1.5, col="darkslateblue")
lines(Pre_HR[,colnames(Pre_HR)=="Triglyceride"],Pre_HR[,colnames(Pre_HR)=="upper"],lty=2,lwd=1.5, col="darkslateblue")
lines(x=range(Pre_HR[,colnames(Pre_HR)=="Triglyceride"]),y=c(1,1),lty=3,col="black",lwd=1.3) 
#points(as.numeric(refvalue),1,pch=16,cex=1.2)
#  Enclosed refvalue Value , The specific location can be modified by yourself 
#text(refvalue + 1, 0.9, paste0("refvalue = ",refvalue)) 


#  Draw a legend , Pay attention to nonlinearity p The value is in the variable p Position in 
legend("topright",
       paste0("P-overall ",
              ifelse(anova(fit)[1,3] < 0.001,"< 0.001",round(anova(fit)[1,3],3)),
              "\nP-non-linear ",
              ifelse(anova(fit)[2,3] < 0.001,"< 0.001",round(anova(fit)[2,3],3))),
       bty="n",cex=0.8)

legend("topleft",lty = c(1,2),col = c("black","darkslateblue"),
       c("Hazard Ratio","95% CI"),
       bty="n",cex=0.8) 
legend("bottomright",
       paste0("N= 127124",
              "\nEvents= 48751"),
       bty="n",cex=0.8) 

####Posthoc####

library(reshape2)
library(dplyr)
library(rstatix)
library(openxlsx)

TG.data$TG.group.label <- factor(TG.data$TG.group.label, levels = TG_list) #Give order the TG group

posthoc_fun <- function(variable){result <- as.data.frame(pairwise_prop_test(table(data$TG.group.label,
                                                                                   data[[variable]])
                                                                             , p.adjust.method = "bonferroni")) %>%
  arrange(factor(group1, levels = TG_list))
return(result)} #function to perform posthoc analysis


posthoc_output <- lapply(var_list,posthoc_fun) #apply the function to all lists


# Create a new workbook
wb <- createWorkbook()

# Loop through the posthoc_output list and add each dataframe to a new sheet
for (i in 1:length(posthoc_output)) {
  addWorksheet(wb, as.character(var_list[i]))
  writeData(wb, as.character(var_list[i]), posthoc_output[[i]])
}

# Save the workbook to a file
saveWorkbook(wb, "posthoc_output.xlsx", overwrite = TRUE)

#####Cox regression and competing risk####

#  
library(survival)
library(cmprsk)
library(dplyr)

 
# Cox regression
cox_model <- coxph(Surv(time, status == 1) ~ 
                   as.factor(triglyceride) +
                     age + sex +smoking+alcohol+
                     hypertension+diabetes+obesity+stroke+
                     anemia+CAD+PVD+dyslipidemia+AF+cirrhosis+
                     CKF+other_CKD+cancer+rheumatism+depression+
                     sleep.apnea+
                     Aspirin+ACEI+ARB+bblocker+CCB+diuretics+
                     statin+insulin+other.anti.DM+
                     TC+LDL+HDL, 
                   data = data)
summary(cox_model)

 
# competing risk regression
cr_model <- crr(ftime = data$time, fstatus = data$status , 
                cov1 = data.frame(triglyceride = factor(data$triglyceride), 
                                  age = data$age, 
                                  sex = data$sex, 
                                  smoking = data$smoking,
                                  alcohol = data$alcohol,
                                  hypertension= data$hypertension,
                                  diabetes= data$diabetes,
                                  obesity= data$obesity,
                                  stroke= data$stroke,
                                  anemia = data$anemia,
                                  CAD= data$CAD,
                                  PVD= data$PVD,
                                  dyslipidemia= data$dyslipidemia,
                                  AF= data$AF,
                                  cirrhosis= data$cirrhosis,
                                  CKF= data$CKF,
                                  other_CKD= data$other_CKD,
                                  cancer= data$cancer,
                                  rheumatism= data$rheumatism,
                                  depression= data$depression,
                                  sleep.apnea= data$sleep.apnea,
                                  Aspirin= data$Aspirin,
                                  ACEI= data$ACEI,
                                  ARB= data$ARB,
                                  bblocker= data$bblocker,
                                  CCB= data$CCB,
                                  diuretics= data$diuretics,
                                  statin= data$statin,
                                  insulin= data$insulin,
                                  other.anti.DM= data$other.anti.DM,
                                  TC= data$TC,
                                  LDL= data$LDL,
                                  HDL= data$HDL 
                                  ))

summary(cr_model)








