#Chicken vs. human choice landing assays: Data Analysis
#Noreuil and Fritz updated 12/20/2020

#GOAL: Investigate variation in host preference across eight AG and BG pipiens populations 
#using a mixed logistic regression model with a binomial error structure. 

#Assess if three types of host response vary according to mosquito population: 
#(1) acceptance of either host during day 1
#(2) selection/preference for a human host
#(3) probability to alternate host across multi-day testing.


# Bootstrapping function --------------------------------------------------
boot.fn <- function(x, N=5000) {
  Int.1 <- replicate(N, mean(sample(x, size= length(x), replace=T)))
  Int.CI <- quantile(Int.1, probs=c(0.025,0.975))
  Int.CI
}


#Loading libraries
library(lme4); library(lmtest)

#Loading datasets
#full_ba_set <- read.csv("~/Desktop/Noreuil_Fritz2020/Data_BehavioralAssay.csv")#Anna
full_ba_set <- read.csv("~/Downloads/Data_BehavioralAssay.csv")#Megan
full_ba_set$time <- as.factor(full_ba_set$time)
str(full_ba_set)
head(full_ba_set)
tail(full_ba_set)

#Getting 'chick' ID to use as a random effect
full_ba_set$chick <- paste(full_ba_set$chick.color,full_ba_set$chick.no)
full_ba_set$chick <- as.factor(full_ba_set$chick)
str(full_ba_set)

#Getting new mosquito IDs because 2 from diff strains have same ID
full_ba_set$new.mosq.ID <- paste(full_ba_set$mosq.id,full_ba_set$strain)
full_ba_set$new.mosq.ID <- as.factor(full_ba_set$new.mosq.ID)

#Dropping old mosq.id
full_ba_set <- full_ba_set[,-1]
str(full_ba_set)

#####################################################################
#(1) acceptance of either host during day 1
#####################################################################
#recode host response - all nr values as 0 and c/h as a 1
full_ba_set$host_resp<- ifelse(full_ba_set$host == "nr", 0, 1)

T1_ba_set <- subset(full_ba_set, time == 1)
table(T1_ba_set$strain)#getting numbers of individs tested at time1

resp_only <- subset(T1_ba_set, host_resp == 1)
table(resp_only$strain)#getting numbers of responders tested at time1


#No'time' in model since all testing for this analysis took place on day 1
#Maximal model - inclusion of all fixed effects and interactions between them 
#reduced model - removal of strain*time interaction term 

model_full <- glmer(host_resp ~ 1 + strain + (1|chick), data = T1_ba_set, 
                    family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

model_red1 <- glmer(host_resp ~ 1 + (1|chick), data = T1_ba_set, 
                    family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

lrtest(model_red1, model_full)


#####Need a quick test for Cal1_2wk v Cal1_3wk
Cals_only <- subset(T1_ba_set, strain == "cal1" | strain == "cal1_2w")

#did response rate differ?
model_cals <- glmer(host_resp ~ 1 + strain + (1|chick), data = Cals_only, 
                                  family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

model_cals_red <- glmer(host_resp ~ 1 + (1|chick), data = Cals_only, 
                    family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

lrtest(model_cals_red, model_cals)#not different - see updated model for T1 below



#####Rerunning model with cal1 com 2 & 3 wk combined
T1_ba_set$strain <- gsub("_2w", "", T1_ba_set$strain)
T1_ba_set$strain <- as.factor(T1_ba_set$strain)
levels(T1_ba_set$strain)#sanity check


###########Anna - use these for final 
model_full <- glmer(host_resp ~ 1 + strain + (1|chick), data = T1_ba_set, 
                    family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

model_red1 <- glmer(host_resp ~ 1 + (1|chick), data = T1_ba_set, 
                    family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

lrtest(model_red1, model_full)


mean_day1_resp <- tapply(T1_ba_set$host_resp, T1_ba_set$strain, mean)
mean_day1_resp

boot_day1_resp <- tapply(T1_ba_set$host_resp, T1_ba_set$strain, boot.fn)
boot_day1_resp

#####################################################################
#(2) selection/preference for a human host
#####################################################################
#NOTE: sweep here then re-run code through setting up chick as a random effect (line 31)

#Setting up the dataset
#Host choice (i.e. chick vs. human) was examined for individuals across test dates.  
#data for individuals that failed to respond on multiple days were removed from the dataset (all 'nr' values removed)
resp_only_full_ba_set <- subset(full_ba_set, host != "nr")
head(resp_only_full_ba_set)
sum(resp_only_full_ba_set$host!= "nr")

#For this model, a host selection of the ‘human’ host was given a value of 1, and selection of a ‘chick’ host a 0.
#recoding the host choices to make binary variable
resp_only_full_ba_set$host_hum <- ifelse(resp_only_full_ba_set$host == "c", 0, 1)
resp_only_full_ba_set$host_chick <- ifelse(resp_only_full_ba_set$host == "c", 1, 0)

#Because host response was measured across the same individuals over the multi-day testing period, 
#‘day’ (time) was also included as a second random effect in this model 
#Necessary to account for autocorrelation, as we are measuring host response across the same individual over time

model_full_final <- glmer(host_hum ~ 1 + strain + time + (1|chick), data = resp_only_full_ba_set, 
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

model_red <- glmer(host_hum ~ 1 + time + (1|chick), data = resp_only_full_ba_set, 
                   family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

lrtest(model_red, model_full_final) 


#####Need a quick test for Cal1_2wk v Cal1_3wk
Cals_only_multi <- subset(resp_only_full_ba_set, strain == "cal1" | strain == "cal1_2w")

#did response rate differ?
model_cals <- glm(host_resp ~ 1 + strain + time, data = Cals_only_multi, 
                    family = binomial, maxit = 500)

model_cals_red <- glm(host_resp ~ 1 + time, data = Cals_only_multi, 
                      family = binomial, maxit = 500)

summary(model_cals)#don't expect problems, but checking to see that the model coefficients differ here
summary(model_cals_red)

lrtest(model_cals_red, model_cals)#not diff - see updated model below



#####Rerunning model with cal1 com 2 & 3 wk combined
resp_only_full_ba_set$strain <- gsub("_2w", "", resp_only_full_ba_set$strain)
resp_only_full_ba_set$strain <- as.factor(resp_only_full_ba_set$strain)
levels(resp_only_full_ba_set$strain)#sanity check

str(resp_only_full_ba_set)


#########Anna - use these for final
model_full_final <- glmer(host_hum ~ 1 + strain + time + (1|chick), data = resp_only_full_ba_set, 
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

model_red <- glmer(host_hum ~ 1 + time + (1|chick), data = resp_only_full_ba_set, 
                   family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

lrtest(model_red, model_full_final) 


mean_host_hum <- tapply(resp_only_full_ba_set$host_hum, resp_only_full_ba_set$strain, mean)
mean_host_hum

boot_host_hum <- tapply(resp_only_full_ba_set$host_hum, resp_only_full_ba_set$strain, boot.fn)
boot_host_hum


#####################################################################
##(3) probability to alternate host across multi-day testing.
#####################################################################

#Setting up dataset:

#getting unique mosquito IDs of multi-day responders
length(unique(resp_only_full_ba_set$new.mosq.ID)) #All unique mosquito IDs
uniq_resp <- data.frame(unique(resp_only_full_ba_set$new.mosq.ID))
str(uniq_resp)
colnames(uniq_resp) <- "new.mosq.ID"

uniq_resp$new.mosq.ID <- as.character(uniq_resp$new.mosq.ID)

#producing merged dataset with only responders, but including "nr" datapoint.
host_switch_data <- full_ba_set[,c(1,5,6,9,10)] #new dataframe with relevant columns
nrow(host_switch_data)

#####Running model with cal1 2 & 3 wk combined
host_switch_data$strain <- gsub("_2w", "", host_switch_data$strain)
host_switch_data$strain <- as.factor(host_switch_data$strain)
levels(host_switch_data$strain)#sanity check

str(host_switch_data)


host_switch_data$new.mosq.ID <- as.character(host_switch_data$new.mosq.ID)

merged_host_switch_data <- merge(host_switch_data, uniq_resp, by = "new.mosq.ID")
str(merged_host_switch_data)

ID_counts <- table(merged_host_switch_data$new.mosq.ID)
ID_counts_log <- ID_counts > 3
which(ID_counts_log == TRUE)#no repeated IDs = good


#Reshaping dataset long to wide
ord_host_switch_data <- merged_host_switch_data[order(merged_host_switch_data[,1], merged_host_switch_data[,4]),]

host_switch_data_reshaped <- reshape(ord_host_switch_data, idvar ="new.mosq.ID", timevar = "time", direction = "wide")
host_switch_data_reshaped <- host_switch_data_reshaped[,c(1,2,3,4,5,8)]

host_switch_data_reshaped$host.1 <- as.character(host_switch_data_reshaped$host.1)#need to be characters for NA function to work
host_switch_data_reshaped$host.2 <- as.character(host_switch_data_reshaped$host.2)
host_switch_data_reshaped$host.3 <- as.character(host_switch_data_reshaped$host.3)
str(host_switch_data_reshaped)
summary(host_switch_data_reshaped)

#Recoding responses as numeric with non-responders and NAs = 0
host_switch_data_reshaped[is.na(host_switch_data_reshaped)] <- 0
head(host_switch_data_reshaped)
tail(host_switch_data_reshaped)

host_switch_data_reshaped$host.1 <- ifelse(host_switch_data_reshaped$host.1 == "c", 1, ifelse(host_switch_data_reshaped$host.1 == "h", 4, 0))#human is 4
host_switch_data_reshaped$host.2 <- ifelse(host_switch_data_reshaped$host.2 == "c", 1, ifelse(host_switch_data_reshaped$host.2 == "h", 4, 0))#human is 4
host_switch_data_reshaped$host.3 <- ifelse(host_switch_data_reshaped$host.3 == "c", 1, ifelse(host_switch_data_reshaped$host.3 == "h", 4, 0))#human is 4

#sanity check for variable recode
head(host_switch_data_reshaped)
tail(host_switch_data_reshaped)#looks good


#Adding host switch column - first finding the sum of 2 or 3 responses per female
#if sum is 1 or 4, exclude
#if sum is 5, 6, 9 it's a switcher 

host_switch_data_reshaped$row_sum <-rowSums(host_switch_data_reshaped[,c("host.1", "host.2", "host.3")], na.rm=TRUE)
str(host_switch_data_reshaped)

#dropping single responders
host_switch_data_dropped <- subset(host_switch_data_reshaped, row_sum != 1 & row_sum != 4)

host_switch_data_dropped$host_switch <- ifelse(host_switch_data_dropped$row_sum == 5 | host_switch_data_dropped$row_sum == 6 | host_switch_data_dropped$row_sum == 9, 1, 0)  

table(host_switch_data_dropped$strain.1)


####now for the host-switching models
####had to drop chick as a random effect because we did not have enough observations - only 244 actually responded multiple times.
model_full <- glm(host_switch ~ 1 + strain.1 , data = host_switch_data_dropped, 
                    family = binomial)

model_red1 <- glm(host_switch ~ 1, data = host_switch_data_dropped, 
                  family = binomial)

lrtest(model_red1, model_full)


#getting proportions that host-switched
mean_host_switch <- tapply(host_switch_data_dropped$host_switch, host_switch_data_dropped$strain.1, mean)
mean_host_switch

boot_host_switch <- tapply(host_switch_data_dropped$host_switch, host_switch_data_dropped$strain.1, boot.fn)
boot_host_switch
