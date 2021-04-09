require(lubridate)
require(ggplot2)
require(cowplot)
require(readxl)
require(dplyr)
require(pscl)
library(glmmTMB)
library(ciTools)
library(penalized)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('../')

moving_average <- function(vec, dist = 3){
  v2 = c()
  for(i in 1:length(vec)){
    v2[i] = mean(vec[(max(1, i - dist)):(min(length(vec), i + dist))])
  }
  return(v2)
}

#
##### Prepping the data - moving average #####
# pull in the data
res_data = data.table::fread('data/comarques_diari.csv', data.table = F)

# reformat data and filter by date
res_data$date<-as.Date(res_data$DATA,"%Y-%m-%d")
res_data <- res_data%>%mutate(year=year(date))
res_data <- res_data %>% filter(date <= as.Date('2021-03-28'))

res_data<-subset(res_data, NOM!="VALL D'ARAN")
comarcas <- unique(res_data$NOM)

#keep all cases in the community or the ones in residences >65 
res_data<-subset(res_data,RESIDENCIA=="No"|GRUP_EDAT!="Entre 15 i 64")
res_data<-subset(res_data,RESIDENCIA=="No"|GRUP_EDAT!="Menors de 15")

# group by week and merge in age cases
df_out<-res_data %>% group_by(date,year,RESIDENCIA,NOM) %>% summarise(cases = sum(CASOS_CONFIRMAT), vacc_1 = sum(VACUNATS_DOSI_1), vacc_2 = sum(VACUNATS_DOSI_2), deaths = sum(EXITUS))

# turn the data to wide with community and nursing home cases per row
df_sub = df_out %>% filter(RESIDENCIA == 'Si') %>%
  mutate(cases_N = cases) %>%  select(-cases)
df_com = df_out %>% filter(RESIDENCIA == 'No') %>%
  mutate(cases_C = cases, vacc_1_C = vacc_1, vacc_2_C = vacc_2) %>% ungroup() %>% select(-cases, - RESIDENCIA, - vacc_1, - vacc_2, -deaths)
df_wide = merge(df_sub, df_com, by = c('date','year','NOM')) %>% select(-RESIDENCIA)
# df_wide = df_wide %>% filter(date >= as.Date('2020-07-06'))

# getting cumulative vaccinations
df_wide = df_wide %>%
  group_by(NOM) %>%
  arrange(date) %>%
  mutate(cum_vacc_1 = cumsum(vacc_1), cum_vacc_2 = cumsum(vacc_2),
         cum_vacc_1_C = cumsum(vacc_1_C), cum_vacc_2_C = cumsum(vacc_2_C))

# getting shifted vaccinations and cases and susceptibility
df_wide2 = NULL
for(n in comarcas){
  tmp = df_wide %>% 
    filter(NOM == n) %>%
    arrange(date)
  
  # get cumulative vaccinations and cases (before doing the moving average)
  tmp$cum_cases_N_m1 = c(0,cumsum(tmp$cases_N)[1:(nrow(tmp)-1)])
  tmp$cum_vacc_1_m1 = c(0, tmp$cum_vacc_1[1:(nrow(tmp)-1)])
  tmp$cum_vacc_2_m1 = c(0, tmp$cum_vacc_2[1:(nrow(tmp)-1)])
  tmp$cum_vacc_1_m2 = c(0, 0, tmp$cum_vacc_1[1:(nrow(tmp)-2)])
  tmp$cum_vacc_2_m2 = c(0, 0, tmp$cum_vacc_2[1:(nrow(tmp)-2)])
  tmp$cum_vacc_1_m3 = c(0, 0, 0, tmp$cum_vacc_1[1:(nrow(tmp)-3)])
  tmp$cum_vacc_2_m3 = c(0, 0, 0, tmp$cum_vacc_2[1:(nrow(tmp)-3)])
  tmp$cum_vacc_1_m4 = c(0, 0, 0, 0, tmp$cum_vacc_1[1:(nrow(tmp)-4)])
  tmp$cum_vacc_2_m4 = c(0, 0, 0, 0, tmp$cum_vacc_2[1:(nrow(tmp)-4)])
  
  # get death moving average
  tmp$deaths_MA = moving_average(tmp$deaths)
  
  # apply moving average to cases and vaccinations
  tmp$cases_N = moving_average(tmp$cases_N)
  tmp$vacc_1 = moving_average(tmp$vacc_1)
  tmp$vacc_2 = moving_average(tmp$vacc_2)
  tmp$cases_C = moving_average(tmp$cases_C)
  tmp$vacc_1_C = moving_average(tmp$vacc_1_C)
  tmp$vacc_2_C = moving_average(tmp$vacc_2_C)
  
  # get lagged vaccinations and casses
  tmp$vacc_1_m1 = c(rep(0,7), tmp$vacc_1[1:(nrow(tmp)-7)])
  tmp$vacc_2_m1 = c(rep(0,7), tmp$vacc_2[1:(nrow(tmp)-7)])
  tmp$vacc_1_m2 = c(rep(0,14), tmp$vacc_1[1:(nrow(tmp)-14)])
  tmp$vacc_2_m2 = c(rep(0,14), tmp$vacc_2[1:(nrow(tmp)-14)])
  tmp$cases_N_m1 = c(rep(0,7), tmp$cases_N[1:(nrow(tmp)-7)])
  tmp$cases_N_m2 = c(rep(0,14), tmp$cases_N[1:(nrow(tmp)-14)])
  tmp$cases_N_m3 = c(rep(0,21), tmp$cases_N[1:(nrow(tmp)-21)])
  tmp$cases_N_m4 = c(rep(0,28), tmp$cases_N[1:(nrow(tmp)-28)])
  tmp$cases_C_m1 = c(rep(0,7),tmp$cases_C[1:(nrow(tmp)-7)])
  tmp$cases_C_m2 = c(rep(0,14), tmp$cases_C[1:(nrow(tmp)-14)])
  
  tmp$max_vacc = tmp$cum_vacc_1[nrow(tmp)]
  
  # a measure of immunity
  tmp$Im = 0.1
  I_old = 0.1
  
  for(i in 2:nrow(tmp)){
    I_new = I_old*(1-0.18/365) - tmp$deaths[i-1]/tmp$max_vacc[1] + tmp$cases_N[i]/tmp$max_vacc[1]
    # bound the value by 0 and 1
    I_old = min(max(I_new, 0),1)
    tmp$Im[i] = I_old
  }
  
  df_wide2 = rbind(df_wide2, tmp)
}

# rename the data frame
df_wide = df_wide2

# filter out old dates
df_wide = df_wide %>% filter(date >= as.Date('2020-07-06'))

# susceptibility
df_wide$S = 1- df_wide$Im

# make training and testing data
train = df_wide %>% filter(date <= as.Date('2020-12-27'))
test = df_wide %>% filter(date > as.Date('2020-12-27'))

# save(df_wide, comarcas, train, test, file = 'data/comarcas_moving_average_03282021.RData')

#
##### Prepping the data - region level moving average #####
#prepare the surveillance database
res_data <- data.table::fread(file = 'data/regio_diari.csv',data.table=F) 

# reformat
res_data$date<-as.Date(res_data$DATA,"%Y-%m-%d")
res_data <- res_data%>%mutate(year=year(date))
res_data = res_data %>% filter(date <= as.Date('2021-03-28'))

res_data<-subset(res_data, NOM!="VALL D'ARAN")
comarcas <- unique(res_data$NOM)

#keep all cases in the community or the ones in residences >65 
res_data<-subset(res_data,RESIDENCIA=="No"|GRUP_EDAT!="Entre 15 i 64")
res_data<-subset(res_data,RESIDENCIA=="No"|GRUP_EDAT!="Menors de 15")

# group by week and merge in age cases
df_out<-res_data %>% group_by(date,year,RESIDENCIA,NOM) %>% summarise(cases = sum(CASOS_CONFIRMAT), vacc_1 = sum(VACUNATS_DOSI_1), vacc_2 = sum(VACUNATS_DOSI_2), deaths = sum(EXITUS))

# turn the data to wide with community and nursing home cases per row
df_sub = df_out %>% filter(RESIDENCIA == 'Si') %>%
  mutate(cases_N = cases) %>%  select(-cases)
df_com = df_out %>% filter(RESIDENCIA == 'No') %>%
  mutate(cases_C = cases, vacc_1_C = vacc_1, vacc_2_C = vacc_2) %>% ungroup() %>% select(-cases, - RESIDENCIA, - vacc_1, - vacc_2, -deaths)
df_wide = merge(df_sub, df_com, by = c('date','year','NOM')) %>% select(-RESIDENCIA)
#df_wide = df_wide %>% filter(date >= as.Date('2020-07-06'))

# getting cumulative vaccinations
df_wide = df_wide %>%
  group_by(NOM) %>%
  arrange(date) %>%
  mutate(cum_vacc_1 = cumsum(vacc_1), cum_vacc_2 = cumsum(vacc_2),
         cum_vacc_1_C = cumsum(vacc_1_C), cum_vacc_2_C = cumsum(vacc_2_C))

# getting shifted vaccinations and cases and "susceptibility"
df_wide2 = NULL
for(n in comarcas){
  tmp = df_wide %>% 
    filter(NOM == n) %>%
    arrange(date)
  
  # get cumulative vaccinations and cases (before doing the moving average)
  tmp$cum_cases_N_m1 = c(0,cumsum(tmp$cases_N)[1:(nrow(tmp)-1)])
  tmp$cum_vacc_1_m1 = c(0, tmp$cum_vacc_1[1:(nrow(tmp)-1)])
  tmp$cum_vacc_2_m1 = c(0, tmp$cum_vacc_2[1:(nrow(tmp)-1)])
  tmp$cum_vacc_1_m2 = c(0, 0, tmp$cum_vacc_1[1:(nrow(tmp)-2)])
  tmp$cum_vacc_2_m2 = c(0, 0, tmp$cum_vacc_2[1:(nrow(tmp)-2)])
  tmp$cum_vacc_1_m3 = c(0, 0, 0, tmp$cum_vacc_1[1:(nrow(tmp)-3)])
  tmp$cum_vacc_2_m3 = c(0, 0, 0, tmp$cum_vacc_2[1:(nrow(tmp)-3)])
  tmp$cum_vacc_1_m4 = c(0, 0, 0, 0, tmp$cum_vacc_1[1:(nrow(tmp)-4)])
  tmp$cum_vacc_2_m4 = c(0, 0, 0, 0, tmp$cum_vacc_2[1:(nrow(tmp)-4)])
  
  # get death moving average
  tmp$deaths_MA = moving_average(tmp$deaths)
  
  # apply moving average to cases and vaccinations
  tmp$cases_N = moving_average(tmp$cases_N)
  tmp$vacc_1 = moving_average(tmp$vacc_1)
  tmp$vacc_2 = moving_average(tmp$vacc_2)
  tmp$cases_C = moving_average(tmp$cases_C)
  tmp$vacc_1_C = moving_average(tmp$vacc_1_C)
  tmp$vacc_2_C = moving_average(tmp$vacc_2_C)
  
  # get lagged values of vaccinations and cases
  tmp$vacc_1_m1 = c(rep(0,7), tmp$vacc_1[1:(nrow(tmp)-7)])
  tmp$vacc_2_m1 = c(rep(0,7), tmp$vacc_2[1:(nrow(tmp)-7)])
  tmp$vacc_1_m2 = c(rep(0,14), tmp$vacc_1[1:(nrow(tmp)-14)])
  tmp$vacc_2_m2 = c(rep(0,14), tmp$vacc_2[1:(nrow(tmp)-14)])
  tmp$cases_N_m1 = c(rep(0,7), tmp$cases_N[1:(nrow(tmp)-7)])
  tmp$cases_N_m2 = c(rep(0,14), tmp$cases_N[1:(nrow(tmp)-14)])
  tmp$cases_N_m3 = c(rep(0,21), tmp$cases_N[1:(nrow(tmp)-21)])
  tmp$cases_N_m4 = c(rep(0,28), tmp$cases_N[1:(nrow(tmp)-28)])
  tmp$cases_C_m1 = c(rep(0,7),tmp$cases_C[1:(nrow(tmp)-7)])
  tmp$cases_C_m2 = c(rep(0,14), tmp$cases_C[1:(nrow(tmp)-14)])
  tmp$cases_C_m3 = c(rep(0,21), tmp$cases_C[1:(nrow(tmp)-21)])
  tmp$cases_C_m4 = c(rep(0,28), tmp$cases_C[1:(nrow(tmp)-28)])
  
  tmp$max_vacc = tmp$cum_vacc_1[nrow(tmp)]
  
  # a measure of immunity
  tmp$Im = 0.1
  I_old = 0.1
  
  for(i in 2:nrow(tmp)){
    I_new = I_old*(1-0.18/365) - tmp$deaths[i-1]/tmp$max_vacc[1] + tmp$cases_N[i]/tmp$max_vacc[1]
    # bound the value by 0 and 1
    I_old = min(max(I_new, 0),1)
    tmp$Im[i] = I_old
  }
  tmp = tmp %>% arrange(date)
  
  df_wide2 = rbind(df_wide2, tmp)
}

# reset the named data frame
df_wide = df_wide2 %>% filter(date >= as.Date('2020-07-06'))

# susceptibility
df_wide$S = 1- df_wide$Im

# make training and testing data
train = df_wide %>% filter(date <= as.Date('2020-12-27'))
test = df_wide %>% filter(date > as.Date('2020-12-27'))

# save(df_wide, comarcas, train, test, file = 'data/region_moving_average_03282021.RData')

#
##### Prepping the data - weekly #####
res_data = data.table::fread('data/comarques_diari.csv', data.table = F)

# reformat date
res_data$date<-as.Date(res_data$DATA,"%Y-%m-%d")
res_data <- res_data%>%mutate(week=isoweek(date),year=year(date))
res_data = res_data %>% filter(date <= as.Date('2021-03-28'))

res_data<-subset(res_data, NOM!="VALL D'ARAN")
comarcas <- unique(res_data$NOM)

#keep all cases in the community or the ones in residences >65 
res_data<-subset(res_data,RESIDENCIA=="No"|GRUP_EDAT!="Entre 15 i 64")
res_data<-subset(res_data,RESIDENCIA=="No"|GRUP_EDAT!="Menors de 15")

# group by week and merge in age cases
df_out<-res_data %>% group_by(week,year,RESIDENCIA,NOM) %>% summarise(cases = sum(CASOS_CONFIRMAT), vacc_1 = sum(VACUNATS_DOSI_1), vacc_2 = sum(VACUNATS_DOSI_2), deaths = sum(EXITUS))
df_out <- merge(df_out, age_cases, by=c('week','year','RESIDENCIA','NOM'), all.x = T)

# correct the problem of duplicated week 53
for (i in comarcas){
  
  # fix cases
  {
  df_out$cases[which(df_out$week==53 & df_out$year==2020
                     &df_out$RESIDENCIA=="No"&
                       df_out$NOM==i)]<-
    sum(df_out$cases[which(df_out$week==53 &
                             df_out$RESIDENCIA=="No"&
                             df_out$NOM==i)])

  df_out$cases[which(df_out$week==53 & df_out$year==2020
                     &df_out$RESIDENCIA=="Si"&
                       df_out$NOM==i)]<-
    sum(df_out$cases[which(df_out$week==53 &
                             df_out$RESIDENCIA=="Si"&
                             df_out$NOM==i)])
  }
  
  # fix vacc_1
  {
    df_out$vacc_1[which(df_out$week==53 & df_out$year==2020
                       &df_out$RESIDENCIA=="No"&
                         df_out$NOM==i)]<-
      sum(df_out$vacc_1[which(df_out$week==53 &
                               df_out$RESIDENCIA=="No"&
                               df_out$NOM==i)])
    
    df_out$vacc_1[which(df_out$week==53 & df_out$year==2020
                       &df_out$RESIDENCIA=="Si"&
                         df_out$NOM==i)]<-
      sum(df_out$vacc_1[which(df_out$week==53 &
                               df_out$RESIDENCIA=="Si"&
                               df_out$NOM==i)])
  }

  # fix vacc_2
  {
    df_out$vacc_2[which(df_out$week==53 & df_out$year==2020
                       &df_out$RESIDENCIA=="No"&
                         df_out$NOM==i)]<-
      sum(df_out$vacc_2[which(df_out$week==53 &
                               df_out$RESIDENCIA=="No"&
                               df_out$NOM==i)])
    
    df_out$vacc_2[which(df_out$week==53 & df_out$year==2020
                       &df_out$RESIDENCIA=="Si"&
                         df_out$NOM==i)]<-
      sum(df_out$vacc_2[which(df_out$week==53 &
                               df_out$RESIDENCIA=="Si"&
                               df_out$NOM==i)])
  }
}
df_out<-subset(df_out,week<=52 | year<=2020)
df_out$week[which(df_out$year>=2021)]<-df_out$week[which(df_out$year>=2021)]+53

# turn the data to wide with community and nursing home cases per row
df_sub = df_out %>% filter(RESIDENCIA == 'Si') %>%
  mutate(cases_N = cases) %>%  select(-cases, - RESIDENCIA)

df_com = df_out %>% filter(RESIDENCIA == 'No') %>%
  mutate(cases_C = cases, vacc_1_C = vacc_1, vacc_2_C = vacc_2) %>% ungroup() %>% select(-cases, - RESIDENCIA, - vacc_1, - vacc_2, -deaths)
df_wide = merge(df_sub, df_com, by = c('week','year','NOM'))
df_wide = df_wide %>% filter(week >= 28)

# getting cumulative vaccinations
df_wide = df_wide %>%
  group_by(NOM) %>%
  arrange(week) %>%
  mutate(cum_vacc_1 = cumsum(vacc_1), cum_vacc_2 = cumsum(vacc_2),
         cum_vacc_1_C = cumsum(vacc_1_C), cum_vacc_2_C = cumsum(vacc_2_C))

# getting shifted vaccinations and cases and "susceptibility"
df_wide2 = NULL
for(n in comarcas){
  tmp = df_wide %>% 
    filter(NOM == n) %>%
    arrange(week)
  
  tmp$cum_vacc_1_m1 = c(0, tmp$cum_vacc_1[1:(nrow(tmp)-1)])
  tmp$cum_vacc_2_m1 = c(0, tmp$cum_vacc_2[1:(nrow(tmp)-1)])
  
  tmp$cum_vacc_1_m2 = c(0, 0, tmp$cum_vacc_1[1:(nrow(tmp)-2)])
  tmp$cum_vacc_2_m2 = c(0, 0, tmp$cum_vacc_2[1:(nrow(tmp)-2)])
  
  tmp$cum_vacc_1_m3 = c(0, 0, 0, tmp$cum_vacc_1[1:(nrow(tmp)-3)])
  tmp$cum_vacc_2_m3 = c(0, 0, 0, tmp$cum_vacc_2[1:(nrow(tmp)-3)])
  
  tmp$cum_vacc_1_m4 = c(0, 0, 0, 0, tmp$cum_vacc_1[1:(nrow(tmp)-4)])
  tmp$cum_vacc_2_m4 = c(0, 0, 0, 0, tmp$cum_vacc_2[1:(nrow(tmp)-4)])
  
  tmp$cases_C_m1 = c(0,tmp$cases_C[1:(nrow(tmp)-1)])
  tmp$cases_C_m2 = c(0,0, tmp$cases_C[1:(nrow(tmp)-2)])
  
  tmp$max_vacc = tmp$cum_vacc_1[nrow(tmp)]

  # a measure of immunity
  tmp$Im = 0.1
  I_old = 0.1
  
  for(i in 2:nrow(tmp)){
    I_new = I_old*(1-0.18/52) - tmp$deaths[i-1]/tmp$max_vacc[1] + tmp$cases_N[i]/tmp$max_vacc[1]
    # bound the value by 0 and 1
    I_old = min(max(I_new, 0),1)
    tmp$Im[i] = I_old
  }
  
  # a scaling factor
  tmp$delta = tmp$cases_C/tmp$cases_C_m1
  tmp$delta[is.infinite(tmp$delta)] = 2 # (a number/0)
  tmp$delta[is.nan(tmp$delta)] = 1 # (0/0)

  df_wide2 = rbind(df_wide2, tmp)
}

# rename
df_wide = df_wide2

# adding in factor for week (when using autocorrelation)
df_wide$week.f = factor(df_wide$week, levels = sort(unique(df_wide$week)))

# susceptibility
df_wide$S = 1- df_wide$Im

# make training and testing data
train = df_wide %>% filter(week <= 52)
test = df_wide %>% filter(week >= 53)

# save(df_wide, comarcas, train, test, file = 'data/comarcas_data_03282021.RData')


##### Prepping the data - region level weekly #####
res_data <- data.table::fread(file = 'data/regio_diari.csv',data.table=F) 

# reformat date
res_data$date<-as.Date(res_data$DATA,"%Y-%m-%d")
res_data <- res_data%>%mutate(week=isoweek(date),year=year(date))
res_data = res_data %>% filter(date <= as.Date('2021-03-28'))

res_data<-subset(res_data, NOM!="VALL D'ARAN")
comarcas <- unique(res_data$NOM)

#keep all cases in the community or the ones in residences >65 
res_data<-subset(res_data,RESIDENCIA=="No"|GRUP_EDAT!="Entre 15 i 64")
res_data<-subset(res_data,RESIDENCIA=="No"|GRUP_EDAT!="Menors de 15")

# group by week and merge in age cases
df_out<-res_data %>% group_by(week,year,RESIDENCIA,NOM) %>% summarise(cases = sum(CASOS_CONFIRMAT), vacc_1 = sum(VACUNATS_DOSI_1), vacc_2 = sum(VACUNATS_DOSI_2), deaths = sum(EXITUS))
df_out <- merge(df_out, age_cases, by=c('week','year','RESIDENCIA','NOM'), all.x = T)

# correct the problem of duplicated week 53
for (i in comarcas){
  # fix cases
  {
    df_out$cases[which(df_out$week==53 & df_out$year==2020
                       &df_out$RESIDENCIA=="No"&
                         df_out$NOM==i)]<-
      sum(df_out$cases[which(df_out$week==53 &
                               df_out$RESIDENCIA=="No"&
                               df_out$NOM==i)])
    
    df_out$cases[which(df_out$week==53 & df_out$year==2020
                       &df_out$RESIDENCIA=="Si"&
                         df_out$NOM==i)]<-
      sum(df_out$cases[which(df_out$week==53 &
                               df_out$RESIDENCIA=="Si"&
                               df_out$NOM==i)])
  }
  
  # fix vacc_1
  {
    df_out$vacc_1[which(df_out$week==53 & df_out$year==2020
                        &df_out$RESIDENCIA=="No"&
                          df_out$NOM==i)]<-
      sum(df_out$vacc_1[which(df_out$week==53 &
                                df_out$RESIDENCIA=="No"&
                                df_out$NOM==i)])
    
    df_out$vacc_1[which(df_out$week==53 & df_out$year==2020
                        &df_out$RESIDENCIA=="Si"&
                          df_out$NOM==i)]<-
      sum(df_out$vacc_1[which(df_out$week==53 &
                                df_out$RESIDENCIA=="Si"&
                                df_out$NOM==i)])
  }
  
  # fix vacc_2
  {
    df_out$vacc_2[which(df_out$week==53 & df_out$year==2020
                        &df_out$RESIDENCIA=="No"&
                          df_out$NOM==i)]<-
      sum(df_out$vacc_2[which(df_out$week==53 &
                                df_out$RESIDENCIA=="No"&
                                df_out$NOM==i)])
    
    df_out$vacc_2[which(df_out$week==53 & df_out$year==2020
                        &df_out$RESIDENCIA=="Si"&
                          df_out$NOM==i)]<-
      sum(df_out$vacc_2[which(df_out$week==53 &
                                df_out$RESIDENCIA=="Si"&
                                df_out$NOM==i)])
  }
}
df_out<-subset(df_out,week<=52 | year<=2020)
df_out$week[which(df_out$year>=2021)]<-df_out$week[which(df_out$year>=2021)]+53

# turn the data to wide with community and nursing home cases per row
df_sub = df_out %>% filter(RESIDENCIA == 'Si') %>%
  mutate(cases_N = cases) %>%  select(-cases, - RESIDENCIA)
df_com = df_out %>% filter(RESIDENCIA == 'No') %>%
  mutate(cases_C = cases, vacc_1_C = vacc_1, vacc_2_C = vacc_2) %>% ungroup() %>% select(-cases, - RESIDENCIA, - vacc_1, - vacc_2, -deaths)
df_wide = merge(df_sub, df_com, by = c('week','year','NOM'))
df_wide = df_wide %>% filter(week >= 28)

# getting cumulative vaccinations
df_wide = df_wide %>%
  group_by(NOM) %>%
  arrange(week) %>%
  mutate(cum_vacc_1 = cumsum(vacc_1), cum_vacc_2 = cumsum(vacc_2),
         cum_vacc_1_C = cumsum(vacc_1_C), cum_vacc_2_C = cumsum(vacc_2_C))

# getting shifted vaccinations and cases and "susceptibility"
df_wide2 = NULL
for(n in comarcas){
  tmp = df_wide %>% 
    filter(NOM == n) %>%
    arrange(week)
  
  tmp$cum_vacc_1_m1 = c(0, tmp$cum_vacc_1[1:(nrow(tmp)-1)])
  tmp$cum_vacc_2_m1 = c(0, tmp$cum_vacc_2[1:(nrow(tmp)-1)])
  
  tmp$cum_vacc_1_m2 = c(0, 0, tmp$cum_vacc_1[1:(nrow(tmp)-2)])
  tmp$cum_vacc_2_m2 = c(0, 0, tmp$cum_vacc_2[1:(nrow(tmp)-2)])
  
  tmp$cum_vacc_1_m3 = c(0, 0, 0, tmp$cum_vacc_1[1:(nrow(tmp)-3)])
  tmp$cum_vacc_2_m3 = c(0, 0, 0, tmp$cum_vacc_2[1:(nrow(tmp)-3)])
  
  tmp$cum_vacc_1_m4 = c(0, 0, 0, 0, tmp$cum_vacc_1[1:(nrow(tmp)-4)])
  tmp$cum_vacc_2_m4 = c(0, 0, 0, 0, tmp$cum_vacc_2[1:(nrow(tmp)-4)])
  
  tmp$cases_C_m1 = c(0,tmp$cases_C[1:(nrow(tmp)-1)])
  tmp$cases_C_m2 = c(0,0, tmp$cases_C[1:(nrow(tmp)-2)])
  
  tmp$cases_N_m1 = c(0,tmp$cases_N[1:(nrow(tmp)-1)])
  tmp$cases_N_m2 = c(0,0, tmp$cases_N[1:(nrow(tmp)-2)])
  tmp$cases_N_m3 = c(0,0,0, tmp$cases_N[1:(nrow(tmp)-3)])
  
  tmp$max_vacc = tmp$cum_vacc_1[nrow(tmp)]
  
  # a measure of immunity
  tmp$Im = 0.1
  I_old = 0.1
  
  for(i in 2:nrow(tmp)){
    I_new = I_old*(1-0.18/52) - tmp$deaths[i-1]/tmp$max_vacc[1] + tmp$cases_N[i]/tmp$max_vacc[1]
    I_old = max(I_new, 0)
    tmp$Im[i] = I_old
  }
  
  # a scaling factor
  tmp$delta = tmp$cases_C/tmp$cases_C_m1
  tmp$delta[is.infinite(tmp$delta)] = 2 # (a number/0)
  tmp$delta[is.nan(tmp$delta)] = 1 # (0/0)
  
  df_wide2 = rbind(df_wide2, tmp)
}

# rename
df_wide = df_wide2

# wack measure of susceptibility
df_wide$S = 1- df_wide$Im

# make training and testing data
train = df_wide %>% filter(week <= 52)
test = df_wide %>% filter(week >= 53)

# save(df_wide, train, test, file = 'data/region_level_df_wide_03282021.RData')


