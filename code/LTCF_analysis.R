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


##### Get aggregated data #####
load('data/comarcas_moving_average_03282021.RData')

full_data = df_wide %>% 
  group_by(date) %>% 
  summarize(n_C = sum(cases_C), n_N = sum(cases_N), vacc_1 = sum(vacc_1), vacc_2 = sum(vacc_2)) %>%
  ungroup() %>%
  mutate(vacc_1 = cumsum(vacc_1), vacc_2 = cumsum(vacc_2))

n = nrow(full_data)
full_data$n_N_d = c(0, full_data$n_N[2:n] - full_data$n_N[1:(n-1)])
full_data$n_C_d = c(0, full_data$n_C[2:n] - full_data$n_C[1:(n-1)])
full_data$n_N_p = c(0, 100*(full_data$n_N[2:n] - full_data$n_N[1:(n-1)])/full_data$n_N[1:(n-1)])
full_data$n_C_p = c(0, 100*(full_data$n_C[2:n] - full_data$n_C[1:(n-1)])/full_data$n_C[1:(n-1)])


### save the data
full_data = full_data %>% 
  mutate(cases_community = n_C, cases_nursing = n_N, prop_vacc_1 = vacc_1/59448, prop_vacc_2 = vacc_2/59448) %>%
  select(date, cases_community, cases_nursing, vacc_1, vacc_2)

# write.csv(full_data, row.names = F, file = 'results/cases_and_vaccinations_03282021.csv')

#
##### All of Catalon: cases, deaths, mortality. #####
load('data/region_moving_average_03282021.RData')

# roll the data up into the whole region
tt = df_wide %>% select(date, cases_N, cases_N_m1, cases_N_m2, cases_N_m3, cases_C, cases_C_m1, cases_C_m2, cases_C_m3, cases_C_m4, S, max_vacc, vacc_1, vacc_2, deaths, deaths_MA) %>%
  group_by(date) %>%
  summarize(deaths_MA = sum(deaths_MA),
            deaths = sum(deaths),
            cases_N = sum(cases_N), 
            cases_N_m1 = sum(cases_N_m1), 
            cases_N_m2 = sum(cases_N_m2), 
            cases_N_m3 = sum(cases_N_m3), 
            cases_C = sum(cases_C),
            cases_C_m1 = sum(cases_C_m1),
            cases_C_m2 = sum(cases_C_m2),
            cases_C_m3 = sum(cases_C_m3),
            cases_C_m4 = sum(cases_C_m4),
            vacc_1 = sum(vacc_1),
            vacc_2 = sum(vacc_2),
            S_N = sum(S*max_vacc),
            N = sum(max_vacc)) %>%
  mutate(S = S_N/N,
         C_log = log1p(cases_C),
         C_m1_log = log1p(cases_C_m1)) %>%
  arrange(date)
tt$cum_vacc_1 = cumsum(tt$vacc_1)
tt$cum_vacc_2 = cumsum(tt$vacc_2)
tt$cases_N2 = tt$cases_N/tt$S

# get the training data
train_tt = tt %>% filter(date <= '2020-12-27')

R = 5000
# fit the susceptibility offset model
# 
if(TRUE){
model_fit = MASS::glm.nb(cases_N ~ C_log + C_m1_log + offset(log(S)), data = train_tt, control = glm.control(maxit = 5000))
AIC(model_fit)

# Get model results for bootstrapping
beta_hat <- model_fit$coefficients
beta_vcov <- vcov(model_fit)
theta = summary(model_fit)$theta

# run the parametric bootstrap
sim.boot <- sapply(1:R, function(r){
  
  #indicator 
  beta_boot <- MASS::mvrnorm(1,beta_hat,beta_vcov)
  pred_boot <- (tt %>% 
                  mutate(intercept=1) %>%
                  dplyr::select(intercept, C_log, C_m1_log) %>%
                  as.matrix())%*%as.matrix(beta_boot)
  
  pred_boot_exp <- exp(pred_boot + log(tt$S)) 
  x = MASS::rnegbin(n = nrow(tt), mu = pred_boot_exp, theta = theta)
  
  x # UPDATE
  
}) 

tt$LPB0.025 <- apply(sim.boot,1,quantile,.025)
tt$LPB0.05 <- apply(sim.boot,1,quantile,.05)
tt$LPB0.25 <- apply(sim.boot,1,quantile,.25)
tt$UPB0.75 <- apply(sim.boot,1,quantile,.75)
tt$UPB0.95 <- apply(sim.boot,1,quantile,.95)
tt$UPB0.975 <- apply(sim.boot,1,quantile,.975)
tt$median <- apply(sim.boot,1,quantile,.50,na.rm=TRUE)
tt$pred = predict(model_fit, tt, type = 'response')

# plot it!
p1 <- ggplot(tt, aes(x = date, y = cases_N)) +
  geom_line() +
  geom_line(aes(y = pred), col = 'red') + 
  geom_line(aes(y = median), col = 'yellow') + 
  geom_line(aes(y = LPB0.25), col = 'blue') +
  geom_line(aes(y = UPB0.75), col = 'blue') +
  geom_line(aes(y = LPB0.05), col = 'purple') +
  geom_line(aes(y = UPB0.95), col = 'purple') +
  #geom_line(aes(y = LPB0.495), col = 'yellow') +
  ggtitle('case predictions') + 
  theme_minimal()
p1
}

# fit the no-offset model
if(FALSE){
model_fit = MASS::glm.nb(cases_N ~ log1p(cases_C) + log1p(cases_C_m1), data = train_tt, control = glm.control(maxit = 5000))

tt <- ciTools::add_pi(df = tt, fit = model_fit)
tt <- ciTools::add_pi(df = tt, fit = model_fit, alpha = 0.1)
tt <- ciTools::add_pi(df = tt, fit = model_fit, alpha = 0.5)

p1 <- ggplot(tt, aes(x = date, y = cases_N)) +
  geom_line() +
  geom_line(aes(y = pred), col = 'red') + 
  #geom_line(aes(y = median), col = 'yellow') + 
  geom_line(aes(y = LPB0.25), col = 'blue') +
  geom_line(aes(y = UPB0.75), col = 'blue') +
  geom_line(aes(y = LPB0.05), col = 'purple') +
  geom_line(aes(y = UPB0.95), col = 'purple') +
  #geom_line(aes(y = LPB0.495), col = 'yellow') +
  ggtitle('case predictions') + 
  theme_minimal()
p1
}

#date_cutoff = tt$date[min(which(tt$cum_vacc_1/tt$N >= 0.7))] # '2021-01-14'
date_cutoff = tt$date[min(which(tt$cum_vacc_2/tt$N >= 0.7))] # '2021-02-06'

# get the cases averted by the date cutoff
res = tt %>%
  filter(date >= date_cutoff) %>%
  mutate(diff = pred - cases_N, diff_U = UPB0.95 - cases_N, diff_L = LPB0.05 - cases_N) %>% 
  select(date, cases_N, diff, diff_L, diff_U, pred, LPB0.05, UPB0.95, LPB0.25, UPB0.75) %>%
  summarize(cases_averted = sum(diff), lower = sum(diff_L), upper = sum(diff_U), pred = sum(pred), pred_lower = sum(LPB0.05), pred_lower25 = sum(LPB0.25), pred_upper = sum(UPB0.95), pred_upper25 = sum(UPB0.75), cases_N = sum(cases_N))

res
# cutoff 1: 1659 (-690, 4817)
# cutoff 2: 1371 (157, 2867)

# get the  % reduction = 100(1-(E-O)/E)
100*(res$pred - res$cases_N)/res$pred
100*(res$pred_upper - res$cases_N)/res$pred_upper
100*(res$pred_lower - res$cases_N)/res$pred_lower
# cutoff 1: 42.35347 (-44.01956, 68.0842)
# cutoff 2: 75.25877 (36.35829, 86.41696)


# save the data
tt2 = tt %>% select(date, cases_N, pred, LPB0.025, LPB0.05, LPB0.25, UPB0.75,UPB0.95, UPB0.975)
# write.csv(tt2, row.names = F, file = 'results/case_predictions_full_model_03282021.csv')


### Now doing deaths
model_fit = lm(deaths_MA ~ cases_C_m1 + cases_C_m2 + cases_C_m3 + cases_C_m4 + 0, data = train_tt)

# 95% PI
tmp <- predict(model_fit, newdata = tt, interval = 'predict', level = .95)
tt$pred = tmp[,1]
tt$LPB0.025 = tmp[,2]
tt$UPB0.975 = tmp[,3]

# 90% PI
tmp <- predict(model_fit, newdata = tt, interval = 'predict', level = .9)
tt$LPB0.05 = tmp[,2]
tt$UPB0.95 = tmp[,3]

# 50% PI
tmp <- predict(model_fit, newdata = tt, interval = 'predict', level = .5)
tt$LPB0.25 = tmp[,2]
tt$UPB0.75 = tmp[,3]

# plot it!
p2 <- ggplot(tt, aes(x = date, y = deaths_MA)) +
  geom_line() +
  geom_line(aes(y = pred), col = 'red') + 
  geom_line(aes(y = LPB0.025), col = 'blue') +
  geom_line(aes(y = UPB0.975), col = 'blue') +
  ggtitle('death predictions') + 
  theme_minimal()
p2

#date_cutoff = tt$date[min(which(tt$cum_vacc_1/tt$N >= 0.7))] # '2021-01-14'
date_cutoff = tt$date[min(which(tt$cum_vacc_2/tt$N >= 0.7))] # '

# get deaths averted by date
res = tt %>%
  filter(date >= date_cutoff) %>%
  summarize(pred = sum(pred), pred_lower = sum(LPB0.05), pred_upper = sum(UPB0.95), deaths = sum(deaths_MA)) %>%
  mutate(deaths_averted = pred - deaths, deaths_lower = pred_lower - deaths, deaths_upper = pred_upper - deaths)

res
# cutoff 1: 382 (54.5, 709)
# cutoff 2: 445 (220, 669)

# get % reduction = 100((E-O)/E) 
100*(res$pred - res$deaths)/res$pred 
100*(res$pred_upper - res$deaths)/res$pred_upper 
100*(res$pred_lower - res$deaths)/res$pred_lower
# cutoff 1: 37.58441 (7.91803, 52.7932)
# cutoff 2: 73.76589 (58.20783, 80.88272)

# save the data
tt2 = tt %>% select(date, deaths_MA, pred, LPB0.025, LPB0.05, LPB0.25, UPB0.75,UPB0.95, UPB0.975)
# write.csv(tt2, row.names = F, file = 'results/death_predictions_full_model_03282021.csv')

### Now doing fatality 
model_fit = lm(deaths_MA ~ cases_N_m2 + cases_N_m3 , data = train_tt)
tmp <- predict(model_fit, newdata = tt, interval = 'predict', level = .95)
tt$pred = tmp[,1]
tt$LPB0.025 = tmp[,2]
tt$UPB0.975 = tmp[,3]

# 90% PI
tmp <- predict(model_fit, newdata = tt, interval = 'predict', level = .9)
tt$LPB0.05 = tmp[,2]
tt$UPB0.95 = tmp[,3]

# 50% PI
tmp <- predict(model_fit, newdata = tt, interval = 'predict', level = .5)
tt$LPB0.25 = tmp[,2]
tt$UPB0.75 = tmp[,3]

p1 <- ggplot(tt, aes(x = date, y = deaths_MA)) +
  geom_line() +
  geom_line(aes(y = pred), col = 'red') + 
  geom_line(aes(y = LPB0.25), col = 'blue') +
  geom_line(aes(y = UPB0.75), col = 'blue') +
  geom_line(aes(y = LPB0.05), col = 'purple') +
  geom_line(aes(y = UPB0.95), col = 'purple') +
  #geom_line(aes(y = LPB0.495), col = 'yellow') +
  ggtitle('case predictions') + 
  theme_minimal()
p1

tt2 = tt %>% select(date, deaths_MA, pred, LPB0.025, LPB0.05, LPB0.25, UPB0.75,UPB0.95, UPB0.975)
# write.csv(tt2, row.names = F, file = 'results/fatality_rate_predictions_full_model_03282021.csv')

#
##### Region-level case prediction #####
load('data/region_moving_average_03282021.RData')

df_wide = df_wide %>% 
  mutate(C_log = log1p(cases_C),
         C_m1_log = log1p(cases_C_m1))

train = train %>% 
  mutate(C_log = log1p(cases_C),
         C_m1_log = log1p(cases_C_m1))

plot_list = list()
iter = 0
full_res = NULL
R = 5000

# run the model and parametric bootstrap for each region
for(n in unique(df_wide$NOM)){
  iter = iter + 1
  
  # get full and training data for this region
  full_tmp = df_wide %>% filter(NOM == n) %>% ungroup()
  train_tmp = train %>% filter(NOM == n) %>% ungroup()
  
  # fit the model
  model_fit = MASS::glm.nb(cases_N ~ C_log + C_m1_log + offset(log(S)), data = train_tmp, control = glm.control(maxit = 5000))
  
  # get the model parameters
  beta_hat <- model_fit$coefficients
  beta_vcov <- vcov(model_fit)
  theta = summary(model_fit)$theta
  
  # parametric bootstrap
  sim.boot <- sapply(1:R, function(r){
    
    #indicator 
    beta_boot <- MASS::mvrnorm(1,beta_hat,beta_vcov)
    pred_boot <- (full_tmp %>% 
                    mutate(intercept=1) %>%
                    dplyr::select(intercept, C_log, C_m1_log) %>%
                    as.matrix())%*%as.matrix(beta_boot)
    
    pred_boot_exp <- exp(pred_boot + log(full_tmp$S)) 
    x = MASS::rnegbin(n = nrow(full_tmp), mu = pred_boot_exp, theta = theta)
    x # UPDATE
    
  }) 
  
  # get predictions and prediction intervals
  full_tmp$pred = predict(model_fit, full_tmp, type = 'response')
  full_tmp$LPB0.025 <- apply(sim.boot,1,quantile,.025)#,na.rm=TRUE)
  full_tmp$LPB0.05 <- apply(sim.boot,1,quantile,.05)#,na.rm=TRUE)
  full_tmp$LPB0.25 <- apply(sim.boot,1,quantile,.25)#,na.rm=TRUE)
  full_tmp$UPB0.75 <- apply(sim.boot,1,quantile,.75)#,na.rm=TRUE)
  full_tmp$UPB0.95 <- apply(sim.boot,1,quantile,.95)#,na.rm=TRUE)
  full_tmp$UPB0.975 <- apply(sim.boot,1,quantile,.975)#,na.rm=TRUE)
  full_tmp$median <- apply(sim.boot,1,quantile,.50,na.rm=TRUE)
  
  # get residuals  
  full_tmp$res = full_tmp$cases_N - full_tmp$pred
  
  # update full results
  full_res = rbind(full_res, full_tmp)
  
  # make the plot
  p1 <- ggplot(full_tmp, aes(x = date, y = cases_N)) +
    geom_line() +
    geom_line(aes(y = pred), col = 'red') + 
    geom_line(aes(y = LPB0.25), col = 'blue') +
    geom_line(aes(y = UPB0.75), col = 'blue') +
    ggtitle(n) + 
    theme_minimal()
  
  plot_list[[iter]] = p1
}

# plot them all
cowplot::plot_grid(plotlist = plot_list)

full_res = full_res %>%
  select(date, year, NOM, cases_N, pred, LPB0.025, LPB0.05, LPB0.25, UPB0.75, UPB0.95, UPB0.975)
# write.csv(full_res, row.names = F, file = 'results/region_epidemic_fits_03282021.csv')

### get the cases averted by the date cutoff - this is for creating Table 1 in the paper
date_cutoff = as.Date('2021-01-14')
a = full_res %>% group_by(NOM) %>%
  filter(date >= date_cutoff) %>%
  mutate(diff = pred - cases_N, diff_U = UPB0.95 - cases_N, diff_L = LPB0.05 - cases_N) %>% 
  summarize(cases_averted = sum(diff), lower = sum(diff_L), upper = sum(diff_U), pred = sum(pred), pred_lower = sum(LPB0.05), pred_upper = sum(UPB0.95), cases_N = sum(cases_N)) %>%
  select(NOM, cases_averted_Jan14 = cases_averted, lower_Jan14 = lower, upper_Jan14 = upper)

date_cutoff = as.Date('2021-02-06')
b = full_res %>% group_by(NOM) %>%
  filter(date >= date_cutoff) %>%
  mutate(diff = pred - cases_N, diff_U = UPB0.95 - cases_N, diff_L = LPB0.05 - cases_N) %>% 
  summarize(cases_averted = sum(diff), lower = sum(diff_L), upper = sum(diff_U), pred = sum(pred), pred_lower = sum(LPB0.05), pred_upper = sum(UPB0.95), cases_N = sum(cases_N)) %>%
  select(NOM, cases_averted_Feb6 = cases_averted, lower_Feb6 = lower, upper_Feb6 = upper)

c = merge(a,b)

# make strings for putting in the paper table
c$str_1 = sprintf('%s (%s, %s)',round(c$cases_averted_Jan14), round(c$lower_Jan14), round(c$upper_Jan14))
c$str_2 = sprintf('%s (%s, %s)',round(c$cases_averted_Feb6), round(c$lower_Feb6), round(c$upper_Feb6))

# write.csv(c, row.names = F, file = 'results/region_cases_averted_03282021.csv')

##### Region-level death prediction #####
load('data/region_moving_average_03282021.RData')

df_wide = df_wide %>% 
  mutate(C_log = log1p(cases_C),
         C_m1_log = log1p(cases_C_m1))

train = train %>% 
  mutate(C_log = log1p(cases_C),
         C_m1_log = log1p(cases_C_m1))

res = NULL
res2 = NULL
res3 = NULL
par(mfrow = c(3,3))
for(n in unique(df_wide$NOM)){
  tt = df_wide %>% filter(NOM == n)
  train_tt = train %>% filter(NOM == n)
  
  # train model
  model_fit = lm(deaths_MA ~ cases_C_m1 + cases_C_m2 + cases_C_m3 + cases_C_m4 + 0, data = train_tt)
  
  tmp <- predict(model_fit, newdata = tt, interval = 'predict', level = .95)
  tt$pred = tmp[,1]
  tt$LPB0.025 = tmp[,2]
  tt$UPB0.975 = tmp[,3]
  
  # 90% PI
  tmp <- predict(model_fit, newdata = tt, interval = 'predict', level = .9)
  tt$LPB0.05 = tmp[,2]
  tt$UPB0.95 = tmp[,3]
  
  # 50% PI
  tmp <- predict(model_fit, newdata = tt, interval = 'predict', level = .5)
  tt$LPB0.25 = tmp[,2]
  tt$UPB0.75 = tmp[,3]
  
  res = rbind(res, tt)
  
  p1 <- ggplot(tt, aes(x = date, y = deaths_MA)) +
    geom_line() +
    geom_line(aes(y = pred), col = 'red') + 
    geom_line(aes(y = LPB0.25), col = 'blue') +
    geom_line(aes(y = UPB0.75), col = 'blue') +
    geom_line(aes(y = LPB0.05), col = 'purple') +
    geom_line(aes(y = UPB0.95), col = 'purple') +
    #geom_line(aes(y = LPB0.495), col = 'yellow') +
    ggtitle('death predictions') + 
    theme_minimal()
  plot(p1)
  
  # get the deaths averted per region
  date_cutoff = as.Date('2021-01-14')

  # get deaths averted by date
  tt2 = tt %>%
    filter(date >= date_cutoff) %>%
    summarize(pred = sum(pred), 
              pred_lower = sum(LPB0.05), 
              pred_upper = sum(UPB0.95), 
              deaths = sum(deaths_MA)) %>%
    mutate(deaths_averted = pred - deaths, 
           deaths_lower = pred_lower - deaths, 
           deaths_upper = pred_upper - deaths,
           deaths_perc_avert = 100*(pred-deaths)/pred,
           deaths_perc_avert_upper = 100*(pred_upper - deaths)/pred_upper,
           deaths_perc_avert_lower = 100*(pred_lower - deaths)/pred_lower)
  
  res2 = rbind(res2, tt2)
  
  date_cutoff = as.Date('2021-02-06')
  
  # get deaths averted by second date
  tt3 = tt %>%
    filter(date >= date_cutoff) %>%
    summarize(pred = sum(pred), 
              pred_lower = sum(LPB0.05), 
              pred_upper = sum(UPB0.95), 
              deaths = sum(deaths_MA)) %>%
    mutate(deaths_averted = pred - deaths, 
           deaths_lower = pred_lower - deaths, 
           deaths_upper = pred_upper - deaths,
           deaths_perc_avert = 100*(pred-deaths)/pred,
           deaths_perc_avert_upper = 100*(pred_upper - deaths)/pred_upper,
           deaths_perc_avert_lower = 100*(pred_lower - deaths)/pred_lower)
  
  res3 = rbind(res3, tt3)
}

tt2 = res %>% select(NOM, date, year, deaths_MA, pred, LPB0.025, LPB0.05, LPB0.25, UPB0.75, UPB0.95, UPB0.975)
# write.csv(tt2, row.names = F, file = 'results/death_predictions_region_model_03282021.csv')


# for Table 1: only select numbers needed for paper and merge them
res2 = res2 %>% 
  select(NOM, deaths_averted_Jan14 = deaths_averted, lower_Jan14 = deaths_lower, upper_Jan14 = deaths_upper)

res3 = res3 %>% 
  select(NOM, deaths_averted_Feb6 = deaths_averted, lower_Feb6 = deaths_lower, upper_Feb6 = deaths_upper)

c = merge(res2, res3)

# make strings for putting in the paper table
c$str_1 = sprintf('%s (%s, %s)',round(c$deaths_averted_Jan14), round(c$lower_Jan14), round(c$upper_Jan14))
c$str_2 = sprintf('%s (%s, %s)',round(c$deaths_averted_Feb6), round(c$lower_Feb6), round(c$upper_Feb6))

# write.csv(c, row.names = F, file = 'results/region_deaths_averted_03282021.csv')
# 
##### Logistic Regression Analysis #####
load('data/comarcas_data_03282021.RData')

# Create the binary outcome
df_wide$y = as.integer(df_wide$cases_N > 0)
train$y = as.integer(train$cases_N > 0) 

# only keep locations with at least one of each class (at least one transmission week and one non-transmission week)
locations = train %>%
  group_by(NOM) %>%
  summarize(y0 = sum(y == 0), y = sum(y)) %>%
  ungroup() %>%
  filter(y >= 1, y0 >= 1) %>%
  pull(NOM)

### Do single-county model fits
# initialize results df
res = NULL 
res_noLOO = NULL # results without leave-one-out prediction (not used in the paper)

# Cycle through each county
for(nn in locations){
  # subset for that county
  tmp = df_wide %>% filter(NOM == nn)
  train_tmp = tmp %>% filter(week <= 52) #
  test_tmp = tmp %>% filter(week >= 53)
  
  # set up predictions
  tmp$y_pred = 0
  
  # do Leave-One-Out predictions
  for(ww in sort(unique(train_tmp$week))){
    train_LOO = train_tmp %>% filter(week != ww)
    test_LOO = train_tmp %>% filter(week == ww)
    
    # log offset
    glm_fit <- tryCatch({
      glm(y ~ cases_C + cases_C_m1 + offset(log(S)), family = binomial(), data = train_LOO, control = list(maxit = 50))
    }, error = function(e){
      print(sprintf('not using offset for %s', nn))
      glm(y ~ cases_C + cases_C_m1 + offset(S), family = binomial(), data = train_LOO, control = list(maxit = 50))
    })
    
    test_LOO$y_pred =  predict(glm_fit, newdata = test_LOO, type = 'response')
    
    res = rbind(res, test_LOO)
  }
  
  # store the LOO model fits on the evaluation period
  test_tmp$y_pred = predict(glm_fit, newdata = test_tmp, type = 'response')
  res = rbind(res, test_tmp)
  
  # store the model fits without LOO
  tmp$y_pred = predict(glm_fit, newdata = tmp, type = 'response')
  res_noLOO = rbind(res_noLOO, tmp)
}

# Get expected numbers and observed by week
tmp = res %>% 
  group_by(week) %>%
  summarize(cases_C = sum(cases_C), obs = sum(y), exp = sum(y_pred)) %>%
  mutate(O_over_E = obs/exp)

# par(mfrow = c(1,1))
# plot(tmp$week, tmp$O_over_E, xlab = 'week', ylab = 'observed/predicted outbreaks', main = 'observed outbreaks over predicted')
# abline(1,0)

tmp = tmp %>% select(week, obs, exp, O_over_E)
# write.csv(tmp, row.names = F, file = 'results/logistic_outbreak_analysis_03282021.csv')
