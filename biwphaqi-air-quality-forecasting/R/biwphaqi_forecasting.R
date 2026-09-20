# ==============================================================================
# Formulation, Characterization, and Time Series Forecasting of a Baseline
# Impact-Weighted Philippine Air Quality Index (BIWPhAQI)  (STT263A / STT163A)
#
# Annotated analysis code reproduced from Appendix A of the final written report.
# Pipeline: pollutant conversion -> RR-based weights -> covariance-adjusted index
#           -> population-weighted monthly series -> SARIMA / ETS / NNAR -> forecasts
#
# Data: see data/README.md. Set your working directory to the folder holding the
# CSV files before running (e.g. setwd("data"); source("../R/...R")).
# ==============================================================================

#Load Packages
library(fpp2)
library(devtools)
library(GGally)
library(patchwork)
library(seasonal)
library(moments)
library(tseries)
library(lubridate)
library(dplyr)
library(tidyr)
library(fracdiff)
library(urca)

#Data Processing
df <- read.csv('climate_air_quality.csv')

#Converting each pollutant to the same metric
conversions <- c(no2 = 1.868,
          o3 = 1.949,
          so2 = 2.602,
          co = 1145)

df <- df %>%
 mutate(
   no2 = no2 * conversions["no2"],
   o3 = o3 * conversions["o3"],
   so2 = so2 * conversions["so2"],
   co = co * conversions["co"]
 )

#Taking the B coefficients from RR as weights
weights <- c(pm25 = log(1.14)/10,
        pm10 = log(1.12)/10,
        no2 = log(1.05)/10,
        o3 = log(1.05)/10,
        so2 = log(1.0067)/10,
        co = log(1.0318)/1000)

#Computing the covariance-adjusted Index
pollutants <- df %>%
 select(pm25, pm10, no2, o3, so2, co) %>%


 drop_na()
Sigma <- cov(pollutants)
W <- diag(weights)
weighted_cov <- W %*% Sigma %*% W
X <- as.matrix(pollutants)
adj_index <- sqrt(rowSums((X %*% weighted_cov) * X))
df_clean <- df %>%
 drop_na(pm25, pm10, no2, o3, so2, co) %>%
 mutate(adj_index = adj_index)

df2 <- df_clean %>% select(adm4_pcode,date,adj_index)

#Importing population data
pop <- read.csv('worldpop_population.csv')

pop_ts <- pop %>%
 select(adm4_pcode, date, pop_count_total) %>%
 arrange(adm4_pcode, date) %>%
 mutate(date = as.Date(date))

# Forecast 2021 and 2022
pop_forecast <- pop_ts %>%
 group_by(adm4_pcode) %>%
 group_modify(~{
  df <- .x %>% filter(!is.na(pop_count_total))
  if(nrow(df) < 2) return(tibble(date = as.Date(character()), pop_count_total = numeric()))

  df <- df %>% mutate(year_num = year(date))
  fit <- lm(pop_count_total ~ year_num, data = df)
  new_dates <- tibble(year_num = 2021:2022)
  preds <- predict(fit, newdata = new_dates)

  tibble(date = as.Date(paste0(2021:2022, "-01-01")),
       pop_count_total = preds)
 }) %>%
 ungroup()

# Combine original and forecasted data
pop_combined <- bind_rows(pop_ts, pop_forecast) %>%
 arrange(adm4_pcode, date)

#Adding population data to air quality data
df2 <- df2 %>% mutate(date = as.Date(date))
pop_combined <- pop_combined %>% mutate(date = as.Date(date))
df2 <- df2 %>% mutate(year = year(date))
pop_combined <- pop_combined %>% mutate(year = year(date))


df3 <- df2 %>%
 left_join(pop_combined %>% select(adm4_pcode, year, pop_count_total),
        by = c("adm4_pcode", "year"))

#Weighted Averaging based on population
daily_index <- df3 %>%
 group_by(date) %>%
 summarise(
   adj_index_pop = sum(`adj_index` * pop_count_total, na.rm = TRUE) /
sum(pop_count_total, na.rm = TRUE),
   .groups = "drop"
 )

#Monthly Averaging
df_monthly <- daily_index %>%
 mutate(month = floor_date(date, unit = "month")) %>%
 group_by(month) %>%
 summarise(
   adj_index_pop = mean(adj_index_pop, na.rm = TRUE),
   .groups = "drop"
 )
write.csv(df_monthly, "monthly_adj_index_fixed.csv", row.names = FALSE)

#Making a time series object
str(df_monthly[,2])
aqi_ts <- ts(unlist(df_monthly[,2]), start = c(2003, 1), freq = 12)

#Characterizing the AQI
autoplot(aqi_ts) + xlab("Year") + ylab("IB-Ph AQI") + theme_minimal()

summary(aqi_ts)

find_date <- function(value) {
  idx <- which.min(abs(aqi_ts - value))
  year <- 2003 + (idx - 1) %/% 12
  month <- 1 + (idx - 1) %% 12
  return(data.frame(value = aqi_ts[idx], idx = idx, year = year, month = month))
}

find_date(7.790); find_date(9.423); find_date(10.010);find_date(10.718);find_date(15.320)

skewness(aqi_ts); kurtosis(aqi_ts); sd(aqi_ts)

tsoutliers(aqi_ts)

plot(decompose(aqi_ts))


plot(decompose(aqi_ts,type='multiplicative'))
autoplot(stl(aqi_ts, s.window='periodic')) + ggtitle("STL decomposition") + theme_minimal()

ggseasonplot(aqi_ts) +
 ggtitle("Seasonal Plot of aqi_ts") +
 ylab("Temperature Anomaly") +
 xlab("Month") +
 theme(legend.position = "none") +
 theme_minimal()

ggAcf(aqi_ts) + theme_minimal()
ggPacf(aqi_ts) + theme_minimal()

summary(ur.kpss(aqi_ts))
lam <- BoxCox.lambda(aqi_ts)
aqi_transformed <- BoxCox(aqi_ts, lambda = BoxCox.lambda(aqi_ts))
summary(ur.kpss(diff(aqi_transformed)))
autoplot(diff(aqi_transformed)) + theme_minimal()

#Model Training
train<- window(aqi_ts,end=c(2020,12))
test <- window(aqi_ts,start=c(2021,1))

#ARIMA Selection
get.best.arima <- function(transformed_train, maxord=c(2,2,2,2,2,2))
{
  best.aic<-1e8
  n<-length(transformed_train)
  for (p in 0:maxord[1]) for (d in 0:maxord[2]) for (q in 0:maxord[3])
    for (P in 0:maxord[4]) for (D in 0:maxord[5]) for (Q in 0:maxord[6])
    {
      fit<-arima (transformed_train, order=c(p,d,q), seas=list(order=c(P,D,Q),
frequency(transformed_train)), method="CSS")
      fit.aic<--2*fit$loglik+(log(n)+1)*length(fit$coef)
      if (fit.aic<best.aic)
      {
        best.aic<-fit.aic
        best.fit<-fit
        best.model<-c(p,d,q,P,D,Q)
      }
    }
  list(best.aic,best.fit,best.model)
}

# NOTE: transformed_train is defined here (the report's appendix defines it later,
# under "NNAR Selection") so that the script runs top-to-bottom without error.
transformed_train <- BoxCox(train, lambda = BoxCox.lambda(train))
get.best.arima(transformed_train)




(model_102_210 <- auto.arima(train,lambda = BoxCox.lambda(train),biasadj = T))
(model_200_101 <- Arima(train,order=c(2,0,0),seas=list(order=c(1,0,1),frequency(train)),
lambda = BoxCox.lambda(train),biasadj = T))
(model_111_100 <- Arima(train,order=c(1,1,1),seas=list(order=c(1,0,0),frequency(train)),
lambda = BoxCox.lambda(train),biasadj = T))
(model_611_100 <- Arima(train,order=c(6,1,1),seas=list(order=c(1,0,0),frequency(train)),
lambda = BoxCox.lambda(train),biasadj = T))
(model_611_101 <- Arima(train,order=c(6,1,1),seas=list(order=c(1,0,1),frequency(train)),
lambda = BoxCox.lambda(train),biasadj = T))

models <- list(
  m1 = model_102_210,
  m2 = model_200_101,
  m3 = model_111_100,
  m4 = model_611_100,
  m5 = model_611_101
)

sapply(models, AIC)



#NNAR Selection
transformed_train <- BoxCox(train, lambda = BoxCox.lambda(train))

AR1 <- Arima(
  transformed_train,
  order = c(1, 0, 0),
  seasonal = list(order = c(1, 0, 0), period = 12)
)

AR2 <- Arima(
  transformed_train,
  order = c(2, 0, 0),
  seasonal = list(order = c(1, 0, 0), period = 12)
)

AR3 <- Arima(
  transformed_train,
  order = c(3, 0, 0),
  seasonal = list(order = c(1, 0, 0), period = 12)
)

AR4 <- Arima(
 transformed_train,
 order = c(4, 0, 0),


    seasonal = list(order = c(1, 0, 0), period = 12)
)


AR5 <- Arima(
  transformed_train,
  order = c(5, 0, 0),
  seasonal = list(order = c(1, 0, 0), period = 12)
)

AR6 <- Arima(
  transformed_train,
  order = c(6, 0, 0),
  seasonal = list(order = c(1, 0, 0), period = 12)
)

AR7 <- Arima(
  transformed_train,
  order = c(7, 0, 0),
  seasonal = list(order = c(1, 0, 0), period = 12)
)

AR_models <- list(
  AR1 = AR1,
  AR2 = AR2,
  AR3 = AR3,
  AR4 = AR4,
  AR5 = AR5,
  AR6 = AR6,
  AR7 = AR7
)

sapply(AR_models, AIC)

nnetar1 <- nnetar(transformed_train, p = 5, P = 1, size = 4)
nnetar1



#Fitting the best models for ARIMA, ETS, and NNAR
(fit1 <- Arima(train,order=c(6,1,1),seas=list(order=c(1,0,1),frequency(train)), lambda =
BoxCox.lambda(train),biasadj = T))
checkresiduals(fit1)
(fit2 <- ets(train,lambda = BoxCox.lambda(train),biasadj = T))
checkresiduals(fit2)
(fit3 <- nnetar(train,lambda = BoxCox.lambda(train),biasadj = T))


checkresiduals(fit3$residuals)

#Testing 2021-2022
fc1 <- forecast(fit1,h=24)
fc2 <- forecast(fit2,h=24)
fc3 <- forecast(fit3,h=24)

accuracy(fc1,test);accuracy(fc2,test);accuracy(fc3,test)

autoplot(train) + autolayer(fc1,series='arima',PI=FALSE) +
autolayer(fc2,series='ets',PI=FALSE) + autolayer(fc3,series='nnet') +
autolayer(test,series='test') + theme_minimal()
autoplot(test,series='test')+autolayer(fc1,series='arima',PI=FALSE) +
autolayer(fc2,series='ets',PI=FALSE) + autolayer(fc3,series='nnet') + theme_minimal()

#Forecasting 2021-2026
fc1.pt <- forecast(fit1,h=72)
fc2.pt <- forecast(fit2,h=72)
fc3.pt <- forecast(fit3,h=72)

accuracy(fc1.pt,test);accuracy(fc2.pt,test);accuracy(fc3.pt,test)

autoplot(train) + autolayer(fc1.pt,series='arima',PI=FALSE) +
autolayer(fc2.pt,series='ets',PI=FALSE) + autolayer(fc3.pt,series='nnet') +
autolayer(test,series='test') + theme_minimal()
autoplot(test,series='test')+autolayer(fc1.pt,series='arima',PI=FALSE) +
autolayer(fc2.pt,series='ets',PI=FALSE) + autolayer(fc3.pt,series='nnet') + theme_minimal()
