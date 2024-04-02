#################################################
# Dataset dissertação - v2
# Feito por: João Vitor Dias
# Orientador: Carlos eugênio
#################################################

library(tidyverse)
library(readxl)
library(lubridate)
library(fredr)
library(openxlsx)
library(lmtest)
library(sandwich)
library(ustyc)
library(xts)
library(zoo)

rm(list = ls())

directory <- "G:/Meu Drive/JOAO/Mestrado/Segundo ano/Dissertacao/Data"
setwd(directory)

#set.seed(123) 

##############################################################################################

# seguindo a dissertação do Jivago Vasconcelos

# Annual US population and consumption expenditures data are from Robert
# Shillerís database between 1947 and 2004. We use the return of 90-day
# Treasury bill deáated by the CPI index as our measure of US real interest
# rate. Stock market excess returns, Treasury bill rates and CPI are from
# CRSP.

#Stock market excess returns 
#Treasury bill rates  
#CPI


# Shiller's data hhttp://www.econ.yale.edu/~shiller/data/chapt26.xlsx

# Consumption: annual average personal consumption expenditures on nondurable goods and
# services series from the National Income and Product Accounts of the United States.

# https://apps.bea.gov/iTable/?reqid=19&step=2&isuri=1&categories=survey#eyJhcHBpZCI6MTksInN0ZXBzIjpbMSwyLDMsM10sImRhdGEiOltbImNhdGVnb3JpZXMiLCJTdXJ2ZXkiXSxbIk5JUEFfVGFibGVfTGlzdCIsIjY2Il0sWyJGaXJzdF9ZZWFyIiwiMjAyMSJdLFsiTGFzdF9ZZWFyIiwiMjAyMyJdLFsiU2NhbGUiLCItOSJdLFsiU2VyaWVzIiwiQSJdLFsiU2VsZWN0X2FsbF95ZWFycyIsIjEiXV19
# Treasuries: https://www.federalreserve.gov/DataDownload/Choose.aspx?rel=H15

# API key: 1661c2b2a9df9279db5603f57bb1d3f4 

# search the ids

#popular_funds_series <- fredr_series_search_text(
#  search_text = "Real personal consumption expenditures per capita",
#  order_by = "popularity",
#  filter_variable =  'units',
#  filter_value = 'Billions of Chained 2012 Dollars',
#  limit = 100)

# NIPA table 2.3.6
# available in https://fred.stlouisfed.org/release/tables?rid=53&eid=43891

# NIPA table 2.4.3
# available in https://fred.stlouisfed.org/release/tables?rid=53&eid=43953

services <- fredr(
  series_id = "A797RX0Q048SBEA",
  observation_start = as.Date("1947-01-01")) %>%
  dplyr::select(c("date", "value")) %>% 
  rename(Date = date) %>% 
  rename(Services = value) 

nondurables <- fredr(
  series_id = "A796RX0Q048SBEA",
  observation_start = as.Date("1947-01-01")) %>%
  dplyr::select(c("date", "value")) %>% 
  rename(Date = date) %>% 
  rename(Nondurables = value)


SP500 <- read_excel("SeP500.xlsx", sheet = 1, skip = 5, col_names = TRUE) %>% 
          dplyr::select(Date, PX_LAST) %>% 
          rename(SP500 = PX_LAST)  %>% 
          arrange(desc(Date))

SP500$Date <- as.Date(ceiling_date(SP500$Date, "quarter"))
# - months(3)

CPI <- read_excel("cpi_ocde.xlsx", sheet = 1, skip = 5, col_names = TRUE) %>% 
        dplyr::select(Date, PX_LAST) %>% 
        rename(CPI = PX_LAST) %>% 
        arrange(desc(Date))

CPI$Date <- as.Date(ceiling_date(CPI$Date, "quarter"))

Tbill90 <- read_excel("tbill_3months.xlsx", sheet = 1, skip = 5, col_names = TRUE) %>% 
            dplyr::select(Date, PX_LAST)%>% 
            rename(Tbill90 = PX_LAST)

Tbill90$Date <- as.Date(ceiling_date(Tbill90$Date, "quarter"))


data <- inner_join(services, nondurables, by = "Date") %>% 
        mutate(Consumption = services$Services + nondurables$Nondurables) %>%
        inner_join(CPI, by = "Date") %>% 
        mutate(Inflation = log(CPI) - lag(log(CPI))) %>% 
        mutate(growth_Consumption = log(Consumption) - lag(log(Consumption))) %>% 
        inner_join(SP500, by = "Date") %>% 
        mutate(e_return = log(SP500) - lag(log(SP500))) %>% 
        inner_join(Tbill90, by = "Date") %>% 
        mutate(real_yield = (1+Tbill90/100)/(1+Inflation*4)-1)


dividend <- read_excel("grid1_w5mougtl.xlsx", sheet = 1, skip = 6, col_names = TRUE) %>% 
  dplyr::select(-PX_LAST)

dividend$Date <- as.Date(ceiling_date(dividend$Date, "quarter"))

names(dividend) <- c("Date", "Dividends")

ref_inf = as.numeric(data[data$Date == max(data$Date), "CPI"])

data <- data %>% inner_join(dividend, by = "Date") 
data <- data %>% mutate(real_dividend = (Dividends/CPI)*ref_inf ) %>%
  mutate(real_SP = (SP500/CPI)*ref_inf ) %>% 
  mutate(PD = real_SP/real_dividend) %>% 
  mutate(total_return = log(real_dividend + SP500)-lag(log(SP500)))


data_wachter <- data %>% filter(Date < as.Date('2004-07-01'))

# tsc = 4  # Intervalo de tempo para simulação (no caso, mensal)
# g = 0.0228 / tsc  # Média da evolução do consumo
# sig = 0.009 / sqrt(tsc)  # Desvio padrão da evolução do consumo
# rf0 = 0.0127 / tsc  # Log da taxa livre de risco
# B = 0.009  # Coeficiente angular de rf
# gamma = 2.5  # Curvatura da função utilidade
# phi = 0.902 ^ (1 / tsc)  # Taxa de persistência


#### Average consumption growth ####

avg_consumption <- mean(data$growth_Consumption, na.rm = TRUE)

sd_consumption <- sd(data$growth_Consumption, na.rm = TRUE)

avg_consumption_a <- avg_consumption*4

sd_consumption_a <- sd_consumption*sqrt(4)

mean(data_wachter$growth_Consumption, na.rm = TRUE)*4
sd(data_wachter$growth_Consumption, na.rm = TRUE)*sqrt(4)

#### Average inflation ####

avg_inflation <- mean(data$Inflation, na.rm = TRUE)*4

sd_inflation <- sd(data$Inflation, na.rm = TRUE)*sqrt(2)

mean(data_wachter$Inflation, na.rm = TRUE)*4
sd(data_wachter$Inflation, na.rm = TRUE)*sqrt(2)

#### Average risk free rate ####

####### checar #######

avg_risk_free_a <- mean(data$real_yield, na.rm = TRUE)

sd_risk_free_a <- sd(data$real_yield, na.rm = TRUE)

avg_risk_free_N_a <- mean(data$Tbill90, na.rm = TRUE)/100

sd_risk_free_N_a <- sd(data$real_yield, na.rm = TRUE)

mean((data_wachter$real_yield), na.rm = TRUE)
mean((data_wachter$Tbill90), na.rm = TRUE)-3.68
sd((data_wachter$real_yield), na.rm = TRUE)

#### Average equity return ####

#### multiplico por 4????? #### acho q sim

avg_equity_return <- mean(data$e_return, na.rm = TRUE)*400

sd_equity_return <- sd(data$e_return, na.rm = TRUE)*sqrt(4)*100

risk_p <- mean(log(data$Dividends + data$SP500) - lag(log(data$SP500)) - data$Tbill90/400, na.rm = TRUE)*400
risk_p_sd <- sd(log(data$Dividends + data$SP500) - lag(log(data$SP500)) - data$Tbill90/400, na.rm = TRUE)*200

sharpe <- risk_p/risk_p_sd 

annual_data <- data %>%
  select(Date, real_SP, real_dividend) %>% 
  mutate(Year = as.integer(format(Date, "%Y"))) %>%
  group_by(Year) %>%
  summarize(Annual_Dividend = sum(real_dividend), 
            SP500_End_Year = last(real_SP))

PD <- mean(log(annual_data$SP500_End_Year/annual_data$Annual_Dividend))
PD_sd <- sd(log(annual_data$SP500_End_Year/annual_data$Annual_Dividend))

(mean(data_wachter$e_return, na.rm = TRUE)*4 - 0.0147)*100

sd(data_wachter$e_return, na.rm = TRUE)

#### Autocorrelation of the dividend yield ####

autocorrel <- acf(data$PD, lag.max = 1, plot = FALSE)

autocorrel_values <- (autocorrel$acf)
print(autocorrel_values)

phi <- autocorrel_values[2]

phi^(4)

### elevar a 4 ou nao? nao sei

#### Wachter's regression for B ####

#phi = 0.89^(1/4)

LW = 40

W_phi <- vector(length = LW)
for (i in 1:LW){
  W_phi[i] <- phi^i 
}

W_consumption <- vector(length = length(data$Consumption)-LW)
W_temp <- matrix(0, nrow = LW, ncol = 2)

reg_c <- data %>% select(Date, growth_Consumption) %>% arrange(desc(Date))

for (i in 1:length(W_consumption)){
  for (j in 1:LW){
    W_temp[j,] <- as.matrix(reg_c[i+j,])
  }
  W_consumption[i] <- (t(as.numeric(W_temp[,2])) %*% W_phi)[1,1] 
}


W_data <- tibble(Date = reg_c$Date[1:length(W_consumption)],
                 rates = rev(data$real_yield[1:length(W_consumption)]),
                 W_consumption = W_consumption)

W_regression <- lm(rates ~ W_consumption, data = W_data)

summary(W_regression)

coeftest(W_regression, vcov = vcovHC(W_regression, type = "HC3"))

########################################################################

phi_w = 0.89^(1/4)

W_phi2 <- vector(length = LW)
for (i in 1:LW){
  W_phi2[i] <- phi_w^i 
}

W_consumption2 <- vector(length = length(data_wachter$Consumption)-LW)
W_temp2 <- matrix(0, nrow = LW, ncol = 2)

reg_c2 <- data_wachter %>% select(Date, growth_Consumption) %>% arrange(desc(Date))

for (i in 1:length(W_consumption2)){
  for (j in 1:LW){
    W_temp2[j,] <- as.matrix(reg_c2[i+j,])
  }
  W_consumption2[i] <- (t(as.numeric(W_temp2[,2])) %*% W_phi2)[1,1] 
}


W_data2 <- tibble(Date = reg_c2$Date[1:length(W_consumption2)],
                 rates = rev(data_wachter$real_yield[1:length(W_consumption2)]),
                 W_consumption = W_consumption2)

W_regression2 <- lm(rates ~ W_consumption, data = W_data2)

coeftest(W_regression2, vcov = vcovHC(W_regression2, type = "HC3"))

########################

#h <- data %>% dplyr::select(Date, CPI, Consumption)
h <- data %>% dplyr::select(Date, Inflation, growth_Consumption)

write.xlsx(h, file = "mle_data.xlsx")
write.csv(h, file = "mle_data.csv", row.names = FALSE)

#### Autocorrelation of the dividend yield ####

yield_original <- read_excel("FRB_H15.xlsx", sheet = 1, col_names = TRUE)

yield_curve <- yield_original %>% 
               select(c("Date", "1_year", "3_years","5_years","10_years","20_years"))



yield_curve$Date <- as.POSIXct(yield_curve$Date)

# Create an xts object from the data frame, excluding the 'Date' column from the data part
yield_curve_xts <- xts(yield_curve[, -which(names(yield_curve) == "Date")], order.by=yield_curve$Date)

# Ensure all days are covered, filling missing days with NA, and then forward fill NA values
all_dates <- seq(start(yield_curve_xts), end(yield_curve_xts), by="day")
yield_curve_xts_all <- merge(yield_curve_xts, xts(, all_dates))
yield_curve_xts_filled <- na.locf(yield_curve_xts_all)

# Convert the daily series into quarterly for each interest rate column
quarterly_rates <- apply.quarterly(yield_curve_xts_filled, last)

index(quarterly_rates) <- as.POSIXct(index(quarterly_rates))

# Adjust the dates to the first day of the next quarter
adjusted_dates <- ceiling_date(index(quarterly_rates), "quarter")

# Apply the adjusted dates back to the xts object index
index(quarterly_rates) <- adjusted_dates

quarterly_rates[quarterly_rates == "ND"] <- NA


quarterly_rates_numeric <- apply(quarterly_rates, 2, function(x) as.numeric(as.character(x)))

# Convert back to xts, using the original index (dates)
quarterly_rates_numeric_xts <- xts(quarterly_rates_numeric, order.by=index(quarterly_rates))

quarterly_real_rates  <- data.frame(Date = as.Date(index(quarterly_rates_numeric_xts)), coredata(quarterly_rates_numeric_xts)) %>% 
                        inner_join(data, by = "Date") %>% 
                        mutate(r_1y = X1_year/100 - Inflation*4) %>% 
                        mutate(r_3y = X3_years/100 - Inflation*4) %>% 
                        mutate(r_5y = X5_years/100 - Inflation*4) %>% 
                        mutate(r_10y = X10_years/100 - Inflation*4) %>% 
                        mutate(r_20y = X20_years/100 - Inflation*4) %>% 
                        select("r_1y","r_3y","r_5y","r_10y","r_20y")

# Calculate the standard deviation for each column, ignoring NA values
std_devs_all_columns <- apply(quarterly_real_rates, 2, sd, na.rm = TRUE)
mean_devs_all_columns <- apply(quarterly_real_rates, 2, mean, na.rm = TRUE)

########################################################################################################################################################33


data_q <- read_excel('REER_database_ver11Feb2024.xls', sheet = 'REER_MONTHLY_51') %>% select(Date, REER_51_US)

data_q$REER_51_US <- as.numeric(as.character(data_q$REER_51_US))

data_q$Date <- as.yearqtr(as.character(data_q$Date), "%YM%m")

data_q_quarterly <- aggregate(zoo(data_q$REER_51_US), data_q$Date, tail, 1)

# Convert 'yearqtr' to the first day of the quarter in 'Date' format
data_q_quarterly <- aggregate(zoo(data_q$REER_51_US), data_q$Date, tail, 1)
index(data_q_quarterly) <- as.Date(index(data_q_quarterly))

data_q_quarterly_df <- data.frame(Date = index(data_q_quarterly), REER_51_US = coredata(data_q_quarterly)) %>% 
  mutate(delta_q = log(REER_51_US) - lag(log(REER_51_US)))

sd(data_q_quarterly_df$delta_q, na.rm = TRUE)
