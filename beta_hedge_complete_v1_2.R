library(quantmod)
library(PerformanceAnalytics)
library(tidyquant)
library(BatchGetSymbols)
library(rollRegres)

data_env <- new.env()
##function to build adjusted price xts object##

# set start date
start_date <- as.Date("2015-01-01")

# set desired period length or rehedging length
period_length <- "weekly"

# set beta window, number of periods in calculation
beta_length <- 78

# retrieve symbols of desired universe
stocks <- tq_index("DOW")
stock_list <- sort(stocks$symbol)

# tidyquamt version of retrieving stock data
#symbol_prices <- BatchGetSymbols(stock_list)

# function to retrieve stock universe prices
symbol_prices <- function(sym, enviro, source = "yahoo", priceType, from = NA){
  for(i in 1:length(sym)){
    getSymbols(sym[i], env = enviro, src = source, from = start_date)
  }
  output <- do.call(merge, eapply(enviro, priceType))
  return(output)
}
stock_prices <- symbol_prices(stock_list, enviro = data_env, source = "yahoo", priceType = "Ad")

# remove any members which have missing data for simplicity
stock_prices <- stock_prices[ , colSums(is.na(stock_prices)) == 0]

# calculate returns
returns <- do.call(merge, lapply(stock_prices, periodReturn, period = period_length))
colnames(returns) <- colnames(stock_prices)

sp500 <- getSymbols("^GSPC", auto.assign = FALSE, from = start_date)
sp500ad <- Ad(sp500)
sp500returns <- periodReturn(sp500ad, period = period_length)
colnames(sp500returns) <- "SP500"

# get 3 month treasury yield
tbill <- getSymbols("^IRX", from = start_date)
tbill <- Cl(IRX)
tbill <- na.locf(tbill)
tbill_yield <-  (1 + (tbill/100))^(1/252)
offset_index <- index(sp500returns) - start_date
tbill_period_yield <- rollapply(tbill_yield, width = offset_index, function(r) cumprod(r))
tbill_period_yield <- tbill_period_yield - 1
tbill_period_yield <- tbill_period_yield[index(sp500returns)]

# calculate standard deviation of market and universe
rolling_sd <- rollapplyr(returns, width = beta_length, FUN = StdDev)
sp500_rolling_sd <- rollapplyr(sp500returns, width = beta_length, FUN = StdDev)

# calculate correlations of market and universe
market <- merge(sp500returns, returns)
roll_cor <- rollapplyr(data = market, width = beta_length, function(x) cor(x[,1], x[,-1]), by.column = FALSE)

# calculate universe betas
betas <- (rolling_sd * roll_cor) / as.vector(sp500_rolling_sd)

# subset betas to include full window
na_betas <- betas[(beta_length + 1):nrow(betas),]

# create matrix for beta ranks
elements <- nrow(na_betas) * ncol(na_betas)
rank <- matrix(rep(1, elements), nrow = nrow(na_betas))

# rank the betas of the universe from smallest to largest
for (k in 1:nrow(na_betas)){
  for (i in 1:ncol(na_betas)){
    for (j in 1:ncol(na_betas)){
      if (na_betas[k, i] > na_betas[k, j]){
        rank[k, i] <- rank[k, i] + 1
      }
    }
  }
}

rank <- xts(rank, order.by = index(na_betas))
colnames(rank) <- colnames(returns)

# create upper beta portfolio weights
length_upper <- round(ncol(betas) / 2)
length_lower <- ncol(betas) - length_upper
weights_upper <- rep(0, length_upper)
for (i in 1:length_upper){
  weights_upper[i] <- (length_upper + 1 - i) / (length_upper * (length_upper + 1) / 2)
}
weights_upper <- sort(weights_upper)
beta_weights_upper <- c(rep(0, length_lower), weights_upper)

# create lower portfolio weight vector
weights_lower <- rep(0,length_lower)
for (i in 1:length_lower){
  weights_lower[i] <- (length_lower + 1 - i) / (length_lower * (length_lower + 1) / 2)
}
beta_weights_lower <- c(weights_lower, rep(0, length_upper))

# create total portfolio weights vector
weights<- c(beta_weights_lower, - 1 * beta_weights_upper)

# calculate individual stock weight by beta rank in lower portfolio
lower_port_betas <- na_betas * beta_weights_lower[rank]

# calculate the lower portfolio beta before leverage
lower_port_beta <- .xts(x = rowSums(lower_port_betas), .index(lower_port_betas))

# calculate the lower porfolios weight after leverage
lower_port_weight <- 1 / lower_port_beta

# calculate the weigthed return of each stock in the lower portfolio
lower_port_returns <- returns[(beta_length + 1):nrow(returns),] * beta_weights_lower[rank]

# calculate the return for each period of the stocks in the lower portfolio with leverage
weighted_lower_port_return <- lower_port_returns * as.vector(lower_port_weight)

# begin process for upper portfolio
# calculate individual stock weight by beta rank in upper portfolio
upper_port_betas <- na_betas * beta_weights_upper[rank]

# calculate the upper portfolio beta before leverage
upper_port_beta <- .xts(x = rowSums(upper_port_betas), .index(upper_port_betas))

# calculate the upper porfolios weight after leverage, and as a short position
upper_port_weight <- -1 / upper_port_beta

# calculate the weigthed return of each stock in the upper portfolio
upper_port_returns <- returns[(beta_length + 1):nrow(returns),] * beta_weights_upper[rank]

#  calculate the return for each period of the stocks in the upper portfolio with leverage, as short positions
weighted_upper_port_return <- upper_port_returns * as.vector(upper_port_weight)

# determine tbill yield weight
tbill_weight <- -1 * (lower_port_weight + upper_port_weight)
tbill_weight_return <- tbill_weight * tbill_period_yield

# calculate total return of upper and lower portfolios each period
total_lower_return <- .xts(x = rowSums(weighted_lower_port_return), .index(weighted_lower_port_return))
total_upper_return <- .xts(x = rowSums(weighted_upper_port_return), .index(weighted_upper_port_return))
colnames(total_lower_return) <- "Low Beta Portfolio"
colnames(total_upper_return) <- "High Beta Portfolio"

# determine total portfolio return
total_return <- total_lower_return + total_upper_return - tbill_weight_return
colnames(total_return) <- "Portfolio Returns"
all_returns <- merge(total_return, total_lower_return, total_upper_return)

# add sp500 returns for comparison
sp500returns_comp <- sp500returns[index(all_returns)]
all_returns <- merge(all_returns, sp500returns_comp)

# determine total cumulative portfolio return
total_cum_returns <- cumprod(1+total_return)
total_cum_lower_returns <- cumprod(1+total_lower_return)
total_cum_upper_returns <- cumprod(1+total_upper_return)

# plot cummulative returns of total portfllio return, high beta port return and low beta port return
plot.zoo(merge(total_cum_returns, total_cum_lower_returns, total_cum_upper_returns), 
         main = "Portfolio Performance", col = (1:3), 
         ylab = c("Portfolio Returns", "Low Beta Port", "High Beta Port", plot.type = "multiple"))

chart.CumReturns(all_returns, main = "Portfolio Performance",legend.loc = "topleft", wealth.index = TRUE)
charts.PerformanceSummary(all_returns, Rf = tbill_period_yield, wealth.index = TRUE, 
                          main = "Portfolio Performance")

# Portfolio Performance Statistics
table.AnnualizedReturns(all_returns, Rf = tbill_period_yield)
VaR(all_returns)
table.Drawdowns(total_return)