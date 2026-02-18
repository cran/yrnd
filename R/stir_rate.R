#' stir_rate
#'
#' @param call_prices a vector of call prices in numeric format
#' @param call_strikes a vector of call strikes attached to the call prices in numeric format
#' @param put_prices a vector of put prices in numeric format
#' @param put_strikes a vector of put strikes attached to the put prices in numeric format
#' @param nb_log a number for the number of lognormal densities to model the futures contracts, either 2 or 3, in numeric format
#' @param r a number for the riskfree discount rate whose maturity is equal to the option's maturity, in numeric format
#' @param r_2 a number for the riskfree discount rate whose maturity is equal to the futures contract's maturity, in numeric format
#' @param day_count_conv a number for the day count convention, either 1 (ACT/ACT), 2 (ACT/360), 3 (ACT/365) or 4 (30/360), in numeric format
#' @param cot a number for the type of quotation of the options, either 1 (European options), 2 (American options quoted as futures) or 3 (American options), in numeric format
#' @param fut_price a number for the futures contract price on calibration date, in numeric format
#' @param fut_matu a date for the maturity date of the futures contract, in Date format
#' @param option_matu a date for the maturity date of the options, in Date format
#' @param start_date a date for the calibration date, in Date format
#' @param ref_rate a character for the name of the STIR, in character format
#' @param currency a character for the currency in which the futures contract and the options are traded, in character format
#'
#' @returns a series of values for STIR rate in numeric format, the probability density attached to each value of the STIR rate in numeric format, the cumulative density attached to each value of the STIR rate in numeric format, the type of convergence in numeric format with 0 indicating successful convergence, the mean, the standard deviation, the skewness and the kurtosis of the STIR rates' distribution at options maturity in numeric format, a plot of the RND of the STIR rates, a plot of the CDF of the STIR rates, quantiles of order 0.1%, 0.5%, 1%, 5%, 10%, 25%, 50%, 75%, 90%, 95%, 99%, 99.5% and 99.9% of the distribution of futures prices at options' maturity, in numeric format
#' @export
#' @importFrom stats approx constrOptim density dlnorm nlminb plnorm pnorm
#' @importFrom utils head tail
#' @import dplyr
#' @import lubridate
#' @import zoo
#' @import ggplot2
#'
#' @examples
#' \donttest{
#' stir_rate( c(11.44500, 1.32000, 1.19750, 1.07500, 0.95750,
#' 0.84250, 0.78750, 0.73250, 0.68000, 0.62750, 0.57750, 0.53000, 0.48500,
#' 0.44000, 0.39750, 0.35750, 0.32000, 0.28500, 0.25250, 0.22250, 0.19500,
#' 0.17000, 0.14750, 0.12750, 0.10750, 0.09250, 0.07750, 0.06500, 0.05500,
#' 0.04500, 0.03750, 0.03000, 0.02500, 0.02000, 0.01500, 0.01250, 0.01000,
#' 0.00750, 0.00500, 0.00500, 0.00250, 0.00250, 0.00250, 0.00250, 0.00024,
#' 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024,
#' 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024,
#' 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024,
#' 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024,
#' 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024,
#' 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024),
#' c(93.2500, 93.3750, 93.5000, 93.6250, 93.7500, 93.8750, 93.9375, 94.0000,
#' 94.0625, 94.1250, 94.1875, 94.2500, 94.3125, 94.3750, 94.4375, 94.5000,
#' 94.5625, 94.6250, 94.6875, 94.7500, 94.8125, 94.8750, 94.9375, 95.0000,
#' 95.0625, 95.1250, 95.1875, 95.2500,95.3125, 95.3750, 95.4375, 95.5000,
#' 95.5625, 95.6250, 95.6875, 95.7500, 95.8125, 95.8750, 95.9375, 96.0000,
#' 96.0625, 96.1250, 96.1875, 96.2500, 96.3125, 96.3750, 96.4375, 96.5000,
#' 96.5625, 96.6250, 96.6875, 96.7500, 96.8125, 96.8750, 96.9375, 97.0000,
#' 97.0625, 97.1250, 97.1875, 97.2500, 97.3125, 97.3750, 97.4375, 97.5000,
#' 97.5625, 97.6250, 97.6875, 97.7500, 97.8125, 97.8750, 97.9375, 98.0000,
#' 98.0625, 98.1250, 98.1875, 98.2500, 98.3125, 98.3750, 98.4375, 98.5000,
#' 98.5625, 98.6250, 98.6875, 98.7500, 98.8125, 98.8750, 99.0000, 99.1250,
#' 99.2500, 99.3750, 99.5000), c(0.0025, 0.0050, 0.0075, 0.0125, 0.0175,
#' 0.0300, 0.0350, 0.0425, 0.0525, 0.0625, 0.0750, 0.0900, 0.1050, 0.1225,
#' 0.1425, 0.1650, 0.1900, 0.2175, 0.2450, 0.2775, 0.3125, 0.3500, 0.3875,
#' 0.4300, 0.4725, 0.5175, 0.5675, 0.6150, 0.6675, 0.7200, 0.7750, 0.8300,
#' 0.8850, 0.9425, 1.0025, 1.0625, 1.1225, 1.1825, 1.2425, 1.3050, 1.3675,
#' 1.4300, 1.4925, 1.5550, 1.6175, 1.6800, 1.7425, 1.8050, 1.8675, 1.9300,
#' 1.9925, 2.0550, 2.1175, 2.1800, 2.2425, 2.3050, 2.3675, 2.4300, 2.4925,
#' 2.5550, 2.6175, 2.6800, 2.7425, 2.8050, 2.8675, 2.9300, 2.9925, 3.0550,
#' 3.1175, 3.1800, 3.2425, 3.3050, 3.3675, 3.4300, 3.4925, 3.5550, 3.6175,
#' 3.6800, 3.7425, 3.8050, 3.8675, 3.9300, 3.9925, 4.0550, 4.1175, 4.1800,
#' 4.3050, 4.4300, 4.5550, 4.6800, 4.8050), c(11.44500, 1.32000, 1.19750,
#' 1.07500, 0.95750,0.84250, 0.78750, 0.73250, 0.68000, 0.62750, 0.57750,
#' 0.53000, 0.48500, 0.44000, 0.39750, 0.35750, 0.32000, 0.28500, 0.25250,
#' 0.22250, 0.19500, 0.17000, 0.14750, 0.12750, 0.10750, 0.09250, 0.07750,
#' 0.06500, 0.05500, 0.04500, 0.03750, 0.03000, 0.02500, 0.02000, 0.01500,
#' 0.01250, 0.01000, 0.00750, 0.00500, 0.00500, 0.00250, 0.00250, 0.00250,
#' 0.00250, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024,
#' 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024,
#' 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024,
#' 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024,
#' 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024,
#' 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024, 0.00024),
#' 2, 0.0537, 0.0539, 1, 3, 94.7, as.Date("2024-02-29"), as.Date("2024-02-25"),
#' as.Date("2023-12-18"), "fed_fund_rate", "USD")
#' }

stir_rate <- function(call_prices, call_strikes, put_prices, put_strikes, nb_log, r, r_2, day_count_conv,
                      cot, fut_price, fut_matu, option_matu, start_date, ref_rate = NA, currency = NA){

  C <- call_prices
  P <- put_prices
  KC <- call_strikes
  KP <- put_strikes
  nb_log <- nb_log
  r <- r
  r_2 <- r_2
  day_count_conv <- day_count_conv
  cot <- cot
  fut_price <- fut_price
  fut_matu <- fut_matu
  option_matu <- option_matu
  start_date <- start_date
  ref_rate <- ref_rate
  currency <- currency

  if(length(nb_log) == 1 & length(r) == 1 & length(r_2) == 1 & length(day_count_conv) == 1 & length(cot) == 1 &
     length(fut_price) == 1 & length(fut_matu) == 1 & length(option_matu) == 1 & length(start_date) == 1 &
     length(ref_rate) == 1 & length(currency) == 1 & length(C) > 1 & length(KC) > 1 & length(P) > 1 & length(KP) > 1){

    contract_fut <- data.frame(fut_price, option_matu, start_date, fut_matu, ref_rate, currency) %>%
      rename_with(~c("fut_price", "option_matu", "start_date", "fut_matu", "name", "currency"))

    if(contract_fut$start_date < contract_fut$option_matu & contract_fut$option_matu <= contract_fut$fut_matu){

      if(day_count_conv == 1){
        contract_fut <- contract_fut %>% mutate(option_term = as.numeric(option_matu - start_date)/
                                                  as.numeric(ceiling_date(option_matu, "year") - floor_date(start_date, "year")),
                                                res_term = as.numeric(fut_matu - option_matu)/
                                                  as.numeric(ceiling_date(fut_matu, "year") - floor_date(option_matu, "year") ))
      } else if(day_count_conv == 2){
        contract_fut <- contract_fut %>% mutate(option_term = as.numeric(option_matu - start_date)/360,
                                                res_term = as.numeric(fut_matu - option_matu)/360)
      } else if(day_count_conv == 3){
        contract_fut <- contract_fut %>% mutate(option_term =as.numeric(option_matu - start_date)/365,
                                                res_term = as.numeric(fut_matu - option_matu)/365)
      } else {
        contract_fut <- contract_fut %>% mutate(stub_1 = max(0, 30 - as.numeric(format(start_date, "%d"))),
                                                stub_2 = min(30, as.numeric(format(option_matu, "%d"))),
                                                plain_months = round((as.numeric(floor_date(option_matu, "months") -
                                                                                   ceiling_date(start_date, "months")))/30),
                                                option_term = (stub_1 + stub_2 + plain_months*30)/360,
                                                stub_1_res = 30 - as.numeric(format(option_matu, "%d")),
                                                stub_2_res = min(30, as.numeric(format(fut_matu, "%d"))),
                                                plain_months_res = round(as.numeric(floor_date(fut_matu, "months") -
                                                                                      ceiling_date(option_matu, "months") )/30),
                                                res_term =  (stub_1_res + stub_2_res + max(0, plain_months_res)*30)/360)}

      rate_table <- data.frame(term = contract_fut$option_term + c(0, contract_fut$res_term), rates = c(r, r_2)) %>%
        mutate(d_fact = term*rates)
      fwd_1 <- diff(rate_table$d_fact)/diff(rate_table$term)

      call <- function(x, KC){
        d1_C <- (x[1] + x[2]^2 - log(KC))/x[2]
        d2_C <- d1_C - x[2]
        if(cot %in%c(1, 2)){call <- exp(-r*T)*(FWD*pnorm(d1_C) - KC*pnorm(d2_C))
        } else(call <- FWD*pnorm(d1_C) - KC*pnorm(d2_C))
      }

      call_mix <- function(x, KC){
        ifelse(length(x) == 7, return(x[5]*call(x[c(1, 3)], KC) + (1 - x[5])*call(x[c(2, 4)], KC) ),
               return(x[7]*call(x[c(1, 4)], KC) + x[8]*call(x[c(2, 5)], KC) + (1 - sum(x[7:8]))*call(x[c(3, 6)], KC))) }

      esp <- function(x){exp(x[1] + (x[2]^2/2))}

      esp_mix <- function(x){
        ifelse(length(x) == 7, x[5]*esp(x[c(1, 3)]) + (1 - x[5])*esp(x[c(2, 4)]),
               x[7]*esp(x[c(1, 4)]) + x[8]*esp(x[c(2, 5)]) + (1 - sum(x[7:8]))*esp(x[c(3, 6)])) }

      put <- function(x, KP){
        d1_C <- (x[1] + x[2]^2 - log(KP))/x[2]
        d2_C <- d1_C - x[2]
        if(cot %in%c(1, 2)){put <- exp(-r*T)*( -FWD*pnorm(-d1_C) + KP*pnorm(-d2_C))
        } else(put <- -FWD*pnorm(-d1_C) + KP*pnorm(-d2_C))
      }

      put_mix <- function(x, KP){
        ifelse(length(x) == 7, return(x[5]*put(x[c(1, 3)], KP) + (1 - x[5])*put(x[c(2, 4)], KP) ),
               return(x[7]*put(x[c(1, 4)], KP) + x[8]*put(x[c(2, 5)], KP) + (1 - sum(x[7:8]))*put(x[c(3, 6)], KP))) }

      if(nb_log == 2) {PR <- seq(0.1, 0.49, 0.01)
      } else {PR <- seq(0.1, 1, 0.01)
      PR <- expand.grid(c(rep(list(PR), 2)))
      PR <- PR[rowSums(PR) < 1, ]}

      suppressWarnings({

        if(cot %in%c(1, 3)){
          MSE_mix <- function(x){
            MSE_mix <- sum((C - call_mix(x, KC))^2, na.rm = T) + sum((P - put_mix(x, KP))^2, na.rm = T) + (FWD - esp_mix(x))^2
            return(MSE_mix) }
        } else { MSE_mix <- function(x){
          C_INF <- pmax(esp_mix(x) - KC, call_mix(x, KC))
          C_SUP <- exp(r*T)*call_mix(x, KC)
          P_INF <- pmax(KP - esp_mix(x), put_mix(x, KP))
          P_SUP <- exp(r*T)*put_mix(x, KP)
          itm_fwd_call <- as.numeric(KC <= esp_mix(x))
          itm_fwd_put <- as.numeric(KP >= esp_mix(x))
          w_call <- itm_fwd_call*first(tail(x, 2)) + (1 - itm_fwd_call)*last(x)
          w_put <- itm_fwd_put*first(tail(x, 2)) + (1 - itm_fwd_put)*last(x)
          CALL <- w_call*C_INF + (1 - w_call)*C_SUP
          PUT <- w_put*P_INF + (1 - w_put)*P_SUP
          MSE_mix <- sum((C - CALL)^2, na.rm = T) + sum((P - PUT)^2, na.rm = T) + (FWD - esp_mix(x))^2
          return(MSE_mix)}
        }
      })

      objective <- function(x){
        ifelse( length(PR) !=2, MSE_mix( c(x[1:4], PR[i])), MSE_mix( c(x[1:6], PR[i, 1], PR[i, 2] ))) }

      T <- contract_fut$option_term
      FWD <- contract_fut$fut_price

      m1 <- m2 <- m3 <- s1 <- s2 <- s3 <- SCE <- NA
      if(nb_log == 2){PARA <- as.matrix(data.frame(m1, m2, s1, s2, pr = PR, w1 = 0.5, w2 = 0.5, SCE))
      } else {PARA <- as.matrix(data.frame(m1, m2, m3, s1, s2, s3, pr1 = PR[, 1], pr2 = PR[, 2], w1 = 0.5, w2 = 0.5,
                                           p1_p2 = rowSums(PR), SCE))}
      start <- rep(c(log(FWD), 0.1), each = nb_log)
      if(FWD != 1){
        lower <- rep(c( (sign(1 - FWD)*0.5 + 1)*log(FWD), 1e-6), each = nb_log)
        upper <- rep(c( (sign(FWD - 1)*0.5 + 1)*log(FWD), 0.8), each = nb_log)
      } else {
        lower <- rep(c( 1, 1e-6), each = nb_log)
        upper <- rep(c( -1, 0.8), each = nb_log) }

      suppressWarnings({
        for (i in 1:length(PR)){
          sol <- nlminb(start = start, objective = objective, lower = lower, upper = upper, control = list(iter.max = 500))
          PARA[i, grep( paste( c("m", "s"), collapse = "|"), colnames(PARA))] <- sol$par
          PARA[i, "SCE"] <- sol$objective
        }
      })

      if(length(which(PARA[, ncol(PARA)]!="Inf")) != 0){
        PARA <- PARA
      } else{
        if(nb_log == 2){
          PARA <- data.frame(m1 = start[1], m2 = start[2], s1 = start[3], s2 = start[4],
                             pr = 0.25, w1 = 0.5, w2 = 0.5, SCE = NA)
        }else{PARA <- data.frame(m1 = start[1], m2 = start[2], m3 = start[3], s1 = start[4],
                                 s2 = start[5], s3 = start[6], pr1 = 0.2, pr2 = 0.2, w1 = 0.5, w2 = 0.5, SCE = NA)  }
      }

      PARA <- PARA[ !is.na(PARA[, "m1"]), ]
      if(nrow(PARA) > 1){
        param <- PARA[which.min(PARA[, "SCE"]), -ncol(PARA)]
      } else {param <- as.matrix(PARA[, -ncol(PARA)], nrow = nrow(PARA))}

      param[param == 0] <- 1e-6
      param <- matrix(param)

      L <- U <- rep(0, length(param))
      L[sign(param) == -1] <- 1.5*param[sign(param) == -1]
      L[sign(param) == 1] <- 1e-2*param[sign(param) == 1]
      U[sign(param) == -1] <- 1e-2*param[sign(param) == -1]
      U[sign(param) == 1] <- 1.5*param[sign(param) == 1]
      CI <- c(L, -U)
      UI <- rbind(diag(length(L)), -diag(length(L)))
      suppressWarnings({
        solu <- constrOptim(param, MSE_mix, NULL, ui = UI, ci = CI, mu = 1e-05, method = "Nelder-Mead")
      })
      params <- solu$par[1:(3*nb_log - 1)]

      range_px <- range(c(KP, KC))
      PX <- Reduce(seq, 1e3*range_px)*1e-3

      sub <- function(x, y){ x[3]*dlnorm(y, meanlog = x[1], sdlog = x[2]) }
      PDF <- function(x, y){
        ifelse(length(params) == 5,
               return(sub(x[c(1, 3, 5)], y) + sub(c(x[c(2, 4)], 1 - x[5]), y) ),
               return(sub(x[c(1, 4, 7)], y) + sub(x[c(2, 5, 8)], y) + sub( c(x[c(3, 6)], 1 - sum(x[7:8])), y))) }

      sub_2 <- function(x, y){ x[3]*plnorm(y, meanlog = x[1], sdlog = x[2]) }
      CDF <- function(x, y){
        ifelse(length(params) == 5,
               return(sub_2(x[c(1, 3, 5)], y) + sub_2(c(x[c(2, 4)], 1 - x[5]), y) ),
               return(sub_2(x[c(1, 4, 7)], y) + sub_2(x[c(2, 5, 8)], y) + sub_2( c(x[c(3, 6)], 1 - sum( x[7:8])), y)) ) }

      DNR <- PDF(params, PX)

      x_axis <- 1e-2
      PX_2 <- PX
      range_px_2 <- range_px
      ratio <- 1
      while (ratio > 1e-5){
        integral <- sum(rollmean(PDF(params, PX_2), 2)*diff(PX_2), na.rm = T)
        range_px_2 <- c(1 - x_axis, 1)*range_px_2
        PX_2 <- Reduce(seq, 1e3*range_px_2)*1e-3
        integral_2 <- sum(rollmean(PDF(params, PX_2), 2)*diff(PX_2), na.rm = T)
        ratio <- integral_2 - integral}

      ratio <- 1
      while (ratio > 1e-5){
        integral <- sum(rollmean(PDF(params, PX_2), 2)*diff(PX_2), na.rm = T)
        range_px_2 <- c(1, 1 + x_axis)*range_px_2
        PX_2 <- Reduce(seq, 1e3*range_px_2)*1e-3
        integral_2 <- sum(rollmean(PDF(params, PX_2), 2)*diff(PX_2), na.rm = T)
        ratio <- integral_2 - integral}

      while (sum(rollmean(PDF(params, PX_2), 2)*diff(PX_2), na.rm = T) < 0.99){
        range_px_2 <- c(1 - x_axis, 1 + x_axis)*range_px_2
        PX_2 <- Reduce(seq, 1e3*range_px_2)*1e-3}

      extension <- diff(range(PX_2))/diff(range(PX))

      if(extension <= 10){
        DNR_2 <- PDF(params, PX_2)

        NCDF <- CDF(params, PX_2)

        if(DNR_2[1] < DNR_2[2] & DNR_2[1] & DNR_2[length(DNR_2) - 1] > DNR_2[length(DNR_2)] & min(DNR_2)%in%DNR_2[c(1, length(DNR_2))]){

          PX_3 <- rev(100*exp(-fwd_1*contract_fut$res_term) - PX_2*exp(-fwd_1*contract_fut$res_term))

          sub_3 <- function(x, y){
            x[3]*(dlnorm( (100 - exp(fwd_1*contract_fut$res_term)*y),
                          meanlog = x[1], sdlog = x[2])*exp(fwd_1*contract_fut$res_term) ) }

          PDF_y <- function(x, y){
            ifelse(length(params) == 5,
                   return(sub_3(x[c(1, 3, 5)], y) + sub_3(c(x[c(2, 4)], 1 - x[5]), y) ),
                   return(sub_3(x[c(1, 4, 7)], y) + sub_3(x[c(2, 5, 8)], y) + sub_3( c(x[c(3, 6)], 1 - sum(x[7:8])), y))) }

          DNR_y <- PDF_y(params, PX_3)

          df_y <- data.frame(price = PX_3, density = DNR_y)

          NCDF_y <- cumsum(rollmean(DNR_y, 2)*diff(PX_3))
          cdf_y <- data.frame(price = PX_3[-1], cdf = NCDF_y)

          E_y <- sum(rollmean(PX_3*DNR_y, 2)*diff(PX_3))
          moments_y <- function(x){ return(sum(rollmean(DNR_y*(PX_3 - E_y)^x , 2)*diff(PX_3)))}
          SD_y <- sqrt(moments_y(2))
          SK_y <- moments_y(3)/SD_y^3
          KU_y <- moments_y(4)/SD_y^4

          thres <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 0.995, 0.999)

          if(length(which(cdf_y$cdf > last(thres))) >0){

            quantiles <- list()

            for (j in 1:length(thres)){
              quantiles[[j]] <- mean(df_y$price[c(min(which(cdf_y$cdf > thres[j] - 1e-3)), max(which(cdf_y$cdf < thres[j] + 1e-3)))])}

            qt <- data.frame(quantiles) %>% rename_with(~paste0("q", 100*thres))

            graph <- PX_3 >= qt$q0.1 & PX_3 <= qt$q99.9
            PX_graph <- PX_3[graph]
            DNR_graph <- DNR_y[graph]
            NCDF_graph <- cdf_y$cdf[graph]

            df_graph <- data.frame(price = PX_graph, density = DNR_graph)
            cdf_graph <- data.frame(price = PX_graph, cdf = NCDF_graph)

            pdf_y <- ggplot() + geom_line(data = df_graph, aes(x = price, y = density)) +
              labs(x = "future yields (%)", y = "probability density") + theme_bw() +
              theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm")) +
              labs(title = paste0(contract_fut$name, " future yields (%) on ", contract_fut$option_matu, " as of ",
                                  contract_fut$start_date),
                   subtitle = paste0("Probability Density for a mixture of ", nb_log, " lognormals"))

            ncdf_y <- ggplot() + geom_line(data = cdf_graph, aes(x = price, y = cdf)) +
              labs(x = "future yields (%)", y = "cumulative probability") + theme_bw() +
              theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm")) +
              labs(title = paste0(contract_fut$name, " future yields (%) on ", contract_fut$option_matu, " as of ",
                                  contract_fut$start_date),
                   subtitle = paste0("Cumulative Probability for a mixture of ", nb_log, " lognormals"))

            stir_rate <- list(df_y$price, df_y$density,  cdf_y$cdf, solu$convergence, E_y, SD_y, SK_y, KU_y,
                              pdf_y, ncdf_y, qt)
            names(stir_rate) <- c("rates", "rnd_r", "cdf_r", "CV", "mean", "stddev", "skew", "kurt",
                                  "rnd_r_plot", "cdf_r_plot", "quantiles")
            return(stir_rate)

          } else {message(paste0("A mixture of ", nb_log, " lognormal distributions is not convient for this data"))}
        } else {message("impossible to retrieve a density")}
      } else {message("impossible to retrieve a density")}
    } else {message("input dates are not consistent")}
  } else {message("inputs have not the required length")}
}
