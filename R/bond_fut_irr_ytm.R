#' bond_fut_irr_ytm
#'
#' @param bond_call_prices a vector of call prices on a bond futures, in numeric format
#' @param bond_call_strikes a vector of call strikes attached to the call prices, in numeric format
#' @param bond_put_prices a vector of put prices on a bond futures, in numeric format
#' @param bond_put_strikes a vector of put strikes attached to the put prices, in numeric format
#' @param stir_call_prices a vector of call prices on a STIR futures, in numeric format
#' @param stir_call_strikes a vector of call strikes attached to the call prices, in numeric format
#' @param stir_put_prices a vector of put prices on a STIR futures, in numeric format
#' @param stir_put_strikes a vector of put strikes attached to the put prices, in numeric format
#' @param r a number for the riskfree spot rate whose maturity is equal to the options' maturity, in numeric format
#' @param r_2 a number for the spot repo funding rate of the underlying bond, with maturity equal to the futures' maturity, in numeric format
#' @param day_count_conv a number for the day count convention, 1 for ACT/ACT, 2 for ACT/360, 3 for ACT/365 and 4 for 30/360, in numeric format
#' @param cot_bond a number for the bond options' style, 1 for European options, 2 for American options and 3 for American options with futures-style margin, in numeric format
#' @param bond_fut_price a number for the bond futures' price on calibration date, in numeric format
#' @param bond_cp a number for the coupon rate of the current Cheapest-to-Deliver Bond, in numeric format
#' @param bond_cp_f a number for the frequency of coupon payment of the current Cheapest-to-Deliver Bond, either 1 if the frequency is annual or 2 if semi-annual
#' @param bond_conv_factor a number for the conversion factor assigned by the futures exchange to the current Cheapest-to-Deliver Bond of the bond futures, in numeric format
#' @param sett a number for the number of days between the ex-coupon date and the coupon payment date of the current Cheapest-to-Deliver Bond, in numeric format
#' @param bond_N a number for the value of the principal of the current Cheapest-to-Deliver Bond, in numeric format (100 by default)
#' @param bond_matu a date for the maturity date of the current Cheapest-to-Deliver Bond, in Date format
#' @param bond_fut_matu a date for the maturity date of the bond futures contract, in Date format
#' @param bond_option_matu a date for the maturity date of the bond options, in Date format
#' @param cot_stir a number for the STIR options' style, 1 for European options, 2 for American options and 3 for American options with futures-style margin, in numeric format
#' @param stir_fut_price a number for the STIR futures' price on calibration date, in numeric format
#' @param stir_fut_matu a date for the maturity date of the STIR futures contract, in Date format
#' @param stir_option_matu a date for the maturity date of the STIR options, in Date format
#' @param start_date a date for the observation date, in Date format
#' @param ref_rate a character for the name of the STIR, in character format (NA by default)
#' @param country a character for the country of the issuer of the bond underlying the futures contract, in character format (NA by default)
#' @param currency a character for the currency in which the futures contract and the options are traded, in character format (NA by default)
#'
#' @returns provided they can be extracted, the discretized domains and RNDs of the bond forward repo rate and of the bond forward yield to maturity, the mean, standard deviation, skewness, kurtosis and mode of the bond forward repo rate and of the bond forward yield to maturity, in numeric format, the plots of the RND and the CDF of the bond forward repo rate and of the bond forward yield to maturity, the parameters of the RND of the bond forward price and of the bond forward repo price, in numeric format, the convergence with 0 indicating successful convergence, in numeric format
#' @export
#' @importFrom stats approx constrOptim density dlnorm nlminb plnorm pnorm
#' @importFrom utils head tail
#' @import dplyr
#' @import lubridate
#' @import zoo
#' @import ggplot2
#' @import tvm
#' @import tibble
#' @import DEoptim
#'
#' @examples
#' \donttest{
#' bond_fut_irr_ytm(c(12.64, 12.14, 11.65, 11.15, 10.65, 10.15, 9.65,
#' 9.16, 8.66, 8.16, 7.67,  7.17, 6.68, 6.19, 5.70, 5.21, 4.73, 4.25,
#' 3.78, 3.32, 2.87, 2.45, 2.05, 1.67, 1.33, 1.03, 0.77, 0.57, 0.41,
#' 0.29, 0.21, 0.15, 0.10, 0.07, 0.06, 0.04, 0.04, 0.03, 0.02, 0.02,
#' 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01),
#' c(seq(114, 136.5, 0.5), 137.5),
#' c(0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.03, 0.03, 0.03, 0.04,
#' 0.04, 0.05, 0.06, 0.07, 0.08, 0.10, 0.12, 0.15, 0.19, 0.24, 0.32,
#' 0.42, 0.54, 0.70, 0.90, 1.14, 1.44, 1.78, 2.16, 2.58, 3.02, 3.47,
#' 3.94, 4.43, 4.91, 5.41, 5.90, 6.39, 6.89, 7.39, 7.89, 8.38, 8.88,
#' 9.38, 9.88, 10.88),
#' c(seq(114, 136.5, 0.5), 137.5),
#' c(1.704999924, 1.579999924, 1.454999924, 1.329999924, 1.204999924,
#' 1.079999924, 0.954999983, 0.892499983, 0.829999983, 0.767499983,
#' 0.704999983, 0.642499983, 0.582499981, 0.519999981, 0.457499981,
#' 0.397499979, 0.334999979, 0.275000006, 0.217500001, 0.162499994,
#' 0.112499997, 0.072499998, 0.039999999, 0.02, 0.0075, 0.0025, 0.0025,
#' 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#' c(seq(95.75, 96.5, 0.125), seq(96.5625, 98.75, 0.0625),
#' seq(98.875, 99.5, 0.125)),
#' c(0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0.0025,  0.0025,
#' 0.0025,  0.005, 0.005,  0.0075,  0.012499999,  0.02, 0.029999999,
#' 0.055, 0.085000001, 0.127499998, 0.177499995, 0.234999999,
#' 0.297499985, 0.357499987, 0.419999987, 0.482499987, 0.544999957,
#' 0.607499957, 0.669999957, 0.732499957, 0.794999957, 0.857499957,
#' 0.919999957, 0.982499957, 1.044999957, 1.107499957, 1.169999957,
#' 1.232499957, 1.294999957, 1.419999957, 1.544999957, 1.669999957,
#' 1.794999957, 1.919999957, 2.044999838),
#' c(seq(95.75, 96.5, 0.125), seq(96.5625, 98.75, 0.0625),
#' seq(98.875, 99.5, 0.125)),
#' 0.0224,
#' 0.02232,
#' 1,
#' 3,
#' 126.66,
#' 0.026,
#' 1,
#' 0.770088,
#' 2,
#' bond_N = 100,
#' as.Date("2035-08-15"),
#' as.Date("2026-09-08"),
#' as.Date("2026-08-21"),
#' 3,
#' 97.45,
#' as.Date("2026-09-14"),
#' as.Date("2026-08-14"),
#' as.Date("2026-06-17"),
#' ref_rate = "3-mth Euribor",
#' country = "Germany",
#' currency = "EUR")
#' }

bond_fut_irr_ytm <- function(bond_call_prices, bond_call_strikes, bond_put_prices, bond_put_strikes,
                             stir_call_prices, stir_call_strikes, stir_put_prices, stir_put_strikes,
                             r, r_2, day_count_conv, cot_bond, bond_fut_price, bond_cp, bond_cp_f,
                             bond_conv_factor, sett, bond_N = 100, bond_matu, bond_fut_matu, bond_option_matu,
                             cot_stir, stir_fut_price, stir_fut_matu, stir_option_matu,
                             start_date, ref_rate = NA, country = NA, currency = NA){

  if(length(r) == 1 & length(r_2) == 1 & length(day_count_conv) == 1 & length(cot_bond) == 1 &
     length(bond_matu) == 1 & length(bond_fut_price) == 1 & length(bond_fut_matu) == 1 &
     length(bond_option_matu) == 1 & length(start_date) == 1 & length(bond_cp) == 1 &
     length(bond_cp_f) == 1 & length(bond_conv_factor) == 1 &  length(sett) == 1 & length(country) == 1 &
     length(currency) == 1 & length(cot_stir) == 1 & length(stir_fut_price) == 1 &
     length(stir_fut_matu) == 1 & length(stir_option_matu) & length(ref_rate) == 1 &
     length(stir_call_prices) > 1 & length(stir_call_strikes) > 1 & length(stir_put_prices) > 1&
     length(stir_put_strikes) > 1 & length(bond_call_prices) > 1 & length(bond_call_strikes) > 1 &
     length(bond_put_prices) > 1 & length(bond_put_strikes) > 1){

    bond_charac_2 <- data.frame(bond_conv_factor, bond_cp, bond_matu, bond_fut_price, bond_option_matu, start_date,
                                bond_cp_f, sett, bond_fut_matu, country, currency, bond_N) %>%
      rename_with(~c("conv_factor", "bond_cp", "bond_matu", "fut_price", "option_matu", "start_date",
                     "cp_f", "sett", "fut_matu", "country", "currency", "Nomi")) %>%
      mutate_at(c("bond_matu", "option_matu", "fut_matu", "start_date"), as.Date)

    stir_fut_matu <- as.Date(stir_fut_matu)
    stir_option_matu <- as.Date(stir_option_matu)

    if( abs(as.numeric(stir_option_matu) - as.numeric(bond_option_matu)) < 30){

      if(bond_charac_2$start_date < bond_charac_2$option_matu & bond_charac_2$option_matu <= bond_charac_2$fut_matu){

        nb_log <- 2

        bond_fut_rnd <- bond_future_price(bond_call_prices, bond_call_strikes, bond_put_prices, bond_put_strikes,
                                          nb_log, r, day_count_conv, cot_bond, bond_matu, bond_fut_price,
                                          bond_fut_matu, bond_option_matu, start_date, country, currency)

        stir_fut_p <- stir_future_price(stir_call_prices, stir_call_strikes, stir_put_prices,
                                        stir_put_strikes, nb_log, r, day_count_conv, cot_stir,
                                        stir_fut_price, stir_fut_matu, stir_option_matu, start_date,
                                        ref_rate, currency)

        marginal_bond <- marginal_repo <- ""

        if(length(bond_fut_rnd) > 0 & length(stir_fut_p) > 0){

          bond_charac_2 <- bond_charac_2 %>%
            mutate(prev_cp_dt = as.Date(paste0(format(option_matu, "%Y"), "-", format(bond_matu, "%m-%d"))))

          if(bond_charac_2$cp_f == 1){bond_fut <- bond_charac_2 %>%
            mutate_at("prev_cp_dt", ~as.Date(ifelse(option_matu < ., . %m-% years(1), .))) %>%
            mutate(curr_cp_dt = prev_cp_dt %m+% years(1), next_cp_dt = curr_cp_dt %m+% years(1))
          } else { bond_fut <- bond_charac_2 %>%
            mutate_at("prev_cp_dt", ~as.Date(ifelse(option_matu - . < - months(6), . %m-% years(1),
                                                    ifelse(option_matu - . < 0, . %m-% months(6), .)))) %>%
            mutate(curr_cp_dt = prev_cp_dt %m+% months(6), next_cp_dt = curr_cp_dt %m+% months(6))}

          if(day_count_conv == 1){
            bond_fut <- bond_fut %>% mutate(option_term = as.numeric(option_matu - start_date)/
                                              as.numeric(ceiling_date(option_matu, "year") - floor_date(start_date, "year")),
                                            res_term = as.numeric(fut_matu - option_matu)/
                                              as.numeric(ceiling_date(fut_matu, "year") - floor_date(option_matu, "year") ))
          } else if(day_count_conv == 2){
            bond_fut <- bond_fut %>% mutate(option_term = as.numeric(option_matu - start_date)/360,
                                            res_term = as.numeric(fut_matu - option_matu)/360)
          } else if(day_count_conv == 3){
            bond_fut <- bond_fut %>% mutate(option_term =as.numeric(option_matu - start_date)/365,
                                            res_term = as.numeric(fut_matu - option_matu)/365)
          } else {
            bond_fut <- bond_fut %>% mutate(stub_1 = max(0, 30 - as.numeric(format(start_date, "%d"))),
                                            stub_2 = min(30, as.numeric(format(option_matu, "%d"))),
                                            plain_months = round((as.numeric(floor_date(option_matu, "months") -
                                                                               ceiling_date(start_date, "months")))/30),
                                            option_term = (stub_1 + stub_2 + plain_months*30)/360,
                                            stub_1_res = max(0, 30 - as.numeric(format(option_matu + sett, "%d"))),
                                            stub_2_res = min(30, as.numeric(format(fut_matu, "%d"))),
                                            plain_months_res = round(as.numeric(floor_date(fut_matu, "months") -
                                                                                  ceiling_date(option_matu + sett, "months") )/30),
                                            res_term = (stub_1_res + stub_2_res + max(0, plain_months_res)*30)/360)}

          bond_fut <- bond_fut %>% mutate(res_term_2 = 0)

          rate_table <- data.frame(term = bond_fut$option_term + c(0, bond_fut$res_term), rates = c(r, r_2)) %>%
            mutate(d_fact = term*rates)
          fwd_1 <- diff(rate_table$d_fact)/diff(rate_table$term)

          if(bond_fut$fut_matu < bond_fut$curr_cp_dt){
            if(day_count_conv == 1){
              bond_fut <- bond_fut %>% mutate(acc_matu = bond_fut$Nomi*bond_cp*as.numeric(fut_matu - prev_cp_dt - sett)/
                                                as.numeric(ceiling_date(fut_matu, "year") - floor_date(fut_matu, "year") ))
            } else if(day_count_conv == 2) {
              bond_fut <- bond_fut %>% mutate(acc_matu = bond_fut$Nomi*bond_cp*as.numeric(fut_matu - prev_cp_dt - sett)/360)
            } else if(day_count_conv == 3){
              bond_fut <- bond_fut %>% mutate(acc_matu = bond_fut$Nomi*bond_cp*as.numeric(fut_matu - prev_cp_dt - sett)/365)
            } else{
              bond_fut <- bond_fut %>% mutate(stub_1 = max(0, 30 - as.numeric(format(prev_cp_dt + sett, "%d"))),
                                              stub_2 = min(30, as.numeric(format(fut_matu, "%d"))),
                                              plain_months = round(as.numeric(floor_date(fut_matu, "months") -
                                                                                ceiling_date(prev_cp_dt + sett, "months") )/30),
                                              acc_matu = bond_fut$Nomi*bond_cp/360*(stub_2 + stub_1 + max(0, plain_months)*30)) }
          } else{
            if(day_count_conv == 1){
              bond_fut <- bond_fut %>% mutate(res_term_2 = as.numeric(fut_matu - curr_cp_dt - sett)/
                                                as.numeric(ceiling_date(fut_matu, "year") - floor_date(fut_matu, "year") ),
                                              acc_matu = bond_fut$Nomi*bond_cp*res_term_2 )
            } else if(day_count_conv == 2) {
              bond_fut <- bond_fut %>% mutate(res_term_2 = as.numeric(fut_matu - curr_cp_dt - sett)/360,
                                              acc_matu = bond_fut$Nomi*bond_cp*res_term_2)
            } else if(day_count_conv == 3){
              bond_fut <- bond_fut %>% mutate(res_term_2 = as.numeric(fut_matu - curr_cp_dt - sett)/365,
                                              acc_matu = bond_fut$Nomi*bond_cp*res_term_2)
            } else{
              bond_fut <- bond_fut %>% mutate(stub_1 = max(0, 30 - as.numeric(format(curr_cp_dt + sett, "%d"))),
                                              stub_2 = min(30, as.numeric(format(fut_matu, "%d"))),
                                              plain_months = round(as.numeric(floor_date(fut_matu, "months") -
                                                                                ceiling_date(curr_cp_dt + sett, "months") )/30),
                                              res_term_2 = (stub_1 + stub_2 + max(0, plain_months)*30)/360,
                                              acc_matu = bond_fut$Nomi*bond_cp*res_term_2)}
            rate_table <- rate_table %>%
              add_row(term = bond_fut$res_term + bond_fut$option_term - bond_fut$res_term_2)
            rate_table$rates[3] <- approx(rate_table$term[1:2], rate_table$rates[1:2], xout = rate_table$term[3],
                                          method = "linear", n = 50, rule = 2, f = 0, ties = "ordered", na.rm = F)$y
            rate_table$d_fact[3] <- rate_table$term[3]*rate_table$rates[3]
            fwd_2 <- diff(rate_table$d_fact[-1])/diff(rate_table$term[-1])
          }

          if(bond_fut$cp_f == 1){ true_cp_dt <- seq(from = bond_fut$curr_cp_dt, to = bond_fut$bond_matu, by = "year")
          } else { true_cp_dt <- seq(from = bond_fut$curr_cp_dt, to = bond_fut$bond_matu, by = "quarter")
          true_cp_dt <- true_cp_dt[seq(1, length(true_cp_dt), by = 2)] }
          true_cp_dt <- c(head(true_cp_dt, -1) + bond_fut$sett, tail(true_cp_dt, 1))
          true_cp_dt <- c(bond_fut$option_matu, true_cp_dt)

          if(day_count_conv == 1){
            if(bond_fut$cp_f == 1){ stub <- as.numeric(diff(head(true_cp_dt, 2)))/
              as.numeric(bond_fut$curr_cp_dt - bond_fut$prev_cp_dt)
            } else {stub <- as.numeric(diff(head(true_cp_dt, 2)))/
              as.numeric(ceiling_date(bond_fut$prev_cp_dt + sett, "year") - floor_date(bond_fut$prev_cp_dt + sett, "year")) }
            cp_dt_2 <- stub + c(0, seq(1, length(true_cp_dt) - 2)/bond_fut$cp_f)
          } else if(day_count_conv == 2){ cp_dt_2 <- tail(as.numeric(true_cp_dt - first(true_cp_dt)), -1)/360
          } else if(day_count_conv == 3){ cp_dt_2 <- tail(as.numeric(true_cp_dt - first(true_cp_dt)), -1)/365
          } else {stub <- (max(0, 30 - as.numeric(format(true_cp_dt[1], "%d"))) +
                             as.numeric(round((floor_date(true_cp_dt[2], "months") - ceiling_date(true_cp_dt[1], "months"))/30))*30 +
                             min(30, as.numeric(format(true_cp_dt[2], "%d"))))/360
          cp_dt_2 <- stub + c(0, seq(1, length(true_cp_dt) - 2)/bond_fut$cp_f)}

          cp_dt_2 <- list(list(cp_dt_2))
          cf_matu <- bond_fut$Nomi*(1 + bond_fut$bond_cp/bond_fut$cp_f)
          cf_other <- split(rep(bond_fut$bond_cp/bond_fut$cp_f*bond_fut$Nomi, sapply(cp_dt_2, lengths) - 1),
                            rep(seq_along(cp_dt_2), sapply(cp_dt_2, lengths) - 1))

          esp <- function(x){exp(x[1] + x[2] + 0.5*(x[3]^2 + x[4]^2 + 2*prod(x[3:5]) )  )}

          esp_mix <- function(x){
            esp_mix <- x[9]*x[10]*esp(x[c(1, 3, 5, 7, 11)]) + x[9]*(1 - x[10])*esp(x[c(1, 4, 5, 8, 12)]) +
              (1 - x[9])*x[10]*esp(x[c(2, 3, 6, 7, 13)]) + (1 - x[9])*(1 - x[10])*esp(x[c(2, 4, 6, 8, 14)])
          }

          put <- function(x, KP){
            sigma <- x[3]^2 + x[4]^2 + 2*prod(x[3:5])
            d1_C <- (x[1] + x[2] + sigma - log(KP))/sqrt(sigma )
            d2_C <- d1_C - sqrt(sigma)
            put <- -esp(x)*pnorm(-d1_C) + KP*pnorm(-d2_C)
            if(cot_bond %in%c(1, 2)){put <- exp(-r*T)*put
            } else{put <- put}
          }

          put_mix <- function(x, KP){
            put_mix <- x[9]*x[10]*put(x[c(1, 3, 5, 7, 11)], KP) + x[9]*(1 - x[10])*put(x[c(1, 4, 5, 8, 12)], KP) +
              (1 - x[9])*x[10]*put(x[c(2, 3, 6, 7, 13)], KP) + (1 - x[9])*(1 - x[10])*put(x[c(2, 4, 6, 8, 14)], KP)
          }

          call <- function(x, KC){
            sigma <- x[3]^2 + x[4]^2 + 2*prod(x[3:5])
            d1_C <- (x[1] + x[2] + sigma - log(KC))/sqrt(sigma)
            d2_C <- d1_C - sqrt(sigma)
            call <- esp(x)*pnorm(d1_C) - KP*pnorm(d2_C)
            if(cot_bond %in%c(1, 2)){call <- exp(-r*T)*call
            } else{call <- call}
          }

          call_mix <- function(x, KC){
            call_mix <- x[9]*x[10]*call(x[c(1, 3, 5, 7, 11)], KC) + x[9]*(1 - x[10])*call(x[c(1, 4, 5, 8, 12)], KC) +
              (1 - x[9])*x[10]*call(x[c(2, 3, 6, 7, 13)], KC) + (1 - x[9])*(1 - x[10])*call(x[c(2, 4, 6, 8, 14)], KC)
          }

          PR <- matrix(seq(0.01, 0.49, 0.01), ncol = 1)

          if(cot_bond %in%c(1, 3)){
            model_prices <- function(x){
              return(list(model_call_price = call_mix(x, KC), model_put_price = put_mix(x, KP)))}
          } else {model_prices <- function(x){
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
            return(list(model_call_price = CALL, model_put_price = PUT))}
          }

          if(bond_fut$fut_matu < bond_fut$curr_cp_dt){
            allc <- bond_fut$acc_matu
          } else { allc <- bond_fut$acc_matu + bond_fut$bond_cp*bond_fut$Nomi/bond_fut$cp_f }
          C <- bond_fut$conv_factor*bond_call_prices
          P <- bond_fut$conv_factor*bond_put_prices
          KC <- bond_fut$conv_factor*bond_call_strikes
          KP <- bond_fut$conv_factor*bond_put_strikes
          T <- bond_fut$option_term
          FWD <- bond_fut$fut_price
          fwd_term <- bond_fut$res_term
          FWD_r <- exp(fwd_1*fwd_term)
          FWD_b <- (bond_fut$conv_factor*FWD + allc)/FWD_r

          esp_a_mix <- function(x){
            esp_a_mix <- x[5]*exp(x[1] + (x[3]^2)/2 ) + (1 - x[5])*exp(x[2] + (x[4]^2)/2 )
          }

          var_mix <- function(x){
            var_mix <- x[5]*exp(2*x[1] + 2*x[3]^2) + (1 - x[5])*exp(2*x[2] + 2*x[4]^2) - esp_a_mix(x)^2
          }

          MSE_mix <- function(x){
            MSE_mix <- sum((C - model_prices(x)$model_call_price)^2, na.rm = T) +
              sum((P - model_prices(x)$model_put_price)^2, na.rm = T) +
              (FWD_b - esp_a_mix(x[c(1, 2, 5, 6, 9)]))^2 +
              (FWD_r - esp_a_mix(x[c(3, 4, 7, 8, 10)]))^2 +
              (FWD_b*FWD_r - esp_mix(x))^2
            return(MSE_mix)}

          if( as.numeric(bond_charac_2$option_matu - bond_charac_2$start_date) > 90){
            volat_stir <- as.numeric(sqrt(stir_fut_p$params[5]*stir_fut_p$params[3]^2 +
                                            (1 - stir_fut_p$params[5])*stir_fut_p$params[4]^2 +
                                            (1 - stir_fut_p$params[5])*stir_fut_p$params[5]*(stir_fut_p$params[1] - stir_fut_p$params[2] )^2 ))
            volat_bond <- as.numeric(sqrt(bond_fut_rnd$params[5]*bond_fut_rnd$params[3]^2 +
                                            (1 - bond_fut_rnd$params[5])*bond_fut_rnd$params[4]^2 +
                                            (1 - bond_fut_rnd$params[5])*bond_fut_rnd$params[5]*(bond_fut_rnd$params[1] - bond_fut_rnd$params[2] )^2 ))
          } else {
            volat_stir <- as.numeric(stir_fut_p$params[5]*stir_fut_p$params[3] +
                                       (1 - stir_fut_p$params[5])*stir_fut_p$params[4])
            volat_bond <- as.numeric(bond_fut_rnd$params[5]*bond_fut_rnd$params[3] +
                                       (1 - bond_fut_rnd$params[5])*bond_fut_rnd$params[4])
          }

          m1 <- m2 <- s1 <- s2 <- rho_1 <- rho_2 <- rho_3 <- rho_4 <- SCE <- NA
          m3 <- m4 <- log(FWD_r)
          s3 <- s4 <- 0.25*sqrt(bond_fut$res_term)*volat_stir
          pi2 <- 0.5

          PARA <- as.matrix(data.frame(m1, m2, m3, m4, s1, s2, s3, s4, pi1 = PR[, 1], pi2,
                                       rho_1, rho_2, rho_3, rho_4, w1 = 0, w2 = 0, SCE))

          if(FWD_b != 1){
            lower <- c( rep((sign(1 - FWD_b)*0.5 + 1)*log(FWD_b), nb_log),
                        rep(0.5*volat_bond, nb_log),  rep(-1 + 1e-6, 4), 1e-6, 1e-6)
            upper <- c( rep((sign(FWD_b - 1)*0.5 + 1)*log(FWD_b), nb_log),
                        rep(1.2*volat_bond, nb_log), rep(1 - 1e-6, 4), 1 - 1e-6, 1 - 1e-6)
          } else {
            lower <- c( rep(0, nb_log), rep(0.5*volat_bond, nb_log), rep(-1 + 1e-6, 4), 1e-6, 1e-6)
            upper <- c( rep(0, nb_log), rep(1.2*volat_bond, nb_log), rep(1 - 1e-6, 4), 1 - 1e-6, 1 - 1e-6)
          }

          objective <- function(x){
            if(cot_bond %in% c(1, 3)){
              MSE_mix( c(x[1:2], m3, m4, x[3:4], s3, s4, PR[i, 1], pi2, x[5:8] ))
            } else{MSE_mix( c(x[1:2], m3, m4, x[3:4], s3, s4, PR[i, 1], pi2, x[5:10] )) }
          }

          set.seed(123)

          suppressWarnings({
            for (i in 1:nrow(PR)){
              if(cot_bond %in%c(1,3)){
                sol <- DEoptim(objective, lower[1:(length(lower) - 2)], upper = upper[1:(length(upper) - 2)],
                               DEoptim.control(trace = FALSE, NP = 80, itermax = 100))
                PARA[i, c("m1", "m2", "s1", "s2", "rho_1", "rho_2", "rho_3", "rho_4")] <- sol$optim$bestmem
              } else{
                sol <- DEoptim(objective, lower, upper = upper,
                               DEoptim.control(trace = FALSE, NP = 80, itermax = 100))
                PARA[i, c("m1", "m2", "s1", "s2", "rho_1", "rho_2", "rho_3", "rho_4", "w1", "w2")] <- sol$optim$bestmem }
              PARA[i, "SCE"] <- sol$optim$bestval
            }
          })

          PARA <- PARA[ !is.na(PARA[, "m1"]), ]
          if(nrow(PARA) > 1){
            param <- PARA[which.min(PARA[, "SCE"]), -ncol(PARA)]
          } else {param <- as.matrix(PARA[, -ncol(PARA)], nrow = nrow(PARA))}

          param[param == 0] <- 1e-4

          L <- U <- rep(0, length(param))

          L[sign(param) == -1] <- 1.2*param[sign(param) == -1]
          L[sign(param) == 1] <- 0.8*param[sign(param) == 1]
          U[sign(param) == -1] <- 0.8*param[sign(param) == -1]
          U[sign(param) == 1] <- 1.2*param[sign(param) == 1]

          if(cot_bond%in%c(1, 3)){
            L <- L[1: (length(L) - 2)]
            U <- U[1: (length(U) - 2)]
          } else{ U[(length(U) - 1):length(U)] <- pmin(U[(length(U) - 1):length(U)], 1)
          L[(length(L) - 1):length(L)] <- pmax(L[(length(L) - 1):length(L)], 0)}

          L[9] <- max(0, L[9])
          U[9] <- min(1, U[9])
          L[10] <- 0
          U[10] <- 1
          L[11:14] <- pmax(-1, L[11:14])
          U[11:14] <- pmin(1, U[11:14])

          CI <- c(L, -U)
          UI <- rbind(diag(length(L)), -diag(length(L)))

          if(cot_bond%in%c(1, 3)){
            param <- param[1: (length(param) - 2)]
          } else{param <- param }

          param <- matrix(param)

          suppressWarnings({
            solu <- constrOptim(param, MSE_mix, NULL, ui = UI, ci = CI, control = list(maxit = 5000),
                                mu = 1e-5, method = "Nelder-Mead")
          })

          model_p <- mapply("/", model_prices(c(solu$par)), bond_fut$conv_factor, SIMPLIFY = F)

          if(cot_bond%in%c(1,3) ){ params <- solu$par
          } else {params <- solu$par[1: (length(solu$par) - 2)]}

          marginal_bond <- params[which(colnames(PARA)%in%c("m1", "m2","s1", "s2", "pi1"))]
          marginal_repo <- params[which(colnames(PARA)%in%c("m3", "m4","s3", "s4", "pi2"))]

          sub <- function(x, y){ x[3]*dlnorm(y, meanlog = x[1], sdlog = x[2]) }
          PDF <- function(x, y){
            return(sub(x[c(1, 3, 5)], y) + sub(c(x[c(2, 4)], 1 - x[5]), y) ) }

          sub_2 <- function(x, y){ x[3]*plnorm(y, meanlog = x[1], sdlog = x[2]) }
          CDF <- function(x, y){
            return(sub_2(x[c(1, 3, 5)], y) + sub_2(c(x[c(2, 4)], 1 - x[5]), y) )}

          range_px_B <- range(c(KP, KC)/FWD_r)
          PX_B <- Reduce(seq, 1e3*range_px_B)*1e-3
          DNR_B <- PDF(marginal_bond, PX_B)

          esp_r <- esp_a_mix(marginal_repo)
          std_r <- sqrt(var_mix(marginal_repo))
          step_r <- std_r/300
          PX_r <- seq(esp_r - std_r, esp_r + std_r, step_r)
          range_px_r <- range(PX_r)
          DNR_r <- PDF(marginal_repo, PX_r)

          repo <- bond <- ""

          if(sum(rollmean(DNR_r, 2)*diff(PX_r), na.rm = T) < 1){

            x_axis <- 1e-2
            PX_r_2 <- PX_r
            range_px_r_2 <- range_px_r
            ratio <- 1
            while (ratio > 1e-5){
              integral <- sum(rollmean(PDF(marginal_repo, PX_r_2), 2)*diff(PX_r_2), na.rm = T)
              range_px_r_2 <- c(1 - x_axis, 1)*range_px_r_2
              PX_r_2 <- seq(min(range_px_r_2), max(range_px_r_2), step_r)
              integral_2 <- sum(rollmean(PDF(marginal_repo, PX_r_2), 2)*diff(PX_r_2), na.rm = T)
              ratio <- integral_2 - integral}

            ratio <- 1
            while (ratio > 1e-5){
              integral <- sum(rollmean(PDF(marginal_repo, PX_r_2), 2)*diff(PX_r_2), na.rm = T)
              range_px_r_2 <- c(1, 1 + x_axis)*range_px_r_2
              PX_r_2 <- seq(min(range_px_r_2), max(range_px_r_2), step_r)
              integral_2 <- sum(rollmean(PDF(marginal_repo, PX_r_2), 2)*diff(PX_r_2), na.rm = T)
              ratio <- integral_2 - integral}

            while (sum(rollmean(PDF(marginal_repo, PX_r_2), 2)*diff(PX_r_2), na.rm = T) < 0.9991){
              range_px_r_2 <- c(1 - x_axis, 1 + x_axis)*range_px_r_2
              PX_r_2 <- seq(min(range_px_r_2), max(range_px_r_2), step_r)}

            extension <- diff(range(PX_r_2))/diff(range(PX_r))

            DNR_2_r <- PDF(marginal_repo, PX_r_2)
            NCDF_r <- CDF(marginal_repo, PX_r_2)

            if(DNR_2_r[1] <= DNR_2_r[2] &
               DNR_2_r[length(DNR_2_r) - 1] >= DNR_2_r[length(DNR_2_r)] &
               min(DNR_2_r)%in%DNR_2_r[c(1, length(DNR_2_r))]){

              repo_rate <- function(x){repo_rate <- exp(x*fwd_term)}

              PX_r_3 <- log(PX_r_2)/fwd_term

              sub_r <- function(x, y){
                x[3]*dlnorm( repo_rate(y[-1]), meanlog = x[1], sdlog = x[2])*diff(repo_rate(y))/diff(y)  }

              PDF_r <- function(x, y){
                return(sub_r(x[c(1, 3, 5)], y) + sub_r(c(x[c(2, 4)], 1 - x[5]), y) ) }

              DNR_repo <- PDF_r(marginal_repo, PX_r_3)

              remove <- which(is.na(DNR_repo))

              if(length(remove) > 0){
                PX_r_4 <- PX_r_3[-remove]
                DNR_repo_2 <- DNR_repo[-remove]
              } else { PX_r_4 <- PX_r_3
              DNR_repo_2 <- DNR_repo}

              df_r <- data.frame(price = PX_r_4[-1], density = DNR_repo_2)
              cdf_r <- data.frame(price = PX_r_4[-c(1, 2)], cdf = cumsum(rollmean(DNR_repo_2, 2)*diff(PX_r_4[-1])))

              thres <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 0.995, 0.999)

              if(length(which(cdf_r$cdf > last(thres))) > 0 & length(which(cdf_r$cdf < first(thres))) > 0){

                quantiles <- list()
                for (j in 1:length(thres)){
                  quantiles[[j]] <- mean(df_r$price[c(min(which(cdf_r$cdf > thres[j] - 1e-3)),
                                                      max(which(cdf_r$cdf < thres[j] + 1e-3)))])}

                qt <- data.frame(quantiles) %>% rename_with(~paste0("q", 100*thres))

                E_r <- sum(rollmean(PX_r_4[-1]*DNR_repo_2, 2)*diff(PX_r_4[-1]))
                moments_r <- function(x){ return(sum(rollmean(DNR_repo_2*(PX_r_4[-1] - E_r)^x , 2)*diff(PX_r_4[-1])))}
                SD_r <- sqrt(moments_r(2))
                SK_r <- moments_r(3)/SD_r^3
                KU_r <- moments_r(4)/SD_r^4
                moments_r <- c(mean = E_r, stddev = SD_r, skewness = SK_r, kurtosis = KU_r)
                mode_r <- PX_r_4[which.max(DNR_repo_2)]

                graph <- PX_r_4 >= qt$q0.1 & PX_r_4 <= qt$q99.9
                PX_graph_repo <- PX_r_4[graph]
                DNR_graph_repo <- DNR_repo_2[graph]
                NCDF_graph_repo <- cdf_r$cdf[graph]
                df_graph_repo <- data.frame(rate = PX_graph_repo, density = DNR_graph_repo)
                cdf_graph_repo <- data.frame(rate = PX_graph_repo, cdf = NCDF_graph_repo)

                pdf_r <- ggplot() + geom_line(data = df_graph_repo, aes(x = rate, y = density)) +
                  labs(x = paste0("repo rate (%) as of ",  bond_charac_2$start_date,
                                  " from ", bond_charac_2$option_matu, " to ", bond_charac_2$fut_matu),
                       y = "probability density") + theme_bw() +
                  theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm")) +
                  labs(title = paste0("Forward repo rate on ", country, " ",100*bond_cp, "% ", bond_matu),
                       subtitle = paste0("Risk Neutral Probability Density for a mixture of ", nb_log, " lognormals")) +
                  scale_x_continuous(labels = scales::percent)

                ncdf_r <- ggplot() + geom_line(data = cdf_graph_repo, aes(x = rate, y = cdf)) +
                  labs(x = paste0("repo rate (%) as of ",  bond_charac_2$start_date,
                                  " from ", bond_charac_2$option_matu, " to ", bond_charac_2$fut_matu),
                       y = "cumulative probability") + theme_bw() +
                  theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm")) +
                  labs(title =  paste0("Forward repo rate on ", country, " ", 100*bond_cp, "% ", bond_matu),
                       subtitle = paste0("Risk Neutral Cumulative Probability for a mixture of ", nb_log, " lognormals")) +
                  scale_x_continuous(labels = scales::percent)

                repo = list(moments_repo = moments_r, mode_repo = mode_r, discretized_rnd_repo = tibble(domain = PX_graph_repo, rnd = DNR_graph_repo),
                            rnd_plot_repo = pdf_r, cdf_plot_repo = ncdf_r)

              }
            }
          }

          if(sum(rollmean(PDF(marginal_bond, PX_B), 2)*diff(PX_B), na.rm = T) < 1){
            x_axis <- 1e-2
            PX_B_2 <- PX_B
            range_px_b_2 <- range_px_B
            ratio <- 1
            while (ratio > 1e-5){
              integral <- sum(rollmean(PDF(marginal_bond, PX_B_2), 2)*diff(PX_B_2), na.rm = T)
              range_px_b_2 <- c(1 - x_axis, 1)*range_px_b_2
              PX_B_2 <- Reduce(seq, 1e3*range_px_b_2)*1e-3
              integral_2 <- sum(rollmean(PDF(marginal_bond, PX_B_2), 2)*diff(PX_B_2), na.rm = T)
              ratio <- integral_2 - integral}

            ratio <- 1
            while (ratio > 1e-5){
              integral <- sum(rollmean(PDF(marginal_bond, PX_B_2), 2)*diff(PX_B_2), na.rm = T)
              range_px_b_2 <- c(1, 1 + x_axis)*range_px_b_2
              PX_B_2 <- Reduce(seq, 1e3*range_px_b_2)*1e-3
              integral_2 <- sum(rollmean(PDF(marginal_bond, PX_B_2), 2)*diff(PX_B_2), na.rm = T)
              ratio <- integral_2 - integral}

            while (sum(rollmean(PDF(marginal_bond, PX_B_2), 2)*diff(PX_B_2), na.rm = T) < 0.9991){
              range_px_b_2 <- c(1 - x_axis, 1 + x_axis)*range_px_b_2
              PX_B_2 <- Reduce(seq, 1e3*range_px_b_2)*1e-3}

            extension <- diff(range(PX_B_2))/diff(range(PX_B))

            DNR_2 <- PDF(marginal_bond, PX_B_2)

            NCDF <- CDF(marginal_bond, PX_B_2)

            if(DNR_2[1] <= DNR_2[2] &
               DNR_2[length(DNR_2) - 1] >= DNR_2[length(DNR_2)] &
               min(DNR_2)%in%DNR_2[c(1, length(DNR_2))]){

              dirty <- function(x){
                dcf <- mapply("/", list(c(unlist(cf_other), cf_matu)),
                              mapply("^", 1 + x, list(unlist(cp_dt_2)), SIMPLIFY = F),
                              SIMPLIFY = F)
                dirty <- unlist(lapply(dcf, sum))}

              tri <- function(x){
                tri <- mapply(xirr, cf = mapply(c, -x, cf_other, cf_matu, SIMPLIFY = F),
                              tau = mapply(c, 0, mapply(unlist, cp_dt_2, SIMPLIFY = F), SIMPLIFY = F) )
              }

              PX_B_3 <- rev(tri(PX_B_2))
              sub_3 <- function(x, y){
                x[3]*dlnorm( dirty(y[-1]), meanlog = x[1], sdlog = x[2])*(-diff(dirty(y)))/diff(y)            }

              PDF_y <- function(x, y){
                return(sub_3(x[c(1, 3, 5)], y) + sub_3(c(x[c(2, 4)], 1 - x[5]), y) ) }

              DNR_y <- PDF_y(params, PX_B_3)
              remove <- which(is.na(DNR_y))

              if(length(remove) > 0){
                PX_B_4 <- PX_B_3[-remove]
                DNR_y_2 <- DNR_y[-remove]
              } else { PX_B_4 <- PX_B_3
              DNR_y_2 <- DNR_y}

              df_y <- data.frame(price = PX_B_4[-1], density = DNR_y_2)
              cdf_y <- data.frame(price = PX_B_4[-c(1,2)], cdf = cumsum(rollmean(DNR_y_2, 2)*diff(PX_B_4[-1])))

              thres <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 0.995, 0.999)

              if(length(which(cdf_y$cdf > last(thres))) > 0 & length(which(cdf_y$cdf < first(thres))) > 0){

                quantiles <- list()
                for (j in 1:length(thres)){
                  quantiles[[j]] <- mean(df_y$price[c(min(which(cdf_y$cdf > thres[j] - 1e-3)), max(which(cdf_y$cdf < thres[j] + 1e-3)))])}

                qt <- data.frame(quantiles) %>% rename_with(~paste0("q", 100*thres))

                E_y <- sum(rollmean(PX_B_4[-1]*DNR_y_2, 2)*diff(PX_B_4[-1]))
                moments_y <- function(x){ return(sum(rollmean(DNR_y_2*(PX_B_4[-1] - E_y)^x , 2)*diff(PX_B_4[-1])))}
                SD_y <- sqrt(moments_y(2))
                SK_y <- moments_y(3)/SD_y^3
                KU_y <- moments_y(4)/SD_y^4
                moments_y <- c(mean = E_y, stddev = SD_y, skewness = SK_y, kurtosis = KU_y)
                mode_y <- PX_B_4[which.max(DNR_y_2)]

                graph <- PX_B_4 >= qt$q0.1 & PX_B_4 <= qt$q99.9
                PX_graph <- PX_B_4[graph]
                DNR_graph <- DNR_y_2[graph]
                NCDF_graph <- cdf_y$cdf[graph]
                df_graph <- data.frame(price = PX_graph, density = DNR_graph)
                cdf_graph <- data.frame(price = PX_graph, cdf = NCDF_graph)

                pdf_y <- ggplot() + geom_line(data = df_graph, aes(x = price, y = density)) +
                  labs(x = paste0("Bond yield to maturity (%) on ", bond_charac_2$option_matu,
                                  " as of ",  bond_charac_2$start_date),
                       y = "probability density") + theme_bw() +
                  theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm")) +
                  labs(title = paste0("Forward Bond yield on ", country, " ", 100*bond_cp, "% ", bond_matu),
                       subtitle = paste0("Risk Neutral Probability Density for a mixture of ", nb_log, " lognormals")) +
                  scale_x_continuous(labels = scales::percent)

                ncdf_y <- ggplot() + geom_line(data = cdf_graph, aes(x = price, y = cdf)) +
                  labs(x = paste0("Bond yield to maturity (%) on ", bond_charac_2$option_matu,
                                  " as of ",  bond_charac_2$start_date),
                       y = "cumulative probability") + theme_bw() +
                  theme(legend.position = "none", plot.margin = margin(.8,.5,.8,.5, "cm")) +
                  labs(title = paste0("Forward Bond yield on ", country, " ", 100*bond_cp, "% ", bond_matu),
                       subtitle = paste0("Risk Neutral Cumulative Probability for a mixture of ", nb_log, " lognormals")) +
                  scale_x_continuous(labels = scales::percent)

                bond = list(moments_ytm = moments_y, mode_ytm = mode_y, discretized_rnd_ytm = tibble(domain = PX_graph, rnd = DNR_graph),
                            rnd_plot_ytm = pdf_y, cdf_plot_ytm = ncdf_y)

              }
            }
          }

          if(length(repo)  == 0 ){ message("impossible to retrieve a density for the forward repo price")}
          if(length(bond)  == 0 ){ message("impossible to retrieve a density for the forward bond price")}

          all <- c(params_bond = data.frame(unlist(marginal_bond)),
                   params_repo = data.frame(unlist(marginal_repo)), repo, bond, CV = solu$convergence)
          return(all)

        }
      } else {message("input dates are not consistent")}
    } else{ message ("please enter STIR options with maturity close to bond options' maturity (distance below 30 days)")}
  } else {message("inputs do not have the required length")}
}
