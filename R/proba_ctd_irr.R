#' proba_ctd_irr
#'
#' @param bond_call_prices a vector of call prices on a bond futures, in numeric format
#' @param bond_call_strikes a vector of call strikes attached to the call prices, in numeric format
#' @param bond_put_prices a vector of put prices on the same bond futures, in numeric format
#' @param bond_put_strikes a vector of put strikes attached to the put prices, in numeric format
#' @param stir_call_prices a vector of call prices on a STIR futures whose maturity date is near the maturity date of the bond futures (maturity distance below 30 days), in numeric format
#' @param stir_call_strikes a vector of call strikes attached to the call prices, in numeric format
#' @param stir_put_prices a vector of put prices on the same STIR futures, in numeric format
#' @param stir_put_strikes a vector of put strikes attached to the put prices, in numeric format
#' @param r a number for the riskfree spot rate whose maturity is equal to the options' maturity, in numeric format
#' @param r_2 a number for the spot repo funding rate of the underlying bond, with maturity equal to the options' maturity, in numeric format
#' @param r_3 a number for the spot repo funding rate of the underlying bond, with maturity equal to the futures' maturity, in numeric format
#' @param day_count_conv a number for the day count convention, 1 for ACT/ACT, 2 for ACT/360, 3 for ACT/365 and 4 for 30/360, in numeric format
#' @param cot_bond a number for the bond options' style, 1 for European options, 2 for American options and 3 for American options with futures-style margin, in numeric format
#' @param bond_fut_price a number for the bond futures' price on calibration date, in numeric format
#' @param bond_cp a vector of the coupon rates for the bonds in the delivery basket, in numeric format
#' @param bond_cp_f a vector of the corresponding frequencies of coupon payment, either 1 for annual payment or 2 for semiannual payment, in numeric format
#' @param bond_conv_factor a vector of the corresponding conversion factors for the bonds in the delivery basket, in numeric format
#' @param bond_ytm a vector of the corresponding yields to maturity at observation date for the bonds in the delivery basket, in numeric format
#' @param sett a number for the number of days between the ex-coupon date and the coupon payment date of the current Cheapest-to-Deliver Bond, in numeric format
#' @param Nomi a single number for the value of the principal of the bonds in the delivery basket, in numeric format (100 by default)
#' @param bond_matu a vector of the corresponding maturity dates for the bonds in the basket of deliverable bonds, in Date format
#' @param bond_ISIN a vector of the corresponding ISIN codes for the bonds in the delivery basket of the bond futures, in character format
#' @param bond_fut_matu a date for the maturity date of the bond futures contract, in Date format
#' @param bond_option_matu a date for the maturity date of the bond options, in Date format
#' @param cot_stir a number for the STIR options' style, 1 for European options, 2 for American options and 3 for American options with futures-style margin, in numeric format
#' @param stir_fut_price a number for the STIR futures' price on calibration date, in numeric format
#' @param stir_fut_matu a date for the maturity date of the STIR futures contract, in Date format
#' @param stir_option_matu a date for the maturity date of the STIR options, in Date format
#' @param start_date a date for the observation date, in Date format
#' @param term_stir a number for the term to maturity, in months, of the STIR at observation date, in numeric format
#'
#' @returns for the bonds in the delivery basket, their ISIN in character format and their probability of being the CtD bond at options' maturity based on each deliverable bond's joint distribution of bond yield and repo rate, in numeric format
#' @export
#'
#' @importFrom stats approx constrOptim density pnorm qlnorm
#' @importFrom utils head tail
#' @importFrom MASS mvrnorm
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
#' proba_ctd_irr(
#' c(19.12, 18.12, 17.12, 16.14, 15.14, 14.16, 13.18, 12.20, 11.22, 10.24, 9.30, 8.34,
#' 7.42, 6.52, 5.66, 4.84, 4.06, 3.32, 2.66, 2.08, 1.64, 1.28, 1.00, 0.78, 0.60, 0.46,
#' 0.34, 0.26, 0.20, 0.14, 0.10, 0.08, 0.06, 0.04, 0.04, 0.02, 0.02, 0.02, 0.02, 0.02,
#' 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02),
#' seq(84, 131),
#' c(0.04, 0.04, 0.04, 0.06, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.22, 0.26, 0.34, 0.44,
#' 0.58, 0.76, 0.98, 1.24, 1.58, 2.00, 2.56, 3.20, 3.92, 4.70, 5.52, 6.38, 7.26, 8.18,
#' 9.12, 10.06, 11.02, 12.00, 12.98, 13.96, 14.96, 15.94, 16.94, 17.94, 18.92, 19.92,
#' 20.92, 21.92, 22.92, 23.92, 24.92, 25.92, 26.92, 27.92),
#' seq(84, 131),
#' c(3.1850, 3.0600, 2.9350, 2.8100, 2.6850, 2.5600, 2.4350, 2.3100, 2.1850, 2.0625,
#' 1.9375, 1.8125, 1.6875, 1.5625, 1.4375, 1.3150, 1.1900, 1.1275, 1.0650, 0.9400,
#' 0.8800, 0.8175, 0.7550, 0.6950, 0.6325, 0.5725, 0.5125, 0.4525, 0.3925, 0.3350,
#' 0.2775, 0.2225, 0.1700, 0.1250, 0.0875, 0.0575, 0.0325, 0.0150, 0.0075, 0.0050,
#' 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025,
#' 0.0025, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
#' 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
#' 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000),
#' c(seq(94, 96, 0.125), 96.0625, 96.125, seq(96.25, 99.25, 0.0625), 99.3750,
#' 99.5000, 99.6250,  99.6875, 99.7500,  99.8750, 100.0000, 101.0000),
#' c(0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0025,
#' 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0050, 0.0050, 0.0050, 0.0050, 0.0050,
#' 0.0075, 0.0075, 0.0075, 0.0100, 0.0100, 0.0125, 0.0150, 0.0175, 0.0200, 0.0250,
#' 0.0300, 0.0375, 0.0475, 0.0650, 0.0900, 0.1225, 0.1600, 0.2050, 0.2600, 0.3200,
#' 0.3800, 0.4425, 0.5050, 0.5675, 0.6300, 0.6925, 0.7550, 0.8175, 0.8800, 0.9425,
#' 1.0050, 1.0650, 1.1275, 1.1900, 1.2525, 1.3150, 1.3775, 1.4400, 1.5025, 1.5650,
#' 1.6275, 1.6900, 1.7525, 1.8150, 1.8775, 1.9400, 2.0025, 2.0650, 2.1900, 2.3150,
#' 2.4400, 2.5025, 2.5650, 2.6900, 2.8150, 3.8150),
#' c(seq(94, 96, 0.125), 96.0625, 96.125, seq(96.25, 99.25, 0.0625), 99.3750,
#' 99.5000, 99.6250,  99.6875, 99.7500,  99.8750, 100.0000, 101.0000),
#' 0.02303,
#' 0.02303,
#' 0.02388,
#' 1,
#' 3,
#' 103.08,
#' c(0.000, 0.018, 0.018, 0.025, 0.029),
#' rep(1, 5),
#' c(0.365252, 0.643086, 0.643086, 0.751530, 0.810737),
#' c(0.03811, 0.03816, 0.03806, 0.03813, 0.03812),
#' 2,
#' Nomi = 100,
#' as.Date(c("2052-08-15", "2053-08-15", "2053-08-15", "2054-08-15", "2056-08-15")),
#' c("DE0001102572", "DE0001102614", "DE0001030757", "DE000BU2D004", "DE000BU2D012"),
#' as.Date("2026-12-08"),
#' as.Date("2026-11-20"),
#' 3,
#' 97.185,
#' as.Date("2026-12-14"),
#' as.Date("2026-12-14"),
#' as.Date("2026-08-31"),
#' 3)
#' }
#'

proba_ctd_irr <- function(bond_call_prices, bond_call_strikes, bond_put_prices, bond_put_strikes,
                          stir_call_prices, stir_call_strikes, stir_put_prices, stir_put_strikes,
                          r, r_2, r_3, day_count_conv, cot_bond, bond_fut_price, bond_cp, bond_cp_f,
                          bond_conv_factor, bond_ytm, sett, Nomi = 100, bond_matu, bond_ISIN, bond_fut_matu,
                          bond_option_matu, cot_stir, stir_fut_price, stir_fut_matu, stir_option_matu,
                          start_date, term_stir){

  if(length(r) == 1 & length(r_2) == 1 & length(r_3) == 1 & length(day_count_conv) == 1 & length(cot_bond) == 1 &
     length(bond_fut_price) == 1 & length(bond_fut_matu) == 1 & length(bond_option_matu) == 1 & length(start_date) == 1 &
     length(sett) == 1 & length(cot_stir) == 1 & length(stir_fut_price) == 1 & length(stir_fut_matu) == 1 &
     length(stir_option_matu) ==1 & length(term_stir) == 1 & length(stir_call_prices) > 1 & length(stir_call_strikes) > 1 &
     length(stir_put_prices) > 1& length(stir_put_strikes) > 1 & length(bond_call_prices) > 1 & length(bond_call_strikes) > 1 &
     length(bond_put_prices) > 1 & length(bond_put_strikes) > 1 & length(bond_ISIN) > 1 &
     identical(length(bond_ISIN), length(bond_cp), length(bond_cp_f), length(bond_matu),
               length(bond_conv_factor), length(bond_ytm))){


    simulate_mixture <- function(x, U) {
      U1 <- pmin(U/x[5], 1)
      U2 <- pmax((U - x[5])/(1 - x[5]), 0)
      ifelse(U < x[5], qlnorm(U1, meanlog = x[1], sdlog = x[3]),
             qlnorm(U2, meanlog = x[2], sdlog = x[4]) )
    }

    simulate_correlated_returns <- function(n, params, Sigma) {
      Z <- mvrnorm(n, mu = rep(0, nrow(params)), Sigma, tol = 1e-06, empirical = F)
      U <- pnorm(Z)
      simR <- data.frame(matrix(nrow = n, ncol = nrow(params)))
      for (j in 1:nrow(params)) {
        p <- params %>% filter(row_number() == j) %>% unlist
        F_T <- simulate_mixture(p, U[ ,j])
        simR[, j] <- F_T }
      simR
    }

    stir_fut_matu <- as.Date(stir_fut_matu)
    stir_option_matu <- as.Date(stir_option_matu)
    bond_fut_matu <- as.Date(bond_fut_matu)
    bond_option_matu <- as.Date(bond_option_matu)
    bond_fut_price <- as.numeric(bond_fut_price)
    stir_fut_price <- as.numeric(stir_fut_price)
    start_date <- as.Date(start_date)
    nb_log <- 2
    sett <- sett
    Nomi <- Nomi
    day_count_conv <- day_count_conv

    if( abs(as.numeric(stir_option_matu) - as.numeric(bond_option_matu)) < 30){

      if(start_date < bond_option_matu & bond_option_matu <= bond_fut_matu &
         start_date < stir_fut_matu & stir_option_matu <= stir_fut_matu &
         length(which(as.Date(bond_fut_matu) - as.Date(bond_matu) > 0)) == 0){

        bond_fut_rnd <- bond_future_price(bond_call_prices, bond_call_strikes, bond_put_prices, bond_put_strikes,
                                          nb_log, r, day_count_conv, cot_bond, bond_matu[1], bond_fut_price,
                                          bond_fut_matu, bond_option_matu, start_date)

        stir_fut_rnd <- stir_future_price(stir_call_prices, stir_call_strikes, stir_put_prices,
                                          stir_put_strikes, nb_log, r, day_count_conv, cot_stir,
                                          stir_fut_price, stir_fut_matu, stir_option_matu, start_date)

        if(length(bond_fut_rnd) > 0 & length(stir_fut_rnd) > 0){

          deliverables <- data.frame(bond_ISIN, bond_cp, bond_matu, bond_conv_factor, bond_cp_f, bond_ytm, Nomi,
                                     bond_fut_price, bond_option_matu, start_date, sett, bond_fut_matu) %>%
            rename_with(~c("ISIN", "coupon", "bond_matu", "conv_factor", "cp_freq", "ytm", "Nomi",
                           "bond_fut_price", "bond_option_matu", "start_date", "sett", "bond_fut_matu")) %>%
            mutate_at(c("bond_matu", "bond_option_matu", "bond_fut_matu", "start_date"), as.Date)

          bond_fut <- marginal_bond <- marginal_repo <- true_cp_dt <- cp_dt_2 <- cf_matu <- cf_other <- params_each <- correl <- corr_matrix <- joint <- list()

          for (k in 1:nrow(deliverables)){

            bond_fut[[k]] <- deliverables[k, ] %>%
              mutate(prev_cp_dt = as.Date(paste0(format(bond_option_matu, "%Y"), "-", format(bond_matu, "%m-%d"))))

            if(bond_fut[[k]]$cp_f == 1){bond_fut[[k]] <- bond_fut[[k]] %>%
              mutate_at("prev_cp_dt", ~as.Date(ifelse(bond_option_matu < ., . %m-% years(1), .))) %>%
              mutate(curr_cp_dt = prev_cp_dt %m+% years(1), next_cp_dt = curr_cp_dt %m+% years(1))
            } else { bond_fut[[k]] <- bond_fut[[k]] %>%
              mutate_at("prev_cp_dt", ~as.Date(ifelse(bond_option_matu - . < - months(6), . %m-% years(1),
                                                      ifelse(bond_option_matu - . < 0, . %m-% months(6), .)))) %>%
              mutate(curr_cp_dt = prev_cp_dt %m+% months(6), next_cp_dt = curr_cp_dt %m+% months(6))}

            if(day_count_conv == 1){
              bond_fut[[k]] <- bond_fut[[k]] %>% mutate(option_term = as.numeric(bond_option_matu - start_date)/
                                                          as.numeric(ceiling_date(bond_option_matu, "year") - floor_date(start_date, "year")),
                                                        res_term = as.numeric(bond_fut_matu - bond_option_matu)/
                                                          as.numeric(ceiling_date(bond_fut_matu, "year") - floor_date(bond_option_matu, "year") ))
            } else if(day_count_conv == 2){
              bond_fut[[k]] <- bond_fut[[k]] %>% mutate(option_term = as.numeric(bond_option_matu - start_date)/360,
                                                        res_term = as.numeric(bond_fut_matu - bond_option_matu)/360)
            } else if(day_count_conv == 3){
              bond_fut[[k]] <- bond_fut[[k]] %>% mutate(option_term =as.numeric(bond_option_matu - start_date)/365,
                                                        res_term = as.numeric(bond_fut_matu - bond_option_matu)/365)
            } else {
              bond_fut[[k]] <- bond_fut[[k]] %>% mutate(stub_1 = max(0, 30 - as.numeric(format(start_date, "%d"))),
                                                        stub_2 = min(30, as.numeric(format(bond_option_matu, "%d"))),
                                                        plain_months = round((as.numeric(floor_date(bond_option_matu, "months") -
                                                                                           ceiling_date(start_date, "months")))/30),
                                                        option_term = (stub_1 + stub_2 + plain_months*30)/360,
                                                        stub_1_res = max(0, 30 - as.numeric(format(bond_option_matu + sett, "%d"))),
                                                        stub_2_res = min(30, as.numeric(format(bond_fut_matu, "%d"))),
                                                        plain_months_res = round(as.numeric(floor_date(bond_fut_matu, "months") -
                                                                                              ceiling_date(bond_option_matu + sett, "months") )/30),
                                                        res_term = (stub_1_res + stub_2_res + max(0, plain_months_res)*30)/360)}

            bond_fut[[k]] <- bond_fut[[k]] %>% mutate(res_term_2 = 0)


            if(bond_fut[[k]]$bond_fut_matu < bond_fut[[k]]$curr_cp_dt){
              if(day_count_conv == 1){
                bond_fut[[k]] <- bond_fut[[k]] %>% mutate(acc_matu = Nomi*coupon*as.numeric(bond_fut_matu - prev_cp_dt - sett)/
                                                            as.numeric(ceiling_date(bond_fut_matu, "year") - floor_date(bond_fut_matu, "year") ))
              } else if(day_count_conv == 2) {
                bond_fut[[k]] <- bond_fut[[k]] %>% mutate(acc_matu = Nomi*coupon*as.numeric(bond_fut_matu - prev_cp_dt - sett)/360)
              } else if(day_count_conv == 3){
                bond_fut[[k]] <- bond_fut[[k]] %>% mutate(acc_matu = Nomi*coupon*as.numeric(bond_fut_matu - prev_cp_dt - sett)/365)
              } else{
                bond_fut[[k]] <- bond_fut[[k]] %>% mutate(stub_1 = max(0, 30 - as.numeric(format(prev_cp_dt + sett, "%d"))),
                                                          stub_2 = min(30, as.numeric(format(bond_fut_matu, "%d"))),
                                                          plain_months = round(as.numeric(floor_date(bond_fut_matu, "months") -
                                                                                            ceiling_date(prev_cp_dt + sett, "months") )/30),
                                                          acc_matu = Nomi*coupon/360*(stub_2 + stub_1 + max(0, plain_months)*30)) }
            } else{
              if(day_count_conv == 1){
                bond_fut[[k]] <- bond_fut[[k]] %>% mutate(res_term_2 = as.numeric(bond_fut_matu - curr_cp_dt - sett)/
                                                            as.numeric(ceiling_date(bond_fut_matu, "year") - floor_date(bond_fut_matu, "year") ),
                                                          acc_matu = Nomi*coupon*res_term_2 )
              } else if(day_count_conv == 2) {
                bond_fut[[k]] <- bond_fut[[k]] %>% mutate(res_term_2 = as.numeric(bond_fut_matu - curr_cp_dt - sett)/360,
                                                          acc_matu = Nomi*coupon*res_term_2)
              } else if(day_count_conv == 3){
                bond_fut[[k]] <- bond_fut[[k]] %>% mutate(res_term_2 = as.numeric(bond_fut_matu - curr_cp_dt - sett)/365,
                                                          acc_matu = Nomi*coupon*res_term_2)
              } else{
                bond_fut[[k]] <- bond_fut[[k]] %>% mutate(stub_1 = max(0, 30 - as.numeric(format(curr_cp_dt + sett, "%d"))),
                                                          stub_2 = min(30, as.numeric(format(bond_fut_matu, "%d"))),
                                                          plain_months = round(as.numeric(floor_date(bond_fut_matu, "months") -
                                                                                            ceiling_date(curr_cp_dt + sett, "months") )/30),
                                                          res_term_2 = (stub_1 + stub_2 + max(0, plain_months)*30)/360,
                                                          acc_matu = Nomi*coupon*res_term_2)}
            }

            if(bond_fut[[k]]$cp_f == 1){ true_cp_dt <- seq(from = bond_fut[[k]]$curr_cp_dt, to = bond_fut[[k]]$bond_matu, by = "year")
            } else { true_cp_dt <- seq(from = bond_fut[[k]]$curr_cp_dt, to = bond_fut[[k]]$bond_matu, by = "quarter")
            true_cp_dt <- true_cp_dt[seq(1, length(true_cp_dt), by = 2)] }
            true_cp_dt <- c(head(true_cp_dt, -1) + bond_fut[[k]]$sett, tail(true_cp_dt, 1))
            true_cp_dt <- c(bond_fut[[k]]$bond_option_matu, true_cp_dt)

            if(day_count_conv == 1){
              if(bond_fut[[k]]$cp_freq == 1){ stub <- as.numeric(diff(head(true_cp_dt, 2)))/
                as.numeric(bond_fut[[k]]$curr_cp_dt - bond_fut[[k]]$prev_cp_dt)
              } else {stub <- as.numeric(diff(head(true_cp_dt, 2)))/
                as.numeric(ceiling_date(bond_fut[[k]]$prev_cp_dt + sett, "year") - floor_date(bond_fut[[k]]$prev_cp_dt + sett, "year")) }
              cp_dt_2[[k]] <- stub + c(0, seq(1, length(true_cp_dt) - 2)/bond_fut[[k]]$cp_freq)
            } else if(day_count_conv == 2){ cp_dt_2[[k]] <- tail(as.numeric(true_cp_dt - first(true_cp_dt)), -1)/360
            } else if(day_count_conv == 3){ cp_dt_2[[k]] <- tail(as.numeric(true_cp_dt - first(true_cp_dt)), -1)/365
            } else {stub <- (max(0, 30 - as.numeric(format(true_cp_dt[1], "%d"))) +
                               as.numeric(round((floor_date(true_cp_dt[2], "months") - ceiling_date(true_cp_dt[1], "months"))/30))*30 +
                               min(30, as.numeric(format(true_cp_dt[2], "%d"))))/360
            cp_dt_2[[k]] <- stub + c(0, seq(1, length(true_cp_dt) - 2)/bond_fut[[k]]$cp_freq)}

            cp_dt_2[[k]] <- list(list(cp_dt_2[[k]]))
            cf_matu[[k]] <- bond_fut[[k]]$Nomi*(1 + bond_fut[[k]]$coupon/bond_fut[[k]]$cp_freq)
            cf_other[[k]] <- split(rep(bond_fut[[k]]$coupon/bond_fut[[k]]$cp_freq*bond_fut[[k]]$Nomi, sapply(cp_dt_2[[k]], lengths) - 1),
                                   rep(seq_along(cp_dt_2[[k]]), sapply(cp_dt_2[[k]], lengths) - 1))


            rate_table <- data.frame(term = bond_fut[[k]]$option_term + c(0, bond_fut[[k]]$res_term), rates = c(r_2, r_3)) %>%
              mutate(d_fact = term*rates)
            fwd_1 <- diff(rate_table$d_fact)/diff(rate_table$term)

            esp <- function(x){exp(x[1] + x[2] + 0.5*(x[3]^2 + x[4]^2 + 2*prod(x[3:5]) )  )}

            esp_mix <- function(x){
              esp_mix <- x[9]*x[10]*esp(x[c(1, 3, 5, 7, 11)]) +
                x[9]*(1 - x[10])*esp(x[c(1, 4, 5, 8, 12)]) +
                (1 - x[9])*x[10]*esp(x[c(2, 3, 6, 7, 13)]) +
                (1 - x[9])*(1 - x[10])*esp(x[c(2, 4, 6, 8, 14)])
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
              put_mix <- x[9]*x[10]*put(x[c(1, 3, 5, 7, 11)], KP) +
                x[9]*(1 - x[10])*put(x[c(1, 4, 5, 8, 12)], KP) +
                (1 - x[9])*x[10]*put(x[c(2, 3, 6, 7, 13)], KP) +
                (1 - x[9])*(1 - x[10])*put(x[c(2, 4, 6, 8, 14)], KP)
            }

            call <- function(x, KC){
              sigma <- x[3]^2 + x[4]^2 + 2*prod(x[3:5])
              d1_C <- (x[1] + x[2] + sigma - log(KC))/sqrt(sigma)
              d2_C <- d1_C - sqrt(sigma)
              call <- esp(x)*pnorm(d1_C) - KC*pnorm(d2_C)
              if(cot_bond %in%c(1, 2)){call <- exp(-r*T)*call
              } else{call <- call}
            }

            call_mix <- function(x, KC){
              call_mix <- x[9]*x[10]*call(x[c(1, 3, 5, 7, 11)], KC) +
                x[9]*(1 - x[10])*call(x[c(1, 4, 5, 8, 12)], KC) +
                (1 - x[9])*x[10]*call(x[c(2, 3, 6, 7, 13)], KC) +
                (1 - x[9])*(1 - x[10])*call(x[c(2, 4, 6, 8, 14)], KC)
            }

            PR <- matrix(seq(0.01, 0.99, 0.03), ncol = 1)

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


            allc <- bond_fut[[k]]$acc_matu + bond_fut[[k]]$coupon*bond_fut[[k]]$Nomi/bond_fut[[k]]$cp_freq*
              max(0, as.numeric(bond_fut[[k]]$bond_fut_matu - bond_fut[[k]]$curr_cp_dt))
            C <- bond_fut[[k]]$conv_factor*bond_call_prices
            P <- bond_fut[[k]]$conv_factor*bond_put_prices
            KC <- bond_fut[[k]]$conv_factor*bond_call_strikes + allc
            KP <- bond_fut[[k]]$conv_factor*bond_put_strikes + allc
            T <- bond_fut[[k]]$option_term
            FWD <- bond_fut[[k]]$bond_fut_price
            fwd_term <- bond_fut[[k]]$res_term
            FWD_r <- exp(fwd_1*fwd_term)
            FWD_b <- (bond_fut[[k]]$conv_factor*FWD + allc)/FWD_r

            esp_a_mix <- function(x){
              esp_a_mix <- x[5]*exp(x[1] + (x[3]^2)/2 ) + (1 - x[5])*exp(x[2] + (x[4]^2)/2 )
            }

            var_mix <- function(x){
              var_mix <- x[5]*exp(2*x[1] + 2*x[3]^2) + (1 - x[5])*exp(2*x[2] + 2*x[4]^2) - esp_a_mix(x)^2
            }

            MSE_mix <- function(x){
              MSE_mix <- sum((C - model_prices(x)$model_call_price)^2, na.rm = T) +
                sum((P - model_prices(x)$model_put_price)^2, na.rm = T) +
                length(P)*(FWD_r - esp_a_mix(x[c(3, 4, 7, 8, 10)]))^2 +
                length(P)*(FWD_b - esp_a_mix(x[c(1, 2, 5, 6, 9)]))^2
              return(MSE_mix)}

            w_b <- bond_fut_rnd$params[5]
            volat_bond <- as.numeric(sqrt(w_b*bond_fut_rnd$params[3]^2 +
                                            (1 - w_b)*bond_fut_rnd$params[4]^2 +
                                            (1 - w_b)*w_b*(diff(bond_fut_rnd$params[1:2]) )^2 ))

            m1 <- m2 <- s1 <- s2 <- rho_1 <- rho_2 <- rho_3 <- rho_4 <- SCE <- NA
            m3 <- m4 <- log(FWD_r)
            s3 <- s4 <- sqrt(12/term_stir)*(unique(bond_fut[[k]]$res_term))^(3/2)/100*
              as.numeric(stir_fut_rnd$moments[2])
            pi1_r <- 0.5

            PARA <- as.matrix(data.frame(m1, m2, m3, m4, s1, s2, s3, s4, pi1_b = PR[, 1], pi1_r,
                                         rho_1, rho_2, rho_3, rho_4, w1 = 0, w2 = 0, SCE))

            if(FWD_b != 1){
              lower <- c( rep(c( (sign(1 - FWD_b)*0.2 + 1)*log(FWD_b), 0.5*volat_bond), each = nb_log),
                          rep(-1 + 1e-6, 2*nb_log), rep(1e-6, 2) )
              upper <- c( rep(c( (sign(FWD_b - 1)*0.2 + 1)*log(FWD_b), 2*volat_bond), each = nb_log),
                          rep(1 - 1e-6, 2*nb_log), rep(1 - 1e-6, 2) )
            } else {
              lower <- c( rep( c(-5*1e-4,  0.5*volat_bond), each = nb_log), rep(-1 + 1e-6, 2*nb_log), rep(1e-6, 2) )
              upper <- c( rep( c(1.5*1e-3, 2*volat_bond), each = nb_log), rep(1 - 1e-6, 2*nb_log), rep(1 - 1e-6, 2) )
            }

            objective <- function(x){
              if(cot_bond %in% c(1, 3)){
                MSE_mix( c(x[1:2], m3, m4, x[3:4], s3, s4, PR[i, 1], pi1_r, x[5:8] ))
              } else{MSE_mix( c(x[1:2], m3, m4, x[3:4], s3, s4, PR[i, 1], pi1_r, x[5:10] )) }
            }

            set.seed(123)

            suppressWarnings({
              for (i in 1:nrow(PR)){
                if(cot_bond %in%c(1,3)){
                  sol <- DEoptim(objective, lower[1:(length(lower) - 2)], upper = upper[1:(length(upper) - 2)],
                                 DEoptim.control(trace = FALSE, NP = 80, itermax = 100))
                  PARA[i, c("m1", "m2", "s1", "s2", "rho_1", "rho_2", "rho_3", "rho_4")] <-
                    sol$optim$bestmem
                } else{
                  sol <- DEoptim(objective, lower, upper = upper,
                                 DEoptim.control(trace = FALSE, NP = 120, itermax = 100))
                  PARA[i, c("m1", "m2", "s1", "s2", "rho_1", "rho_2", "rho_3", "rho_4",
                            "w1", "w2")] <- sol$optim$bestmem }
                PARA[i, "SCE"] <- sol$optim$bestval
              }
            })

            PARA <- PARA[ !is.na(PARA[, "m1"]), ]
            if(nrow(PARA) > 1){
              param <- PARA[which.min(PARA[, "SCE"]), -ncol(PARA)]
            } else {param <- as.matrix(PARA[, -ncol(PARA)], nrow = nrow(PARA))}

            param[param == 0] <- 1e-4

            L <- U <- rep(0, length(param))

            L[sign(param) == -1] <- 1.1*param[sign(param) == -1]
            L[sign(param) == 1] <- 0.9*param[sign(param) == 1]
            U[sign(param) == -1] <- 0.9*param[sign(param) == -1]
            U[sign(param) == 1] <- 1.1*param[sign(param) == 1]

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

            if(cot_bond%in%c(1,3) ){ params <- solu$par
            } else {params <- solu$par[1: (length(solu$par) - 2)]}

            marginal_bond[[k]] <- params[which(colnames(PARA)%in%c("m1", "m2","s1", "s2", "pi1_b"))]
            marginal_repo[[k]] <- params[which(colnames(PARA)%in%c("m3", "m4","s3", "s4", "pi1_r"))]

            params_each[[k]] <- cbind(marginal_repo[[k]], marginal_bond[[k]]) %>% t %>%
              data.frame %>% mutate_all(as.numeric)

            correl[[k]] <- (esp_mix(params) - esp_a_mix(marginal_repo[[k]])*esp_a_mix(marginal_bond[[k]]))/
              sqrt(var_mix(marginal_repo[[k]])*var_mix(marginal_bond[[k]]))

            corr_matrix[[k]] <- matrix(c(1, rep(correl[[k]], 2), 1), nrow = 2)

            n <- 10000
            joint[[k]] <- simulate_correlated_returns(n, params_each[[k]], corr_matrix[[k]]) %>%
              rename_with(~c("repo", "bond")) %>%
              mutate(fut_p = (repo*bond - allc)/bond_fut[[k]]$conv_factor) %>%
              rename_at(1, ~paste0(., k)) %>% arrange(bond)

          }

          bond_fut <- do.call(rbind, bond_fut)

          tri <- function(x){
            tri <- mapply(xirr, cf = mapply(c, -x, cf_other[[i]], cf_matu[[i]], SIMPLIFY = F),
                          tau = mapply(c, 0, mapply(unlist, cp_dt_2[[i]], SIMPLIFY = F), SIMPLIFY = F) )
          }

          dirty <- function(x){
            dcf <- mapply("/", list(c(unlist(cf_other[[j]]), cf_matu[[j]])),
                          mapply("^", 1 + x, list(unlist(cp_dt_2[[j]])), SIMPLIFY = F), SIMPLIFY = F)
            dirty <- unlist(lapply(dcf, sum))}

          ytm <- dy <- ctd <- px <- list()
          for (i in 1:length(joint)){
            ytm[[i]] <- tri(joint[[i]]$bond)
            dy[[i]] <- ytm[[i]] - bond_fut$ytm[i]
            ctd[[i]] <- px[[i]] <- list()
            for (j in 1:length(joint)){
              ctd[[i]][[j]] <- dy[[i]] + bond_fut$ytm[j]
              px[[i]][[j]] <- dirty(ctd[[i]][[j]])
            }
            px[[i]] <- do.call(cbind, px[[i]])
          }

          all_c <- bond_fut$acc_matu + bond_fut$coupon*bond_fut$Nomi/bond_fut$cp_freq*
            pmax(0, as.numeric(bond_fut$bond_fut_matu - bond_fut$curr_cp_dt))

          fut <- list()
          for (i in 1:length(joint)){
            fut[[i]] <- mapply("+", mapply("*", list(joint[[i]]$fut_p), bond_fut$conv_factor, SIMPLIFY = F),
                               all_c, SIMPLIFY = F)
            fut[[i]] <- do.call(cbind, fut[[i]])
          }

          rp <- mapply("/", fut, px, SIMPLIFY = F)

          prob <- list()
          for (i in 1:length(joint)){
            prob[[i]] <- apply(rp[[i]], 1, which.max)}

          proba_all <- do.call(cbind, prob)

          for (i in 1:ncol(proba_all)){
            proba_all[proba_all[, i] == i, -i] <- 0}

          prb <- rowSums(proba_all)

          doub <- doub_2 <- count <- list()
          for (i in 1:ncol(proba_all)){
            doub[[i]] <- rowSums(proba_all == i)
            doub_2[[i]] <- which(!doub[[i]]%in%c(0, 1))
          }

          if( length(sort(unique(unlist(doub_2)))) > 0){
            add_on <- matrix(proba_all[sort(unique(unlist(doub_2))) ,], ncol = nrow(deliverables))
            for (i in 1:nrow(add_on)){
              count[[i]] <- list()
              for (j in 1:length(joint)){
                count[[i]][[j]] <- length(which(add_on[i, ] == j))
              }
              count[[i]] <- unlist(count[[i]])
            }
            prb[sort(unique(unlist(doub_2)))] <- unlist(lapply(count, which.max))
          } else{ prb <- prb }

          proba_ctd <- list()
          for (k in 1:nrow(deliverables)){
            proba_ctd[[k]] <- length(which(prb == k))/length(prb)
          }

          proba_all <- data.frame(ISIN = deliverables$ISIN,
                                  proba = round(1000*unlist(proba_ctd))/1000) %>%
            arrange(desc(proba))

          return(proba_all)

        }
      } else {message("input dates are not consistent")}
    } else{ message ("please enter STIR options with maturity close to bond options' maturity (distance below 30 days)")}
  } else {message("inputs do not have the required length")}
}
