#' sc_repo_spread
#'
#' @param bond_call_prices_1 a vector of call prices on the first bond futures, in numeric format
#' @param bond_call_strikes_1 a vector of call strikes attached to the call prices, in numeric format
#' @param bond_put_prices_1 a vector of put prices on the first bond futures, in numeric format
#' @param bond_put_strikes_1 a vector of put strikes attached to the put prices, in numeric format
#' @param bond_call_prices_2 a vector of call prices on the second bond futures, in numeric format
#' @param bond_call_strikes_2 a vector of call strikes attached to the call prices, in numeric format
#' @param bond_put_prices_2 a vector of put prices on the second bond futures, in numeric format
#' @param bond_put_strikes_2 put_strikes_2 a vector of put strikes attached to the put prices, in numeric format
#' @param stir_call_prices a vector of call prices on a STIR futures whose maturity date is near the maturity date of the bond futures (maturity distance below 30 days), in numeric format
#' @param stir_call_strikes a vector of call strikes attached to the call prices, in numeric format
#' @param stir_put_prices a vector of put prices on the same STIR futures, in numeric format
#' @param stir_put_strikes a vector of put strikes attached to the put prices, in numeric format
#' @param r a number for the riskfree spot rate whose maturity is equal to the options' maturity, in numeric format
#' @param r_2 a number for the spot repo funding rate of the underlying bond, with maturity equal to the options' maturity (by default, the General Collateral repo rate of the jurisdiction), in numeric format
#' @param r_3 a number for the spot repo funding rate of the underlying bond, with maturity equal to the futures' maturity (by default, the General Collateral repo rate of the jurisdiction), in numeric format
#' @param day_count_conv a number for the day count convention of the bond futures, 1 for ACT/ACT, 2 for ACT/360, 3 for ACT/365 and 4 for 30/360, in numeric format
#' @param cot_bond_1 a number for style of the options on the first bond futures, 1 for European options, 2 for American options and 3 for American options with futures-style margin, in numeric format
#' @param cot_bond_2 a number for style of the options on the second bond futures, 1 for European options, 2 for American options and 3 for American options with futures-style margin, in numeric format
#' @param bond_fut_price_1 a number for the first bond futures' price on calibration date, in numeric format
#' @param bond_fut_price_2 a number for the second bond futures' price on calibration date, in numeric format
#' @param bond_cp_1 a number for the coupon rate of the current Cheapest-to-Deliver bond of the first futures contract, in numeric format
#' @param bond_cp_f_1 a number for the frequency of coupon payment of the current Cheapest-to-Deliver bond in the first futures contract, either 1 for annual payment or 2 for semiannual payment, in numeric format
#' @param bond_conv_factor_1 a number for the conversion factor of the current Cheapest-to-Deliver bond of the first futures contract, in numeric format
#' @param bond_matu_1 a date for the maturity date of the current Cheapest-to-Deliver Bond of the first bond futures, in Date format
#' @param bond_cp_2 a number for the coupon rate of the current Cheapest-to-Deliver bond of the second futures contract, in numeric format
#' @param bond_cp_f_2 a number for the frequency of coupon payment of the current Cheapest-to-Deliver bond in the second futures contract, either 1 for annual payment or 2 for semiannual payment, in numeric format
#' @param bond_conv_factor_2 a number for the conversion factor of the current Cheapest-to-Deliver bond of the second futures contract, in numeric format
#' @param bond_matu_2 a date for the maturity date of the current Cheapest-to-Deliver Bond of the second bond futures, in Date format
#' @param sett a number for the number of days between the ex-coupon date and the coupon payment date of the bonds, in numeric format
#' @param bond_N a single number for the nominal of the bonds (100 by default) in numeric format
#' @param bond_fut_matu a date for the bond futures' maturity date, in Date format
#' @param bond_option_matu a date for the maturity date of the options, in Date format
#' @param cot_stir a number for the STIR options' style, 1 for European options, 2 for American options and 3 for American options with futures-style margin, in numeric format
#' @param stir_fut_price a number for the STIR futures' price on calibration date, in numeric format
#' @param stir_fut_matu a date for the maturity date of the STIR futures contract, in Date format
#' @param stir_option_matu a date for the maturity date of the STIR options, in Date format
#' @param start_date a date for the observation date, in Date format
#' @param country_1 a character for the nationality of the issuer of the bond in the first bond futures, in character format
#' @param country_2 a character for the nationality of the issuer of the bond in the second bond futures, in character format
#' @param term_stir a number for the term to maturity, in months, of the STIR at observation date, in numeric format
#' @param corr_rp a number for the correlation coefficient between the two bond repo rates, in numeric format
#'
#' @returns based on the user's correlation coefficient, a density plot of the spread, in bps, of 10,000 realizations of the repo rate spread between two cheapest-to-deliver bonds of either different issuers or of the same issuer with different maturities, using options on two bond futures of identical maturity and options of a STIR contract with a maturity close to the maturity of the bond futures options
#' @export
#'
#' @importFrom stats approx constrOptim density dlnorm nlminb plnorm pnorm
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
#' sc_repo_spread(
#' c(20.71, 20.21, 19.71, 19.21, 18.71, 18.21, 17.71,
#' 17.21, 16.71, 16.22, 15.72, 15.22, 14.72, 14.22, 13.72,
#' 13.22, 12.73, 12.23, 11.73, 11.24, 10.74, 10.24, 9.75,
#' 9.25, 8.76, 8.27, 7.78, 7.29, 6.80, 6.32, 5.84, 5.36,
#' 4.89, 4.42, 3.97, 3.53, 3.10, 2.69, 2.30, 1.93, 1.59,
#' 1.28, 1.01, 0.77, 0.58, 0.43, 0.31, 0.23, 0.17, 0.13,
#' 0.10, 0.08, 0.06, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02,
#' 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
#' seq(105, 138.5, 0.5),
#' c(0.01,  0.01,  0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
#' 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.03, 0.03, 0.03,
#' 0.04, 0.04, 0.04, 0.05, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
#' 0.12, 0.14, 0.16, 0.19, 0.22, 0.27, 0.33, 0.40, 0.49,
#' 0.60, 0.73, 0.89, 1.08, 1.31, 1.57, 1.88, 2.23, 2.61, 3.03,
#' 3.47, 3.93, 4.40, 4.88, 5.36, 5.85, 6.34, 6.83, 7.33, 7.82,
#' 8.32, 8.82, 9.32, 9.81, 10.31, 10.81, 11.31, 11.81, 12.31, 12.81),
#' seq(105, 138.5, 0.5),
#' c(2.205, 2.105, 2.010, 1.910, 1.810, 1.715, 1.615, 1.520,
#' 1.425, 1.330, 1.235, 1.140, 1.045, 0.955, 0.865, 0.780,
#' 0.695, 0.615, 0.540, 0.470, 0.400, 0.335, 0.280, 0.225,
#' 0.180, 0.140, 0.110, 0.085, 0.065, 0.050, 0.040, 0.030, 0.025,
#' 0.020, 0.015, 0.015, 0.015, 0.010, 0.010, 0.010, 0.010,
#' 0.010, 0.010, 0.005, 0.005, 0.005, 0.005, 0.005),
#' c(seq(103.6, 107.9, 0.1), 108.3, 108.5, 108.6, 108.7),
#' c(0.010, 0.010, 0.015, 0.015, 0.015, 0.020, 0.020, 0.025,
#' 0.030, 0.035, 0.040, 0.045, 0.050, 0.060, 0.070, 0.085,
#' 0.100, 0.120, 0.145, 0.175, 0.205, 0.240, 0.285, 0.330,
#' 0.385, 0.445, 0.515, 0.590, 0.670, 0.755, 0.845, 0.935,
#' 1.030, 1.125, 1.220, 1.320, 1.420, 1.515, 1.615, 1.715,
#' 1.815, 1.915, 2.015, 2.110, 2.510, 2.710, 2.810, 2.910),
#' c(seq(103.6, 107.9, 0.1), 108.3, 108.5, 108.6, 108.7),
#' c(2.5300, 2.4050, 2.2800, 2.1550, 2.0300, 1.9050, 1.7800,
#' 1.6550, 1.5300, 1.4700, 1.4075, 1.3450, 1.2825, 1.1575,
#' 1.0350, 0.9125, 0.8525, 0.7900, 0.7300, 0.6700,
#' 0.6125, 0.5525, 0.4950, 0.4400, 0.3850, 0.3325, 0.2800,
#' 0.2325, 0.1875, 0.1475, 0.1100, 0.0775, 0.0550, 0.0375,
#' 0.0250, 0.0150, 0.0125, 0.0100, 0.0075, 0.0075, 0.0050,
#' 0.0050, 0.0050, 0.0050, 0.0050, 0.0050, 0.0025, 0.0025,
#' 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025,
#' 0.0025, 0.0025, 0.0025, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000),
#' c(seq(95.0000, 96.000, 0.125), 96.0625, 96.1250, 96.1875 ,
#' 96.2500, 96.3750, 96.5000,
#' seq(96.6250, 99.2500, 0.0625),
#' 100.3750, 100.6250, 100.7500, 100.8750, 101.0000),
#' c(0.000, 0.0000, 0.0000, 0.0000, 0.000, 0.0000, 0.0000,
#' 0.0000, 0.000, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025,
#' 0.0050, 0.0075, 0.0100, 0.0100, 0.0125, 0.0150, 0.0200,
#' 0.0225, 0.0275, 0.0350, 0.0425, 0.0525, 0.0625, 0.0775,
#' 0.0950, 0.1175, 0.1425, 0.1725, 0.2125, 0.2575, 0.3075,
#' 0.3600, 0.4200, 0.4800, 0.5400, 0.6025, 0.6625, 0.7250,
#' 0.7875, 0.8500, 0.9125, 0.9750, 1.0350, 1.0975, 1.1600,
#' 1.2225, 1.2850, 1.3475, 1.4100, 1.4725, 1.5350, 1.5975,
#' 1.6600, 1.7225, 2.8450, 3.0950, 3.2200, 3.3450, 3.4700),
#' c(seq(95.0000, 96.000, 0.125), 96.0625, 96.1250,
#' 96.1875, 96.2500, 96.3750, 96.5000,
#' seq(96.6250, 99.2500, 0.0625), 100.3750,
#' 100.6250, 100.7500, 100.8750, 101.0000),
#' 0.02052,
#' 0.02052,
#' 0.02059,
#' 1,
#' 3,
#' 3,
#' 125.7,
#' 105.795,
#' 0.025,
#' 1,
#' 0.768181,
#' as.Date("2035-02-15"),
#' 0.021,
#' 1,
#' 0.936365,
#' as.Date("2028-03-15"),
#' 2,
#' bond_N = 100,
#' as.Date("2026-06-08"),
#' as.Date("2026-05-22"),
#' 3,
#' 97.53,
#' as.Date("2026-06-15"),
#' as.Date("2026-06-15"),
#' as.Date("2026-04-01"),
#' "Germany",
#' "Germany",
#' 3,
#' 0.9)
#' }
#

sc_repo_spread <- function(bond_call_prices_1, bond_call_strikes_1, bond_put_prices_1, bond_put_strikes_1,
                           bond_call_prices_2, bond_call_strikes_2, bond_put_prices_2, bond_put_strikes_2,
                           stir_call_prices, stir_call_strikes, stir_put_prices, stir_put_strikes,
                           r, r_2, r_3, day_count_conv, cot_bond_1, cot_bond_2, bond_fut_price_1, bond_fut_price_2,
                           bond_cp_1, bond_cp_f_1, bond_conv_factor_1, bond_matu_1,
                           bond_cp_2, bond_cp_f_2, bond_conv_factor_2, bond_matu_2, sett, bond_N = 100,
                           bond_fut_matu, bond_option_matu, cot_stir, stir_fut_price,
                           stir_fut_matu, stir_option_matu, start_date, country_1,
                           country_2, term_stir, corr_rp){

  if(length(corr_rp) == 1){

    bond_1 <- bond_fut_irr_ytm(bond_call_prices_1, bond_call_strikes_1, bond_put_prices_1,
                               bond_put_strikes_1, stir_call_prices, stir_call_strikes, stir_put_prices, stir_put_strikes, r, r_2, r_3,
                               day_count_conv, cot_bond_1, bond_fut_price_1, bond_cp_1, bond_cp_f_1, bond_conv_factor_1, sett, bond_N = 100,
                               bond_matu_1, bond_fut_matu, bond_option_matu, cot_stir, stir_fut_price, stir_fut_matu, stir_option_matu,
                               start_date, country = country_1, currency = NA, term_stir)

    bond_2 <- bond_fut_irr_ytm(bond_call_prices_2, bond_call_strikes_2, bond_put_prices_2,
                               bond_put_strikes_2, stir_call_prices, stir_call_strikes, stir_put_prices, stir_put_strikes, r, r_2, r_3,
                               day_count_conv, cot_bond_2, bond_fut_price_2, bond_cp_2, bond_cp_f_2, bond_conv_factor_2, sett, bond_N = 100,
                               bond_matu_2, bond_fut_matu, bond_option_matu, cot_stir, stir_fut_price, stir_fut_matu, stir_option_matu,
                               start_date, country = country_2, currency = NA, term_stir)

    bond_fut <- data.frame(bond_option_matu, start_date, sett, bond_fut_matu, currency = NA, bond_N) %>%
      rename_with(~c("option_matu", "start_date", "sett", "fut_matu", "country", "Nomi")) %>%
      mutate_at(c("option_matu", "fut_matu", "start_date"), as.Date)

    if(length(bond_1)  == 18 & length(bond_2)  == 18 & length(bond_fut) > 0){

      if(day_count_conv == 1){
        bond_fut <- bond_fut %>% mutate(res_term = as.numeric(fut_matu - option_matu)/
                                          as.numeric(ceiling_date(fut_matu, "year") - floor_date(option_matu, "year") ))
      } else if(day_count_conv == 2){
        bond_fut <- bond_fut %>% mutate(res_term = as.numeric(fut_matu - option_matu)/360)
      } else if(day_count_conv == 3){
        bond_fut <- bond_fut %>% mutate(res_term = as.numeric(fut_matu - option_matu)/365)
      } else {
        bond_fut <- bond_fut %>% mutate(stub_1_res = max(0, 30 - as.numeric(format(option_matu + sett, "%d"))),
                                        stub_2_res = min(30, as.numeric(format(fut_matu, "%d"))),
                                        plain_months_res = round(as.numeric(floor_date(fut_matu, "months") -
                                                                              ceiling_date(option_matu + sett, "months") )/30),
                                        res_term = (stub_1_res + stub_2_res + max(0, plain_months_res)*30)/360)}

      fwd_term <- bond_fut$res_term

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

      params_repo <- rbind(bond_1$params_repo.unlist.marginal_repo.,
                           bond_2$params_repo.unlist.marginal_repo.) %>%
        data.frame

      corr_matrix <- matrix(c(1, rep(corr_rp, 2), 1), nrow = 2)

      n <- 10000
      set.seed(123)
      joint <- simulate_correlated_returns(n, params_repo, corr_matrix) %>%
        rename_with(~c("repo_1", "repo_2")) %>% mutate_all(~log(.)/fwd_term) %>%
        mutate(spread = 10000*(repo_2 - repo_1))

      spread_plot <- ggplot(joint, aes(x = spread)) + geom_density() +
        labs(title = paste0("spread (bp) between ", country_2, " ", 100*bond_cp_2, "% ", bond_matu_2, " repo rate and ",
                            country_1, " ", 100*bond_cp_1, "% ", bond_matu_1, " repo rate as of ", start_date),
             subtitle = "Risk Neutral Probability Density",
             x = "spread", y = "probability density") +
        theme(legend.position = "bottom", legend.title = element_blank(), plot.margin = margin(.8,.5,.8,.5, "cm")) +
        theme_bw()

      return(spread_plot = spread_plot)

    } else {message("impossible to retrieve a density for the spread between the bond repo rates") }
  } else {message("inputs do not have the required length")}
}
