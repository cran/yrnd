#' bond_future_charac_bbg
#'
#' @param bbg_tickers a vector of Bloomberg tickers, in character format
#'
#' @returns for each bond future contract, provided it has a Cheapest-to-Deliver Bond, the maturity date of the futures contract, in Date format, its currency in character format, the conversion factor assigned by the exchange to the Cheapest-to-Deliver Bond of the futures contract in numeric format, the maturity date of the CtD Bond in Date format, the coupon rate of the CtD Bond in numeric format, the frequency of coupon payment of the CtD Bond in character format, the day count convention of the CtD Bond in character format
#' @export
#'
#' @import dplyr
#' @import Rblpapi
#' @import tibble
#' @examples
#' \dontrun{
#' bond_future_charac_bbg(c("TYU24", "CNM6", "KAAH6", "IKM6", "OATU6"))
#' }

bond_future_charac_bbg <- function(bbg_tickers){
  blpConnect()
  flds <- c("FUT_HAS_CTD", "LAST_TRADEABLE_DT", "CRNCY", "FUT_CNVS_FACTOR", "FUT_CTD_ISIN",
            "FUT_CTD_MTY", "FUT_CTD_FREQ", "FUT_CTD_DAY_CNT")

  extract <- bdp(paste0(bbg_tickers, " Comdty"), flds) %>%
    rename_with(~c("ctd_exists", "fut_matu", "fut_currency", "ctd_conv_factor", "ctd_isin",
                   "ctd_matu", "ctd_cp_fq", "day_count_conv")) %>% rownames_to_column("fut_bbg_ticker") %>%
    mutate_at("fut_bbg_ticker", ~gsub(" Comdty", "", .))
  no_ctd <- extract %>% filter(ctd_exists == F) %>% dplyr::select(fut_bbg_ticker)
  extract <- extract %>% filter(ctd_exists == T)

  if(nrow(extract) > 0){
    extract <- extract %>% mutate_at("ctd_isin", ~paste0(., " Corp"))
    temp <- bdp(extract$ctd_isin, "CPN") %>% rename_with(~"ctd_cpn_rate") %>%
      rownames_to_column("ctd_isin") %>% mutate_at("ctd_cpn_rate", ~./100)
    extract <- extract %>% left_join(temp, by="ctd_isin") %>% relocate(ctd_cpn_rate, .before = ctd_cp_fq) %>%
      dplyr::select(-c(ctd_exists, ctd_isin))

    if( length(no_ctd) == 0){
      return(extract)
    } else{
      message(paste0("there is no CtD Bond for ",  no_ctd))
      return(extract)}
  } else {message(paste0("there is no CtD Bond for ",  no_ctd))}
}
