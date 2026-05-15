#' bond_future_charac_bbg
#'
#' @param bbg_tickers a vector of Bloomberg tickers of bond futures contracts, in character format
#'
#' @returns for each bond future contract, its maturity date in Date format, its currency, its listing place and its underlying asset in character format. Provided the contract has a Cheapest-to-Deliver Bond, its ISIN code in character format, its conversion factor in numeric format, its maturity date in Date format, its (annualized) coupon rate in numeric format, its frequency of coupon payment and its day count convention in character format
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

  underlying <- bdp(paste0(bbg_tickers, " Comdty"), "FUT_NOTL_BOND")

  bbg_tickers <- bbg_tickers[underlying$FUT_NOTL_BOND != ""]

  underlying <- underlying %>% filter(FUT_NOTL_BOND != "")

  if(nrow(underlying) == 0){message("please enter a bond future Bloomberg ticker")
  } else{

    flds <- c("FUT_HAS_CTD", "LAST_TRADEABLE_DT", "CRNCY", "FUT_EXCH_NAME_LONG", "FUT_NOTL_BOND",
              "FUT_CTD_ISIN", "FUT_CNVS_FACTOR", "FUT_CTD_MTY", "FUT_CTD_FREQ", "FUT_CTD_DAY_CNT")

    extract <- bdp(paste0(bbg_tickers, " Comdty"), flds) %>%
      rename_with(~c("ctd_exists", "fut_matu", "fut_currency", "fut_exchange", "fut_notional",
                     "ctd_isin", "ctd_conv_factor", "ctd_matu", "ctd_cp_fq", "ctd_day_count_conv")) %>%
      rownames_to_column("fut_bbg_ticker") %>% mutate_at("fut_bbg_ticker", ~gsub(" Comdty", "", .))
    no_ctd <- extract %>% filter(ctd_exists == F) %>%
      dplyr::select(c(fut_bbg_ticker, fut_matu, fut_currency, fut_exchange, fut_notional))
    extract <- extract %>% filter(ctd_exists == T) %>% dplyr::select(-ctd_exists)

    if(nrow(extract) > 0){
      extract <- extract %>% mutate_at("ctd_isin", ~paste0(., " Corp"))
      temp <- bdp(extract$ctd_isin, "CPN") %>% rename_with(~"ctd_cpn_rate") %>%
        rownames_to_column("ctd_isin") %>% mutate_at("ctd_cpn_rate", ~./100)
      extract <- extract %>% left_join(temp, by="ctd_isin") %>%
        relocate(ctd_cpn_rate, .before = ctd_cp_fq) %>% mutate_at("ctd_isin", ~gsub(" Corp", "", .))

      if( nrow(no_ctd) == 0){
        return(extract)
      } else{
        message(paste0("there is no CtD Bond for ",  no_ctd$fut_bbg_ticker))
        return(list(bond_future_wo_ctd_charac = no_ctd, bond_future_with_ctd_charac = extract))
        return(extract)}
    } else {
      message(paste0("there is no CtD Bond for ",  no_ctd$fut_bbg_ticker))
      return(no_ctd)}
  }
}
