#' bond_future_charac_bbg
#'
#' @param bbg_tickers a vector of Bloomberg tickers of bond futures contracts, in character format
#'
#' @returns for each bond futures contract, its currency, its listing place and its underlying asset all in character format, its first delivery date and its maturity date, all in Date format
#' @export
#'
#' @import dplyr
#' @import Rblpapi
#' @import tibble
#' @examples
#' \dontrun{
#' bond_future_charac_bbg(c("TYU24", "CNM6", "KAAH6", "IKM6", "OATU6"))
#' }
#'
bond_future_charac_bbg <- function(bbg_tickers){

  blpConnect()

  underlying <- bdp(paste0(bbg_tickers, " Comdty"), "FUT_NOTL_BOND")

  bbg_tickers <- bbg_tickers[underlying$FUT_NOTL_BOND != ""]

  underlying <- underlying %>% filter(FUT_NOTL_BOND != "")

  if(nrow(underlying) == 0){message("please enter a bond futures Bloomberg ticker")
  } else{

    flds <- c("CRNCY", "FUT_EXCH_NAME_LONG", "FUT_NOTL_BOND", "FUT_FIRST_TRADE_DT",
              "FUT_DLV_DT_FIRST", "LAST_TRADEABLE_DT")

    extract <- bdp(paste0(bbg_tickers, " Comdty"), flds) %>%
      rename_with(~c("fut_currency", "fut_exchange", "fut_notional", "first_trade_dt",
                     "first_delivery_dt", "fut_matu")) %>%
      rownames_to_column("fut_bbg_ticker") %>% mutate_at("fut_bbg_ticker", ~gsub(" Comdty", "", .))
    return(extract)
  }
}
