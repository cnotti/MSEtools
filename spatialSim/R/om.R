#' @export
om = function(a = NULL, data_om, seed = 105) {
  set.seed(seed)

  if (is.null(a) || a == 1) {
    # simulate species-space-time and species-space-time-size innovations    
    sim = TMB::MakeADFun(data = data_om[!names(data_om) %in% "extra"],
                         parameters = data_om["dummyParameter"], 
                         DLL = "SIM")$simulate()
    
    data_om["E_csb"] = sim["E_csb"]
    data_om["E_cst"] = sim["E_cst"]
    data_om["M_csbl"] = sim["M_csbl"]
    data_om["M_cstl"] = sim["M_cstl"]
  }
  
  if (!is.null(a)) {
    data_om[["tstop"]] = data_om[["extra"]][["tstop_a"]][a]
    data_om[["tstart"]] = data_om[["extra"]][["tstart_a"]][a]
  }
  
  t1 = Sys.time()
  om = TMB::MakeADFun(data = data_om[!names(data_om) %in% "extra"],
                      parameters = data_om["dummyParameter"], 
                      type = "Fun", DLL = "OM")
  om_rep = om$report(unlist(data_om["dummyParameter"]))
  run_time = Sys.time() - t1
  cat(paste("harvest completed in", 
            round(run_time, 3), attr(run_time, "units"),
            "\n"))
  
  return(list(data_om = data_om,
              om_rep  = om_rep))
}
