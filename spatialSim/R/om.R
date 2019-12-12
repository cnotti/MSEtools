#' @export
om = function(a = NULL, data_om, seed = NULL, verbose = TRUE) {
  set.seed(seed)

  if (is.null(a) || a == 1) {
    t1 = Sys.time()
    # simulate species-space-time and species-space-time-size innovations    
    sim = TMB::MakeADFun(data = data_om[!names(data_om) %in% "extra"],
                         parameters = data_om["dummyParameter"], 
                         DLL = "SIM")$simulate()
    data_om["I_b"] = sim["I_b"]
    data_om["I_t"] = sim["I_t"]
    data_om["H_sbl"] = sim["H_sbl"]
    data_om["H_stl"] = sim["H_stl"]
    data_om["E_csb"] = sim["E_csb"]
    data_om["E_cst"] = sim["E_cst"]
    data_om["M_csbl"] = sim["M_csbl"]
    data_om["M_cstl"] = sim["M_cstl"]
    allSim = c(data_om[["M_csbl"]], data_om[["M_cstl"]], 
               data_om[["E_csb"]], data_om[["E_cst"]])
    if (any(is.na(allSim)) | any(is.null(allSim)) | any(is.infinite(allSim))){
      stop("simulated random processes are not all real numbers")
    }
    run_time = Sys.time() - t1
    if (verbose == TRUE) {
      cat(paste("simulation of random processes completed in", 
                round(run_time, 3), attr(run_time, "units"),
                "\n"))
    }
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
  if (verbose == TRUE) {
    cat(paste("harvest completed in", 
              round(run_time, 3), attr(run_time, "units"),
              "\n"))
  }
  
  return(list(data_om = data_om,
              om_rep  = om_rep,
              harvest_time = run_time))
}
