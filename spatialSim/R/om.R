#' Run operating model
#'
#' \code{om} runs the \code{operating model}
#'
#' @param a \code{NULL} or the number of time periods to run OM before updating management plan
#' @param data_om tagged list containing model parameters
#' @param seed seed for random number generation
#' @param verbose if \code{TRUE} then runtime feedback is provided
#' 
#' @return object of class \code{om}, containing om output:
#' 
#' \describe{
#'   \item{data_om}{updated list of model parameters}
#'   \item{om_rep}{reported values from \code{C++} file}
#'   \itemize{
#'     \item{SSB0_c}{unfished spawning-stock-biomass for each species}
#'     \item{N_csl}{numbers of individuals-at-size from last simulated time period. This variable is used for model re-entry when \code{a != NULL}.}
#'     \item{B_cst}{}
#'     \item{ncatch_cstl}{number of harvested individuals-at-size attributed to each mesh loation}
#'     \item{c_i}{species associated with the ith harvest tow}
#'     \item{f_i}{fishing grid cell associated with the ith harvest tow}
#'     \item{t_i}{time period associated with the ith harvest tow}
#'     \item{catch_i}{harvested biomass from ith harvest tow}
#'     \item{ncatch_il}{harvested numbers-at-size from ith harvest tow}
#'   }
#'   \item{harvest_time}{information on model execution time}
#' }
#' 
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
  
  om_out = list(data_om = data_om,
                om_rep  = om_rep,
                harvest_time = run_time)
  class(om_out) = "om"
  return(om_out)
}


#' Print output from \code{\link{om}}
#'
#' @title Print data
#' @param x Output from \code{\link{om}}
#' @param ... Not used
#' @return NULL
#' @method print om
#' @export
print.om <- function(x, ...) {
  cat(paste("om output for period up to t=", 
            data_om[["tstop"]], 
            "\n"))
  cat(paste("harvest for previous period completed in", 
            round(run_time, 3), attr(run_time, "units"),
            "\n"))
}