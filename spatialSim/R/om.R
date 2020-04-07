#' Run operating model
#'
#' \code{om} runs the \code{operating model}
#'
#' @param a \code{NULL} or the number of time periods to run OM before updating management plan
#' @param data_om tagged list containing model parameters
#' @param seed seed for random number generation
#' @param verbose if \code{TRUE} then runtime feedback is provided
#' 
#' @return object of class \code{om}, containing the om output:
#' 
#' \describe{
#'    \itemize{
#'      \item \code{data_om}: updated list of model parameters
#'      \itemize{
#'        \item ...
#'      }
#'      \item \code{om_rep}: reported values from \code{C++} file
#'      \itemize{
#'        \item \code{SSB0_c}: unfished spawning-stock-biomass for each species
#'        \item \code{N_csl}: numbers of individuals-at-size from last simulated time period. This variable is used for model re-entry when \code{a != NULL}
#'        \item \code{B_cst}: biomass by species, time and spatial location
#'        \item \code{ncatch_cstl}: number of harvested individuals-at-size attributed to each mesh loation
#'        \item \code{c_i}: zero-indexed species associated with the ith harvest tow
#'        \item \code{f_i}: zero-indexed fishing grid cell associated with the ith harvest tow
#'        \item \code{t_i}: zero-indexed time period associated with the ith harvest tow
#'        \item \code{catch_i}: harvested biomass from ith harvest tow
#'        \item \code{ncatch_il}: harvested numbers-at-size from ith harvest tow
#'      }
#'      \item \code{harvest_time}: information on model execution time
#'      \itemize{
#'        \item ...
#'      }
#'    }
#' }
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
  harvest_time = paste("harvest completed in", 
                       round(run_time, 3), attr(run_time, "units"),
                       "\n")
  if (verbose == TRUE) {
    cat(harvest_time)
  }
  
  om_out = list(data_om = data_om,
                om_rep  = om_rep,
                harvest_time = harvest_time)
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
print.om <- function(x, units = "kg", ...) {
  Brange = round(apply(
    apply(x$om_rep[["B_cst"]], c(1,3), "sum"), 1, "range")
    )
  
  Brange.d = apply(
    with(x$data_om, x$om_rep[["B_cst"]]/(areaOmegastar/ns)), 1, "range"
  )
  
  cat("==================== OM summary ====================\n")
  cat(paste("Output for period up to t =", 
            x$data_om[["tstop"]], 
            "\n\n"))
  
  for (c in 1:x$data_om$nc) {
    cat(paste("Output for species c =", 
              x$data_om$extra$species_names[c], ":\n"))
    cat(paste("\t* total biomass range", 
               paste(Brange[,c], collapse = " -- "), units, "\n"))
    cat(paste("\t* density range per unit area sq:", 
              paste(round(Brange.d[,c], 2), collapse = " -- "), units, "\n\n"))
  }
  
  cat(x$harvest_time)
  cat("====================================================\n")
}

