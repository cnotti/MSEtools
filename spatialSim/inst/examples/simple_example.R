library("spatialSim")

seed = 50
set.seed(seed)

n_cases = 20L
n_mult = 4L
om_out = vector("list", n_cases)
harvest_times = vector("list", n_cases)

for (j in 1:5) {
  for (i in 1:n_cases) {
    # number of species
    ny = 5
    np = 12
    nc = 2
    nl = 10
    fishp1 = 11
    
    # define spatial objects
    areadredge = 100
    xx = c(174.0, 174.025)
    multiplier = seq(1, n_mult, len = n_cases)[i]
    yy = c(-41.52, -41.52 + 0.24*multiplier) 
  
    loc_s = expand.grid(seq(xx[1], xx[2], len = 10), 
                        seq(yy[1], yy[2], len = 50 * multiplier))
    ns = nrow(loc_s)
    mesh_om = inla.mesh.2d(loc_s, max.n.strict = 3000, 
                           offset = c(0.005, 0.025))
    
    coords = cbind(c(xx[1], xx[2], xx[2], xx[1], xx[1]),
                   c(yy[1], yy[1], yy[2], yy[2], yy[1]))
    domain = SpatialPolygons( list( Polygons(list(Polygon(coords)), "Domain") ) )
    areaOmegastar = geosphere::areaPolygon(domain)
    
    loc_f = coordinates(spsample(domain, 
                                 floor(areaOmegastar/areadredge), 
                                 type = "regular"))
    nf = nrow(loc_f)
    areadredge = areaOmegastar/nf
    
    # calculate psi_cs from bathymetry
    depths = -c(1:7)
    nbounds = length(depths)
    dbounds = c(xx[1], xx[1] + cumsum(rep((xx[2] - xx[1])/nbounds, nbounds)))
    depth_s = rep(0, ns)
    depth_f = rep(0, nf)
    for(d in 1:nbounds) {
      depth_s[dbounds[d] <= loc_s[,1] & dbounds[d+1] >= loc_s[,1]] = depths[d]
      depth_f[dbounds[d] <= loc_f[,1] & dbounds[d+1] >= loc_f[,1]] = depths[d]
    }
    
    lnmud_c = c(1.317284, 1.468462)
    sigmad_c = c(0.3074893, 0.2669919)
    psi_cs = matrix(0, nc, ns)
    for (c in 1:nc) {
      psi_cs[c,] = plnorm(-depth_s, lnmud_c[c], sigmad_c[c]) - 
        plnorm(-depth_s-1, lnmud_c[c], sigmad_c[c])
    }
    
    # setup species specific fishing domain book-keeping vars 
    f_c_fstar = vector("list", nc)
    dbounds_c2 = matrix(-sapply(c(0.1,0.9), 
                                function(x) qlnorm(x, lnmud_c, sigmad_c)), 
                        nrow = nc, ncol = 2)
    for (c in 1:nc) {
      # note depths are stored as negative numbers (i.e., deaper < shallower)
      f_c_fstar[[c]] = which(depth_f < dbounds_c2[c,1] & 
                               depth_f > dbounds_c2[c,2]) - 1
    }
    
    # set sampling sites 
    f_ct_fsurv = as.vector(rep(-1, nc*ny*np), "list")
    ct = 1
    p = 1
    for(c in 1:nc) {
      for(t in 1:(ny*np)){
        if (p == fishp1+1) {
          f_ct_fsurv[[ct]] = sample(f_c_fstar[[c]], 50, replace = FALSE)
        }
        p = ifelse(p < np, p + 1, 1)
        ct = ct + 1
      }
    }
    
    # define growth, selectivity & maturity functions
    linf_c = c(mean(c(60.8, 60.3, 57.6, 52.1, 54.6)),
               mean(c(75.2, 88, 80.6, 72.3, 72.4)))
    beta_c = c(mean(c(0.48, 1.01, 1.74, 0.8, 1.44)),
               mean(c(0.35, 0.57, 0.58, 0.6, 1.84)))/np
    fn_growth = function(lmid_cl) {
      delta_cl = matrix(0, nrow = nrow(lmid_cl), ncol = ncol(lmid_cl))
      for (c in 1:nrow(lmid_cl)){
        delta_cl[c,] = (1 - exp(-beta_c[c])) * (linf_c[c] - lmid_cl[c,])
      }
      delta_cl
    }
    fn_mature = function(lmid_cl) {
      l50_c = 0.25 * linf_c
      l95_c = 0.5 * linf_c
      pmat_cl = matrix(0, nrow = nrow(lmid_cl), ncol = ncol(lmid_cl))
      for (c in 1:nrow(lmid_cl)){
        pmat_cl[c,] = 1/(1+exp(-log(19)*((lmid_cl[c,]-l50_c[c])/(l95_c[c]-l50_c[c]))))
      }
      pmat_cl
    }
    fn_select = function(lmid_cl) {
      l50_c = 0.025 * linf_c
      l95_c = 0.05 * linf_c
      select_cl = matrix(0, nrow = nrow(lmid_cl), ncol = ncol(lmid_cl))
      for (c in 1:nrow(lmid_cl)) {
        select_cl[c,] = 1/(1+exp(-log(19)*((lmid_cl[c,]-l50_c[c])/(l95_c[c]-l50_c[c]))))
      }
      select_cl
    }
    fn_weight = function(lmid_cl, beta_ci) {
      weight_cl = matrix(0, nrow = nrow(lmid_cl), ncol = ncol(lmid_cl))
      beta_ci = cbind("beta0" = c(-8.585949,
                                  -7.5989853),
                      "beta1" = c(0.167854,
                                  0.1181136),
                      "beta2" = c(-0.001375,
                                  -0.0006332))
                                  
      for (c in 1:nc) {
        weight_cl[c,] = exp(beta_ci[c,1] + 
                              beta_ci[c,2]*lmid_cl[c,] + 
                              beta_ci[c,3]*lmid_cl[c,]^2)
      }
      weight_cl
    }
    
    # create input data list
    data_om = setup_om(
      # indexing/population structure
      nb = 600, 
      ny = ny, 
      np = np, 
      nc = nc,
      nx = 50,
      nl = nl, 
      nfishp = 1, 
      fishp1 = fishp1,
      
      # growth 
      fn_growth = fn_growth,
      fn_weight = fn_weight,
      sigmaG_c = c(0.15, 0.15),
      
      # natural mortality
      alpha_c = c(0.79, 0.42),
      mugI = -2,
      mugH = 0,
      sigmagI = 1,
      taugH = 0.1,
      kappaH = exp(3),
      phiIt = 0.75,
      phiHt = 0.75,
      phiHl = 0.75,
      
      # fishing mortality
      fn_select = fn_select,
      limit_c = c(1000, 1000)*1e3,
      areadredge = areadredge,
      ptarget = 0.99,
      F_intensity = 2.5,
      F_settings = 0,
      f_c_fstar = f_c_fstar,
      f_ct_fsurv = f_ct_fsurv,
      
      # recruitment
      fn_mature = fn_mature,
      R0_c = exp(c(5, 3)),
      h_c = c(0.5, 0.5),
      psi_p = c(rep(0, np/2), rep(1/(np/2), np/2)),
      psi_l = matrix(c(1, rep(0, nl - 1)), ncol = 1), 
      psi_cs = psi_cs,
      mugE = -5,
      taugE = 0.00025, #0.00025
      kappaE = exp(4.5),
      phiEt = 0.75,
      Rrange_c = rep(0.05, nc),
      
      # spatial objects
      projection = "ll",
      loc_s = loc_s,
      loc_f = loc_f,
      mesh = mesh_om,
      lnaniso = rep(0, 4),
      
      # population size structure
      lmin_c = rep(0, nc), 
      lmax_c = linf_c,
      
      # extra
      seed = seed,
      species_names = c("SAE", "MMI"))
    
    
    # run model
    cat(paste("running case where domain area =", round(areaOmegastar/1e6, 2), "km sq\n"))
    om_out[[i]] = om(data = data_om, seed = seed)
    harvest_times[[i]][j] = om_out[[i]][["harvest_time"]]
  }
}

if (FALSE) {
  # decorrelation distance
  geosphere::distm(c(174.025, -41.52), c(174.025 + sqrt(8)/exp(3), -41.52))
  geosphere::distm(c(174.025, -41.52), c(174.025 + sqrt(8)/exp(4.5), -41.52))
  
  # biomass map
  data_om = om_out[[1]]$data_om
  om_rep = om_out[[1]]$om_rep
  B_cst = with(data_om, array(0, dim = c(nc, mesh_om$n, nt)))
  y_t = NULL
  for (c in 1:data_om[["nc"]]) {
    for (t in 1:data_om[["nt"]]) {
      y_t[t] = floor((t - 1)/12) + 1
      B_cst[c,data_om[["splus_s"]]+1,t] = om_rep[["B_cst"]][c,,t]/areaf*data_om$ns
    }
  }
  proj_om = inla.mesh.projector(mesh_om, dims = c(200,200))
  pdf("C:/Users/Chris/OneDrive - The University of Auckland/PhD/PhD_Write_up/article1/casestudy.pdf", 
      height = 8, width = 5.5)
  plot_map_cy(z_csy = B_cst[,,(1:5)*12-11], inla_proj = proj_om,
              lab_y = seq(1, 5, 1), leglim_c = rbind(c(0,1), c(0,1)), ncols = 500,
              height = 8, p_width = 0.2, mar = rep(0.25, 4),
              legend_type = 2, species_names = c("S. aequilatera", "M. murchisoni"),
              xlim = xx, ylim = yy) #235
  dev.off()
}

