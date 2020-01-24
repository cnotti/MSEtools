# unload TMB dlls
unloadTMBdll = function(dll) {
  dlls = getLoadedDLLs()
  isTMBdll = function(dll) !is(try(getNativeSymbolInfo("MakeADFunObject",
                                                       dll), TRUE), "try-error")
  TMBdll = sapply(dlls, isTMBdll)
  loaded = names(TMBdll[TMBdll == TRUE])
  #if (dll %in% loaded) dyn.unload(paste0("src/", dll))
  if (dll %in% loaded) dyn.unload(paste0("src/", dynlib(dll)))
}

# parameter transformations
logit = function(p) log(p/(1 - p))
invlogit = function(p) exp(p)/(1 + exp(p))
logitST = function(p, a, b) logit((p - a)/(b - a))
invlogitST =  function(p, a, b) a + (b - a)*invlogit(p)

# logistic curves
logistic = function(l, l50, l95)
  1/(1 + exp(-log(19)*((l - l50)/(l95 - l50))))
dec_logistic = function(l, l50, l95, delta_max)
  delta_max/(1 + exp(log(19)*((l - l50)/l95)))

# calculate area for each mesh node
inla.mesh.dual = function(mesh) {
  # estimate the area of each Voronoi polygon
  if (mesh$manifold == 'R2') {
    ce = t(sapply(1:nrow(mesh$graph$tv), function(i)
      colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
    library(parallel)
    pls = mclapply(1:mesh$n, function(i) {
      p = unique(Reduce('rbind', lapply(1:3, function(k) {
        j = which(mesh$graph$tv[,k] == i)
        if (length(j) > 0)
          return(rbind(ce[j, , drop = FALSE],
                       cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1],
                             mesh$loc[mesh$graph$tv[j, k], 2] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
        else return(ce[j, , drop = FALSE])
      })))
      j1 = which(mesh$segm$bnd$idx[,1] == i)
      j2 = which(mesh$segm$bnd$idx[,2] == i)
      if ((length(j1) > 0) | (length(j2) > 0)) {
        p = unique(rbind(mesh$loc[i, 1:2], p,
                         mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                           mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2,
                         mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                           mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
        yy = p[,2] - mean(p[,2])/2 - mesh$loc[i, 2]/2
        xx = p[,1] - mean(p[,1])/2 - mesh$loc[i, 1]/2
      }
      else {
        yy = p[,2] - mesh$loc[i, 2]
        xx = p[,1] - mesh$loc[i, 1]
      }
      Polygon(p[order(atan2(yy,xx)), ])
    })
    return(SpatialPolygons(lapply(1:mesh$n, function(i)
      Polygons(list(pls[[i]]), i))))
  }
  else stop("It only works for R2!")
}

# calculate tstart and tstop using t_catch
calc_t = function(t_catch) {
  runs = rle(t_catch)
  if (runs$values[1] == 0) {
    t0 = 1:(runs$lengths[1] + 1)
  } else {
    t0 = c(1, runs$lengths[1] + 1)
  }
  for (t in 2:length(runs$lengths)) {
    if (runs$values[t] == 0) {
      t0 = c(t0[-length(t0)], tail(t0,1):(cumsum(runs$lengths)[t] + 1))
    } else {
      t0 = c(t0, cumsum(runs$lengths)[t] + 1)
    }
  }
  list(tstart = t0[-length(t0)]-1, tstop = t0[-1]-1)
}


# calculate depths from bath matrix contour lines
#' @export
get_depth = function(lonlat, bath, levels = 0:(-20)) {
  # bath much be matrix of depths where row/col names = lon/lat coords
  lon = unique(as.numeric(rownames(bath)))
  lat = unique(as.numeric(colnames(bath)))
  cont = contourLines(lon, lat, bath, levels = levels)
  c.lines = maptools::ContourLines2SLDF(cont)
  distmat = rgeos::gDistance(c.lines, sp::SpatialPoints(lonlat), byid = c(T, T))
  invdepth_s = apply(distmat, 1, function(y) which(y == min(y))[1]) - 1
  depth_s = invdepth_s - max(invdepth_s)
  depth_s
}

# create depth polygon by joining depth contour lines
#' @export
make_depth_poly = function(bath, depth_range = c(0, -10), 
                           domain_outer, reverseLine2 = TRUE) {
  lon = unique(as.numeric(rownames(bath)))
  lat = unique(as.numeric(colnames(bath)))
  cont = contourLines(lon, lat, bath, levels = depth_range)
  c.lines = maptools::ContourLines2SLDF(cont)
  depth_lines = gIntersection(c.lines, domain_outer)
  # select which 2 lines to join
  coords = coordinates(depth_lines)[[1]]
  n = length(coords)
  if (n != 2) {
    plot(depth_lines)
    for (line in 1:length(coords)) {
      lines(coords[[line]][,1], coords[[line]][,2],
            col = rainbow(n)[line])
    }
    legend("topright", legend = 1:n, col = rainbow(n), lty = 1)
    drop = as.numeric(readline(prompt = "lines to drop: "))
    Lines = coords[-drop]
  }
  else {
    Lines = coords
  }
  if (reverseLine2 == FALSE) {
    domain_lines = do.call("rbind", list(Lines[[1]], 
                                         Lines[[2]], 
                                         Lines[[1]][1,]))
  }
  else if (reverseLine2 == TRUE) {
    domain_lines = do.call("rbind", list(Lines[[1]], 
                                         Lines[[2]][nrow(Lines[[2]]):1,], 
                                         Lines[[1]][1,]))
  }
  SpatialPolygons(list(Polygons(list(Polygon(domain_lines)), "Domain")))
}


# plot results on map
#' @export
plot_map = function(z_s, inla_proj, land, bath, legend_type = 1,
                    p, y, leglim = c(min(z_s), max(z_s)),
                    xlim, ylim, np = 12, ncols = 100, pretty = TRUE, 
                    cex.text = 1, col.text = "black",
                    bg = "dark blue", ...) {
  # plot image
  cols = colorRampPalette(c(bg, 
                            "light yellow", "yellow",
                            "red", "maroon"))(ncols)

  lon = inla_proj$x
  lat = inla_proj$y
  #if (pretty == TRUE) z_s = ifelse(z_s <= leglim[2], z_s, leglim[2])
  zmat = inla.mesh.project(inla_proj, z_s)
  if (pretty == TRUE) {
    zmat[zmat < leglim[1]] = leglim[1]
    zmat[zmat > leglim[2]] = leglim[2]
  }
  
  if (missing(xlim)) {
    xlim = c(min(lon), max(lon))
  }
  if (missing(ylim)) {
    ylim = c(min(lat), max(lat))
  }
  if (!missing(land)) {
    plot(land, col = "dark green", bg = bg,
         xlim = xlim, ylim = ylim, xpd = FALSE)
    image(lon, lat, zmat,
          col =  cols, breaks = seq(leglim[1], leglim[2], len = ncols + 1),
          add = TRUE)
    plot(land, col = "dark green", add = TRUE, xpd = FALSE)
  } else {
    image(lon, lat, zmat, xaxt = "n", yaxt = "n",
          col =  cols, breaks = seq(leglim[1], leglim[2], len = ncols + 1),
          asp = 1, xlim = xlim, ylim = ylim)
    plot(domain, add = TRUE, bg = "bg")
    }
  if (!missing(bath)) {
    plot(bath, deep = -22, shallow = -2, step=10, drawlabels = TRUE,
         col = "grey70", add = T)
    scaleBathy(bath, deg = 0.05, x = "bottomright", inset = 8, col = "black")
  }
  
  
  usr = par("usr")
  xl = usr[1]; xr = usr[2]; yb = usr[3]; yt = usr[4]

  domain_coords = cbind(c(xlim[1], xlim[2], xlim[2], xlim[1], xlim[1]),
                        c(ylim[1], ylim[1], ylim[2], ylim[2], ylim[1]))
  outer_coords = cbind(c(xl, xr, xr, xl, xl),
                       c(yb, yb, yt, yt, yb))
  domain_sp = SpatialPolygons(list(Polygons(list(Polygon(domain_coords, hole = TRUE),
                                                 Polygon(outer_coords)),
                                            "Domain")))
  
  #plot(domain_sp, add = TRUE, col = bg, border = bg)
  #box()
  
  # add legend
  if (legend_type == 1) {
    h = yt - yb; w = xr - xl
    rl = xl + 0.85*w; rr = xl + 0.95*w
    rb = head(seq(yb + 0.3*h, yt - 0.3*h, len = ncols), -1)
    rt = tail(seq(yb + 0.3*h, yt - 0.3*h, len = ncols), -1)
    rect(rl, rb, rr, rt, col = cols, border = NA)
    rect(rl, head(seq(min(rb), max(rt), len = 5), -2), rr,
         tail(seq(min(rb), max(rt), len = 5), -2))
    text(rl - (rr-rl)*1, seq(min(rb), max(rt), len = 5)[2:4],
         sprintf("%.2f", round(seq(leglim[1], leglim[2], len = 5), 2))[2:4],
         cex = cex.text, col = col.text)
  }
  
  
  if (!missing(p) & !missing(y)) {
    text(xl + 0.25*w, yt - 0.05*h, cex = 2,
         labels = paste0(sprintf(paste0("%0",nchar(np),"d"),
                                (1:np)[p]), "/", y), ...)
  } else if (missing(p) & !missing(y)) {
    text(xl + 0.25*w, yt - 0.05*h, cex = 2,
         labels = paste0("year: ", y), ...)
  }
}

# function to plot spatio-temporal map as gif
#' @export
plot_map_gif = function(z_i_st, inla_proj, land, bath, p1 = 1, y1 = 2018,
                        leglim_i = lapply(z_i_st, function(z_st) c(0, max(z_st))),
                        nt = ncol(z_i_st[[1]]), np = 12,
                        ncols = 100, filename, fps = 10, ...) {
  ni = length(z_i_st)
  y = y1
  p = p1
  fig = magick::image_graph(width = 400*ni, height = 580, res = 96)
  par(mfrow = c(1, ni), mar = c(0,0,0,0), oma = c(0,0,0,0))

  for (t in 1:nt) {
    for (i in 1:ni) {
        plot_map(z_s = z_i_st[[i]][,t], inla_proj = inla_proj, land = land,
               bath = bath, leglim = leglim_i[[i]], p = p, y = y, np = np,
               ncols = ncols, ...)
    }

    if (p < np) {
      p = p + 1
    } else {
      y = y + 1
      p = 1
    }
  }
  dev.off()
  anim = magick::image_animate(fig, fps = fps)
  if (missing(filename)) {
    return(anim)
  } else {
    magick::image_write(anim, filename)
  }
}


#' @export
plot_map_cy = function(z_csy, inla_proj, land, bath, lab_y,
                      p, y, leglim_c = NULL, cex.text = 1,
                      xlim, ylim, col.text = "black",
                      np = 12, ncols = 100, pretty = TRUE,
                      species_names = NULL, height = 5,
                      p_width = 0.3, mar = rep(0.2, 4),
                      legend_type = 1, legend_size = 0.4,
                      bg = "dark blue") {
  
  old.par = par(no.readonly = TRUE)
  
  # get dims
  nc = dim(z_csy)[1]
  ny = dim(z_csy)[3]
  # set up plot region
  par(mfrow = c(ny, nc), mar = mar, oma = c(0,4,0,0), xpd = F)
  
  if (legend_type == 1) {
    Lmat = matrix(0, ncol = nc+1, nrow = ny+1)
    Lmat[1,2:(nc+1)] = 1:nc
    Lmat[2:(ny+1),1] = (nc+1):(nc+ny)
    Lmat[2:(ny+1),2:(nc+1)] = (nc+ny+1):(((nc+1)*(ny+1)) - 1)
    layout(mat = t(Lmat), 
           heights = c(lcm(0.65), rep(lcm(height), nc)), 
           widths = c(lcm(0.65), rep(lcm(height*p_width), ny)))
  } else if (legend_type == 2) {
    Lmat = matrix(0, ncol = nc+1, nrow = ny+2)
    Lmat[1,2:(nc+1)] = 1:nc
    Lmat[2:(ny+1),1] = (nc+1):(nc+ny)
    Lmat[2:(ny+1),2:(nc+1)] = (nc+ny+1):(((nc+1)*(ny+1)) - 1)
    
    Lmat[ny+2,2:(nc+1)] = (nc+1)*(ny+1)
    
    layout(mat = t(Lmat), 
           heights = c(lcm(0.65), rep(lcm(height), nc)), 
           widths = c(lcm(0.65), rep(lcm(height*p_width), ny), lcm(2)))
  }
  
  
  # plot labels
  if (is.null(species_names)) species_names = 1:nc
  for (c in 1:nc) {
    plot.new()
    rect(-2, -2, 2, 2, border = "black", col = "light gray")
    text(0.5, 0.5, species_names[c], srt = 90)
    box()
  }

  for (y in 1:ny) {
    plot.new()
    rect(-2, -2, 2, 2, border = "black", col = "light gray")
    text(0.5, 0.5, paste("y =", lab_y[y]))
    box()
  }
  if (is.null(leglim_c)) {
    leglim_c = matrix(0, nrow = nc, ncol = 2)
    for (c in 1:nc) { 
      leglim_c[c,] = c(min(z_csy[c,,]), max(z_csy[c,,])*0.9)
    }
  }
  # plot maps
  for (c in 1:nc) {  
    for (y in 1:ny) {
      plot_map(z_csy[c,,y], inla_proj = inla_proj, leglim = leglim_c[c,],
               land = land, pretty = pretty, cex.text = cex.text, 
               xlim = xlim, ylim = ylim, col.text = col.text,
               legend_type = legend_type, bg = bg)
      #box()
    }
  }
  if (legend_type == 2) {
    par(xpd = NA)
    
    cols = colorRampPalette(c(bg, 
                              "light yellow", "yellow",
                              "red", "maroon"))(ncols)
    
    plot.new()
    usr = par("usr")
    xl = usr[1]; xr = usr[2]; yb = usr[3]; yt = usr[4]

    h = yt - yb; w = xr - xl
    rl = xl + 0.6*w; rr = xl + 0.8*w
    rb = head(seq(yb + legend_size*h, yt - legend_size*h, len = ncols), -1)
    rt = tail(seq(yb + legend_size*h, yt - legend_size*h, len = ncols), -1)
    rect(rl, rb, rr, rt, col = cols, border = NA)
    rect(rl, head(seq(min(rb), max(rt), len = 5), -2), rr,
         tail(seq(min(rb), max(rt), len = 5), -2))
    text(rl - (rr-rl)*1, seq(min(rb), max(rt), len = 5)[2:4],
         sprintf("%.2f", round(seq(leglim_c[1,1], leglim_c[1,2], len = 5), 2))[2:4],
         cex = cex.text, col = col.text)
  }
  par(old.par)
}




# harvest function version: -1
# initialize obs
if (FALSE) {
  f_i = c_i = catch_i = NULL
  limit_c = data_om$limit_c
  selectivityF_cl = t(do.call("cbind", om_rep$selectivityF_c_l))
  weight_cl = data_om$weight_cl
  areaS = data_om$areaOmegastar/data_om$ns
  targetfprob = 0.75  # maximum harvest rate of a given area
  selweight_cfl = N_cfl = with(data_om, array(0, dim = c(nc, nf, nl)))
  for (c in 1:data_om$nc) {
    selweight_cfl[c,,] = matrix(selectivityF_cl[c,] * weight_cl[c,], 
                                nrow = data_om$nf, ncol = data_om$nl, byrow = TRUE)
  }
  system.time({
    for (c in 1:data_om$nc) {
      N_cfl[c,,] = as.matrix(exp(data_om$A_fs %*% log(data_om$N_csl[c,,] / areaS + 1e-06)))
    }  
    B_cf = apply(N_cfl * selweight_cfl, 1:2, "sum")
    p_cf = B_cf/rowSums(B_cf)
    
    
    for (c in 1:data_om$nc) {
      fprob_f = ifelse(p_cf[c,] < mean(p_cf[c,]), 1 - targetfprob, targetfprob)
      C_n = 0
      n = floor(limit_c[c]/mean(B_cf[c,]))
      while (sum(C_n) < limit_c[c]) {
        if (c == 1) {
          # all f sites available
          C_n = sample(B_cf[c,], size = n, 
                       replace = F, prob = fprob_f)
        } else {
          # remove previously sampled sites f_i
          unfished = -unique(f_i)
          C_n = sample(B_cf[c, unfished], size = n, 
                       replace = F, prob = fprob_f[unfished])
        }
        n = n * 1.2
      }
      # save vector of catch and book-keeping vars
      ni = which.min(abs(cumsum(C_n) - limit_c[c]))
      C_i = C_n[1:ni]
      f_i = c(f_i, match(C_i, B_cf[c,]))
      C_ci = B_cf[,f_i]
      c_i = c(c_i, rep(c, ni))
      catch_i = c(catch_i, C_i) 
      # throw back catch > limit
      for (cstar in (1:data_om$nc)[-c]) {
        if (cstar > c) {
          nretain = which.min(abs(cumsum(C_ci[cstar,]) - limit_c[cstar]))
          C_i = C_ci[cstar,1:nretain]
          # book-keeping
          f_i = c(f_i, match(C_i, B_cf[cstar,]))
          c_i = c(c_i, rep(cstar, nretain))
          catch_i = c(catch_i, C_i)  
        }
      }
    }
  })
}
