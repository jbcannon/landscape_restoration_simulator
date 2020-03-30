#############################################################################################
#### Script for simulating heterogenous restoration treatments using LANDFIRE data
#### Author: Jeffery Cannon (Jeffery.Cannon@jonesctr.org)
#### Institution: The Jones Center at Ichauway
#### Date Created:  21 November 2017
#### Last Modified: 30 March 2020
#############################################################################################
# The goal of this script is to simulate restoration treatments in individual catchments by
# specifying target canopy cover (mean, sd) for wet and dry regions (based on mois), and
# reduce canopy cover while considering contraints of the current forest vegetation structure
#############################################################################################

#######################################START SET UP##########################################
packages <- c('raster', 'rgdal', 'rgeos', 'sp', 'doParallel', 'foreach')
for (package in packages) {
  if (suppressMessages(!require(package, character.only = T))) {
    install.packages(package, repos = 'https://cran.mtu.edu/')
    suppressMessages(library(package, character.only = T))
  }
}

#---> Load data (cite sources)
# Canopy cover from LANDFIRE, coverted to proprtion (x/100)
CC  <-  raster('data/canopy_cover_raster.tif') / 100

# Soil moisture index (Cacluated from solar radiation and topographic wetness, See Cannon et al. 2020)
mois   <-  raster('data/moisture_raster.tif') 

# Treatment boundaries of study area
catches <- readOGR('data/boundary.shp') #catchments

# Table of historic canopy fuel distribution parameters
gamma_pars <- read.csv('data/desired_cc_dist_pars.csv')

#---> Set parameters for simulation
# Moisture index threshold, values > threshold are considered 'wet', see reclass table
mois_thresh <- 0.2 # division between 'wet' and 'dry' (0.2, study area median/mean)
reclass_tab = matrix(c(0, mois_thresh, mois_thresh, 1000, 0, 1), ncol = 3)

# Directory for output rasters (temporary)
trt_results_dir <- 'data/output'
raster_results_filename = 'restoration_simulation_example.tif'
fig_folder = 'data/figures'
scratch_dir <- 'data/scratch'

#########################
# FUNCTION: rank_raster()
#
# This function takes a raster and returns the ranks of the values, used for ranking wettest to driest
# portions of a catchment and matches them for the lowest/highest canopy cover
#########################
rank_raster <- function(r) {
  mask = !is.na(r)
  tmp_r = r
  tmp_r[] = rank(values(r))
  max = cellStats(mask(tmp_r, mask), stat = 'mean')
  tmp_r[!mask] = NA
  return(tmp_r)
}

#########################
# FUNCTION: mois_goal()
#
# This function resassigns canopy cover in a catchment based on treatment desired conditions and 
# soil moisture rankings.
#########################
mois_goal <-
  function(cover_mean_target,
           cc_a,
           cc_b,
           mois,
           cc,
           maxtry = 50) {
    try = 1
    treated = FALSE
    goal_mean = cover_mean_target
    if (cellStats(cc, stat = 'mean') < cover_mean_target) {
      print('Initial canopy cover lower than target. Area left untreated.')
      return(cc)
      break
    }
    while (treated == FALSE) {
      #print(paste('Attempting treatment simulation', try, 'SD =', round(cover_sd)))
      n = cellStats(!is.na(mois), stat = 'sum')
      cc_dist = rgamma(n, shape = cc_a, rate = cc_b)
      cc_dist[cc_dist < 0] <- 0
      cc_dist[cc_dist > 1] <- 1
      cc_dist = cc_dist[order(cc_dist)]
      
      rank_mois = rank_raster(mois)
      mois_vals = values(rank_mois)
      mois_vals = mois_vals[!is.na(mois_vals)]
      
      cc_out = rank_mois
      cc_out[!is.na(cc_out)] = cc_dist[mois_vals]
      
      #sim treatment
      mult = cc_out / cc
      mult[is.na(mult)] <- 1
      mult[mult > 1] = 1
      tent_mean = cellStats(cc * mult, stat = 'mean')
      if (tent_mean < goal_mean & try < maxtry) {
        cc_a = cc_a * 1.02
        try = try + 1
        #print(paste('CC mean =', round(tent_mean,3)))
      } else {
        treated = TRUE
      }
    }
    if (tent_mean < goal_mean)
      print(paste('***Treatment aborted after', try, 'attempts***'))
    else {
      print(paste('        Treatment simulation successful after', try, 'attempts'))
    }
    return(cc * mult)
  }

#########################
# FUNCTION: goal_output()
#########################
goal_output <-
  function(mois,
           mois_rcl,
           mask_value,
           CC,
           cc_mean,
           cc_a,
           cc_b) {
    mois_mask = mois_rcl
    mois_mask[mois_mask == (1 - mask_value)] <- NA
    mois_mask[is.na(mois_mask) == FALSE] <- 1
    mois = mask(mois, mois_mask)
    mois_cc = mask(CC, mois_mask)
    
    # check to ensure there are areas included
    mois_cells = sum(values(mois_mask), na.rm = TRUE)
    if (mois_cells > 0) {
      goal = mois_goal(
        cover_mean_target = cc_mean,
        cc_a = cc_a,
        cc_b = cc_b,
        mois = mois,
        cc = mois_cc
      )
    } else {
      goal = mois_mask
    } #If no cells, then return mask (multiplier = 1)
    return(goal)
  }

# Create output directories
if(!dir.exists(trt_results_dir)) dir.create(trt_results_dir)
if(!dir.exists(scratch_dir)) dir.create(scratch_dir)
if(!dir.exists(fig_folder)) dir.create(fig_folder)
raster_results_location = paste(trt_results_dir, '/',raster_results_filename, sep = '')

# Raster processing and parallel processing options
rasterOptions(maxmemory = 10 ^ 9)
useCores <- detectCores() - 2
cl       <- makeCluster(useCores)
registerDoParallel(cl)

# Create raster of study areas
catches_ras = CC
catches_ras[] = NA
catches_ras = rasterize(catches, catches_ras, 'FEATURE')
catches_ras[!is.na(catches_ras)]=1

#--> Begin clock on analysis
start_time = Sys.time()
########################################END SET UP###########################################

######################################START ANALYSIS#########################################
#---> Clean up, project, and crop data
# Set default projection and project all boundaries to match.
proj    <- proj4string(CC)
bnd <- spTransform(catches, proj)

# Crop rasters to study area
CC  <- crop(CC,  bnd)
mois <- crop(mois, spTransform(bnd, proj4string(mois)))

#---> Loop through all catchments (in parallel), simulate treatment, create figure, update landscape cover
REST_SIM <- foreach(
  i = catches$FEATURE,
  .combine = rbind,
  .init = data.frame(FEATURE = numeric(),
                     result = character())
) %dopar% {
  #load libraries for each core
  for (package in packages) {
    if (suppressMessages(!require(package, character.only = T))) {
      install.packages(package, repos = 'https://cran.mtu.edu/')
      suppressMessages(library(package, character.only = T))
    }
  }
  
  # subset to single catchment, crop mois, CC, and reclass mois, and mask out non-forested areas
  tmp_bnd = subset(catches, FEATURE == i)
  tmp_mask = crop(catches_ras, tmp_bnd)
  tmp_mask = mask(tmp_mask, tmp_bnd)
  tmp_CC  = mask(crop(CC, tmp_bnd), tmp_mask)
  tmp_mois = crop(mois, spTransform(tmp_bnd,proj4string(mois)))
  tmp_mois = projectRaster(tmp_mois, tmp_CC)
  tmp_mois = mask(tmp_mois, tmp_mask)
  tmp_mois_rcl = reclassify(tmp_mois, rcl = reclass_tab)
  
  # Extract CC gamma distribution pars based on zone/mois
  wet = subset(gamma_pars,
               gamma_pars$zone == as.factor(tmp_bnd$zone) &
                 mois == 'wetter')
  dry = subset(gamma_pars, zone == tmp_bnd$zone & mois == 'drier')
  
  # If catchment contains no pixels, then skip
  cells <- cellStats(na.omit(tmp_mois), stat = 'sum', na.rm = TRUE)
  if (cells == 0) {
    df_out <- data.frame(FEATURE = i, result = 'ERROR, NO CELLS')
    return(df_out)
    
  } else {
    goal_dry = goal_output(
      mois = tmp_mois,
      mois_rcl = tmp_mois_rcl,
      mask_value = 0,
      CC = tmp_CC,
      cc_mean = wet$cc_mn,
      cc_a = wet$cc_a,
      cc_b = wet$cc_b
    )
    
    goal_wet = goal_output(
      mois = tmp_mois,
      mois_rcl = tmp_mois_rcl,
      mask_value = 1,
      CC = tmp_CC,
      cc_mean = wet$cc_mn,
      cc_a = wet$cc_a,
      cc_b = wet$cc_b
    )
    
    #combine wet and dry goals in to one
    goal = merge(goal_wet, goal_dry)
    
    # Calculate multipliers necessary to reach goal and calculate post-treatment
    cc_mult <- goal / tmp_CC
    cc_mult[is.na(cc_mult)] <- 1
    cc_mult[cc_mult < 0] <- 0
    cc_mult[cc_mult > 1] <- 1
    outcome <- tmp_CC * cc_mult
    
    res_out <- paste(scratch_dir, '/restsim_', i, '.tif', sep = '')
    writeRaster(outcome, res_out, overwrite = TRUE)
    df_out = data.frame(FEATURE = i, result = 'written to disk')
    
    {pdf(paste0(fig_folder, '/', i, '.pdf'), width = 7.5, height = 7.5)
    par(mfrow = c(2,2), mar = c(2,2,3,0), oma = c(0,0,2.5,0))
    # Untreated
    plot(tmp_CC, main = 'Untreatead\nCanopy cover', col = rev(terrain.colors(11)), breaks = seq(0,1,by=0.1))
    plot(tmp_bnd, add = TRUE, lwd = 2)
    scalebar(1000, type = 'bar', below = 'km', label = c(0,'',1))
    cc_mean = round(cellStats(tmp_CC, stat = 'mean'),3)
    cc_sd   = round(cellStats(tmp_CC, stat = 'sd'),3)
    legend('topleft', paste('Mean CC:', cc_mean, '±', cc_sd), bty = 'n')
    
    #MOIS
    plot(tmp_mois, col = rev(topo.colors(100)), main = 'Moisture index')
    plot(tmp_bnd, add = TRUE, lwd = 2)
    scalebar(1000, type = 'bar', below = 'km', label = c(0,'',1))
    
    #RES
    plot(outcome, main = 'Restoration treatment\nCanopy cover', col = rev(terrain.colors(11)), breaks = seq(0,1,by=0.1))
    plot(tmp_bnd, add = TRUE, lwd = 2)
    scalebar(1000, type = 'bar', below = 'km', label = c(0,'',1))
    cc_mean = round(cellStats(outcome, stat = 'mean'),3)
    cc_sd   = round(cellStats(outcome, stat = 'sd'),3)
    legend('topleft', paste('Mean CC:', cc_mean, '±', cc_sd), bty = 'n')
    
    #RES DIFF
    plot(outcome - tmp_CC, main = 'Restoration treatment\nCanopy reduction', col = heat.colors(11), breaks = seq(0,-1,by=-0.1))
    plot(tmp_bnd, add = TRUE, lwd = 2)
    scalebar(1000, type = 'bar', below = 'km', label = c(0,'',1))
    
    mtext(paste('Catchment #', i), outer = TRUE)
    
    dev.off()}
    return(df_out)
  }
}

# Get list of all rasters in raster folder
ras_list = list()
ras_out = grep('restsim_\\d', list.files(scratch_dir), value = TRUE)
for (i in 1:length(ras_out)) {
  tmp_ras = raster(paste(scratch_dir, '/', ras_out[i], sep = ''))
  ras_list[[i]] = tmp_ras
}
# Combine all rasters into single raster
RES_out = do.call(merge, ras_list)

# Burn in treated catchments to untreated landscape
RES = CC
RES_out = extend(RES_out, CC)
treated_cells = Which(!is.na(RES_out), cells = TRUE)
RES[treated_cells] <- RES_out[treated_cells]
RES=RES*100
#######################################END ANALYSIS##########################################

######################################START OUTPUTS##########################################
end_time = Sys.time()
paste('Total analysis time:', round(difftime(end_time, start_time, units = 'mins'), 2), 'minutes')

# Write output raster to disk
writeRaster(RES, filename = raster_results_location, overwrite = TRUE)

#Delete scratch directory
unlink(scratch_dir, recursive = TRUE)
#######################################END OUTPUTS###########################################

