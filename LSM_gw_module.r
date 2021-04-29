library(HBV.IANIGLA) # for running HBV model
#library(prism) # for downloading PRISM climate data (to 1981); might replace by just usign mopex data
library(hddtools)		# for mopex data
library(hydroGOF)		# for nse calculations
library(dataRetrieval)	# for streamflow data (I think)
library(data.table)
library(sf)


# because NOAA had a facilities issue on 9MAR2021 and MOPEX data is temporarily unavailable via webservice,
# I am inserting a temporary "fix" to allow me to read MOPEX data locally
# the tag TEMPFIX denotes where I am doing this
# ultimately, TEMPFIX should be removed when the data is again available via webservice
  #!!!! TEMPFIX
#mopey_catalogue = catalogueMOPEX()							# reading in the mopex catalogue
mopey_catalogue = gsub(".dly",'', list.files("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\Daily")) 
  #!!!! TEMPFIX





####################################################################
### part 1: monte carlo of w/ and w/o new groundwater module


      #######################################################
      ### step 1: reading in data

# reading in the empirical data calculated for all basins
setwd("C:\\Users\\arik\\Documents\\PhD Research\\D3")
aq_chars = st_read("basin_and_aq_chars_with_polygons_jan2020.gpkg")

# identify the watersheds where we have MOPEX data and are recession data from "reference" watersheds 
ws_ihave = NULL
#!!! TEMPFIX
#for(gg in 1:nrow(mopey_catalogue)){
#if(mopey_catalogue[gg, "USGS_ID"] %in% aq_chars$gage) 	ws_ihave = c(ws_ihave, gg)
#}
#head(mopey_catalogue[ws_ihave, ])
for(gg in mopey_catalogue){
  if(gg %in% aq_chars$gage) 	ws_ihave = c(ws_ihave, gg)
}
#!!! TEMPFIX


# initializing the dataframe for capturing model performance
model_perf = data.frame(
  gage_id = rep(NA,length(ws_ihave)),
  hbv_nse = rep(NA,length(ws_ihave)),
  new_nse = rep(NA,length(ws_ihave)),
#  hbv_log_nse = rep(NA,length(ws_ihave)),
#  new_log_nse = rep(NA,length(ws_ihave)),
  hbv_Q5_nse = rep(NA,length(ws_ihave)),
  new_Q5_nse = rep(NA,length(ws_ihave)),
  hbv_Q10_nse = rep(NA,length(ws_ihave)),
  new_Q10_nse = rep(NA,length(ws_ihave)),
  hbv_Q20_nse = rep(NA,length(ws_ihave)),
  new_Q20_nse = rep(NA,length(ws_ihave)),
  hbv_kge = rep(NA,length(ws_ihave)),
  new_kge = rep(NA,length(ws_ihave)),
  hbv_Q5_kge = rep(NA,length(ws_ihave)),
  new_Q5_kge = rep(NA,length(ws_ihave)),
  hbv_Q10_kge = rep(NA,length(ws_ihave)),
  new_Q10_kge = rep(NA,length(ws_ihave)),
  hbv_Q20_kge = rep(NA,length(ws_ihave)),
  new_Q20_kge = rep(NA,length(ws_ihave)),
  hbv_rmse = rep(NA,length(ws_ihave)),
  new_rmse = rep(NA,length(ws_ihave)),
  hbv_Q5_rmse = rep(NA,length(ws_ihave)),
  new_Q5_rmse = rep(NA,length(ws_ihave)),
  hbv_Q10_rmse = rep(NA,length(ws_ihave)),
  new_Q10_rmse = rep(NA,length(ws_ihave)),
  hbv_Q20_rmse = rep(NA,length(ws_ihave)),
  new_Q20_rmse = rep(NA,length(ws_ihave)),
  hbv_mae = rep(NA,length(ws_ihave)),
  new_mae = rep(NA,length(ws_ihave)),
  hbv_Q5_mae = rep(NA,length(ws_ihave)),
  new_Q5_mae = rep(NA,length(ws_ihave)),
  hbv_Q10_mae = rep(NA,length(ws_ihave)),
  new_Q10_mae = rep(NA,length(ws_ihave)),
  hbv_Q20_mae = rep(NA,length(ws_ihave)),
  new_Q20_mae = rep(NA,length(ws_ihave)),
  hbv_Q5_pctB = rep(NA, length(ws_ihave)),
  new_Q5_pctB = rep(NA, length(ws_ihave)),
  hbv_10_pctB = rep(NA,length(ws_ihave)),
  new_10_pctB = rep(NA,length(ws_ihave)),
  hbv_20_pctB = rep(NA,length(ws_ihave)),
  new_20_pctB = rep(NA,length(ws_ihave)),
  # hbv modules
  hbv_sfcf = rep(NA,length(ws_ihave)),	#snowfall correction factor [-]
  hbv_tr   = rep(NA,length(ws_ihave)),	#solid and liquid precipitation threshold temperature [C]
  hbv_tt   = rep(NA,length(ws_ihave)),	#melt temperature [C]
  hbv_fm   = rep(NA,length(ws_ihave)),	#snowmelt factor [mm/C]
  hbv_fi   = rep(NA,length(ws_ihave)),	#icemelt factor [mm/C]
  hbv_fic  = rep(NA,length(ws_ihave)),	#debris-covered icemelt factor [mm/C]
  # soil module
  hbv_fc   = rep(NA,length(ws_ihave)),
  hbv_lp   = rep(NA,length(ws_ihave)),
  hbv_beta_soils = rep(NA,length(ws_ihave)),
  # routing module
  hbv_k0   = rep(NA,length(ws_ihave)),
  hbv_k1   = rep(NA,length(ws_ihave)),
  hbv_k2   = rep(NA,length(ws_ihave)),	
  hbv_uz1  = rep(NA,length(ws_ihave)),
  hbv_perc = rep(NA,length(ws_ihave)),
  # hbv modules
  new_sfcf = rep(NA,length(ws_ihave)),	#snowfall correction factor [-]
  new_tr   = rep(NA,length(ws_ihave)),	#solid and liquid precipitation threshold temperature [C]
  new_tt   = rep(NA,length(ws_ihave)),	#melt temperature [C]
  new_fm   = rep(NA,length(ws_ihave)),	#snowmelt factor [mm/C]
  new_fi   = rep(NA,length(ws_ihave)),	#icemelt factor [mm/C]
  new_fic  = rep(NA,length(ws_ihave)),	#debris-covered icemelt factor [mm/C]
  # soil module
  new_fc   = rep(NA,length(ws_ihave)),
  new_lp   = rep(NA,length(ws_ihave)),
  new_beta_soils = rep(NA,length(ws_ihave)),
  # routing module
  new_k0   = rep(NA,length(ws_ihave)),
  new_k1   = rep(NA,length(ws_ihave)),
  new_k2   = rep(NA,length(ws_ihave)),	
  new_uz1  = rep(NA,length(ws_ihave)),
  new_perc = rep(NA,length(ws_ihave))
)
  


### model 1: setting the parameter ranges for the HBV runoff modules
# snow and glacier module
  sfcf = c(runif(500, .2, 1), runif(500, 1, 3)) #1, 2)	#snowfall correction factor [-]
  tr   = runif(1000, -6, 5)#-1, 3)	#solid and liquid precipitation threshold temperature [C]
  tt   = runif(1000, -5, 6)#0, 3)	#melt temperature [C]
  fm   = c(runif(500, .2, 1.5), (runif(500, 1.5, 8)))#1, 4)	#snowmelt factor [mm/C]
  fi   = c(runif(500, .2, 1.5), (runif(500, 1.5, 10)))#4, 8)	#icemelt factor [mm/C]
  fic  = runif(1000, 2, 10)#6, 6)	#debris-covered icemelt factor [mm/C]

# soil module
#  fc   = runif(1000, 30, 3000)#100, 200)
  fc   = c(runif(500, 25, 150), (runif(500, 150, 1200)))
  lp   = runif(1000, .2, 1)#0.5, 1)   # parameter to actual ET
  beta_soils = runif(1000, 1, 3)

# routing module
  k0   = c(runif(500, .05, .5), (runif(500, .5, 1)))#runif(1000, .01, 1)#0.09, 0.1)
  k1   = c(runif(500, .005, .09), (runif(500, .09, .5)))#runif(1000, .001, .5)#0.05, 0.07)
  k2   = c(runif(500, .0001, .01), (runif(500, .01, .1)))#runif(1000, .0001, .1)#0.05)	
  uz1  = c(runif(500, .22, 10), (runif(500, 10, 40)))#runif(1000, .1, 100)#1, 5) #max flux rate from STZ to SUZ in mm/d
#  uz1  = runif(1000, 10, 100)
  perc = c(runif(500, .1, .5), (runif(500, .5, 20)))#runif(1000, .001, 100)#.8, 2)  # max flux rate from SUZ to SLZ in m/d
#  perc = runif(1000, 10, 100)




  
for(allws in c(1:59, 60:length(ws_ihave))){
  
  #################################################################################
  # initializing data for our gw module
  
  #!!! TEMPFIX
    which_ws = allws
    #this_ws = ws_ihave[which_ws]				# identifying the ws in mopex
    #this_ws_props = which(aq_chars$gage == mopey_catalogue[this_ws, "USGS_ID"])
    this_ws_props = which(aq_chars$gage == ws_ihave[which_ws])
    model_perf$gage_id[allws] = ws_ihave[which_ws]
    #!!! TEMPFIX
    
if(!is.na(aq_chars$B_depth_dry[this_ws_props])){
    D_dry = aq_chars$B_depth_dry[this_ws_props]  # in m?
    D_wet = aq_chars$B_depth_wet[this_ws_props] # in m?
    K_dry = aq_chars$B_K2_dry[this_ws_props] * (60 * 60 * 24) / 100 # convert to m/d
    K_wet = aq_chars$B_K2_wet[this_ws_props] * (60 * 60 * 24) / 100 # convert to m/d
    Wo = aq_chars$TOT_STREAM_LENGTH[this_ws_props]  * 1000	# stream length in meters
    i_slp =	atan(aq_chars$TOT_BASIN_SLOPE[this_ws_props] / 100) 					# this is i in brut's equation, but I'll probably use i for something else later
    p =	1																# constant / fitting parameter; set to 1 since that what I did when I originally solved
    f = 0.1																# set to 0.1 since that's what I did when originally solved
    A = aq_chars$TOT_BASIN_AREA[this_ws_props] * 1000 *1000								# basin area in square meters
#    B = (A / (2*Wo)) / cos(i_slp)										# aquifer breadth
    B = (A / (2*Wo)) #/ cos(i_slp)										# aquifer breadth
    sin_slp = sin(i_slp)
    strm_crs_sxn = 2*Wo/A   # accounting for streams having to banks
    f = 0.1 # the drainable porosoity we set for recession analysis
    B_half = B/2
    
    
    plot(x=c(D_dry, D_wet), y=c(log10(K_dry), log10(K_wet)))
    con_mod = lm(c(log10(K_dry), log10(K_wet)) ~ 	c(D_dry, D_wet))
#    con_mod = lm(c(K_dry, K_wet) ~ 	c(D_dry, D_wet))
    
    # funscion for variable conductivity
    k_fof_D = function(D_thickness) {
#      max(c(con_mod$coef[1] + con_mod$coef[2] * D_thickness, 0.01))	#why did I put this max here?? oh yeah, the linear model sometimes goes negative
     10^(con_mod$coef[1] + con_mod$coef[2] * D_thickness)	# now trying a log mod??
    }
  
    
    
    #########################################################################################
    # initializing data for HBV
    
    #!!! TEMPFIX
    #mopey_df = tsMOPEX(mopey_catalogue[this_ws, "USGS_ID"])	# identifying the first ts data from the first idnetified mopex gage
    this_ws = ws_ihave[which_ws]
    #mopey_df = tsMOPEX(mopey_catalogue[this_ws, "USGS_ID"])	# identifying the first ts data from the first idnetified mopex gage
    #this_ws_props = which(aq_chars$gage == mopey_catalogue[this_ws, "USGS_ID"])
    idiotic_file = readLines(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\Daily\\", this_ws, ".dly"))
    less_idiotic_file = substring(idiotic_file, 9) # whoever saved the file in this format should be shot... 
    write(less_idiotic_file, paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", this_ws, ".txt"))
    mopey_df = fread(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", this_ws, ".txt"))
    names(mopey_df) = c('P','E','Q','T_max','T_min')
    for(aa in c(1,2,4,5)){
      if(any(mopey_df[ , ..aa] < -90))
        mopey_df[which(mopey_df[ , ..aa] < -90), aa] = 0
    }
    if(any(mopey_df$Q < -90)){ mopey_df$Q[which(mopey_df$Q < -90)] = NA}
    
    #!!! TEMPFIX
    mopey_df$T_mean = apply(mopey_df[,c("T_max","T_min")],1,mean)
    
    
    n_runs = 5000
    
    
    cal_out = data.frame(
      sfcf=rep(NA,n_runs),
      tr=rep(NA,n_runs), 
      tt=rep(NA,n_runs),
      fm=rep(NA,n_runs),
      fi=rep(NA,n_runs),
      fic=rep(NA,n_runs),
      fc=rep(NA,n_runs),
      lp=rep(NA,n_runs),
      beta_soils = rep(NA,n_runs),
      k0 = rep(NA,n_runs),
      k1 = rep(NA,n_runs),
      k2 = rep(NA,n_runs),
      uz1 = rep(NA,n_runs),
      perc = rep(NA,n_runs),
      hbv_nse = rep(NA,n_runs),
      new_nse = rep(NA,n_runs)
    )
    
    # Start the clock
    ptm = proc.time()
    for(jj in 1:n_runs){
        cal_out$sfcf[jj] = sample(sfcf,1)
        cal_out$tr[jj] = sample(tr,1)
        cal_out$tt[jj] = sample(tt,1)
        cal_out$fm[jj] = sample(fm,1)
        cal_out$fi[jj] = sample(fi,1)
        cal_out$fic[jj] = sample(fic,1)
        cal_out$fc[jj] = sample(fc,1)
        cal_out$lp[jj] = sample(lp,1)
        cal_out$beta_soils[jj] = sample(beta_soils,1)
        # since k0>k1>k2 and uz1>perc or an error is thrown, we need a routine to ensure this is true while still allowing random sampling
        cal_out$k0[jj] = sample(k0,1)
        cal_out$k1[jj] = min(sample(k1,1), cal_out$k0[jj]*.99)
        cal_out$k2[jj] = min(sample(k2,1), cal_out$k1[jj]*.99)
        cal_out$uz1[jj] = sample(uz1,1)
        cal_out$perc[jj] = min(sample(perc,1), cal_out$uz1[jj]*.99)
        
        mopey_df_ppt = cbind(mopey_df,	
                             SnowGlacier_HBV(model = 1,
                                             inputData = cbind(mopey_df$T_mean, mopey_df$P),
                                             initCond = c(0,2,0),	#SWE0, surface type (2=soil), surface area of glaciers as ratio [-]
                                             param = c(	cal_out$sfcf[jj],		#SFCF - snowfall correction factor [-]
                                                        cal_out$tr[jj],		#Tr - solid and liquid precipitation threshold temperature [C]
                                                        cal_out$tt[jj],		#Tt - melt temperature [C]
                                                        cal_out$fm[jj],		#fm - snowmelt factor [mm/C]
                                                        cal_out$fi[jj],		#fi - icemelt factor [mm/C]
                                                        cal_out$fic[jj])	 	#fic - debris-covered icemelt factor [mm/C]
                             )	)
        
        # soils model
        mopey_df_rech = cbind(mopey_df_ppt,
                              Soil_HBV(
                                model = 1,
                                inputData = cbind(mopey_df_ppt$Total, mopey_df_ppt$E),
                                initCond = c(50,1),	# initial soil moisture in mm, then relative ratio of soil over teh whole the whole basin
                                param = c(	cal_out$fc[jj],			#FC - soil field capacity [mm]
                                           cal_out$lp[jj],			#LP - parameter ot get actual ET [-]
                                           cal_out$beta_soils[jj])	#beta - exponential value for nonlinear relations between soil storage and runoff
                              )	)
        
        
        mopey_df_disch = cbind(mopey_df_rech,
                               Routing_HBV(
                                 model = 1,	# model=1 gives three stores with routing for each
                                 lake = FALSE,
                                 inputData = cbind(mopey_df_rech$Rech),	# recharge time series
                                 initCond = c(10,10,10),	# initial storage in each reservoir in mm
                                 param = c(	cal_out$k0[jj],	#KO - top bucket (STZ) storage constant [1/t]
                                            cal_out$k1[jj],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
                                            cal_out$k2[jj],	#K2 - lower bucket (SLZ) storage constant [1/t]
                                            cal_out$uz1[jj],	#UZL - max flux rate between STZ and SUZ [mm/t] 
                                            cal_out$perc[jj])	#PERC - max flux rate between SUZ and SLZ [mm/t]
                               )	)
        
        
        ####################################################################
        # new gw module introduced
        # calculating static values so I don't have to recalculate at every loop
        gw_out = NULL
        deep_store = D_dry * f * 1000 #as.numeric(mopey_df_disch$SLZ[1])
        depth = D_dry
        deep_rech = as.numeric(c(0,diff(mopey_df_disch$SLZ)) + mopey_df_disch$Q2)
        min_ridge = 1 / (B_half/2) # setting the minimum "ridge" i.e. dx in dh/dx to xx meter
        #        recent_rech = seq(1/182,1,length=182) / 100 # * drainable porosity and converted to meters
        accum_rech = 0
        riparian_fraction = 0#15 * (strm_crs_sxn * 2)
        evap_gw = riparian_fraction * (mopey_df_rech$E - mopey_df_rech$Eac) 
        for(gg in 1:(nrow(mopey_df_disch)))	{
          #          if(gg > length(recent_rech))	{
          #            accum_rech = sum(recent_rech * deep_rech[(gg - length(recent_rech) + 1):gg]) 
          #          }
          accum_rech = min(accum_rech * 0.95 + max(deep_rech[gg] / (f * 1000), 0), D_wet - depth / 100)
          k_arb = as.numeric(k_fof_D(D_thickness = depth))
          gw_ridge_calc = ifelse(depth < D_wet, ((D_wet - depth) / D_wet)^3, min_ridge)
          gw_ridge = max(gw_ridge_calc, min_ridge)
          gw_out_iter =
            #           f * 				# transmissivity = f * 2 * k * dh/dx * h_outlet * Length_of_streams / Area
            k_arb *			# 
            #          (depth) / (B * gw_ridge) *	# dh / dx, without accounting for slope of landscape
            (((sin_slp * B_half * gw_ridge / 2) + depth ) / (B_half * gw_ridge / 2)) *	# dh / dx, while accounting for slope of imperm layer
            #            (sin_slp * sqrt(B_half) + depth) / sqrt(B_half) *	# dh / dx, while accounting for slope of imperm layer
            (depth / 100  + accum_rech) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
            #           (depth / 4) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
            strm_crs_sxn	* 2 *
            1000                    # converting back to mm 
          if(gw_out_iter > deep_store[gg]) {gw_out_iter = deep_store[gg]}
          
          new_deep_store = deep_store[gg] +
            deep_rech[gg] -
            gw_out_iter
          deep_store = c(deep_store, new_deep_store)
          depth = new_deep_store / (f * 1000) # converting back to m and normalizing for f
          
          gw_stream = max(gw_out_iter - evap_gw[gg], 0) 
          gw_out = c(gw_out, gw_stream)
          
        }
        print(summary(deep_store[-(1:1090)] / (f * 1000)))
        print(c(D_dry, D_wet))
        print(accum_rech)
        print(gw_ridge)
        print(k_arb)
        print(c(K_dry, K_wet))
        print(riparian_fraction)
        mopey_df_disch$deep_store = deep_store[-length(deep_store)]
        mopey_df_disch$deep_rech = deep_rech
        mopey_df_disch$gw_out = gw_out
        
        ##############################################################################
        # plotting and calibrating
        
        
        these_days = which(!is.na(mopey_df_disch$Q))[1095:2190]
#        old_nse = NSE(mopey_df_disch$Qg[-(1:1095)], mopey_df_disch$Q[-(1:1095)])
#        new_nse = NSE((mopey_df_disch$Q0[-(1:1095)] + mopey_df_disch$Q1[-(1:1095)] + 	mopey_df_disch$gw_out[-(1:1095)]), 
#                      mopey_df_disch$Q[-(1:1095)])
        old_nse = KGE(mopey_df_disch$Qg[-(1:1095)], mopey_df_disch$Q[-(1:1095)])
        new_nse = KGE((mopey_df_disch$Q0[-(1:1095)] + mopey_df_disch$Q1[-(1:1095)] + 	mopey_df_disch$gw_out[-(1:1095)]), 
                      mopey_df_disch$Q[-(1:1095)])
        plot(mopey_df_disch$Q[these_days], lwd=2, type='l',
#             ylim=c(0,2),
             ylim=c(0,max(mopey_df_disch$Q[these_days]/2, na.rm=TRUE)),
#             log='y',
             ylab = "Q (mm)", xlab = "Date",
             #!!! TEMPFIX
             #		main=paste(mopey_catalogue[this_ws,"Name"],
             main=paste(this_ws,
                        #!!! TEMPFIX	           
                        ": \n old NSE =", round(old_nse,2), "new =", round(new_nse,2)))
        lines(mopey_df_disch$Qg[these_days], lwd=2, col='red', lty=1)
        lines(mopey_df_disch$Q0[these_days] + mopey_df_disch$Q1[these_days] + mopey_df_disch$gw_out[these_days],
              col='blue', lwd=2)
        lines(gw_out[these_days],col='purple1',lwd=2)
        
        # saving objective function data
        cal_out$hbv_nse[jj] = old_nse
        cal_out$new_nse[jj] = new_nse
      # calibrating by recommendation on https://www.smhi.se/en/research/research-departments/hydrology/hbv-1.90007
      #  actual_nonas = which(!is.na(mopey_df_disch$Q))[-(1:1090)]
      #  actual_flow = sum(mopey_df_disch$Q[actual_nonas])
      #  old_flow = sum(mopey_df_disch$Qg[actual_nonas])
      #  new_flow = sum(mopey_df_disch$Q0[actual_nonas]) + sum(mopey_df_disch$Q1[actual_nonas]) + 	sum(mopey_df_disch$gw_out[actual_nonas])
      #  old_nse_p = old_nse - .1 * abs(old_flow - actual_flow) / actual_flow
      #  new_nse_p = new_nse - .1 * abs(new_flow - actual_flow) / actual_flow
      #  cal_out$hbv_nse_p[jj] = old_nse_p
      #  cal_out$new_nse_p[jj] = new_nse_p
        
    }
      
    hbv_best = which.max(cal_out$hbv_nse)
    new_best = which.max(cal_out$new_nse)
    #hbv_best = which.max(cal_out$hbv_nse_p)
    #new_best = which.max(cal_out$new_nse_p)
    
    model_perf$hbv_nse[allws] = cal_out$hbv_nse[hbv_best]
    model_perf$new_nse[allws] = cal_out$new_nse[new_best]
       # hbv modules
    model_perf$hbv_sfcf[allws] = cal_out$sfcf[hbv_best]
    model_perf$hbv_tr[allws] = cal_out$tr[hbv_best]
    model_perf$hbv_tt[allws] = cal_out$tt[hbv_best]
    model_perf$hbv_fm[allws] = cal_out$fm[hbv_best]
    model_perf$hbv_fi[allws] = cal_out$fi[hbv_best]
    model_perf$hbv_fic[allws] = cal_out$fic[hbv_best]
    # soil module
    model_perf$hbv_fc[allws] = cal_out$fc[hbv_best]
    model_perf$hbv_lp[allws] = cal_out$lp[hbv_best]
    model_perf$hbv_beta_soils[allws] = cal_out$beta_soils[hbv_best]
    # routing module
    model_perf$hbv_k0[allws] = cal_out$k0[hbv_best]
    model_perf$hbv_k1[allws] = cal_out$k1[hbv_best]
    model_perf$hbv_k2[allws] = cal_out$k2[hbv_best]
    model_perf$hbv_uz1[allws] = cal_out$uz1[hbv_best]
    model_perf$hbv_perc[allws] = cal_out$perc[hbv_best]
    # hbv modules
    model_perf$new_sfcf[allws] = cal_out$sfcf[new_best]
    model_perf$new_tr[allws] = cal_out$tr[new_best]
    model_perf$new_tt[allws] = cal_out$tt[new_best]
    model_perf$new_fm[allws] = cal_out$fm[new_best]
    model_perf$new_fi[allws] = cal_out$fi[new_best]
    model_perf$new_fic[allws] = cal_out$fic[new_best]
    # soil module
    model_perf$new_fc[allws] = cal_out$fc[new_best]
    model_perf$new_lp[allws] = cal_out$lp[new_best]
    model_perf$new_beta_soils[allws] = cal_out$beta_soils[new_best]
    # routing module
    model_perf$new_k0[allws] = cal_out$k0[new_best]
    model_perf$new_k1[allws] = cal_out$k1[new_best]
    model_perf$new_k2[allws] = cal_out$k2[new_best]
    model_perf$new_uz1[allws] = cal_out$uz1[new_best]
    model_perf$new_perc[allws] = cal_out$perc[new_best]

  write.csv(model_perf, "C:\\Users\\arik\\Documents\\LSM Research\\model_perf_26apr2021_exp_2x_100div_3up_95_b_hlfB_0evap_KGE.csv")
  print(allws)
}
}
    
#model_perf = fread("C:\\Users\\arik\\Documents\\LSM Research\\model_perf_26apr2021_exp_2x_100div_3up_95_b_hlfB_0evap_KGE.csv")
  
  

  
  
  
  
  
  
################################################################
##### step 2: rerunning each hbv_ and new_ optimized model for each watershed to capture model performance
  
  
  
  for(allws in c(45:74,76:length(ws_ihave))){
    
    #################################################################################
    # initializing data for our gw module
    
    #!!! TEMPFIX
    which_ws = allws
    #this_ws = ws_ihave[which_ws]				# identifying the ws in mopex
    #this_ws_props = which(aq_chars$gage == mopey_catalogue[this_ws, "USGS_ID"])
    this_ws_props = which(aq_chars$gage == ws_ihave[which_ws])
    model_perf$gage_id[allws] = ws_ihave[which_ws]
    #!!! TEMPFIX
    
    D_dry = aq_chars$B_depth_dry[this_ws_props]  # in m?
    D_wet = aq_chars$B_depth_wet[this_ws_props] # in m?
    K_dry = aq_chars$B_K2_dry[this_ws_props] * (60 * 60 * 24) / 100 # convert to m/d
    K_wet = aq_chars$B_K2_wet[this_ws_props] * (60 * 60 * 24) / 100 # convert to m/d
    Wo = aq_chars$TOT_STREAM_LENGTH[this_ws_props]  * 1000	# stream length in meters
    i_slp =	atan(aq_chars$TOT_BASIN_SLOPE[this_ws_props] / 100) 					# this is i in brut's equation, but I'll probably use i for something else later
    p =	1																# constant / fitting parameter; set to 1 since that what I did when I originally solved
    f = 0.1																# set to 0.1 since that's what I did when originally solved
    A = aq_chars$TOT_BASIN_AREA[this_ws_props] * 1000 *1000								# basin area in square meters
    #    B = (A / (2*Wo)) / cos(i_slp)										# aquifer breadth
    B = (A / (2*Wo)) #/ cos(i_slp)										# aquifer breadth
    sin_slp = sin(i_slp)
    strm_crs_sxn = 2*Wo/A   # accounting for streams having to banks
    f = 0.1 # the drainable porosoity we set for recession analysis
    B_half = B/2
    
    
#    plot(x=c(D_dry, D_wet), y=c(log10(K_dry), log10(K_wet)))
    con_mod = lm(c(log10(K_dry), log10(K_wet)) ~ 	c(D_dry, D_wet))
#    con_mod = lm(c(K_dry, K_wet) ~ 	c(D_dry, D_wet))
    
    # funscion for variable conductivity
    k_fof_D = function(D_thickness) {
#      max(c(con_mod$coef[1] + con_mod$coef[2] * D_thickness, 0.001))	#why did I put this max here?? oh yeah, the linear model sometimes goes negative
      10^(con_mod$coef[1] + con_mod$coef[2] * D_thickness)	# now trying a log mod??
    }

    
    
    #########################################################################################
    # initializing data for HBV
    
    #!!! TEMPFIX
    #mopey_df = tsMOPEX(mopey_catalogue[this_ws, "USGS_ID"])	# identifying the first ts data from the first idnetified mopex gage
    this_ws = ws_ihave[which_ws]
    #mopey_df = tsMOPEX(mopey_catalogue[this_ws, "USGS_ID"])	# identifying the first ts data from the first idnetified mopex gage
    #this_ws_props = which(aq_chars$gage == mopey_catalogue[this_ws, "USGS_ID"])
    idiotic_file = readLines(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\Daily\\", this_ws, ".dly"))
    less_idiotic_file = substring(idiotic_file, 9) # whoever saved the file in this format should be shot... 
    write(less_idiotic_file, paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", this_ws, ".txt"))
    mopey_df = fread(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", this_ws, ".txt"))
    names(mopey_df) = c('P','E','Q','T_max','T_min')
    for(aa in c(1,2,4,5)){
      if(any(mopey_df[ , ..aa] < -90))
        mopey_df[which(mopey_df[ , ..aa] < -90), aa] = 0
    }
    if(any(mopey_df$Q < -90)){ mopey_df$Q[which(mopey_df$Q < -90)] = NA}
    
    #!!! TEMPFIX
    mopey_df$T_mean = apply(mopey_df[,c("T_max","T_min")],1,mean)
    
    
    n_runs = 2  # on run for each optimized model
    
    
     
    # Start the clock
    ptm = proc.time()
    for(jj in 1:n_runs){
    
      # snow-glacier model
      mopey_df_ppt = cbind(mopey_df,	
                           SnowGlacier_HBV(model = 1,
                                           inputData = cbind(mopey_df$T_mean, mopey_df$P),
                                           initCond = c(0,2,0),	#SWE0, surface type (2=soil), surface area of glaciers as ratio [-]
                                          if(jj == 1){
                                           param = c(	model_perf$hbv_sfcf[allws],		#SFCF - snowfall correction factor [-]
                                                      model_perf$hbv_tr[allws],		#Tr - solid and liquid precipitation threshold temperature [C]
                                                      model_perf$hbv_tt[allws],		#Tt - melt temperature [C]
                                                      model_perf$hbv_fm[allws],		#fm - snowmelt factor [mm/C]
                                                      model_perf$hbv_fi[allws],		#fi - icemelt factor [mm/C]
                                                      model_perf$hbv_fic[allws])	 	#fic - debris-covered icemelt factor [mm/C]
                                          } else  {
                                           param = c( model_perf$new_sfcf[allws],		#SFCF - snowfall correction factor [-]
                                                      model_perf$new_tr[allws],		#Tr - solid and liquid precipitation threshold temperature [C]
                                                      model_perf$new_tt[allws],		#Tt - melt temperature [C]
                                                      model_perf$new_fm[allws],		#fm - snowmelt factor [mm/C]
                                                      model_perf$new_fi[allws],		#fi - icemelt factor [mm/C]
                                                      model_perf$new_fic[allws])	 	#fic - debris-covered icemelt factor [mm/C]
                                          }
                           )	)
      
      # soils model
      mopey_df_rech = cbind(mopey_df_ppt,
                            Soil_HBV(
                              model = 1,
                              inputData = cbind(mopey_df_ppt$Total, mopey_df_ppt$E),
                              initCond = c(50,1),	# initial soil moisture in mm, then relative ratio of soil over teh whole the whole basin
                            if(jj == 1){
                              param = c( model_perf$hbv_fc[allws],			#FC - soil field capacity [mm]
                                         model_perf$hbv_lp[allws],			#LP - parameter ot get actual ET [-]
                                         model_perf$hbv_beta_soils[allws])	#beta - exponential value for nonlinear relations between soil storage and runoff
                            } else {
                              param = c( model_perf$new_fc[allws],			#FC - soil field capacity [mm]
                                         model_perf$new_lp[allws],			#LP - parameter ot get actual ET [-]
                                         model_perf$new_beta_soils[allws])	#beta - exponential value for nonlinear relations between soil storage and runoff
                              
                            }
                            )	)
      
      
      mopey_df_disch = cbind(mopey_df_rech,
                             Routing_HBV(
                               model = 1,	# model=1 gives three stores with routing for each
                               lake = FALSE,
                               inputData = cbind(mopey_df_rech$Rech),	# recharge time series
                               initCond = c(10,10,10),	# initial storage in each reservoir in mm
                             if(jj == 1){
                               param = c(	model_perf$hbv_k0[allws],	#KO - top bucket (STZ) storage constant [1/t]
                                          model_perf$hbv_k1[allws],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
                                          model_perf$hbv_k2[allws],	#K2 - lower bucket (SLZ) storage constant [1/t]
                                          model_perf$hbv_uz1[allws],	#UZL - max flux rate between STZ and SUZ [mm/t] 
                                          model_perf$hbv_perc[allws])	#PERC - max flux rate between SUZ and SLZ [mm/t]
                             } else {
                               param = c(	model_perf$new_k0[allws],	#KO - top bucket (STZ) storage constant [1/t]
                                          model_perf$new_k1[allws],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
                                         # ifelse(model_perf$hbv_k2[allws] < model_perf$new_k1[allws], model_perf$hbv_k2[allws], model_perf$new_k1[allws]-.0001),	#K2 - lower bucket (SLZ) storage constant [1/t]
                                          model_perf$new_k2[allws],	#K2 - lower bucket (SLZ) storage constant [1/t]
                                          model_perf$new_uz1[allws],	#UZL - max flux rate between STZ and SUZ [mm/t] 
                                          model_perf$new_perc[allws])	#PERC - max flux rate between SUZ and SLZ [mm/t]
                             }
                             )	)
      
      
      ####################################################################
      # new gw module introduced
      # calculating static values so I don't have to recalculate at every loop
      gw_out = NULL
      deep_store = D_dry * f * 1000 #as.numeric(mopey_df_disch$SLZ[1])
      depth = D_dry
      deep_rech = as.numeric(c(0,diff(mopey_df_disch$SLZ)) + mopey_df_disch$Q2)
      min_ridge = 1 / (B_half/2) # setting the minimum "ridge" i.e. dx in dh/dx to xx meter
      #        recent_rech = seq(1/182,1,length=182) / 100 # * drainable porosity and converted to meters
      accum_rech = 0
      riparian_fraction = 15 * (strm_crs_sxn * 2)
      evap_gw = riparian_fraction * (mopey_df_rech$E - mopey_df_rech$Eac) 
      for(gg in 1:(nrow(mopey_df_disch)))	{
        #          if(gg > length(recent_rech))	{
        #            accum_rech = sum(recent_rech * deep_rech[(gg - length(recent_rech) + 1):gg]) 
        #          }
        accum_rech = min(accum_rech * 0.95 + max(deep_rech[gg] / (f * 1000), 0), D_wet - depth / 100)
        k_arb = as.numeric(k_fof_D(D_thickness = depth))
        gw_ridge_calc = ifelse(depth < D_wet, ((D_wet - depth) / D_wet)^3, min_ridge)
        gw_ridge = max(gw_ridge_calc, min_ridge)
        gw_out_iter =
          #           f * 				# transmissivity = f * 2 * k * dh/dx * h_outlet * Length_of_streams / Area
          k_arb *			# 
          #          (depth) / (B * gw_ridge) *	# dh / dx, without accounting for slope of landscape
          (((sin_slp * B_half * gw_ridge / 2) + depth ) / (B_half * gw_ridge / 2)) *	# dh / dx, while accounting for slope of imperm layer
          #            (sin_slp * sqrt(B_half) + depth) / sqrt(B_half) *	# dh / dx, while accounting for slope of imperm layer
          (depth / 100  + accum_rech) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
          #           (depth / 4) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
          strm_crs_sxn	* 2 *
          1000                    # converting back to mm 
        if(gw_out_iter > deep_store[gg]) {gw_out_iter = deep_store[gg]}
        
        new_deep_store = deep_store[gg] +
          deep_rech[gg] -
          gw_out_iter
        deep_store = c(deep_store, new_deep_store)
        depth = new_deep_store / (f * 1000) # converting back to m and normalizing for f
        
        gw_stream = max(gw_out_iter - evap_gw[gg], 0) 
        gw_out = c(gw_out, gw_stream)
        
      }
      print(summary(deep_store[-(1:1090)] / (f * 1000)))
      print(c(D_dry, D_wet))
      print(accum_rech)
      print(gw_ridge)
      print(k_arb)
      print(c(K_dry, K_wet))
      print(riparian_fraction)
      mopey_df_disch$deep_store = deep_store[-length(deep_store)]
      mopey_df_disch$deep_rech = deep_rech
      mopey_df_disch$gw_out = gw_out
      
      ##############################################################################
      # plotting and calibrating
      
      
      these_days = which(!is.na(mopey_df_disch$Q))[2000:3090] # for illustrative plotting
      Q5_days = which(mopey_df_disch$Q <= quantile(mopey_df_disch$Q, .05, na.rm=TRUE))
      Q10_days = which(mopey_df_disch$Q <= quantile(mopey_df_disch$Q, .10, na.rm=TRUE))
      Q20_days = which(mopey_df_disch$Q <= quantile(mopey_df_disch$Q, .20, na.rm=TRUE))
      
      
      if(jj == 1){
        old_Q = mopey_df_disch$Qg
        act_Q = mopey_df_disch$Q + 10^(-5)
        old_nse = NSE(old_Q[-(1:1090)], act_Q[-(1:1090)])
        old_kge = KGE(old_Q[-(1:1090)], act_Q[-(1:1090)])
        #        old_nse_log = NSE(mopey_df_disch$Qg[-(1:365)], mopey_df_disch$Q[-(1:365)], FUN=log)
        old_rmse = rmse(old_Q[-(1:1090)], act_Q[-(1:1090)])
        old_mae = mae(old_Q[-(1:1090)], act_Q[-(1:1090)])
        Q5_old_nse = NSE(old_Q[Q5_days], act_Q[Q5_days])
        Q5_old_kge = KGE(old_Q[Q5_days], act_Q[Q5_days])
        Q5_old_rmse = rmse(old_Q[Q5_days], act_Q[Q5_days])
        Q5_old_mae = mae(old_Q[Q5_days], act_Q[Q5_days])
        Q5_old_pctB = pbias(old_Q[Q5_days], act_Q[Q5_days])
        Q10_old_nse = NSE(old_Q[Q10_days], act_Q[Q10_days])
        Q10_old_kge = KGE(old_Q[Q10_days], act_Q[Q10_days])
        Q10_old_rmse = rmse(old_Q[Q10_days], act_Q[Q10_days])
        Q10_old_mae = mae(old_Q[Q10_days], act_Q[Q10_days])
        Q10_old_pctB = pbias(old_Q[Q10_days], act_Q[Q10_days])  
        Q20_old_nse = NSE(old_Q[Q20_days], act_Q[Q20_days])
        Q20_old_kge = KGE(old_Q[Q20_days], act_Q[Q20_days])
        Q20_old_rmse = rmse(old_Q[Q20_days], act_Q[Q20_days])
        Q20_old_mae = mae(old_Q[Q20_days], act_Q[Q20_days])
        Q20_old_pctB = pbias(old_Q[Q20_days], act_Q[Q20_days])

        model_perf[allws,'hbv_nse'] = old_nse
        model_perf[allws,'hbv_kge'] = old_kge
        #        model_perf[allws,'hbv_log_nse'] = old_nse_log
        model_perf[allws,'hbv_rmse'] = old_rmse
        model_perf[allws,'hbv_mae'] = old_mae
        model_perf[allws,'hbv_Q5_nse'] = Q5_old_nse
        model_perf[allws,'hbv_Q5_kge'] = Q5_old_kge
        model_perf[allws,'hbv_Q5_rmse'] = Q5_old_rmse
        model_perf[allws,'hbv_Q5_mae'] = Q5_old_mae
        model_perf[allws,'hbv_Q5_pctB'] = Q5_old_pctB
        model_perf[allws,'hbv_Q10_nse'] = Q10_old_nse
        model_perf[allws,'hbv_Q10_kge'] = Q10_old_kge
        model_perf[allws,'hbv_Q10_rmse'] = Q10_old_rmse
        model_perf[allws,'hbv_Q10_mae'] = Q10_old_mae
        model_perf[allws,'hbv_Q10_pctB'] = Q10_old_pctB
        model_perf[allws,'hbv_Q20_nse'] = Q20_old_nse
        model_perf[allws,'hbv_Q20_kge'] = Q20_old_kge
        model_perf[allws,'hbv_Q20_rmse'] = Q20_old_rmse
        model_perf[allws,'hbv_Q20_mae'] = Q20_old_mae
        model_perf[allws,'hbv_Q20_pctB'] = Q20_old_pctB
      } else {
        new_Q = mopey_df_disch$Q0 + mopey_df_disch$Q1 + mopey_df_disch$gw_out
        act_Q = mopey_df_disch$Q + 10^(-5)
        
        new_nse = NSE(new_Q[-(1:1090)], act_Q[-(1:1090)])
        new_kge = KGE(new_Q[-(1:1090)], act_Q[-(1:1090)])
        #        new_nse_log = NSE((mopey_df_disch$Q0[-(1:365)] + mopey_df_disch$Q1[-(1:365)] + 	mopey_df_disch$gw_out[-(1:365)]), 
        #                      mopey_df_disch$Q[-(1:365)], FUN=log)
        new_rmse = rmse(new_Q[-(1:1090)], act_Q[-(1:1090)])
        new_mae = mae(new_Q[-(1:1090)], act_Q[-(1:1090)])
        Q5_new_nse = NSE(new_Q[Q5_days], act_Q[Q5_days])
        Q5_new_kge = KGE(new_Q[Q5_days], act_Q[Q5_days])
        Q5_new_rmse = rmse(new_Q[Q5_days], act_Q[Q5_days])
        Q5_new_mae = mae(new_Q[Q5_days], act_Q[Q5_days])
        Q5_new_pctB = pbias(new_Q[Q5_days], act_Q[Q5_days])
        Q10_new_nse = NSE(new_Q[Q10_days], act_Q[Q10_days])
        Q10_new_kge = KGE(new_Q[Q10_days], act_Q[Q10_days])
        Q10_new_rmse = rmse(new_Q[Q10_days], act_Q[Q10_days])
        Q10_new_mae = mae(new_Q[Q10_days], act_Q[Q10_days])
        Q10_new_pctB = pbias(new_Q[Q10_days], act_Q[Q10_days])
        Q20_new_nse = NSE(new_Q[Q20_days], act_Q[Q20_days])
        Q20_new_kge = KGE(new_Q[Q20_days], act_Q[Q20_days])
        Q20_new_rmse = rmse(new_Q[Q20_days], act_Q[Q20_days])
        Q20_new_mae = mae(new_Q[Q20_days], act_Q[Q20_days])
        Q20_new_pctB = pbias(new_Q[Q20_days], act_Q[Q20_days])
        
        model_perf[allws,'new_nse'] = new_nse
        model_perf[allws,'new_kge'] = new_kge
        #        model_perf[allws,'new_log_nse'] = new_nse_log
        model_perf[allws,'new_rmse'] = new_rmse
        model_perf[allws,'new_mae'] = new_mae
        model_perf[allws,'new_Q5_nse'] = Q5_new_nse
        model_perf[allws,'new_Q5_kge'] = Q5_new_kge
        model_perf[allws,'new_Q5_rmse'] = Q5_new_rmse
        model_perf[allws,'new_Q5_mae'] = Q5_new_mae
        model_perf[allws,'new_Q5_pctB'] = Q5_new_pctB
        model_perf[allws,'new_Q10_nse'] = Q10_new_nse
        model_perf[allws,'new_Q10_kge'] = Q10_new_kge
        model_perf[allws,'new_Q10_rmse'] = Q10_new_rmse
        model_perf[allws,'new_Q10_mae'] = Q10_new_mae
        model_perf[allws,'new_Q10_pctB'] = Q10_new_pctB
        model_perf[allws,'new_Q20_nse'] = Q20_new_nse
        model_perf[allws,'new_Q20_kge'] = Q20_new_kge
        model_perf[allws,'new_Q20_rmse'] = Q20_new_rmse
        model_perf[allws,'new_Q20_mae'] = Q20_new_mae
        model_perf[allws,'new_Q20_pctB'] = Q20_new_pctB
       }
     
      plot(mopey_df_disch$Q[Q5_days], lwd=2, type='l',
           ylim=c(0,max(mopey_df_disch$Q[Q5_days], na.rm=TRUE)),
           #log='y',
           ylab = "Q (mm)", xlab = "Date",
           #!!! TEMPFIX
           #		main=paste(mopey_catalogue[this_ws,"Name"],
           main=paste(this_ws,
                      #!!! TEMPFIX	           
                      ": \n old NSE =", round(old_nse,2), "new =", round(new_nse,2)))
      lines(mopey_df_disch$Qg[Q5_days], lwd=2, col='red', lty=1)
      lines(mopey_df_disch$Q0[Q5_days] + mopey_df_disch$Q1[Q5_days] + mopey_df_disch$gw_out[Q5_days],
            col='blue', lwd=2)
      lines(gw_out[Q5_days],col='purple1',lwd=2)
      
  }
}
  
  
  
  
  
  
basic = fread("C:/Users/arik/Documents/LSM Research/model_perf_14apr2021_exp_1x_100div_3up_982_b.csv")
switchy = fread("C:/Users/arik/Documents/LSM Research/model_perf_14apr2021_exp_1x_100div_3up_982_b_switch.csv")


newy = model_perf

complexy = fread("C:/Users/arik/Documents/LSM Research/model_perf_2apr2021_all.csv")  
newy = fread("C:/Users/arik/Documents/LSM Research/model_perf_14apr2021_lin_2x_newrech.csv")
oldy = fread("C:/Users/arik/Documents/LSM Research/model_perf_14apr2021_exp_1x_100div_3up_982_b_switch.csv")


complexy$version = "v.a1"
oldy$version = "v.c1"
newy$version = "v.c2"


allmods = rbind(complexy[,-c(1)], oldy, newy)
  
  # NSE boxplots
boxplot(hbv_nse ~ version, data = allmods,
          boxwex=0.25, at=1:3 - 0.2,
        xlim = c(0.5, 3.5), ylim = c(0, 1), yaxs = "i",
        col='tomato',
        ylab = "NSE")
boxplot(new_nse ~ version, data = allmods, add=TRUE,
        boxwex=0.25, at=1:3 + 0.2,
        col='royalblue3')
legend(.41,.985,c("HBV","new"), fil=c("tomato","royalblue3"))

  # MAE Q10
boxplot(hbv_Q10_mae ~ version, data = allmods,
        boxwex=0.25, at=1:3 - 0.2,
        xlim = c(0.5, 3.5), ylim = c(0, 1), yaxs = "i",
        col='tomato',
        ylab = "Q10 MAE (mm/d)")
boxplot(new_Q10_mae ~ version, data = allmods, add=TRUE,
        boxwex=0.25, at=1:3 + 0.2,
        col='royalblue3')
legend(.41,.985,c("HBV","new"), fil=c("tomato","royalblue3"))

# RMSE Q10
boxplot(hbv_Q10_rmse ~ version, data = allmods,
        boxwex=0.25, at=1:3 - 0.2,
        xlim = c(0.5, 3.5), ylim = c(0, 1), yaxs = "i",
        col='tomato',
        ylab = "Q10 RMSE (mm/d)")
boxplot(new_Q10_rmse ~ version, data = allmods, add=TRUE,
        boxwex=0.25, at=1:3 + 0.2,
        col='royalblue3')
legend(.41,.985,c("HBV","new"), fil=c("tomato","royalblue3"))

# MAE Q5
boxplot(hbv_Q5_mae ~ version, data = allmods,
        boxwex=0.25, at=1:3 - 0.2,
        xlim = c(0.5, 3.5), ylim = c(0, 1), yaxs = "i",
        col='tomato',
        ylab = "Q5 MAE (mm/d)")
boxplot(new_Q5_mae ~ version, data = allmods, add=TRUE,
        boxwex=0.25, at=1:3 + 0.2,
        col='royalblue3')
legend(.41,.985,c("HBV","new"), fil=c("tomato","royalblue3"))

# RMSE Q5
boxplot(hbv_Q5_rmse ~ version, data = allmods,
        boxwex=0.25, at=1:3 - 0.2,
        xlim = c(0.5, 3.5), ylim = c(0, 1), yaxs = "i",
        col='tomato',
        ylab = "Q5 RMSE (mm/d)")
boxplot(new_Q5_rmse ~ version, data = allmods, add=TRUE,
        boxwex=0.25, at=1:3 + 0.2,
        col='royalblue3')
legend(.41,.985,c("HBV","new"), fil=c("tomato","royalblue3"))

  
boxplot(hbv_Q5_rmse ~ version, data = allmods,
        boxwex=0.25, at=1:3 - 0.2,
        xlim = c(0.5, 3.5), ylim = c(0, 1), yaxs = "i",
        col='tomato',
        ylab = "Q5 RMSE (mm/d)")
boxplot(new_Q5_rmse ~ version, data = allmods, add=TRUE,
        boxwex=0.25, at=1:3 + 0.2,
        col='royalblue3')

  
  
boxplot(newy$hbv_Q20_pctB, data = newy,
        boxwex=0.25, at=1 - 0.2,
        xlim = c(0.5, 1.5), ylim = c(-50, 10), yaxs = "i",
        col='tomato',
        ylab = "% bias")
boxplot(newy$new_Q20_pctB , data = newy, add=TRUE,
        boxwex=0.25, at=1 + 0.2,
        col='royalblue3')
abline(h=0)
legend(.5,.985,c("HBV","new"), fil=c("tomato","royalblue3"))

  
  ### data to save for the final run
  

  
  



































































################################################################
### model 1: the basic 3 component HBV runoff module
# snow and glacier module
  sfcf = runif(100, .5, 4)#1, 2)	#snowfall correction factor [-]
  tr   = runif(100, -2, 4)#-1, 3)	#solid and liquid precipitation threshold temperature [C]
  tt   = runif(100, -1, 4)#0, 3)	#melt temperature [C]
  fm   = runif(100, .5, 6)#1, 4)	#snowmelt factor [mm/C]
  fi   = runif(100, 2, 10)#4, 8)	#icemelt factor [mm/C]
  fic  = runif(100, 4, 10)#6, 6)	#debris-covered icemelt factor [mm/C]

  # soil module
  fc   = runif(100, 50, 400)#100, 200)
  lp   = runif(100, .2, 1)#0.5, 1)
  beta_soils = runif(100, 1, 3)

  # routing module
  k0   = runif(100, .01, .9)#0.09, 0.1)
  k1   = runif(100, .001, .5)#0.05, 0.07)
  k2   = runif(100, .0001, .1)#0.05)
  uz1  = runif(100, .1, 20)#1, 5)
  perc = runif(100, .0001, 20)#.8, 2)


this_ws = 8

  #!!!! TEMPFIX
#mopey_df = tsMOPEX(mopey_catalogue$USGS_ID[this_ws])	# identifying the first ts data from the first idnetified mopex gage
idiotic_file = readLines(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\Daily\\", mopey_catalogue[this_ws], ".dly"))
less_idiotic_file = substring(idiotic_file, 9) # whoever saved the file in this format should be shot... 
write(less_idiotic_file, paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", mopey_catalogue[this_ws], ".txt"))
mopey_df = fread(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", mopey_catalogue[this_ws], ".txt"))
names(mopey_df) = c('P','E','Q','T_max','T_min')
for(aa in 1:ncol(mopey_df)){
  if(any(mopey_df[ , ..aa] < -90))
  mopey_df[which(mopey_df[ , ..aa] < -90), aa] = NA
  }
summary(mopey_df$P)
summary(mopey_df$E)
summary(mopey_df$Q)
summary(mopey_df$T_max)
summary(mopey_df$T_min)
#!!!! TEMPFIX


#mopey_df_nonas = subset(mopey_df, !is.na(Q) & !is.na(P))
#mopey_df_nonas$T_mean = apply(mopey_df_nonas[,c("T_max","T_min")],1,mean)
mopey_df$T_mean = apply(mopey_df[,c("T_max","T_min")],1,mean)


n_loops = 10
n_runs = 1000

cal_out = data.frame(
  sfcf=rep(NA,n_runs),
  tr=rep(NA,n_runs), 
  tt=rep(NA,n_runs),
  fm=rep(NA,n_runs),
  fi=rep(NA,n_runs),
  fic=rep(NA,n_runs),
  fc=rep(NA,n_runs),
  lp=rep(NA,n_runs),
  beta_soils = rep(NA,n_runs),
  k0 = rep(NA,n_runs),
  k1 = rep(NA,n_runs),
  k2 = rep(NA,n_runs),
  uz1 = rep(NA,n_runs),
  perc = rep(NA,n_runs),
  nse = rep(NA,n_runs))


for(ii in 1:n_loops){
	for(jj in 1:n_runs){
		cal_out$sfcf[jj] = sample(sfcf,1)
		cal_out$tr[jj] = sample(tr,1)
		cal_out$tt[jj] = sample(tt,1)
		cal_out$fm[jj] = sample(fm,1)
		cal_out$fi[jj] = sample(fi,1)
		cal_out$fic[jj] = sample(fic,1)
		cal_out$fc[jj] = sample(fc,1)
		cal_out$lp[jj] = sample(lp,1)
		cal_out$beta_soils[jj] = sample(beta_soils,1)
		# since k0>k1>k2 and uz1>perc or an error is thrown, we need a routine to ensure this is true while still allowing random sampling
		cal_out$k0[jj] = sample(k0,1)
		cal_out$k1[jj] = min(sample(k1,1), cal_out$k0[jj]*.99)
		cal_out$k2[jj] = min(sample(k2,1), cal_out$k1[jj]*.99)
		cal_out$uz1[jj] = sample(uz1,1)
		cal_out$perc[jj] = min(sample(perc,1), cal_out$uz1[jj]*.99)

		mopey_df_ppt = cbind(mopey_df,	
			SnowGlacier_HBV(model = 1,
				inputData = cbind(mopey_df$T_mean, mopey_df$P),
				initCond = c(0,2,0),	#SWE0, surface type (2=soil), surface area of glaciers as ratio [-]
				param = c(	cal_out$sfcf[jj],		#SFCF - snowfall correction factor [-]
							cal_out$tr[jj],		#Tr - solid and liquid precipitation threshold temperature [C]
							cal_out$tt[jj],		#Tt - melt temperature [C]
							cal_out$fm[jj],		#fm - snowmelt factor [mm/C]
							cal_out$fi[jj],		#fi - icemelt factor [mm/C]
							cal_out$fic[jj])	 	#fic - debris-covered icemelt factor [mm/C]
		)	)

			# soils model
		mopey_df_rech = cbind(mopey_df_ppt,
			Soil_HBV(
			   model = 1,
			   inputData = cbind(mopey_df_ppt$Total, mopey_df_ppt$E),
			   initCond = c(50,1),	# initial soil moisture in mm, then relative ratio of soil over teh whole the whole basin
			   param = c(	cal_out$fc[jj],			#FC - soil field capacity [mm]
							cal_out$lp[jj],			#LP - parameter ot get actual ET [-]
							cal_out$beta_soils[jj])	#beta - exponential value for nonlinear relations between soil storage and runoff
		)	)


		mopey_df_disch = cbind(mopey_df_rech,
			Routing_HBV(
				model = 1,	# model=1 gives three stores with routing for each
				lake = FALSE,
				inputData = cbind(mopey_df_rech$Rech),	# recharge time series
				initCond = c(10,10,10),	# initial storage in each reservoir in mm
				param = c(	cal_out$k0[jj],	#KO - top bucket (STZ) storage constant [1/t]
							cal_out$k1[jj],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
							cal_out$k2[jj],	#K2 - lower bucket (SLZ) storage constant [1/t]
							cal_out$uz1[jj],	#UZL - max flux rate between STZ and SUZ [mm/t] 
							cal_out$perc[jj])	#PERC - max flux rate between SUZ and SLZ [mm/t]
		)	)


		the_nse = NSE(mopey_df_disch$Qg, mopey_df_disch$Q)
		plot(mopey_df_disch$Q[365:730], lwd=2, type='l', ylim=c(0,max(mopey_df_disch$Q[365:730],na.rm=TRUE)),
			#log='y',
			ylab = "Q (mm)", xlab = "Date",
  #!!! TEMPFFIX
#			main=paste(mopey_catalogue[this_ws,"Name"], ": NSE =", round(the_nse,2)))
      main=paste(mopey_catalogue[this_ws], ": NSE =", round(the_nse,2)))
		#!!! TEMPFFIX
    lines(mopey_df_disch$Qg[365:730], lwd=2, col='blue', lty=1)

		# calibrating by NSE
		cal_out$nse[jj][-(1:365)] = the_nse
		# calibrating by recommendation on https://www.smhi.se/en/research/research-departments/hydrology/hbv-1.90007
			# basically NSE - .1 * abs(relative volume error)
		

	}
	print(ii)
	cal_sort = head(cal_out[rev(order(cal_out$nse)),], 100)
	# snow and glacier module
	  sfcf = runif(100, min(cal_sort$sfcf), max(cal_sort$sfcf))#1, 2)	#snowfall correction factor [-]
	  tr   = runif(100, min(cal_sort$tr), max(cal_sort$tr))#-1, 3)	#solid and liquid precipitation threshold temperature [C]
	  tt   = runif(100, min(cal_sort$tt), max(cal_sort$tt))#0, 3)	#melt temperature [C]
	  fm   = runif(100, min(cal_sort$fm), max(cal_sort$fm))#1, 4)	#snowmelt factor [mm/C]
	  fi   = runif(100, min(cal_sort$fi), max(cal_sort$fi))#4, 8)	#icemelt factor [mm/C]
	  fic  = runif(100, min(cal_sort$fic), max(cal_sort$fic))#6, 6)	#debris-covered icemelt factor [mm/C]
	# soil module
	  fc   = runif(100, min(cal_sort$fc), max(cal_sort$fc))#100, 200)
	  lp   = runif(100, min(cal_sort$lp), max(cal_sort$lp))#0.5, 1)
	  beta_soils = runif(100, min(cal_sort$beta_soils), max(cal_sort$beta_soils))
	# routing module
	  k0   = runif(100, min(cal_sort$k0), max(cal_sort$k0))#0.09, 0.1)
	  k1   = runif(100, min(cal_sort$k1), max(cal_sort$k1))#0.05, 0.07)
	  k2   = runif(100, min(cal_sort$k2), max(cal_sort$k2))#0.05)
	  uz1  = runif(100, min(cal_sort$uz1), max(cal_sort$uz1))#1, 5)
	  perc = runif(100, min(cal_sort$perc), max(cal_sort$perc))#.8, 2)
}



#########################################################################################
#### model 2: adding our gw module
library(sf)
	# reading in the data calculated for all basins
setwd("C:\\Users\\arik\\Documents\\PhD Research\\D3")
aq_chars = st_read("basin_and_aq_chars_with_polygons_jan2020.gpkg")

	# choose the watersheds
  #!!! TEMPFIX
#for(gg in 1:nrow(mopey_catalogue)){
for(gg in mopey_catalogue){
  #!!! TEMPFIX
  
  this_ws = gg
    #!!! TEMPFIX
  #if(mopey_catalogue[this_ws, "USGS_ID"] %in% aq_chars$gage) 	print(this_ws)
  if(this_ws %in% aq_chars$gage) 	print(this_ws)
    #!!! TEMPFIX
}

  #!!! TEMPFIX
this_ws = "14232500"
#mopey_df = tsMOPEX(mopey_catalogue[this_ws, "USGS_ID"])	# identifying the first ts data from the first idnetified mopex gage
#this_ws_props = which(aq_chars$gage == mopey_catalogue[this_ws, "USGS_ID"])
idiotic_file = readLines(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\Daily\\", this_ws, ".dly"))
less_idiotic_file = substring(idiotic_file, 9) # whoever saved the file in this format should be shot... 
write(less_idiotic_file, paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", this_ws, ".txt"))
mopey_df = fread(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", this_ws, ".txt"))
names(mopey_df) = c('P','E','Q','T_max','T_min')
for(aa in 1:ncol(mopey_df)){
  if(any(mopey_df[ , ..aa] < -90))
    mopey_df[which(mopey_df[ , ..aa] < -90), aa] = NA
}
summary(mopey_df$P)
summary(mopey_df$E)
summary(mopey_df$Q)
summary(mopey_df$T_max)
summary(mopey_df$T_min)
this_ws_props = which(aq_chars$gage == this_ws)
  #!!! TEMPFIX


D_dry = aq_chars$B_depth_dry[this_ws_props]
D_wet = aq_chars$B_depth_wet[this_ws_props]
K_dry = aq_chars$B_K2_dry[this_ws_props]
K_wet = aq_chars$B_K2_wet[this_ws_props]
Wo = aq_chars$TOT_STREAM_LENGTH[this_ws_props]  * 1000	# stream length in meters
i_slp =	atan((1 + aq_chars$TOT_BASIN_SLOPE[this_ws_props]) / 100) 					# this is i in brut's equation, but I'll probably use i for something else later
p =	1																# constant / fitting parameter; set to 1 since that what I did when I originally solved
f = 0.1																# set to 0.1 since that's what I did when originally solved
A = aq_chars$TOT_BASIN_AREA[this_ws_props] * 1000 *1000								# basin area in square meters
B = (A / (2*Wo)) / cos(i_slp)										# aquifer breadth


plot(x=c(D_dry, D_wet), y=c(K_dry, K_wet))
con_mod = lm(c(K_dry, K_wet) ~ 	c(D_dry, D_wet))

	# funscion for variable conductivity
k_fof_D = function(D_thickness) {
	#max(c(con_mod$coef[1] + con_mod$coef[2] * D_thickness, 0.001))	#why did I put this max here??
	con_mod$coef[1] + con_mod$coef[2] * D_thickness
}


mopey_df_disch$deep_rech = c(0,diff(mopey_df_disch$SLZ)) + mopey_df_disch$Q2
f = 0.1 # the drainable porosoity we set for recession analysis
mopey_df_disch$gw_out = NA
mopey_df_disch$deep_store = NA
mopey_df_disch$deep_store[1] =  mopey_df_disch$SLZ[1]
depth = mopey_df_disch$deep_store[1] / f
B_half = B/2

	# Start the clock
	prm = proc.time()

	# calculating values so I don't have to recalculate at every loop
sin_slp = sin(i_slp)
strm_crs_sxn = 2*Wo/A
		# Start the clock!
		ptm <- proc.time()
for(gg in 1:(nrow(mopey_df_disch)-1))	{
	k_arb = k_fof_D(D_thickness = depth)
	gw_ridge = ifelse(D_wet > depth, ((D_wet - depth) / D_wet)^2, .0001)
	mopey_df_disch$gw_out[gg] = f * 				# transmissivity = f * 2 * k * dh/dx * h_outlet * Length_of_streams / Area
		k_arb *			# 2k to account for 2 hillslopes per reach of stream
#		(((sin(i_slp) * (B/2)) + 2 * depth ) / sqrt(B/2)) *	# dh / dx, while accounting for slope of imperm layer
		(((sin_slp * B_half) + 2 * depth ) / (B * gw_ridge)) *	# dh / dx, while accounting for slope of imperm layer
		(depth / 2) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
		strm_crs_sxn					# 
	mopey_df_disch$deep_store[gg+1] = mopey_df_disch$deep_store[gg] + mopey_df_disch$deep_rech[gg] - mopey_df_disch$gw_out[gg]
	depth = mopey_df_disch$deep_store[gg+1] / f
}
	# Stop the clock
	proc.time() - ptm


		# Start the clock!
		ptm <- proc.time()
f = 0.1 # the drainable porosoity we set for recession analysis
gw_out = NULL
deep_store = as.numeric(mopey_df_disch$SLZ[1])
depth = deep_store / f
deep_rech = as.numeric(c(0,diff(mopey_df_disch$SLZ)) + mopey_df_disch$Q2)
for(gg in 1:(nrow(mopey_df_disch)))	{
	k_arb = as.numeric(k_fof_D(D_thickness = depth))
	gw_ridge = ifelse(depth < D_wet, ((D_wet - depth) / D_wet)^2, .001)
	gw_out_iter = f * 				# transmissivity = f * 2 * k * dh/dx * h_outlet * Length_of_streams / Area
		k_arb *			# 2k to account for 2 hillslopes per reach of stream
#		(((sin(i_slp) * (B/2)) + 2 * depth ) / sqrt(B/2)) *	# dh / dx, while accounting for slope of imperm layer
		(((sin_slp * B_half) + 2 * depth ) / (B * gw_ridge)) *	# dh / dx, while accounting for slope of imperm layer
		(depth / 2) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
		strm_crs_sxn					# 
	
	new_deep_store = deep_store[gg] +
			deep_rech[gg] -
			gw_out_iter
	deep_store = c(deep_store, new_deep_store)
	depth = new_deep_store / f

	gw_out = c(gw_out, gw_out_iter)
	
}
	
mopey_df_disch$deep_store = deep_store[-length(deep_store)]
mopey_df_disch$deep_rech = deep_rech
mopey_df_disch$gw_out = gw_out
	# Stop the clock
	proc.time() - ptm



these_days = 2600:2965
old_nse = NSE(mopey_df_disch$Qg, mopey_df_disch$Q)
new_nse = NSE((mopey_df_disch$Q0 + mopey_df_disch$Q1 + mopey_df_disch$gw_out), mopey_df_disch$Q)
		plot(mopey_df_disch$Q[these_days], lwd=2, type='l', ylim=c(0,max(mopey_df_disch$Q[2600:2965],na.rm=TRUE)),
			#log='y',
			ylab = "Q (mm)", xlab = "Date",
  #!!! TEMPFIX
#			main=paste(mopey_catalogue[this_ws,"Name"],
      main=paste(this_ws,
  #!!! TEMPFIX
			": \n old NSE =", round(old_nse,2), "new =", round(new_nse,2)))
		lines(mopey_df_disch$Qg[these_days], lwd=2, col='red', lty=1)
		lines(mopey_df_disch$Q0[these_days] + mopey_df_disch$Q1[these_days] + mopey_df_disch$gw_out[these_days],
			col='blue', lwd=2)










#######################################################
### step 3: calibrating HBV with our model included
	# reading in the data calculated for all basins
setwd("C:\\Users\\arik\\Documents\\PhD Research\\D3")
aq_chars = st_read("basin_and_aq_chars_with_polygons_jan2020.gpkg")

	# choose the watersheds
ws_ihave = NULL
  #!!! TEMPFIX
#for(gg in 1:nrow(mopey_catalogue)){
  #if(mopey_catalogue[gg, "USGS_ID"] %in% aq_chars$gage) 	ws_ihave = c(ws_ihave, gg)
#}
#head(mopey_catalogue[ws_ihave, ])
for(gg in mopey_catalogue){
  if(gg %in% aq_chars$gage) 	ws_ihave = c(ws_ihave, gg)
}
head(ws_ihave)
  #!!! TEMPFIX




#################################################################################
	# initializing data for our gw module

  #!!! TEMPFIX
which_ws = all_ws
#this_ws = ws_ihave[which_ws]				# identifying the ws in mopex
#this_ws_props = which(aq_chars$gage == mopey_catalogue[this_ws, "USGS_ID"])
this_ws_props = which(aq_chars$gage == ws_ihave[which_ws])
  #!!! TEMPFIX

D_dry = aq_chars$B_depth_dry[this_ws_props]
D_wet = aq_chars$B_depth_wet[this_ws_props]
K_dry = aq_chars$B_K2_dry[this_ws_props]
K_wet = aq_chars$B_K2_wet[this_ws_props]
Wo = aq_chars$TOT_STREAM_LENGTH[this_ws_props]  * 1000	# stream length in meters
i_slp =	atan((1 + aq_chars$TOT_BASIN_SLOPE[this_ws_props]) / 100) 					# this is i in brut's equation, but I'll probably use i for something else later
p =	1																# constant / fitting parameter; set to 1 since that what I did when I originally solved
f = 0.1																# set to 0.1 since that's what I did when originally solved
A = aq_chars$TOT_BASIN_AREA[this_ws_props] * 1000 *1000								# basin area in square meters
B = (A / (2*Wo)) / cos(i_slp)										# aquifer breadth


plot(x=c(D_dry, D_wet), y=c(K_dry, K_wet))
con_mod = lm(c(K_dry, K_wet) ~ 	c(D_dry, D_wet))

	# funscion for variable conductivity
k_fof_D = function(D_thickness) {
	#max(c(con_mod$coef[1] + con_mod$coef[2] * D_thickness, 0.001))	#why did I put this max here??
	con_mod$coef[1] + con_mod$coef[2] * D_thickness
}



#########################################################################################
	# initializing data for HBV

  #!!! TEMPFIX
#mopey_df = tsMOPEX(mopey_catalogue[this_ws, "USGS_ID"])	# identifying the first ts data from the first idnetified mopex gage
this_ws = ws_ihave[which_ws]
#mopey_df = tsMOPEX(mopey_catalogue[this_ws, "USGS_ID"])	# identifying the first ts data from the first idnetified mopex gage
#this_ws_props = which(aq_chars$gage == mopey_catalogue[this_ws, "USGS_ID"])
idiotic_file = readLines(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\Daily\\", this_ws, ".dly"))
less_idiotic_file = substring(idiotic_file, 9) # whoever saved the file in this format should be shot... 
write(less_idiotic_file, paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", this_ws, ".txt"))
mopey_df = fread(paste0("C:\\Users\\arik\\Documents\\LSM Research\\MOPEX_data\\", this_ws, ".txt"))
names(mopey_df) = c('P','E','Q','T_max','T_min')
for(aa in 1:ncol(mopey_df)){
  if(any(mopey_df[ , ..aa] < -90))
    mopey_df[which(mopey_df[ , ..aa] < -90), aa] = NA
}
summary(mopey_df$P)
summary(mopey_df$E)
summary(mopey_df$Q)
summary(mopey_df$T_max)
summary(mopey_df$T_min)
#!!! TEMPFIX




################################################################
### model 1: the basic 3 component HBV runoff module
# snow and glacier module
  sfcf = runif(100, .3, 5)#1, 2)	#snowfall correction factor [-]
  tr   = runif(100, -2, 5)#-1, 3)	#solid and liquid precipitation threshold temperature [C]
  tt   = runif(100, -2, 5)#0, 3)	#melt temperature [C]
  fm   = runif(100, .2, 8)#1, 4)	#snowmelt factor [mm/C]
  fi   = runif(100, 2, 10)#4, 8)	#icemelt factor [mm/C]
  fic  = runif(100, 4, 10)#6, 6)	#debris-covered icemelt factor [mm/C]

  # soil module
  fc   = runif(100, 50, 2000)#100, 200)
  lp   = runif(100, .1, 1)#0.5, 1)
  beta_soils = runif(100, 1, 3)

  # routing module
  k0   = runif(100, .01, 1)#0.09, 0.1)
  k1   = runif(100, .0001, .5)#0.05, 0.07)
  k2   = 0.01#runif(100, .0001, .1)#0.05)	# since we won't be using this component of hbv
  uz1  = runif(100, .1, 20)#1, 5)
  perc = runif(100, .0001, 20)#.8, 2)




mopey_df$T_mean = apply(mopey_df[,c("T_max","T_min")],1,mean)


n_runs = 100
n_loops = 10


cal_out = data.frame(
	sfcf=rep(NA,n_runs),
	tr=rep(NA,n_runs), 
	tt=rep(NA,n_runs),
	fm=rep(NA,n_runs),
	fi=rep(NA,n_runs),
	fic=rep(NA,n_runs),
	fc=rep(NA,n_runs),
	lp=rep(NA,n_runs),
	beta_soils = rep(NA,n_runs),
	k0 = rep(NA,n_runs),
	k1 = rep(NA,n_runs),
	k2 = rep(NA,n_runs),
	uz1 = rep(NA,n_runs),
	perc = rep(NA,n_runs),
	nse = rep(NA,n_runs))

	# Start the clock
ptm = proc.time()
for(ii in 1:n_loops){
	for(jj in 1:n_runs){
		cal_out$sfcf[jj] = sample(sfcf,1)
		cal_out$tr[jj] = sample(tr,1)
		cal_out$tt[jj] = sample(tt,1)
		cal_out$fm[jj] = sample(fm,1)
		cal_out$fi[jj] = sample(fi,1)
		cal_out$fic[jj] = sample(fic,1)
		cal_out$fc[jj] = sample(fc,1)
		cal_out$lp[jj] = sample(lp,1)
		cal_out$beta_soils[jj] = sample(beta_soils,1)
		# since k0>k1>k2 and uz1>perc or an error is thrown, we need a routine to ensure this is true while still allowing random sampling
		cal_out$k0[jj] = sample(k0,1)
		cal_out$k1[jj] = min(sample(k1,1), cal_out$k0[jj]*.99)
		cal_out$k2[jj] = min(sample(k2,1), cal_out$k1[jj]*.99)
		cal_out$uz1[jj] = sample(uz1,1)
		cal_out$perc[jj] = min(sample(perc,1), cal_out$uz1[jj]*.99)

		mopey_df_ppt = cbind(mopey_df,	
			SnowGlacier_HBV(model = 1,
				inputData = cbind(mopey_df$T_mean, mopey_df$P),
				initCond = c(0,2,0),	#SWE0, surface type (2=soil), surface area of glaciers as ratio [-]
				param = c(	cal_out$sfcf[jj],		#SFCF - snowfall correction factor [-]
							cal_out$tr[jj],		#Tr - solid and liquid precipitation threshold temperature [C]
							cal_out$tt[jj],		#Tt - melt temperature [C]
							cal_out$fm[jj],		#fm - snowmelt factor [mm/C]
							cal_out$fi[jj],		#fi - icemelt factor [mm/C]
							cal_out$fic[jj])	 	#fic - debris-covered icemelt factor [mm/C]
		)	)

			# soils model
		mopey_df_rech = cbind(mopey_df_ppt,
			Soil_HBV(
			   model = 1,
			   inputData = cbind(mopey_df_ppt$Total, mopey_df_ppt$E),
			   initCond = c(50,1),	# initial soil moisture in mm, then relative ratio of soil over teh whole the whole basin
			   param = c(	cal_out$fc[jj],			#FC - soil field capacity [mm]
							cal_out$lp[jj],			#LP - parameter ot get actual ET [-]
							cal_out$beta_soils[jj])	#beta - exponential value for nonlinear relations between soil storage and runoff
		)	)


		mopey_df_disch = cbind(mopey_df_rech,
			Routing_HBV(
				model = 1,	# model=1 gives three stores with routing for each
				lake = FALSE,
				inputData = cbind(mopey_df_rech$Rech),	# recharge time series
				initCond = c(10,10,10),	# initial storage in each reservoir in mm
				param = c(	cal_out$k0[jj],	#KO - top bucket (STZ) storage constant [1/t]
							cal_out$k1[jj],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
							cal_out$k2[jj],	#K2 - lower bucket (SLZ) storage constant [1/t]
							cal_out$uz1[jj],	#UZL - max flux rate between STZ and SUZ [mm/t] 
							cal_out$perc[jj])	#PERC - max flux rate between SUZ and SLZ [mm/t]
		)	)


####################################################################
	# new gw module introduced

			f = 0.1 # the drainable porosoity we set for recession analysis
			gw_out = NULL
			deep_store = as.numeric(mopey_df_disch$SLZ[1])
			depth = deep_store / f
			deep_rech = as.numeric(c(0,diff(mopey_df_disch$SLZ)) + mopey_df_disch$Q2)
			recent_rech = seq(1/180,1,length=180)^2
			accum_rech = 0
			for(gg in 1:(nrow(mopey_df_disch)))	{
				if(gg > length(recent_rech))	{
					accum_rech = sum(recent_rech * deep_rech[(gg - length(recent_rech) + 1):gg]) / f
				}
				k_arb = as.numeric(k_fof_D(D_thickness = depth))
				gw_ridge = ifelse(depth < D_wet, ((D_wet - depth) / D_wet)^2, .002)
				gw_out_iter =
				  f * 				# transmissivity = f * 2 * k * dh/dx * h_outlet * Length_of_streams / Area
					k_arb *			# 2k to account for 2 hillslopes per reach of stream
			#		(((sin(i_slp) * (B/2)) + 2 * depth ) / sqrt(B/2)) *	# dh / dx, while accounting for slope of imperm layer
					(((sin_slp * B_half * gw_ridge) + depth ) / (B * gw_ridge)) *	# dh / dx, while accounting for slope of imperm layer
					(depth / 4 + accum_rech) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
					strm_crs_sxn					# 
				if(gw_out_iter > deep_store[gg]) {gw_out_iter = deep_store[gg]}
				
				new_deep_store = deep_store[gg] +
						deep_rech[gg] -
						gw_out_iter
				deep_store = c(deep_store, new_deep_store)
				depth = new_deep_store / f

				gw_out = c(gw_out, gw_out_iter)
				
			}
				
			mopey_df_disch$deep_store = deep_store[-length(deep_store)]
			mopey_df_disch$deep_rech = deep_rech
			mopey_df_disch$gw_out = gw_out

##############################################################################
	# plotting and calibrating
	
	
	these_days = 3300:3665
	old_nse = NSE(mopey_df_disch$Qg[-(1:365)], mopey_df_disch$Q[-(1:365)])
	new_nse = NSE((mopey_df_disch$Q0[-(1:365)] + mopey_df_disch$Q1[-(1:365)] + 	mopey_df_disch$gw_out[-(1:365)]), 
		mopey_df_disch$Q[-(1:365)])
	plot(mopey_df_disch$Q[these_days], lwd=2, type='l',
		ylim=c(0,max(mopey_df_disch$Q[these_days], na.rm=TRUE)),
		#log='y',
		ylab = "Q (mm)", xlab = "Date",
  #!!! TEMPFIX
#		main=paste(mopey_catalogue[this_ws,"Name"],
	  main=paste(this_ws,
	#!!! TEMPFIX	           
		": \n old NSE =", round(old_nse,2), "new =", round(new_nse,2)))
	lines(mopey_df_disch$Qg[these_days], lwd=2, col='red', lty=1)
	lines(mopey_df_disch$Q0[these_days] + mopey_df_disch$Q1[these_days] + mopey_df_disch$gw_out[these_days],
		col='blue', lwd=2)

	# calibrating by NSE
	cal_out$nse[jj] = new_nse
	# calibrating by recommendation on https://www.smhi.se/en/research/research-departments/hydrology/hbv-1.90007
		# basically NSE - .1 * abs(relative volume error)

	}
  
	print(ii)
	cal_sort = head(cal_out[rev(order(cal_out$nse)),], 10)
	print(head(cal_sort, 20))
	# snow and glacier module
	  sfcf = runif(100, min(cal_sort$sfcf), max(cal_sort$sfcf))#1, 2)	#snowfall correction factor [-]
	  tr   = runif(100, min(cal_sort$tr), max(cal_sort$tr))#-1, 3)	#solid and liquid precipitation threshold temperature [C]
	  tt   = runif(100, min(cal_sort$tt), max(cal_sort$tt))#0, 3)	#melt temperature [C]
	  fm   = runif(100, min(cal_sort$fm), max(cal_sort$fm))#1, 4)	#snowmelt factor [mm/C]
	  fi   = runif(100, min(cal_sort$fi), max(cal_sort$fi))#4, 8)	#icemelt factor [mm/C]
	  fic  = runif(100, min(cal_sort$fic), max(cal_sort$fic))#6, 6)	#debris-covered icemelt factor [mm/C]
	# soil module
	  fc   = runif(100, min(cal_sort$fc), max(cal_sort$fc))#100, 200)
	  lp   = runif(100, min(cal_sort$lp), max(cal_sort$lp))#0.5, 1)
	  beta_soils = runif(100, min(cal_sort$beta_soils), max(cal_sort$beta_soils))
	# routing module
	  k0   = runif(100, min(cal_sort$k0), max(cal_sort$k0))#0.09, 0.1)
	  k1   = runif(100, min(cal_sort$k1), max(cal_sort$k1))#0.05, 0.07)
	  #k2   = runif(100, min(cal_sort$k2), max(cal_sort$k2))#0.05)
	  uz1  = runif(100, min(cal_sort$uz1), max(cal_sort$uz1))#1, 5)
	  perc = runif(100, min(cal_sort$perc), max(cal_sort$perc))#.8, 2)
}
	#Stop clock
proc.time() - ptm








# now assessing variations in our new gw module

mopey_df_ppt = cbind(mopey_df,	
			SnowGlacier_HBV(model = 1,
				inputData = cbind(mopey_df$T_mean, mopey_df$P),
				initCond = c(0,2,0),	#SWE0, surface type (2=soil), surface area of glaciers as ratio [-]
				param = c(	cal_sort$sfcf[1],		#SFCF - snowfall correction factor [-]
							cal_sort$tr[1],		#Tr - solid and liquid precipitation threshold temperature [C]
							cal_sort$tt[1],		#Tt - melt temperature [C]
							cal_sort$fm[1],		#fm - snowmelt factor [mm/C]
							cal_sort$fi[1],		#fi - icemelt factor [mm/C]
							cal_sort$fic[1])	 	#fic - debris-covered icemelt factor [mm/C]
		)	)

			# soils model
		mopey_df_rech = cbind(mopey_df_ppt,
			Soil_HBV(
			   model = 1,
			   inputData = cbind(mopey_df_ppt$Total, mopey_df_ppt$E),
			   initCond = c(50,1),	# initial soil moisture in mm, then relative ratio of soil over teh whole the whole basin
			   param = c(	cal_sort$fc[1],			#FC - soil field capacity [mm]
							cal_sort$lp[1],			#LP - parameter ot get actual ET [-]
							cal_sort$beta_soils[1])	#beta - exponential value for nonlinear relations between soil storage and runoff
		)	)


		mopey_df_disch = cbind(mopey_df_rech,
			Routing_HBV(
				model = 1,	# model=1 gives three stores with routing for each
				lake = FALSE,
				inputData = cbind(mopey_df_rech$Rech),	# recharge time series
				initCond = c(10,10,10),	# initial storage in each reservoir in mm
				param = c(	cal_sort$k0[1],	#KO - top bucket (STZ) storage constant [1/t]
							cal_sort$k1[1],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
							cal_sort$k2[1],	#K2 - lower bucket (SLZ) storage constant [1/t]
							cal_sort$uz1[1],	#UZL - max flux rate between STZ and SUZ [mm/t] 
							cal_sort$perc[1])	#PERC - max flux rate between SUZ and SLZ [mm/t]
		)	)


####################################################################
	# new gw module introduced
			length_rech_dist = 20:190			# opt based on 3 runs
			exp_rech_dist = runif(100,1,3) # opt based on 2 runs: smaller is better
			depth_sclr_dist = runif(100,1,8) # opt based on 2 runs: smaller is better (1)
			ridge_exp_dist = runif(100,2,4)	# opt based on 2 runs: doesn't seem to matter
			ridge_max_dist = runif(100,0.0001,0.05) # optimal based on 2 runs: btw .0001 and .002
			
	
		

n_runs = 1000

cal_out_3 = data.frame(
  length_rech=rep(NA,n_runs),
  exp_rech=rep(NA,n_runs), 
  depth_sclr=rep(NA,n_runs),
  ridge_exp=rep(NA,n_runs),
  ridge_max=rep(NA,n_runs),
  nse=rep(NA,n_runs)
)

	# Stop the clock
ptm = proc.time()
	for(yy in 1:n_runs){
		cal_out_3$length_rech[yy] = sample(length_rech_dist,1)
		cal_out_3$exp_rech[yy] = sample(exp_rech_dist,1)
		cal_out_3$depth_sclr[yy] = sample(depth_sclr_dist,1)
		cal_out_3$ridge_exp[yy] = sample(ridge_exp_dist,1)
		cal_out_3$ridge_max[yy] = sample(ridge_max_dist,1)
		
			f = 0.1 # the drainable porosoity we set for recession analysis
			gw_out = NULL
			deep_store = as.numeric(mopey_df_disch$SLZ[1])
			depth = deep_store / f
			deep_rech = as.numeric(c(0,diff(mopey_df_disch$SLZ)) + mopey_df_disch$Q2)
			recent_rech = seq(.01, 1, length=sample(cal_out_3$length_rech[yy],1))^cal_out_3$exp_rech[yy]
			accum_rech = 0
			for(gg in 1:(nrow(mopey_df_disch)))	{
				if(gg > length(recent_rech))	{
					accum_rech = sum(recent_rech * deep_rech[(gg - length(recent_rech) + 1):gg]) / f
				}
				k_arb = as.numeric(k_fof_D(D_thickness = depth))
				gw_ridge = ifelse(depth < D_wet, ((D_wet - depth) / D_wet)^cal_out_3$ridge_exp[yy], cal_out_3$ridge_max[yy])
				gw_out_iter =
				  f * 				# transmissivity = f * 2 * k * dh/dx * h_outlet * Length_of_streams / Area
				  k_arb *			# 2k to account for 2 hillslopes per reach of stream
				  #		(((sin(i_slp) * (B/2)) + 2 * depth ) / sqrt(B/2)) *	# dh / dx, while accounting for slope of imperm layer
				  (((sin_slp * B_half * gw_ridge) + depth ) / (B * gw_ridge)) *	# dh / dx, while accounting for slope of imperm layer
				  (depth / 4 + accum_rech) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
				  strm_crs_sxn					# 
				if(gw_out_iter > deep_store[gg]) {gw_out_iter = deep_store[gg]}
				
				new_deep_store = deep_store[gg] +
						deep_rech[gg] -
						gw_out_iter
				deep_store = c(deep_store, new_deep_store)
				depth = new_deep_store / f

				gw_out = c(gw_out, gw_out_iter)
				
			}
				
			mopey_df_disch$deep_store = deep_store[-length(deep_store)]
			mopey_df_disch$deep_rech = deep_rech
			mopey_df_disch$gw_out = gw_out
	
	these_days = 3300:3965
	old_nse = NSE(mopey_df_disch$Qg[-(1:365)], mopey_df_disch$Q[-(1:365)])
	new_nse = NSE((mopey_df_disch$Q0[-(1:365)] + mopey_df_disch$Q1[-(1:365)] + 	mopey_df_disch$gw_out[-(1:365)]), 
		mopey_df_disch$Q[-(1:365)])
	plot(mopey_df_disch$Q[these_days], lwd=2, type='l',
		ylim=c(0,max(mopey_df_disch$Q[these_days], na.rm=TRUE)),
		#log='y',
		ylab = "Q (mm)", xlab = "Date",
  #!!! TEMPFIX
#		main=paste(mopey_catalogue[this_ws,"Name"],
		main=paste(this_ws,
	#!!! TEMPFIX	           
		": \n old NSE =", round(old_nse,2), "new =", round(new_nse,2)))
	lines(mopey_df_disch$Qg[these_days], lwd=2, col='red', lty=1)
	lines(mopey_df_disch$Q0[these_days] + mopey_df_disch$Q1[these_days] + mopey_df_disch$gw_out[these_days],
		col='blue', lwd=2)

	cal_out_3$nse[yy] = new_nse
	
	}

	cal_sort_3 = cal_out_3[order(cal_out_3$nse),]
	#cal_sort_ny for ny

 plot(cal_out_3$length_rech, cal_out_3$nse)
 plot(cal_out_3$exp_rech, cal_out_3$nse)
 plot(cal_out_3$depth_sclr, cal_out_3$nse)
 plot(cal_out_3$ridge_exp, cal_out_3$nse)
 plot(cal_out_3$ridge_max, cal_out_3$nse)


##############################################################################
	# plotting and calibrating
	
	
	













	









































this_ws_props = which(aq_chars$gage == mopey_catalogue[this_ws, "USGS_ID"])

D_dry = aq_chars$B_depth_dry[this_ws_props]
D_wet = aq_chars$B_depth_wet[this_ws_props]
K_dry = aq_chars$B_K2_dry[this_ws_props]
K_wet = aq_chars$B_K2_wet[this_ws_props]
Wo = aq_chars$TOT_STREAM_LENGTH[this_ws_props]  * 1000	# stream length in meters
i_slp =	atan((1 + aq_chars$TOT_BASIN_SLOPE[this_ws_props]) / 100) 					# this is i in brut's equation, but I'll probably use i for something else later
p =	1																# constant / fitting parameter; set to 1 since that what I did when I originally solved
f = 0.1																# set to 0.1 since that's what I did when originally solved
A = aq_chars$TOT_BASIN_AREA[this_ws_props] * 1000 *1000								# basin area in square meters
B = (A / (2*Wo)) / cos(i_slp)										# aquifer breadth


plot(x=c(D_dry, D_wet), y=c(K_dry, K_wet))
con_mod = lm(c(K_dry, K_wet) ~ 	c(D_dry, D_wet))

	# funscion for variable conductivity
k_fof_D = function(D_thickness) {
	#max(c(con_mod$coef[1] + con_mod$coef[2] * D_thickness, 0.001))	#why did I put this max here??
	con_mod$coef[1] + con_mod$coef[2] * D_thickness
}


mopey_df_disch$deep_rech = c(0,diff(mopey_df_disch$SLZ) + mopey_df_disch$Q2)
f = 0.1 # the drainable porosoity we set for recession analysis
mopey_df_disch$gw_out = NA
mopey_df_disch$deep_store = NA
mopey_df_disch$deep_store[1] =  mopey_df_disch$SLZ[1]
depth = mopey_df_disch$deep_store[1] / f

	













































































################################################################
### model 1: the basic 3 component HBV runoff module
# snow and glacier module
  sfcf = runif(100, .5, 4)#1, 2)	#snowfall correction factor [-]
  tr   = runif(100, -2, 4)#-1, 3)	#solid and liquid precipitation threshold temperature [C]
  tt   = runif(100, -1, 4)#0, 3)	#melt temperature [C]
  fm   = runif(100, .5, 6)#1, 4)	#snowmelt factor [mm/C]
  fi   = runif(100, 2, 10)#4, 8)	#icemelt factor [mm/C]
  fic  = runif(100, 4, 10)#6, 6)	#debris-covered icemelt factor [mm/C]

  # soil module
  fc   = runif(100, 50, 400)#100, 200)
  lp   = runif(100, .2, 1)#0.5, 1)
  beta_soils = runif(100, 1, 3)

  # routing module
  k0   = runif(100, .01, .9)#0.09, 0.1)
  k1   = runif(100, .001, .5)#0.05, 0.07)
  k2   = runif(100, .0001, .1)#0.05)
  uz1  = runif(100, .1, 20)#1, 5)
  perc = runif(100, .0001, 20)#.8, 2)
















mopey_df_ppt = cbind(mopey_df,	
			SnowGlacier_HBV(model = 1,
				inputData = cbind(mopey_df$T_mean, mopey_df$P),
				initCond = c(0,2,0),	#SWE0, surface type (2=soil), surface area of glaciers as ratio [-]
				param = c(	cal_out$sfcf[jj],		#SFCF - snowfall correction factor [-]
							cal_out$tr[jj],		#Tr - solid and liquid precipitation threshold temperature [C]
							cal_out$tt[jj],		#Tt - melt temperature [C]
							cal_out$fm[jj],		#fm - snowmelt factor [mm/C]
							cal_out$fi[jj],		#fi - icemelt factor [mm/C]
							cal_out$fic[jj])	 	#fic - debris-covered icemelt factor [mm/C]
		)	)

			# soils model
		mopey_df_rech = cbind(mopey_df_ppt,
			Soil_HBV(
			   model = 1,
			   inputData = cbind(mopey_df_ppt$Total, mopey_df_ppt$E),
			   initCond = c(50,1),	# initial soil moisture in mm, then relative ratio of soil over teh whole the whole basin
			   param = c(	cal_out$fc[jj],			#FC - soil field capacity [mm]
							cal_out$lp[jj],			#LP - parameter ot get actual ET [-]
							cal_out$beta_soils[jj])	#beta - exponential value for nonlinear relations between soil storage and runoff
		)	)


		mopey_df_disch = cbind(mopey_df_rech,
			Routing_HBV(
				model = 1,	# model=1 gives three stores with routing for each
				lake = FALSE,
				inputData = cbind(mopey_df_rech$Rech),	# recharge time series
				initCond = c(10,10,10),	# initial storage in each reservoir in mm
				param = c(	cal_out$k0[jj],	#KO - top bucket (STZ) storage constant [1/t]
							cal_out$k1[jj],	#K1 - intermediate bucket (SUZ) storage constant [1/t]
							cal_out$k2[jj],	#K2 - lower bucket (SLZ) storage constant [1/t]
							cal_out$uz1[jj],	#UZL - max flux rate between STZ and SUZ [mm/t] 
							cal_out$perc[jj])	#PERC - max flux rate between SUZ and SLZ [mm/t]
		)	)


		the_nse = NSE(mopey_df_disch$Qg, mopey_df_disch$Q)
		plot(mopey_df_disch$Q[635:1000], lwd=2, type='l', ylim=c(0,max(mopey_df_disch$Q[635:1000,],na.rm=TRUE)),
			#log='y',
			ylab = "Q (mm)", xlab = "Date",
			main=paste(mopey_catalogue[this_ws,"Name"], ": NSE =", round(the_nse,2)))
		lines(mopey_df_disch$Qg[635:1000], lwd=2, col='blue', lty=1)

		cal_out$nse[jj] = the_nse
	}






































# snow and glacier module
glacier_range <- rbind(
  sfcf = c(1, 2),	#snowfall correction factor [-]
  tr   = c(-1, 3),	#solid and liquid precipitation threshold temperature [C]
  tt   = c(0, 3),	#melt temperature [C]
  fm   = c(1, 4),	#snowmelt factor [mm/C]
  fi   = c(4, 8),	#icemelt factor [mm/C]
  fic  = c(6, 6)	#debris-covered icemelt factor [mm/C]
)


  # soil module
soil_range <- rbind(
  fc   = c(100, 200),
  lp   = c(0.5, 1),
  beta = c(1, 3)
)

  # routing module (here I will give you the correct values)
routing_range <- rbind(
  k0   = c(0.09, 0.09),
  k1   = c(0.07, 0.07),
  k2   = c(0.05, 0.05),
  uzl  = c(5, 5),
  perc = c(2, 2)
)

  # transfer function module (I will give the correct value)
tf_range <- rbind(
  bmax = c(2.25, 2.25)
)



    param_snow = c(1.1, 0, 0, 2.5)
    param_soil = c(150, 0.90, 1.5)
    param_route = c(0.09, 0.07, 0.05, 5, 2)
    param_tf = c(3.00)





# paramaters I don't need to calibrate
	
	 # precip model
	precip_range <- rbind(
	  p_grad  = c(5, 25)
	)

	 # air temperature model
	tair_range <- rbind(
	  t_grad  = c(-9.8, -2)
	)



## Case example with the first model
inputMatrix <- cbind(
                     runif(n = 200, max = 100, min = 0),
                     runif(n = 200, max = 50, min = 5),
                     runif(n = 100, max = 3, min = 1)
                     )

routeMod1   <- Routing_HBV(model = 1, lake = TRUE, inputData = inputMatrix,
                     initCond = c(10, 15, 20), param = c(0.1, 0.05, 0.001, 1, 0.8))





############################################################
# part 1
# Precip (from gage and possibly include elevation gradient)

	# we don't need to use this unless we want to introduce a topo gradient for ppt
Precip_model(
       model, # =1 for linear ppt gradient; =2 for linear ppt gradient with upper threshold
       inputData, # ppt (mm/t) from gage (vector)
       zmeteo,	# altitude of ppt gage (masl) (numeric)
       ztopo,	# target height (masl) (numberic)
       param # ppt gradient (%/100m); optional second value (for model=2) that sets max threshold elevation above which ppt doesn't increase
)


############################################################
# part 2
# Timeseries of rainfall and snowmelt (from part 1 and also bring in Temp)

SnowGlacier_HBV(
       model, # =1 for temeprature index model; =2 or =3 for variable snow cover and variable glacier cover
       inputData, # two columns, column_1= air temperature time series, column_2= ppt series
       initCond, # 3 values, 1st is SWE initial condition (mm), 2nd is 1/2/3 for clean ice/soil/dirty ice, 3rd set to 0 since we have no glaciers
	   param # snowfall correction factor (SFCF) [-], solid / liquid ppt threshold temperature (Tr), melt temperature (Tt), snowmelt factor (fm) in mm/T>0, icemelt factor (fi), debris-covered icemelt factor
)
	# outputs 1) rainfall; 2) snowfall; 3) SWE; 4) melted snow; rainfall + melted snow

############################################################
# part 3
# Soil routing model

Soil_HBV(
       model,
       inputData,
       initCond,
       param
       )


############################################################
# part 4
# streamflow partitioning model






# we first consider the SnowGlacier module to take into account 
# precipitation partioning and the snow accumulation/melting. 
snow_module <-
  SnowGlacier_HBV(model = 1, 
                  inputData = as.matrix( lumped_hbv[ , c('T(ºC)', 'P(mm/d)')] ),
                  initCond = c(20, 2), 
                  param = c(1.20, 1.00, 0.00, 2.5) )
