library(HBV.IANIGLA) # for running HBV model
#library(prism) # for downloading PRISM climate data (to 1981); might replace by just usign mopex data
library(hddtools)		# for mopex data
library(hydroGOF)		# for nse calculations
library(dataRetrieval)	# for streamflow data (I think)


mopey_catalogue = catalogueMOPEX()							# reading in the mopex catalogue












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


this_ws = 7
mopey_df = tsMOPEX(mopey_catalogue$USGS_ID[this_ws])	# identifying the first ts data from the first idnetified mopex gage

#mopey_df_nonas = subset(mopey_df, !is.na(Q) & !is.na(P))
#mopey_df_nonas$T_mean = apply(mopey_df_nonas[,c("T_max","T_min")],1,mean)
mopey_df$T_mean = apply(mopey_df[,c("T_max","T_min")],1,mean)

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


n_runs = 1000
n_loops = 10

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
		plot(mopey_df_disch$Q[1000:1365], lwd=2, type='l', ylim=c(0,max(mopey_df_disch$Q[635:1000,],na.rm=TRUE)),
			#log='y',
			ylab = "Q (mm)", xlab = "Date",
			main=paste(mopey_catalogue[this_ws,"Name"], ": NSE =", round(the_nse,2)))
		lines(mopey_df_disch$Qg[1000:1365], lwd=2, col='blue', lty=1)

		# calibrating by NSE
		cal_out$nse[jj][-(1:365] = the_nse
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
for(gg in 1:nrow(mopey_catalogue)){
this_ws = gg
if(mopey_catalogue[this_ws, "USGS_ID"] %in% aq_chars$gage) 	print(this_ws)
}
head(mopey_catalogue[7, ])

this_ws = 7
mopey_df = tsMOPEX(mopey_catalogue[this_ws, "USGS_ID"])	# identifying the first ts data from the first idnetified mopex gage


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
deep_rech = as.numeric(c(0,diff(mopey_df_disch$SLZ) + mopey_df_disch$Q2))
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
		plot(mopey_df_disch$Q[these_days], lwd=2, type='l', ylim=c(0,max(mopey_df_disch$Q[2600:2965,],na.rm=TRUE)),
			#log='y',
			ylab = "Q (mm)", xlab = "Date",
			main=paste(mopey_catalogue[this_ws,"Name"],
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
for(gg in 1:nrow(mopey_catalogue)){
	if(mopey_catalogue[gg, "USGS_ID"] %in% aq_chars$gage) 	ws_ihave = c(ws_ihave, gg)
}
head(mopey_catalogue[ws_ihave, ])

#################################################################################
	# initializing data for our gw module
this_ws = ws_ihave[76]				# identifying the ws in mopex
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



#########################################################################################
	# initializing data for HBV
mopey_df = tsMOPEX(mopey_catalogue[this_ws, "USGS_ID"])	# identifying the first ts data from the first idnetified mopex gage



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


n_runs = 1000
n_loops = 10

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
			deep_rech = as.numeric(c(0,diff(mopey_df_disch$SLZ) + mopey_df_disch$Q2))
			recent_rech = seq(.01,1,length=180)^2
			accum_rech = 0
			for(gg in 1:(nrow(mopey_df_disch)))	{
				if(gg > length(recent_rech))	{
					accum_rech = sum(recent_rech * deep_rech[(gg - length(recent_rech) + 1):gg]) / f
				}
				k_arb = as.numeric(k_fof_D(D_thickness = depth))
				gw_ridge = ifelse(depth < D_wet, ((D_wet - depth) / D_wet)^2, .001)
				gw_out_iter = f * 				# transmissivity = f * 2 * k * dh/dx * h_outlet * Length_of_streams / Area
					k_arb *			# 2k to account for 2 hillslopes per reach of stream
			#		(((sin(i_slp) * (B/2)) + 2 * depth ) / sqrt(B/2)) *	# dh / dx, while accounting for slope of imperm layer
					(((sin_slp * B_half * gw_ridge) + depth ) / (B * gw_ridge)) *	# dh / dx, while accounting for slope of imperm layer
					(depth / 8 + accum_rech) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
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

##############################################################################
	# plotting and calibrating
	
	
	these_days = 2300:2665
	old_nse = NSE(mopey_df_disch$Qg[-(1:365)], mopey_df_disch$Q[-(1:365)])
	new_nse = NSE((mopey_df_disch$Q0[-(1:365)] + mopey_df_disch$Q1[-(1:365)] + 	mopey_df_disch$gw_out[-(1:365)]), 
		mopey_df_disch$Q[-(1:365)])
	plot(mopey_df_disch$Q[these_days], lwd=2, type='l',
		ylim=c(0,max(mopey_df_disch$Q[these_days,], na.rm=TRUE)),
		#log='y',
		ylab = "Q (mm)", xlab = "Date",
		main=paste(mopey_catalogue[this_ws,"Name"],
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
	cal_sort = head(cal_out[rev(order(cal_out$nse)),], 100)
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
			length_rech_dist = 1:190			# opt based on 2 runs
			exp_rech_dist = runif(10000,0,4) # opt based on 2 runs: smaller is better
			depth_sclr_dist = 1#runif(10000,1,9) # opt based on 2 runs: smaller is better (1)
			ridge_exp_dist = runif(10000,1,4)	# opt based on 2 runs: doesn't seem to matter
			ridge_max_dist = runif(10000,0.0001,0.02) # optimal based on 2 runs: btw .0001 and .002
			
	
			cal_out_3 = data.frame(
				length_rech=rep(NA,n_runs),
				exp_rech=rep(NA,n_runs), 
				depth_sclr=rep(NA,n_runs),
				ridge_exp=rep(NA,n_runs),
				ridge_max=rep(NA,n_runs),
				nse=rep(NA,n_runs)
				)


n_runs = 10000

	# Stop the clock
ptm = proc.time()
	for(yy in 337:n_runs){
		cal_out_3$length_rech[yy] = sample(length_rech_dist,1)
		cal_out_3$exp_rech[yy] = sample(exp_rech_dist,1)
		cal_out_3$depth_sclr[yy] = sample(depth_sclr_dist,1)
		cal_out_3$ridge_exp[yy] = sample(ridge_exp_dist,1)
		cal_out_3$ridge_max[yy] = sample(ridge_max_dist,1)
		
			f = 0.1 # the drainable porosoity we set for recession analysis
			gw_out = NULL
			deep_store = as.numeric(mopey_df_disch$SLZ[1])
			depth = deep_store / f
			deep_rech = as.numeric(c(0,diff(mopey_df_disch$SLZ) + mopey_df_disch$Q2))
			recent_rech = seq(.01, 1, length=sample(cal_out_3$length_rech[yy],1))^cal_out_3$exp_rech[yy]
			accum_rech = 0
			for(gg in 1:(nrow(mopey_df_disch)))	{
				if(gg > length(recent_rech))	{
					accum_rech = sum(recent_rech * deep_rech[(gg - length(recent_rech) + 1):gg]) / f
				}
				k_arb = as.numeric(k_fof_D(D_thickness = depth))
				gw_ridge = ifelse(depth < D_wet, ((D_wet - depth) / D_wet)^cal_out_3$ridge_exp[yy], cal_out_3$ridge_max[yy])
				gw_out_iter = f * 				# transmissivity = f * 2 * k * dh/dx * h_outlet * Length_of_streams / Area
					k_arb *			# 2k to account for 2 hillslopes per reach of stream
			#		(((sin(i_slp) * (B/2)) + 2 * depth ) / sqrt(B/2)) *	# dh / dx, while accounting for slope of imperm layer
					(((sin_slp * B_half * gw_ridge) + depth ) / (B * gw_ridge)) *	# dh / dx, while accounting for slope of imperm layer
					(depth / cal_out_3$depth_sclr[yy] + accum_rech) * 					# elevation of wt at outlet; for converting to transmissivity; can change later
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
	
	these_days = 3300:3965
	old_nse = NSE(mopey_df_disch$Qg[-(1:365)], mopey_df_disch$Q[-(1:365)])
	new_nse = NSE((mopey_df_disch$Q0[-(1:365)] + mopey_df_disch$Q1[-(1:365)] + 	mopey_df_disch$gw_out[-(1:365)]), 
		mopey_df_disch$Q[-(1:365)])
	plot(mopey_df_disch$Q[these_days], lwd=2, type='l',
		ylim=c(0,max(mopey_df_disch$Q[these_days,], na.rm=TRUE)),
		#log='y',
		ylab = "Q (mm)", xlab = "Date",
		main=paste(mopey_catalogue[this_ws,"Name"],
		": \n old NSE =", round(old_nse,2), "new =", round(new_nse,2)))
	lines(mopey_df_disch$Qg[these_days], lwd=2, col='red', lty=1)
	lines(mopey_df_disch$Q0[these_days] + mopey_df_disch$Q1[these_days] + mopey_df_disch$gw_out[these_days],
		col='blue', lwd=2)

	cal_out_3$nse[yy] = new_nse
	
	}

	cal_sort_3 = cal_out_3[order(cal_out_3$nse),]
	#cal_sort_ny for ny

 plot(cal_sort_3$length_rech, cal_sort_3$nse)
 plot(cal_sort_3$exp_rech, cal_sort_3$nse)
 plot(cal_sort_3$depth_sclr, cal_sort_3$nse)
 plot(cal_sort_3$ridge_exp, cal_sort_3$nse)
 plot(cal_sort_3$ridge_max, cal_sort_3$nse)


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
