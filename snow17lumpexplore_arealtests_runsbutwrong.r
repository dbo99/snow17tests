{rm(list = ls()) 
setwd("~/Documents/Rscripts/snow17")

library(tidyverse)
library(lubridate)

df <- read.csv("stcroix.csv") 
df <- df %>% mutate(date = mdy(date), tavg_c = ((tmax_c + tmin_c)/2), doy = yday(date))

prcp <- df$p_mm
tavg <- df$tavg_c
doy  <- df$doy
elev <- 3000

par <- c(1.1, 0.0, 1.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
snow17 <- function(par, prcp, tavg, elev, doy, ini.states = c(100, 0, 20, 0, 130, 1000)) {
  
 # set parameters (major or minor parameter as assigned by model creator E. Anderson) [units]
                      # guesses from ranges given in papers
 SCF    <-   par[1]   #1.1  (MAjor) correction for snow gauge deficiency, eg under reporting snow depth from wind [unitless]
 PXTEMP <-   par[2]   #0.0  (minor) snow/rain threshold temp [C]
 MFMAX  <-   par[3]   #1.0  (MAjor) max melt factor during non-rain periods [mm/C/timestep] (varies with forest type/aspect, etc)
 MFMIN  <-   par[4]   #0.1  (MAjor) min melt factor during non-rain periods [mm/C/timestep] (varies with forest type/aspect, etc)
 UADJ   <-   par[5]   #0.1  (MAjor) avg wind function during rain on snow events [mm/mb/C]
 MBASE  <-   par[6]   #0.1  (minor) base temperature for snowmelt calcs [C]
 TIPM   <-   par[7]   #0.1  (minor) antecedent temperature index [unitless], intended to represent temperature inside the snow cover but near the surface, for a gradient
 PLWHC  <-   par[8]   #0.1  (minor) max amount of liquid water able to be held by snowpack  (percent of liquid-water capacity) [unitless]
 NMF    <-   par[9]   #0.1  (minor) maximum negative melt factor [mm/C/timestep]
 DAYGM  <-   par[10]  #0.1  (minor) constant melt rate at snow-soil interface [mm/timestep]
 
 # Define constants
 dtt <- 24  # T data timestep hours
 dtp <- 24  # P data timestep hours
 
 meltNrain <- vector(mode = "numeric", length = length(prcp))
 
 # loop through each day (daily timestep for now)
 for (i in 1:length(prcp)) {
   
   # set initial states
   we_solids          <- ini.states[1]
   ATI                <- ini.states[2]
   we_liquid          <- ini.states[3]
   deficit            <- ini.states[4]
   swe                <- ini.states[5]
   snowdepth_tot      <- ini.states[6]
   
   # set current temperature and precipitation
   t_i <- tavg[i]  # mean air temperature of time step [C]
   p_i <- prcp[i]  # precipitation over this time step [mm]
   
   # binary precipitation form
if (t_i <= PXTEMP) {
   # snowing
   swe_newsnow <- p_i
   rain <- 0
 } else {
   # raining
   swe_newsnow <- 0
   rain <- p_i
 }
 
 # new snow swe [mm]
 swe_newsnow_gadj <- swe_newsnow * SCF  # (gage bias correction)
 
 # new snow density [g/cm^3]
 if (t_i <= -15) {
   den_newsnow <- 0.05 
 }
 else {
   den_newsnow <- 0.05 + 0.0017* t_i^1.5
 }
 
 # new snow depth [cm]
 depth_newsnow  = (0.1 * swe_newsnow_gadj)/den_newsnow
 
 # new snow temperature
 if (t_i < 0) 
   t_newsnow <- t_i 
 else 
   t_newsnow <- 0
 
 ## areal extent curve
 # reduces melt factor rate as snow recedes
 
 max_we_ap   <- swe       # [mm] #max water equivalent during accumulation period
 SI          <- 600 # [mm] #threshold above which 100% snow cover always exists, ie >6m, basin's always covered
 arealindex  <- min(max_we_ap, SI)
 meanarealwe <- swe # [mm] #swe before any melt below  ##same as max_we_ap, right?
 meanarealwe_to_ai <- meanarealwe/arealindex
 
 ### areal depletion curve - basin specific!
 meanarealwe_to_ai_x <- c(0, 0.07, .1, .13, .17, .2, .24, .27, .32, .49, 1)
 percentarealsnow_y <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
 df_arealdeplete <- data.frame( meanarealwe_to_ai_x, percentarealsnow_y)
 #prcntarealsnow <- with(df_arealdeplete, approx(meanarealwe_to_ai_x, percentarealsnow_y,  xout = meanarealwe_to_ai)) %>% tail(1) %>% unlist() %>% unname()
 #prcntarealsnow <- paste0(prcntarealsnow) #as.numeric(prcntarealsnow)
 #prcntarealsnow <- as.numeric(prcntarealsnow)
 prcntarealsnow <- 0.7
 
# energy exchange at snow/air surface when no surface melt
#..
# change (increase) in the heat deficit due to new snowfall [mm] (heat amount needed to heat new snow to 0C)
# 80 cal/g: latent heat of fusion
# 0.5 cal/g/C: specific heat of ice

delta_HD_newsnow <- - (t_newsnow * swe_newsnow_gadj)/(80/0.5)

# define/update Antecedent Temperature Index (represents near surface snow temp from past snow & temp history),
# most recent air temps weighed by decreasing amounts
# if new snow w/>6 mm/hr rate, use new snow's temp, otherwise update from last ATI

if (swe_newsnow_gadj > 1.5 * dtp) {
 ATI <- t_newsnow  
} else {
 TIPM_dtt <- 1 - ((1 - TIPM)^(dtt/6))
   ATI <- ATI + TIPM_dtt * (t_i - ATI) 
    }

ATI <- min(ATI, 0) #ATI shouldn't exceed 0

## change (incr. or decr.) in heat deficit due to temperature gradient [mm]

# first define seasonal variation in the non-rain melt factor (assume 365 day year)
N_Mar21 <- doy[i] - 80
Sv <- (0.5 * sin((N_Mar21 * 2 * pi)/365)) + 0.5  # the seasonal pattern assumed
Av <- 1  # Seasonal variation adjustment, none when lat < ~54N
Mf <- dtt/6 * ((Sv * Av * (MFMAX - MFMIN)) + MFMIN)  # the seasonally varying non-rain melt factor
t_snowsurf <- min(0, t_i)
# now with Mf, get the change in heat deficit
delta_HD_T <- NMF * dtp/6 * Mf/MFMAX * (ATI - t_snowsurf) * prcntarealsnow
              

# solar & atmospheric snow melt
t_rain <- max(t_i, 0)  # Temperature of rain (deg C), t_i or 0C, whichever greater
if (rain > 0.25 * dtp) {
 # rain-on-snow melt
 stefan <- 6.12 * (10^(-10))  # Stefan-Boltzman constant (mm/K/hr)
 e_sat <- 2.7489 * (10^8) * exp((-4278.63/(t_i + 242.792)))  # Saturated vapor pressure at t_i (mb)
 P_atm <- 33.86 * (29.9 - (0.335 * (elev/100)) + (0.00022 * ((elev/100)^2.4)))  # elevation is in HUNDREDS of meters (this is incorrectly stated in the manual)
 term1 <- stefan * dtp * (((t_i + 273)^4) - (273^4))
 term2 <- 0.0125 * rain * t_rain
 term3 <- 8.5 * UADJ * (dtp/6) * ((0.9 * e_sat - 6.11) + (0.00057 * P_atm * t_i))
 melt_satmos <- term1 + term2 + term3
 melt_satmos <- max(melt_satmos, 0)
 melt_satmos <- melt_satmos * prcntarealsnow
 
} 

 else if ((rain <= 0.25 * dtp) && (t_i > MBASE)) {
 # non-rain melt
 melt_satmos <- (Mf * (t_i - MBASE) * (dtp/dtt)) + (0.0125 * rain * t_rain)
 melt_satmos <- max(melt_satmos, 0)
 melt_satmos <- melt_satmos * prcntarealsnow
} else {
 melt_satmos <- 0
}

# ripeness of the snow cover
# we_solids: water equivalent of the ice portion of the snow cover
# we_liquid: liquid water in the snow
    # we_liquidmaxlim: liquid water storage capacity
    # Qw: Amount of available water due to melt and rain

we_solids <- we_solids + swe_newsnow_gadj   # water equivalent of the ice portion of the snow cover [mm]


deficit <- max(deficit + delta_HD_newsnow + delta_HD_T, 0)  # deficit <- heat deficit [mm]
if (deficit > (0.33 * we_solids)) {      #check manual for 1/3 clarification
  # limits of heat deficit
  deficit <- 0.33 * we_solids
}

if (melt_satmos < we_solids) {
  we_solids <- we_solids - melt_satmos
  Qw <- melt_satmos + rain
  we_liquidmaxlim <- PLWHC * we_solids
  
  if ((Qw + we_liquid) > (deficit + deficit * PLWHC + we_liquidmaxlim)) {
    # snow is ripe!
    
    E <- Qw + we_liquid - we_liquidmaxlim - deficit - (deficit * PLWHC)  # Excess liquid water [mm]
    we_solids <- we_solids + deficit  # we_solids increases because water refreezes as heat deficit is decreased
    we_liquid <- we_liquidmaxlim + PLWHC * deficit  # fills liquid water capacity
    deficit <- 0
    
  } else if ((Qw + we_liquid) >= deficit) {
    
    
    E <- 0
    we_solids <- we_solids + deficit  # we_solids increases because water refreezes as heat deficit is decreased
    we_liquid <- we_liquid + Qw - deficit
    deficit <- 0
    
  } else if ((Qw + we_liquid) < deficit) {
 
    # THEN the snow is NOT yet ripe
    E <- 0
    we_solids <- we_solids + Qw + we_liquid  # we_solids increases because water refreezes as heat deficit is decreased
    deficit <- deficit - Qw - we_liquid
  }
  
} else {
  
  melt_satmos <- we_solids + we_liquid  
  we_solids <- 0
  we_liquid <- 0
  Qw <- melt_satmos + rain
  E <- Qw
  
  
}

if (deficit == 0) {
  ATI = 0
}

# lithospheric snow melt - constant daily amount of melt that takes place at the snow-terra interface
# (if there's more snow than the daily ground melt rate assumed)
if (we_solids > DAYGM) {
  
  melt_litho_liqloss <- (DAYGM/we_solids) * we_liquid * prcntarealsnow
  melt_litho_solidloss <- DAYGM * prcntarealsnow
  melt_litho <- melt_litho_liqloss + melt_litho_solidloss
  we_solids <- we_solids - melt_litho_solidloss
  we_liquid <- we_liquid - melt_litho_liqloss
  
  E <- E + melt_litho

  swe <- we_solids + we_liquid
  } 
  else {
      melt_litho <- 0
      E <- E + melt_litho
      swe <- 0
      
    }
  
#E <- E * prcntarealsnow + (1 - prcntarealsnow) * p_i
  
    meltNrain[i] <- E
    ini.states <- c(we_solids, ATI, we_liquid, deficit, swe, snowdepth_tot)
    
  }
 

  return(meltNrain)
}

snowmelt <- snow17(par, prcp, tavg, elev, doy)
df <- data.frame(df, snowmelt)
ggplot(df, aes(date, snowmelt)) + geom_line()
}




#liqlag <- 5.33 * (1 - exp  ((-0.03*(dtp/6)* we_solids )/E )   )