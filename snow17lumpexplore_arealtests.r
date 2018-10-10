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
snow17 <- function(par, prcp, tavg, elev, doy, ini.states = c(150, 0, 20, 0, 0.5)) {
  
 # set parameters (major or minor parameter as assigned by model creator E. Anderson) [units]
                      # guesses from ranges given in papers
 SCF    <-   par[1]   #1.1  (major) correction for snow gauge deficiency, eg under reporting snow depth from wind [unitless]
 PXTEMP <-   par[2]   #0.0  (minor) snow/rain threshold temp [C]
 MFMAX  <-   par[3]   #1.0  (major) max melt factor during non-rain periods [mm/C/timestep]
 MFMIN  <-   par[4]   #0.1  (major) min melt factor during non-rain periods [mm/C/timestep]
 UADJ   <-   par[5]   #0.1  (major) avg wind function during rain on snow events [mm/mb/C]
 MBASE  <-   par[6]   #0.1  (minor) base temperature for snowmelt calcs [C]
 TIPM   <-   par[7]   #0.1  (minor) antecedent temperature index [unitless], intended to represent the temperature inside the snow cover but near the surface
 PLWHC  <-   par[8]   #0.1  (minor) max amount of liquid water able to be held by snowpack  (percent of liquid-water capacity) [unitless]
 NMF    <-   par[9]   #0.1  (minor) maximum negative melt factor [mm/C/timestep]
 DAYGM  <-   par[10]  #0.1  (minor) constant melt rate at snow-soil interface [mm/timestep]
 
 # Define constants
 dtt <- 24  # T data timestep hours
 dtp <- 24  # P data timestep hours
 
 meltNrain <- vector(mode = "numeric", length = length(prcp))
 
 # Loop through each day (daily timestep for now)
 for (i in 1:length(prcp)) {
   
   # Set initial states
   swe_solids <- ini.states[1]
   ATI        <- ini.states[2]
   swe_liquid <- ini.states[3]
   Deficit    <- ini.states[4]
   percentarealsnow <- ini.states[5]
   
   # Set current temperature and precipitation
   t_i <- tavg[i]  # mean air temperature of time step [C]
   p_i <- prcp[i]  # precipitation over this time step [mm]
   
   # FORM OF PRECIPITATION
if (t_i <= PXTEMP) {
   # it snows
   swe_newsnow <- p_i
   RAIN <- 0
 } else {
   # it rains
   swe_newsnow <- 0
   RAIN <- p_i
 }
 
 # NEW SNOW SWE [mm]
 swe_newsnow_gadj <- swe_newsnow * SCF  # (gauge bias correction)
 
 # NEW SNOW DENSITY [g/cm^3]
 if (t_i <= -15) {
   den_newsnow <- 0.05 
 }
 else {
   den_newsnow <- 0.05 + 0.0017* t_i^1.5
 }
 
 # NEW SNOW DEPTH [cm]
 depth_newsnow  = (0.1 * swe_newsnow_gadj)/den_newsnow
 
 # NEW SNOW TEMPERATURE
 if (t_i < 0) 
   t_newsnow <- t_i 
 else 
   t_newsnow <- 0
 
# ENERGY EXCHANGE AT SNOW/AIR SURFACE WHEN NO SURFACE MELT

 # change (increase) in the heat deficit due to new snowfall [mm] (heat amount needed to heat new snow to 0C)
# 80 cal/g: latent heat of fusion
# 0.5 cal/g/C: specific heat of ice

delta_HD_newsnow <- -(t_newsnow * swe_newsnow_gadj)/(80/0.5)

# define Antecedent Temperature Index (represents near surface snow temp from past snow & temp history),
# most recent air temps weighed by decreasing amounts
# if new snow w/>6 mm/hr rate, use new snow's temp, otherwise update from last ATI

if (swe_newsnow_gadj > 1.5 * dtp) {
 ATI <- t_newsnow  
} else {
 TIPM_dtt <- 1 - ((1 - TIPM)^(dtt/6))
   ATI <- ATI + TIPM_dtt * (t_i - ATI) 
    }

ATI <- min(ATI, 0) #can't exceed 0

## change (incr. or decr.) in heat deficit due to temperature gradient [mm]

# first define seasonal variation in the non-rain melt factor (assume 365 day year)
N_Mar21 <- doy[i] - 80
Sv <- (0.5 * sin((N_Mar21 * 2 * pi)/365)) + 0.5  # seasonal pattern assumed
Av <- 1  # Seasonal variation adjustment, none when lat < ~54N
Mf <- dtt/6 * ((Sv * Av * (MFMAX - MFMIN)) + MFMIN)  # seasonally varying non-rain melt factor

t_snowsurf <- min(0, t_i)
delta_HD_T <- NMF * dtp/6 * Mf/MFMAX * (ATI - t_snowsurf)

# Snow Melt from the atmosphere
t_rain <- max(t_i, 0)  # Temperature of rain (deg C), t_i or 0C, whichever greater
if (RAIN > 0.25 * dtp) {
 # Rain-on-Snow Melt
 stefan <- 6.12 * (10^(-10))  # Stefan-Boltzman constant (mm/K/hr)
 e_sat <- 2.7489 * (10^8) * exp((-4278.63/(t_i + 242.792)))  # Saturated vapor pressure at t_i (mb)
 P_atm <- 33.86 * (29.9 - (0.335 * (elev/100)) + (0.00022 * ((elev/100)^2.4)))  # Atmospheric pressure (mb) where elevation is in HUNDREDS of meters (this is incorrectly stated in the manual)
 term1 <- stefan * dtp * (((t_i + 273)^4) - (273^4))
 term2 <- 0.0125 * RAIN * t_rain
 term3 <- 8.5 * UADJ * (dtp/6) * ((0.9 * e_sat - 6.11) + (0.00057 * P_atm * t_i))
 melt_atmos <- term1 + term2 + term3
 melt_atmos <- max(melt_atmos, 0)
 
} else if ((RAIN <= 0.25 * dtp) && (t_i > MBASE)) {
 # Non-Rain Melt
 melt_atmos <- (Mf * (t_i - MBASE) * (dtp/dtt)) + (0.0125 * RAIN * t_rain)
 melt_atmos <- max(melt_atmos, 0)
 
} else {
 melt_atmos <- 0
}

# Ripeness of the snow cover
# swe_solids: water equivalent of the ice portion of the snow cover
# swe_liquid: liquid water in the snow
    # swe_liquidmaxlim: liquid water storage capacity
    # Qw: Amount of available water due to melt and rain

swe_solids <- swe_solids + swe_newsnow_gadj   # water equivalent of the ice portion of the snow cover [mm]
#E <- 0  # excess liquid water in the snow cover   

Deficit <- max(Deficit + delta_HD_newsnow + delta_HD_T, 0)  # Deficit <- heat deficit [mm]
if (Deficit > (0.33 * swe_solids)) {      #check manual for 1/3 clarification
  # limits of heat deficit
  Deficit <- 0.33 * swe_solids
}

if (melt_atmos < swe_solids) {
  swe_solids <- swe_solids - melt_atmos
  Qw <- melt_atmos + RAIN
  swe_liquidmaxlim <- PLWHC * swe_solids
  
  if ((Qw + swe_liquid) > (Deficit + Deficit * PLWHC + swe_liquidmaxlim)) {
    # THEN the snow is RIPE
    
    E <- Qw + swe_liquid - swe_liquidmaxlim - Deficit - (Deficit * PLWHC)  # Excess liquid water [mm]
    swe_solids <- swe_solids + Deficit  # swe_solids increases because water refreezes as heat deficit is decreased
    swe_liquid <- swe_liquidmaxlim + PLWHC * Deficit  # fills liquid water capacity
    Deficit <- 0
    
  } else if ((Qw + swe_liquid) >= Deficit) {
    
    # & [[Qw + swe_liquid] <= [[Deficit*[1+PLWHC]] + swe_liquidmaxlim]] THEN the snow is NOT yet ripe, but ice is being
    # melted
    
    E <- 0
    swe_solids <- swe_solids + Deficit  # swe_solids increases because water refreezes as heat deficit is decreased
    swe_liquid <- swe_liquid + Qw - Deficit
    Deficit <- 0
    
  } else if ((Qw + swe_liquid) < Deficit) {
 
    # THEN the snow is NOT yet ripe
    E <- 0
    swe_solids <- swe_solids + Qw + swe_liquid  # swe_solids increases because water refreezes as heat deficit is decreased
    Deficit <- Deficit - Qw - swe_liquid
  }
  
} else {
  
  melt_atmos <- swe_solids + swe_liquid  # melt_atmos >= swe_solids
  swe_solids <- 0
  swe_liquid <- 0
  Qw <- melt_atmos + RAIN
  E <- Qw
  
  
}

if (Deficit == 0) {
  ATI = 0
}

# Snow melt from the ground - constant daily amount of melt that takes place at the snow-soil interface
# (if there's more snow than the daily ground melt assumed)
if (swe_solids > DAYGM) {
  
  melt_terra_liqloss <- (DAYGM/swe_solids) * swe_liquid
  melt_terra_solidloss <- DAYGM
  melt_terra <- melt_terra_liqloss + melt_terra_solidloss
  swe_solids <- swe_solids - melt_terra_solidloss
  swe_liquid <- swe_liquid - melt_terra_liqloss
  
  E <- E + melt_terra
  SWE <- swe_solids + swe_liquid
  
} else {

      melt_terra <- 0
      E <- E + melt_terra
      SWE <- 0
      
    }
    
    meltNrain[i] <- E
    ini.states <- c(swe_solids, ATI, swe_liquid, Deficit)
    
  }
 

  
  return(meltNrain)
}

snowmelt <- snow17(par, prcp, tavg, elev, doy)
df <- data.frame(df, snowmelt)
ggplot(df, aes(date, snowmelt)) + geom_line()
}




#liqlag <- 5.33 * (1 - exp  ((-0.03*(dtp/6)* swe_solids )/E )   )