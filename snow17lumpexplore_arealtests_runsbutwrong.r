{rm(list = ls()) 
setwd("~/Documents/Rscripts/snow17")

library(tidyverse)
library(lubridate) #doesn't always automatically read-in with tidyverse
library(data.table)

df <- read.csv("stcroix.csv") 
df <- df %>% mutate(date = mdy(date), tavg_c = ((tmax_c + tmin_c)/2), doy = yday(date)) %>% filter(date > "1980-8-16")

prcp <- df$p_mm
tavg <- df$tavg_c
doy  <- df$doy

# input timestep interval in hours, current setup requires each to be the same
ts_t <- 24  # temperature [C]
ts_p <- 24  # precipitation [mm]

snow17 <- function(par, prcp, tavg, elev, doy,
          ini.tstep.state = c(0,           0,              0,    0,                 0,                 0))        {
                            #we_solid[mm]  #we_liquid[mm]  #ati  #heatdeficit[mm]   #snowdepth_tot[mm] #swe
 # set parameters (major or minor as assigned by model creator E. Anderson) 
 elev   <-   5000 #representative mean areal elevation [m]             
 scf    <-   1.1  #(major) correction for snow gage deficiency, eg under reporting snow depth from wind/sublimation [unitless]
 mfmax  <-   1.0  #(major) max melt factor during non-rain periods [mm/C/timestep] (varies with forest type/aspect, prevailing wind, etc)
 mfmin  <-   0.1  #(major) min melt factor during non-rain periods [mm/C/timestep] (varies with forest type/aspect, etc)
 uadj   <-   0.1  #(major) avg wind function during rain on snow events [mm/mb/C]
 si     <-   600  #(major) #threshold above which there's always 100% snow cover [mm]
 pxtemp <-   0.0  #(minor) snow/rain threshold temp [C]
 mbase  <-   0.1  #(minor) base temperature for snowmelt calcs [C]
 tipm   <-   0.1  #(minor) antecedent temperature index [unitless], intended to represent temperature inside the snow cover but near the surface, for a gradient
 plwhc  <-   0.1  #(minor) max amount of liquid water able to be held by snowpack  (percent of liquid water holding capacity) [unitless]
 nmf    <-   0.1  #(minor) maximum negative melt factor [mm/C/timestep]
 daygm  <-   0.1  #(minor) constant melt rate at snow-soil interface [mm/timestep]

par <- c(scf, mfmax, mfmin, uadj, si, pxtemp, mbase, tipm, plwhc, nmf, daygm)
 
 ### basin-specific areal depletion curve
 #meanarealwe_to_ai_x <- c(0, 0.07, .1, .13, .17, .2, .24, .27, .32, .49, 1)
 #percentarealsnow_y <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
 #dt_arealdeplete <- data.table(meanarealwe_to_ai_x, percentarealsnow_y)
 
 meltandrain <- vector(mode = "numeric", length = length(prcp))
 #aesc <- vector(mode = "numeric", length = length(meanarealwe_to_ai_x))
 
 # loop through each timestep 
 for (i in 1:length(prcp)) {
   
   # set initial states (update at loop end)
   we_solid           <- ini.tstep.state[1]
   we_liquid          <- ini.tstep.state[2]
   ati                <- ini.tstep.state[3]
   heatdef            <- ini.tstep.state[4]
   snowdepth_tot      <- ini.tstep.state[5]
   swe                <- ini.tstep.state[6]
   
   # set current ground temperature and precipitation
   grndairtemp <- tavg[i]  # mean air temperature of time step [C]
   precip <- prcp[i]  # total precipitation for time step [mm]
   
   # set binary precipitation form & amounts (original fortran version has other choices)
if (grndairtemp <= pxtemp) {
   # snowing
   we_newsnow <- precip
   rain <- 0
 } else {
   # raining
   we_newsnow <- 0
   rain <- precip
 }
 
 # new snow swe [mm]
 we_newsnow_gadj <- we_newsnow * scf  # (bias correct the gage(s))
 
 # new snow density [g/cm^3]
 if (grndairtemp <= -15) {    #snow density assumed about the same below ~15 C 
   den_newsnow <- 0.05 
 }
 else {
   den_newsnow <- 0.05 + 0.0017  * grndairtemp^1.5  #manual's snow density scaling increase with temperature
 }
 
 # new snow depth [cm]  #centimeters for output convenience
 depth_newsnow  = (0.1 * we_newsnow_gadj)/den_newsnow
 
 # new snow temperature
 if (grndairtemp < 0)       #if air temp is below 0C
   t_newsnow <- grndairtemp #new snow temp = air temp
 else  
   t_newsnow <- 0   #otherwise it's 0C
 
 ## areal extent curve
 # to apply heat exchange only to snow covered area, also implicitly reduces melt factor rate as snow recedes
 #swe              <- 300#we_solid + we_liquid #[mm]
#max_we_ap         <- swe #[mm] #max water equivalent during accumulation period - revisit
#arealindex        <- min(max_we_ap, si) #[mm]
#meanarealwe       <- swe # [mm] #swe before any melt below  ##same as max_we_ap, right?
#meanarealwe_to_ai <- meanarealwe/arealindex
#fracarealcover    <- with(dt_arealdeplete, approx(meanarealwe_to_ai_x, percentarealsnow_y,  xout = meanarealwe_to_ai)) %>% tail(1) %>% unlist() %>% unname()

 
# energy exchange at snow/air surface when no surface melt
#..
# change (increase) in the heat deficit due to new snowfall [mm] (heat amount needed to heat new snow to 0C)
# 80 cal/g: latent heat of fusion
# 0.5 cal/g/C: specific heat of ice

heatdefincreasefromnewsnow <- - (t_newsnow * we_newsnow_gadj)/(80/0.5)

# define/update antecedent temperature index (represents near surface snow temp from past snow & temp history),
# most recent air temps weighed by decreasing amounts
# if 'significant' new snow (>1.5 mm/hr), use new snow's temp, otherwise compute ati from last ati as variable going into representation of shallow snow temp evolution

if (we_newsnow_gadj > 1.5 * ts_p) {
 ati <- t_newsnow  
} else {
 tipm_ts_t <- 1 - ((1 - tipm)^(ts_t/6))
   ati <- ati + tipm_ts_t * (grndairtemp - ati) 
    }

ati <- min(ati, 0) #ati can't exceed 0

# define sinusoidal seasonal variation in the non-rain melt factor (assume 365 day year)
N_Mar21 <- doy[i] - 80  
sp <- (0.5 * sin((N_Mar21 * 2 * pi)/365)) + 0.5  # the seasonal pattern assumed
sp_adj <- 1  # Seasonal variation adjustment, none when latitude is below ~54N
mf <- ts_t/6 * ((sp * sp_adj * (mfmax - mfmin)) + mfmin)  # the seasonally varying non-rain melt factor
t_snowsurf <- min(0, grndairtemp)

# now with melt factor for the day of year and the snow surface temperature, get the temperature gradient (increase or decrease)
## driven change in heat deficit due to the temperature gradient within the upper part of the snowpack [mm], for only the snow covered fraction:
heatdefchangefromprofilegradient <- nmf * ts_p/6 * mf/mfmax * (ati - t_snowsurf) #* fracarealcover
              

# solar & atmospheric snow melt
t_rain <- max(grndairtemp, 0)  # rain temp is air temp as long as it's above 0C
if (rain > 0.25 * ts_p) {  # if 'significant' rate of 0.25 mm/hr
 # rain-on-snow melt   #assumed overcast, high humidity (>90%), emissivity of 1 at cloud elev temp, which is close to ground temp
 stefan_bolt <- 6.12 * (10^(-10))  # Stefan-Boltzman constant [mm/C/hr]
 e_sat <- 2.7489 * (10^8) * exp((-4278.63/(grndairtemp + 242.792)))  # saturated vapor pressure at grndairtemp [mb]
 p_atm <- 33.86 * (29.9 - (0.335 * (elev/100)) + (0.00022 * ((elev/100)^2.4)))  # elevation is in hundreds of meters (incorrect in snow17 manual)
 term1 <- stefan_bolt * ts_p * (((grndairtemp + 273)^4) - (273^4))
 term2 <- 0.0125 * rain * t_rain
 term3 <- 8.5 * uadj * (ts_p/6) * ((0.9 * e_sat - 6.11) + (0.00057 * p_atm * grndairtemp))
 melt_satmos <- term1 + term2 + term3
 melt_satmos <- max(melt_satmos, 0) #melt can't be negative
 melt_satmos <- melt_satmos #* fracarealcover  #only snow covered fraction can melt
} 
 else if ((rain <= 0.25 * ts_p) && (grndairtemp > mbase)) { #if insignificant rain and air temp is above snowmelt threshold (usually 0C)
 # non-rain or very little rain melt - melt factor driven and accomodates heat from small rain amounts
 melt_satmos <- (mf * (grndairtemp - mbase) * (ts_p/ts_t)) + (0.0125 * rain * t_rain)
 melt_satmos <- max(melt_satmos, 0) #melt can't be negative
 melt_satmos <- melt_satmos #* fracarealcover #only snow covered fraction can melt
} else {
 melt_satmos <- 0  #otherwise, no solar/atmospheric melt without significant rain or temps
}

# ripeness of the snow cover #
# we_solid: water equivalent of the ice portion of the snow cover [mm]
# we_liquid: liquid water in the snow [mm]
# packliquidmaxmm: liquid water storage capacity [mm]
# Qw: Amount of available water due to melt and rain [mm]

we_solid <- we_solid + we_newsnow_gadj   # water equivalent of the whole ice portion of the snow cover [mm]


heatdef <- max(heatdef + heatdefincreasefromnewsnow + heatdefchangefromprofilegradient, 0)  # [mm] #heat deficit can't exceed zero (when melting starts)
if (heatdef > (1/3 * we_solid)) {   #limit a deep snowpack's ability to keep its surface from melting
  # set heat deficit limit
  heatdef <- 1/3 * we_solid
}

if (melt_satmos < we_solid) {
  we_solid <- we_solid - melt_satmos
  Qw <- melt_satmos + rain
  packliquidmaxmm <- plwhc * we_solid
  
  if ((Qw + we_liquid) > (heatdef + heatdef * plwhc + packliquidmaxmm)) {
    # snow is ripe!
    
    E <- Qw + we_liquid - packliquidmaxmm - heatdef - (heatdef * plwhc)  # excess liquid water [mm]
    we_solid <- we_solid + heatdef  # we_solid increases because water refreezes as heat deficit is decreased
    we_liquid <- packliquidmaxmm + plwhc * heatdef  # fills liquid water capacity
    heatdef <- 0
    
  } else if ((Qw + we_liquid) >= heatdef) {
    
    E <- 0
    we_solid <- we_solid + heatdef  # we_solid increase because water refreezes as heat heatdef is decreased
    we_liquid <- we_liquid + Qw - heatdef
    heatdef <- 0
    
  } else if ((Qw + we_liquid) < heatdef) {
 
    # snow isn't ripe
    E <- 0
    we_solid <- we_solid + Qw + we_liquid  # we_solid increases because water refreezes as heat deficit is decreased
    heatdef <- heatdef - Qw - we_liquid
  }
  
} else {
  
  melt_satmos <- we_solid + we_liquid  
  we_solid <- 0
  we_liquid <- 0
  Qw <- melt_satmos + rain
  E <- Qw
  
}

if (heatdef == 0) {
  ati = 0
}

# lithospheric snow melt - constant daily amount of melt that takes place at the snow-ground interface

if (we_solid > daygm) { # if more swe than the daily ground melt rate assumed
  
  melt_litho_liqloss    <- (daygm/we_solid) * we_liquid #* fracarealcover  
  melt_litho_solidloss  <- daygm #* fracarealcover
  melt_litho            <- melt_litho_liqloss + melt_litho_solidloss
  we_solid              <- we_solid - melt_litho_solidloss
  we_liquid             <- we_liquid - melt_litho_liqloss
  
  E <- E + melt_litho

  swe <- we_solid + we_liquid
  } 
  else {
      melt_litho <- 0 # if bare or less swe than daily ground melt, no ground melt
      E <- E + melt_litho
      
    }
  
    meltandrain[i] <- E
  #  aesc[i] <- fracarealcover
    ini.tstep.state <- c(we_solid, we_liquid, ati, heatdef, snowdepth_tot, swe)
   # swe <- ini.tstep.state[1] +ini.step.state[2]
  }
  return(meltandrain)
  #return(aesc)
}

snowmelt <- snow17(par, prcp, tavg, elev, doy)
df <- data.frame(df, snowmelt)
ggplot(df, aes(date, snowmelt)) + geom_line()
}




#liqlag <- 5.33 * (1 - exp  ((-0.03*(ts_p/6)* we_solid )/E )   )