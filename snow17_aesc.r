{rm(list = ls()) 
setwd("~/Documents/Rscripts/snow17")

library(tidyverse)
library(lubridate) #doesn't always automatically read-in with tidyverse
library(data.table)
#library(plotly)

df <- read.csv("stcroix.csv") 
df <- df %>% mutate(date = mdy(date), tavg_c = ((tmax_c + tmin_c)/2), doy = yday(date)) %>% 
      filter(date > "1988-8-16", date < "1990-8-16")

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
 scf    <-   0.8  #(major) correction for snow gage deficiency, eg under reporting snow depth from wind/sublimation [unitless]
 mfmax  <-   0.2  #(major) max melt factor during non-rain periods [mm/C/timestep] (varies with forest type/aspect, prevailing wind, etc)
 mfmin  <-   0.02 #(major) min melt factor during non-rain periods [mm/C/timestep] (varies with forest type/aspect, etc)
 uadj   <-   0.2  #(major) avg wind function during rain on snow events [mm/mb/C]
 si     <-   500  #(major) #threshold above which there's always 100% snow cover [mm]
 pxtemp <-   0.0  #(minor) snow/rain threshold temp [C]
 mbase  <-   0.1  #(minor) base temperature for snowmelt calcs [C]
 tipm   <-   0.1  #(minor) antecedent temperature index [unitless], intended to represent temperature inside the snow cover but near the surface, for a gradient
 plwhc  <-   0.05 #(minor) max amount of liquid water able to be held by snowpack  (percent of liquid water holding capacity) [unitless]
 nmf    <-   0.1  #(minor) maximum negative melt factor [mm/C/timestep]
 daygm  <-   0.02 #(minor) constant melt rate at snow-soil interface [mm/timestep]

 par <- c(scf, mfmax, mfmin, uadj, si, pxtemp, mbase, tipm, plwhc, nmf, daygm)
 
 ## basin-specific areal depletion curve
 meanarealwe_to_ai_x <- c(0, 0.07, .1, .13, .17, .2, .24, .27, .32, .49, 1)
 percentarealsnow_y <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
 dt_arealdeplete <- data.table(meanarealwe_to_ai_x, percentarealsnow_y)
 
 meltandrain <- vector(mode = "numeric", length = length(prcp))
 aesc <- vector(mode = "numeric", length = length(prcp))
 
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
#si <- 600
#swe              <- 300#we_solid + we_liquid #[mm]
max_we_ap         <- swe #[mm] #max water equivalent during accumulation period - revisit
arealindex        <- max(1.0e-100, min(max_we_ap, si)) #[mm]
meanarealwe       <- swe # [mm] #swe before any melt below  ##same as max_we_ap, right?
meanarealwe_to_ai <- min(1, meanarealwe/si) #max(0, min(1, meanarealwe/arealindex))
fracarealcover    <- with(dt_arealdeplete, approx(meanarealwe_to_ai_x, percentarealsnow_y,  xout = meanarealwe_to_ai)) %>% tail(1) %>% unlist() %>% unname()
 aesc[i] <- fracarealcover
 
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
heatdefchangefromprofilegradient <- nmf * ts_p/6 * mf/mfmax * (ati - t_snowsurf) * fracarealcover
              

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
 melt_satmos <- max(melt_satmos, 0) # melt can't be negative
 melt_satmos <- melt_satmos * fracarealcover  # only snow covered fraction can melt
} 
 else if ((rain <= 0.25 * ts_p) && (grndairtemp > mbase)) { # if insignificant rain and air temp is above snowmelt threshold (usually 0C)
 # non-rain or very little rain melt - melt factor driven and accomodates heat from small rain amounts
 melt_satmos <- (mf * (grndairtemp - mbase) * (ts_p/ts_t)) + (0.0125 * rain * t_rain)
 melt_satmos <- max(melt_satmos, 0) # melt can't be negative
 melt_satmos <- melt_satmos * fracarealcover # only snow covered area can melt
} else {
 melt_satmos <- 0  #otherwise, no solar/atmospheric melt without significant rain or temps (listed above)
}

# update ice water equivalent with bias-corrected new snow amount
we_solid <- we_solid + we_newsnow_gadj   # water equivalent of total ice portion of the snow cover [mm]


# adjust heat deficit from any new snow and the evolving profile gradient from ground air temp & last new snow's temp history 
heatdef <- max(heatdef + heatdefincreasefromnewsnow + heatdefchangefromprofilegradient, 0)  # [mm] 

# but limit a deep snowpack's ability to keep its surface from melting
if (heatdef >= (1/3 * we_solid)) {   
  # set heat deficit limit
  heatdef <- 1/3 * we_solid  #not in 2006 documentation.., check whether in more recent & whether in ops version 
}

if (melt_satmos < we_solid) { # if solar+atmos  melt is less than the ice's water equivalent
  we_solid <- we_solid - melt_satmos # reduce ice water equivalent by the solar+atmos melt amount
  snowsurfavailliq <- melt_satmos + rain # surface liquid content is sum of solar/atmospheric melt and rain
  liquidstorcap <- plwhc * we_solid # but the pack can retain up to this much, the plwhc % of the solid water equivalent
  
  if ((snowsurfavailliq + we_liquid) > (heatdef + (heatdef * plwhc) + liquidstorcap)) {
    # if the solar+atmos melt + rain + existing liquid > heat deficit + deficit's liq stor cap + solid's liq stor cap
    # ie if there's sufficient available liquid water to overcome the total deficit & liquid storage capacity,
    # the snow is ripe:
    
    # excess liquid water is solar+atmos melt + rain + existing liquid - total deficit - liquid held by the pack
    excessliquid <- snowsurfavailliq + we_liquid - heatdef - (heatdef * plwhc) - liquidstorcap 
    
    # increase the just-reduced we_solid, as water 'refreezes' returning pack up to 0C 
    we_solid <- we_solid + heatdef  #written in manual but seems wrong
    # liquid water in the pack is equal to the % maximum
    we_liquid <- liquidstorcap 
    heatdef <- 0
    }
    else if ((snowsurfavailliq + we_liquid) >= heatdef) {
    # if the solar+atmos melt + rain + existing liquid [mm] > heat deficit, but not its and liquidstorcap's extra capacity
    # the snow's not ripe 
    excessliquid <- 0 # there's no excess liquid (aka no melt and rain)
    # still increase the just-reduced we_solid, as water 'refreezes' returning pack up to 0C 
    we_solid <- we_solid + heatdef  
    # the new amount of liquid water is adjusted by the difference of solar+atmos melt + any rain and the heat deficit
    we_liquid <- we_liquid + snowsurfavailliq - heatdef
    # and the pack is at equilibrium - not cold enough for a deficit or warm enough to melt and run off
    heatdef <- 0
    
  } else if ((snowsurfavailliq + we_liquid) < heatdef) {
    # if solar+atmos melt + rain is less than the heat deficit
    # the snow's again not ripe
    excessliquid <- 0 # there's no excess liquid (aka no melt and rain)
    # the solid equivalent is adjusted by and increase of the new solar+atmos melt and rain input and existing liquid (cold pack freezes everything)
    we_solid <- we_solid + snowsurfavailliq + we_liquid  # we_solid increases because water refreezes as heat deficit is decreased
    # the heat deficit
    heatdef <- heatdef - snowsurfavailliq - we_liquid
  }
  }

  else { #if solar+atmos melt > we_solid
  
  melt_satmos <- we_solid + we_liquid   #solar+atmos melt is simply the swe 
  we_solid <- 0  # any ice is now gone, zero
  we_liquid <- 0 # and interstitial liquid, too
  excessliquid <- melt_satmos + rain # add rain for total excessliquid (aka meltandrain)
  
}
# update antecedent temp index - if there's no deficit there can't be a profile gradient
if (heatdef == 0) {
  ati = 0
}

# for existing snow, reduce swe further slightly, slightly increase excess liquid
# lithospheric snow melt - constant daily amount of melt that takes place at the snow-ground interface

if (we_solid > daygm) { # if more swe than the daily ground melt rate assumed
  
  melt_litho_liqloss    <- (daygm/we_solid) * we_liquid * fracarealcover  
  melt_litho_solidloss  <- daygm * fracarealcover
  melt_litho            <- melt_litho_liqloss + melt_litho_solidloss
  we_solid              <- we_solid - melt_litho_solidloss
  we_liquid             <- we_liquid - melt_litho_liqloss
  
  excessliquid <- excessliquid + melt_litho

  swe <- we_solid + we_liquid
  } 
  else {
      melt_litho <- 0 # if bare or less swe than daily ground melt, no ground melt
      excessliquid <- excessliquid + melt_litho

    }
    # save excess liquid as "meltandrain" output - input to rainfall/runoff model (eg sac-sma)
    meltandrain[i] <- excessliquid
    
    # update states for next time step
    ini.tstep.state <- c(we_solid, we_liquid, ati, heatdef, snowdepth_tot, swe = we_solid + we_liquid)
  
  }
  #return(meltandrain)
  return(aesc)
}

arealextent <- snow17(par, prcp, tavg, elev, doy)
df <- data.frame(df, arealextent)
ggplot(df, aes(date, arealextent)) + geom_line()
#ggplotly(p)


}


#df2 <- data.frame(meanarealwe, arealindex, meanarealwe_to_ai)
#liqlag <- 5.33 * (1 - exp  ((-0.03*(ts_p/6)* we_solid )/excessliquid )   )