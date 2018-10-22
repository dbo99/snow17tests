{rm(list = ls()) 
setwd("~/Documents/Rscripts/snow17")

library(tidyverse)
library(lubridate) #should, but doesn't always load with tidyverse
library(data.table)
#library(plotly)

tstep         <- seq(from=as.POSIXct("1984-10-1 0:00", tz="UTC"), to=as.POSIXct("2010-9-30 18:00", tz="UTC"), by="6 hours")  
csv_mat       <- read_csv("sbrc1up_mat.csv") #read-in csv temp timeseries
mat6hr_degf   <- c(na.omit(c(t(csv_mat))))[1:length(tstep)] #convert multi-column data card format to row format
mat6hr_degc   <- (mat6hr_degf - 32) * (5/9)  #convert fahrenheit to celsius
csv_map       <- read_csv("sbrc1up_map.csv") #read-in csv precip timeseries
map6hr_in     <- c(na.omit(c(t(csv_map))))[1:length(tstep)] #convert multi-column data card format to row format
map6hr_mm     <- map6hr_in * 25.4 # convert inches to millimeters

df <- data.frame(tstep, mat6hr_degc, map6hr_mm) %>% mutate(doy = yday(tstep)) #combine in data.frame, add day of year attribute

# time step hourly frequency - current setup requires t's & p's to match
tshrs_t <- 6  # 6-hour input temp data   [C]
tshrs_p <- 6  # 6-hour input precip data [mm]

map <- df$map6hr_mm
mat <- df$mat6hr_degc
doy  <- df$doy

snow17 <- function(par, map, mat, elev, doy,
          ini.tstep.state =  # initialize. assumes timeseries starts late summer
          c(0,   # (1)   we_solid [mm]      
            0,   # (2)   we_liquid [mm]     
            0,   # (3)   atip [-]           
            0,   # (4)   heatdeficit [mm]   
            0,   # (5)   swe [mm]    
            0,   # (6)   swe_b4suffsnewsnow [mm]
            0,   # (7)   swe_aescfdb100_pb [mm] 
            0,   # (8)   aesc [%]  
            0,   # (9)   aesc_b4suffnewsnow [%]
            0,   # (10)  meanarealwe_to_ai [%]                             
            0,   # (11)  integer switch indicating whether last time step had qualifying new snow to leave depletion curve  
            0 )) # (12)  accmax [mm]
          { 
  # set parameters (major or minor as assigned by model creator E. Anderson) 
elev   <-   1768 #representative mean areal elevation [m]     
scf    <-   0.97 #(major) correction for snow gage deficiency, eg under reporting snow depth from wind/sublimation [unitless]
mfmax  <-   0.68 #(major) max melt factor during non-rain periods [mm/C/timestep] (varies with forest type/aspect, prevailing wind, etc)
mfmin  <-   0.15 #(major) min melt factor during non-rain periods [mm/C/timestep] (varies with forest type/aspect, prevailing wind, etc)
uadj   <-   0.09 #(major) avg wind function during rain on snow events [mm/mb/C]
si     <-   750  #(major) threshold above which there's always 100% snow cover [mm]
pxtemp <-   1.0  #(minor) snow/rain threshold temp [C]
mbase  <-   0.5  #(minor) base temperature assumed for snowmelt calcs [C]
tipm   <-   0.1  #(minor) for antecedent temperature index [unitless], intended to represent temperature inside the snow cover but near the surface, for a gradient
plwhc  <-   0.05 #(minor) max amount of liquid water able to be held by snowpack  (percent of liquid water holding capacity) [unitless]
nmf    <-   0.3  #(minor) maximum negative melt factor [mm/C/timestep]
daygm  <-   0.3  #(minor) constant melt rate at snow-soil interface [mm/timestep]
hsnof  <-   0.2  #(unassigned) minimum qualifying hourly snow fall rate for leaving depletion curve [mm/hr]
hsnof2 <-   1.5  #(unassigned) minimum qualifying hourly snow fall rate to set ATI equal to temperature of new snow [mm/hr]
par <- c(elev, scf, mfmax, mfmin, uadj, si, pxtemp, mbase, tipm, plwhc, nmf, daygm)

## set basin/sub-basin specific areal extent depletion curve 
# to apply heat exchange only to snow covered area, also implicitly reduces melt factor rate as snow recedes
meanarealwe_to_ai_x <- c(0.0, 0.05, .09, .15, .23, .37, .56, .72, .86, .93, 1) #remember set/fix lowest to 0.05
percentarealsnow_y  <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1) #base value recommended 0.05 per manual, helps ensure complete seasonal melt
dt_arealdeplete     <- data.table(meanarealwe_to_ai_x, percentarealsnow_y)     #look up table for linear interpolation

meltandrain         <- vector(mode = "numeric", length = length(map))
aesc                <- vector(mode = "numeric", length = length(map))

# loop through each timestep 
for (i in 1:length(map)) {
   
# set initial states (and update at loop end)
we_solid                      <- ini.tstep.state[1]  # water equivalent of ice [mm]  (w_i in manual)
we_liquid                     <- ini.tstep.state[2]  # liquid water [mm] (w_q in manual)
ati                           <- ini.tstep.state[3]  # antecedent temperature index [-]
heatdef                       <- ini.tstep.state[4]  # heat deficit [mm]
swe                           <- ini.tstep.state[5]  # snow water equivalent [mm] (sum we_solid & we_liquid) [mm]
swe_b4suffsnewsnow            <- ini.tstep.state[6]  # snow water equivalent before accumualtion period [mm] (sb in doc)
swe_aescfdb100_pb             <- ini.tstep.state[6]  # swe at time when fresh areal snow extent first drops below 100% for partially bare conditions [mm] (limited by hsnof) (w100 in doc.)
aesc                          <- ini.tstep.state[7]  # areal extent of snow cover [%]
aesc_b4suffnewsnow            <- ini.tstep.state[8]  # areal extent at departure time from depletion curve [%] (sbaesc in doc)
meanarealwe_to_ai             <- ini.tstep.state[9]  # normalized swe index (swe/(min(accmax, si))) [-] (ainorm in manual)
meanarealwe_to_ai_b4newsnow   <- ini.tstep.state[10] # normalized swe index (swe/(min(accmax, si))) [-], at timestep before accumulation period
suffsnowprevtsteptolvdc_int   <- ini.tstep.state[11] # sufficient snow last time step to leave depletion curve? [integer switch] 1 = yes, 0 = no
accmax                        <- ini.tstep.state[12] # maximum swe during 'accumulation period' (defined as swe > 3*swe_b4suffsnewsnow)



# set mean temp and precip for the time step
mat_i               <- mat[i]  # mean areal temperature for time step [C]
map_i               <- map[i]  # mean areal precipitation for time step [mm]

# set binary precip form & amounts (fortran/nws ops versions have other choices)
if (mat_i <= pxtemp) { # mean areal temperature is at or below snow temp threshold
   # it's snowing (all precip assumed as snow)
   we_newsnow <- map_i
   rain <- 0
   } 
else { # otherwise
   # it's raining
   we_newsnow <- 0
   rain <- map_i
   }
 
# new snow swe [mm]
we_newsnow_gadj <- we_newsnow * scf  # (bias correct the gage(s) snow we just defined based on the snow threshold temp)

#

# if 'significant' snow fell (0.2 mm/hr, and it's the start of a new accumulation period, save last time step's swe & aesc - go off the depletion curve
# this represents ability to reset depletion trend for snow on partially bare areas (and the 0.2 mm/hr rate prevents aesc increasing to 100% 
# for very small snow amounts). last time step's values give data for setting linear equation to return to depletion curve, when swe falls below its departure value
# 

#if bias-corrected water equivalent of new snowfall is enough to depart depletion curve, and there wasn't qualifying snow last time step,
if (we_newsnow_gadj >= hsnof * tshrs_p && suffsnowprevtsteptolvdc_int == 0) { 
          swe_b4suffsnewsnow             <- swe      # remember last time step's swe (sb in fortran), and
          aesc_b4newsnow                 <- aesc     # rewmember last time step's aesc (oldaesc in fortran)
          meanarealwe_to_ai_b4newsnow    <- meanarealwe_to_ai
          aesc                           <- 1        # reset aesc to full extent
          swe                            <- swe + we_newsnow_gadj # update pre-energy exchange swe with new snowfall
}
#if bias-corrected water equivalent of new snowfall is enough to depart depletion curve, but there was qualifying snow last time step,
if (we_newsnow_gadj >= hsnof * tshrs_p && suffsnowprevtsteptolvdc_int == 1)   {      
          swe                            <- swe + we_newsnow_gadj # update pre-energy exchange swe with new snowfall
          # add new snow to last time step's swe - new swe before any heat is exchanged
          #accmax                         <- swe + we_newsnow_gadj

          
          }


## define line equation for new depletion line
#y2 <- accmax/si
#y1 <- meanarealwe_to_ai_b4newsnow
#x2 <- aesc
#x1 <- aesc_b4newsnow
#
#m   <- (y2-y1)/(x2-x1)


#
## define temp cover index
#if (we_newsnow_gadj >= hsnof * tshrs_p) {
#       swe_aescfdb100_pb  <- swe + 0.75 * we_newsnow_gadj  #[mm]
#     } 
#else {
#  swe_aescfdb100_pb  <- swe_aescfdb100_pb + 0.75 * we_newsnow_gadj } #[mm]

 
#swe <- swe + we_newsnow_gadj
max_we_ap         <- swe #[mm] #max water equivalent during accumulation period - revisit
arealindex        <- max(1.0e-100, min(max_we_ap, si)) #[mm]
meanarealwe       <- swe # [mm] #swe before any melt below  ##same as max_we_ap, right?
meanarealwe_to_ai <- min(1, meanarealwe/si) #max(0, min(1, meanarealwe/arealindex))
aesc              <- with(dt_arealdeplete, approx(meanarealwe_to_ai_x, percentarealsnow_y,  xout = meanarealwe_to_ai)) %>% 
                     tail(1) %>% unlist() %>% unname()


#aesc[i] <- aesc
 
# energy exchange at snow/air surface when no surface melt
#..
# change (increase) in the heat deficit due to new snowfall [mm] (heat amount needed to heat new snow to 0C)
# 80 cal/g: latent heat of fusion
# 0.5 cal/g/C: specific heat of ice
# new snow temperature
t_newsnow <- min(0, mat_i)
heatdefincreasefromnewsnow <- - (t_newsnow * we_newsnow_gadj)/(80/0.5)

# define/update antecedent temperature index (represents near surface snow temp from past snow & temp history),
# most recent air temps weighed by decreasing amounts
# if 'significant' new snow (>1.5 mm/hr), use new snow's temp, otherwise compute ati from last ati as variable going into representation of shallow snow temp evolution

if (we_newsnow_gadj > hsnof2 * tshrs_p) {
 ati <- t_newsnow  
} else {
 tipm_tshrs_t <- 1 - ((1 - tipm)^(tshrs_t/6))  #formula based on 6 hr timestep, accommodates any factor of 6
   ati <- ati + tipm_tshrs_t * (mat_i - ati) 
    }

ati <- min(ati, 0) #ati can't exceed 0

# define sinusoidal seasonal variation in the non-rain melt factor (assume 365 day year)
N_Mar21 <- doy[i] - 80  
sp <- (0.5 * sin((N_Mar21 * 2 * pi)/365)) + 0.5  # the seasonal pattern assumed
sp_adj <- 1  # latitudinal seasonal variation adjustment (per 'E. Anderson'06 manual, none when below ~54N)
mf <- tshrs_t/6 * ((sp * sp_adj * (mfmax - mfmin)) + mfmin)  # the seasonally varying non-rain melt factor
t_snowsurf <- min(0, mat_i)

# now with melt factor for the day of year and the snow surface temperature, get the temperature gradient (increase or decrease)
## driven change in heat deficit due to the temperature gradient within the upper part of the snowpack [mm], for only the snow covered fraction:
heatdefchangefromprofilegradient <- nmf * tshrs_p/6 * mf/mfmax * (ati - t_snowsurf) * aesc
              

# solar & atmospheric snow melt
t_rain <- max(mat_i, 0)  # rain temp is air temp as long as it's above 0C
if (rain > 0.25 * tshrs_p) {  # if 'significant' rate of 0.25 mm/hr
 # rain-on-snow melt   #assumed overcast, high humidity (>90%), emissivity of 1 at cloud elev temp, which is assumed equal to ground temp for this condition
 stefan_bolt <- 6.12 * (10^(-10))  # Stefan-Boltzman constant [mm/K/hr]
 e_sat <- 2.7489 * (10^8) * exp((-4278.63/(mat_i + 242.792)))  # saturated vapor pressure at mat_i [mb]
 p_atm <- 33.86 * (29.9 - (0.335 * (elev/100)) + (0.00022 * ((elev/100)^2.4)))  # elevation is in hundreds of meters (incorrect in snow17 manual)
 term1 <- stefan_bolt * tshrs_p * (((mat_i + 273)^4) - (273^4))
 term2 <- 0.0125 * rain * t_rain
 term3 <- 8.5 * uadj * (tshrs_p/6) * ((0.9 * e_sat - 6.11) + (0.00057 * p_atm * mat_i))
 melt_satmos <- term1 + term2 + term3
 melt_satmos <- max(melt_satmos, 0) # enforces positive melt
 melt_satmos <- melt_satmos * aesc  # only snow covered fraction can melt
} 
 else if ((rain <= 0.25 * tshrs_p) && (mat_i > mbase)) { # if insignificant rain and air temp is above snowmelt threshold (usually 0C)
 # non-rain or very little rain melt - melt factor driven and accomodates heat from small rain amounts
 melt_satmos <- (mf * (mat_i - mbase) * (tshrs_p/tshrs_t)) + (0.0125 * rain * t_rain)
 melt_satmos <- max(melt_satmos, 0) # melt can't be negative
 melt_satmos <- melt_satmos * aesc # only snow covered area can melt
} else {
 melt_satmos <- 0  #otherwise, no solar/atmospheric melt without significant rain or temps (listed above)
}

# update ice water equivalent with bias-corrected new snow amount
we_solid <- we_solid + we_newsnow_gadj   # water equivalent of total ice portion of the snow cover [mm]


# adjust heat deficit from any new snow and the evolving profile gradient from air temp & last new snow's temp history 
heatdef <- max(heatdef + heatdefincreasefromnewsnow + heatdefchangefromprofilegradient, 0)  # [mm] 

## but limit a deep snowpack's ability to keep its surface from melting
#if (heatdef >= (1/3 * we_solid)) {   
#  # set heat deficit limit
#  heatdef <- 1/3 * we_solid  #not in 2006 documentation.., check whether in more recent & whether in ops version 
#}

if (melt_satmos < we_solid) { # if solar+atmos  melt is less than the ice's water equivalent
  we_solid <- we_solid - melt_satmos # reduce ice water equivalent by the solar+atmos melt amount
  snowsurfavailliq <- melt_satmos + rain # surface liquid content is sum of solar/atmospheric melt and any rain
  liquidstorcap <- plwhc * we_solid # but the pack can retain up this much, the plwhc % of the solid water equivalent
  
  if ((snowsurfavailliq + we_liquid) > (heatdef + (heatdef * plwhc) + liquidstorcap)) {
    # if the solar+atmos melt + rain + existing liquid > heat deficit + deficit's liq stor cap + solid's liq stor cap
    # ie if there's sufficient available liquid water to overcome the total deficit & liquid storage capacity,
    # the snow is ripe:
    
    # excess liquid water is solar+atmos melt + rain + existing liquid - total deficit - liquid held by the pack
    excessliquid <- snowsurfavailliq + we_liquid - heatdef - (heatdef * plwhc) - liquidstorcap 
    
    # complete energy balance by adding in heatdef
    we_solid <- we_solid + heatdef  
    # liquid water in the pack is equal to the % maximum
    we_liquid <- liquidstorcap 
    heatdef <- 0
    }
    else if ((snowsurfavailliq + we_liquid) >= heatdef) {
    # if the solar+atmos melt + rain + existing liquid [mm] > heat deficit, but not its and liquidstorcap's extra capacity
    # the snow's not ripe 
    excessliquid <- 0 # there's no excess liquid (aka no melt and rain)
    # still increase the just-reduced we_solid, to complete heat accounting (water 'refreezes' returning pack up to 0C )
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
  
  melt_litho_liqloss    <- (daygm/we_solid) * we_liquid * aesc  
  melt_litho_solidloss  <- daygm * aesc
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
    
    
if (we_newsnow_gadj >= hsnof * tshrs_p) {
      suffsnowprevtsteptolvdc_int <- 1
    }
else {
      suffsnowprevtsteptolvdc_int <- 0
    }
    
    # update states for next time step
    ini.tstep.state <- c(we_solid, we_liquid, ati, heatdef, swe = we_solid + we_liquid, swe_b4suffsnewsnow, swe_aescfdb100_pb,
                         aesc, aesc_b4suffnewsnow, meanarealwe_to_ai,  suffsnowprevtsteptolvdc_int, accmax)
  
  }
  return(meltandrain)
  #return(aesc)

 }

{
raim_r  <- snow17(par, map, mat, elev, doy) 
raim_r  <- data.frame(raim_r, tstep) %>% transmute(tstep, raim_mm = raim_r, source = "r")
raim_icp <- read_csv("sbrc1up_icpraim.csv") %>% transmute(tstep, raim_mm = icpraim_mm, source = "icp")
df <- rbind(raim_r, raim_icp)
#p <- ggplot(df, aes(tstep, raim_mm, color = source, linetype  = source)) + geom_line()

ggplot(df, aes(tstep, raim_mm, color = source, linetype  = source)) + geom_line()
#ggplotly(p)
}
}






#Density/metamorphism stuff for later (alpmg with lag)
# new snow density [g/cm^3]
#if (mat_i <= -15) {    # snow density assumed same 'dryness'/'lightness' below ~15 C 
#  den_newsnow <- 0.05 
#}
#else {
#  den_newsnow <- 0.05 + 0.0017  * mat_i^1.5  # manual's snow density scaling increase with temperature
#}
#
## new snow depth [cm]  # centimeters for output convenience
#depth_newsnow  = (0.1 * we_newsnow_gadj)/den_newsnow

#df2 <- data.frame(meanarealwe, arealindex, meanarealwe_to_ai)
#liqlag <- 5.33 * (1 - exp  ((-0.03*(tshrs_p/6)* we_solid )/excessliquid )   )
