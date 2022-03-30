# Will Wieder [wwieder@ucar.edu]
# sets up lab incubation experiment, 
# Feb 2022

# !! this code base has not yet been published !!
# Please check with W. Wieder before distributing this code
# Exlolores how we might replicate lab incubation for Macrosystem Biology project.


rm(list=ls())
library(rootSolve)

#----------------Define function for fluxes-----------------------
# fPHYS and (fCHEM?) are set to zero using parameter file
# TODO add water scalar function, only uses temperature and litter quality
# Currently model assumes historic litter quality for fMET,
# Could also consider how historic mositure or temperature could constrain microbial physiology
# Extra pools and fluxes retained here for capatability.

# Define the REVERSE MODEL 
RXEQ <- function(t, y, pars) {
  with (as.list(c(y, pars)),{
    
    #Flows to and from MIC_1
    LITmin[1] = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + MIC_1)   #MIC_1 decomp of MET lit
    LITmin[2] = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + MIC_1)   #MIC_1 decomp of STRUC lit
    MICtrn[1] = MIC_1 * tau[1]  * fPHYS[1]                  #MIC_1 turnover to PHYSICAL SOM 
    MICtrn[2] = MIC_1 * tau[1]  * fCHEM[1]                  #MIC_1 turnover to CHEMICAL SOM  
    MICtrn[3] = MIC_1 * tau[1]  * fAVAI[1]                  #MIC_1 turnover to AVAILABLE SOM  
    SOMmin[1] = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + MIC_1)   #decomp of SOMa by MIC_1
    
    #Flows to and from MIC_2
    LITmin[3] = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + MIC_2)   #decomp of MET litter
    LITmin[4] = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + MIC_2)   #decomp of SRUCTURAL litter
    MICtrn[4] = MIC_2 * tau[2]  * fPHYS[2]                  #MIC_2 turnover to PHYSICAL  SOM 
    MICtrn[5] = MIC_2 * tau[2]  * fCHEM[2]                  #MIC_2 turnover to CHEMICAL  SOM  
    MICtrn[6] = MIC_2 * tau[2]  * fAVAI[2]                  #MIC_2 turnover to AVAILABLE SOM  
    SOMmin[2] = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + MIC_2)   #decomp of SOMa by MIC_2
    
    DEsorb    = SOM_1 * desorb  #* (MIC_1 + MIC_2)  	#desorbtion of PHYS to AVAIL (function of fCLAY)
    OXIDAT    = ((MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + MIC_2)) +
                   (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + MIC_1)))  #oxidation of C to A

    dLIT_1 = I[1]*(1-FI[1]) - LITmin[1] - LITmin[3]
    dMIC_1 = CUE[1]*(LITmin[1]+ SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3])
    dSOM_1 = I[1]*FI[1] + MICtrn[1] + MICtrn[4]- DEsorb 
    CO2_1  = (1-CUE[1])*(LITmin[1]+ SOMmin[1]) + (1-CUE[2])*(LITmin[2])  # CO2 fluxes from MICr

    dLIT_2 = I[2] * (1-FI[2]) - LITmin[2] - LITmin[4]
    dMIC_2 = CUE[3]*(LITmin[3]+ SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6])  
    dSOM_2 = I[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT
    CO2_2  = (1-CUE[3])*(LITmin[3]+ SOMmin[2]) + (1-CUE[4])*(LITmin[4])  # CO2 fluxes from MICk
    
    dSOM_3 = MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]
    
    list(c(dLIT_1, dLIT_2, dMIC_1, dMIC_2, dSOM_1, dSOM_2, dSOM_3,CO2_1,CO2_2))
  })
}

# Define the REVERSE MODEL, here for a daily timestep

RXEQ_day <- function(t, y, pars) {
  with (as.list(c(y, pars)),{
    
    #Flows to and from MIC_1
    LITmin[1] = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + MIC_1)   #MIC_1 decomp of MET lit
    LITmin[2] = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + MIC_1)   #MIC_1 decomp of STRUC lit
    MICtrn[1] = MIC_1 * tau[1]  * fPHYS[1]                  #MIC_1 turnover to PHYSICAL SOM 
    MICtrn[2] = MIC_1 * tau[1]  * fCHEM[1]                  #MIC_1 turnover to CHEMICAL SOM  
    MICtrn[3] = MIC_1 * tau[1]  * fAVAI[1]                  #MIC_1 turnover to AVAILABLE SOM  
    SOMmin[1] = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + MIC_1)   #decomp of SOMa by MIC_1
    
    #Flows to and from MIC_2
    LITmin[3] = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + MIC_2)   #decomp of MET litter
    LITmin[4] = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + MIC_2)   #decomp of SRUCTURAL litter
    MICtrn[4] = MIC_2 * tau[2]  * fPHYS[2]                  #MIC_2 turnover to PHYSICAL  SOM 
    MICtrn[5] = MIC_2 * tau[2]  * fCHEM[2]                  #MIC_2 turnover to CHEMICAL  SOM  
    MICtrn[6] = MIC_2 * tau[2]  * fAVAI[2]                  #MIC_2 turnover to AVAILABLE SOM  
    SOMmin[2] = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + MIC_2)   #decomp of SOMa by MIC_2
    
    DEsorb    = SOM_1 * desorb  #* (MIC_1 + MIC_2)  	#desorbtion of PHYS to AVAIL (function of fCLAY)
    OXIDAT    = ((MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + MIC_2)) +
                   (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + MIC_1)))  #oxidation of C to A
    
    dLIT_1 = 24*(I[1]*(1-FI[1]) - LITmin[1] - LITmin[3] )
    dMIC_1 = 24*(CUE[1]*(LITmin[1]+ SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3]))
    dSOM_1 = 24*(I[1]*FI[1] + MICtrn[1] + MICtrn[4]- DEsorb )
    CO2_1  = 24*((1-CUE[1])*(LITmin[1]+ SOMmin[1]) + (1-CUE[2])*(LITmin[2]) ) # CO2 fluxes from MICr
    
    dLIT_2 = 24*(I[2] * (1-FI[2]) - LITmin[2] - LITmin[4])
    dMIC_2 = 24*(CUE[3]*(LITmin[3]+ SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6])  )
    dSOM_2 = 24*(I[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT)
    CO2_2  = 24*((1-CUE[3])*(LITmin[3]+ SOMmin[2]) + (1-CUE[4])*(LITmin[4])  )# CO2 fluxes from MICk
    
    dSOM_3 = 24*(MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2])
    
    list(c(dLIT_1, dLIT_2, dMIC_1, dMIC_2, dSOM_1, dSOM_2, dSOM_3,CO2_1,CO2_2))
  })
}

#---------------------------------------------------------
# (A)       Read in parameters and site level data
#---------------------------------------------------------
para_file <- paste("parameters_LIDET-MIM-REV_LabInc.csv", sep = "") 
parameters <- read.csv(para_file)
names(parameters)
attach(parameters)

# multiplied default vMOD by 0.3 
data <- read.csv("NEON_SITE_1.csv") #site level forcing variables
names(data)
attach(data)

ANPP    <- ANPP / 2       			# convert to gC/m2/y from g/m2/y
clay   <- CLAY2/100  				    # convert from clay fraction to %
nsites <- length(Site)

lig    <- LIG/100
Nnew   <- 1/CN/2.5                  	                       # N in litter additions
fMET1  <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig / Nnew)   # as partitioned in Daycent

## Create table for summary output
exp = c('control')
nexp = length(exp)
strSite <- as.character(data$Site)  #convert site names to string
cols  <- c('site','LITm', 'LITs', 'MICr', 'MICK', 'SOMp', 'SOMc', 'SOMa','CO2')
ncols = length(cols)
table  <- array(NA, dim=c(nexp, nsites,ncols), dimnames=list(exp,as.character(Site),cols)) 
table[,,1]       <- as.character(Site)

table
# ----------------------------------------------------------
# (B)     start lab incubation experiment
# ----------------------------------------------------------

## read in and calculate parameters 
## Not all of these will end up being used for this lab experiment

#TODO, loop over sites and experiments
i <- 1
fMET    <- fMET1[i] 
fCLAY   <- clay[i]
mat     <- MAT[i]
map     <- MAP[i]
## set temperature for controled lab experiment
TSOI  <- 25
theta <- 40

EST_LIT_in <- ANPP[i] / (365*24) # gC/m2/h (from gC/m2/y)
depth      <- parameters$Depth[1] 
h2y        <- 24*365
MICROtoECO <- depth * 1e4 * 1e6 / 1e6    #mgC/cm3 to kgC/km2

#-----------------caclulate parameters---------------------------
Vmax     <- exp(Vint  + TSOI * Vslope ) * aV
Km       <- exp(Kint  + TSOI * Kslope ) * aK

# density edependent microbial turnover, based on NPP 
Tau_MOD1 <- sqrt(ANPP[i]/Tau_MOD[1])          
Tau_MOD2 <- Tau_MOD[4]                         
Tau_MOD1[Tau_MOD1 < Tau_MOD[2]] <- Tau_MOD[2] 
Tau_MOD1[Tau_MOD1 > Tau_MOD[3]] <- Tau_MOD[3] 
tau <- c(tau_r[1]*exp(tau_r[2]*fMET), 
         tau_K[1]*exp(tau_K[2]*fMET))   
tau <- tau * Tau_MOD1 * Tau_MOD2

# -- allocation of microbial residues to SOM pools --- 
fPHYS    <- c(fPHYS_r[1] * exp(fPHYS_r[2]*fCLAY), 
              fPHYS_K[1] * exp(fPHYS_K[2]*fCLAY)) 	            #fraction to SOMp
fCHEM    <- c(fCHEM_r[1] * exp(fCHEM_r[2]*fMET) * fCHEM_r[3], 
              fCHEM_K[1] * exp(fCHEM_K[2]*fMET) * fCHEM_K[3]) 	#fraction to SOMc
fAVAI    <- 1- (fPHYS + fCHEM)
desorb   <- fSOM_p[1] * exp(fSOM_p[2]*(fCLAY))                 #CHANGED FOR GLOBAL RUN!!!     
pSCALAR  <- PHYS_scalar[1] * exp(PHYS_scalar[2]*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp

#------------MODIFY fluxes from SOM pools as function of clay content---------
v_MOD    <- vMOD
k_MOD    <- kMOD # to avoid writing over orig. parameters
k_MOD[3] <- k_MOD[3] * pSCALAR    
k_MOD[6] <- k_MOD[6] * pSCALAR    

VMAX     <- Vmax * v_MOD 
KM       <- Km / k_MOD

LITmin  <- rep(NA, dim=4)
MICtrn  <- rep(NA, dim=6)
SOMmin  <- rep(NA, dim=2)
DEsorb  <- rep(NA, dim=1)
OXIDAT  <- rep(NA, dim=1)


nday   <- 200
day    <- 1


#initialize arrays to store daily output data
r.names = c('LITm','LITs')
LIT    <- array(NA, dim = c(2,nday), dimnames = list(r.names,NULL))
r.names = c('MICr','MICk')
MIC    <- array(NA, dim = c(2,nday), dimnames = list(r.names,NULL))
CO2    <- array(NA, dim = c(2,nday), dimnames = list(r.names,NULL))
r.names = c('SOMp','SOMc','SOMa')
SOM    <- array(NA, dim = c(3,nday), dimnames = list(r.names,NULL))

#initialize pools and fluxes
I        <- rep(0,2)
LIT_1    <- 100   
LIT_2    <- 100
MIC_1    <- 0.01
MIC_2    <- 0.01
SOM_1    <- 0
SOM_2    <- 0
SOM_3    <- 0
CO2_1    <- 0
CO2_2    <- 0

# Loop over days
for (d in 1:nday)  {
  for (h in 1:24)   {
    
    UPpars  <- c( I = I, VMAX = VMAX, KM = KM, CUE = CUE, 
                  fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
                  tau   = tau, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
                  desorb= desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO)
    UPy     <- c( LIT_1 = LIT_1, LIT_2 = LIT_2, 
                  MIC_1 = MIC_1, MIC_2 = MIC_2, 
                  SOM_1 = SOM_1, SOM_2 = SOM_2, SOM_3 = SOM_3 )
    update <- RXEQ(y = UPy, pars = UPpars)

    LIT_1  <- LIT_1 + update[[1]][1]
    LIT_2  <- LIT_2 + update[[1]][2]
    MIC_1  <- MIC_1 + update[[1]][3]
    MIC_2  <- MIC_2 + update[[1]][4]
    SOM_1  <- SOM_1 + update[[1]][5]
    SOM_2  <- SOM_2 + update[[1]][6]
    SOM_3  <- SOM_3 + update[[1]][7]
    CO2_1  <- CO2_1 + update[[1]][8]
    CO2_2  <- CO2_2 + update[[1]][9]
    remove(UPpars, UPy, update)
    
    #write out daily results
    if (h == 24) {
      LIT[1,d] <- LIT_1
      LIT[2,d] <- LIT_2
      MIC[1,d] <- MIC_1
      MIC[2,d] <- MIC_2
      SOM[1,d] <- SOM_1
      SOM[2,d] <- SOM_2
      SOM[3,d] <- SOM_3
      CO2[1,d] <- CO2_1
      CO2[2,d] <- CO2_2
    }	   						#close daily results counter
  }							#close hour loop
}								#close daily loop


par(mfrow=c(3,1),mar=c(0,4,3,2))
maxy = max(LIT)*1.01
miny = min(LIT)*0.9
x = seq(1,nday)
plot(x,LIT[1,],col=1,type='l',lwd=2,
     xlab=NA,ylab='Litter C Remaining', xaxt='n',
     ylim=c(miny,maxy) ,main='Hourly function' )
at1 <- seq(0, 200, 50)
axis(side =1, at1, labels = F)
lines(x,LIT[2,],lwd=2,col=2)
legend(-8,(maxy+miny)/2,c('LITs','LITm'),
       text.col=c(2,1),text.font = (cex=2),
       bty = "n")


par(mar=c(1.5,4,1.5,2))
maxy = max(c(MIC))*1.1
plot(x,MIC[1,],col=1,lwd=2,type='l',
     xlab=NA,ylab='Microbial and Soil C',
     ylim=c(1e-2,maxy), xaxt='n')
at1 <- seq(0, 200, 50)
axis(side =1, at1, labels = F)
lines(x,MIC[2,],lwd=2,col=2)
lines(x,SOM[3,],lwd=2,col=4)
lines(x,SOM[2,],lwd=2,col=5)
legend(-8,maxy,c('SOMc','SOMa','MICk','MICr'),
       text.col=c(5,4,2,1),text.font = (cex=2),
       bty = "n")

par(mar=c(3,4,0,2))
totalC = colSums(LIT)+colSums(MIC)+colSums(SOM)
initC = totalC[1]
CO2_calc = (initC - totalC) / initC
CO2_mod = colSums(CO2) / initC
plot(x,CO2_calc,col=1,lwd=2,type='l',
     xlab='day',ylab='CO2 (fraction of initial)' )
lines(x,CO2_mod,lwd=2,col=4)

normRate = CO2_mod/colSums(LIT)[1]/(x*24)
plot(x,normRate)


# -----------------
# Repeat with daily function
# -----------------


#initialize arrays to store daily output data
r.names = c('LITm','LITs')
LIT    <- array(NA, dim = c(2,nday), dimnames = list(r.names,NULL))
r.names = c('MICr','MICk')
MIC    <- array(NA, dim = c(2,nday), dimnames = list(r.names,NULL))
CO2    <- array(NA, dim = c(2,nday), dimnames = list(r.names,NULL))
r.names = c('SOMp','SOMc','SOMa')
SOM    <- array(NA, dim = c(3,nday), dimnames = list(r.names,NULL))

#initialize pools and fluxes
I        <- rep(0,2)
LIT_1    <- 100   
LIT_2    <- 100
MIC_1    <- 0.01
MIC_2    <- 0.01
SOM_1    <- 0
SOM_2    <- 0
SOM_3    <- 0
CO2_1    <- 0
CO2_2    <- 0

# Loop over days
for (d in 1:nday)  {

    UPpars  <- c( I = I, VMAX = VMAX, KM = KM, CUE = CUE, 
                  fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
                  tau   = tau, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
                  desorb= desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO)
    UPy     <- c( LIT_1 = LIT_1, LIT_2 = LIT_2, 
                  MIC_1 = MIC_1, MIC_2 = MIC_2, 
                  SOM_1 = SOM_1, SOM_2 = SOM_2, SOM_3 = SOM_3 )
    update <- RXEQ_day(y = UPy, pars = UPpars)
    
    LIT_1  <- LIT_1 + update[[1]][1]
    LIT_2  <- LIT_2 + update[[1]][2]
    MIC_1  <- MIC_1 + update[[1]][3]
    MIC_2  <- MIC_2 + update[[1]][4]
    SOM_1  <- SOM_1 + update[[1]][5]
    SOM_2  <- SOM_2 + update[[1]][6]
    SOM_3  <- SOM_3 + update[[1]][7]
    CO2_1  <- CO2_1 + update[[1]][8]
    CO2_2  <- CO2_2 + update[[1]][9]
    remove(UPpars, UPy, update)
    
    #write out daily results
    LIT[1,d] <- LIT_1
    LIT[2,d] <- LIT_2
    MIC[1,d] <- MIC_1
    MIC[2,d] <- MIC_2
    SOM[1,d] <- SOM_1
    SOM[2,d] <- SOM_2
    SOM[3,d] <- SOM_3
    CO2[1,d] <- CO2_1
    CO2[2,d] <- CO2_2
}								#close daily loop


par(mfrow=c(3,1),mar=c(0,4,3,2))
maxy = max(LIT)*1.01
miny = min(LIT)*0.9
x = seq(1,nday)
plot(x,LIT[1,],col=1,type='l',lwd=2,
     xlab=NA,ylab='Litter C Remaining', xaxt='n',
     ylim=c(miny,maxy),main='Daily function' )
at1 <- seq(0, 200, 50)
axis(side =1, at1, labels = F)
lines(x,LIT[2,],lwd=2,col=2)
legend(-8,(maxy+miny)/2,c('LITs','LITm'),
       text.col=c(2,1),text.font = (cex=2),
       bty = "n")


par(mar=c(1.5,4,1.5,2))
maxy = max(c(MIC))*1.1
plot(x,MIC[1,],col=1,lwd=2,type='l',
     xlab=NA,ylab='Microbial and Soil C',
     ylim=c(1e-2,maxy), xaxt='n')
at1 <- seq(0, 200, 50)
axis(side =1, at1, labels = F)
lines(x,MIC[2,],lwd=2,col=2)
lines(x,SOM[3,],lwd=2,col=4)
lines(x,SOM[2,],lwd=2,col=5)
legend(-8,maxy,c('SOMc','SOMa','MICk','MICr'),
       text.col=c(5,4,2,1),text.font = (cex=2),
       bty = "n")

par(mar=c(3,4,0,2))
totalC = colSums(LIT)+colSums(MIC)+colSums(SOM)
initC = totalC[1]
CO2_calc = (initC - totalC) / initC
CO2_mod = colSums(CO2) / initC
plot(x,CO2_calc,col=1,lwd=2,type='l',
     xlab='day',ylab='CO2 (fraction of initial)' )
lines(x,CO2_mod,lwd=2,col=4)

normRate = CO2_mod/colSums(LIT)[1]/(x*24)
plot(x,normRate)


