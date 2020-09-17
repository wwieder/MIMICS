# -----------------------------------------------------------
# MIMICS2_testing_reverseSTODEonly
# -----------------------------------------------------------
# Will Wieder
# Oct 23, 2013
# Modified Aug, 2015; April 2016; Sept 2020
# Solve a system of non-linear equations for enzyme SOC solution
# uses packages deSolve / rootSolve

# Single layer model with: 
# 2 litter C pools (LIT), corresponding to metabolic and structural litter.
# 2 microbial pool (MIC; i.e. r vs. K strategists)
# 3 SOM pools corresponding to physically & chemically protected & available pools

# -----------------------------------------------------------
# (A) Reads in site level data from LTER sites
# (B) Calculates steady state C pools using RXEQ & STODE function
# -----------------------------------------------------------

rm(list=ls())
dir <- '/Users/wwieder/Will/git_repos_local/MIMICS/'
setwd(dir)

library(rootSolve)
library(boot)

#REVERSE MODEL 
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
    OXIDAT    = ((MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + MIC_1)) +
                 (MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + MIC_2)) )  #oxidation of C to A

    dLIT_1 = I[1]*(1-FI[1]) - LITmin[1] - LITmin[3]
    dMIC_1 = CUE[1]*(LITmin[1]+ SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3])
    dSOM_1 = I[1]*FI[1] + MICtrn[1] + MICtrn[4]- DEsorb 
    
    dLIT_2 = I[2] * (1-FI[2]) - LITmin[2] - LITmin[4]
    dMIC_2 = CUE[3]*(LITmin[3]+ SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6])  
    dSOM_2 = I[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT
    
    dSOM_3 = MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]
    
    list(c(dLIT_1, dLIT_2, dMIC_1, dMIC_2, dSOM_1, dSOM_2, dSOM_3))
  })
}

#---------------------------------------------------------
# ---- Main program -------
#---------------------------------------------------------

#---------------------------------------------------------
# (A)       Read in parameters and site level data
#---------------------------------------------------------
  para_file <- paste("parameters_LIDET-MIM-REV.csv", sep = "") 
  parameters <- read.csv(para_file)
  names(parameters)
  attach(parameters)
  depth <- parameters$depth[1]
  
  litter <- read.csv("LIDET_SITE_obs/LitterCharacteristics.csv")
  litter <- litter[1:6,]   #subset foliar litter only 
  attach(litter)
  calcN    <- (1 / litCN) / 2.5 * 100    
  lit_fMET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * litLIG / calcN)   
  
  data <- read.csv("LIDET_SITE_obs/LTER_SITE_1.csv") #site level forcing variables
  names(data)
  attach(data)
  
  KO <- parameters$KO
  ANPP  <- data$ANPP / 2           		# if needed convert to gC/m2/y from g/m2/y
  clay  <- data$CLAY2/100  				    # if needed, convert from % clay to fraction
  tsoi  <- MAT
  nsites<- length(Site)
  
  lig    <- LIG #/ 100
  Nnew   <- N                                                  #N in litter additions
  fMET1  <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig / Nnew)   #as partitioned in Daycent
  MIMLIT <- rep(NA, nsites)           	                       #Vector for results
  MIMMIC <- rep(NA, nsites)           
  MIM_CO <- rep(NA, nsites)           
  MIMSOC <- rep(NA, nsites)           
  
  strSite <- as.character(data$Site)                           #convert site names to string
  pools  <- c('site','LITm', 'LITs', 'MICr', 'MICK', 'SOMp', 'SOMc', 'SOMa')
  table  <- array(NA, dim=c(nsites,8), dimnames=list(as.character(Site),pools)) 
  table[,1]       <- as.character(Site)
  
  POOLS  <- c("LIT","MIC","SOC")
  npools <- length(POOLS)  
  LITpool<- c('LITm', 'LITs') 
  MICpool<- c('MICr', 'MICk') 
  SOMpool<- c('SOMp', 'SOMa','SOMc') 
  
  MIMLIT <- array(NA, dim=c(nsites))           	             # array for results
  MIMMIC <- array(NA, dim=c(nsites))           
  MIM_CO <- array(NA, dim=c(nsites))           
  MIMSOC <- array(NA, dim=c(nsites))           
  
  strSite <- as.character(data$Site)  #convert site names to string
  
  #Make vectors to store model results that match obs temporal resolution
  npts   <- 6*10*14   				#6 litter * 10 years * 14 sites
  xyLIT  <- rep(NA, npts) 
  xyTIME <- rep(NA, npts) 
  xySITE <- rep(NA, npts) 
  xyOBS  <- rep(NA, npts)
  xyMIM  <- rep(NA, npts)
  xyCount<- 1
  
  #-----------------------------------------------------------
  # (B)       RXEQ for site using STODE function
  # Starts Big loop over all sites (i) & all litter types (j)
  #-----------------------------------------------------------
  
  for (i in 1:nsites) {         #speeds up debugging     	
    # Read in site characteristics -----------
    print(paste("-------- starting ", data$Site[i], " --------") )
    
    fMET       <- mean(lit_fMET)           # uses mean litter fmet from LIDET
    # fMET       <- fMET1[i]                 # uses site estimate for fmet
    fCLAY      <- clay[i]
    TSOI       <- tsoi[i]
    EST_LIT_in <- ANPP[i] / (365*24)         # gC/m2/h (from gC/m2/y)
    h2y        <- 24*365
    MICROtoECO <- depth * 1e4 * 1e-3         # mgC/cm3 to g/m2
    EST_LIT    <- EST_LIT_in  * 1e3 / 1e4    # mgC/cm2/h(from gC/m2/h) 
    
    # ------------ caclulate parameters ---------------
    Vmax     <- exp(TSOI * Vslope + Vint) * aV
    Km       <- exp(TSOI * Kslope + Kint) * aK
    
    #ANPP strongly correlated with MAP
    Tau_MOD1 <- sqrt(ANPP[i]/Tau_MOD[1])          # basicaily standardize against NWT
    Tau_MOD2 <- Tau_MOD[4]                        # increased 3-fold for SS SOC pools
    Tau_MOD1[Tau_MOD1 < Tau_MOD[2]] <- Tau_MOD[2] # correction not used in LIDET resutls 
    Tau_MOD1[Tau_MOD1 > Tau_MOD[3]] <- Tau_MOD[3] 
    tau <- c(tau_r[1]*exp(tau_r[2]*fMET), 
             tau_K[1]*exp(tau_K[2]*fMET))   
    tau <- tau * Tau_MOD1 * Tau_MOD2
    
    fPHYS    <- c(fPHYS_r[1] * exp(fPHYS_r[2]*fCLAY), 
                  fPHYS_K[1] * exp(fPHYS_K[2]*fCLAY)) 	            #fraction to SOMp
    fCHEM    <- c(fCHEM_r[1] * exp(fCHEM_r[2]*fMET) * fCHEM_r[3], 
                  fCHEM_K[1] * exp(fCHEM_K[2]*fMET) * fCHEM_K[3]) 	#fraction to SOMc
    fAVAI    <- 1- (fPHYS + fCHEM)
    desorb   <- fSOM_p[1] * exp(fSOM_p[2]*(fCLAY))                  #CHANGED FOR GLOBAL RUN!!!     
    
    #desorb   <- desorb/10 # modified as in MIMdef from Zhang et al 2020
    #fPHYS    <- fPHYS/5  # to reduce allocation to physically protected pool 5x
    
    pSCALAR  <- PHYS_scalar[1] * exp(PHYS_scalar[2]*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
    v_MOD    <- vMOD  # to avoid writing over orig. parameters
    k_MOD    <- kMOD 
    k_MOD[3] <- k_MOD[3] * pSCALAR    
    k_MOD[6] <- k_MOD[6] * pSCALAR    
    
    VMAX     <- Vmax * v_MOD 
    KM       <- Km / k_MOD
    
    #initialize pools
    I       <- array(NA, dim=2)              #Litter inputs to MET/STR
    I[1]    <- (EST_LIT / depth) * fMET     
    I[2]    <- (EST_LIT / depth) * (1-fMET)
    lit     <- I   
    mic     <- I  
    som     <- rep(NA, 3) 
    som[1]  <- I[1]
    som[2]  <- I[2]
    som[3]  <- I[1] 
    LITmin  <- rep(NA, dim=4)
    MICtrn  <- rep(NA, dim=6)
    SOMmin  <- rep(NA, dim=2)
    DEsorb  <- rep(NA, dim=1)
    OXIDAT  <- rep(NA, dim=1)
    
    #Calculate RXEQ pools  
    Tpars <- c( I = I, VMAX = VMAX, KM = KM, CUE = CUE, 
                fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
                tau = tau, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
                desorb = desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO)
    Ty    <- c( LIT_1 = lit[1], LIT_2 = lit[2], 
                MIC_1 = mic[1], MIC_2 = mic[2], 
                SOM_1 = som[1], SOM_2 = som[2], SOM_3 = som[3] )
    test  <- stode(y = Ty, time = 1e6, fun = RXEQ, parms = Tpars, positive = TRUE)
    if (i == 2) {print(test[1])}
    remove(lit, mic, som)
  
    table[i,2:8] <- as.numeric(test[[1]])
    MIMLIT[i]    <- (test[[1]][[1]]+test[[1]][[2]])  * depth *1e4 / 1e6#convert kgC/m2 from mgC/cm3 (0-30 cm) 
    MIMMIC[i]    <- (test[[1]][[3]]+test[[1]][[4]])  * depth *1e4 / 1e6
    MIM_CO[i]    <-  test[[1]][[3]]/test[[1]][[4]]
    MIMSOC[i]    <- sum(test[[1]])  * depth *1e4 / 1e6   
    
    remove(test, Ty, Tpars, LITmin, MICtrn, SOMmin)
  
  } # close i loop (sites)
  
  r <- cor.test(data$SOC, MIMSOC)
  cor1text <- paste("r = ", signif(r[[4]][[1]], digits=2))
  par(mar=c(5,5,2,1))
  plot(MIMSOC, data$SOC,ylim=c(0,10), xlim=c(0,10), col=0,
       pch=16, cex.lab=1.2, cex.axis = 1.2, cex=1.8,
       main=paste0('KO = ',KO[1]),
       ylab=expression(paste("Observed SOC (kg C ", m^-2, ")")),
       xlab=expression(paste("Predicted SOC (kg C ", m^-2, ")")))
  legend("topleft", pch = 16, cex=1.1, col=0, pt.cex=0.9, legend=cor1text,bty="n")

abline(0,1, lty=2)
text(10,9.5,"1:1")
text(MIMSOC, data$SOC, labels = Site,col=1)

