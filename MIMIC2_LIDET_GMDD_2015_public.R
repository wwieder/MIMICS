# Will Wieder
# Generates data for Fig. 2 of 
# Wieder et al. Geosci. Model Dev. Discuss., 8, 2011–2052, 2015  doi:10.5194/gmdd-8-2011-2015
# Data or findings that are based on this model, please cite the manuscript above

# Created Aug 16 2013 (Modified May 2014)
# Solve a system of non-linear equations for enzyme SOC solution
# uses packages deSolve / rootSolve
# compares results from LIDET / LTER sites
# MODIFIED TO USE MEAN LITTER CHEMISTRY FROM LIDET LITTER
# seperate file ending in "_test.R" USES SITE LEVEL LITTER CHEMISTIRY TO SPIN UP
# Here tao modified by site level productivity (sqrt(ANPP/100))


# Single layer model with: 
  # 2 litter C pools (LIT), corresponding to metabolic and structural litter.
  # 2 microbial pool (MIC; i.e. R vs. K strategists)
  # 3 SOM pools corresponding to physically & chemically protected & available pools
# ADDS COMPLEXITY TO LITTER DECOMPOSITION

# (A) Reads in site level data from LTER sites (Table C1; data in seperate .csv file)
# (B) Calculates steady state C pools using XEQ & STODE function
# (C) runs for another 20 years with mean soil temperature 
# (D) Replicates LIDET experiment, adding litter Oct 1 of year 0

library(rootSolve)

#----------------analytical solutin using stode function-----------------------

XEQ <- function(t, y, pars) {
    with (as.list(c(y, pars)),{

    #Flows to and from MIC_1
    LITmin[1] = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + LIT_1)   #MIC_1 decomp of MET lit
    LITmin[2] = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + LIT_2)   #MIC_1 decomp of STRUC lit
    MICtrn[1] = MIC_1 * tao[1]  * fPHYS[1]                  #MIC_1 turnover to PHYSICAL SOM 
    MICtrn[2] = MIC_1 * tao[1]  * fCHEM[1]                  #MIC_1 turnover to CHEMICAL SOM  
    MICtrn[3] = MIC_1 * tao[1]  * fAVAI[1]                  #MIC_1 turnover to AVAILABLE SOM  
    SOMmin[1] = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + SOM_3)   #decomp of SOMa by MIC_1

    #Flows to and from MIC_2
    LITmin[3] = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + LIT_1)   #decomp of MET litter
    LITmin[4] = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + LIT_2)   #decomp of SRUCTURAL litter
    MICtrn[4] = MIC_2 * tao[2]  * fPHYS[2]                  #MIC_2 turnover to PHYSICAL  SOM 
    MICtrn[5] = MIC_2 * tao[2]  * fCHEM[2]                  #MIC_2 turnover to CHEMICAL  SOM  
    MICtrn[6] = MIC_2 * tao[2]  * fAVAI[2]                  #MIC_2 turnover to AVAILABLE SOM  
    SOMmin[2] = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + SOM_3)   #decomp of SOMa by MIC_2
	
	DEsorb    = SOM_1 * desorb  #* (MIC_1 + MIC_2)		#desorbtion of PHYS to AVAIL (function of fCLAY)
	OXIDAT    = ((MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + SOM_2)) +
		          (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + SOM_2)))  #oxidation of C to A
		#can make fluxes from CHEM a function of microbial biomass size?

    dLIT_1 = I[1]*(1-FI[1]) - LITmin[1] - LITmin[3]
    dMIC_1 = CUE[1]*(LITmin[1]+ SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3])
    dSOM_1 = I[1]*FI[1] + MICtrn[1] + MICtrn[4]- DEsorb 
	    
    dLIT_2 = I[2] * (1-FI[2]) - LITmin[2] - LITmin[4]
    dMIC_2 = CUE[3]*(LITmin[3]+ SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6])  
    dSOM_2 = I[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT

	dSOM_3  = MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]

    list(c(dLIT_1, dLIT_2, dMIC_1, dMIC_2, dSOM_1, dSOM_2, dSOM_3))
    })
  }

#---------------------------------------------------------
# (A)       Read in site level data
#---------------------------------------------------------
data <- read.csv("LTER_SITE_1.csv") #site level forcing variables
names(data)
attach(data)

ANPP  	<- ANPP / 2       			# convert to gC/m2/y from g/m2/y
strSite <- as.character(data$Site)  #convert site names to string
nsites  <- length(strSite)
npts    <- 6*10*14 					#6 litter * 10 years * 14 sites

#Make vectors to store model results that match obs temporal resolution
xyLIT  <- rep(NA, npts) 
xyTIME <- rep(NA, npts) 
xySITE <- rep(NA, npts) 
xyOBS  <- rep(NA, npts)
xyMIM  <- rep(NA, npts)
xyCount<- 1

#for (s in 1:nsites) {
for (s in 1:2) {
	
  clay   <- CLAY2/100					#convert from clay fraction to %
  tsoi   <- MAT
  nsites <- length(Site)

  lig    <- LIG/100
  Nnew   <- 1/CN/2.5                  	#N in litter additions
  fMET1  <- 0.85 - 0.013 * lig / Nnew   #as partitioned in Daycent
  MIMLIT <- rep(NA, nsites)           	#Vector for results
  MIMMIC <- rep(NA, nsites)           
  MIMSOC <- rep(NA, nsites)           

  #---------------------------------------------------------
  # (B)       XEQ for site using STODE function
  #---------------------------------------------------------

  LITtype  <- c('TRAEf', 'PIREf','THPLf','ACSAf','QUPRf','DRGLf')
  bagMET   <- c(10.6, 36.2, 37.4, 56.8, 37.1, 49.3) #from Gordon's LitterCharacteristics.txt
  bagLIG   <- c(16.2, 19.2, 26.7, 15.9, 23.5, 10.9) # % from Gordon's LitterCharacteristics.txt
  bagN     <- c(0.38, 0.59, 0.62, 0.81, 1.03, 1.97) # %N 
  bagCN    <- c(133.3,92.7, 83.1, 61.8, 50.5, 24.2)
  calcN    <- (1 / bagCN) / 2.5 * 100    
  calcMET  <- 0.85 - 0.013 * bagLIG/calcN 			#as calculated in DAYCENT
  bagMET   <- bagMET / 100
  bagMET   <- calcMET

  #bagMET[1]<- mean(bagMET)  #SPEEDS search

  fMET        <- mean(calcMET)
  TSOI        <- tsoi[s]   
  fCLAY       <- clay[s]
  EST_LIT_in  <- ANPP[s] / (365*24)   		#gC/m2/h (from g/m2/y, Knapp et al. Science 2001)
  BAG_LIT_in  <- 100      					#gC/m2/h
  depth       <- 30
  h2y         <- 24*365
  MICROtoECO  <- depth * 1e4 * 1e6 / 1e6   	#mgC/cm3 to kgC/km2
  EST_LIT     <- EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/h 
  BAG_LIT     <- BAG_LIT_in  * 1e3 / 1e4    #mgC/cm2/h

  #-----------------caclulate parameters---------------------------
  #Calculate Vmax & (using parameters from German 2012, as in Wieder et al. 2013 Nature Climate Change)
  Vslope   <- array(0.063,dim=6)
  Vint     <- 5.47
  aV       <- 8e-6
  Vmax     <- exp(TSOI * Vslope + Vint) * aV

  Kslope   <- array(NA,dim=6)
  Kslope[1]<- 0.017 #META LIT to MIC_1
  Kslope[2]<- 0.027 #STRU LIT to MIC_1 
  Kslope[3]<- 0.017 #AVAI SOM to MIC_1 
  Kslope[4]<- 0.017 #META LIT to MIC_2
  Kslope[5]<- 0.027 #STRU LIT to MIC_2
  Kslope[6]<- 0.017 #AVAI SOM to MIC_2
  Kint     <- 3.19
  aK       <- 10
  Km       <- exp(Kslope * TSOI + Kint) * aK

  	CUE        <- c(0.55, 0.25, 0.75, 0.35)  #for LITm and LITs entering MICr and MICK, respectively
	#ANPP strongly correlated with MAP
	Tao_MOD1 <- sqrt(ANPP[s]/100)  #basicaily standardize against NWT
	tao      <- c(5.2e-4*exp(0.3*fMET), 2.4e-4*exp(0.1*fMET))	
    tao      <- tao * Tao_MOD1
    
	#------NEW Parameters--------------
	fPHYS    <- c(0.3 * exp(1.3*fCLAY), 0.2 * exp(0.8*fCLAY)) 	#fraction to SOMp
	fCHEM    <- c(0.1 * exp(-3*fMET)  , 0.3 * exp(-3*fMET)  ) 	#fraction to SOMc
	fAVAI    <- 1- (fPHYS + fCHEM)
	desorb   <- 9e-4 * exp(-3*(sqrt(fCLAY))) #if modified by MIC!
	desorb   <- 3e-4 * exp(-4*(sqrt(fCLAY))) #if stand alone rate constant
    desorb   <- 1.5e-5 * exp(-1.5*(fCLAY))      #CHANGED FOR GLOBAL RUN!!!   
	
	k        <- 2.0    #2.0			#REDUCED FROM 3 TO 1, REDUCES TEXTURE EFFECTS ON SOMa decay
	a        <- 2.0    #2.2			#increased from 4.0 to 4.5
	
	cMAX     <- 1.4                    #ORIG 1.4 Maximum CHEM SOM scalar w/   0% Clay 
	cMIN     <- 1.2                    #ORIG 1.4 Minimum CHEM SOM scalar w/ 100% Clay 
	cSLOPE   <- cMIN - cMAX            #Slope of linear function of cSCALAR for CHEM SOM  

	pSCALAR  <- a * exp(-k*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp

#------------!!MODIFIERS AS IN MIMICS2_b!!---------------
	MOD1     <- c(10, 2, 10, 3, 3, 2) 
	MOD2     <- c( 8, 2 ,4 * pSCALAR, 2, 4, 6 * pSCALAR) 	

	VMAX     <- Vmax * MOD1 
	KM       <- Km / MOD2
	KO       <- c(4,4)      #scalar modifies Km of Oxidat	
	I        <- array(NA, dim=2)              #Litter inputs to MET/STR
	I[1]     <- (EST_LIT / depth) * fMET      #partitioned to layers
	I[2]     <- (EST_LIT / depth) * (1-fMET)

	BAG      <- array(NA, dim=c(6,2))              #litter BAG inputs to MET/STR

	for (i in 1:6) {
	  BAG[i,1]   <- (BAG_LIT / depth) * bagMET[i]      #partitioned to layers
	  BAG[i,2]   <- (BAG_LIT / depth) * (1-bagMET[i])
	}

	FI       <- c(0.05, 0.05)
	#initialize pools
	LIT     <- I   # * 1e3
	MIC     <- I   # * 1e2
	SOM     <- rep(NA, 3) 
	SOM[1]  <- I[1]
	SOM[2]  <- I[2]
	SOM[3]  <- I[1] 
	LITmin  <- rep(NA, dim=4)
	MICtrn  <- rep(NA, dim=6)
	SOMmin  <- rep(NA, dim=2)
	DEsorb  <- rep(NA, dim=1)
	OXIDAT  <- rep(NA, dim=1)

    #Calculate XEQ pools

	Tpars <- c( I = I, VMAX = VMAX, KM = KM, CUE = CUE, 
				fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
    	        tao = tao, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
    	        desorb = desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO)
	Ty    <- c( LIT_1 = LIT[1], LIT_2 = LIT[2], 
    	        MIC_1 = MIC[1], MIC_2 = MIC[2], 
        	    SOM_1 = SOM[1], SOM_2 = SOM[2], SOM_3 = SOM[3] )
	test  <- stode(y = Ty, time = 1e6, fun = XEQ, parms = Tpars, positive = TRUE)
	test[[1]]

	remove(LIT, MIC, SOM)

# ----------------------------------------------------------
# (C)     Run to steady state using daily TSOI
# ----------------------------------------------------------

nday   <- 365 * 30			#SPEEDS up exploration
day    <- seq(1,nday,1)
year   <- day/365
doy    <- 1

#initialize arrays to store daily output data
LIT    <- array(NA, dim = c(2,nday), dimnames = c("LITpool", "cumDAY"))
MIC    <- array(NA, dim = c(2,nday), dimnames = c("MICpool", "cumDAY"))
SOM    <- array(NA, dim = c(3,nday), dimnames = c("SOMpool", "cumDAY"))


  #initialize pools
  LIT_1    <- test[[1]][[1]]    
  LIT_2    <- test[[1]][[2]]
  MIC_1    <- test[[1]][[3]]
  MIC_2    <- test[[1]][[4]]
  SOM_1    <- test[[1]][[5]]
  SOM_2    <- test[[1]][[6]]
  SOM_3    <- test[[1]][[7]]

  for (d in 1:nday)  {
    for (h in 1:24)   {
    #Fluxes at each time step
	LITmin  <- rep(NA, dim=4)
	MICtrn  <- rep(NA, dim=6)
	SOMmin  <- rep(NA, dim=2)
	DEsorb  <- rep(NA, dim=1)
	OXIDAT  <- rep(NA, dim=1)
	
    #Flows to and from MIC_1
    LITmin[1] = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + LIT_1)   #MIC_1 decomp of MET lit
    LITmin[2] = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + LIT_2)   #MIC_1 decomp of STRUC lit
    SOMmin[1] = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + SOM_3)   #Decomp of SOMa by MIC_1

    MICtrn[1] = MIC_1 * tao[1]  * fPHYS[1]                  #MIC_1 turnover to SOMp
    MICtrn[2] = MIC_1 * tao[1]  * fCHEM[1]                  #MIC_1 turnover to SOMc  
    MICtrn[3] = MIC_1 * tao[1]  * fAVAI[1]                  #MIC_1 turnover to SOMa 
    
    #Flows to and from MIC_2
    LITmin[3] = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + LIT_1)   #decomp of MET litter
    LITmin[4] = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + LIT_2)   #decomp of SRUCTURAL litter
    SOMmin[2] = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + SOM_3)   #decomp of PHYSICAL SOM by MIC_1

    MICtrn[4] = MIC_2 * tao[2]  * fPHYS[2]                  #MIC_2 turnover to SOMp 
    MICtrn[5] = MIC_2 * tao[2]  * fCHEM[2]                  #MIC_2 turnover to SOMc 
    MICtrn[6] = MIC_2 * tao[2]  * fAVAI[2]                  #MIC_2 turnover to SOMa  


	DEsorb    = SOM_1 * desorb  #* (MIC_1 + MIC_2)  #desorbtion of PHYS to AVAIL (function of fCLAY)
	OXIDAT    = ((MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + SOM_2)) +
		         (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + SOM_2)))  #oxidation of C to A


    LIT_1 = LIT_1 + I[1]*(1-FI[1]) - LITmin[1] - LITmin[3]
    MIC_1 = MIC_1 + CUE[1]*(LITmin[1]+ SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3])
    SOM_1 = SOM_1 + I[1]*FI[1] + MICtrn[1] + MICtrn[4]- DEsorb 
	    
    LIT_2 = LIT_2 + I[2] * (1-FI[2]) - LITmin[2] - LITmin[4]
    MIC_2 = MIC_2 + CUE[3]*(LITmin[3]+ SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6])  
    SOM_2 = SOM_2 + I[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT

	SOM_3 = SOM_3 + MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]	

	#write out daily results
	if (h == 24) {
	  LIT[1,d] <- LIT_1
	  LIT[2,d] <- LIT_2
	  MIC[1,d] <- MIC_1
	  MIC[2,d] <- MIC_2
	  SOM[1,d] <- SOM_1
	  SOM[2,d] <- SOM_2
	  SOM[3,d] <- SOM_3
	  
	#advancy day of year counter
	  if (doy == 365) {
	    doy <- 1 
	    print(paste(c(strSite[s], "finished initial year", year[d])))
	  } else {
	  	  doy <- doy + 1
	  }                         #close day of year counter
	}	   						#close daily results counter

#    remove(Vmax, VMAX, Km, KM)
    }							#close hour loop
  }								#close daily loop

  sLIT1 <- LIT_1
  sLIT2 <- LIT_2
  sMIC1 <- MIC_1
  sMIC2 <- MIC_2
  sSOM1 <- SOM_1
  sSOM2 <- SOM_2
  sSOM3 <- SOM_3
test[[1]]

quartz()
par(mfrow=c(3,1), mar=c(4,4,1,1))
  plot(year,  LIT[1,], lwd=3, ylim=c(min(LIT)*0.7, max(LIT))*1.15, type ="l", xlab="", main = paste(strSite[s]))
    lines(year, LIT[2,], lwd=3, col = 2)
	abline(h=test[[1]][[1]], col=1, lty=2)
	abline(h=test[[1]][[2]], col=2, lty=2)
    legend("topright", legend=c("Met","Struc"), col=c(1,2), lty = 1, 
    	lwd = 3, cex=1.3, bty="n")
    	
  plot(year,  MIC[1,], lwd=3, type ="l", xlab="", ylim=c(min(MIC)*0.7, max(MIC))*1.15, )
    lines(year, MIC[2,], lwd=3, col = 2) 
 	abline(h=test[[1]][[3]], col=1, lty=2)
 	abline(h=test[[1]][[4]], col=2, lty=2)
    legend("topright", legend=c("Mic_r","Mic_K"), col=c(1,2), lty = 1, 
    	lwd = 3, cex=1.3, bty="n")

  plot(year,  SOM[1,], lwd=3, type ="l", ylim=c(min(SOM)*0.7, max(SOM))*1.15, )
    lines(year, SOM[2,], lwd=3, col = 2)
    lines(year, SOM[3,], lwd=3, col = 4)
 	abline(h=test[[1]][[5]], col=1, lty=2)
 	abline(h=test[[1]][[6]], col=2, lty=2)
 	abline(h=test[[1]][[7]], col=4, lty=2)
    legend("topright", legend=c("Phys","Chem","Avail"), col=c(1,2,4), lty = 1, 
    	lwd = 3, cex=1.3, bty="n")


remove(nday, day, year, LIT, MIC, SOM)

# ----------------------------------------------------------
# (D)     start litter bag experiment
#            add litter Oct 1, d=144
# ----------------------------------------------------------
nday   <- 365 * 10 + 200
day    <- seq(1,nday,1)
year   <- (day-143)/365
doy    <- 1


#initialize arrays to store daily output data
LIT    <- array(NA, dim = c(2,nday), dimnames = c("LITpool", "cumDAY"))
LITBAG <- array(NA, dim = c(6,2,nday), dimnames = c("BAGpool", "cumDAY")) #for litter bag study
MIC    <- array(NA, dim = c(2,nday), dimnames = c("MICpool", "cumDAY"))
SOM    <- array(NA, dim = c(3,nday), dimnames = c("SOMpool", "cumDAY"))

for (i in 1:6)     {	#SPEEDS search
  doy <- 1
  #initialize pools
  LIT_1    <- sLIT1    
  LIT_2    <- sLIT2
  LITbag_1 <- sLIT1
  LITbag_2 <- sLIT2
  MIC_1    <- sMIC1
  MIC_2    <- sMIC2
  MICbag_1 <- sMIC1
  MICbag_2 <- sMIC2
  SOM_1    <- sSOM1
  SOM_2    <- sSOM2
  SOM_3    <- sSOM3

  LIT    <- array(NA, dim = c(2,nday), dimnames = c("LITpool", "cumDAY"))
  MIC    <- array(NA, dim = c(2,nday), dimnames = c("MICpool", "cumDAY"))
  SOM    <- array(NA, dim = c(3,nday), dimnames = c("SOMpool", "cumDAY"))

  for (d in 1:nday)  {
    for (h in 1:24)   {
   #Fluxes at each time step
	LITmin  <- rep(NA, dim=4)
	LITbag  <- rep(NA, dim=4)
	MICtrn  <- rep(NA, dim=6)
	SOMmin  <- rep(NA, dim=2)
	DEsorb  <- rep(NA, dim=1)
	OXIDAT  <- rep(NA, dim=1)
	
    #Flows to and from MIC_1
    LITmin[1] = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + LIT_1)   #MIC_1 decomp of MET lit
    LITmin[2] = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + LIT_2)   #MIC_1 decomp of STRUC lit
    LITbag[1] = MIC_1 * VMAX[1] * LITbag_1 / (KM[1] + LITbag_1)   #MIC_1 mineralization of METABOLIC litter
    LITbag[2] = MIC_1 * VMAX[2] * LITbag_2 / (KM[2] + LITbag_2)   #MIC_1 mineralization of STRUC litter
    SOMmin[1] = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + SOM_3)   #Decomp of SOMa by MIC_1

    MICtrn[1] = MIC_1 * tao[1]  * fPHYS[1]                  #MIC_1 turnover to SOMp
    MICtrn[2] = MIC_1 * tao[1]  * fCHEM[1]                  #MIC_1 turnover to SOMc  
    MICtrn[3] = MIC_1 * tao[1]  * fAVAI[1]                  #MIC_1 turnover to SOMa 
    
    #Flows to and from MIC_2
    LITmin[3] = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + LIT_1)   #decomp of MET litter
    LITmin[4] = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + LIT_2)   #decomp of SRUCTURAL litter
    LITbag[3] = MIC_2 * VMAX[4] * LITbag_1 / (KM[4] + LITbag_1)   #mineralization of MET litter
    LITbag[4] = MIC_2 * VMAX[5] * LITbag_2 / (KM[5] + LITbag_2)   #mineralization of SRUCTURAL litter
    SOMmin[2] = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + SOM_3)   #decomp of PHYSICAL SOM by MIC_1

    MICtrn[4] = MIC_2 * tao[2]  * fPHYS[2]                  #MIC_2 turnover to SOMp 
    MICtrn[5] = MIC_2 * tao[2]  * fCHEM[2]                  #MIC_2 turnover to SOMc 
    MICtrn[6] = MIC_2 * tao[2]  * fAVAI[2]                  #MIC_2 turnover to SOMa  


	DEsorb    = SOM_1 * desorb  #* (MIC_1 + MIC_2)  #desorbtion of PHYS to AVAIL (function of fCLAY)
	OXIDAT    = (MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + SOM_2)) 
		      + (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + SOM_2))  #oxidation of C to A


    LIT_1 = LIT_1 + I[1]*(1-FI[1]) - LITmin[1] - LITmin[3]
    LITbag_1 <- LITbag_1 + I[1]*(1-FI[1]) - LITbag[1] - LITbag[3]
    MIC_1 = MIC_1 + CUE[1]*(LITmin[1]+ SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3])
    SOM_1 = SOM_1 + I[1]*FI[1] + MICtrn[1] + MICtrn[4]- DEsorb 
	    
    LIT_2 = LIT_2 + I[2] * (1-FI[2]) - LITmin[2] - LITmin[4]
    LITbag_2 <- LITbag_2 + I[2] * (1-FI[2]) - LITbag[2] - LITbag[4]
    MIC_2 = MIC_2 + CUE[3]*(LITmin[3]+ SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6])  
    SOM_2 = SOM_2 + I[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT

	SOM_3 = SOM_3 + MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]	
	#add litter bag on Oct 1
	if (d == 143)   {
	  if (h == 24)  {
		LITbag_1 <- LITbag_1 + BAG[i,1]
		LITbag_2 <- LITbag_2 + BAG[i,2]
		print(paste("------added litter",LITtype[i]))
	  }
	}

	#write out daily results
	if (h == 24) {
	  LIT[1,d] <- LIT_1
	  LIT[2,d] <- LIT_2
	  LITBAG[i,1,d]  <- LITbag_1
	  LITBAG[i,2,d]  <- LITbag_2
	  MIC[1,d] <- MIC_1
	  MIC[2,d] <- MIC_2
	  SOM[1,d] <- SOM_1
	  SOM[2,d] <- SOM_2
	  SOM[3,d] <- SOM_3
	  
	#advancy day of year counter
	  if (doy == 365) {
	    doy <- 1 
#	    print(paste(c(strSite[s], "finished year", year[d])))
	  } else {
	  	  doy <- doy + 1
	  }                         #close day of year counter
	}	   						#close daily results counter

#    remove(Vmax, VMAX, Km, KM)
    }							#close hour loop
  }								#close daily loop

	
  print(paste(c("finished litter", LITtype[i])))
}

#--------------------------------------------------
#               Calculate averages
#--------------------------------------------------

allLIT  <- colSums(LIT,     dims = 1)
allBAG  <- array(NA, dim=c(6,nday))
difBAG  <- array(NA, dim=c(6,nday))
maxBAG  <- array(NA, dim=c(6,nday))
BAGleft <- array(NA, dim=c(6,nday), dimnames=list(LITtype,c(as.character(year))))

  for (i in 1:6) {
	allBAG[i,]  <- colSums(LITBAG[i,,])  #SPEED UP ANALYSIS change LITBAG[i ] to 1
	difBAG[i,]  <- allBAG[i,] - allLIT
	maxBAG      <- max(difBAG[i,])
	BAGleft[i,] <-  100* difBAG[i,] / maxBAG        #by taking LIT + BAG mineralization at each step
  }

BAGleft[1:6,1:143] <- NA
meanBAG            <- colMeans(BAGleft)
sdBAG              <- apply(BAGleft, 2, sd)

#quartz()
xx <- c(year, rev(year))
yy <- c(meanBAG + sdBAG, rev(meanBAG - sdBAG))
par(mar=c(5,5.2,0,1))

fout          <- c('results/',as.character(data$Site[s]), '.pdf')
fout          <- paste(fout, collapse="")
pdf(fout, width=5, height=5)

plot(xx,yy, type="n", main=paste(strSite[s], "leaf decomp"), ylab="Mass remaining (%)", 
	xlab="time (y)", ylim=c(0,100),cex.lab = 1.3, cex.axis = 1.2)
  axis(side = 4, at = seq(0,100,20), labels = FALSE)
  polygon(xx, yy, col="grey", border = NA)
  lines(year, meanBAG, lwd=2)
 dev.off()

BagLeft <- BAGleft
BagLeft[,1] <- LITtype
    fout          <- c('results/MIMICS_',as.character(data$Site[s]), '_LIDET.csv')
    fout          <- paste(fout, collapse="")
	write.table(BagLeft, file=fout, sep=",", col.names=TRUE, row.names=FALSE)
    print(paste("wrote out", fout))


  remove(allLIT, allBAG, difBAG, maxBAG, BAGleft, LITBAG, meanBAG, BagLeft, sdBAG, fout, xx, yy)
  
}  #finish site loop

