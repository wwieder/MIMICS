# Will Wieder
# Oct 23, 2013
# Modified Aug, 2015; April 2016
# March 2017, added litter removal, 70% clay,50 y simulation
# Solve a system of non-linear equations for enzyme SOC solution
# uses packages deSolve / rootSolve

rm(list=ls())
dir <- "/Users/wwieder/Desktop/Working_files/soilCN/enzyme/Theory_2/Sulman_microbialModels/"
setwd(dir)

library(rootSolve)
library(boot)
#----------------analytical solutin using stode function-----------------------

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
    OXIDAT    = ((MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + MIC_2)) +
                 (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + MIC_1)))  #oxidation of C to A
    
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
# (A)       Read in parameters and site level data
#---------------------------------------------------------
para_file <- paste("parameters_LIDET-MIM-REV_test_lowKm.csv", sep = "") 
parameters <- read.csv(para_file)
names(parameters)
attach(parameters)

data <- read.csv("Test_SITES_1.csv", sep = ",") #site level forcing variables
names(data)
attach(data)

ANPP   <- ANPP              		# if needed, convert to gC/m2/y from g/m2/y
clay   <- CLAY/100  				    # if needed, convert from % clay to fraction
tsoi   <- MAT
nsites <- length(Site)
FI     <- FI * 0.1           #reduce FI by fraction of 10, as fPHYS
lig    <- LIG 
Nnew   <- N                                                  #N in litter additions
fMET1  <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig / Nnew)   #as partitioned in Daycent
MIMLIT <- rep(NA, nsites)           	                       #Vector for results
MIMMIC <- rep(NA, nsites)           
MIM_CO <- rep(NA, nsites)           
MIMSOC <- rep(NA, nsites)           

strSite <- as.character(data$Site)  #convert site names to string
pools  <- c('site','LITm', 'LITs', 'MICr', 'MICK', 'SOMp', 'SOMc', 'SOMa')

POOLS  <- c("LIT","MIC","SOC")
exper  <- c("Control","Warm2","Warm5", "1.3xLitter","2.0xLitter",
            "Prime_1.3","Warm2_1.3xLitter",'Warm2_1.3xprime',"Prime_1.3_fMET",
            "0xLitter","0xLitter_lowMicTurn")
npools <- length(POOLS)  
nexper <- length(exper)
nday   <- 365 * 60      #SPEEDS up exploration, ultimately want to be 50 years
LITpool<- c('LITm', 'LITs') 
MICpool<- c('MICr', 'MICk') 
SOMpool<- c('SOMp', 'SOMa','SOMc') 
ALLpools <- c(LITpool, MICpool, SOMpool,'CO2','INPUTS')

LIT    <- array(NA, dim = c(nsites,nexper,2,nday), dimnames = list(strSite,exper,LITpool,rep(NA, nday)))
MIC    <- array(NA, dim = c(nsites,nexper,2,nday), dimnames = list(strSite,exper,MICpool,rep(NA, nday)))
SOM    <- array(NA, dim = c(nsites,nexper,3,nday), dimnames = list(strSite,exper,SOMpool,rep(NA, nday)))
dim(LIT)

#---------------------------------------------------------
#---------------------------------------------------------
# (B)       RXEQ for site using STODE function
# (C)       Generate times series w/o inputs
#---------------------------------------------------------
#---------------------------------------------------------

#for (i in 1:nsites) {         #speeds up debugging     	
for (i in 1:2) {       	#speeds up debugging     	
  fMET       <- fMET1[i] 
  fCLAY      <- clay[i]
  TSOI       <- tsoi[i]
  EST_LIT_in <- ANPP[i] / (365*24) # gC/m2/h (from gC/m2/y)
  depth      <- parameters$Depth[1]
  h2y        <- 24*365
  MICROtoECO <- depth * 1e4 * 1e-3         #mgC/cm3 to g/m2
  EST_LIT    <- EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/h(from gC/m2/h) 
  
  #-----------------caclulate parameters---------------------------
  Vmax     <- exp(TSOI * Vslope + Vint) * aV
  Km       <- exp(Kslope * TSOI + Kint) * aK
  
  #ANPP strongly correlated with MAP
  Tau_MOD1 <- sqrt(ANPP[i]/Tau_MOD[1])          # basicaily standardize against NWT
  Tau_MOD2 <- Tau_MOD[4]                        # increased 3-fold for SS SOC pools
  Tau_MOD1[Tau_MOD1 < Tau_MOD[2]] <- Tau_MOD[2] # correction not used in LIDET resutls 
  Tau_MOD1[Tau_MOD1 > Tau_MOD[3]] <- Tau_MOD[3] 
  tau <- c(tau_r[1]*exp(tau_r[2]*fMET), 
           tau_K[1]*exp(tau_K[2]*fMET))   
  tau <- tau * Tau_MOD1 * Tau_MOD2
  
  #------NEW Parameters--------------
  fPHYS    <- c(fPHYS_r[1] * exp(fPHYS_r[2]*fCLAY), 
                fPHYS_K[1] * exp(fPHYS_K[2]*fCLAY)) 	            #fraction to SOMp
  fPHYS    <- fPHYS  * 0.1                                        #reduce fraction to fPHYS
  fCHEM    <- c(fCHEM_r[1] * exp(fCHEM_r[2]*fMET) * fCHEM_r[3], 
                fCHEM_K[1] * exp(fCHEM_K[2]*fMET) * fCHEM_K[3]) 	#fraction to SOMc
  fAVAI    <- 1- (fPHYS + fCHEM)
  desorb   <- fSOM_p[1] * exp(fSOM_p[2]*(fCLAY))                  #CHANGED FOR GLOBAL RUN!!!     
  pSCALAR  <- PHYS_scalar[1] * exp(PHYS_scalar[2]*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
  
  desorb   <- desorb * 0.1
  #------------MODIFY fluxes from SOM pools as function of clay content---------
  v_MOD    <- vMOD
  k_MOD    <- kMOD # to avoid writing over orig. parameters
  k_MOD[3] <- k_MOD[3] * pSCALAR    
  k_MOD[6] <- k_MOD[6] * pSCALAR    
  
  VMAX     <- Vmax * v_MOD 
  KM       <- Km / k_MOD
  
  I        <- array(NA, dim=2)              #Litter inputs to MET/STR
  I[1]     <- (EST_LIT / depth) * fMET      #partitioned to layers
  I[2]     <- (EST_LIT / depth) * (1-fMET)
  
  #initialize pools
  lit     <- I   # * 1e3
  mic     <- I   # * 1e2
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
              SOM_1 = som[1], SOM_2 = som[2], SOM_3 = som[3])
  test  <- stode(y = Ty, time = 1e6, fun = RXEQ, parms = Tpars, positive = TRUE)
  remove(lit, mic, som, I)
  
  # ----------------------------------------------------------
  # Generate plots of transient results
  # ----------------------------------------------------------
  #initialize arrays to store daily output data
  day    <- seq(1,nday,1)
  year   <- day/365

  for (j in 1:nexper) {
#  for (j in 10:10) {         #speeds up debugging
    dataOUT  <- array(NA, dim =c(nday,length(ALLpools)), dimnames = list(seq(1,nday,1),ALLpools))
    doy    <- 1
    LIT_1    <- test[[1]][[1]]               #initialize pools    
    LIT_2    <- test[[1]][[2]]
    MIC_1    <- test[[1]][[3]]
    MIC_2    <- test[[1]][[4]]
    SOM_1    <- test[[1]][[5]]
    SOM_2    <- test[[1]][[6]]
    SOM_3    <- test[[1]][[7]]
    TSOI     <- tsoi[i]
    for (d in 1:nday)  {
      I        <- array(NA, dim=2)              #Litter inputs to MET/STR
      I[1]     <- (EST_LIT / depth) * fMET      
      I[2]     <- (EST_LIT / depth) * (1-fMET)
      fMET     <- fMET1[i]
      # define each manipulation (j) here
      if (d > 3650) {                          #speed up simulation 
        if (j == 2) { TSOI = tsoi[i] + 2    }
        if (j == 3) { TSOI = tsoi[i] + 5    }
        if (j == 4) { I    = I    * 1.3     }
        if (j == 5) { I    = I    * 2.      }
        if (j == 6) { I[1] = I[1] * 1.3     }      
        if (j == 7) { TSOI = tsoi[i] + 2 
                      I = I * 1.3           }
        if (j == 8) { TSOI = tsoi[i] + 2 
                      I[1] = I[1] * 1.3     }
        if (j == 9) { I[1] = I[1] * 1.3     
                      fMET = I[1] / sum(I)  }        
        if (j == 10){ I    = I    * 0.      } 
        if (j == 11){ I    = I    * 0.       
          Tau_MOD1 <- sqrt(ANPP[i]*0/Tau_MOD[1])          # modify tau accordingly
          Tau_MOD2 <- Tau_MOD[4]                       
          Tau_MOD1[Tau_MOD1 < Tau_MOD[2]] <- Tau_MOD[2] 
          Tau_MOD1[Tau_MOD1 > Tau_MOD[3]] <- Tau_MOD[3] 
          tau <- c(tau_r[1]*exp(tau_r[2]*fMET), 
                   tau_K[1]*exp(tau_K[2]*fMET))   
          tau <- tau * Tau_MOD1 * Tau_MOD2
        }
      } # close manipulation loop
      # re-calculate parameters
      Vmax  <- exp(Vslope * TSOI + Vint) * aV
      Km    <- exp(Kslope * TSOI + Kint) * aK
      VMAX  <- Vmax * v_MOD 
      KM    <- Km / k_MOD

      tau   <- c(tau_r[1]*exp(tau_r[2]*fMET), 
               tau_K[1]*exp(tau_K[2]*fMET))   
      tau   <- tau * Tau_MOD1 * Tau_MOD2
      
      fCHEM <- c(fCHEM_r[1] * exp(fCHEM_r[2]*fMET) * fCHEM_r[3], 
                 fCHEM_K[1] * exp(fCHEM_K[2]*fMET) * fCHEM_K[3]) 	#fraction to SOMc
      fAVAI <- 1- (fPHYS + fCHEM)
      
      for (h in 1:24) {
        UPpars  <- c( I = I, VMAX = VMAX, KM = KM, CUE = CUE, 
                      fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
                      tau   = tau, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
                      desorb= desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO)
        UPy     <- c( LIT_1 = LIT_1, LIT_2 = LIT_2, 
                      MIC_1 = MIC_1, MIC_2 = MIC_2, 
                      SOM_1 = SOM_1, SOM_2 = SOM_2, SOM_3 = SOM_3 )

        # seperately calculate CO2 fluxes
        LITmin1 = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + MIC_1)   #MIC_1 decomp of MET lit
        LITmin2 = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + MIC_1)   #MIC_1 decomp of STRUC lit
        SOMmin1 = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + MIC_1)   #decomp of SOMa by MIC_1

        LITmin3 = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + MIC_2)   #decomp of MET litter
        LITmin4 = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + MIC_2)   #decomp of SRUCTURAL litter
        SOMmin2 = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + MIC_2)   #decomp of SOMa by MIC_2
            
        CO2     = (1-CUE[1])*(LITmin1+ SOMmin1) + (1-CUE[2])*(LITmin2)  +
                  (1-CUE[3])*(LITmin3+ SOMmin2) + (1-CUE[4])*(LITmin4)  #mgC/cm3/h
                    
        update <- RXEQ(y = UPy, pars = UPpars)
        LIT_1  <- LIT_1 + update[[1]][1]
        LIT_2  <- LIT_2 + update[[1]][2]
        MIC_1  <- MIC_1 + update[[1]][3]
        MIC_2  <- MIC_2 + update[[1]][4]
        SOM_1  <- SOM_1 + update[[1]][5]
        SOM_2  <- SOM_2 + update[[1]][6]
        SOM_3  <- SOM_3 + update[[1]][7]
        if (h == 24) {             #write out daily results
          LIT[i,j,1,d] <- LIT_1
          LIT[i,j,2,d] <- LIT_2
          MIC[i,j,1,d] <- MIC_1
          MIC[i,j,2,d] <- MIC_2
          SOM[i,j,1,d] <- SOM_1
          SOM[i,j,2,d] <- SOM_2
          SOM[i,j,3,d] <- SOM_3
          dataOUT[d,1] <- LIT_1
          dataOUT[d,2] <- LIT_2
          dataOUT[d,3] <- MIC_1
          dataOUT[d,4] <- MIC_2
          dataOUT[d,5] <- SOM_1
          dataOUT[d,6] <- SOM_2
          dataOUT[d,7] <- SOM_3
          dataOUT[d,8] <- CO2     # mgC/cm3/h
          dataOUT[d,9] <- sum(I)
          
          if (doy == 365) {       #advancy day of year counter
            doy <- 1 
            print(paste(strSite[i], " ",exper[j], " finished year ", year[d],sep=""))
          } else {
            doy <- doy + 1
          }                         #close day of year counter
          remove(UPpars, UPy, update)
        }	   						            #close daily results counter
      }							                #close hour loop
    }					 			                #close daily loop    

    # write out results for each experiment
    dout <- paste('Time_series/data/MIMICS_',exper[j],"_CLAY",CLAY[i],"_LIG",lig[i],'.csv', sep='')
    write.csv(dataOUT, file=dout)
  
    # plot change in stocks over time
    fout <- paste("Time_series/MOD_fPHYS",exper[j],"_",strSite[i],"_Reverse.pdf", sep="")
    pdf(fout)
    par(mfrow=c(3,1), mar=c(4,4,1,1))
    plot(year,  LIT[i,j,1,], lwd=3, 
         main = paste(strSite[i]," ",exper[j]),
         ylim = c(min(LIT, na.rm=T)*0.7, max(LIT,na.rm=T))*1.15, 
         type ="l", xlab="" )
    lines(year, LIT[i,j,2,], lwd=3, col = 2)
    abline(h=test[[1]][[1]], col=1, lty=2)
    abline(h=test[[1]][[2]], col=2, lty=2)
    legend("topright", legend=c("Met","Struc"), col=c(1,2), lty = 1, 
           lwd = 3, cex=1.3, bty="n")
  
    plot(year,  MIC[i,j,1,], lwd=3, type ="l", xlab="", 
         ylim=c(min(MIC, na.rm=T)*0.7, max(MIC, na.rm=T))*1.15)
    lines(year, MIC[i,j,2,], lwd=3, col = 2) 
    abline(h=test[[1]][[3]], col=1, lty=2)
    abline(h=test[[1]][[4]], col=2, lty=2)
    legend("topright", legend=c("Mic_r","Mic_K"), col=c(1,2), lty = 1, 
           lwd = 3, cex=1.3, bty="n")
  
    plot(year,  SOM[i,j,1,], lwd=3, type ="l", 
         ylim=c(min(SOM, na.rm=T)*0.7, max(SOM, na.rm=T))*1.15)
    lines(year, SOM[i,j,2,], lwd=3, col = 2)
    lines(year, SOM[i,j,3,], lwd=3, col = 4)
    abline(h=test[[1]][[5]], col=1, lty=2)
    abline(h=test[[1]][[6]], col=2, lty=2)
    abline(h=test[[1]][[7]], col=4, lty=2)
    legend("topright", legend=c("Phys","Chem","Avail"), col=c(1,2,4), lty = 1, 
           lwd = 3, cex=1.3, bty="n")
    dev.next()
    
    # -- changes w/ perturbations --
    if (j > 1) {
      LITdiff_m <- 100 * (LIT[i,j,,] / LIT[i,1,,] - 1)
      MICdiff_m <- 100 * (MIC[i,j,,] / MIC[i,1,,] - 1)
      SOMdiff_m <- 100 * (SOM[i,j,,] / SOM[i,1,,] - 1)

      plot( year, LITdiff_m[1,], col=1, type ="l", lwd = 3,
            ylim = c(min(LITdiff_m, na.rm=T)*0.85, max(LITdiff_m, na.rm=T))*1.15,  
            xlab="", 
            main = paste(strSite[i]," ",exper[j]))
      lines(year, LITdiff_m[2,], col=2,lwd = 3)
      abline(h=0, col=1, lty=2)
      legend("topleft", legend=c("Met","Struc"), col=c(1,2), lty = 1, 
             lwd = 3, cex=1.3, bty="n")
  
      plot( year, MICdiff_m[1,], col=1, type ="l", lwd = 3,
            xlab="", 
            ylim=c(min(MICdiff_m, na.rm=T)*0.85, max(MICdiff_m, na.rm=T))*1.15)
      lines(year, MICdiff_m[2,], col=2, lwd = 3)
      abline(h=0, col=1, lty=2)
      legend("topleft", legend=c("Mic_r","Mic_K"), col=c(1,2), lty = 1, 
             lwd = 3, cex=1.3, bty="n")
    
      plot( year, SOMdiff_m[1,], col=1, type ="l", lwd = 3,
            ylim=c(min(SOMdiff_m, na.rm=T)*0.85, max(SOMdiff_m, na.rm=T))*1.15)
      lines(year, SOMdiff_m[2,], col=2, lwd = 3)
      lines(year, SOMdiff_m[3,], col=4, lwd = 3)
      abline(h=0, col=1, lty=2)
      legend("topleft", legend=c("Phys","Chem","Avail"), col=c(1,2,4), lty = 1, 
              lwd = 3, cex=1.3, bty="n")
      dev.next()
      
    }             #close difference plots

    plot(year,  dataOUT[,8], lwd=3, type ="l",
         main = paste(strSite[i]," ",exper[j]))
    lines(year, dataOUT[,9], lwd=3, col=3)
    legend("topright", legend=c("CO2","Inputs"), col=c(1,3), lty = 1, 
           lwd = 3, cex=1.3, bty="n")
    
    dev.off()
    print(paste('wrote ',fout))
  }                 # close j loop (for each experiment)
  
  remove(test, Ty, Tpars, fout,dout, fMET, I)
  remove(LITmin, MICtrn, SOMmin)
  
}      # close i loop (sites)


#------------ FINISHED MAIN LOOP------------------------



LITall <- (LIT[,,1,] + LIT[,,2,]) * MICROtoECO / 1e3   #kg/m2
MICall <- (MIC[,,2,] + MIC[,,2,]) * MICROtoECO / 1e3
SOMall <- (SOM[,,1,] + SOM[,,2,] + SOM[,,3,]) * MICROtoECO / 1e3

# change for sandy soils, low fmet
dLIT <- array(NA, dim = c(nsites,nexper,nday), dimnames = list(strSite,exper,rep(NA, nday)))
dMIC <- array(NA, dim = c(nsites,nexper,nday), dimnames = list(strSite,exper,rep(NA, nday)))
dSOM <- array(NA, dim = c(nsites,nexper,nday), dimnames = list(strSite,exper,rep(NA, nday)))

for (j in 1:nexper) {
  dLIT[,j,] <- 100 * (LITall[,j,] / LITall[,1,] - 1)
  dMIC[,j,] <- 100 * (MICall[,j,] / MICall[,1,] - 1)
  dSOM[,j,] <- 100 * (SOMall[,j,] / SOMall[,1,] - 1)
}

dim(dLIT)
Site

fout <-  ('Time_series/MOD_fPHYS_MIMICS_summary.pdf')
pdf(fout)
par(mfrow=c(2,2), mar=c(1,5,2,0), cex = 1.3)
plot( year, SOMall[1,1,], col = 1, lwd = 3, type='l',
      main= "Clayey soils",  xaxt="n",
      ylim=c(min(SOMall[c(1,3),,], na.rm=T),max(SOMall[c(1,3),,], na.rm=T)),
      ylab=expression(paste('SOM (kg C ',m^-2,')')))
for (i in 1:8) {
  lines(year, SOMall[1,i,], col=i, lwd = 3, lty=1 )
  lines(year, SOMall[3,i,], col=i, lwd = 2, lty=2 )
}

par(mar=c(1,3,2,2))
plot( year, SOMall[2,1,], col = 1, lwd = 3, type='l',
      main = 'Sandy soils',  xaxt="n",
      ylim=c(min(SOMall[c(2,4),,], na.rm=T),max(SOMall[c(2,4),,], na.rm=T)))
for (i in 1:8) {
  lines(year, SOMall[2,i,], col=i, lwd = 3, lty=1 )
  lines(year, SOMall[4,i,], col=i, lwd = 2, lty=2 )
}

# -------------relative chagnes----------------------
par(mar=c(2,5,1,0))
plot( year, dSOM[1,1,], col = 1, lwd = 3, type='l',
      ylim=c(-15,15),
      ylab=expression(paste(Delta,' SOM (%)')))
for (i in 1:8) {
  lines(year, dSOM[1,i,], col=i, lwd = 3, lty=1 )
  lines(year, dSOM[3,i,], col=i, lwd = 2, lty=2 )
}

par(mar=c(2,3,1,2))
plot( year, dSOM[2,1,], col = 1, lwd = 3, type='l',
      ylim=c(-15,15))
for (i in 1:8) {
  lines(year, dSOM[2,i,], col=i, lwd = 3, lty=1 )
  lines(year, dSOM[4,i,], col=i, lwd = 2, lty=2 )
}

dev.off()

