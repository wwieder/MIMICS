# Will Wieder [wwieder@ucar.edu]
# Oct 23, 2013
# Modified Aug, 2015; April 2016; Dec 2018
# Solve a system of non-linear equations for enzyme SOC solution
# uses packages deSolve / rootSolve
# here uses a quadratic for Vmax based on data from Machuuller & Crowther
# Vmax = e^(Vint + Vs1 * T + Vs2 * T^2) * av 

# !! this code base has not yet been published !!
# Please check with W. Wieder before distributing this code
# Here parameter values were manually tuned to match steady states SOC pools across LTER / LIDET sites
# Checked Jan 2016, for handoff to Wally and Steve for data assimilation. 
# Solve a system of non-linear equations for enzyme SOC solution
# uses packages deSolve / rootSolve

rm(list=ls())
library(rootSolve)
#----------------analytical solutiuon using stode function-----------------------

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
    #can make fluxes from CHEM a function of microbial biomass size?
    
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
para_file <- paste("parameters_LIDET-MIM-REV_test_lowKm_10cm.csv", sep = "") 
parameters <- read.csv(para_file)
names(parameters)
attach(parameters)

data <- read.csv("LTER_SITE_1.csv") #site level forcing variables
names(data)
attach(data)

ANPP    <- ANPP / 2       			# convert to gC/m2/y from g/m2/y
clay   <- CLAY2/100  				    # convert from clay fraction to %
tsoi   <- MAT
nsites <- length(Site)

lig    <- LIG/100
Nnew   <- 1/CN/2.5                  	                       # N in litter additions
fMET1  <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig / Nnew)   # as partitioned in Daycent
MIMLIT <- array(NA, dim=c(3,nsites))           	             # array for results
MIMMIC <- array(NA, dim=c(3,nsites))           
MIM_CO <- array(NA, dim=c(3,nsites))           
MIMSOC <- array(NA, dim=c(3,nsites))           

strSite <- as.character(data$Site)  #convert site names to string
pools  <- c('site','LITm', 'LITs', 'MICr', 'MICK', 'SOMp', 'SOMc', 'SOMa')
exp    <- c('Vmax1', 'Vmax2','Vmax2*')
table  <- array(NA, dim=c(3, nsites,8), dimnames=list(exp,as.character(Site),pools)) 
table[,,1]       <- as.character(Site)
table_doub      <- table

#---------------------------------------------------------
# (B)       RXEQ for site using STODE function
#---------------------------------------------------------
for (h in 1:length(exp)) {            # for default & quadratic Vmax formualtion
  for (i in 1:nsites)      {    			# speeds up debugging     	
    fMET       <- fMET1[i] 
    fCLAY      <- clay[i]
    TSOI       <- tsoi[i]
    EST_LIT_in <- ANPP[i] / (365*24) # gC/m2/h (from gC/m2/y)
    depth      <- parameters$Depth[1] 
    h2y        <- 24*365
    MICROtoECO <- depth * 1e4 * 1e6 / 1e6    #mgC/cm3 to kgC/km2
    EST_LIT    <- EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/h(from gC/m2/h) 
    EST_LIT2   <- EST_LIT * 2	               #double ANPP   
  
    #-----------------caclulate parameters---------------------------
    Vmax     <- exp(Vint  + TSOI * Vslope ) * aV
    Km       <- exp(Kint  + TSOI * Kslope ) * aK
    Vmax2    <- exp(Vint2 + TSOI * Vs1 + TSOI^2 * Vs2) * aV2
    # try reducing biases in stocks
    Vmax3    <- exp(Vint2 + TSOI * Vs1 + TSOI^2 * Vs2 * 0.25) * aV2
    if (h == 1) {Vmax <- Vmax}            #default
    if (h == 2) {Vmax <- Vmax2}           #quadratic
    if (h == 3) {Vmax <- Vmax3}           #quadratic, less curve from Vs2
    
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
    
    I        <- array(NA, dim=2)              #Litter inputs to MET/STR
    I[1]     <- (EST_LIT / depth) * fMET      #partitioned to layers
    I[2]     <- (EST_LIT / depth) * (1-fMET)
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
    
    #Calculate RXEQ pools
    Tpars <- c( I = I, VMAX = VMAX, KM = KM, CUE = CUE, 
                fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
                tau = tau, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
                desorb = desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO)
    Ty    <- c( LIT_1 = LIT[1], LIT_2 = LIT[2], 
                MIC_1 = MIC[1], MIC_2 = MIC[2], 
                SOM_1 = SOM[1], SOM_2 = SOM[2], SOM_3 = SOM[3] )
    test  <- stode(y = Ty, time = 1e6, fun = RXEQ, parms = Tpars, positive = TRUE)
  
    table[h, i,2:8] <- as.numeric(test[[1]])
    MIMLIT[h, i]    <- (test[[1]][[1]]+test[[1]][[2]])  * depth *1e4 / 1e6#convert kgC/m2 from mgC/cm3 (0-30 cm) 
    MIMMIC[h, i]    <- (test[[1]][[3]]+test[[1]][[4]])  * depth *1e4 / 1e6
    MIM_CO[h, i]    <-  test[[1]][[3]]/test[[1]][[4]]
    MIMSOC[h, i]    <- sum(test[[1]])  * depth *1e4 / 1e6   
  
    print(paste(exp[h],Site[[i]]))
    print(test[[1]])
  
    remove(LIT, MIC, SOM, test, Ty, Tpars)
    remove(LITmin, MICtrn, SOMmin)
  }      # close i loop (sites)
}        # chose h loop (vmax)



# ----------------------------------------
# ---------- now plot results ------------
# ----------------------------------------

# results from DAYCENT run across sites, sensu Wieder et al. 2014, Global Biogeochemical Cycles
DAYLIT <- c(1321,2880,1401,2066,708,2322,3131,685,1166,2918,969,
            1431,1412,826)
DAYSOC <- c(4227.32, 7549.74, 6003.994, 5385.7, 1908.9, 6799.8, 9569.7,
            2733.78, 4014.1, 9832.8, 5708.1, 4296.9, 4501.7, 3828.4)
DAYSOC      <- DAYSOC/1000 #kgC/m2

summary(lm(SOC ~ MIMSOC[1,]))
summary(lm(SOC ~ MIMSOC[2,]))
summary(lm(SOC ~ MIMSOC[3,]))
summary(lm(SOC ~ DAYSOC))
T1a <- cor.test(SOC, MIMSOC[1,])
T1b <- cor.test(SOC, MIMSOC[2,])
T1c <- cor.test(SOC, MIMSOC[3,])
T2  <- cor.test(SOC, DAYSOC)

cor1texta <- paste(exp[1], "r = ", signif(T1a[[4]][[1]], digits=2))
cor1textb <- paste(exp[2], "r = ", signif(T1b[[4]][[1]], digits=2))
cor1textc <- paste(exp[3], "r = ", signif(T1c[[4]][[1]], digits=2))
cor2text  <- paste("DAYCENT r = ", signif(T2[[4]][[1]],  digits=2))

fout <- 'Steady_state_SOC_test_lowKm_10cm.pdf'
pdf(fout,width=6.5, height=4)
par(mar=c(5,5,2,1))
plot(MIMSOC[1,], SOC,ylim=c(0,10), xlim=c(0,10), col=0,
      pch=16, cex.lab=1.3, cex.axis = 1.2, cex=1.8,
      main = 'default parameterizaton',
      ylab=expression(paste("Observed SOC (kg C ", m^-2, ")")),
      xlab=expression(paste("Predicted SOC (kg C ", m^-2, ")")))
legend("topleft", pch = 16, cex=1.3, col=0, pt.cex=1.8, 
       legend=cor1texta,bty="n")
abline(0,1, lty=2)
text(10,9.5,"1:1")
text(MIMSOC[1,], SOC, labels = Site,col=1)
dev.next()

plot(MIMSOC[2,], SOC,ylim=c(0,10), xlim=c(0,10), col=0,
     pch=16, cex.lab=1.3, cex.axis = 1.2, cex=1.8,
     main = 'Machmuller Vmax',
     ylab=expression(paste("Observed SOC (kg C ", m^-2, ")")),
     xlab=expression(paste("Predicted SOC (kg C ", m^-2, ")")))
legend("topleft", pch = 16, cex=1.3, col=0, pt.cex=1.8, 
       legend=cor1textb,bty="n")
abline(0,1, lty=2)
text(10,9.5,"1:1")
text(MIMSOC[2,], SOC, labels = Site,col=2)
dev.next()

plot(MIMSOC[3,], SOC,ylim=c(0,10), xlim=c(0,10), col=0,
     pch=16, cex.lab=1.3, cex.axis = 1.2, cex=1.8,
     main = 'Machmuller Vmax*',
     ylab=expression(paste("Observed SOC (kg C ", m^-2, ")")),
     xlab=expression(paste("Predicted SOC (kg C ", m^-2, ")")))
legend("topleft", pch = 16, cex=1.3, col=0, pt.cex=1.8, 
       legend=cor1textc,bty="n")
abline(0,1, lty=2)
text(10,9.5,"1:1")
text(MIMSOC[3,], SOC, labels = Site,col=4)
dev.next()

bias_1 <- MIMSOC[1,]-SOC
bias_2 <- MIMSOC[2,]-SOC
bias_3 <- MIMSOC[3,]-SOC

plot(bias_1~MAT, pch=16, col=1, cex=1.8,
     ylab="MIMICS bias", ylim=range(bias_2))
points(bias_2~MAT, pch=16, col=2, cex=1.8)
points(bias_3~MAT, pch=16, col=4, cex=1.8)
abline(h=0, lty=2)
legend("topleft", pch = 16, cex=1.3, col=c(1,2,4), pt.cex=1.8, 
         legend=exp, bty="n")
dev.next()  


plot(MIMSOC[1,], SOC,ylim=c(0,11), xlim=c(0,11), col=1,
     pch=16, cex.lab=1.3, cex.axis = 1.2, cex=1.8,
     ylab=expression(paste("Observed SOC (kg C ", m^-2, ")")),
     xlab=expression(paste("Predicted SOC (kg C ", m^-2, ")")))
points(MIMSOC[2,],SOC, col=2, pch=16, cex=1.8)
points(DAYSOC,SOC, col=4,     pch=17, cex=1.8)
legend("topleft", pch = c(16,16,17), cex=1.3, col=c(1,2,4), pt.cex=1.8, 
       legend=c(cor1texta, cor1textb, cor2text), bty="n")
abline(0,1, lty=2)
text(11,10.5,"1:1")
#text(0.1, 9.75, "(a)", cex=1.4)
dev.off()


# ------------------------------------------


fout          <- 'LTER_STEADY_control_test_lowKm_10cm.csv'
write.table(table, file=fout, sep=",", col.names=TRUE, row.names=FALSE)
print(paste("wrote out", fout))

header <- c('site','OBS_SOC', 'MIM_SOC')

table_summary     <- array(NA, dim=c(nsites,3), dimnames=list(as.character(Site),header)) 
table_summary[,1] <- as.character(Site)
table_summary[,2] <- SOC
table_summary[,3] <- MIMSOC
fout          <- 'LTER_Summary_Table_10cm.csv'
write.table(table_summary, file=fout, sep=",", col.names=TRUE, row.names=FALSE)
print(paste("wrote out", fout))
