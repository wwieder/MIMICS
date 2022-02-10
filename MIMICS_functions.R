
rm(list=ls())

para_file <- paste("parameters_LIDET-MIM-REV_test_lowKm_10cm.csv", sep = "") 
parameters <- read.csv(para_file)
names(parameters)
attach(parameters)

# Range of environemtal conditions
fmet  <- seq(0.1,1,0.01)
fclay <- fmet
TSOI  <- seq(-20,35, 0.1)

Vmax  <- exp(Vint[1]  + TSOI * Vslope[1] ) * aV[1]
Vmax1  <- exp(Vint[1]*1.1  + TSOI * Vslope[1] ) * aV[1] #increased Vint
Vmax1b <- exp(Vint[1]*1.1  + TSOI * Vslope[1]*0.7) * aV[1] #increased Vint, decreased Vslowpe
Vmax1c <- exp(6  + TSOI * 0.045) * aV[1] #increased Vint, decreased Vslowpe
Vmax2 <- exp(Vint2[1] + TSOI * Vs1[1] + TSOI^2 * Vs2[1]) 
Vmax3 <- exp(Vint2[1] + TSOI * Vs1[1] + TSOI^2 * Vs2[1] * 0.25) 
# use this modifier to tune aV2 in parameter file 
Vmax2 <- Vmax2 * aV[1] * 5   
Vmax3 <- Vmax3 * aV[1] * 5   

plot(TSOI,Vmax, type='l', lwd=2)#,ylim=c(2e-4,2e-2))
lines(TSOI,Vmax1,  col=2, lwd=2)
lines(TSOI,Vmax1b, col=4, lwd=2)
lines(TSOI,Vmax1c, col=3, lwd=2)

lines(TSOI,Vmax2,  col=4, lwd=2, lty=2)

Km    <- exp(Kint[1]  + TSOI * Kslope[1] ) * aK[1]
Km1   <- exp(Kint[1]*0.9  + TSOI * Kslope[1] ) * aK[1]
Km1b  <- exp(Kint[1]*0.93  + TSOI * Kslope[1]*1.3 ) * aK[1]

Km    <- exp(Kint[2]  + TSOI * Kslope[2] ) * aK[2]
Km1   <- exp(Kint[2]*0.9  + TSOI * Kslope[2] ) * aK[2]
Km1b  <- exp(Kint[2]*0.93  + TSOI * Kslope[2]*1.3 ) * aK[2]
Km1c  <- exp(2.7  + TSOI * 0.035) * aK[2]

plot(TSOI,Km, type='l',  lwd=2,ylim=c(1,8))
lines(TSOI,Km1,  col=2, lwd=2)
lines(TSOI,Km1b, col=4, lwd=2)
lines(TSOI,Km1c, col=4, lwd=2)

# --------------------------------------------------
# try hard coding parameter here to play with them
fPHYS_1 <- 0.03* exp(1.3*fclay) 
fPHYS_2 <- 0.02* exp(0.8*fclay) 
fPHYS_3 <- 0.05* exp(0.8*fclay) 
fPHYS_4 <- 0.03* exp(0.4*fclay) 
plot(fPHYS_1, ylim=c(0.,0.15), 
     ylab='Fraction MIC turnover to SOMp',
     xlab='% clay') 
lines(fPHYS_2, col=1,lw=2)
points(fPHYS_3, col=4)
lines(fPHYS_4, col=4,lw=2)
legend(1,0.15,c('MIC_r','MIC_K'),pch=c(1,NA),lty=c(0,1),bty='n')

fCHEM_1 <- 0.1* exp(-3 *fmet ) 
fCHEM_2 <- 0.3* exp(-3 *fmet ) 

fAVAIL_1 <- 1-(fPHYS_1 + fCHEM_1)
fAVAIL_2 <- 1-(fPHYS_2 + fCHEM_2)
fAVAIL_3 <- 1-(fPHYS_3 + fCHEM_1)
fAVAIL_4 <- 1-(fPHYS_4 + fCHEM_2)

plot(fAVAIL_1, ylim=c(0.6,1)) 
points(fAVAIL_2, col=2)
lines(fAVAIL_3, col=1)
lines(fAVAIL_4, col=2)
abline(h=1)

fCHEM_1 <- 0.1* exp(-3 *fmet ) *3
fCHEM_2 <- 0.3* exp(-3 *fmet ) *3
fAVAIL_1 <- 1-(fPHYS_1 + fCHEM_1)
fAVAIL_2 <- 1-(fPHYS_2 + fCHEM_2)


lines(fAVAIL_1, col=1, lwd=2)
lines(fAVAIL_2, col=2, lwd=2)
abline(h=0)


Vslope1 <- 0.063
Vslope2 <- 0.077
Vslope3 <- 0.073
Vint    <- 5.47
Vint2   <- 5.2
av      <- 1.25E-08
Tsoil   <- seq(-20,30,1)

Vmax1   <- exp(Vslope1 * Tsoil + Vint) * av     # default
Vmax2   <- exp(Vslope2 * Tsoil + Vint2) * av    # modify slope & intercept
Vmax3   <- exp(Vslope3 * Tsoil + Vint) * av     # change slope only
plot(Tsoil, Vmax1, type='l', log='y')
lines(Tsoil, Vmax2, col=2, lwd=2)
lines(Tsoil, Vmax3, col=3, lwd=2, lty=2)

summary(lm(log10(Vmax1)~Tsoil))
summary(lm(log10(Vmax3)~Tsoil))

#desorption function, flux / h
fSOM_p <- c(1.00E-05, -2)
desorb <- fSOM_p[1] * exp(fSOM_p[2] * fclay) 
desorb <- 1/(desorb * 365 * 24)
plot(fclay,desorb)

q10Desorp = desorb[30] * 1.1^((Tsoil-15)/10)
plot(Tsoil,q10Desorp,type='l',lwd=2,col=4)
lines(Tsoil,q10Desorp, lwd=2, lty=2, col=4)
abline(h=desorb[30],lw=2)
# ----- Flexible microbial stoichiometry--------
# Parameterized as a function of litter quality
fmet = seq(0.3,0.6,0.01)
cnScale = sqrt(.45/fmet)
cn_K = 10*cnScale
cn_r = 6*cnScale
plot(fmet,cn_K, ylim=c(3,13), ylab='microbial C:N')
points(fmet,cn_r, col=2)

points(fmet,cn_K)

