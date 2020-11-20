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
Km    <- exp(Kint[1]  + TSOI * Kslope[1] ) * aK[1]
Vmax2 <- exp(Vint2[1] + TSOI * Vs1[1] + TSOI^2 * Vs2[1]) 
Vmax3 <- exp(Vint2[1] + TSOI * Vs1[1] + TSOI^2 * Vs2[1] * 0.25) 
# use this modifier to tune aV2 in parameter file 
Vmax2 <- Vmax2 * aV[1] * 5   
Vmax3 <- Vmax3 * aV[1] * 5   

plot(TSOI,Vmax, type='l', log='y', lwd=2)
lines(TSOI,Vmax2, col=2, lwd=2)
lines(TSOI,Vmax3, col=4, lwd=2, lty=2)
# --------------------------------------------------
# try hard coding parameter here to play with them
fPHYS_1 <- 0.3* exp(1.3*fclay) 
fPHYS_2 <- 0.2* exp(0.8*fclay) 
fCHEM_1 <- 0.1* exp(-3 *fmet ) 
fCHEM_2 <- 0.3* exp(-3 *fmet ) 

fAVAIL_1 <- 1-(fPHYS_1 + fCHEM_1)
fAVAIL_2 <- 1-(fPHYS_2 + fCHEM_2)

plot(fAVAIL_1, ylim=c(-0.1,0.8)) 
points(fAVAIL_2, col=2)
abline(h=0)

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

# ----- Flexible microbial stoichiometry--------
# Parameterized as a function of litter quality
fmet = seq(0.3,0.6,0.01)
cnScale = sqrt(.45/fmet)
cn_K = 10*cnScale
cn_r = 6*cnScale
plot(fmet,cn_K, ylim=c(3,13), ylab='microbial C:N')
points(fmet,cn_r, col=2)

points(fmet,cn_K)

