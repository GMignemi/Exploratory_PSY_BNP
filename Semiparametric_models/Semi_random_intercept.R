setwd("C:/Users/giuse/Desktop/Folders/BNP Explorative Psy")
library(haven)
library(nimble)
################################################################################
# ----------------------------- READ DATA -------------------------------------- 
data <- read_sav("PTSD_data_V1.sav")

naniar::vis_miss(Data)

# Carica il pacchetto dplyr
library(dplyr)

# Filtra le righe con casi completi nelle tre colonne specifiche
Data <- data_full %>%
  filter(!is.na(SUM_PM_PCL5) &
           !is.na(PCS) &
           !is.na(MCS)&
           !is.na(Sum_GAD2)&
           !is.na(Sum_PHQ2)
  )

# ------------------------------  MODEL  --------------------------------------- 
code <- nimbleCode({  
  for(i in 1:I) {
    y[i] ~ dnorm(beta_0 + beta_1*x1[i] + beta_2*x2[i] + beta_3*x3[i]+ beta_4*x4[i], sigma  ) # devo dirgli che distrib ha nei miei dati la variabile outcome 
  }
  
  beta_1 ~ dnorm(0, var = 30)
  beta_2 ~ dnorm(0, var = 30)
  beta_3 ~ dnorm(0, var = 30)
  beta_4 ~ dnorm(0, var = 30)
  
  ## Individual effects
  
  ## CRP for clustering individual effects
  zi[1:I] ~ dCRP(alpha, size = I)
  alpha   ~ dgamma(a, b)  
  ## Mixture component parameter drawn from the base measure
  for(i in 1:I) {
    beta_0[i] ~ dnorm(mu[i], var = s2[i])  
    mu[i]     <- muTilde[zi[i]]                 
    s2[i]     <- s2Tilde[zi[i]]   
  }
  
  for(m in 1:M) {
    muTilde[m] ~ dnorm(0, var = s2_mu)
    s2Tilde[m] ~ dinvgamma(nu1, nu2)
  }
  
  sigma ~ dgamma(0.5,0.5)  # prior for variance components based on Gelman (2006)
})


# --------------- VARIABLE PREPARATION -----------------------------------------
SUM_PM_PCL5 = Data$SUM_PM_PCL5
PCS         = Data$PCS
MCS         = Data$MCS
Sum_GAD2    = Data$Sum_GAD2
Sum_PHQ2    = Data$Sum_PHQ2

# ------------------------- CONSTANT -------------------------------------------

constants <- list(  I = length(Data),
                    x1 = PCS,
                    x2 = MCS,
                    x3 = Sum_GAD2,
                    x4 = Sum_PHQ2,
                    M=50)

# ------------------- OUTCOME VARIABLE -----------------------------------------

data = list(y = SUM_PM_PCL5)

# ------------------- INITIAL VALUES-- -----------------------------------------

inits <- list(beta_1 = 0, beta_2 = 0,beta_3 = 0,beta_4 = 0, sigma = 1,
              zi     = sample(1:constants$M,constants$I,replace = TRUE),
              alpha  = 1,
              muTilde= rep(0,constants$M),
              s2Tilde= rep(1,constants$M),
              mu     = rep(0,constants$I),
              s2     = rep(1,constants$I),
              nu1    = 2.01, 
              nu2    = 1.01,
              s2_mu  = 2,
              a      = 1,
              b      = 3)

# ----------------------- PUT TOGETHER -----------------------------------------

model <- nimbleModel(code, constants = constants, data = data, inits = inits)
monitors = c("beta_1","beta_2", "beta_3","beta_4", "beta_0","eta","zi", "muTilde", "s2Tilde", "alpha" ) # parameters

model               <- nimbleModel(code, constants, data, inits)

cmodel              <- compileNimble(model)

conf                <- configureMCMC(model, monitors = monitors)

modelMCMC           <- buildMCMC(conf)
cModelMCMC          <- compileNimble(modelMCMC, project = model)

# --------------------------- RUN MCMC -----------------------------------------

system.time(samples <- runMCMC(cModelMCMC, niter=55000, nburnin = 5000, thin=10 ))

################################################################################
################################################################################
################################################################################

# ----------------------- CALL THE VARIABLES -----------------------------------

# TO CHECK
betaCols1 <- grep("beta_1", colnames(samples))
betaCols2 <- grep("beta_2", colnames(samples))
betaCols3 <- grep("beta_3", colnames(samples))
betaCols4 <- grep("beta_4", colnames(samples))

betaCols0 <- grep("beta_1", colnames(samples))

# -------------------------- SUMARY OUTPUT -------------------------------------

samplesSummary(samples[, c(betaCols1)])
samplesSummary(samples[, c(betaCols2)])
samplesSummary(samples[, c(betaCols3)])
samplesSummary(samples[, c(betaCols4)])
samplesSummary(samples[, c(betaCols0)])


#------------------------ TRACE PLOT -------------------------------------------
par(mfrow = c(2, 2), cex = 1.1)

ts.plot(samples[ , betaCols1], xlab = 'iteration', ylab = colnames(samples)[ betaCols1])
ts.plot(samples[ , betaCols2], xlab = 'iteration', ylab = colnames(samples)[ betaCols2])
ts.plot(samples[ , betaCols3], xlab = 'iteration', ylab = colnames(samples)[ betaCols3])
ts.plot(samples[ , betaCols4], xlab = 'iteration', ylab = colnames(samples)[ betaCols4])

#----------------------- from DP prior -----------------------------------------
# https://r-nimble.org/bayesian-nonparametric-models-in-nimble-general-multivariate-models
# https://www.r-bloggers.com/2018/12/bayesian-nonparametric-models-in-nimble-part-2-nonparametric-random-effects/
# https://hal.science/hal-03367099/document

samplesG <- getSamplesDPmeasure(cModelMCMC)









################################################################################
#------------------------ MODEL FIT WAIC ---------------------------------------

calculateWAIC(samples, model)

















