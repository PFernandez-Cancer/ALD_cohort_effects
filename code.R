rm(list=ls())
library(OpenMx)

# Code based on:
# Hunter, M. D. (2018). State space modeling in an open source, modular, structural equation modeling environment. Structural Equation Modeling: A Multidisciplinary Journal, 25(2), 307-324.


##################################
### SSM-CT WITH COHORT EFFECTS ###
##################################

dataW <- readRDS("example_SSM.rds")

## Generate data in a format adequate for SSM ----
n <- nrow(dataW)
dataL <- list()

for(i in 1:n){
  dataL[[i]] <- data.frame(coh = dataW[i,"coh"], 
                           age = dataW[i,grep("age", colnames(dataW), value=TRUE)], 
                           y = dataW[i, grep("Y", colnames(dataW), value=TRUE)])
  colnames(dataL[[i]]) <- c("coh", "age", "y")
  rownames(dataL[[i]]) <- NULL }



## Specify State-Space Model ----
opmxL <- list()

# Autoregressive dynamics
opmxL$Acoef <- mxMatrix(name="Acoef", nrow=2, ncol=3, free = c(T,F,T,F,F,F),
                        values = c(-.45, 0, -.1818, 0, 1, 0),
                        labels = c("beta", NA, "reg_beta", NA,NA,NA) )
opmxL$Acoh <- mxMatrix(name="Acoh", nrow=3, ncol=2, free = FALSE,
                       values = c(1, NA, 0, 0, 0, 1),
                       labels = c(NA, "data.coh", NA,NA,NA,NA))
opmxL$amat_ct <- with(opmxL, mxAlgebra(name = "A", expression = Acoef%*%Acoh, 
                                       dimnames = list(c("y0", "yA"), c("y0", "yA") )))

# Input effects on the latent variables
opmxL$bmat <- mxMatrix(name = "B", "Zero", 2, 1)

#  Factor loadings in the measurement model
opmxL$cmat <- mxMatrix(name = "C", "Full", 1, 2, free = FALSE, values = c(1,0), 
                       dimnames = list(c("y"), c("y0", "yA")) )

# Input effects on the observed variables
opmxL$dmat <- mxMatrix("Zero", 1, 1, name = "D")

# Dynamic error (i.e., innovations)
opmxL$qmat <- mxMatrix("Zero", 2, 2, name = "Q")

# Measurement error
opmxL$rmat <- mxMatrix("Diag", 1, 1, TRUE, 2, name = "R", labels = "mery")

## Mean vector
opmxL$inicoh <- mxMatrix(name='inicoh', nrow=2, ncol=1, values=c(1, NA), 
                         labels=c(NA, 'data.coh'))
opmxL$inicoef <- mxMatrix(name='iniCoef', nrow=2, ncol=2, 
                          values=c(10, 13.5, 0, -1.0909),
                          # free = TRUE,
                          free = c(T,T,F,T),
                          labels=c('y0mn', 'yAmn', 'reg_y0mn', 'reg_yAmn') )
opmxL$xmat <- with(opmxL, mxAlgebra(name='x0', iniCoef %*% inicoh) )


### Covariance matrix
# Standard deviation of y0
opmxL$sd_p011 <- mxMatrix(name = "sd_p011", type = "Full", nrow=1, ncol=1, 
                          free=TRUE, values=5, labels="y0sd")

# Standard deviation of yA
opmxL$sd_coef <- mxMatrix(name='sd_coef', nrow=1, ncol=2,
                          values=c(2.25, -.1),
                          free = TRUE,
                          labels=c('yAsd', 'reg_yAsd') )
opmxL$sd_p022 <- with(opmxL, mxAlgebra(name="sd_p022", sd_coef %*% inicoh))

# Diagonal of standard deviations
opmxL$diag_sd <- with(opmxL, mxAlgebra(name="diag_sd", 
                                       expression = rbind(cbind(sd_p011, 0),
                                                          cbind(0, sd_p022))))

# Correlation matrix
opmxL$corrs <- mxMatrix(name="corrs", nrow=2, ncol=2, type="Symm",
                        values=c(1, .7, 1),
                        free = c(F, T, F),
                        labels = c(NA, "y0Ar", NA))

# Covariance matrix P0
opmxL$pmat <- with(opmxL, mxAlgebra(name="P0", expression = diag_sd%*%corrs%*%diag_sd))

# covariates
opmxL$umat <- mxMatrix("Zero", 1, 1, name = "u")

# Specification of the time index
opmxL$tmat <- mxMatrix('Full', 1, 1, name='time', labels='data.age')

# Store matrices and algebra
opmxL$modataL_ct <- with(opmxL, list(opmxL$amat_ct, opmxL$Acoef, opmxL$Acoh,
                                     opmxL$bmat, opmxL$cmat, opmxL$dmat,
                                     opmxL$qmat, opmxL$rmat, opmxL$umat, opmxL$tmat,
                                     opmxL$xmat, opmxL$inicoef, opmxL$inicoh, 
                                     opmxL$sd_p011, opmxL$sd_coef, opmxL$sd_p022, 
                                     opmxL$diag_sd, opmxL$corrs, opmxL$pmat))
# Select expectation
opmxL$expSSCT <- mxExpectationStateSpaceContinuousTime(A = "A", B = "B",
                                                       C = "C", D = "D",
                                                       Q = "Q", R = "R",
                                                       x0 = "x0", P0 = "P0",
                                                       u = "u", t = "time")

## Create SSM-CT multisubject model
indivmodels <- list()
modNames <- paste0("indiv", 1:n)
  
for(k in 1:n){
  DataSetForSubjectK <- dataL[[k]]
  indivmodels[[k]] <- mxModel(name=modNames[k],
                              opmxL$modataL_ct,
                              opmxL$expSSCT,
                              mxFitFunctionML(),
                              mxData(DataSetForSubjectK, type='raw')) 
}
  
LCS.model <- mxModel(name="LCS_CT", indivmodels,
                     mxFitFunctionMultigroup(modNames))
  
## Run model
LCS.fit <- mxRun(LCS.model)
summary(LCS.fit)





##################################
### SEM-DT WITH COHORT EFFECTS ###
##################################

dataW <- readRDS("example_SSM.rds")

# Transform data for SEM-DT
dataW <- data.frame(readRDS("example_SSM.rds")[,c(1, 6:9)])

d.0 <- dataW[dataW$coh == 0,2:5]
d.1 <- dataW[dataW$coh == 1,2:5]
d.2 <- dataW[dataW$coh == 2,2:5]
d.3 <- dataW[dataW$coh == 3,2:5]
d.4 <- dataW[dataW$coh == 4,2:5]
d.5 <- dataW[dataW$coh == 5,2:5]
d.6 <- dataW[dataW$coh == 6,2:5]
d.7 <- dataW[dataW$coh == 7,2:5]
d.8 <- dataW[dataW$coh == 8,2:5]
d.9 <- dataW[dataW$coh == 9,2:5]
d.10 <- dataW[dataW$coh == 10,2:5]
d.11 <- dataW[dataW$coh == 11,2:5]

d.0NA <- cbind(d.0, 
               matrix(NA, ncol=11, nrow=nrow(d.0)))

d.1NA <- cbind(matrix(NA, ncol=1, nrow=nrow(d.1)),  d.1, 
               matrix(NA, ncol=10, nrow=nrow(d.1)))

d.2NA <- cbind(matrix(NA, ncol=2, nrow=nrow(d.2)),  d.2, 
               matrix(NA, ncol=9, nrow=nrow(d.2)))

d.3NA <- cbind(matrix(NA, ncol=3, nrow=nrow(d.3)),  d.3, 
               matrix(NA, ncol=8, nrow=nrow(d.3)))

d.4NA <- cbind(matrix(NA, ncol=4, nrow=nrow(d.4)),  d.4, 
               matrix(NA, ncol=7, nrow=nrow(d.4)))

d.5NA <- cbind(matrix(NA, ncol=5, nrow=nrow(d.5)),  d.5, 
               matrix(NA, ncol=6, nrow=nrow(d.5)))

d.6NA <- cbind(matrix(NA, ncol=6, nrow=nrow(d.6)),  d.6, 
               matrix(NA, ncol=5, nrow=nrow(d.6)))

d.7NA <- cbind(matrix(NA, ncol=7, nrow=nrow(d.7)),  d.7, 
               matrix(NA, ncol=4, nrow=nrow(d.7)))

d.8NA <- cbind(matrix(NA, ncol=8, nrow=nrow(d.8)),  d.8, 
               matrix(NA, ncol=3, nrow=nrow(d.8)))

d.9NA <- cbind(matrix(NA, ncol=9, nrow=nrow(d.9)),  d.9, 
               matrix(NA, ncol=2, nrow=nrow(d.9)))

d.10NA <- cbind(matrix(NA, ncol=10, nrow=nrow(d.10)),  d.10, 
                matrix(NA, ncol=1, nrow=nrow(d.10)))

d.11NA <- cbind(matrix(NA, ncol=11, nrow=nrow(d.11)),  d.11)

colnames(d.0NA) <- paste0("Y", 0:14)
colnames(d.1NA) <- paste0("Y", 0:14)
colnames(d.2NA) <- paste0("Y", 0:14)
colnames(d.3NA) <- paste0("Y", 0:14)
colnames(d.4NA) <- paste0("Y", 0:14)
colnames(d.5NA) <- paste0("Y", 0:14)
colnames(d.6NA) <- paste0("Y", 0:14)
colnames(d.7NA) <- paste0("Y", 0:14)
colnames(d.8NA) <- paste0("Y", 0:14)
colnames(d.9NA) <- paste0("Y", 0:14)
colnames(d.10NA) <- paste0("Y", 0:14)
colnames(d.11NA) <- paste0("Y", 0:14)

data <- rbind(d.0NA, d.1NA, d.2NA, d.3NA, d.4NA, d.5NA, 
              d.6NA, d.7NA, d.8NA, d.9NA, d.10NA, d.11NA)
data$coh <- dataW$coh


## Build model
Tmax <- 15
Y_manif <- paste0("Y", 0:(Tmax-1))
y_lat <- paste0("y", 0:(Tmax-1))

LCS <- mxModel("LCSxcoh", 
               type="RAM",
               mxData(observed = data, type="raw"),
               manifestVars=Y_manif,
               latentVars=c("y00", "yA", y_lat),
               # From latent to manifest
               mxPath(from=y_lat, to=Y_manif, arrows=1, free=FALSE, values=1),
               # Measurement error variance
               mxPath(from=Y_manif, arrows=2, free=TRUE, values=2, labels="mery"),
               # From initial condition to first measurement occasion
               mxPath(from="y00", to=y_lat[1], arrows=1, free=FALSE, values=1),
               # From additive component to latent
               mxPath(from="yA", to=y_lat[2:Tmax], arrows=1, 
                      free=FALSE, values=1),
               
               # Self-feedback effect
               mxMatrix("Full", 1, 2, free=TRUE, values=c(.55, .03636), 
                        labels=c("beta", "lambda_beta"), name="b_coeff"),
               mxMatrix("Full", 2, 1, free=FALSE, values=c(1, NA), 
                        labels=c(NA, "data.coh"), name="coh_coeff"),
               mxAlgebra(b_coeff%*%coh_coeff, name="regbeta"),
               mxPath(from=y_lat[1:(Tmax-1)], to=y_lat[2:Tmax], arrows=1, values=NA,
                      free=F, labels="regbeta[1,1]"),
               
               # Latent variances
               mxMatrix("Full", 1, 2, free=TRUE, values=c(2.25, -.1818),
                        labels=c("sdyA", "lambda_yAv"), name="v_trans1"),
               mxMatrix("Full", 2, 2, free=FALSE, values=c(1, 0, 0, NA), 
                        labels=c(NA, NA, NA, "data.coh"), name="v_trans2"),
               mxMatrix("Full", 2, 2, free=FALSE, values=c(1,1,1,1), name="unity"),
               mxAlgebra(v_trans1%*%v_trans2, name="v_coeff"),
               mxAlgebra(v_coeff%*%unity%*%t(v_coeff), name="regvyA"),
               
               mxPath(from="y00", arrows=2, free=T, labels="y00v", values=25),
               mxPath(from="yA", arrows=2, free=F, labels="regvyA[1,1]"),
               
               # Latent covariance
               mxMatrix("Full", 1, 2, free=TRUE, values=c(7.875, -.6363), 
                        labels=c("y0Acv", "lambda_y0Acv"), name="cov_coeff"),
               mxAlgebra(cov_coeff%*%coh_coeff, name="regCovar"),
               mxPath(from = "yA", to = "y00", arrows = 2, free=FALSE,
                      labels = "regCovar[1,1]"),
               
               # Mean structure
               mxMatrix("Full", nrow=2, ncol=2, 
                        # free=TRUE,
                        free=c(T,T,F,T),
                        values=c(10, 13.5, 0, -1.0909),
                        labels=c('y0mn', 'yAmn',
                                 'lambda_y0mn', 'lambda_yAmn'), name='iniCoef' ),
               mxAlgebra(name='initmeans', iniCoef %*% coh_coeff),
               mxPath(from="one", to=c("y00", "yA"), arrows=1, free=F, 
                      labels=c("initmeans[1,1]", "initmeans[2,1]")) 
)

LCS.fit <- mxRun(LCS)

s <- summary(LCS.fit)
print(s)













