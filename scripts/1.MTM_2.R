#' Fits a (Bayesian) Multivariate Gaussian Mixed Effects Model using a Gibbs
#' Sampler.
#'
#' Data equation: Y = 1mu' + XB + U1 + ... + Uq + E \cr\cr where: \itemize{
#' \item{Y (numeric, nxp) is a matrix of phenotypes (individuals in rows, traits
#' in columns, NAs accepted),} \item{mu is a vector of (p) intercepts (included
#' by default),} \item{X (nxq) is an incidence matrix for q fixed effects,}
#' \item{B is a matrix of fixed effects,} \item{U1, ..., Uq are (nxq)
#' matrices of random effects with vec(Uj)~N(0, kronecker(Gj, Kj)), where Gj is a
#' pxp (unknown) co-variance matrix, Kj is an nxn user-defined covariance
#' matrix (kernel) which must by numeric, symmetric positive semi-definite,}
#' \item{E (nxp) is a matrix of model residuals, assumed to follow a MVN
#' distribution vec(E)~N(0, kronecker(R0, I)).} }
#'
#' Intercepts are included by default, fixed and random effects are optional.
#' The same fixed effects are applied to all traits. For each random effect the
#' user must provide a kernel (K). By default the residual co-variance matrix
#' (R0) and the co-variance matrices of random effects are un-structured;
#' however users can specify other models (e.g., DIAG=diagonal, FA=factor
#' analysis, and REC=recursive).
#'
#' @param Y Phenotype matrix (nxp numeric, traits (p) in columns, individuals
#'   (n) in rows).
#' @param K A 2-level list, 1st level defines random effects, inside each level
#'   a list is used to provide the kernel (K), a covariance structure (type,
#'   'UN', 'DIAG', 'FA' supported) and hyper-parameters (degree of freedom, df0,
#'   and scale, S0).
#' @param resCov A list used to define the co-variance matrix for model
#'   residuals (R0). Example: resCov=list(type='UN', df0=x, S0=V) specifies an
#'   un-structured covariance matrix, with an Inverse Whishart prior with degree
#'   of freedom df0 (scalar) and scale matrix (pxp) V.
#' @param nIter The number of iterations (integer).
#' @param thin Thinin interval (integer).
#' @param burnIn The number of iterations to be discarded as burn-in (integer).
#' @param XF A numeric design matrix (nxq) for fixed effects. For factors use
#'   XF=as.matrix(model.matrix(~x+y...))[, -1].
#' @param saveAt A character path and a prefix used to define where to store
#'   samples (e.g., saveAt='c:/mtmFit/test_'. By default samples are saved in the
#'   current directory and filenames have no prefix.
#' @param tolD A numeric parameter used to define the minimum eigenvalue to be
#'   maintained in the model. Eigenvectors of kernels smaller than tolD are
#'   removed. The default value is tolD=1e-6.
#' @example examples/MTM.R
#' @return List containing estimated posterior means and estimated posterior
#'   standard deviations, including: $yHat.
#' @export
MTM_2 <- function(Y, XF = NULL, K = NULL,
                resCov = list(type = "UN", df0 = 0, S0 = diag(0,ncol(as.matrix(Y)))),
                nIter = 110, burnIn = 10, thin = 2, saveAt = "", tolD = 1e-05) {
    
    if ((nIter - burnIn - thin) < 0) {
        stop("nIter must be greater than thin+burnIn")
    }
    iter <- 0
    
    Y <- as.matrix(Y)
    traits <- ncol(Y)
    n <- nrow(Y)
    hasXF <- !is.null(XF)
    hasK <- !is.null(K)
    
    ## Initializing overall mean
    mu <- colMeans(Y, na.rm = TRUE)
    post_mu <- rep(0, traits)
    post_YHat <- matrix(nrow = n, ncol = traits, 0)
    post_logLik <- 0
    
    ## Initializing missing values
    YStar <- Y
    whichNa <- list()
    whichNa$subjects <- which(apply(FUN = any, X = is.na(Y), MARGIN = 1))
    nNa <- length(whichNa$subjects)
    whichNa$traits <- list()
    if (nNa > 0) {
        for (k in 1:nNa) {
            whichNa$traits[[k]] <- which(is.na(Y[whichNa$subjects[[k]], ]))
            tmpSubject <- whichNa$subject[[k]]
            tmpTraits <- whichNa$traits[[k]]
            YStar[tmpSubject, tmpTraits] <- mu[tmpTraits]
        }
    }
    hasNa <- rowSums(is.na(Y)) > 0
    
    ## Initialization residuals
    E <- t(t(YStar) - mu)
    
    ## Initializing Fixed effects
    if (hasXF) {
        XF <- as.matrix(XF)
        dimX.f <- ncol(XF)
        tmp <- eigen(crossprod(XF))
        if (any(tmp$values < 0)) {
            Stop("XF is not full-column rank")
        } 
        if (length(tmp$values) == 1) {
            Tb.f <- tmp$vectors %*% matrix(1/sqrt(tmp$values))
        } else {
            Tb.f <- tmp$vectors %*% diag(1/sqrt(tmp$values))
        }
        XTb.f <- XF %*% Tb.f
        B.f <- matrix(nrow = dimX.f, ncol = traits, 0)
        for (i in 1:traits) {
            B.f[, i] <- lm(E[, i] ~ XF - 1)$coef
        }
        post_B.f <- B.f
        E <- E - XF %*% B.f
    }
    
    ## Initialization R0
    resCov$R <- var(E)/2
    resCov$L <- chol(resCov$R)
    resCov$RInv <- chol2inv(resCov$L)
    resCov$post_R <- matrix(nrow = traits, ncol = traits, 0)
    
    if (resCov$type == "REC") {
        resCov$B <- matrix(nrow = traits, ncol = traits, 0)
        resCov$PSI <- diag(resCov$R)
        resCov$post_PSI <- rep(0, traits)
        resCov$post_B <- matrix(nrow = traits, ncol = traits, 0)
    }
    
    if (resCov$type == "FA") {
        resCov$nF <- ncol(resCov$M)
        sdU <- apply(FUN = sd, MARGIN = 2, X = E/2)
        FA <- factanal(E/2, factors = resCov$nF)
        resCov$B <- matrix(nrow = traits, ncol = resCov$nF, 0)
        resCov$B[resCov$M] <- (diag(sdU) %*% FA$loadings)[resCov$M]
        resCov$PSI <- (sdU^2) * FA$uniquenesses + 1e-04
        resCov$R <- tcrossprod(resCov$B) + diag(resCov$PSI)
        resCov$L <- chol(resCov$R)
        resCov$RInv <- chol2inv(resCov$L)
        resCov$F <- matrix(nrow = n, ncol = resCov$nF, 0)
        resCov$post_PSI <- rep(0, traits)
        resCov$post_B <- matrix(nrow = traits, ncol = resCov$nF, 0)
    }
    
    if (resCov$type == "DIAG") {
        resCov$R <- diag(apply(FUN = var, X = E, MARGIN = 2))/2
        resCov$RInv <- diag(1/diag(resCov$R))
        resCov$post_R <- matrix(nrow = traits, ncol = traits, 0)
    }
    
    ## Initializing Kernel-components
    if (hasK) {
        post_K <- list()
        nK <- length(K)
        
        for (k in 1:nK) {
            K[[k]]$G <- resCov$R/nK
            if (is.null(K[[k]]$EVD)) {
                tmp <- eigen(K[[k]]$K)
            } else {
                tmp <- K[[k]]$EVD
            }
            
            K[[k]]$V <- tmp$vectors[, tmp$values > tolD]
            K[[k]]$d <- tmp$values[tmp$values > tolD]
            K[[k]]$nD <- length(K[[k]]$d)
            K[[k]]$U <- matrix(nrow = n, ncol = traits, 0)
            post_K[[k]] <- list()
            post_K[[k]]$G <- matrix(nrow = traits, ncol = traits, 0)
            post_K[[k]]$U <- matrix(nrow = n, ncol = traits, 0)
            
            if (K[[k]]$COV$type == "REC") {
                K[[k]]$B <- diag(traits)
                K[[k]]$PSI <- diag(K[[k]]$G)
                post_K[[k]]$B <- matrix(nrow = traits, ncol = traits, 0)
                post_K[[k]]$PSI <- rep(0, traits)
            }
            
            if (K[[k]]$COV$type == "FA") {
                K[[k]]$COV$nF <- ncol(K[[k]]$COV$M)
                sdU <- sqrt(diag(K[[k]]$G))
                FA <- factanal(covmat = K[[k]]$G, factors = K[[k]]$COV$nF)
                K[[k]]$B <- matrix(nrow = traits, ncol = K[[k]]$COV$nF, 0)
                K[[k]]$B[K[[k]]$COV$M] <- (diag(sdU) %*% FA$loadings)[K[[k]]$COV$M]
                K[[k]]$PSI <- (sdU^2) * FA$uniquenesses + 1e-04
                K[[k]]$G <- tcrossprod(K[[k]]$B) + diag(K[[k]]$PSI)
                K[[k]]$F <- matrix(nrow = K[[k]]$nD, ncol = K[[k]]$COV$nF, 0)
                post_K[[k]]$PSI <- rep(0, traits)
                post_K[[k]]$B <- matrix(nrow = traits, ncol = K[[k]]$COV$nF, 0)
            }
        }
    }
    
    ### Gibbs Sampler
    time <- proc.time()[3]
    for (i in 1:nIter) {
        
        logLik <- 0
        ## Fixed effects
        if (hasXF) {
            E <- E + XF %*% B.f
            B.f <- sampleBf(Y = E, XTb = XTb.f, Tb = Tb.f, R = resCov$R, L = resCov$L,
                            RInv = resCov$RInv, traits = traits, dimX = dimX.f)
            E <- E - XF %*% B.f
        }
        
        ## Kernels
        if (hasK)
        {
            for (k in 1:nK) {
                E <- E + K[[k]]$U
                GInv <- chol2inv(chol(K[[k]]$G))
                tmp <- sampleU(Y = E, RInv = resCov$RInv, GInv = GInv, V = K[[k]]$V,
                               d = K[[k]]$d, traits = traits, df0 = K[[k]]$df0, S0 = K[[k]]$S0)
                E <- E - tmp$U
                K[[k]]$U <- tmp$U
                
                if (K[[k]]$COV$type == "UN") {
                    tmp <- tmp$U0/sqrt(K[[k]]$d)
                    SS <- crossprod(tmp) + K[[k]]$COV$S0
                    df <- K[[k]]$nD + K[[k]]$COV$df0
                    K[[k]]$G <- MCMCpack::riwish(S = SS, v = df)
                }
                
                if (K[[k]]$COV$type == "REC") {
                    tmp <- tmp$U0/sqrt(K[[k]]$d)
                    
                    tmp <- sampleG.REC(U = tmp, M = K[[k]]$COV$M, PSI = K[[k]]$PSI,
                                       traits = traits, df0 = K[[k]]$COV$df0, S0 = K[[k]]$COV$S0,
                                       priorVar = K[[k]]$COV$var)
                    K[[k]]$G <- tmp$G
                    K[[k]]$B <- tmp$B
                    K[[k]]$PSI <- tmp$PSI
                }
                
                if (K[[k]]$COV$type == "FA") {
                    tmp <- tmp$U0/sqrt(K[[k]]$d)
                    tmp <- sampleG.FA(U = tmp, F = K[[k]]$F, M = K[[k]]$COV$M, B = K[[k]]$B,
                                      PSI = K[[k]]$PSI, G = K[[k]]$G, traits = traits, nF = K[[k]]$COV$nF,
                                      df0 = K[[k]]$COV$df0, S0 = K[[k]]$COV$S0, priorVar = K[[k]]$COV$var,
                                      n = K[[k]]$nD)
                    K[[k]]$G <- tmp$G
                    K[[k]]$PSI <- tmp$PSI
                    K[[k]]$B <- tmp$B
                    K[[k]]$F <- tmp$F
                }
                if (K[[k]]$COV$type == "DIAG") {
                    tmp <- tmp$U0/sqrt(K[[k]]$d)
                    K[[k]]$G <- sampleG.DIAG(U = tmp, traits = traits, df0 = K[[k]]$COV$df0,
                                             S0 = K[[k]]$COV$S0, n = n)
                }
            }
            
        }  #ends has(K)
        
        ## Overal Mean
        E <- t(t(E) + mu)
        mu <- sampleMu(Y = E, L = resCov$L, n = n, traits = traits)
        E <- t(t(E) - mu)
        ## Residual Variance
        if (resCov$type == "UN") {
            SS <- crossprod(E) + resCov$S0
            df <- n + resCov$df0
            resCov$R <- MCMCpack::riwish(v = df, S = SS)
            resCov$L <- chol(resCov$R)
            resCov$RInv <- chol2inv(resCov$L)
        }
        if (resCov$type == "REC") {
            tmp <- sampleG.REC(U = E, M = resCov$M, PSI = resCov$PSI, traits = traits,
                               df0 = resCov$df0, S0 = resCov$S0, priorVar = resCov$var)
            resCov$R <- tmp$G
            resCov$L <- chol(resCov$R)
            resCov$RInv <- chol2inv(resCov$L)
            resCov$B <- tmp$B
            resCov$PSI <- tmp$PSI
        }
        if (resCov$type == "FA") {
            
            tmp <- sampleG.FA(U = E, F = resCov$F, M = resCov$M, B = resCov$B, PSI = resCov$PSI,
                              G = resCov$R, traits = traits, nF = resCov$nF, df0 = resCov$df0,
                              S0 = resCov$S0, priorVar = resCov$var, n = n)
            resCov$R <- tmp$G
            resCov$L <- chol(resCov$R)
            resCov$RInv <- chol2inv(resCov$L)
            resCov$PSI <- tmp$PSI
            resCov$B <- tmp$B
            resCov$F <- tmp$F
        }
        
        if (resCov$type == "DIAG") {
            resCov$R <- sampleG.DIAG(U = E, traits = traits, df0 = resCov$df0, S0 = resCov$S0,
                                     n = n)
        }
        resCov$RInv <- chol2inv(chol(resCov$R))
        
        ## Imputing missing values
        
        YHat <- YStar - E
        if ((nNa > 0)) {
            for (j in 1:nNa) {
                subject <- whichNa$subject[[j]]
                missing <- whichNa$traits[[j]]
                observed <- (1:traits)[-missing]
                tmp <- sampleY(R = resCov$R, y = Y[subject, ], yHat = YHat[subject,
                ], e = E[subject, ], missing = missing, observed = observed, traits = traits)
                YStar[subject, ] <- tmp$y
                E[subject, ] <- YStar[subject, ] - YHat[subject, ]
                logLik <- logLik + tmp$logLik
            }
        }
        
        
        ## Completing the logLik computation
        if (nNa < n) {
            logLik <- logLik + dMVNorm(X = E[!hasNa, ], Sigma = resCov$R, mu = rep(0,
                                                                                   traits), log = TRUE)
        }
        
        ## Running means
        if ((i > burnIn) & (i%%thin == 0)) {
            iter <- iter + 1
            k <- (iter - 1)/(iter)
            
            
            post_mu <- post_mu * k + mu/iter
            resCov$post_R <- resCov$post_R * k + resCov$R/iter
            post_YHat <- post_YHat * k + YHat/iter
            post_logLik <- post_logLik * k + logLik/iter
            
            if (resCov$type %in% c("REC", "FA")) {
                resCov$post_B <- resCov$post_B * k + resCov$B/iter
                resCov$post_PSI <- resCov$post_PSI * k + resCov$PSI/iter
            }
            
            if (hasXF) {
                post_B.f <- post_B.f * k + B.f/iter
            }
            
            if (hasK) {
                for (j in 1:nK) {
                    post_K[[j]]$U <- post_K[[j]]$U * k + K[[j]]$U/iter
                    post_K[[j]]$G <- post_K[[j]]$G * k + K[[j]]$G/iter
                    
                    if (K[[j]]$COV$type %in% c("REC", "FA")) {
                        post_K[[j]]$B <- post_K[[j]]$B * k + K[[j]]$B/iter
                        post_K[[j]]$PSI <- post_K[[j]]$PSI * k + K[[j]]$PSI/iter
                    }
                }
            }
            
        }
        
        ## Saving Sammples
        tmp <- i%%thin == 0
        if ((tmp)) {
            
            tmp <- logLik
            fileName <- paste(saveAt, "logLik.dat", sep = "")
            write(tmp, ncol = length(tmp), file = fileName, append = T, sep = " ")
            
            
            tmp <- MCMCpack::vech(resCov$R)
            fileName <- paste(saveAt, "R.dat", sep = "")
            write(tmp, ncol = length(tmp), file = fileName, append = T, sep = " ")
            
            if (resCov$type %in% c("REC", "FA")) {
                if (sum(resCov$M) > 0) {
                    tmp <- resCov$B[resCov$M]
                    fileName <- paste(saveAt, "B_R.dat", sep = "")
                    write(tmp, ncol = length(tmp), file = fileName, append = T, sep = " ")
                }
                tmp <- resCov$PSI
                fileName <- paste(saveAt, "PSI_R.dat", sep = "")
                write(tmp, ncol = length(tmp), file = fileName, append = T, sep = " ")
                
            }
            
            tmp <- c(mu)
            fileName <- paste(saveAt, "mu", ".dat", sep = "")
            write(tmp, ncol = length(tmp), file = fileName, append = T, sep = " ")
            
            if (hasXF) {
                tmp <- as.numeric(B.f)
                fileName <- paste(saveAt, "mu", ".dat", sep = "")
                write(tmp, ncol = length(tmp), file = fileName, append = T, sep = " ")
            }
            
            if (hasK) {
                for (k in 1:nK) {
                    tmp <- MCMCpack::vech(K[[k]]$G)
                    fileName <- paste(saveAt, "G_", k, ".dat", sep = "")
                    write(tmp, ncol = length(tmp), file = fileName, append = T, sep = " ")
                    
                    if ((K[[k]]$COV$type %in% c("REC", "FA"))) {
                        tmp <- K[[k]]$PSI
                        fileName <- paste(saveAt, "PSI_G_", k, ".dat", sep = "")
                        write(tmp, ncol = length(tmp), file = fileName, append = T, sep = " ")
                        if (sum(resCov$M) > 0) {
                            tmp <- t(K[[k]]$B)[t(K[[k]]$COV$M)]
                            fileName <- paste(saveAt, "B_G_", k, ".dat", sep = "")
                            write(tmp, ncol = length(tmp), file = fileName, append = T,
                                  sep = " ")
                        }
                    }
                }
            }
        }
        
        tmp <- proc.time()[3]
        cat(paste("Iter: ", i, "time: ", (round(tmp - time, 4))))
        cat("\n")
        cat("\n")
        time <- tmp
    }
    
    tmp <- list()
    tmp$R <- resCov$post_R
    if (resCov$type %in% c("REC", "FA")) {
        tmp$B <- resCov$post_B
        tmp$PSI <- resCov$post_PSI
    }
    out <- list(mu = post_mu, YHat = post_YHat, resCov = tmp, E = E)
    
    if (hasXF) {
        out$B.f = post_B.f
    }
    
    if (hasK) {
        out$K <- post_K
    }
    
    out$DIC <- getDIC(Y = Y, YHat = post_YHat, R = resCov$post_R, meanLogLik = post_logLik)
    
    return(out)
}


### FIXED EFFECTS
### #############################################################################
sampleBf <- function(Y, XTb, Tb, R, L, RInv, traits, dimX) {
    SOL <- crossprod(XTb, Y)
    Z <- matrix(nrow = dimX, ncol = traits, rnorm(dimX * traits))
    E <- tcrossprod(Z, L)
    B <- SOL + E
    B <- Tb %*% B
    return(B)
}

sampleMu <- function(Y, L, n, traits) {
    sol <- colMeans(Y)
    L <- L/sqrt(n)
    mu <- as.numeric(crossprod(L, rnorm(traits)) + sol)
    return(mu)
}

### MISSING RECORDS
### #############################################################################
sampleY <- function(R, y, yHat, e, missing, observed, traits) {
    if (length(missing) == traits) {
        e <- crossprod(chol(R), rnorm(traits))
        y <- yHat + e
        logLik <- 0
    } else {
        Roo <- matrix(R[observed, observed], nrow = length(observed), ncol = length(observed))
        RooInv <- chol2inv(chol(Roo))
        Rmm <- matrix(R[missing, missing], nrow = length(missing), ncol = length(missing))
        Rom <- matrix(R[observed, missing], nrow = length(observed), ncol = length(missing))
        
        Bmo <- crossprod(Rom, RooInv)
        
        yHat2 <- as.numeric(Bmo %*% e[observed])
        CondVar <- Rmm - Bmo %*% Rom
        L <- chol(CondVar)
        e <- crossprod(L, rnorm(length(missing)))
        y[missing] <- yHat[missing] + yHat2 + e
        tmp <- (y - yHat)[observed]
        logLik <- dMVNorm_i(x_i = tmp, SigmaInv = RooInv, mu = rep(0, length(observed)),
                            log = TRUE)
    }
    out <- list(y = y, logLik = logLik)
    return(out)
}


### Sample U
### #############################################################################
sampleUj <- function(j, Y, RInv, GInv, V, d, traits) {
    CInv <- chol2inv(chol(RInv + GInv/d[j]))
    T <- RInv %*% CInv
    YStar <- Y %*% T
    sol <- as.numeric(crossprod(V[, j], YStar))
    L <- chol(CInv)
    uj <- as.numeric(crossprod(L, rnorm(traits))) + sol
    return(uj)
}

sampleU <- function(Y, RInv, GInv, V, d, traits, df0 = 0, S0 = 0) {
    tmp <- matrix(unlist(lapply(FUN = sampleUj, Y = Y, RInv = RInv, GInv = GInv,
                                V = V, d = d, traits = traits, X = 1:length(d))), byrow = TRUE, ncol = traits)
    
    U <- V %*% tmp
    return(list(U = U, U0 = tmp))
}



### Sample G
### #############################################################################


sampleG.REC <- function(U, M, PSI, traits, priorVar = 100, df0 = rep(0, traits),
                        S0 = rep(0, traits)) {
    ### Model: U=UB+D (ONLY RECURSIVE ALLOWED!)  Current sample of random effects
    ### ('data') T a pxp matrix with TRUE/FALSE indicating possition of non-null
    ### recursive effects (FALSE in diagonal!)  PSI px1 the variance of the orthogonal
    ### shocks ...
    B <- matrix(nrow = traits, ncol = traits, 0)
    for (i in 1:traits) {
        
        dimX <- sum(M[i, ])
        
        if (dimX > 0) {
            tmpX <- U[, M[i, ]]
            tmpY <- U[, i]
            
            C <- crossprod(tmpX)/PSI[i] + 1/priorVar
            CInv <- chol2inv(chol(C))
            rhs <- crossprod(tmpX, tmpY)/PSI[i]
            
            sol <- crossprod(CInv, rhs)
            L <- chol(CInv)
            shock <- crossprod(L, rnorm(dimX))
            tmpB <- as.numeric(sol + shock)
            B[i, M[i, ]] <- tmpB
            uStar <- tmpY - matrix(tmpX, ncol = dimX) %*% (tmpB)
            SS <- as.numeric(crossprod(uStar)) + S0[i]
            df <- nrow(U) + df0[i]
            PSI[i] <- SS/rchisq(n = 1, df = df)
        } else {
            SS <- as.numeric(crossprod(U[, i])) + S0[i]
            df <- nrow(U) + df0
            PSI[i] <- SS/rchisq(n = 1, df = df)
        }
    }
    
    tmp <- solve(diag(traits) - B)
    G <- tmp %*% diag(PSI) %*% t(tmp)
    
    out <- list(B = B, PSI = PSI, G = G)
    
    return(out)
}



############################################################################################
sampleG.FA <- function(U, F, M, B, PSI, G, traits, nF, n, df0 = rep(1, traits), S0 = rep(1/100,
                                                                                         traits), priorVar = 100) {
    ### Gibbs sampler for FA model Model: U=BF+D
    
    ## sampling common factors LOOP OVER FACTORS
    for (i in 1:nF) {
        tmpY <- U - F[, -i] %*% matrix((B[, -i]), ncol = traits)
        rhs <- tmpY %*% matrix(B[, i]/PSI, ncol = 1)
        CInv <- 1/(sum((B[, i]^2)/PSI) + 1)
        sol <- CInv * rhs
        SD <- sqrt(CInv)
        F[, i] <- rnorm(n = n, sd = SD, mean = sol)
    }
    
    # sampling loadings LOOP OVER TRAITS LOOP OVER FACTORS
    for (i in 1:traits) {
        for (j in 1:nF) {
            if (M[i, j]) {
                tmpY <- U[, i] - F[, -j] %*% matrix(B[i, -j], ncol = 1)
                CInv <- 1/as.numeric(crossprod(F[, j])/PSI[i] + 1/priorVar)
                rhs <- as.numeric(crossprod(F[, j], tmpY)/PSI[i])
                sol <- CInv * rhs
                SD <- sqrt(CInv)
                B[i, j] <- rnorm(n = 1, mean = sol, sd = SD)
            }
        }
        D <- U[, i] - F %*% B[i, ]
        df <- df0[i] + n
        SS <- S0[i] + crossprod(D)
        PSI[i] <- SS/rchisq(df = df, n = 1)
    }
    if (nF > 1) {
        B <- varimax(B)$loadings[]
    }
    G <- tcrossprod(B) + diag(PSI)
    out <- list(F = F, PSI = PSI, B = B, G = G)
}

########################################################################################################################

############################################################################################
sampleG.DIAG <- function(U, traits = ncol(U), n = nrow(U), df0 = rep(0, traits),
                         S0 = diag(0, traits)) {
    ### Gibbs sampler for DIAG covaraince matrices
    G <- matrix(nrow = traits, ncol = traits, 0)
    ## sampling common factors LOOP OVER FACTORS
    for (i in 1:traits) {
        tmp_SS <- sum(U[, i]^2) + S0[i]
        tmp_df <- n + df0[i]
        G[i, i] <- tmp_SS/rchisq(df = tmp_df, n = 1)
    }
    
    return(G)
}


if (FALSE) {
    B <- cbind(rnorm(1000, sd = 1), rnorm(1000, sd = 2))
    sampleG.DIAG(B)
    
}


########################################################################################################################
dScaleInvChisq <- function(x, df, S, log = FALSE) {
    y <- S/x
    out <- dchisq(y, df = df, log = TRUE) + log(S/(x^2))
    if (!log) {
        out <- exp(out)
    }
    return(out)
}

dMVNorm <- function(X, Sigma, mu, log = TRUE) {
    ## X may be a matrix, rows are IID draws
    
    TInv <- solve(chol(Sigma))
    Z <- crossprod((t(X) - mu), TInv)
    out <- sum(dnorm(x = as.vector(Z), log = TRUE)) + nrow(X) * sum(log(diag(TInv)))
    if (!log) {
        out <- exp(out)
    }
    return(out)
    
}


dMVNorm_i <- function(x_i, SigmaInv, mu, log = TRUE) {
    
    ## works for a single random draw and requires SigmaInv
    
    e <- as.matrix(x_i - mu)
    out <- -(length(e)/2 * log(2 * pi)) + log(det(SigmaInv))/2 - (crossprod(e, SigmaInv) %*%
                                                                      e)/2
    if (!log) {
        out <- exp(out)
    }
    return(out)
}



getDIC <- function(Y, YHat, R, meanLogLik = NULL, logLik = NULL) {
    #### Returns Deviance Information Criterion and Effective Number of Parameters Y
    #### (nxp) the data matrix (NA's for misssing) YHat the posterior mean of the
    #### conditional expectation logLik the samples of the log-likelihood genreated by
    #### MTM() meanLogLik the posterior mean of logLik
    if (!is.null(logLik)) {
        meanLogLik <- mean(logLik)
    }
    E <- Y - YHat
    logLikAtPostMean <- 0
    for (i in 1:nrow(Y)) {
        observed <- !is.na(Y[i, ])
        e <- matrix(nrow = 1, E[i, observed])
        mu <- rep(0, sum(observed))
        Rtmp <- R[observed, observed]
        if (length(Rtmp) > 0) {
            logLikAtPostMean <- logLikAtPostMean + dMVNorm(X = e, Sigma = Rtmp, mu = mu,
                                                           log = TRUE)
        }
    }
    
    pD <- -2 * (meanLogLik - logLikAtPostMean)
    DIC <- pD - 2 * meanLogLik
    return(list(meanLogLik = meanLogLik, logLikAtPostMean = logLikAtPostMean, DIC = DIC,
                pD = pD))
    
}
