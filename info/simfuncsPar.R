dSorensen2p11 <- function(d, p01, p10) {
  p11 <- 0.5 * (p10 + p01) * (1 - d) / d
  invalid.p11 <- ((p11 + p10 + p01) > 1) | (p11 < 0)
  if (any(invalid.p11)) {
    warning("Some resulting p11 not compatible with p10 or p01 and the dissimilarities range")
    p11 <- p11[!invalid.p11]
    d <- d[!invalid.p11[length(invalid.p11):1]]
  }
  attr(p11, "d") <- d
  return(p11)
}


getP11ToSim <- function(dRange, p01, p10, d0, nPoints = 21) {
  if (dRange[1] < dRange[2]) {
    dmin <- dRange[1]
    dmax <- dRange[2]
  } else {
    dmin <- dRange[2]
    dmax <- dRange[1]
  }
  # in fact, p11min > p11max, but they correspond to min and max dissimilarity

  if (!missing(d0)) {
    p11.0 <- dSorensen2p11(d0, p01, p10)
    if (length(p11.0) == 0) {
      stop("Equivalence threshold dissimilarity d0 to be simulated not compatible with p10 and p01")
    }
    nPoints <- nPoints - 1
    dToSim <- seq(from = dmin, to = dmax, length.out = nPoints)
    id0 <- findInterval(d0, dToSim)
    if (id0 == 0) {
      dToSim <- c(d0, dToSim)
    } else if (id0 == length(dToSim)) {
      dToSim <- c(dToSim, d0)
    } else {
      dToSim <- c(dToSim[1:id0], d0, dToSim[(id0+1):length(dToSim)])
    }
    result <- dSorensen2p11(dToSim, p01, p10)
    attr(result, "d0") <- d0
  } else {
    dToSim <- seq(from = dmin, to = dmax, length.out = nPoints)
    result <- dSorensen2p11(dToSim, p01, p10)
    attr(result, "d0") <- NULL
  }
  attr(result, "p10") <- p10
  attr(result, "p01") <- p01
  # attr(result, "dSorensen") <- dToSim
  return(result)
}

buildScenarios <- function(d0ToSim, nToSim, p01p10ToSim, dRange, nPoints = 21) {
  result <- lapply(d0ToSim, function(d0) {
    # cat("d0 = ", d0, "\n")
    result <- lapply(nToSim, function(n, d0) {
      # cat("d0 = ", d0, "  n = ", n, "\n")
      apply(p01p10ToSim, 1, function(p01p10, n, d0) {
        # cat("d0 = ", d0, "  n = ", n, "p01p10 = ", p01p10, "\n")
        return(list(d0 = d0, n = n, p11 = getP11ToSim(dRange, p01p10[1], p01p10[2], d0, nPoints)))
      }, n = n, d0 = d0)
    }, d0 = d0)
    names(result) <- paste("n_", nToSim, sep = "")
    return(result)
  })
  names(result) <- paste("d0_", round(d0ToSim, 4), sep = "")
  return(result)
}

simScenario <- function(scenario, nsim = 100000, 
                        seed = NULL, boot.seed = NULL, resultsFolder = ".\\") {
  fileId <- paste(resultsFolder,
                  paste(paste("d0=", round(scenario$d0, 4), "n=", scenario$n, 
                              "p01=", attr(scenario$p11, "p01"), "p10=", attr(scenario$p11, "p10"), 
                              "nsim=", nsim, sep = ""),
                        "Rda", sep = "."), sep = "")
  cat(fileId, "\n")
  
  cl <- makeCluster(4)
  clusterExport(cl, "equivTestSorensen")
  
  len.d <- length(attr(scenario$p11, "d"))
  prRej <- rbind(attr(scenario$p11, "d"), rep(NA, len.d), rep(NA, len.d))
  rownames(prRej) <- c("d", "pr{Rej}", "prB{Rej}")
  ip11 <- 0
  for (p11 in scenario$p11) {
    ip11 <- ip11 + 1
    cat("d=", attr(scenario$p11, "d")[ip11], "\n")
    if (!is.null(seed)) {
      set.seed(seed)
    }
    simRes <- simSorensen(nsim = nsim, n = scenario$n, 
                          p = c(p11, attr(scenario$p11, "p01"), attr(scenario$p11, "p10")),
                          d0 = scenario$d0, boot.seed = boot.seed,
                          cl = cl)
    print(formatC(simRes, format = "g", digits = 10))
    prRej[2, ip11] <- simRes["pr(Rej)"]
    prRej[3, ip11] <- simRes["prB(Rej)"]
  }
  
  stopCluster(cl)
  
  save(prRej, file = fileId)
  return(prRej)
}

iterScenarios <- function(scenarios, nsim = 100000, seed = NULL, boot.seed = NULL,
                          resultsFolder = ".\\",
                          simId = "results") {
  prRejScenarios <- lapply(scenarios, function(this.d0) {
    prRejd0 <- lapply(this.d0, function(this.n) {
      prRejn <- lapply(this.n, function(this.p01p10) {
        return(simScenario(this.p01p10, nsim = nsim, 
                           seed = seed, boot.seed = boot.seed,
                           resultsFolder = resultsFolder))
      })
      names(prRejn) <- names(this.n)
      return(prRejn)
    })
    names(prRejd0) <- names(this.d0)
    return(prRejd0)
  })
  names(prRejScenarios) <- names(scenarios)
  save(prRejScenarios, file = paste(paste(resultsFolder, simId, sep = ""), "Rda", sep = "."))
  return(prRejScenarios)
}

simSorensen <- function(nsim = 100000, n, p, 
                        d0 = 1/(1 + 2*1.25), 
                        conf.level = 0.95, boot.seed = NULL,
                        cl) {
  simRun <- function(isim, freqs, d0, conf.level, boot.seed) {
    fr <- matrix(freqs[,isim], ncol = 2)
    equiv <- equivTestSorensen(fr, d0 = d0, conf.level = conf.level,
                               check.table = FALSE)
    if (any(!is.finite(c(equiv$estimate, equiv$stderr))) || equiv$stderr == 0) {
      return(c(NA, NA, NA, NA))
    } 
    if (!is.null(boot.seed)) {
      set.seed(boot.seed)
    }
    equivBoot <- equivTestSorensen(fr, d0 = d0, conf.level = conf.level,
                                   boot = TRUE,
                                   check.table = FALSE)
    return(c(equiv$estimate, equiv$stderr, equiv$p.value, equivBoot$p.value))
  }
  if (length(p) != 3) {
    stop("Argument p must have length 3")
  }
  if ((sum(p) > 1) || any(p < 0)) {
    stop("Invalid probability vector to simulate")
  }
  # pEnrich <- sum(p)
  nu <- n * p
  alpha <- 1 - conf.level
  dS <- dSorensen(p)
  result <- c(NA, alpha, dS, d0, n, p, 0, 0, #NA, 
              0, 0, 0, nu, NA, NA)
  names(result) <- c("nsim", "alpha", "dS", "d0", "n", "p11", "p01", "p10", 
                     "E(^dS)", "sd(^dS)", #"sqrt(4p11(p10+...", 
                     "E(^se)", "pr(Rej)", "prB(Rej)", 
                     "E(nu11)", "E(nu01)", "E(nu10)", 
                     "bias(^dS)", "bias(^se)")
  print(date())
  freqs <- rmultinom(nsim, n, c(p, 1 - sum(p)))
  # cl <- makeCluster(4)
  # clusterExport(cl, "equivTestSorensen")
  testResults <- parSapply(cl, 1:nsim, simRun, 
                           freqs = freqs, d0 = d0, conf.level = conf.level,
                           boot.seed = boot.seed)
  # stopCluster(cl)
  print(date())
  result["nsim"] <- nsim - sum(apply(testResults, 2, function(testRes) {
    any(is.na(testRes))
  }))
  result["E(^dS)"] <- mean(testResults[1,], na.rm = TRUE)
  result["sd(^dS)"] <- sd(testResults[1,], na.rm = TRUE)
  result["E(^se)"] <- mean(testResults[2,], na.rm = TRUE)
  result["pr(Rej)"] <- mean(testResults[3,] <= alpha, na.rm = TRUE)
  result["prB(Rej)"] <- mean(testResults[4,] <= alpha, na.rm = TRUE)
  result["bias(^dS)"] <- result["E(^dS)"] - result["dS"]
  result["bias(^se)"] <- round(result["E(^se)"] - result["sd(^dS)"], digits = 8)
  return(result)
}

bias.dS <- function(p, n) {
  dS <- dSorensen(p)
  return(dS * dS * p[1] / (n * (sum(p) + p[1])))
}

