############################
## Main HTKmeans function ##
############################

HTKmeans <- function(X, k, lambdas,
                     standardize = TRUE,
                     iter.max = 100,
                     nstart = 100) {
  # Regularized k-means for a fixed value of k over a grid of lambda
  # values 
  #
  # args:
  #   X: a n x p data matrix
  #   lambdas: a (sequence of) value(s) for the regularization parameter lambda
  #   iter.max: the maximum number of iterations
  #   nstart: number of starts used by the regular k-means algorithm
  # 
  # returns:
  #   HTKmeans.out: a list with components:
  #   centers: the centers 
  #   cluster: the labels of the data points
  # 
  
  inputargs <- list(k = k,
                    lambdas = lambdas,
                    iter.max = iter.max,
                    nstart = nstart,
                    standardize = standardize,
                    X = X)
  if (standardize) {
    X <- scale(X)
  }
  HTKmeans.out <- list()
  lambdas     <- sort(lambdas)
  
  for (i in 1:length(lambdas)) {
    startval.out <- getStartvalues(X = X, k = k,
                                   lambda = lambdas[i],
                                   iter.max = iter.max,
                                   nstart = nstart)
    centers    <- startval.out$centers
    IDs        <- startval.out$IDs
    
    HTKmeans.out[[i]] <- HTKmeans_inner(X = X,
                                        centers = centers,
                                        lambda =  lambdas[i],
                                        iter.max = iter.max, 
                                        oldIDs = IDs)
    if (max(abs(HTKmeans.out[[i]]$centers)) < 1e-10) {
      HTKmeans.out[[i]]$cluster <- rep(0, dim(X)[1])
      if (i < length(lambdas)) {
        for (j in (i + 1):length(lambdas)) {
          HTKmeans.out[[j]] <- HTKmeans.out[[i]]
        }
        break
      }
    }
  }
  
  HTKmeans.out$inputargs <- inputargs
  return(HTKmeans.out)
}


####################################
## Support functions for HTKmeans ##
####################################

getIDs <- function(X, lambda, centers, oldIDs) {
  # calculate the new cluster assignments
  #
  
  n        <- dim(X)[1]
  IDs      <- rep(0, n)
  
  # assign new IDs
  for (i in 1:n) {
    x      <- X[i, ]
    cost   <- apply(centers, 1, function(y) sum((y - x)^2)) / n
    IDs[i] <- which.min(cost)
  }
  return(IDs)
}

getObj <- function(X, IDs, centers, lambda) {
  # calculates the value of the objective function
  # 
  
  n            <- dim(X)[1]
  obj          <- 0 # full objective function
  obj_penalty  <- 0 # penalty part of the objective function
  WCSS         <- 0 # within cluster sums of squares, scaled by n
  WCSS_nonZero <- 0 # WCSS for the non-zero variables, scaled by n
  
  nonzeroCenters <- which(apply(centers, 2, function(y) sum(abs(y))) > 0)
  
  if (length(nonzeroCenters) == 0) {
    obj <- WCSS <- sum(X^2) / n
    WCSS_nonZero <- 0
    obj_penalty <- 0
  } else {
    
    for (i in 1:n) {
      x     <- X[i, ]
      id    <- IDs[i]
      WCSS_ctb     <- sum((x - centers[id, ])^2) / n
      WCSS_nZ_ctb  <- sum((x[nonzeroCenters] - centers[id, nonzeroCenters])^2) / n
      WCSS         <- WCSS + WCSS_ctb
      WCSS_nonZero <- WCSS_nonZero + WCSS_nZ_ctb
    }
    
    lambda <- rep(lambda, dim(centers)[2])
    obj_penalty <- sum(lambda * (apply(centers, 2, function(y) sum(y^2)) > 0))
    obj <- obj_penalty + WCSS
  }
  return(list(obj = obj,
              obj_penalty = obj_penalty,
              WCSS = WCSS,
              WCSS_nonZero = WCSS_nonZero, 
              nBactive = length(nonzeroCenters), 
              activeVars = nonzeroCenters))
  
}

getClassicCenters <- function(X, clusID, k) {
  # calculate classical centers
  #
  
  muhat <- sapply(1:k, FUN = function(y) colMeans(X[which(clusID == y),
                                                    , drop = FALSE]))
  return(t(muhat))
}

getCenters <- function(X, clusID, k,
                       lambda,
                       iter.max = 100) {
  # Updates the centers for the regularized k-means, assuming fixed
  # assignment of the observations to k clusters
  # args:
  #   X: n x p data matrix
  #   clusID: vector of length n with cluster memberships in 1, ..., k
  #   k: number of clusters
  #   lambda: penalization parameter
  # returns:
  #   centers: k x p matrix of updated cluster centers
  #
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  centers <- matrix(0, k, p)
  
  centers_clas <- getClassicCenters(X = X, clusID = clusID, k = k)
  lambdaMat    <- matrix(lambda, k, p)
  
  for (j in 1:p) {
    muhat        <- centers_clas[, j]
    withinSS     <- sum(sapply(1:k, FUN = function(y)
      sum((X[which(clusID == y), j] - muhat[y])^2)))
    centers[, j] <- (sum(X[, j]^2)  > (withinSS + n * lambdaMat[1, j] )) * muhat
  }
  
  return(centers)
}

getStartvalues <- function(X, k, lambda, 
                           iter.max,
                           nstart) {
  # returns initial centers
  #
  
  if (lambda == 0) {
    kmeans.out <- kmeans(X, k, nstart = nstart)
    centers    <- kmeans.out$centers
    IDs <- kmeans.out$cluster
  } else {
    km.out   <- kmeans(X, centers = k, nstart = nstart)
    meansize <- apply(km.out$centers, 2, function(y) sum(y^2))
    
    pervars  <- c(1, 2, 5, 10, 25, 50, 100) / 100
    objvals  <- rep(Inf, length(pervars))
    ordering <- order(meansize, decreasing = TRUE)
    centers  <- IDs <- list()
    for (i in 1:length(pervars)) {
      Inds <- ordering[1:(pervars[i] * dim(X)[2])]
      if (length(Inds) > 0) {
        if (i == length(pervars)) {
          kmeans.out <- km.out
        } else {
          kmeans.out <- kmeans(X[, Inds], k, nstart = nstart)
        }
        centers[[i]]  <- getCenters(X = X, clusID = kmeans.out$cluster, k = k,
                                    lambda = lambda,
                                    iter.max = iter.max)
        IDs[[i]]      <- kmeans.out$cluster
        objvals[i]    <- getObj(X = X, IDs = kmeans.out$cluster,
                                centers = centers[[i]], lambda = lambda)$obj
      }
    }
    centers <- centers[[which.min(objvals)]]
    IDs     <- IDs[[which.min(objvals)]]
  }
  
  return(list(centers = centers, IDs = IDs))
}

HTKmeans_inner <- function(X, centers, lambda,
                           iter.max, oldIDs) {
  # iterates the two steps of Lloyds algorithm
  #
  
  n         <- dim(X)[1]
  converged <- FALSE
  clusIDOld <- rep(0, n)
  itnb      <- 0
  k         <- dim(centers)[1]
  clusID    <- oldIDs
  
  # start iteration
  while ((itnb < iter.max) & !converged ) {
    # Get new cluster assignments
    clusID <- getIDs(X, lambda, centers,
                     oldIDs = clusID)
    
    if (length(unique(clusID)) != k) {
      centers <- matrix(0, k, dim(X)[2])
      converged <- FALSE
      break
    }
    
    # get new centers
    centers  <- getCenters(X = X, clusID = clusID, k = k,
                           lambda = lambda)
    converged <- sum(abs(clusID - clusIDOld)) == 0
    clusIDOld <- clusID
    itnb      <- itnb + 1
  }
  return(list(centers = centers,
              cluster = clusID,
              itnb = itnb,
              converged = converged))
}



#############################################
## Functions for analyzing HTKmeans output ##
#############################################

diagPlot = function(HTKmeans.out, type = 1)  {
  # Make the diagnostic plot based on the output of HTKmeans
  #
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  if (type == 1) {
  inputargs <- HTKmeans.out$inputargs
  lambdas <- inputargs$lambdas
  
  plotdata <- sapply(1:length(lambdas),
                     function(z)
                       colSums(abs(HTKmeans.out[[z]]$centers)))
  xlims <- rev(range(inputargs$lambdas))
  
  matplot(lambdas, t(plotdata), type = "l",
          xlab = expression(lambda), ylab = "norm of center vector",
          col = 1:ncol(plotdata), xlim = xlims, lwd = 4, lty = 1,
          cex.lab = 2, cex.axis = 2)
  
  } else {
  
  
  clusterInfo <- extractClusterInfo(HTKmeans.out,  summarize = TRUE)
  IDmat <- matrix(unlist(clusterInfo$IDs),
                  ncol = length(clusterInfo$IDs), byrow = FALSE)
  
  deltaARI <- rep(0, dim(IDmat)[2] - 1)
  for (i in 1:(dim(IDmat)[2] - 1)) {
    deltaARI[i] <-  1 - mclust::adjustedRandIndex(IDmat[, i+1], IDmat[,i])
  }
  
  par(mar=c(3.5, 5, 1.5, 6))
  
  plot(rev(clusterInfo$nbActiveVars)[-1],
       diff(rev(clusterInfo$WCSS_active)) /
         diff(rev(clusterInfo$nbActiveVars)),
       ylab = "",
       xlab = "",
       type = "l", lwd = 3,axes = FALSE,
       ylim = c(0, 1))
  axis(2, ylim=c(0,1),col="black",las=1, lwd = 3,
       lwd.ticks = 2, cex.axis = 2)  ## las=1 makes horizontal labels
  mtext(expression(paste(Delta," WCSS")), side=2, line=3.5, cex = 2)
  
  
  ## Allow a second plot on the same graph
  par(new=TRUE)
  
  # plot delta ARI
  plot(clusterInfo$nbActiveVars[-c(length(deltaARI), length(deltaARI)-1)],
       deltaARI[-c(length(deltaARI))],
       xlab="", ylab="", ylim = c(0, 1),
       axes=FALSE, type="l", xlim = c(1,max(clusterInfo$nbActiveVars)),
       col="firebrick", lwd = 3)
  mtext(expression(paste(Delta," ARI")),side=4,col="firebrick",line=4, cex = 2) 
  axis(4, ylim=c(0, 1), col="firebrick",col.axis="firebrick",
       las=1, cex.axis = 2, 
       lwd = 3, lwd.ticks = 2, at = seq(0,1, by = 0.2),
       labels = seq(0,1, by = 0.2))
  
  ## Draw the x axis
  axis(1,seq(max(clusterInfo$nbActiveVars), 1, -1),lwd = 3,
       lwd.ticks = 2, cex.axis = 2)
  mtext("Number of active variables",
        side=1,col="black",line=2.5, cex = 2)  
  }
}

extractClusterInfo <- function(HTKmeans.out,
                               y = NULL,
                               summarize = FALSE) {
  # Extract and organize cluster info from the output of HTKmeans 
  # can be used to make diagnostic plots
  #
  
  lambdas    <- HTKmeans.out$inputargs$lambdas
  objective  <- penalty <- WCSS <- WCSS_nonZero  <- nbActiveVars <- ARIs <- rep(0, length(lambdas))
  activeVars <- list()
  centers    <- list()
  IDs        <- list()
  
  if (HTKmeans.out$inputargs$standardize) {
    Z <- scale(HTKmeans.out$inputargs$X)
  } else {
    Z <- HTKmeans.out$inputargs$X
  }
  for (i in 1:length(lambdas)) {
    temp <- getObj(X = Z, IDs = HTKmeans.out[[i]]$cluster, 
                   centers = HTKmeans.out[[i]]$centers,
                   lambda = lambdas[i])
    WCSS[i]         <- temp$WCSS
    objective[i]    <- temp$obj
    penalty[i]      <- temp$obj_penalty
    WCSS_nonZero[i] <- temp$WCSS_nonZero
    nbActiveVars[i] <- temp$nBactive
    activeVars[[i]] <- temp$activeVars
    IDs[[i]]        <- HTKmeans.out[[i]]$cluster
    centers[[i]] <-  HTKmeans.out[[i]]$centers
    if (!is.null(y)) {
      ARIs[i] <- mclust::adjustedRandIndex(HTKmeans.out[[i]]$cluster, y)
    }
  }
  
  if (summarize) {
    uniqNbVars <- sort(unique(nbActiveVars), decreasing = TRUE)
    objective_sum <- penalty_sum <- WCSS_sum <- WCSS_nonZero_sum  <- nbActiveVars_sum <- ARIs_sum <- rep(0, length(uniqNbVars))
    activeVars_sum <- centers_sum <- IDs_sum <- list()
    lambdas_sum <- rep(0, length(uniqNbVars))
    for (j in 1:length(uniqNbVars)) {
      simInds <- which(nbActiveVars == uniqNbVars[j])
      
      best.one            <- simInds[which.min(objective[simInds])]
      WCSS_sum[j]         <- WCSS[best.one]
      objective_sum[j]    <- objective[best.one]
      penalty_sum[j]      <- penalty[best.one]
      WCSS_nonZero_sum[j] <- WCSS_nonZero[best.one]
      nbActiveVars_sum[j] <- nbActiveVars[best.one]
      activeVars_sum[[j]] <- activeVars[[best.one]]
      centers_sum[[j]]    <- centers[[best.one]]
      lambdas_sum[j]      <- lambdas[best.one]
      IDs_sum[[j]]        <- IDs[[best.one]]
      if (!is.null(y)) {
        ARIs_sum[j] <- ARIs[best.one]
      }
    }
    
    objective    <- objective_sum
    penalty      <- penalty_sum
    WCSS         <- WCSS_sum
    WCSS_nonZero <- WCSS_nonZero_sum
    nbActiveVars <- nbActiveVars_sum
    activeVars   <- activeVars_sum
    ARIs         <- ARIs_sum
    lambdas      <- lambdas_sum
    centers      <- centers_sum
    IDs          <- IDs_sum
  }
  
  
  return(list(objective = objective,
              penalty = penalty,
              WCSS = WCSS,
              WCSS_active = WCSS_nonZero, 
              nbActiveVars = nbActiveVars,
              activeVars = activeVars,
              ARIs = ARIs,
              lambdas = lambdas,
              centers = centers,
              IDs = IDs))
}

##################################################
## Functions for selecting the tuning parameter ##
##################################################

getLambda <- function(HTKmeans.out, type  = "AIC") {
  # Select the value of the regularization parameter lambda
  # based on the output of HTKmeans with hard thresholding
  #
  
  clusterInfo <- extractClusterInfo(HTKmeans.out,
                                    y = NULL,
                                    summarize = TRUE)
  k       <- HTKmeans.out$inputargs$k
  lambdas <- clusterInfo$lambdas
  n       <- nrow(HTKmeans.out$inputargs$X)
  if (type == "AIC") {
    AICvals <- clusterInfo$WCSS * n +
      2 * k * clusterInfo$nbActiveVars
    lambda <- lambdas[which.min(AICvals)]
  } else if (type == "BIC") {
    BICvals <- clusterInfo$WCSS * n +
      k * log(n) * clusterInfo$nbActiveVars
    lambda <- lambdas[which.min(BICvals)]
  }
  
  return(lambda)
}

