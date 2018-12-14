#' Generate the structure/penalty matrix for sarim formulae
#' 
#' @description Function used in definition of smooth terms within sarim model formulae.
#' 
#' @param formula usual formula, y ~ x1 + sx(x2), be aware an intercept can only be
#'                included by setting intercept = "TRUE".
#' @param data a data.frame.
#' @param intercept decide if an intercept should be included or not. 
#' @param nIter number of MCMC-Iteration.
#' @param burnin number of Burn-In iterations (currently NOT available).
#' @param family family can currently only be "gaussian", "binomial" or "poisson".
#' @param link link-function, currently only be "identity", "logit" or "log" implemented,
#'             if family = "binomial", please choose "logit", if family = "poisson", please choose "log".
#' @param sigma start value for the dispersion parameter, it is assumed that sigma ~ IG(sigma_a, sigma_b),
#'              i.e. it is assumed that sigma is Inverse Gamma distributed with parameters sigma_a and sigma_b.
#' @param sigma_a prior-value, see sigma.
#' @param sigma_b prior-value, see sigma.
#' @param Ntrials Number of trails, especially for binomial distribution.
#' @param m number of maximal Lanczos-iteration for start iteration, after first iteration, it use the 
#'              number of the first + 10, in the third iteration it use the number of the second + 10; 
#'              for saving memory.
#' @param thr threshold, when the Lanczos-algorithm or conjugate gradient-algorithm should break.
#' @import Matrix parallel
#' @export
#' @useDynLib Sarim
sarim <- function(formula, data = list(), intercept = "FALSE", nIter = 1000L, burnin = 100L,
                  family = "gaussian", link = "identity",
                  sigma = 0.1, sigma_a = 0.0001, sigma_b = 0.0001, Ntrials = 1L,
                  m = 250L, thr = 0.0001) {
  
    nIter =50L
    burnin = 50L
    mf <- stats::model.frame(formula = formula, data = data)
    y <- as.numeric(stats::model.response(mf))
    
    # initialise parameterlist for MCMC-sampling
    Z <- list() # design matrix
    K <- list() # structure matrix
    K_rk <- list() # rank of structure matrix
    kappaList <- list() # list for kappa-parameters "kappa ~ Ga(kappa_a, kappa_b)"
    gammaList <- list() 
    solverList <- list() # list solving method, i.e. "rue" or "lanczos"
    kappa_startList <- list() # start values for kappa
    constraintList <- list() # if linear constraint is required
     
    # check correct possible situations
    check_intercept <- list("TRUE", "FALSE")
    if (!is.element(intercept, check_intercept))
        stop('Intercept can only be "TRUE" or "FALSE" ')
    
    check_family <- list("gaussian", "binomial", "poisson")
    if (!is.element(family, check_family))
        stop('"family" can only be "gaussian", "poisson" or "binomial"')
    
    check_link <- list("identity", "log", "logit")
    if (!is.element(link, check_link))
        stop('"link" can only be "identity" (gaussian), "log" (poisson) or "logit" (binomial) ')
    
    if ( (family != "gaussian") && (link == "identity") )
        stop('Gaussian family only possible if link = "identity" ')
    
    
    # calculate Z_k, K_k, K_rk and kappa/sigma from sx-function
    for (i in 2:length(mf)) {
        # distinguish between WITH intercept and WITHOUT intercept
        if (intercept == "TRUE") {
            # assume that first column is intercept
            kappa_startList[[1]] <- as.numeric(1)
            solverList[[1]] <- "rue"
            constraintList[[1]] <- "FALSE"
            gammaList[[1]] <- c(stats::rnorm(n = 1, sd = 0.1))
            kappaList[[1]] <- c(1, 0.00001)
            K[[1]] <- as(as.matrix(1), "dgCMatrix")
            K_rk[[1]] <- Matrix::rankMatrix(K[[1]])[1] 
            Z[[1]] <- as(as.matrix(rep(1, length(y)), nrow = length(y)), "dgCMatrix")
            
            # differ between vector and matrix as input feature
            if (!is.vector(mf[[i]])) {
                kappa_startList[[i]] <- attr(mf[[i]], "ka_start")
                kappaList[[i]] <- c(attr(mf[[i]], "ka_a"), attr(mf[[i]], "ka_b"))
                solverList[[i]] <- attr(mf[[i]], "solver")
                constraintList[[i]] <- attr(mf[[i]], "linear_const")
                
                # distinguish between given Z
                if (is.null(attr(mf[[i]], "Z"))) {
                    if (is.null(attr(mf[[i]], "gamma"))) {
                        gammaList[[i]] <- c(stats::rnorm(ncol(mf[[i]]), sd = 0.1))
                    } else {
                        gammaList[[i]] <- c(attr(mf[[i]], "gamma"))
                    }
                    K[[i]] <- as(K_matrix(dimension = dim(mf[[i]])[2], penalty = attr(mf[[i]], "penalty"),
                                          data = attr(mf[[i]], "K")), "dgCMatrix")
                    K_rk[[i]] <- Matrix::rankMatrix(K[[i]], method = "qr")[1]
                    attributes(mf[[i]]) <- attributes(mf[[i]])["dim"]
                    Z[[i]] <- as(as.matrix(mf[[i]]), "dgCMatrix")
                }
                
                if (!is.null(attr(mf[[i]], "Z"))) {
                    if (is.null(attr(mf[[i]], "gamma"))) {
                        gammaList[[i]] <- c( stats::rnorm(ncol(attr(mf[[i]], "Z")), sd = 0.1) )
                    } else {
                        gammaList[[i]] <- c(attr(mf[[i]], "gamma"))
                    }
                    K[[i]] <- as(attr(mf[[i]], "K"), "dgCMatrix")
                    K_rk[[i]] <- Matrix::rankMatrix(K[[i]], method = "qr")[1]
                    Z[[i]] <- as(attr(mf[[i]], "Z"), "dgCMatrix")
                }
                
            }
            if (is.vector(mf[[i]])) {
                kappa_startList[[i]] <- attr(mf[[i]], "ka_start")
                solverList[[i]] <- "rue"
                gammaList[[i]] <- c(stats::rnorm(1, sd = 0.1))
                kappaList[[i]] <- c(1, 0.00001)
                K[[i]] <- as(as.matrix(1), "dgCMatrix")
                K_rk[[i]] <- Matrix::rankMatrix(K[[i]], method = "qr")[1]
                Z[[i]] <- as(as.matrix(mf[[i]]), "dgCMatrix")
            }
        } else {
            if (!is.vector(mf[[i]])) {
                kappa_startList[[i - 1]] <- attr(mf[[i]], "ka_start")
                kappaList[[i - 1]] <- c(attr(mf[[i]], "ka_a"), attr(mf[[i]], "ka_b"))
                solverList[[i - 1]] <- attr(mf[[i]], "solver")
                constraintList[[i - 1]] <- attr(mf[[i]], "linear_const")
                gammaList[[i - 1]] <- c(stats::rnorm(ncol(mf[[i]])))
                
                # distinguish between given Z
                if (is.null(attr(mf[[i]], "Z"))) {
                    if (is.null(attr(mf[[i]], "gamma"))) {
                        gammaList[[i - 1]] <- c( stats::rnorm( ncol(mf[[i]]), sd = 0.1 ) )
                    } else {
                        gammaList[[i - 1]] <- c(attr(mf[[i]], "gamma"))
                    }
                    K[[i - 1]] <- as(K_matrix(dimension = dim(mf[[i]])[2], penalty = attr(mf[[i]], "penalty"),
                                              data = attr(mf[[i]], "K")), "dgCMatrix")
                    K_rk[[i - 1]] <- Matrix::rankMatrix(K[[i - 1]], method = "qr")[1]
                    
                    attributes(mf[[i]]) <- attributes(mf[[i]])["dim"]
                    Z[[i - 1]] <- as(as.matrix(mf[[i]]), "dgCMatrix")
                }
                
                if (!is.null(attr(mf[[i]], "Z"))) {
                    if (is.null(attr(mf[[i]], "gamma"))) {
                        gammaList[[i - 1]] <- c( stats::rnorm( ncol(attr( mf[[i]], "Z" )), sd = 0.1 ) )
                    } else {
                        gammaList[[i - 1]] <- c(attr(mf[[i]], "gamma"))
                    }
                    K[[i - 1]] <- as(attr(mf[[i]], "K"), "dgCMatrix")
                    K_rk[[i - 1]] <- Matrix::rankMatrix(K[[i - 1]], method = "qr")[1]
                    Z[[i - 1]] <- as(attr(mf[[i]], "Z"), "dgCMatrix")
                }
                
            }
            if (is.vector(mf[[i]])) {
                kappa_startList[[i - 1]] <- as.numeric(1)
                kappaList[[i - 1]] <- c(1, 0.00001)
                solverList[[i - 1]] <- "rue"
                constraintList[[i - 1]] <- "FALSE"
                gammaList[[i - 1]] <- c(stats::rnorm(1, sd = 0.1))
                K[[i - 1]] <- as(as.matrix(1), "dgCMatrix")
                K_rk[[i - 1]] <- Matrix::rankMatrix(K[[i - 1]], method = "qr")
                Z[[i - 1]] <- as(as.matrix(mf[[i]]), "dgCMatrix")
            }
        }
        
        
    } 
    
    sigmavalues <- c(sigma_a, sigma_b)
    sigma <- as.numeric(sigma)
    
    # Using Gibbs or Gibbs-with-MH depending on family
    if (family == "gaussian") {
        out <- Sarim::sarim_gibbs(y = y, Z = Z, K = K, K_rank = K_rk, gamma = gammaList,
                                  ka_start = kappa_startList, ka_values = kappaList,
                                  solver = solverList, lin_constraint = constraintList,
                                  sigma = sigma, sigma_values = sigmavalues,
                                  nIter = nIter,
                                  m = m, thr = thr)
        list_out <- list("coef_results" = out$coef_results,
                         "kappa_results" = out$kappa_results,
                         "sigma_results" = out$sigma_results,
                         "lanzcos_iterations" = out$lanzcos_iterations)
    }
    
    
    if (family != "gaussian") {
      # out <- Sarim::sarim_mcmc(y=y, Z = Z, K = K, K_rank = K_rk, gamma = gammaList,
      #                          ka_start = kappa_startList, ka_values = kappaList,
      #                          solver = solverList, lin_constraint = constraintList,
      #                          family = family, link = link,
      #                          nIter = nIter, Ntrials = Ntrials,
      #                          m = m, thr = thr)
      initialburnin=burnin
      burnin = FALSE
      
      gammaListList <- lapply(1:4,function(i,x)return(x),x=gammaList)
      kappaListList <- lapply(1:4,function(i,x)return(x),x=kappa_startList)
      kappa_mean<-gamma_mean<-list()
      gamma_mean<-lapply(1:length(gammaList),function(i)rep(0.0,length(gammaList[[i]])))
      kappa_mean<-lapply(1:length(kappa_startList),function(x)return(0))
  

      gamma_mean<-gamma_mean2<-gamma_mean3<-lapply(1:4,function(i,x)return(x),x=gamma_mean)
      kappa_mean<-kappa_mean2<-kappa_mean3<-lapply(1:4,function(i,x)return(x),x=kappa_mean)
      
      itercounter<-lapply(1:4,function(i)return(0))
      cat("\n")
      while(!burnin)                     
      {
       out <- parallel::mclapply(1:4,function(i,y, Z, K, K_rk, gammaListList,
       #out <- lapply(1:4,function(i,y, Z, K, K_rk, gammaListList,
                                kappaListList, ka_values,
                                solverList, lin_constraint,
                                family, link,
                                nIter, Ntrials,
                                m, thr, 
                                gamma_mean, gamma_mean2, gamma_mean3,
                                kappa_mean, kappa_mean2, kappa_mean3,
                                itercounter){
                                Sarim::sarim_mcmc_nosamples(y=y, Z = Z, K = K, K_rank = K_rk, gamma = gammaListList[[i]],
                                          ka_start = kappaListList[[i]], ka_values = kappaList,
                                          solver = solverList, lin_constraint = constraintList,
                                          family = family, link = link,
                                          nIter = nIter, burnin = initialburnin, Ntrials = Ntrials,
                                          m = m, thr = thr, 
                                          gammamean = gamma_mean[[i]], 
                                          gamma2mean = gamma_mean2[[i]], 
                                          gamma3mean = gamma_mean3[[i]],
                                          kappamean = kappa_mean[[i]], 
                                          kappa2mean = kappa_mean2[[i]], 
                                          kappa3mean = kappa_mean3[[i]],
                                          iterationcounter = itercounter[[i]]
                                )},
                                y, Z, K, K_rk, gammaListList, kappaListList, kappaList,
                                solverList, constraintList,
                                family, link, nIter, Ntrials,
                                m, thr, 
                                gamma_mean, gamma_mean2, gamma_mean3,
                                kappa_mean, kappa_mean2, kappa_mean3,
                                itercounter 
                     ,mc.cores = 4, mc.set.seed = FALSE
                     )
     # print(out)
      p<-psrf(out)
      burnin <- (p<1.05)
      cat(paste0("Computing burnin: ",out[[1]]$iterationcounter[[1]], " iterations done. R = ",round(p,2),"\n"))
      
      for (i in 1:length(out))
      {
        if(is.null(out[[i]]))
        {
          print(paste("Chain",i,"failed :("))
          gammaListList[[i]]<-gammaList  
          kappaListList[[i]]<-kappaList
          gamma_mean[[i]]<-gamma_mean2[[i]]<-gamma_mean3[[i]]<-lapply(1:length(gammaList),function(i)rep(0.0,length(gammaList[[i]])))
          kappa_mean[[i]]<-kappa_mean2[[i]]<-kappa_mean3[[i]]<-lapply(1:length(kappa_startList),function(x)return(0))
          itercounter[[i]]<-list(0)
        }
        else
        {
        for (j in 1:length(out[[i]]$coef_results))
        {
          #last<-dim(out[[i]]$coef_results[[j]])[2]
          gammaListList[[i]][[j]]<-out[[i]]$coef_results[[j]][,1]
          kappaListList[[i]][[j]]<-out[[i]]$kappa_results[[j]][1]
        }
          gamma_mean[[i]]=out[[i]]$gamma_mean
          gamma_mean2[[i]]=out[[i]]$gamma_mean2
          gamma_mean3[[i]]=out[[i]]$gamma_mean3
          kappa_mean[[i]]=out[[i]]$kappa_mean
          kappa_mean2[[i]]=out[[i]]$kappa_mean2
          kappa_mean3[[i]]=out[[i]]$kappa_mean3
          itercounter[[i]]=out[[i]]$iterationcounter[[1]]
        }
        }

      initalburnin <- 0
      }
      
      cat("Burnin done, sampling...\n")
      itercounter<-lapply(1:4,function(i)return(0))
      kappa_mean<-gamma_mean<-list()
      gamma_mean<-lapply(1:length(gammaList),function(i)rep(0.0,length(gammaList[[i]])))
      kappa_mean<-lapply(1:length(kappa_startList),function(x)return(0))
      
      
      gamma_mean<-gamma_mean2<-gamma_mean3<-lapply(1:4,function(i,x)return(x),x=gamma_mean)
      kappa_mean<-kappa_mean2<-kappa_mean3<-lapply(1:4,function(i,x)return(x),x=kappa_mean)
      
      out <- parallel::mclapply(1:4,function(i,y, Z, K, K_rk, gammaListList,
                                             #out <- lapply(1:4,function(i,y, Z, K, K_rk, gammaListList,
                                             kappaListList, ka_values,
                                             solverList, lin_constraint,
                                             family, link,
                                             nIter, Ntrials,
                                             m, thr, 
                                             gamma_mean, gamma_mean2, gamma_mean3,
                                             kappa_mean, kappa_mean2, kappa_mean3,
                                             itercounter){
        Sarim::sarim_mcmc(y=y, Z = Z, K = K, K_rank = K_rk, gamma = gammaListList[[i]],
                          ka_start = kappaListList[[i]], ka_values = kappaList,
                          solver = solverList, lin_constraint = constraintList,
                          family = family, link = link,
                          nIter = 1000, burnin = 0, Ntrials = Ntrials,
                          m = m, thr = thr, 
                          gammamean = gamma_mean[[i]], 
                          gamma2mean = gamma_mean2[[i]], 
                          gamma3mean = gamma_mean3[[i]],
                          kappamean = kappa_mean[[i]], 
                          kappa2mean = kappa_mean2[[i]], 
                          kappa3mean = kappa_mean3[[i]],
                          iterationcounter = itercounter[[i]]
        )},
        y, Z, K, K_rk, gammaListList, kappaListList, kappaList,
        solverList, constraintList,
        family, link, nIter, Ntrials,
        m, thr, 
        gamma_mean, gamma_mean2, gamma_mean3,
        kappa_mean, kappa_mean2, kappa_mean3,
        itercounter, mc.cores = 4)
         list_out <- out               
        #list_out <- list("coef_results" = out$coef_results,
        #                 "kappa_results" = out$kappa_results,
        #                 "accept_rate" = out$accept_rate,
        #                 "lanzcos_iterations" = out$lanzcos_iterations)
    }
    cat("\n")
    
    return(list_out)
}