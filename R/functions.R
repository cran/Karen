#' @keywords internal
eval.constraints <- function(theta, constr.lst){
  val.constr <- sapply(constr.lst, function(cnstr){
    eval(parse(text = gsub(pattern = "\\",
                           replacement = "",
                           x = unlist(strsplit(cnstr, split = "=", fixed = TRUE))[2],
                           fixed = TRUE)))
  })
  fixed.params <- get.fixedParams(constr.lst)
  names(val.constr) <- fixed.params
  return(val.constr)
}

#' @keywords internal
get.fixedParams <- function(constr.lst){
  fixedParams <- sapply(constr.lst, function(cnstr){
    gsub(pattern = "\\'\\]",
         replacement = "",
         x = gsub(pattern = "theta\\[\\'",
                  replacement = "",
                  x = unlist(strsplit(cnstr, split = "=", fixed = TRUE))[1],
                  fixed = TRUE),
         fixed = TRUE)
  })
  return(fixedParams)
}

#' Get the clone-average of the first two-order smoothing moments from a fitted Kalman Reaction Network.
#'
#' This function returns the clone-average of the first two-order smoothing moments from a Kalman Reaction Network previously fitted on a clonal tracking dataset.
#' @param res.fit A list returned by get.fit() containing the information of a fitted Kalman Reaction Network.
#' @param X Stochastic process. A 3-dimensional array whose dimensions are the time, the cell type and the clone respectively.
#' @param cell.cols Color legend for the cell types. Defaults to NULL, in which case no color legend for the cell types is provided.
#' @return No return value.
#' @examples
#' rcts <- c("HSC->T", ## reactions
#'           "HSC->M",
#'           "T->0",
#'           "M->0")
#'
#' cnstr <- c("theta\\[\\'HSC->T\\'\\]=(theta\\[\\'T->0\\'\\])",
#'            "theta\\[\\'HSC->M\\'\\]=(theta\\[\\'M->0\\'\\])")
#' latsts <- "HSC" ## latent cell types
#'
#' ctps <- unique(setdiff(c(sapply(rcts, function(r){ ## all cell types
#'   as.vector(unlist(strsplit(r, split = "->", fixed = TRUE)))
#' }, simplify = "array")), c("0", "1")))
#'
#' Y0 <- Y_CT$WAS[,setdiff(ctps,"HSC"),]
#' topClones <- 2
#' Y0 <- Y0[,,names(head(sort(apply(Y0!=0, 3, sum), decreasing = TRUE), topClones)),drop=FALSE]
#'
#' ## cluster parameters:
#' cl <- parallel::makeCluster(2, type = "PSOCK")
#'
#' ## initial condition:
#' X0 <- rep(0, length(ctps))
#' names(X0) <- ctps
#' X0["HSC"] <- 1
#'
#' ## mean vector and covariance matrix of X0:
#' m_0 <- replicate(dim(Y0)[3], X0, simplify = "array")
#' colnames(m_0) <- dimnames(Y0)[[3]]
#' P_0 <- Matrix::Diagonal(length(ctps) * dim(Y0)[3], 10)
#' rownames(P_0) <- colnames(P_0) <- rep(dimnames(Y0)[[3]], each = length(ctps))
#'
#' ## fit Karen on data:
#' res.fit <- get.fit(rct.lst = rcts,
#'                    constr.lst = cnstr,
#'                    latSts.lst = latsts,
#'                    ct.lst = ctps,
#'                    Y = Y0,
#'                    m0 = m_0,
#'                    P0 = P_0,
#'                    cl = cl,
#'                    list(nLQR = 1,
#'                         lmm = 0, ## needs to be >=5 for real applications
#'                         pgtol = 0,
#'                         relErrfct = 1e-5,
#'                         tol = 1e-3,
#'                         maxit = 0, ## needs to be increased for real applications
#'                         maxitEM = 1, ## needs to be increased for real applications
#'                         trace = 1,
#'                         verbose = TRUE,
#'                         FORCEP = FALSE))
#' parallel::stopCluster(cl)
#' get.sMoments.avg(res.fit)
##' @export
get.sMoments.avg <- function(res.fit, X = NULL, cell.cols = NULL){
  V <- res.fit$V # net-effect matrix
  nProc <- length(res.fit$cloneChunks) # nrow(summary(cl)) # number of cores
  Y <- res.fit$Y # simulated measurements
  Y_NA <- Y
  Y_NA[Y_NA == 0] <- NA
  nCL <- dim(Y)[3] # number of clones

  if(!is.null(cell.cols)){
    cols <- cell.cols[rownames(V)]
  }else{
    cols <- palette.colors(nrow(V), palette = "Classic Tableau")
  }

  if(!is.null(X)){
    tps <- as.numeric(rownames(X))
    rownames(X) <- (tps - min(tps))/(max(tps) - min(tps))

    X_avg <- apply(X, c(1,2), mean, na.rm = T)
    X_avg[is.nan(X_avg)] <- NA
  }

  Y_avg <- apply(Y_NA, c(1,2), mean, na.rm = T)
  Y_avg[is.nan(Y_avg)] <- NA
  mean_smooth_avg <- abind::abind(sapply(1:nProc, function(cnk){
    sapply(1:length(res.fit$cloneChunks[[cnk]]), function(cl){
      m_xt <- array(data = NA, dim = dim(Y[,,1]) + c(1,0))
      rownames(m_xt) <- c("0", rownames(Y))
      m_xt[rownames(t(res.fit$bwd.res$m_xt_Yn[[cnk]][,cl,])),] <- t(res.fit$bwd.res$m_xt_Yn[[cnk]][,cl,])
      return(m_xt)
    }, simplify = "array")
  }, simplify = FALSE))

  mean_smooth_avg <- apply(mean_smooth_avg, c(1,2), mean, na.rm = T)
  sd_smooth_avg <- abind::abind(sapply(1:nProc, function(cnk){
    idx.clones <- matrix(data = 1:(length(res.fit$cloneChunks[[cnk]])*nrow(res.fit$V)), nrow = nrow(res.fit$V), ncol = length(res.fit$cloneChunks[[cnk]]))
    sapply(1:length(res.fit$cloneChunks[[cnk]]), function(cl){
      sd_xt <- array(data = NA, dim = dim(Y[,,1]) + c(1,0))
      rownames(sd_xt) <- c("0", rownames(Y))
      sd_xt[dimnames(res.fit$bwd.res$V_xt_Yn[[cnk]][res.fit$idx.clones[,cl],res.fit$idx.clones[,cl],])[[3]],] <- sqrt(t(sapply(1:(dim(res.fit$bwd.res$V_xt_Yn[[cnk]][idx.clones[,cl],idx.clones[,cl],])[3]),
                                                                                                                               FUN = function(t){
                                                                                                                                 diag(nearestPD(res.fit$bwd.res$V_xt_Yn[[cnk]][idx.clones[,cl],idx.clones[,cl],t]))})))
      return(sd_xt)
    }, simplify = "array")
  }, simplify = FALSE))
  sd_smooth_avg <- apply(sd_smooth_avg, c(1,2), mean, na.rm = TRUE)

  matplot(as.numeric(rownames(Y_NA)), Y_avg, lty = 1, pch = 20, add = F, col = alpha(cols, alpha = .8), cex = 7,
          cex.axis = 2, cex.lab = 2, xlab = "t", ylab = expression("Y"[t]), main = paste("clone average ", sep = ""), cex.main = 2,
          xlim = c(0, max(as.numeric(rownames(Y)))), ylim = range(c(Y_avg, mean_smooth_avg, mean_smooth_avg - 1.96*sd_smooth_avg, mean_smooth_avg + 1.96*sd_smooth_avg), na.rm = T))
  if(!is.null(X)){
    matplot(as.numeric(rownames(X_avg)), X_avg, add = T, pch = 1, cex = 6, lwd = 2, col = alpha(cols, alpha = .8))
  }
  matplot(as.numeric(rownames(mean_smooth_avg)), mean_smooth_avg, lwd = 3, lty = 1, type = 'l', add = TRUE, col = cols)
  matplot(as.numeric(rownames(mean_smooth_avg)), mean_smooth_avg - 1.96*sd_smooth_avg, lwd = 3, lty = 3, type = 'l', add = TRUE, col = cols)
  matplot(as.numeric(rownames(mean_smooth_avg)), mean_smooth_avg + 1.96*sd_smooth_avg, lwd = 3, lty = 3, type = 'l', add = TRUE, col = cols)
}

#' Get the first two-order smoothing moments from a fitted Kalman Reaction Network.
#'
#' This function returns the first two-order smoothing moments from a Kalman Reaction Network previously fitted on a clonal tracking dataset.
#' @param res.fit A list returned by get.fit() containing the information of a fitted Kalman Reaction Network.
#' @param X Stochastic process. A 3-dimensional array whose dimensions are the time, the cell type and the clone respectively.
#' @param cell.cols Color legend for the cell types. Defaults to NULL, in which case no color legend for the cell types is provided.
#' @return No return value.
#' @examples
#' rcts <- c("HSC->T", ## reactions
#'           "HSC->M",
#'           "T->0",
#'           "M->0")
#'
#' cnstr <- c("theta\\[\\'HSC->T\\'\\]=(theta\\[\\'T->0\\'\\])",
#'            "theta\\[\\'HSC->M\\'\\]=(theta\\[\\'M->0\\'\\])")
#' latsts <- "HSC" ## latent cell types
#'
#' ctps <- unique(setdiff(c(sapply(rcts, function(r){ ## all cell types
#'   as.vector(unlist(strsplit(r, split = "->", fixed = TRUE)))
#' }, simplify = "array")), c("0", "1")))
#'
#'
#'
#' Y0 <- Y_CT$WAS[,setdiff(ctps,"HSC"),]
#' topClones <- 2
#' Y0 <- Y0[,,names(head(sort(apply(Y0!=0, 3, sum), decreasing = TRUE), topClones)),drop=FALSE]
#'
#'
#' ## cluster parameters:
#' cl <- parallel::makeCluster(2, type = "PSOCK")
#'
#' ## initial condition:
#' X0 <- rep(0, length(ctps))
#' names(X0) <- ctps
#' X0["HSC"] <- 1
#'
#' ## mean vector and covariance matrix of X0:
#' m_0 <- replicate(dim(Y0)[3], X0, simplify = "array")
#' colnames(m_0) <- dimnames(Y0)[[3]]
#' P_0 <- Matrix::Diagonal(length(ctps) * dim(Y0)[3], 10)
#' rownames(P_0) <- colnames(P_0) <- rep(dimnames(Y0)[[3]], each = length(ctps))
#'
#' ## fit Karen on data:
#' res.fit <- get.fit(rct.lst = rcts,
#'                    constr.lst = cnstr,
#'                    latSts.lst = latsts,
#'                    ct.lst = ctps,
#'                    Y = Y0,
#'                    m0 = m_0,
#'                    P0 = P_0,
#'                    cl = cl,
#'                    list(nLQR = 1,
#'                         lmm = 0, ## needs to be >=5 for real applications
#'                         pgtol = 0,
#'                         relErrfct = 1e-5,
#'                         tol = 1e-3,
#'                         maxit = 0, ## needs to be increased for real applications
#'                         maxitEM = 1, ## needs to be increased for real applications
#'                         trace = 1,
#'                         verbose = TRUE,
#'                         FORCEP = FALSE))
#' parallel::stopCluster(cl)
#' oldpar <- par(no.readonly = TRUE)
#' par(mar = c(5,5,2,2), mfrow = c(1,3))
#' get.sMoments(res.fit)
#' par(oldpar)
##' @export
get.sMoments <- function(res.fit, X = NULL, cell.cols = NULL){
  V <- res.fit$V # net-effect matrix
  nProc <- length(res.fit$cloneChunks) # nrow(summary(cl)) # number of cores
  Y <- res.fit$Y # simulated measurements
  Y_NA <- Y
  Y_NA[Y_NA == 0] <- NA

  if(!is.null(cell.cols)){
    cols <- cell.cols[rownames(V)]
  }else{
    cols <- palette.colors(nrow(V), palette = "Classic Tableau")
  }

  if(!is.null(X)){
    tps <- as.numeric(rownames(X))
    rownames(X) <- (tps - min(tps))/(max(tps) - min(tps))
  }

  lapply(1:nProc, function(cnk){
    lapply(1:length(res.fit$cloneChunks[[cnk]]), function(cl){
      idx.clones <- matrix(data = 1:(length(res.fit$cloneChunks[[cnk]])*nrow(V)), nrow = nrow(V), ncol = length(res.fit$cloneChunks[[cnk]]))
      mean_smooth <- t(res.fit$bwd.res$m_xt_Yn[[cnk]][,cl,])
      sd_smooth <- sqrt(t(sapply(1:(dim(res.fit$bwd.res$V_xt_Yn[[cnk]][idx.clones[,cl],idx.clones[,cl],])[3]), FUN = function(t){diag(nearestPD(res.fit$bwd.res$V_xt_Yn[[cnk]][idx.clones[,cl],idx.clones[,cl],t]))})))
      rownames(sd_smooth) <- rownames(mean_smooth)

      matplot(as.numeric(rownames(Y_NA)), Y_NA[,,res.fit$cloneChunks[[cnk]][cl]], lty = 1, pch = 20, type = 'p', add = F, col = alpha(palette.colors(nrow(V), palette = "Classic Tableau"), alpha = .8), cex = 2,
              cex.axis = 2, cex.lab = 2, xlab = "t", ylab = expression("Y"[t]), main = paste("clone ", cl, sep = ""), cex.main = 2,
              xlim = c(0, max(as.numeric(rownames(Y)))), ylim = range(c(Y_NA[,,res.fit$cloneChunks[[cnk]][cl]], mean_smooth, mean_smooth - 1.96*sd_smooth, mean_smooth + 1.96*sd_smooth), na.rm = T))
      if(!is.null(X)){
        matplot(as.numeric(rownames(X)), X[,,res.fit$cloneChunks[[cnk]][cl]], add = T, pch = 1, cex = 1.5, lwd = 2, col = alpha(palette.colors(nrow(V), palette = "Classic Tableau"), alpha = .8))
      }
      matplot(as.numeric(rownames(mean_smooth)), mean_smooth, lwd = 2, lty = 1, type = 'l', add = TRUE, col = palette.colors(nrow(V), palette = "Classic Tableau"))
      matplot(as.numeric(rownames(mean_smooth)), mean_smooth - 1.96*sd_smooth, lwd = 2, lty = 3, type = 'l', add = TRUE, col = palette.colors(nrow(V), palette = "Classic Tableau"))
      matplot(as.numeric(rownames(mean_smooth)), mean_smooth + 1.96*sd_smooth, lwd = 2, lty = 3, type = 'l', add = TRUE, col = palette.colors(nrow(V), palette = "Classic Tableau"))
      # legend(x = "topright", legend = rownames(V), col = palette.colors(nrow(V), palette = "Classic Tableau"), pch = 20, lwd = 5, cex = 1.5)
    })
  })
  plot.new()
  legend(x = "center", legend = rownames(V), col = palette.colors(nrow(V), palette = "Classic Tableau"), pch = 20, lwd = 5, cex = 1.5)
}


#' Get the cell differentiation network from a fitted Kalman Reaction Network.
#'
#' This function returns the cell differentiation network from a Kalman Reaction Network previously fitted on a clonal tracking dataset.
#' @param res.fit A list returned by get.fit() containing the information of a fitted Kalman Reaction Network.
#' @param edges.lab (logical) Defaults to FALSE, in which case the labels (weights) will not be printed on the network edges.
#' @param AIC (logical) Defaults to FALSE, in which case the Akaike Information Criterion is not reported.
#' @param cell.cols Color legend for the cell types. Defaults to NULL, in which case no color legend for the cell types is provided.
#' @return No return value.
#' @examples
#' rcts <- c("HSC->T", ## reactions
#'           "HSC->M",
#'           "T->0",
#'           "M->0")
#'
#' cnstr <- c("theta\\[\\'HSC->T\\'\\]=(theta\\[\\'T->0\\'\\])",
#'            "theta\\[\\'HSC->M\\'\\]=(theta\\[\\'M->0\\'\\])")
#' latsts <- "HSC" ## latent cell types
#'
#' ctps <- unique(setdiff(c(sapply(rcts, function(r){ ## all cell types
#'   as.vector(unlist(strsplit(r, split = "->", fixed = TRUE)))
#' }, simplify = "array")), c("0", "1")))
#'
#'
#'
#' Y0 <- Y_CT$WAS[,setdiff(ctps,"HSC"),]
#' topClones <- 2
#' Y0 <- Y0[,,names(head(sort(apply(Y0!=0, 3, sum), decreasing = TRUE), topClones)),drop=FALSE]
#'
#'
#' ## cluster parameters:
#' cl <- parallel::makeCluster(2, type = "PSOCK")
#'
#' ## initial condition:
#' X0 <- rep(0, length(ctps))
#' names(X0) <- ctps
#' X0["HSC"] <- 1
#'
#' ## mean vector and covariance matrix of X0:
#' m_0 <- replicate(dim(Y0)[3], X0, simplify = "array")
#' colnames(m_0) <- dimnames(Y0)[[3]]
#' P_0 <- Matrix::Diagonal(length(ctps) * dim(Y0)[3], 10)
#' rownames(P_0) <- colnames(P_0) <- rep(dimnames(Y0)[[3]], each = length(ctps))
#'
#' ## fit Karen on data:
#' res.fit <- get.fit(rct.lst = rcts,
#'                    constr.lst = cnstr,
#'                    latSts.lst = latsts,
#'                    ct.lst = ctps,
#'                    Y = Y0,
#'                    m0 = m_0,
#'                    P0 = P_0,
#'                    cl = cl,
#'                    list(nLQR = 1,
#'                         lmm = 0, ## needs to be >=5 for real applications
#'                         pgtol = 0,
#'                         relErrfct = 1e-5,
#'                         tol = 1e-3,
#'                         maxit = 0, ## needs to be increased for real applications
#'                         maxitEM = 1, ## needs to be increased for real applications
#'                         trace = 1,
#'                         verbose = TRUE,
#'                         FORCEP = FALSE))
#' parallel::stopCluster(cl)
#'
#' get.cdn(res.fit)
##' @export
get.cdn <- function(res.fit, edges.lab = FALSE, AIC = FALSE, cell.cols = NULL){

  if(!is.null(cell.cols)){
    cols <- cell.cols[rownames(res.fit$V)]
  }else{
    cols <- palette.colors(nrow(res.fit$V), palette = "Classic Tableau")
    names(cols) <- rownames(res.fit$V)
  }

  mycircle <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v]
    }
    vertex.size  <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v]
    }
    vertex.frame.color <- params("vertex", "frame.color")
    if (length(vertex.frame.color) != 1 && !is.null(v)) {
      vertex.frame.color <- vertex.frame.color[v]
    }
    vertex.frame.width <- params("vertex", "frame.width")
    if (length(vertex.frame.width) != 1 && !is.null(v)) {
      vertex.frame.width <- vertex.frame.width[v]
    }

    mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
           vertex.size, vertex.frame.width,
           FUN=function(x, y, bg, fg, size, lwd) {
             symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                     circles=size, add=TRUE, inches=FALSE)
           })
  }

  add.vertex.shape("fcircle", clip=igraph.shape.noclip,
                   plot=mycircle, parameters=list(vertex.frame.color=1,
                                                  vertex.frame.width=1))

  phi.curr <- res.fit$fit$par
  phi.curr <- c(eval.constraints(theta = head(phi.curr,-2), constr.lst = res.fit$constr.lst), phi.curr)

  adjNet <- matrix(data = NA, nrow = nrow(res.fit$V), ncol = nrow(res.fit$V))
  rownames(adjNet) <- colnames(adjNet) <- rownames(res.fit$V)
  for (j in 1:length(head(phi.curr,-2))) {
    r <- as.vector(unlist(strsplit(names(phi.curr)[j], split = "->", fixed = TRUE)))[1]
    p <- as.vector(unlist(strsplit(names(phi.curr)[j], split = "->", fixed = TRUE)))[2]
    if(p != "0" & p != "1"){
      adjNet[p,r] <- phi.curr[j]
    }else{
      if(p == "1"){
        adjNet[r,r] <- phi.curr[j]
      }
    }
  }
  # for (p in setdiff(res.fit$latSts.lst, "HSC")) {
  #   adjNet[p,"HSC"] <- sum(adjNet[setdiff(rownames(res.fit$V), p),p], na.rm = T)
  # }
  adjNet[adjNet == 0] <- 1e-8

  adjNet <- t(t(adjNet)/colSums(adjNet, na.rm = T))
  adjNet[is.nan(adjNet)] <- 0
  adjNet[is.na(adjNet)] <- 0
  net <- graph_from_adjacency_matrix(t(adjNet), mode = "directed", weighted = TRUE)

  E(net)$width <- 15
  V(net)[names(cols)]$color <- cols

  edge.start <- ends(net, es=E(net), names=F)[,1]
  edge.loops <- which(apply(ends(net, es=E(net), names=F), 1, function(e){e[1]==e[2]}))
  V(net)$label.cex <- c(1.5, rep(2.1, nrow(res.fit$V) - 1))
  arrow.width <- rep(2.5, length(E(net)))
  arrow.width[edge.loops] <- 1
  E(net)[edge.loops]$width <- 7

  col.bins <- cut(sort(c(unique(c(adjNet[!is.na(adjNet)])),1)), breaks = seq(0, 1, len = 100),
                  include.lowest = TRUE)
  col.pal <- colorRampPalette(c("lightgray", "red", "black"))(99)[col.bins]
  E(net)[order(E(net)$weight, decreasing = FALSE)]$color <- tail(head(col.pal, -1), -1)

  latSts.lst <- setdiff(res.fit$latSts.lst, "HSC")

  if(length(latSts.lst) > 0){
    if(length(latSts.lst) == 1){
      intProgCoords <- 0
    }else{
      intProgCoords <- seq(-length(latSts.lst), length(latSts.lst), length.out = length(latSts.lst))
      if(length(latSts.lst) == 3){
        intProgCoords <- c(intProgCoords[1], intProgCoords[3], intProgCoords[2])
      }
    }
    coords <- rbind(c(0,.3), # HSC
                    cbind(intProgCoords*1.2, .25), # INTERMEDIATE PROGENITORS
                    cbind(seq(-length(setdiff(rownames(res.fit$V), c("HSC", latSts.lst))),
                              length(setdiff(rownames(res.fit$V), c("HSC", latSts.lst))),
                              length.out = length(setdiff(rownames(res.fit$V), c("HSC", latSts.lst)))), .2))
  }else{
    coords <- rbind(c(0,.3), # HSC
                    cbind(seq(-length(setdiff(rownames(res.fit$V), c("HSC", latSts.lst))),
                              length(setdiff(rownames(res.fit$V), c("HSC", latSts.lst))),
                              length.out = length(setdiff(rownames(res.fit$V), c("HSC", latSts.lst)))), .2))
  }
  coords <- tail(coords, length(rownames(res.fit$V)))

  plot.network(net,
               edge.color = E(net)$color,
               vertex.size = c(35, rep(35, nrow(res.fit$V) - 1)),
               edge.label = ifelse(rep(edges.lab, length(E(net))), round(E(net)$weight,2), ""),
               edge.label.cex = 2,
               edge.label.font = 2,
               layout = coords,
               vertex.label.family = "Helvetica",
               edge.label.color = "black",
               vertex.frame.color = "black",
               vertex.frame.width = 4,
               vertex.frame.cex = 2,
               vertex.label.font = 2,
               vertex.label.color = "white",
               # edge.curved = .05,
               vertex.shape="fcircle",
               edge.arrow.size = arrow.width
  )
  if(AIC){
    mtext(side = 3, text = paste("AIC = ", round(res.fit$AIC,3), sep = ""), cex = 2, font = 2)
  }
}

#' @keywords internal
check.Y <- function(Y){
  if(!is.array(Y)){
    stop("Y must be a 3-dimensional array or a matrix.\n")
  }else{
    if(length(dim(Y)) < 2 | length(dim(Y)) > 3){
      stop("Y must be a 3-dimensional array or a matrix.\n")
    }else{
      if(length(dim(Y)) == 2){
        if(is.null(rownames(Y))){
          stop("Please specify row names in Y as the time-points.\n")
        }else{
          if(sum(is.na(as.numeric(rownames(Y)))) > 0){
            stop("Please specify row names in Y as the time-points.\n")
          }
        }
        if(is.null(colnames(Y))){
          stop("Please specify column names in Y as the cell types.\n")
        }
        Y_arr <- array(Y, dim = c(dim(Y),1))
        dimnames(Y_arr) <- list(rownames(Y),
                                colnames(Y),
                                "clone_1")
        Y <- Y_arr
        rm(Y_arr)
      }else{
        if(is.null(rownames(Y))){
          stop("Please specify row names in Y as the time-points.\n")
        }else{
          if(sum(is.na(as.numeric(rownames(Y)))) > 0){
            stop("Please specify row names in Y as the time-points.\n")
          }
        }
        if(is.null(colnames(Y))){
          stop("Please specify column names in Y as the cell types.\n")
        }
        if(is.null(dimnames(Y)[[3]])){
          stop("Please specify third-dimension names of Y as the unique barcode/clone identifiers.\n")
        }
      }
    }
  }
}

#' @keywords internal
check.rcts <- function(Y,
                       rct.lst,
                       latSts.lst,
                       ct.lst){
  pattern <- "[A-Z0-9]{1,}->[A-Z0-9]{0,}[0-1]{0,}"

  if(sum(grepl(pattern,
               rct.lst,
               ignore.case = FALSE)) < length(rct.lst)){
    stop(paste("Unvalid reactions provided in 'rct.lst': ",
              paste(rct.lst[which(!grepl(pattern,
                                         c(rct.lst),
                                         ignore.case = FALSE))], collapse = ", "), "\n", sep = ""))
  }

  if(sum(unlist(lapply(str_split(string = rct.lst, pattern = "->"), function(rct){rct[1]==rct[2]}))) > 0){
    stop(paste("Unvalid reactions provided in 'rct.lst': ",
              paste(rct.lst[which(unlist(lapply(str_split(string = rct.lst, pattern = "->"),
                                                function(rct){rct[1]==rct[2]})))], collapse = ", "), "\n", sep = ""))
  }

  ct.rct.lst <- setdiff(as.vector(unlist(c(sapply(rct.lst, function(r){
    as.vector(unlist(str_split(r, "->")))
  }, simplify = "array")))), c("0","1"))

  rgts.lst <- as.vector(unlist(lapply(str_split(string = rct.lst, pattern = "->"), function(s){s[1]})))
  rgts_prod.lst <- unique(c(ct.rct.lst, rgts.lst))

  if(sum(!(setdiff(rgts_prod.lst, latSts.lst) %in% colnames(Y))) > 0){
    stop(paste("Observed cell types in 'rct.lst' not present in 'Y': ",
              paste(as.vector(unlist(lapply(setdiff(rgts_prod.lst, latSts.lst)[which(!(setdiff(rgts_prod.lst, latSts.lst) %in% colnames(Y)))],
                                            function(l){paste("'", l, "'", sep = "")}))), collapse = ", "),
              "\nPlease provide compatible (observed) cell types in 'Y' and 'rct.lst'.\n",
              sep = ""))
  }

  if(sum(!(colnames(Y) %in% rgts_prod.lst)) > 0){
    stop(paste("Observed cell types in 'Y' not present in 'rct.lst': ",
              paste(as.vector(unlist(lapply(colnames(Y)[which(!(colnames(Y) %in% rgts_prod.lst))],
                                            function(l){paste("'", l, "'", sep = "")}))), collapse = ", "),
              "\nPlease provide compatible (observed) cell types in 'Y' and 'rct.lst'.\n",
              sep = ""))
  }

  if(sum(!(setdiff(rgts_prod.lst, colnames(Y)) %in% latSts.lst)) > 0){
    stop(paste("Latent cell types in 'rct.lst' not present in 'latSts.lst': ",
              paste(as.vector(unlist(lapply(setdiff(rgts_prod.lst, colnames(Y))[which(!(setdiff(rgts_prod.lst, colnames(Y)) %in% latSts.lst))],
                                            function(l){paste("'", l, "'", sep = "")}))), collapse = ", "),
              "\nPlease provide compatible (latent) cell types in 'rct.lst' and 'latSts.lst'.\n",
              sep = ""))
  }

  if(sum(!(latSts.lst %in% rgts_prod.lst)) > 0){
    stop(paste("Latent cell types in 'latSts.lst' not present in 'rct.lst': ",
              paste(as.vector(unlist(lapply(latSts.lst[which(!(latSts.lst %in% rgts_prod.lst))],
                                            function(l){paste("'", l, "'", sep = "")}))), collapse = ", "),
              "\nPlease provide compatible (latent) cell types in 'rct.lst' and 'latSts.lst'.\n",
              sep = ""))
  }

  if(sum(!(ct.lst %in% rgts_prod.lst)) > 0){
    stop(paste("Cell types in 'ct.lst' not present in 'rct.lst': ",
              paste(as.vector(unlist(lapply(ct.lst[which(!(ct.lst %in% rgts_prod.lst))],
                                            function(l){paste("'", l, "'", sep = "")}))), collapse = ", "),
              "\nPlease provide compatible cell types in 'rct.lst' and 'ct.lst'.\n",
              sep = ""))
  }

  if(sum(!(rgts_prod.lst %in% ct.lst)) > 0){
    stop(paste("Cell types in 'rct.lst' not present in 'ct.lst': ",
              paste(as.vector(unlist(lapply(rgts_prod.lst[which(!(rgts_prod.lst %in% ct.lst))],
                                            function(l){paste("'", l, "'", sep = "")}))), collapse = ", "),
              "\nPlease provide compatible cell types in 'rct.lst' and 'ct.lst'.\n",
              sep = ""))
  }
}

#' @keywords internal
check.cnstrs <- function(rct.lst,
                         constr.lst){
  pattern <- "^(theta\\\\\\[\\\\\\'[A-Z0-9]{1,}\\->[A-Z0-9]{0,}[0-1]{0,}\\\\\\'\\\\\\]){1,1}(\\=){1,1}[\\(]{1,1}(theta\\\\\\[\\\\\\'[A-Z0-9]{1,}\\->[A-Z0-9]{0,}[0-1]{0,}\\\\\\'\\\\\\]){1,1}(+[[:space:]]+\\++[[:space:]]+theta\\\\\\[\\\\\\'[A-Z0-9]{1,}\\->[A-Z0-9]{0,}[0-1]{0,}\\\\\\'\\\\\\]){0,}[\\)]{1,1}$"

  if(sum(grepl(pattern,
               constr.lst,
               ignore.case = FALSE)) < length(constr.lst)){
    stop(paste("Unvalid reactions provided in 'constr.lst':\n",
              paste(constr.lst[which(!grepl(pattern,
                                            c(constr.lst),
                                            ignore.case = FALSE))], collapse = ",\n"), "\n", sep = ""))
  }

  rcts.constr <- strsplits(constr.lst, splits = c("theta\\[\\'",
                                                  "\\'\\]=(theta\\[\\'",
                                                  "\\'\\]=(",
                                                  "\\'\\] + ",
                                                  "\\'\\])"), fixed = TRUE)

  if(sum(!(rcts.constr %in% rct.lst)) > 0){
    stop(paste("Unvalid reactions provided in 'constr.lst':\n",
              paste(as.vector(unlist(lapply(rcts.constr[which(!(rcts.constr %in% rct.lst))],
                                            function(l){paste("'", l, "'", sep = "")}))), collapse = ", "),
              "\nPlease provide compatible reactions in 'rct.lst' and 'constr.lst'.\n",
              sep = ""))
  }

  fixed.constr <- get.fixedParams(constr.lst)
  if(sum(duplicated(fixed.constr)) > 0){
    stop(paste("There is one or more parameter with multiple constraints in 'constr.lst':\n",
               paste(fixed.constr[duplicated(fixed.constr)], collapse = ", "),
               "\nPlease provide only one constraint per parameter.\n",
               sep = ""))
  }
}

#' @keywords internal
check.m0 <- function(m0, ct.lst){

  if(sum(is.na(rownames(m0))) > 0){
    stop(paste("NA values in rownames(m0) are not allowed!\n", sep = ""))
  }

  if(sum(!(ct.lst %in% rownames(m0)))){
    stop(paste("Cell types of 'ct.lst' not present in the initial condition m0:\n",
              paste(ct.lst[which(!(ct.lst %in% rownames(m0)))], collapse = ", "),
              sep = ""))
  }

  if(sum(!(rownames(m0) %in% ct.lst))){
    warning(paste("Cell types of the IC m0 not present in 'ct.lst':\n",
              paste(rownames(m0)[which(!(rownames(m0) %in% ct.lst))], collapse = ", "),
              "\nThese will be discarded.\n",
              sep = ""))
  }
}

#' @keywords internal
check.input <- function(Y,
                        rct.lst,
                        latSts.lst,
                        constr.lst,
                        m0,
                        ct.lst,
                        cl){

  check.cluster(Y, cl)
  check.Y(Y)
  check.rcts(Y,
             rct.lst,
             latSts.lst,
             ct.lst)
  if(length(constr.lst) > 0){
    check.cnstrs(rct.lst,
                 constr.lst)
  }
  check.m0(m0, ct.lst)
}

#' @keywords internal
check.cluster <- function(Y, cl){
  nCL <- dim(Y)[3]
  nProc <- nrow(summary(cl))
  if(nProc > nCL){
    stop("More processors than clones!\n Please provide a cluster with a number of processors lower than the number of clones.")
  }
}

#' @keywords internal
check.rcts.sim <- function(rct.lst,
                       latSts.lst,
                       ct.lst){

  pattern <- "[A-Z0-9]{1,}->[A-Z0-9]{0,}[0-1]{0,}"

  if(sum(grepl(pattern,
               rct.lst,
               ignore.case = FALSE)) < length(rct.lst)){
    stop(paste("Unvalid reactions provided in 'rct.lst': ",
              paste(rct.lst[which(!grepl(pattern,
                                         c(rct.lst),
                                         ignore.case = FALSE))], collapse = ", "), "\n", sep = ""))
  }

  if(sum(unlist(lapply(str_split(string = rct.lst, pattern = "->"), function(rct){rct[1]==rct[2]}))) > 0){
    stop(paste("Unvalid reactions provided in 'rct.lst': ",
              paste(rct.lst[which(unlist(lapply(str_split(string = rct.lst, pattern = "->"),
                                                function(rct){rct[1]==rct[2]})))], collapse = ", "), "\n", sep = ""))
  }

  ct.rct.lst <- setdiff(as.vector(unlist(c(sapply(rct.lst, function(r){
    as.vector(unlist(str_split(r, "->")))
  }, simplify = "array")))), c("0","1"))

  rgts.lst <- as.vector(unlist(lapply(str_split(string = rct.lst, pattern = "->"), function(s){s[1]})))
  rgts_prod.lst <- unique(c(ct.rct.lst, rgts.lst))

  if(sum(!(latSts.lst %in% rgts_prod.lst)) > 0){
    stop(paste("Latent cell types in 'latSts.lst' not present in 'rct.lst': ",
              paste(as.vector(unlist(lapply(latSts.lst[which(!(latSts.lst %in% rgts_prod.lst))],
                                            function(l){paste("'", l, "'", sep = "")}))), collapse = ", "),
              "\nPlease provide compatible (latent) cell types in 'rct.lst' and 'latSts.lst'.\n",
              sep = ""))
  }

  if(sum(!(ct.lst %in% rgts_prod.lst)) > 0){
    stop(paste("Cell types in 'ct.lst' not present in 'rct.lst': ",
              paste(as.vector(unlist(lapply(ct.lst[which(!(ct.lst %in% rgts_prod.lst))],
                                            function(l){paste("'", l, "'", sep = "")}))), collapse = ", "),
              "\nPlease provide compatible cell types in 'rct.lst' and 'ct.lst'.\n",
              sep = ""))
  }

  if(sum(!(rgts_prod.lst %in% ct.lst)) > 0){
    stop(paste("Cell types in 'rct.lst' not present in 'ct.lst': ",
              paste(as.vector(unlist(lapply(rgts_prod.lst[which(!(rgts_prod.lst %in% ct.lst))],
                                            function(l){paste("'", l, "'", sep = "")}))), collapse = ", "),
              "\nPlease provide compatible cell types in 'rct.lst' and 'ct.lst'.\n",
              sep = ""))
  }
}

#' @keywords internal
check.X0 <- function(X0, ct.lst){

  if(sum(is.na(names(X0))) > 0){
    stop(paste("NA values in names(X0) are not allowed!\n", sep = ""))
  }

  if(sum(!(ct.lst %in% names(X0)))){
    stop(paste("Cell types of 'ct.lst' not present in the initial condition X0:\n",
              paste(ct.lst[which(!(ct.lst %in% names(X0)))], collapse = ", "),
              sep = ""))
  }

  if(sum(!(names(X0) %in% ct.lst))){
    warning(paste("Cell types of the IC X0 not present in 'ct.lst':\n",
              paste(names(X0)[which(!(names(X0) %in% ct.lst))], collapse = ", "),
              "\nThese will be discarded.\n",
              sep = ""))
  }
}

#' @keywords internal
check.input.sim <- function(rct.lst,
                        latSts.lst,
                        constr.lst,
                        ct.lst,
                        X0){
  check.rcts.sim(rct.lst,
                 latSts.lst,
                 ct.lst)

  if(length(constr.lst) > 0){
    check.cnstrs(rct.lst,
                 constr.lst)
  }

  check.X0(X0, ct.lst)
}

#' @keywords internal
strsplits <- function(x, splits, ...)
{
  for (split in splits)
  {
    x <- unlist(strsplit(x, split, ...))
  }
  return(x[!x == ""]) # Remove empty values
}

#' Fit the state-space model to a clonal tracking dataset
#'
#' This function fits a state-space model to a clonal tracking dataset using an extended Kalman filter approach.
#' @param rct.lst A list of biochemical reactions defining the cell differentiation network.
#' A differentiation move from cell type "A" to cell type "B" must be coded as "A->B"
#' Duplication of cell "A" must be coded as "A->1"
#' Death of cell "A" must be coded as "A->0".
#' @param constr.lst (defaults to NULL, when no constraints are needed) List of linear constraints that must be applied to the biochemical reactions.
#' For example, if we need the constraint "A->B = B->C + B->D", this must be coded using the following syntax
#' c("theta\\[\\'A->B\\'\\]=(theta\\[\\'B->C\\'\\] + theta\\[\\'B->D\\'\\])").
#' @param latSts.lst List of the latent cell types. If for example counts are not available for cell types "A" and "B", then latSts.lst = c("A", "B").
#' @param ct.lst List of all the cell types involved in the network formulation.
#' For example, if the network is defined by the biochemical reactions are A->B" and "A->C", then ct.lst = c("A", "B", "C").
#' @param Y A 3-dimensional array whose dimensions are the time, the cell type and the clone respectively.
#' @param m0 mean vector of the initial condition \eqn{x_0}{x0}
#' @param P0 covariance matrix of the initial condition \eqn{x_0}{x0}
#' @param cl An object of class "cluster" specifying the cluster to be used for parallel execution. See makeCluster for more information.
#' If the argument is not specified, the default cluster is used. See setDefaultCluster for information on how to set up a default cluster.
#' @param control A a list of control parameters for the optimization routine:
#' \itemize{
#'  \item{"nLQR"}{(defaults to 3) is an integer giving the order of the Gauss-Legendre approximation for integrals.}
#'  \item{"lmm"}{(defaults to 25) is an integer giving the number of BFGS updates retained in the "L-BFGS-B" method.}
#'  \item{"pgtol"}{(defaults to 0 when check is suppressed) is a tolerance on the projected gradient
#'  in the current search direction of the "L-BFGS-B" method.}
#'  \item{"relErrfct"}{(defaults to 1e-5) is the relative error on the function value for the "L-BFGS-B" optimization.
#'  That is, the parameter "factr" of the optim() function is set to relErrfct/.Machine$double.eps.}
#'  \item{"tol"}{(defaults to 1e-9) is the relative error tolerance for the expectation-maximization algorithm
#'  of the extended Kalman filter optimization. That is, the optimization is run until the relative error of the function
#'  and of the parameter vector are lower than tol.}
#'  \item{"maxit"}{(defaults to 1000) The maximum number of iterations for the "L-BFGS-B" optimization.}
#'  \item{"maxitEM"}{(defaults to 10) The maximum number of iterations for the expectation-maximization algorithm.}
#'  \item{"trace"}{(defaults to 1) Non-negative integer. If positive, tracing information on the progress of the optimization is produced.
#' This parameter is also passed to the optim() function.
#' Higher values may produce more tracing information: for method "L-BFGS-B" there are six levels of tracing.
#' (To understand exactly what these do see the source code: higher levels give more detail.)}
#'  \item{"verbose"}{(defaults to TRUE) Logical value. If TRUE, then information messages on the progress of the
#'  filtering/smoothing algorithm are printed to the console.}
#'  \item{"FORCEP"}{(defaults to TRUE) Logical value. If TRUE, then all the covariance matrices involved in the algorithm
#'  are forced to be positive-definite and it helps the convergence of the optimization.}
#' }
#' @return A list containing the following:
#' \itemize{
#'  \item{"fit"}{The output list returned by the optim() function (See documenttion of optim() for more details).}
#'  \item{"bwd.res"}{First two-order moments of the estimated smoothing distribution.}
#'  \item{"m0.res"}{Mean vector of the smoothing distribution at time t = 0.}
#'  \item{"P0.res"}{Covariance matrix of the smoothing distribution at time t = 0.}
#'  \item{"AIC"}{Akaike Information Criterion (AIC) of the fitted model.}
#'  \item{"cloneChunks"}{List containing the chunks of clones that have been defined for parallel-computing.}
#'  \item{"V"}{The net-effect matrix associated to the differentiation network.}
#'  \item{"Y"}{The complete clonal tracking dataset that includes also the missing cell types.}
#'  \item{"rct.lst"}{The list of biochemical reactions.}
#'  \item{"constr.lst"}{The linear constraints applied on the reactions.}
#'  \item{"latSts.lst"}{The missing/latent cell types.}
#' }
#' @examples
#' rcts <- c("HSC->T", ## reactions
#'           "HSC->M",
#'           "T->0",
#'           "M->0")
#'
#' cnstr <- c("theta\\[\\'HSC->T\\'\\]=(theta\\[\\'T->0\\'\\])",
#'            "theta\\[\\'HSC->M\\'\\]=(theta\\[\\'M->0\\'\\])")
#' latsts <- "HSC" ## latent cell types
#'
#' ctps <- unique(setdiff(c(sapply(rcts, function(r){ ## all cell types
#'   as.vector(unlist(strsplit(r, split = "->", fixed = TRUE)))
#' }, simplify = "array")), c("0", "1")))
#'
#'
#'
#' Y0 <- Y_CT$WAS[,setdiff(ctps,"HSC"),]
#' topClones <- 2
#' Y0 <- Y0[,,names(head(sort(apply(Y0!=0, 3, sum), decreasing = TRUE), topClones)),drop=FALSE]
#'
#'
#' ## cluster parameters:
#' cl <- parallel::makeCluster(2, type = "PSOCK")
#'
#' ## initial condition:
#' X0 <- rep(0, length(ctps))
#' names(X0) <- ctps
#' X0["HSC"] <- 1
#'
#' ## mean vector and covariance matrix of X0:
#' m_0 <- replicate(dim(Y0)[3], X0, simplify = "array")
#' colnames(m_0) <- dimnames(Y0)[[3]]
#' P_0 <- Matrix::Diagonal(length(ctps) * dim(Y0)[3], 10)
#' rownames(P_0) <- colnames(P_0) <- rep(dimnames(Y0)[[3]], each = length(ctps))
#'
#' ## fit Karen on data:
#' res.fit <- get.fit(rct.lst = rcts,
#'                    constr.lst = cnstr,
#'                    latSts.lst = latsts,
#'                    ct.lst = ctps,
#'                    Y = Y0,
#'                    m0 = m_0,
#'                    P0 = P_0,
#'                    cl = cl,
#'                    list(nLQR = 1,
#'                         lmm = 0, ## needs to be >=5 for real applications
#'                         pgtol = 0,
#'                         relErrfct = 1e-5,
#'                         tol = 1e-3,
#'                         maxit = 0, ## needs to be increased for real applications
#'                         maxitEM = 1, ## needs to be increased for real applications
#'                         trace = 1,
#'                         verbose = TRUE,
#'                         FORCEP = FALSE))
##' @export
get.fit <- function(rct.lst,
                    constr.lst = NULL,
                    latSts.lst,
                    ct.lst,
                    Y,
                    m0,
                    P0,
                    cl = getDefaultCluster(),
                    control = list(nLQR = 3,
                                   lmm = 25,
                                   pgtol = 0,
                                   relErrfct = 1e-5, # 1e-9
                                   tol = 1e-9, # 1e-4
                                   maxit = 1000,
                                   maxitEM = 10,
                                   trace = 1,
                                   verbose = TRUE,
                                   FORCEP = TRUE)){

  check.input(Y, rct.lst, latSts.lst, constr.lst, m0, ct.lst, cl)

  control_params <- c("nLQR",
                      "lmm",
                      "pgtol",
                      "relErrfct",
                      "tol",
                      "maxit",
                      "maxitEM",
                      "trace",
                      "verbose",
                      "FORCEP")

  if(sum(!(names(control) %in% control_params))){
    control <- control[-which(!(names(control) %in% control_params))]
  }

  if(is.null(control$nLQR) |
     is.null(control$lmm) |
     is.null(control$pgtol) |
     is.null(control$relErrfct) |
     is.null(control$tol) |
     is.null(control$maxit) |
     is.null(control$maxitEM) |
     is.null(control$trace) |
     is.null(control$verbose) |
     is.null(control$FORCEP)){

    control_def = list(nLQR = 3,
                       lmm = 25,
                       pgtol = 0,
                       relErrfct = 1e-5, # 1e-9
                       tol = 1e-9, # 1e-4
                       maxit = 1000,
                       maxitEM = 10,
                       trace = 1,
                       verbose = TRUE,
                       FORCEP = TRUE)

    NULL_params <- control_params[which(!(control_params %in% names(control)))]
    warning(paste(paste(NULL_params, collapse = ", "), " not provided. Default values will be used.\n", sep = ""))

    control[NULL_params] <- control_def[NULL_params]
    control <- control[control_params]
    rm(control_def)
    rm(NULL_params)
  }
  rm(control_params)

  generate.h(rct.lst = rct.lst, constr.lst = constr.lst, envir = environment())
  generate.H(rct.lst = rct.lst, constr.lst = constr.lst, envir = environment())
  generate.get.Vth(ct.lst, rct.lst, constr.lst, envir = environment())

  V <- get.V(ct.lst = ct.lst, rct.lst = rct.lst)
  m0 <- m0[rownames(V),,drop=FALSE]

  if(control$verbose){cat(paste("n. of clones: ", dim(Y)[3], "\n", sep = ""))}

  tps <- c(0, as.numeric(rownames(Y)))
  tps <- (tps - min(tps))/(max(tps) - min(tps))
  rownames(Y) <- tail(tps,-1)
  rm(tps)

  Y.wl <- array(data = 0,
                dim = c(nrow(Y), ncol(Y) + length(latSts.lst), dim(Y)[3]),
                dimnames = list(rownames(Y),
                                rownames(V),
                                dimnames(Y)[[3]]))
  Y.wl[,colnames(Y),] <- Y
  Y <- Y.wl
  rm(Y.wl)

  Y <- Y[which(apply(Y != 0, 1, sum) != 0),,,drop=FALSE]
  Y <- Y[,rownames(V),,drop=FALSE]
  Y_NA <- Y
  Y_NA[Y_NA == 0] <- NA

  nCL <- dim(Y)[3]

  m0.curr <- m0
  P0.curr <- P0
  rm(m0)
  rm(P0)

  if(is.null(colnames(m0.curr))){
    colnames(m0.curr) <- dimnames(Y)[[3]]
  }

  if(is.null(rownames(P0.curr)) | is.null(colnames(P0.curr))){
    rownames(P0.curr) <- colnames(P0.curr) <- rep(dimnames(Y)[[3]], each = nrow(V))
  }

  if(length(constr.lst) > 0){
    fixed.params <- get.fixedParams(constr.lst)
    phi.curr <- c(rep(1, ncol(V) - length(fixed.params)), 1, 1)
    names(phi.curr) <- c(setdiff(colnames(V), fixed.params), "r0", "r1")
    rm(fixed.params)
  }else{
    phi.curr <- c(rep(1, ncol(V)), 1, 1)
    names(phi.curr) <- c(colnames(V), "r0", "r1")
  }

  Vjs <- bdiag(lapply(1:(length(phi.curr) - 2), function(j){V}))

  nLQR <- control$nLQR
  if(nLQR == 1){
    LQD.rule <- legendre.quadrature.rules(n = 2)[[1]]
  }else{
    LQD.rule <- tail(legendre.quadrature.rules(n = nLQR), 1)[[1]]
    LQD.rule <- LQD.rule[nrow(LQD.rule):1,]
    rownames(LQD.rule) <- 1:nrow(LQD.rule)
    if((nLQR %% 2) != 0){
      LQD.rule[ceiling(nLQR/2),1] <- 0
    }
  }

  relErrfct <- control$relErrfct
  eps <- .Machine$double.eps
  FACTR <- (relErrfct/eps)
  tol <- control$tol
  FORCEP <- control$FORCEP

  ################################
  nProc <- nrow(summary(cl))
  if(nCL %% nProc == 0){
    chunkSizes <- rep(nCL/nProc, nProc)
  }else{
    chunkSizes <- c(rep(x = floor(nCL/nProc), times = nProc - 1), floor(nCL/nProc) +  (nCL %% nProc))
  }
  cloneChunks <- split(dimnames(Y)[[3]], rep(1:nProc, chunkSizes))

  clusterEvalQ(cl, library("Matrix"))
  clusterEvalQ(cl, library("igraph"))
  clusterExport(cl, varlist = as.list(c("get.symmetric",
                                        "tr",
                                        "ldet",
                                        "get.mkt",
                                        "get.nl",
                                        "get.gnl",
                                        "get.forward.l",
                                        "get.forward",
                                        "get.dPkt_ds2",
                                        "get.dPkt",
                                        "ubdiag",
                                        "ubdiag2",
                                        "bdiag_m",
                                        "get.dmkt",
                                        "get.Pkt",
                                        "LQR.mat",
                                        "expm",
                                        "solve",
                                        "get.V",
                                        "get.h",
                                        "get.H",
                                        "get.Vth",
                                        "get.truncated",
                                        "nearestPD",
                                        "cloneChunks",
                                        "Y",
                                        "phi.curr",
                                        "ct.lst",
                                        "rct.lst",
                                        "m0.curr",
                                        "P0.curr",
                                        "V",
                                        "FORCEP",
                                        "LQD.rule",
                                        "get.backward")), envir = environment())
  ################################
  f.old <- Inf
  f.new <- get.nl.p(phi.curr, Y, V, m0.curr, P0.curr, cl = cl, cloneChunks = cloneChunks, forceP = FORCEP, LQD.rule = LQD.rule)#/nCL
  nSmoothingSteps <- 0
  totTime <- 0
  relErr <- relErr_f <- 1e3
  maxIT <- control$maxitEM
  curr.opt <- NULL

  while (nSmoothingSteps < maxIT & (relErr > tol | relErr_f > tol)) {
    curr.old <- curr.opt
    if(control$verbose){cat(paste("\tnSmooth = ", nSmoothingSteps,
              ";\terror = ", relErr,
              ";\tr0 = ", tail(phi.curr,2)[1],
              ";\tr1 = ", tail(phi.curr,2)[2],
              "\n",
              sep = ""))}
    phi.old <- phi.curr

    m0.curr_old <- m0.curr
    P0.curr_old <- P0.curr

    if(control$verbose){cat(paste("Complete log-likelihood optimization...\n", sep = ""))}
    time.s <- proc.time()
    curr.opt <- try(optim(par = phi.curr,
                          fn = get.nl.p,
                          gr = get.gnl.p,
                          Y = Y, V = V, m0 = m0.curr, P0 = P0.curr, cl = cl, cloneChunks = cloneChunks, forceP = FORCEP, LQD.rule = LQD.rule,
                          method = "L-BFGS-B",
                          lower = c(rep(0, length(phi.curr) - 2), 1e-3, 0),
                          upper = c(rep(+Inf, length(phi.curr) - 2), +Inf, +Inf),
                          control = list(lmm = control$lmm,
                                         # fnscale = FNSCALE,
                                         factr = FACTR,
                                         trace = control$trace,
                                         pgtol = control$pgtol,
                                         maxit = control$maxit)), silent = TRUE)

    time.e <- proc.time()
    diff.MLE <- time.e - time.s
    totTime <- totTime + diff.MLE[3]
    if(control$verbose){cat("\tDONE\n")}

    if(inherits(curr.opt,'try-error')){
      curr.opt <- curr.old
      break
    }

    f.new <- curr.opt$value
    if(f.old - f.new < 0 | (nSmoothingSteps > 0 & curr.opt$convergence %in% c(10, 51, 52))){
      break
    }
    f.old <- f.new

    phi.curr <- phi.0 <- curr.opt$par

    relErr_f <- as.numeric(abs(f.new - f.old)/(f.new))
    relErr <- as.numeric(sqrt(t(phi.curr - phi.old) %*% (phi.curr - phi.old))/sqrt(t(phi.old)%*%phi.old))

    if(control$verbose){cat("\tSmoothing step...")}

    clusterExport(cl, varlist = as.list("phi.curr"), envir = environment())

    bwd.res <- get.backward.p(Y, V, phi.curr, m0.curr, P0.curr, cloneChunks = cloneChunks, cl = cl, forceP = FORCEP, LQD.rule = LQD.rule)
    m0.curr <- Reduce(cbind, lapply(bwd.res$m_xt_Yn, function(A){A[,,1]}))
    if(nCL == 1){
      m0.curr <- matrix(m0.curr, nrow = nrow(V))
    }
    P0.curr <- bdiag(lapply(bwd.res$V_xt_Yn, function(A){A[,,1]}))
    if(control$verbose){cat("\tDONE\n")}

    if(!FORCEP){
      m0.curr <- get.truncated(m0.curr)
      P0.curr <- nearestPD(P0.curr)
    }

    rownames(m0.curr) <- rownames(V)
    colnames(m0.curr) <- dimnames(Y)[[3]]
    rownames(P0.curr) <- colnames(P0.curr) <- rep(dimnames(Y)[[3]], each = nrow(V))

    if(control$verbose){cat("\tExporting new initial conditions to the cluster...")}
    clusterExport(cl, varlist = as.list("m0.curr", "P0.curr"), envir = environment())
    if(control$verbose){cat("\tDONE\n")}

    gc()

    nSmoothingSteps <- nSmoothingSteps + 1
  }
  AIC.res <- get.AIC.p(phi.curr, Y, V, m0.curr, P0.curr, cl = cl, cloneChunks = cloneChunks, forceP = FORCEP, LQD.rule = LQD.rule)

  res <- list()
  res$fit <- curr.opt
  res$bwd.res <- bwd.res
  res$m0.res <- m0.curr
  res$P0.res <- P0.curr
  res$AIC <- AIC.res
  res$cloneChunks <- cloneChunks
  res$V <- V
  res$Y <- Y
  res$rct.lst <- rct.lst
  res$constr.lst <- constr.lst
  res$latSts.lst <- latSts.lst

  return(res)
}


#' @keywords internal
get.t <- function(x, a, b){
  ((b - a)*(x - min(x)))/(max(x) - min(x)) + a
}

#' @keywords internal
get.AIC <- function(phi, Y, V, m0, P0, forceP = FALSE){
  p <- length(phi)
  n <- nrow(Y)
  nl <- get.nl(phi, Y, V, m0, P0, forceP = forceP)
  AICc <- 2*p + 2*nl

  return(AICc)
}

#' @keywords internal
get.AIC.p <- function(phi, Y, V, m0, P0, cl, cloneChunks, forceP = FALSE, LQD.rule){
  p <- length(phi)
  n <- nrow(Y)
  nl <- get.nl.p(phi, Y, V, m0, P0, cl = cl, cloneChunks = cloneChunks, forceP = forceP, LQD.rule = LQD.rule)

  AICc <- 2*p + 2*nl

  return(AICc)
}

#' @keywords internal
sde.sim.Euler <- function(rct.lst, constr.lst, latSts.lst, V, t0=0, T=1, X0=1, N=100, delta = 1, th, s2, trunc = F){
  generate.h(rct.lst = rct.lst, constr.lst = constr.lst, envir = environment())
  generate.H(rct.lst = rct.lst, constr.lst = constr.lst, envir = environment())

  t <- c(t0,t0+cumsum(rep(delta,N)))
  Z <- matrix(data = NA, nrow = N + 1, ncol = nrow(V))
  X <- matrix(data = NA, nrow = N + 1, ncol = nrow(V))
  colnames(X) <- rownames(V)

  Dt <- (T-t0)/(N + 1)
  sDt <- sqrt(Dt)
  X[1,] <- X0

  for(i in 2:(N+1)){
    d1 <- V %*% eval(parse(text="get.h(x = X[i-1,], theta = th)"), envir = environment())
    s1 <- V %*% diag(eval(parse(text="get.h(x = X[i-1,], theta = th)"), envir = environment())) %*% t(V)
    diag(s1) <- diag(s1) + s2
    if(trunc){
      Z[i-1,] <- as.numeric(rtmvnorm(1, mean = rep(0, nrow(V)), sigma = s1,
                                     lower=rep(0, length = nrow(V))))
    }else{
      Z[i-1,] <- as.numeric(rmvnorm(1, mean = rep(0, nrow(V)), sigma = s1))
    }

    X[i,] <- X[i-1,] + d1*Dt + sDt*Z[i-1,]
    X[i,X[i,] < 0] <- 0
  }
  return(X)
}

#' Simulate a clonal tracking dataset from a given cell differentiation network.
#'
#' This function simulates clone-specific trajectories for a cell differentiation network associated to a set of (constrained) biochemical reactions,
#' cell types, and missing/latent cell types.
#' @param rct.lst A list of biochemical reactions defining the cell differentiation network.
#' A differentiation move from cell type "A" to cell type "B" must be coded as "A->B"
#' Duplication of cell "A" must be coded as "A->1"
#' Death of cell "A" must be coded as "A->0".
#' @param constr.lst (defaults to NULL, when no constraints are needed) List of linear constraints that must be applied to the biochemical reactions.
#' For example, if we need the constraint "A->B = B->C + B->D", this must be coded using the following syntax
#' c("theta\\[\\'A->B\\'\\]=(theta\\[\\'B->C\\'\\] + theta\\[\\'B->D\\'\\])").
#' @param latSts.lst List of the latent cell types. If for example counts are not available for cell types "A" and "B", then latSts.lst = c("A", "B").
#' @param ct.lst List of all the cell types involved in the network formulation.
#' For example, if the network is defined by the biochemical reactions are A->B" and "A->C", then ct.lst = c("A", "B", "C").
#' @param th The vector parameter that must be used for simulation. The length of th equals the number of unconstrained reactions plus 2
#' (for the noise parameters \eqn{(\rho_0, \rho_1)}{(r0, r1)}). Only positive parameters can be provided.
#' @param S The length of each trajectory.
#' @param nCL An integer defining the number of distinct clones.
#' @param X0 A p-dimensional vector for the initial condition of the cell types,
#' where \eqn{p}{p} is the number of distinct cell types provided in ct.lst.
#' @param s2 (defaults to 1e-8) A positive value for the overall noise variance.
#' @param r0 (defaults to 0) A positive value for the intercept defining the
#' noise covariance matrix \eqn{R_k = \rho_0 + \rho_1G_kX_k}{Rk = r0 + r1GkXk}).
#' @param r1 (defaults to 0) A positive value for the slope defining the
#' noise covariance matrix \eqn{R_k = \rho_0 + \rho_1G_kX_k}{Rk = r0 + r1GkXk}).
#' @param f (defaults to 0) The fraction of measurements that must be considered as missing/latent.
#' @param ntps Number of time points to consider from the whole simulated clonal tracking dataset.
#' @param trunc (defaults to FALSE) Logical, indicating whether sampling from a truncated multivariate normal must be performed.
#' @return A list containing the following:
#' \itemize{
#'  \item{"X"}{The simulated process.}
#'  \item{"Y"}{The simulated noisy-corrupted measurements.}
#' }
#' @examples
#' rcts <- c("HSC->T", ## reactions
#'           "HSC->M",
#'           "T->0",
#'           "M->0")
#'
#' cnstr <- NULL
#' latsts <- "HSC" ## latent cell types
#'
#' ctps <- unique(setdiff(c(sapply(rcts, function(r){ ## all cell types
#'   as.vector(unlist(strsplit(r, split = "->", fixed = TRUE)))
#' }, simplify = "array")), c("0", "1")))
#'
#' ## simulation parameters:
#' S <- 100 ## trajectories length
#' nCL <- 2 ## number of clones
#' X0 <- rep(0, length(ctps)) ## initial condition
#' names(X0) <- ctps
#' X0["HSC"] <- 1
#' ntps <- 5 ## number of time-points
#' f_NA <- 0 ## fraction of observed data
#'
#' th.true <- c(1.9538674, 1.0559815, 0.7232172, 0.7324133) ## dynamic parameters
#' names(th.true) <- rcts
#' s2.true <- 1e-8 ## additonal noise
#' r0.true <- .1 ## intercept noise parameter
#' r1.true <- .01 ## slope noise parameter
#'
#' ## simulate trajectories:
#' XY <- get.sim.trajectories(rct.lst = rcts,
#'                            constr.lst = cnstr,
#'                            latSts.lst = latsts,
#'                            ct.lst = ctps,
#'                            th = th.true,
#'                            S = S,
#'                            nCL = nCL,
#'                            X0 = X0,
#'                            s2 = s2.true,
#'                            r0 = r0.true,
#'                            r1 = r1.true,
#'                            f = f_NA,
#'                            ntps = ntps,
#'                            trunc = FALSE)
#'
#' XY$X ## process
#' XY$Y ## measurements
##' @export
get.sim.trajectories <- function(rct.lst,
                                 constr.lst = NULL,
                                 latSts.lst,
                                 ct.lst,
                                 th, S, nCL, X0, s2 = 1e-8, r0 = 0, r1 = 0, f = 0, ntps, trunc = FALSE){

  check.input.sim(rct.lst,
                  latSts.lst,
                  constr.lst,
                  ct.lst,
                  X0)

  V <- get.V(ct.lst = ct.lst, rct.lst = rct.lst)

  X <- array(data = NA, dim = c(S, nrow(V), nCL))
  dimnames(X)[[1]] <- 1:nrow(X)
  dimnames(X)[[2]] <- rownames(V)
  dimnames(X)[[3]] <- 1:nCL
  X[1,,] <- X0[rownames(V)]
  for (cl in 1:nCL) {
    X[,,cl] <- sde.sim.Euler(rct.lst = rct.lst,
                             constr.lst = constr.lst,
                             latSts.lst = latSts.lst,
                             V = V,
                             t0 = 0,
                             T = 1,
                             X0 = X0,
                             N = S - 1,
                             th = th,
                             s2 = s2,
                             trunc = trunc)
  }

  Y <- X
  Y <- Y + array(data = rnorm(n = prod(dim(Y)), mean = 0, sd = sqrt(s2 + r0 + r1*c(X))),
                 dim = dim(X))
  Y[Y < 0] <- 0

  tps <- round(seq(2,S, length.out = ntps))
  Y <- Y[tps,,]
  X <- X[tps,,]

  Y[,latSts.lst,] <- 0

  if(f > 0){
    Y[as.logical(rbinom(n = prod(dim(Y)), size = 1, prob = f))] <- 0
    Y <- Y[which(apply(Y != 0, 1, sum) != 0),,]
  }

  res <- list(X = X,
              Y = Y)

  return(res)
}

#' @keywords internal
get.symmetric <- function(X){
  return((X + t(X))/2)
}

#' @keywords internal
tr <- function(A){
  return(as.numeric(sum(diag(A))))
}

#' @keywords internal
ldet <- function(A, log = T){
  return(as.numeric(determinant(A, logarithm = log)$modulus))
}

#' @keywords internal
get.mkt <- function(t0, t, Vth, m0){
  mkt <- Matrix(expm::expm(Vth*(t - t0))) %*%m0
  return(mkt)
}

#' @keywords internal
get.nl <- function(phi, Y, V, m0, P0, forceP = FALSE, LQD.rule){
  th <- head(phi, -2)
  r0 <- tail(phi, 2)[1]
  r1 <- tail(phi, 2)[2]
  K <- nrow(Y)

  res.frwd <- get.forward.l(Y,V,phi, m0, P0, forceP, LQD.rule = LQD.rule)
  S <- bdiag(res.frwd$Sks)
  Si <- bdiag(lapply(seq(K), function(k){return(solve(res.frwd$Sks[[k]]))}))
  y <- as.vector(unlist(res.frwd$yks))
  mu <- as.vector(unlist(c(res.frwd$muks)))

  nl <- as.numeric(ldet(S) + t(y - mu) %*% Si %*% (y - mu))
  nl <- ifelse((is.infinite(nl) | is.nan(nl)), 1e50, nl)

  return(nl)
}

#' @keywords internal
get.nl.p <- function(phi, Y, V, m0, P0, cl, cloneChunks, forceP = FALSE, LQD.rule){
  nProc <- nrow(summary(cl))
  nl.p <- sum(parSapply(cl = cl, X = 1:nProc, FUN =  function(cnk){
    Y_chunk <- Y[,,cloneChunks[[cnk]],drop=FALSE]
    Y_chunk <- Y_chunk[names(which(apply(Y_chunk != 0, 1, sum) != 0)),,,drop = FALSE]
    return(  get.nl(phi = phi,
                    Y = Y_chunk,
                    V = V,
                    m0 = m0[,cloneChunks[[cnk]],drop=FALSE],
                    P0 =  P0[rownames(P0) %in% cloneChunks[[cnk]],colnames(P0) %in% cloneChunks[[cnk]]],
                    forceP = forceP,
                    LQD.rule = LQD.rule))
  }))

  return(nl.p)
}

#' @keywords internal
get.gnl <- function(phi, Y, V, m0, P0, forceP = FALSE, LQD.rule){
  th <- head(phi, -2)
  r0 <- tail(phi, 2)[1]
  r1 <- tail(phi, 2)[2]
  K <- nrow(Y)
  nCL <- dim(Y)[[3]]

  res.frwd <- get.forward(Y,V,phi, m0, P0, forceP = forceP, LQD.rule = LQD.rule)
  S <- bdiag(res.frwd$Sks)
  Si <- bdiag(lapply(seq(K), function(k){return(solve(res.frwd$Sks[[k]]))}))
  y <- as.vector(unlist(res.frwd$yks))
  mu <- as.vector(unlist(c(res.frwd$muks)))

  dSs <- res.frwd$dSks
  dmus <- res.frwd$dmuks

  dS_dr0s <- res.frwd$dSks_dr0
  dmu_dr0s <- res.frwd$dmuks_dr0
  dS_dr1s <- res.frwd$dSks_dr1
  dmu_dr1s <- res.frwd$dmuks_dr1

  p <- length(th)

  gnl_th <- sapply(1:length(th), function(j){
    dS_thj <- bdiag(lapply(1:K, function(k){

      sizes <- sapply(1:nCL, function(cl){sum(Y[k,,cl] > 0)})
      names(sizes) <- 1:nCL
      lags <- c(0, cumsum(head(sizes, -1)*p))
      sel.idx <- as.vector(unlist(lapply(which(sizes != 0), function(s){
        d <- sizes[s]
        return(((j-1)*d + 1):(j*d) + lags[s])
      })))

      dSkj <- dSs[[k]][sel.idx,sel.idx]
      return(dSkj)
    }))

    dmu_thj <- as.vector(unlist(lapply(1:K, function(k){
      sizes <- sapply(1:nCL, function(cl){sum(Y[k,,cl] > 0)})
      lags <- c(0, cumsum(head(sizes, -1)*p))
      sel.idx <- as.vector(unlist(lapply(which(sizes != 0), function(s){
        d <- sizes[s]
        return(((j-1)*d + 1):(j*d) + lags[s])
      })))

      return(dmus[[k]][sel.idx])
    })))
    gnl_thj <- as.numeric(tr(Si %*% dS_thj)  -2* t(dmu_thj) %*% (Si %*% (y - mu))   - t(y - mu) %*% Si %*% (dS_thj %*% Si %*% (y - mu)))
    return(gnl_thj)
  })

  dS_r0 <- bdiag(dS_dr0s)
  dmu_r0 <- as.vector(unlist(lapply(1:K, function(k){as.numeric(dmu_dr0s[[k]])})))
  gnl_r0 <- as.numeric(tr(Si %*% dS_r0)  -2* t(dmu_r0) %*% (Si %*% (y - mu))   - t(y - mu) %*% Si %*% (dS_r0 %*% Si %*% (y - mu)))

  dS_r1 <- bdiag(dS_dr1s)
  dmu_r1 <- as.vector(unlist(lapply(1:K, function(k){as.numeric(dmu_dr1s[[k]])})))
  gnl_r1 <- as.numeric(tr(Si %*% dS_r1)  -2* t(dmu_r1) %*% (Si %*% (y - mu))   - t(y - mu) %*% Si %*% (dS_r1 %*% Si %*% (y - mu)))

  return(c(gnl_th, gnl_r0, gnl_r1))
}

#' @keywords internal
get.gnl.p <- function(phi, Y, V, m0, P0, cl, cloneChunks, forceP = FALSE, LQD.rule){
  nProc <- nrow(summary(cl))
  gnl.p <- rowSums(parSapply(cl = cl, X = 1:nProc, FUN =  function(cnk){
    Y_chunk <- Y[,,cloneChunks[[cnk]], drop=FALSE]
    Y_chunk <- Y_chunk[names(which(apply(Y_chunk != 0, 1, sum) != 0)),,,drop = FALSE]
    return(  get.gnl(phi = phi,
                     Y = Y_chunk,
                     V = V,
                     m0 = m0[,cloneChunks[[cnk]], drop=FALSE],
                     P0 =  P0[rownames(P0) %in% cloneChunks[[cnk]],colnames(P0) %in% cloneChunks[[cnk]]],
                     forceP = forceP,
                     LQD.rule = LQD.rule))
  }))

  return(gnl.p)
}

#' @keywords internal
get.forward.l <- function(Y, V, phi, m0, P0, forceP, LQD.rule){

  th <- head(phi, -2)
  r0 <- tail(phi, 2)[1]
  r1 <- tail(phi, 2)[2]
  K <- nrow(Y)
  p <- length(th)
  nCL <- dim(Y)[[3]]
  tps <- as.numeric(rownames(Y))
  dT <- diff(c(0,tps))

  G <- diag(1, nrow(V))
  rownames(G) <- colnames(G) <- rownames(V)

  Vth <- eval(parse(text="get.Vth(th)"), envir = environment())
  Vthjs <- bdiag(lapply(1:length(th), function(j){Vth}))


  Sks <- muks <- vector("list", length = K)

  yks <- vector("list", length = K)
  m0s <- mks <- array(data = NA, dim = c(nrow(V), nCL, K))
  P0s <- Pks <- vector("list", length = K)

  for (k in 1:K) {
    obsSts.lst.k <- apply(which(!(Y[k,,,drop=FALSE] == 0), arr.ind = T), 1, function(r){paste(rownames(V)[r[2]], r[3], sep = "_")})

    ####### PREDICTION STEPS #######

    mk <- get.mkt(t0 = tps[k],
                  t = tps[k] + dT[k],
                  Vth = Vth,
                  m0 = m0)
    if(forceP){
      mk <- get.truncated(mk)
    }

    rownames(mk) <- rownames(V)

    Pk <- get.Pkt(t0 = tps[k],
                  t = tps[k] + dT[k],
                  th = th,
                  V = V,
                  Vth = Vth,
                  P0.tp = P0,
                  m0.tp = m0,
                  dt = dT[k],
                  forceP = forceP,
                  nCL = nCL,
                  LQD.rule = LQD.rule)

    ####### UPDATE STEPS #######

    G.bd <- bdiag(lapply(1:nCL, function(cl){G}))
    rownames(G.bd) <- colnames(G.bd) <- ave(rep(rownames(G), nCL) , rep(rownames(G), nCL) , FUN = function(i) paste0(i, '_', seq_along(i)))
    if(length(obsSts.lst.k) > 1){
      G.bd <- G.bd[obsSts.lst.k, ]
    }else{
      G.bd <- Matrix(G.bd[obsSts.lst.k, ], ncol = ncol(G.bd))
    }

    ####### update of ICs:
    muk <- G.bd %*% as.numeric(mk)
    Rk <- Diagonal(nrow(G.bd), r0) + r1*Diagonal(nrow(G.bd), as.numeric(muk))
    yk <- G.bd %*% as.numeric(Y[k,,])
    Sk <- G.bd %*% (Pk %*% t(G.bd))  + Rk
    Ski <- solve(Sk)
    Kk <- Pk%*% (t(G.bd) %*% Ski)
    m0 <- mk + matrix(data = as.numeric(Kk %*% as.numeric(yk - muk)), nrow(V), ncol = nCL)
    P0 <- Pk - Kk %*% (Sk %*% t(Kk))
    P0 <- get.symmetric(P0)
    if(forceP){
      m0 <- get.truncated(m0)
      P0 <- nearestPD(P0)
    }

    # store current results
    Sks[[k]] <- Sk
    muks[[k]] <- as.matrix(muk)
    yks[[k]] <- as.matrix(yk)
    m0s[,,k] <- as.matrix(m0)
    P0s[[k]] <- P0
    mks[,,k] <- as.matrix(mk)
    Pks[[k]] <- Pk
  }

  res <- list()
  res$Sks <- Sks
  res$muks <- muks
  res$yks <- yks
  res$m0s <- m0s
  res$P0s <- P0s
  res$mks <- mks
  res$Pks <- Pks

  return(res)
}

#' @keywords internal
get.forward <- function(Y, V, phi, m0, P0, forceP, LQD.rule){

  th <- head(phi, -2)
  r0 <- tail(phi, 2)[1]
  r1 <- tail(phi, 2)[2]
  K <- nrow(Y)
  p <- length(th)
  nCL <- dim(Y)[[3]]
  tps <- as.numeric(rownames(Y))
  dT <- diff(c(0,tps))

  G <- diag(1, nrow(V))
  rownames(G) <- colnames(G) <- rownames(V)

  Vth <- eval(parse(text="get.Vth(th)"), envir = environment())
  Vthjs <- bdiag(lapply(1:length(th), function(j){Vth}))

  Sks <- muks <- dSks <- dmuks <- dSks_dr0 <- dmuks_dr0 <- dSks_dr1 <- dmuks_dr1 <- vector("list", length = K)
  yks <- vector("list", length = K)

  m0s <- mks <- array(data = NA, dim = c(nrow(V), nCL, K))
  P0s <- Pks <- vector("list", length = K)

  dm0 <- Matrix(0, nrow(V)*length(th), ncol = nCL, sparse = T)
  rownames(dm0) <- rep(rownames(V), times = length(th))
  dP0 <- Diagonal(nrow(V)*length(th)*nCL, 0)

  #######################################################  NEW DERIVATIVES
  dm0r0 <- Matrix(0, nrow(V), ncol = nCL, sparse = T)
  rownames(dm0r0) <- rownames(V)
  dP0r0 <- Diagonal(nrow(V)*nCL, 0)
  rownames(dP0r0) <- colnames(dP0r0) <- rep(rownames(V), nCL)
  dm0r1 <- Matrix(0, nrow(V), ncol = nCL, sparse = T)
  rownames(dm0r1) <- rownames(V)
  dP0r1 <- Diagonal(nrow(V)*nCL, 0)
  rownames(dP0r1) <- colnames(dP0r1) <- rep(rownames(V), nCL)
  #######################################################

  as <- diag(1, length(th))
  rownames(as) <- names(th)
  dVth <- bdiag(lapply(1:ncol(as), function(j){
    eval(parse(text="get.Vth(as[,j])"), envir = environment())
  }))

  for (k in 1:K) {
    obsSts.lst.k <- apply(which(!(Y[k,,,drop=FALSE] == 0), arr.ind = T), 1, function(r){paste(rownames(V)[r[2]], r[3], sep = "_")})

    ####### PREDICTION STEPS #######
    mk <- get.mkt(t0 = tps[k],
                  t = tps[k] + dT[k],
                  Vth = Vth,
                  m0 = m0)
    if(forceP){mk <- get.truncated(mk)}
    rownames(mk) <- rownames(V)

    Pk <- get.Pkt(t0 = tps[k],
                  t = tps[k] + dT[k],
                  th = th,
                  V = V,
                  Vth = Vth,
                  P0.tp = P0,
                  m0.tp = m0,
                  dt = dT[k],
                  forceP = forceP,
                  nCL = nCL,
                  LQD.rule = LQD.rule)
    rownames(Pk) <- colnames(Pk) <- apply(cbind(rep(rownames(V), times = nCL), rep(1:nCL, each = nrow(V))),
                                          1, paste, collapse = "-")

    dmkjs <- get.dmkt(t0 = tps[k],
                      t = tps[k] + dT[k],
                      th = th,
                      Vthjs = Vthjs,
                      dVth = dVth,
                      dm0 = dm0,
                      m0.tp = m0,
                      dt = dT[k],
                      LQD.rule = LQD.rule)
    rownames(dmkjs) <- rep(rownames(V), length(th))

    dPks <- get.dPkt(t0 = tps[k],
                     t = tps[k] + dT[k],
                     th = th,
                     V = V,
                     Vth = Vth,
                     Vthjs = Vthjs,
                     dVth = dVth,
                     dm0 = dm0,
                     dP0 = dP0,
                     m0.tp = m0,
                     P0.tp = P0,
                     dt = dT[k],
                     forceP = forceP,
                     nCL = nCL,
                     LQD.rule = LQD.rule)
    rownames(dPks) <- colnames(dPks) <- rep(rownames(V), length(th)*nCL)

    dmkr0 <- get.mkt(t0 = tps[k],
                     t = tps[k] + dT[k],
                     Vth = Vth,
                     m0 = dm0r0)
    rownames(dmkr0) <- rownames(V)

    dPk_dr0 <- get.dPkt_ds2(t0 = tps[k],
                            t = tps[k] + dT[k],
                            th = th,
                            V = V,
                            Vth = Vth,
                            dm0s2 = dm0r0,
                            dP0s2 = dP0r0,
                            dt = dT[k],
                            nCL = nCL,
                            LQD.rule = LQD.rule)
    rownames(dPk_dr0) <- colnames(dPk_dr0) <- rownames(Pk)

    dmkr1 <- get.mkt(t0 = tps[k],
                     t = tps[k] + dT[k],
                     Vth = Vth,
                     m0 = dm0r1)
    rownames(dmkr1) <- rownames(V)

    dPk_dr1 <- get.dPkt_ds2(t0 = tps[k],
                            t = tps[k] + dT[k],
                            th = th,
                            V = V,
                            Vth = Vth,
                            dm0s2 = dm0r1,
                            dP0s2 = dP0r1,
                            dt = dT[k],
                            nCL = nCL,
                            LQD.rule = LQD.rule)
    rownames(dPk_dr1) <- colnames(dPk_dr1) <- rownames(Pk)

    ####### UPDATE STEPS #######
    G.bd <- bdiag(lapply(1:nCL, function(cl){G}))
    rownames(G.bd) <- colnames(G.bd) <- ave(rep(rownames(G), nCL) , rep(rownames(G), nCL) , FUN = function(i) paste0(i, '_', seq_along(i)))
    if(length(obsSts.lst.k) > 1){
      G.bd <- G.bd[obsSts.lst.k, ]
    }else{
      G.bd <- Matrix(G.bd[obsSts.lst.k, ], ncol = ncol(G.bd))
    }

    ####### update of ICs:
    muk <- G.bd %*% as.numeric(mk)
    Rk <- Diagonal(nrow(G.bd), r0) + r1*Diagonal(nrow(G.bd), as.numeric(muk))
    yk <- G.bd %*% as.numeric(Y[k,,])
    Sk <- G.bd %*% (Pk %*% t(G.bd))  + Rk
    Ski <- Matrix(solve(Sk), sparse = T)
    Kk <- Pk%*% (t(G.bd) %*% Ski)
    m0 <- mk + matrix(data = as.numeric(Kk %*% as.numeric(yk - muk)), nrow(V), ncol = nCL)
    P0 <- Pk - Kk %*% (Sk %*% t(Kk))
    P0 <- get.symmetric(P0)
    if(forceP){
      m0 <- get.truncated(m0)
      P0 <- nearestPD(P0)
    }

    Gjs.bd <- bdiag(lapply(1:nCL, function(cl){bdiag(replicate(length(th), matrix(G[names(which(Y[k,,cl] != 0)),], ncol = nrow(V)), simplify = F))}))
    ####### update of derivative's ICs w.r.t. theta:
    dmukjs <- Gjs.bd %*% as.numeric(dmkjs)
    dSkjs <- Gjs.bd %*% dPks %*% t(Gjs.bd) + r1 * Diagonal(nrow(Gjs.bd), as.numeric(dmukjs))

    Sk_blocks <- ubdiag(Sk, plot.graph = FALSE)
    Sk_js <- bdiag(lapply(Sk_blocks, function(b){bdiag_m(lapply(1:length(th), function(j){b}))}))
    Ski_blocks <- ubdiag(Ski, plot.graph = FALSE)
    Ski_js <- bdiag(lapply(Ski_blocks, function(b){bdiag_m(lapply(1:length(th), function(j){b}))}))
    Pk_blocks <- ubdiag2(Pk, plot.graph = FALSE, nCL = nCL)
    Pk_js <- bdiag_m(lapply(Pk_blocks, function(b){bdiag_m(lapply(1:length(th), function(j){b}))}))

    dKkjs <- dPks %*% (t(Gjs.bd) %*% Ski_js) - Pk_js %*% t(Gjs.bd) %*% (Ski_js %*% (dSkjs %*%  Ski_js))
    sizes <- sapply(1:nCL, function(cl){sum(Y[k,,cl] > 0)})
    yk_js <- as.vector(unlist(lapply(split(as.numeric(yk), rep(1:nCL, sizes)), function(s){rep(s, times = length(th))})))
    muk_js <- as.vector(unlist(lapply(split(as.numeric(muk), rep(1:nCL, sizes)), function(s){rep(s, times = length(th))})))

    Kk_js <- Pk_js %*% (t(Gjs.bd) %*% Ski_js)

    dm0 <- dmkjs + Matrix(data = as.numeric(dKkjs %*% as.numeric(yk_js - muk_js) - Kk_js %*% as.numeric(dmukjs)), nrow(V)*length(th), ncol = nCL)
    rownames(dm0) <- rep(rownames(V), times = length(th))

    dP0 <- get.symmetric(dPks - dKkjs %*% (Sk_js %*% t(Kk_js))
                         - Kk_js %*% (dSkjs %*% t(Kk_js))
                         - Kk_js %*% (Sk_js %*% t(dKkjs)))

    ####### update of derivative's ICs w.r.t r0:
    Id <- Diagonal(nrow(G.bd), 1)
    dmuk_dr0 <- G.bd %*% as.numeric(dmkr0)
    dSk_dr0 <- G.bd %*% dPk_dr0 %*% t(G.bd)  + r1*Diagonal(nrow(G.bd), as.numeric(dmuk_dr0)) + Id
    dKk_dr0 <- dPk_dr0 %*% (t(G.bd) %*% Ski) - Pk %*% t(G.bd) %*% (Ski %*% (dSk_dr0 %*% Ski))

    dm0r0 <- dmkr0 + matrix(data = as.numeric(dKk_dr0 %*% as.numeric(yk - muk) - Kk %*% as.numeric(dmuk_dr0)), nrow(V), ncol = nCL)
    rownames(dm0r0) <- rownames(V)
    dP0r0 <- dPk_dr0 - dKk_dr0 %*% (Sk %*% t(Kk)) - Kk %*% (dSk_dr0 %*% t(Kk)) - Kk %*% (Sk %*% t(dKk_dr0))
    dP0r0 <- get.symmetric(dP0r0)

    ####### update of derivative's ICs w.r.t r1:
    dmuk_dr1 <- G.bd %*% as.numeric(dmkr1)
    dSk_dr1 <- G.bd %*% dPk_dr1 %*% t(G.bd)  + Diagonal(nrow(G.bd), as.numeric(muk)) + r1*Diagonal(nrow(G.bd), as.numeric(dmuk_dr1))
    dKk_dr1 <- dPk_dr1 %*% (t(G.bd) %*% Ski) - Pk %*% t(G.bd) %*% (Ski %*% (dSk_dr1 %*% Ski))

    dm0r1 <- dmkr1 + matrix(data = as.numeric(dKk_dr1 %*% as.numeric(yk - muk) - Kk %*% as.numeric(dmuk_dr1)), nrow(V), ncol = nCL)
    rownames(dm0r1) <- rownames(V)
    dP0r1 <- dPk_dr1 - dKk_dr1 %*% (Sk %*% t(Kk)) - Kk %*% (dSk_dr1 %*% t(Kk)) - Kk %*% (Sk %*% t(dKk_dr1))
    dP0r1 <- get.symmetric(dP0r1)


    # store current results

    Sks[[k]] <- Sk
    muks[[k]] <- as.matrix(muk)
    dSks[[k]] <- dSkjs
    dmuks[[k]] <- dmukjs

    dSks_dr0[[k]] <- dSk_dr0
    dmuks_dr0[[k]] <- dmuk_dr0
    dSks_dr1[[k]] <- dSk_dr1
    dmuks_dr1[[k]] <- dmuk_dr1

    yks[[k]] <- as.matrix(yk)
    m0s[,,k] <- as.matrix(m0)
    P0s[[k]] <- P0
    mks[,,k] <- as.matrix(mk)
    Pks[[k]] <- Pk
  }

  res <- list()
  res$Sks <- Sks
  res$muks <- muks
  res$dSks <- dSks
  res$dmuks <- dmuks
  res$yks <- yks
  res$m0s <- m0s
  res$P0s <- P0s
  res$mks <- mks
  res$Pks <- Pks

  res$dSks_dr0 <- dSks_dr0
  res$dmuks_dr0 <- dmuks_dr0
  res$dSks_dr1 <- dSks_dr1
  res$dmuks_dr1 <- dmuks_dr1

  return(res)
}

#' @keywords internal
get.dPkt_ds2 <- function(t0, t, th, V, Vth, dm0s2, dP0s2, dt = 1, nCL, LQD.rule){
  X01 <- diag(1,nrow(V))
  rownames(X01) <- colnames(X01) <- rownames(dm0s2)

  dBK_mk <- lapply(1:nrow(V), function(i){
    dBK_mki <- V %*% (diag(eval(parse(text="get.h(x = X01[,i], theta = th)"), envir = environment())) %*% t(V))
    dBK_mki <- bdiag(replicate(nCL, dBK_mki, simplify = FALSE))
    return(dBK_mki)
  })

  fct <- function(x, t0, t, V, Vth, dt = dt, nCL) {

    dms2x <- get.mkt(t0 = t0,
                     t = x,
                     Vth = Vth,
                     m0 = dm0s2)
    rownames(dms2x) <- rownames(dm0s2)

    Qx <- dt*Reduce("+", lapply(1:nrow(V), function(i){dBK_mk[[i]] * rep(c(dms2x[i,]), each = nrow(V))}))


    EVth_t_x <- expm::expm((t - x)*Vth)
    EVth_t_x.bd <- bdiag(lapply(1:nCL, function(cl){EVth_t_x}))

    fx <- EVth_t_x.bd %*% (Qx %*% t(EVth_t_x.bd))
    return(fx)
  }


  I <- LQR.mat(fct = fct,
                lower = t0,
                upper = t,
                rule = LQD.rule,
                t0 = t0, t = t, V = V, Vth = Vth, dt = dt, nCL = nCL)

  I <- I + t(I)
  I <- I/2

  EVth_t_t0 <- expm::expm((t - t0)*Vth)
  EVth_t_t0.bd <- bdiag(lapply(1:nCL, function(cl){EVth_t_t0}))

  dPkt_ds2 <- EVth_t_t0.bd %*% (dP0s2 %*% t(EVth_t_t0.bd))  + I
  return(as.matrix(dPkt_ds2))
}

#' @keywords internal
get.dPkt <- function(t, t0, th, V, Vth, Vthjs, dVth, dm0, dP0, m0.tp, P0.tp, dt = 1, forceP, nCL, LQD.rule){
  p <- length(th)
  X01 <- diag(1,nrow(V))
  rownames(X01) <- colnames(X01) <- rownames(m0.tp)

  Vjs <- bdiag(lapply(1:length(th), function(j){V}))
  dBK_mk <- lapply(1:nrow(V), function(i){
    dBK_mki <- Vjs %*% (Diagonal(ncol(V)*length(th), rep(eval(parse(text="get.h(x = X01[,i], theta = th)"), envir = environment()), length(th))) %*% t(Vjs))
    dBK_mki <- bdiag(replicate(nCL, dBK_mki, simplify = FALSE))
    return(dBK_mki)
  })

  th01 <- diag(1,length(th))
  rownames(th01) <- colnames(th01) <- names(th)

  dVth.bd <- bdiag(lapply(1:nCL, function(cl){dVth}))

  idx.mat <- matrix(data = as.vector(unlist(split(1:(p*nrow(V)*nCL), ceiling(seq_along(1:(p*nrow(V)*nCL))/(p))))), nrow = p, ncol = nrow(V)*nCL, byrow = T)


  fct <- function(x, t0, t, th, V, Vth, Vthjs, m0.tp, P0.tp, dVth, dVth.bd, dt = dt, forceP = forceP, nCL = nCL, LQD.rule) {
    mx <- get.mkt(t0 = t0, t = x, Vth = Vth, m0 = m0.tp)
    if(forceP){mx <- get.truncated(mx)}
    rownames(mx) <- rownames(m0.tp)

    dBk_th <- bdiag(lapply(1:ncol(th01), function(j){
      VV <- bdiag(lapply(1:nCL, function(cl){V}))
      H <- Diagonal(ncol(V)*nCL, x = c(t(eval(parse(text="get.H(mx, th01[,j], nCL)"), envir = environment()))))
      return(VV %*% (H %*% t(VV)))
    }))

    dBk_th <- bdiag(sapply(1:nCL, function(cl){
      return(dBk_th[c(t(idx.mat[,((cl - 1)*nrow(V) + 1):(cl*nrow(V))])), c(t(idx.mat[,((cl - 1)*nrow(V) + 1):(cl*nrow(V))]))])
    }))

    Px <- get.Pkt(x, t0, th, V = V, Vth, P0.tp, m0.tp, dt = dt, forceP = forceP, nCL = nCL, LQD.rule = LQD.rule)
    Px_blocks <- ubdiag(Px, plot.graph = FALSE)
    Px <- bdiag(lapply(Px_blocks, function(b){bdiag(lapply(1:length(th), function(j){b}))}))

    dmx <- get.dmkt(x, t0, th, Vthjs, dVth, dm0, m0.tp, dt = dt, LQD.rule = LQD.rule)
    rownames(dmx) <- rownames(dm0)

    Qx <- (dVth.bd %*% Px + Px %*% t(dVth.bd)
           + dt*(Reduce("+", lapply(1:nrow(V),
                                    function(i){
                                      dmxi <- dmx[which(rownames(dmx) %in% rownames(V)[i]),,drop=FALSE]
                                      dmxi <- dmxi[rep(1:nrow(dmxi), each = nrow(V)),]
                                      return(dBK_mk[[i]] * as.numeric(dmxi))
                                    })) + dBk_th))

    EVthjs_t_x <- expm::expm((t - x)*Vthjs)
    EVthjs_t_x.bd <- bdiag(lapply(1:nCL, function(cl){EVthjs_t_x}))

    fx <- EVthjs_t_x.bd %*% (Qx %*% t(EVthjs_t_x.bd))
    return(fx)
  }


  I <- LQR.mat(fct = fct,
                lower = t0,
                upper = t,
                rule = LQD.rule,
                t0 = t0, t = t, th = th, V = V, Vth = Vth, Vthjs = Vthjs, m0.tp = m0.tp, P0.tp = P0.tp, dVth = dVth, dVth.bd, dt = dt, forceP = forceP, nCL = nCL, LQD.rule = LQD.rule)

  I <- I + t(I)
  I <- I/2

  EVthjs_t_t0 <- expm::expm((t - t0)*Vthjs)
  EVthjs_t_t0.bd <- bdiag(lapply(1:nCL, function(cl){EVthjs_t_t0}))


  dPkt <- EVthjs_t_t0.bd %*% (dP0 %*% t(EVthjs_t_t0.bd))  + I

  return(dPkt)
}

#' @keywords internal
ubdiag <- function(mat, plot.graph = FALSE) {
  stopifnot(nrow(mat) == ncol(mat))
  x <- mat
  diag(x) <- 1
  if(isDiagonal(x)){
    edges <- cbind(1:nrow(x), 1:nrow(x))
    colnames(edges) <- c("i", "j")
  }else{
    x <- Matrix(x, sparse = T)
    edges <- as.matrix(summary(x)[c("i", "j")])
  }
  g <- graph.edgelist(edges, directed = FALSE)
  groups <- unique(lapply(Map(sort, neighborhood(g, nrow(mat))), function(l){as.vector(unlist(l))}))
  sub.Mat <- Map(`[`, list(mat), groups, groups, drop = FALSE)
  sub.mat <- Map(as.matrix, sub.Mat)
  return(sub.mat)
}

#' @keywords internal
ubdiag2 <- function(A, plot.graph = FALSE, nCL = nCL){
  bl.idxs <- matrix(data = 1:nrow(A), nrow = nrow(A)/nCL, ncol = nCL)
  A.l <- apply(bl.idxs, 2,  function(l){
    A[l,l]
  }, simplify = FALSE)
  return(A.l)
}

#' @keywords internal
bdiag_m <- function(lmat) {
  k <- nrow(lmat[[1]])
  N <- length(lmat)
  if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- N * k
  ## result: an   M x M  matrix
  return(new("dgCMatrix", Dim = c(M,M),
             i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
             p = k * 0L:M,
             x = as.numeric(unlist(lapply(lmat, function(l){as.numeric(unlist(l))})))))
}

#' @keywords internal
get.dmkt <- function(t, t0, th, Vthjs, dVth, dm0, m0.tp, dt = 1, LQD.rule){

  m0.tpjs <- do.call(rbind, replicate(length(th), m0.tp, simplify = F))

  fct <- function(x, t0, t, Vthjs, m0.tpjs, dVth) {
    return(expm::expm(as.numeric(t - x)*Vthjs) %*% dVth %*% (expm::expm((x - t0)*Vthjs) %*% m0.tpjs))
  }

  I <- Matrix(LQR.mat(fct = fct,
                       lower = t0,
                       upper = t,
                       rule = LQD.rule,
                       t0 = t0, t = t, Vthjs = Vthjs, m0.tpjs = m0.tpjs, dVth = dVth))

  dmkt <- I + expm::expm((t - t0)*Vthjs) %*% dm0

  return(dmkt)
}

#' @keywords internal
get.Pkt <- function(t, t0, th, V, Vth, P0.tp, m0.tp, dt = 1, forceP, nCL, LQD.rule){

  fct <- function(x, Vth, dt, V, m0.tp, th, t, forceP, nCL) {
    mx <- expm::expm(Vth*(x - t0)) %*%m0.tp
    if(forceP){mx <- get.truncated(mx)}
    rownames(mx) <- rownames(m0.tp)

    VV <- bdiag(lapply(1:nCL, function(cl){V}))
    H <- Diagonal(ncol(V)*nCL, x = c(t(eval(parse(text="get.H(mx, th, nCL)"), envir = environment()))))
    EVth <- expm::expm(as.numeric(t - x)*Vth)
    EVth <- bdiag_m(lapply(1:nCL, function(cl){EVth}))
    M <- EVth %*% (dt * VV %*% H %*% t(VV)) %*% t(EVth)
    return(M)
  }


  I <- Matrix(LQR.mat(fct = fct,
                       lower = t0,
                       upper = t,
                       rule = LQD.rule,
                       Vth = Vth, dt = dt, V = V, m0.tp = m0.tp, th = th, t = t, forceP = forceP, nCL = nCL))

  I <- I + t(I)
  I <- I/2

  EVth <- expm::expm((t - t0)*Vth)
  EVth <- bdiag_m(lapply(1:nCL, function(cl){EVth}))

  Pkt <- EVth %*% (P0.tp %*% t(EVth)) + I

  if(forceP){
    return(nearestPD(Pkt))
  }else{
    return(Pkt)
  }
}

#' @keywords internal
LQR.mat <- function(fct, lower, upper, rule, ...){
  lambda <- (upper - lower)/(2)
  mu <- (lower + upper)/(2)
  y <- lambda * rule$x + mu
  w <- rule$w
  return(lambda * Reduce("+", lapply(1:length(y), function(s){w[s]*fct(y[s], ...)})))
}

#' @keywords internal
solve <- function(A){
  Ai <- try(chol2inv(chol(A)),silent = TRUE)
  if(inherits(Ai,'try-error')){
    Ai <- try(base::solve(A),silent = TRUE)
    if(inherits(Ai,'try-error')){
      Ai <- ginv(as.matrix(A))
    }
  }
  return(Ai)
}

#' @keywords internal
get.backward <- function(Y, V, phi, m0, P0, forceP = FALSE, LQD.rule){
  nCL <- dim(Y)[3]
  dT <- diff(c(0, as.numeric(rownames(Y))))

  frwd.res <- get.forward.l(Y = Y,
                            V = V,
                            phi = phi,
                            m0 = m0,
                            P0 = P0,
                            forceP = forceP,
                            LQD.rule = LQD.rule)

  K <- nrow(Y)
  m0s <- frwd.res$m0s
  P0s <- sapply(frwd.res$P0s, function(l){as.matrix(l)}, simplify = "array")
  mks <- frwd.res$mks
  Pks <- sapply(frwd.res$Pks, function(l){as.matrix(l)}, simplify = "array")

  m_xt_Yn <- array(data = NA, dim = c(nrow(mks), ncol(mks), dim(mks)[3] + 1))
  V_xt_Yn <- array(data = NA, dim = c(nrow(Pks), ncol(Pks), dim(Pks)[3] + 1))

  m_xt_Yn[,,K + 1] <- m0s[,,K]
  V_xt_Yn[,,K + 1] <- P0s[,,K]

  for (k in (K-1):0) {
    V.th <- eval(parse(text="get.Vth(head(phi,-1))"), envir = environment())
    eV.tht <- expm::expm(t(V.th) * dT[k + 1])
    eV.tht.bd <- bdiag(lapply(1:nCL, function(cl){eV.tht}))
    if(k > 0){
      Bk <- Matrix(P0s[,,k]) %*% eV.tht.bd %*% solve(Matrix(Pks[,,k + 1]))
      m_xt_Yn[,,k + 1] <- m0s[,,k] + matrix(data = Bk %*% c(m_xt_Yn[,,k+2] - mks[,,k + 1]), nrow = nrow(V), ncol = nCL)
      if(forceP){m_xt_Yn[,,k + 1] <- get.truncated(m_xt_Yn[,,k + 1])}
      V_xt_Yn[,,k + 1] <- as.matrix(P0s[,,k] + Bk %*% (Matrix(V_xt_Yn[,,k+2]) - Matrix(Pks[,,k + 1])) %*% t(Bk))
      if(forceP){V_xt_Yn[,,k + 1] <- as.matrix(nearestPD(Matrix(V_xt_Yn[,,k + 1])))}
      V_xt_Yn[,,k + 1] <- get.symmetric(V_xt_Yn[,,k + 1])
    }else{
      Bk <- P0 %*% eV.tht.bd %*% solve(Matrix(Pks[,,k + 1]))
      m_xt_Yn[,,k + 1] <- m0 + matrix(data = Bk %*% c(m_xt_Yn[,,k+2] - mks[,,k + 1]), nrow = nrow(V), ncol = nCL)
      if(forceP){m_xt_Yn[,,k + 1] <- get.truncated(m_xt_Yn[,,k + 1])}
      V_xt_Yn[,,k + 1] <- as.matrix(P0 + Bk %*% (Matrix(V_xt_Yn[,,k+2]) - Matrix(Pks[,,k + 1])) %*% t(Bk))
      if(forceP){V_xt_Yn[,,k + 1] <- as.matrix(nearestPD(Matrix(V_xt_Yn[,,k + 1])))}
      V_xt_Yn[,,k + 1] <- get.symmetric(V_xt_Yn[,,k + 1])
    }
  }

  res <- list()
  res$m_xt_Yn <- m_xt_Yn
  res$V_xt_Yn <- V_xt_Yn

  return(res)
}

#' @keywords internal
get.backward.p <- function(Y, V, phi, m0, P0, cloneChunks, cl, forceP = FALSE, LQD.rule){
  nProc <- nrow(summary(cl))
  bkw.p <- parLapply(cl = cl, X = 1:nProc, fun = function(cnk){
    Y_chunk <- Y[,,cloneChunks[[cnk]], drop = FALSE]
    Y_chunk <- Y_chunk[names(which(apply(Y_chunk != 0, 1, sum) != 0)),,,drop = FALSE]
    return(get.backward(Y = Y_chunk,
                        V = V,
                        phi = phi,
                        m0 = m0[,cloneChunks[[cnk]], drop = FALSE],
                        P0 =  P0[rownames(P0) %in% cloneChunks[[cnk]],colnames(P0) %in% cloneChunks[[cnk]]],
                        forceP = forceP,
                        LQD.rule = LQD.rule))
  })

  tps_chunks <- lapply(1:nProc, function(cnk){
    Y_chunk <- Y[,,cloneChunks[[cnk]], drop = FALSE]
    return(names(which(apply(Y_chunk != 0, 1, sum) != 0)))
  })

  res <- list()
  res$m_xt_Yn <- lapply(1:nProc, function(cnk){
    temp <- bkw.p[[cnk]]$m_xt_Yn
    dimnames(temp)[[3]] <- c(0, tps_chunks[[cnk]])
    return(temp)
  })

  res$V_xt_Yn <- lapply(1:nProc, function(cnk){
    temp <- bkw.p[[cnk]]$V_xt_Yn
    dimnames(temp)[[3]] <- c(0, tps_chunks[[cnk]])
    return(temp)
  })

  return(res)
}

#' @keywords internal
get.V <- function(ct.lst = ct.lst, rct.lst = rct.lst){
  V <- matrix(data = 0, nrow = length(ct.lst), ncol = length(rct.lst))
  rownames(V) <- ct.lst
  colnames(V) <- rct.lst

  for (r in rct.lst) {
    rgts_prod <- as.vector(unlist(strsplit(r, split = "->", fixed = TRUE)))
    V[rgts_prod[1],r] <- -1
    if(rgts_prod[2] != "0"){
      if(rgts_prod[2] == "1"){
        V[rgts_prod[1],r] <- 1
      }else{
        V[rgts_prod[2],r] <- 2
      }
    }
  }
  return(V)
}

#' @keywords internal
generate.h <- function(rct.lst, constr.lst, envir){
  get.h.string <- paste("get.h <- function(x, theta){
                h <- c(",
                        paste(sapply(rct.lst, function(r){
                          rgts <- as.vector(unlist(strsplit(r, split = "->", fixed = TRUE)))[1]
                          hx <- paste("x['", rgts, "']*theta['", r, "']", sep = "")
                        }, simplify = "array"), collapse = ",\n"),
                        ")
return(h)
}", sep = "")

  for (constr in constr.lst) {
    constr <- as.vector(unlist(strsplit(constr, split = "=", fixed = TRUE)))
    pattern <- constr[1]
    replacement <- constr[2]

    get.h.string <- gsub(pattern = pattern,
                         replacement = replacement,
                         x = get.h.string)
  }

  eval(parse(text=get.h.string), envir = envir) # envir = .GlobalEnv
}

#' @keywords internal
generate.H <- function(rct.lst, constr.lst, envir){
  get.H.string <- paste("get.H <- function(x, theta, nCL){
  h <- matrix(data = c(",
                        paste(sapply(rct.lst, function(r){
                          rgts <- as.vector(unlist(strsplit(r, split = "->", fixed = TRUE)))[1]
                          hx <- paste("x['", rgts, "',]*theta['", r, "']", sep = "")
                        }, simplify = "array"), collapse = ",\n"),
                        "), nrow = nCL, ncol = ncol(V))
return(h)
}", sep = "")

  for (constr in constr.lst) {
    constr <- as.vector(unlist(strsplit(constr, split = "=", fixed = TRUE)))
    pattern <- constr[1]
    replacement <- constr[2]

    get.H.string <- gsub(pattern = pattern,
                         replacement = replacement,
                         x = get.H.string)
  }

  eval(parse(text=get.H.string), envir = envir) # envir = .GlobalEnv
}

#' @keywords internal
generate.get.Vth <- function(ct.lst, rct.lst, constr.lst, envir){
  V <- get.V(ct.lst = ct.lst, rct.lst = rct.lst)
  Vth <- matrix(data = "0", nrow(V), nrow(V))
  rownames(Vth) <- colnames(Vth) <- ct.lst
  for (ct in ct.lst) {
    for (r in names(V[ct,V[ct,] != 0])) {
      ct.rgt <- as.vector(unlist(strsplit(r, split = "->", fixed = TRUE)))[1]
      V_ct_r <- eval(parse(text=paste("V['", ct, "','", r, "']", sep = "")))
      if(Vth[ct, ct.rgt] == "0"){
        Vth[ct, ct.rgt] <- paste("theta['", r, "']*(", V_ct_r, ")", sep = "")
      }else{
        Vth[ct, ct.rgt] <- paste(Vth[ct, ct.rgt], " + ", paste("theta['", r, "']*(", V_ct_r, ")", sep = ""), sep = "")
      }
    }
  }

  paste(as.vector(unlist(lapply(split(as.vector(t(Vth)), ceiling(seq_along(1:(nrow(V)^2))/nrow(V))), function(lst){
    paste(lst, collapse = ",")
  }))), collapse = ",\n")

  get.Vth.string <- paste("get.Vth <- function(theta){
  Vth <- matrix(data = c(",
                          paste(as.vector(unlist(lapply(split(as.vector(t(Vth)), ceiling(seq_along(1:(nrow(V)^2))/nrow(V))), function(lst){
                            paste(lst, collapse = ",")
                          }))), collapse = ",\n"),
                          "), byrow = T, nrow(V), nrow(V))
  return(Vth)
}", sep = "")

  for (constr in constr.lst) {
    constr <- as.vector(unlist(strsplit(constr, split = "=", fixed = TRUE)))
    pattern <- constr[1]
    replacement <- constr[2]

    get.Vth.string <- gsub(pattern = pattern,
                           replacement = replacement,
                           x = get.Vth.string)
  }

  eval(parse(text=get.Vth.string), envir = envir)
}

#' @keywords internal
get.truncated <- function(x){
  x.tr <- x
  x.tr[x.tr < 0] <- 0
  return(x.tr)
}


#' Nearest Positive Definite Matrix
#'
#' This function first check if a matrix A is positive definite, typically a correlation or variance-covariance matrix.
#' If A is not positive definite, this function computes the nearest positive definite matrix of A using the function nearPD from package Matrix.
#' @param A numeric \eqn{n \times n}{n by n} approximately positive definite matrix, typically an approximation to a correlation or covariance matrix.
#' If A is not symmetric (and ensureSymmetry is not false), symmpart(A) is used.
#' @param ... Further arguments to be passed to nearPD (see package Matrix for details).
#' @return The nearest positive definite matrix of A.
#' @examples
#' nearestPD(diag(c(1,0,1)))
##' @export
nearestPD <- function(A,...){
  if(inherits(try(chol2inv(chol(A)),silent = TRUE), 'try-error')){
    A.pd <- try(nearPD(A,...),silent = TRUE)
    if(!inherits(A.pd, 'try-error')){
      A <- A.pd$mat
    }
  }
  return(A)
}

#' @keywords internal
plot.network <- function(graph,
                      layout,
                      edge.arrow.size,
                      vertex.shape="circle",
                      edge.lty=rep(1, ecount(graph)),
                      edge.curved=autocurve.edges(graph), ...) {

  plot.igraph(graph,
              edge.arrow.size=0,
              edge.lty=0,
              layout=layout,
              vertex.shape="none")
  for (e in seq_len(ecount(graph))) {
    graph2 <- delete.edges(graph, E(graph)[(1:ecount(graph))[-e]])
    plot.igraph(graph2,
                edge.lty=edge.lty[e],
                edge.arrow.size=edge.arrow.size[e],
                edge.curved=edge.curved[e],
                layout=layout,
                vertex.shape="none",
                vertex.label=NA,
                add=TRUE)
  }
  plot.igraph(graph,
              edge.lty=0,
              edge.arrow.size=0,
              layout=layout,
              vertex.shape=vertex.shape,
              add=TRUE, ...)
  invisible(NULL)
}
