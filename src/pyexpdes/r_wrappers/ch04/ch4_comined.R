args <- commandArgs(trailingOnly = TRUE)
in_dir <- args[1]
out_dir <- args[2]

# This wrapper reproduces the "Candidate generation / CoMinED" section in ch4.R
# and exports candidate sets for Python plotting.
#
# Note: The original script calls `Lattice(...)` without explicitly attaching a
# package that defines it. To avoid fragile package dependencies, we define a
# small deterministic `Lattice()` helper here.

Lattice <- function(N, p){
  m <- ceiling(N^(1/p))
  grid <- lapply(1:p, function(i) (0:(m-1) + 0.5) / m)
  pts <- as.matrix(expand.grid(grid))
  if (nrow(pts) < N) stop("Lattice(): not enough points generated")
  return(pts[1:N, , drop=FALSE])
}

constraint <- function(x){
  c1 <- (x[1] - sqrt(50 * (x[2] - 0.52)^2 + 2) + 1)
  c2 <- (sqrt(120 * (x[2] - 0.48)^2 + 1) - 0.75 - x[1])
  c3 <- (0.65^2 - x[1]^2 - x[2]^2)
  return (c(c1,c2,c3))
}

# define constraint contour
x1 <- x2 <- matrix(NA, nrow = 3, ncol = 1001)
x2.seq <- seq(0,1,length.out = 1001)
x2[1,] <- x2.seq
x1[1,] <- sqrt(50 * (x2.seq - 0.52)^2 + 2) - 1
x2[2,] <- x2.seq
x1[2,] <- sqrt(120 * (x2.seq - 0.48)^2 + 1) - 0.75
x2[3,] <- x2.seq
x1[3,] <- sqrt(pmax(0, 0.65^2 - x2.seq^2))

contour <- list(x1 = x1, x2 = x2)

logdis <- function(x, s=2){
  if (s > 0) length(x)/s*log(sum(abs(x)^s)) else sum(log(abs(x)))
}

mined.seq <- function(n, cand, cand.lf, s = 2, return.obj = F){
  N <- nrow(cand)
  if (n > N) stop(message("Not enough candidate points!"))
  p <- ncol(cand)
  idx <- which.max(cand.lf)
  xi <- c(idx)
  val <- NULL
  for (i in 2:n){
    cand.dist <- cand - rep(1,N) %*% t(cand[idx,])
    val <- cbind(val, 0.5*cand.lf[idx] + 0.5*cand.lf + apply(cand.dist, 1, logdis, s=s))
    idx <- which.max(apply(val, 1, min))
    xi <- c(xi, idx)
  }
  if (return.obj){
    val <- val[xi,]
    val <- cbind(val, c(val[n,],-Inf))
    obj <- min(val[upper.tri(val)])
    return (list(idx = xi, obj = obj))
  }
  return (xi)
}

comined <- function(n, p, tau, constraint, n.aug, auto.scale = F,
                    s = 2){
  # initialize sample
  samp <- Lattice(n*n.aug, p)
  min.dist <- min(dist(samp))
  min.ddist <- Inf # dimensional minimum
  for (i in 1:p){
    ddist <- min(dist(samp[,i]))
    if (ddist < min.ddist) min.ddist <- ddist
  }
  samp.gval <- matrix(t(apply(samp,1,constraint)), nrow=nrow(samp))
  scale <- rep(1, ncol(samp.gval))
  if (auto.scale) scale <- apply(samp.gval,2,mad,center=0)
  samp.lf <- apply(samp.gval,1,function(x) sum(pnorm(-tau[2]*x/scale,log.p=T)))
  med.op <- mined.seq(n, samp, samp.lf, s = s, return.obj = TRUE)
  samp.med <- samp[med.op$idx,]

  for (k in 3:length(tau)){
    min.dist <- min.dist / 2
    min.ddist <- min.ddist / 2
    no.decimal <- attr(regexpr("(?<=\\.)0+", format(min.ddist, scientific = F),
                               perl = TRUE), "match.length") + 1
    samp.med.dist <- as.matrix(dist(samp.med))
    samp.aug <- NULL
    for (i in 1:n){
      nn.idx <- order(samp.med.dist[i,])[2:(n.aug+1)]
      samp.aug <- rbind(samp.aug, 0.5*(samp.med[nn.idx,]+rep(1,n.aug)%*%t(samp.med[i,])))
      samp.aug <- rbind(samp.aug, 0.5*(3*samp.med[nn.idx,]-rep(1,n.aug)%*%t(samp.med[i,])))
    }
    samp.aug.rep <- round(samp.aug, digits = no.decimal)
    samp.aug.dp <- duplicated(samp.aug.rep)
    samp.aug <- samp.aug[!samp.aug.dp,]
    samp.aug.out <- apply(samp.aug, 1, function(x) (any(x<0)|any(x>1)))
    samp.aug <- samp.aug[!samp.aug.out,]
    no.aug <- nrow(samp.aug)
    samp.rep <- rbind(samp.aug, samp)
    samp.rep <- round(samp.rep, digits = no.decimal)
    samp.rep.dp <- duplicated(samp.rep, fromLast = TRUE)
    samp.aug <- samp.aug[!(samp.rep.dp[1:no.aug]),]
    if (nrow(samp.aug) > 0){
      samp.aug.gval <- matrix(t(apply(samp.aug,1,constraint)),
                              nrow=nrow(samp.aug))
      samp <- rbind(samp, samp.aug)
      samp.gval <- rbind(samp.gval, samp.aug.gval)
    }
    if (auto.scale) scale <- apply(samp.gval,2,mad,center=0)
    samp.lf <- apply(samp.gval,1,function(x) sum(pnorm(-tau[k]*x/scale,log.p=T)))
    if (s == 0){
      hp <- floor(p/2)
      nl <- sqrt((min.dist^2-min.ddist^2*hp)/(p-hp))
      min.logdis <- logdis(c(rep(min.ddist,hp),rep(nl,(p-hp))),s=s)
    } else {
      min.logdis <- logdis(rep(min.dist/sqrt(p),p),s=s)
    }
    samp.lf.cv <- sort(samp.lf, decreasing = T)[n] + 2.5 * min.logdis
    samp.cand.idx <- (samp.lf > samp.lf.cv)
    samp.cand <- samp[samp.cand.idx,]
    samp.cand.lf <- samp.lf[samp.cand.idx]
    med.op <- mined.seq(n, samp.cand, samp.cand.lf, s = s, return.obj = TRUE)
    samp.med <- samp.cand[med.op$idx,]
  }

  feasible.idx <- !(apply(samp.gval, 1, function(x) any(x>0)))
  return (list(med = samp.med,
               cand = samp,
               feasible.idx = feasible.idx))
}

set.seed(20210329)
tau <- c(0,exp(c(1:7)),1e6)
out <- comined(n = 53, p = 2, tau = tau, constraint = constraint,
              n.aug = 5, auto.scale = FALSE, s = 2)

cand <- out$cand
# IMPORTANT: keep matrix shape even if only 1 row is selected
feas <- out$cand[out$feasible.idx, , drop = FALSE]
med <- out$med

# Ensure 2D outputs are always written as two numeric columns (x1,x2)
as_xy <- function(M){
  # Handle the common R pitfall where subsetting drops to a vector
  if (is.null(dim(M))) {
    if (length(M) != 2) stop(paste0("Expected length-2 vector but got length=", length(M)))
    M <- matrix(as.numeric(M), nrow = 1, ncol = 2, byrow = TRUE)
  }
  if (!is.matrix(M)) M <- as.matrix(M)
  if (ncol(M) != 2) stop(paste0("Expected 2 columns (x1,x2) but got ncol=", ncol(M)))
  df <- data.frame(x1 = as.numeric(M[,1]), x2 = as.numeric(M[,2]))
  return(df)
}

write.csv(as_xy(cand), file.path(out_dir, "comined_cand.csv"), row.names = FALSE)
write.csv(as_xy(feas), file.path(out_dir, "comined_feasible.csv"), row.names = FALSE)
write.csv(as_xy(med), file.path(out_dir, "comined_med.csv"), row.names = FALSE)
write.csv(data.frame(x1=as.vector(contour$x1), x2=as.vector(contour$x2), row=rep(1:3, each=1001)), file.path(out_dir, "constraint_contour.csv"), row.names = FALSE)


