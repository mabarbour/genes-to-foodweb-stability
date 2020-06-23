## Code for critical boundaries and plotting feasibility domains ----

# Author: Matthew A. Barbour

## Load required libraries ----
library(rgl)

## Code for a perfectly round sphere ----
# note that I don't use this for any visualizations in the manuscript
# https://stackoverflow.com/questions/39778093/how-to-increase-smoothness-of-spheres3d-in-rgl
sphere1.f <- function(x0 = 0, y0 = 0, z0 = 0, r = 1, n = 101, ...){
  f <- function(s,t){
    cbind(   r * cos(t)*cos(s) + x0,
             r *        sin(s) + y0,
             r * sin(t)*cos(s) + z0)
  }
  persp3d(f, slim = c(-pi/2,pi/2), tlim = c(0, 2*pi), n = n, add = T, ...)
}

## Initialize rgl device for 3D visualization of feasibility domain for 3 species ----
# from: http://www.sthda.com/english/wiki/a-complete-guide-to-3d-visualization-device-system-in-r-r-software-and-data-visualization
rgl_init <- function(new.device = FALSE, bg = "white", width = 640) {
  if( new.device | rgl.cur() == 0 ) {
    rgl.open()
    par3d(windowRect = 50 + c( 0, 0, width, width ) )
    rgl.bg(color = bg )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
  rgl.viewpoint(theta = 15, phi = 20, zoom = 0.7)
}

## Calculate cross product between two vectors ----
CrossProduct <- function(a, b){
  # https://www.mathsisfun.com/algebra/vectors-cross-product.html
  x <- 1
  y <- 2
  z <- 3
  matrix(c((a[y]*b[z]-a[z]*b[y]),
           (a[z]*b[x]-a[x]*b[z]),
           (a[x]*b[y]-a[y]*b[x])),
         ncol = 1)
}

## Calculate the normalized angles from critical boundary for initial food web ----
FeasibilityBoundary <- function(A, r){
  origin <- matrix(c(0,0,0), ncol = 1)
  # normalize each intrinsic growth rate to have length = 1
  r <- r / sqrt(r[1]^2 + r[2]^2 + r[3]^2)
  # normalize each column vector to have length = 1
  for(j in 1:ncol(A)){
    A[,j] <- A[,j] / sqrt(A[1,j]^2 + A[2,j]^2 + A[3,j]^2)
  }
  # calculate normal vector (https://web.ma.utexas.edu/users/m408m/Display12-5-4.shtml)
  NV.A12 <- CrossProduct(-A[,1], -A[,2])
  # we multiply this by -1 so that normal vector points into the
  # feasible parameter space. It also makes it so that negative
  # alphas indicate how far away the growth rate vector is from the
  # feasibility boundary. It also gives insight to which feasibility
  # boundary defines the structural stability of the system.
  NV.A13 <- CrossProduct(-A[,1], -A[,3])
  NV.A23 <- CrossProduct(-A[,2], -A[,3])
  # calculate angle between normal vector and growth rate vector
  theta.A12 <- angle(c(r), c(NV.A12), degree = T)
  theta.A13 <- angle(c(r), c(NV.A13), degree = T)
  theta.A23 <- angle(c(r), c(NV.A23), degree = T)
  # calculate angle between growth rate vector and plane
  # note that it is not exactly possible to determine where
  # the "correct" angle is. This has to be checked manually
  alpha.A12 <- 90 - theta.A12
  alpha.A13 <- 90 - theta.A13
  alpha.A23 <- 90 - theta.A23
  if(all(-1*inv(A) %*% r > 0) == TRUE){
    feasibility <- 1
  } else{
    feasibility <- 0
  }
  # output
  return(c(alpha.A12 = alpha.A12, alpha.A13 = alpha.A13, alpha.A23 = alpha.A23, feasibility = feasibility))
}

## Calculate normalized angle from critical boundary for LYER and Ptoid ----
BoundaryLYER.Ptoid <- function(A, r){
  # normalize each intrinsic growth rate to have length = 1
  r <- r / sqrt(r[1]^2 + r[2]^2)
  # normalize each column vector to have length = 1
  for(j in 1:ncol(A)){
    A[,j] <- A[,j] / sqrt(A[1,j]^2 + A[2,j]^2)
  }
  if(all(-1*inv(A) %*% r > 0) == TRUE){
    feasibility <- 1
    boundary <- angle(-1*A[,1], r)
  } else{
    feasibility <- 0
    boundary <- -1*angle(-1*A[,1], r)
  }
  return(c(boundary = boundary, feasibility = feasibility))
}

## Plot feasibility domain for 3 species ----
PlotFeasibilityDomain3sp <- function(A, # list of Interaction Matrices
                                    r, # list of intrinsic growth rates
                                    A.color, # vector of colors for Interaction Matrices
                                    r.color, # vector of colors for intrinsic growth rates
                                    sphere.color = "grey",
                                    normalize = TRUE, # normalize vectors to length = 1?
                                    species.labels = 1:ncol(A[[1]]), # default to column order, otherwise put in column order
                                    constraints = NULL, # specific additional constraints on feasibility domain
                                    sphere.alpha = 0.1, # alpha for sphere
                                    barb.n = 100,
                                    arc.width = rep(2, ncol(A[[1]])) # default to size 2 for all
                                    ){

  # plot origin
  rgl_init()
  # plot sphere
  sphere1.f(col = "grey", alpha = sphere.alpha)

  if(normalize == TRUE){
    # normalize each intrinsic growth rate to have length = 1
    for(i in 1:length(r)){
      r[[i]] <- r[[i]] / sqrt(r[[i]][1]^2 + r[[i]][2]^2 + r[[i]][3]^2)

      # normalize each column vector to have length = 1
      for(j in 1:ncol(A[[i]])){
        A[[i]][,j] <- A[[i]][,j] / sqrt(A[[i]][1,j]^2 + A[[i]][2,j]^2 + A[[i]][3,j]^2)
      }
      # label column vectors (i.e. critical boundaries)
      rgl.texts(x = -A[[i]][1,1],
                y = -A[[i]][2,1],
                z = -A[[i]][3,1], col = A.color[i], text = species.labels[1], pos = 3)
      rgl.texts(x = -A[[i]][1,2],
                y = -A[[i]][2,2],
                z = -A[[i]][3,2], col = A.color[i], text = species.labels[2], pos = 3)
      rgl.texts(x = -A[[i]][1,3],
                y = -A[[i]][2,3],
                z = -A[[i]][3,3], col = A.color[i], text = species.labels[3], pos = 3)
      # region of compatible intrinsic growth rates
      arc3d(from = c(-A[[i]][1,1], -A[[i]][2,1], -A[[i]][3,1]),
            to = c(-A[[i]][1,2], -A[[i]][2,2], -A[[i]][3,2]),
            center = c(0,0,0),
            col = A.color[i], radius = 1, lwd = arc.width[i])
      arc3d(from = c(-A[[i]][1,1], -A[[i]][2,1], -A[[i]][3,1]),
            to = c(-A[[i]][1,3], -A[[i]][2,3], -A[[i]][3,3]),
            center = c(0,0,0),
            col = A.color[i], radius = 1, lwd = arc.width[i])
      arc3d(from = c(-A[[i]][1,2], -A[[i]][2,2], -A[[i]][3,2]),
            to = c(-A[[i]][1,3], -A[[i]][2,3], -A[[i]][3,3]),
            center = c(0,0,0),
            col = A.color[i], radius = 1, lwd = arc.width[i])
      # doubling and reversing direction, this fixes the fact
      # that the arc wasn't extending all of the way
      # for some unknown reason...
      arc3d(from = c(-A[[i]][1,2], -A[[i]][2,2], -A[[i]][3,2]),
            to = c(-A[[i]][1,1], -A[[i]][2,1], -A[[i]][3,1]),
            center = c(0,0,0),
            col = A.color[i], radius = 1, lwd = arc.width[i])
      arc3d(from = c(-A[[i]][1,3], -A[[i]][2,3], -A[[i]][3,3]),
            to =  c(-A[[i]][1,1], -A[[i]][2,1], -A[[i]][3,1]),
            center = c(0,0,0),
            col = A.color[i], radius = 1, lwd = arc.width[i])
      arc3d(from = c(-A[[i]][1,3], -A[[i]][2,3], -A[[i]][3,3]),
            to =  c(-A[[i]][1,2], -A[[i]][2,2], -A[[i]][3,2]),
            center = c(0,0,0),
            col = A.color[i], radius = 1, lwd = arc.width[i])

      # plot vector of intrinsic growth rates
      lines3d(x = c(0,r[[i]][1]),
              y = c(0,r[[i]][2]),
              z = c(0,r[[i]][3]), col = r.color[i])
      points3d(x = c(r[[i]][1]),
               y = c(r[[i]][2]),
               z = c(r[[i]][3]), col = r.color[i], size = 10, pch = 1)
    }}
  else{
    for(i in 1:length(r)){
    # region of compatible intrinsic growth rates
    arc3d(from = c(-A[[i]][1,1], -A[[i]][2,1], -A[[i]][3,1]),
          to = c(-A[[i]][1,2], -A[[i]][2,2], -A[[i]][3,2]),
          center = c(0,0,0),
          col = A.color[i], radius = 1)
    arc3d(from = c(-A[[i]][1,1], -A[[i]][2,1], -A[[i]][3,1]),
          to = c(-A[[i]][1,3], -A[[i]][2,3], -A[[i]][3,3]),
          center = c(0,0,0),
          col = A.color[i], radius = 1)
    arc3d(from = c(-A[[i]][1,2], -A[[i]][2,2], -A[[i]][3,2]),
          to = c(-A[[i]][1,3], -A[[i]][2,3], -A[[i]][3,3]),
          center = c(0,0,0),
          col = A.color[i], radius = 1)
    # doubling and reversing direction, this fixes the fact
    # that the arc wasn't extending all of the way
    # for some unknown reason...
    arc3d(from = c(-A[[i]][1,2], -A[[i]][2,2], -A[[i]][3,2]),
          to = c(-A[[i]][1,1], -A[[i]][2,1], -A[[i]][3,1]),
          center = c(0,0,0),
          col = A.color[i], radius = 1)
    arc3d(from = c(-A[[i]][1,3], -A[[i]][2,3], -A[[i]][3,3]),
          to =  c(-A[[i]][1,1], -A[[i]][2,1], -A[[i]][3,1]),
          center = c(0,0,0),
          col = A.color[i], radius = 1)
    arc3d(from = c(-A[[i]][1,3], -A[[i]][2,3], -A[[i]][3,3]),
          to =  c(-A[[i]][1,2], -A[[i]][2,2], -A[[i]][3,2]),
          center = c(0,0,0),
          col = A.color[i], radius = 1)
    # plot vector of intrinsic growth rates
    points3d(x = c(r[[i]][1]),
             y = c(r[[i]][2]),
             z = c(r[[i]][3]), col = r.color[i], size = 7)
    }
  }
  if(length(constraints) > 0){
    # cut sphere into constrained region of intrinsic growth rates compatible with the biology of the system
    clipplanes3d(a = constraints[1], b = 0, c = 0, d = 0)
    clipplanes3d(a = 0, b = constraints[2], c = 0, d = 0)
    clipplanes3d(a = 0, b = 0, c = constraints[3], d = 0)
  }
}

## Plot feasibility domain for 2 species ----
PlotFeasibilityDomain2sp <- function(A, # list of Interaction Matrices
                                     r, # list of intrinsic growth rates
                                     labels, # labels for Interaction Matrices and growth rates
                                     normalize # normalize vectors to length = 1?
){

  if(normalize == TRUE){
    # normalize each intrinsic growth rate to have length = 1
    for(i in 1:length(r)){
     r[[i]] <- r[[i]] / sqrt(r[[i]][1]^2 + r[[i]][2]^2)

     # normalize each column vector to have length = 1
     for(j in 1:ncol(A[[i]])){
       A[[i]][,j] <- A[[i]][,j] / sqrt(A[[i]][1,j]^2 + A[[i]][2,j]^2)
     }
    }}
  else{
      r <- r
      A <- A
  }

  # data frame for plotting feasibility domain
  # of each Interaction Matrix with corresponding
  # intrinsic growth rates
  df.A <- list()
  for(i in 1:length(A)){
    df.A[[i]] <- data.frame(
      A_ID = rep(labels[i], 3),
      Type = c("A","A","r"),
      Sp_1 = c(-A[[i]][1,1], -A[[i]][1,2], r[[i]][1]),
      Sp_2 = c(-A[[i]][2,1], -A[[i]][2,2], r[[i]][2])
    )
  }
  plyr::ldply(df.A) %>%
    ggplot(., aes(x = Sp_1, y = Sp_2, color = factor(A_ID))) +
    geom_segment(aes(x = 0, y = 0, xend = Sp_1, yend = Sp_2, linetype = Type),
                 arrow = arrow(length = unit(0.1,"cm")))
}

## Get feasibility domain for 2 species for more tailored plotting ----
FeasibilityDomain2sp <- function(A, # list of Interaction Matrices
                                 r, # list of intrinsic growth rates
                                 labels, # labels for Interaction Matrices and growth rates
                                 normalize # normalize vectors to length = 1?
){

  if(normalize == TRUE){
    # normalize each intrinsic growth rate to have length = 1
    for(i in 1:length(r)){
      r[[i]] <- r[[i]] / sqrt(r[[i]][1]^2 + r[[i]][2]^2)

      # normalize each column vector to have length = 1
      for(j in 1:ncol(A[[i]])){
        A[[i]][,j] <- A[[i]][,j] / sqrt(A[[i]][1,j]^2 + A[[i]][2,j]^2)
      }
    }}
  else{
    r <- r
    A <- A
  }

  # data frame for plotting feasibility domain
  # of each Interaction Matrix with corresponding
  # intrinsic growth rates
  df.A <- list()
  for(i in 1:length(A)){
    df.A[[i]] <- data.frame(
      A_ID = rep(labels[i], 3),
      Type = c("A","A","r"),
      Sp_1 = c(-A[[i]][1,1], -A[[i]][1,2], r[[i]][1]),
      Sp_2 = c(-A[[i]][2,1], -A[[i]][2,2], r[[i]][2])
    )
  }
  plyr::ldply(df.A)
}
