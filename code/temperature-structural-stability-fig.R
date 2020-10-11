# run locally rather on my Rstudio server
load('output/structural-stability-keystone.RData')
source('code/plot-feasibility-domain.R')
library(rgl)

# plot structural stability ----
constraints.matrix <- matrix(c(-1,0,0,
                               0,-1,0,
                               0,0,1),
                             ncol = 3, byrow = T)
PlotFeasibilityDomain3sp(A = list(temp20.mat, temp23.mat, constraints.matrix),
                         r = list(temp20.IGR, temp23.IGR, temp23.IGR), # IGR doesn't matter for constraints
                         A.color = c("steelblue","firebrick1","grey"),
                         r.color = c("steelblue","firebrick1","grey"),
                         normalize = TRUE,
                         sphere.alpha = 0,
                         arc.width = c(2,2,2),
                         barb.n = 2,
                         species.labels = c("","",""))
# add axes
rgl.lines(x = c(0,1.1), y = c(0,0), z = c(0,0), color = "black", lwd = 3)
rgl.lines(x = c(0,0), y = c(1.1,0), z = c(0,0), color = "black", lwd = 3)
rgl.lines(x = c(0,0), y = c(0,0), z = c(-1.1,0), color = "black", lwd = 3)

scene3d()
rglwidget()

# print snapshot
# manually rotated into position I liked
rgl.snapshot("figures/initial-foodweb-structural-stability.png")
