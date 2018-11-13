library(RColorBrewer)
library(Rcpp)
library(Matrix)
library(Sarim)

# functions for different smoothing pictures
f1 <- function(x, y, n_x, n_y) {
    value <- (x - n_x/2) * (y - n_y/2)
}
f2 <- function(x, y, n_x, n_y) {
    value <- x - n_x/2 + n_y * sin(y / n_y)
}
f3 <- function(x, y, n_x, n_y) {
    value <- sqrt( (x - n_x / 2)**2 + (y - n_y / 2)**2 )
}


# scale-function to scale piture values form [-0.5, 0.5]
scalefun <- function(z) {
    z_min <- min(z)
    z_max <- max(z)
    zx <- ncol(z)
    zy <- nrow(z)
    y <- matrix(ncol = zx, nrow = zy)
    for (i in 1:zx) {
        for (j in 1:zy) {
            y[i,j] <- -0.5 + (z[i,j] - z_min) / (z_max - z_min)
        }
    }
    return(y)
}

# number of pixels and generate "picture" with scaling form [-0.5, 0.5]
nx <- 120
x <- seq(1, nx) 
mat <- data.frame("x" = rep(x, each = length(x)), "y" = rep(x, length(x)))
im1 <- matrix(f1(mat$x, mat$y, nx, nx), nrow = length(x))
im1 <- scalefun(im1)
im2 <- matrix(f2(mat$x, mat$y, nx, nx), nrow = length(x))
im2 <- scalefun(im2)
im3 <- matrix(f3(mat$x, mat$y, nx, nx), nrow = length(x))
im3 <- scalefun(im3)


# see original image
filled.contour(x = x, y = x, z = im1, nlevels = 40, 
               color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))) 
filled.contour(x = x, y = x, z = im2, nlevels = 40, 
               color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))) 
filled.contour(x = x, y = x, z = im3, nlevels = 40, 
               color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))) 

# generate z values
m <- 100
set.seed(042018)
z1 <- runif(m, -1, 1)
z2 <- runif(m, -1, 1)
z3 <- runif(m, -1, 1)
X1 <- matrix(c(rep(1, m*nx*nx)))
X2 <- matrix(c(sample(c(0,1), m*nx*nx, replace = TRUE)))

Z1 <- as(kronecker(z1, .symDiagonal(nx*nx)), "dgCMatrix")
Z2 <- as(kronecker(z2, .symDiagonal(nx*nx)), "dgCMatrix")
Z3 <- as(kronecker(z3, .symDiagonal(nx*nx)), "dgCMatrix")

im1v <- as.vector(t(im1))
im2v <- as.vector(t(im2))
im3v <- as.vector(t(im3))

# calculate eta
eta = Z1 %*% im1v + Z2 %*% im2v + Z3 %*% im3v + X1 * (-0.2) + X2 * (-0.2) 
eta <- as.numeric(eta)

# apply Cpp-function due to better speed
cppFunction('NumericVector rngCpp(NumericVector eta) {
                int n = eta.size();
                NumericVector y(n);
                
                for (int i = 0; i < n; ++i) {
                    // return binomial samples from Bin(1, size = 5, prob = exp(eta)/(1+exp(eta)))
                    y(i) = R::rbinom( 5, (exp(eta[i])/(1 + exp(eta[i]))) ) ;
                }
                
                return(y);
            }') 
y <- rngCpp(eta)


# create structure matrices
Ps1 <- as(diff(diag(nx), differences = 1), "dgCMatrix")
Ps <- as(crossprod(Ps1, Ps1), "dgCMatrix")
Is <- .symDiagonal(nx)
K1 <- as(kronecker(Ps, Is) + kronecker(Is, Ps), "dgCMatrix")
K2 <- as(kronecker(Ps, Is) + kronecker(Is, Ps), "dgCMatrix")
K3 <- as(kronecker(Ps, Is) + kronecker(Is, Ps), "dgCMatrix")

df <- data.frame("y" = y)


start.time <- Sys.time()
mod_binomial <- Sarim::sarim(y ~ sx(Z = Z1, K = K1, penalty = "gmrf", solver = "lanczos", 
                         ka_start = 50, ka_a = 1, ka_b = 0.00005, linear_constraint = "TRUE") + 
                      sx(Z = Z2, K = K2, penalty = "gmrf", solver = "lanczos", 
                         ka_start = 50, ka_a = 1, ka_b = 0.00005, linear_constraint = "TRUE") + 
                      sx(Z = Z3, K = K3, penalty = "gmrf", solver = "lanczos", 
                         ka_start = 50, ka_a = 1, ka_b = 0.00005, linear_constraint = "TRUE") +
                      sx(x = X1, knots = 1, penalty = "identity", solver = "rue", 
                         ka_start = 50, ka_a = 1, ka_b = 0.00005) +
                      sx(x = X2, knots = 1, penalty = "identity", solver = "rue", 
                         ka_start = 50, ka_a = 1, ka_b = 0.00005), 
                  family = "binomial", link = "logit", Ntrials = 5, 
                  data = df, nIter = 3500, intercept = "FALSE") 

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
# Time difference of 2.161646 hours

# save results for later inspection
setwd("~/Schreibtisch/Master/C++/Abgabe/simulation_binomial/")
saveRDS(mod_binomial, "simulation_binomial_output.rds")
#mod_binomial <- readRDS("simulation_binomial_output.rds")

################################################################################
# take a look at the results
mod_binomial$accept_rate
mod_binomial$lanzcos_iterations


ga1 <- apply(mod_binomial$coef_results[[1]][, 500:3501], 1, FUN = mean)
gamat1 <- matrix(ga1, nrow = nx, byrow = TRUE)
filled.contour(x = x, y = x, z = gamat1, nlevels = 40, 
               color = colorRampPalette(rev(brewer.pal(11, "RdYlBu"))))
gadif1 <- gamat1 - im1
filled.contour(x = x, y = x, z = gadif1, nlevels = 40, zlim = c(-0.5, 0.5),
               color = colorRampPalette(rev(brewer.pal(11, "RdYlBu"))))


ga2 <- apply(mod_binomial$coef_results[[2]][, 500:3501], 1, FUN = mean)
gamat2 <- matrix(ga2, nrow = nx, byrow = TRUE)
filled.contour(x = x, y = x, z = gamat2, nlevels = 40, 
               color = colorRampPalette(rev(brewer.pal(11, "RdYlBu"))))
gadif2 <- gamat2 - im2
filled.contour(x = x, y = x, z = gadif2, nlevels = 40, zlim = c(-0.5, 0.5),
               color = colorRampPalette(rev(brewer.pal(11, "RdYlBu"))))


ga3 <- apply(mod_binomial$coef_results[[3]][, 500:3501], 1, FUN = mean)
gamat3 <- matrix(ga3, nrow = nx, byrow = TRUE)
filled.contour(x = x, y = x, z = gamat3, nlevels = 40, 
               color = colorRampPalette(rev(brewer.pal(11, "RdYlBu"))))
gadif3 <- gamat3 - im3
filled.contour(x = x, y = x, z = gadif3, nlevels = 40, zlim = c(-0.5, 0.5),
               color = colorRampPalette(rev(brewer.pal(11, "RdYlBu"))))


hist(mod_binomial$coef_results[[4]][500:3501], breaks = 50)
ga4 <- mean(mod_binomial$coef_results[[4]][500:3501])
ga4


hist(mod_binomial$coef_results[[5]][500:3501], breaks = 50)
ga5 <- mean(mod_binomial$coef_results[[5]][500:3501])
ga5
