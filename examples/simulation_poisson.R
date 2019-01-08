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


# scale-function to scale picture values form [-0.5, 0.5]
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
nx <- 128
x <- seq(1, nx) 
mat <- data.frame("x" = rep(x, each = length(x)), "y" = rep(x, length(x)))
im1 <- matrix(f1(mat$x, mat$y, nx, nx), nrow = length(x))
im1 <- scalefun(im1)
im2 <- matrix(f2(mat$x, mat$y, nx, nx), nrow = length(x))
im2 <- scalefun(im2)
im3 <- matrix(f3(mat$x, mat$y, nx, nx), nrow = length(x))
im3 <- scalefun(im3)


# see original smooth image
filled.contour(x = x, y = x, z = im1, nlevels = 40, 
               color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))) 
filled.contour(x = x, y = x, z = im2, nlevels = 40, 
               color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))) 
filled.contour(x = x, y = x, z = im3, nlevels = 40, 
               color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))) 

# generate z values
m <- 100
#set.seed(042018)
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
eta = Z1 %*% im1v + Z2 %*% im2v + Z3 %*% im3v + X1 * 2.8 + X2 * (-0.2) 
eta <- as.numeric(eta)

# apply Cpp-function due to better speed
cppFunction('NumericVector rngCpp(NumericVector eta) {
                int n = eta.size();
                NumericVector y(n);
                
                for (int i = 0; i < n; ++i) {
                    // return poisson samples from Po(1, lambda = exp(eta))
                    y(i) = R::rpois(exp(eta(i)));
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

system.time({
out <- sarim(y ~ 
              sx(x = X1, knots = 1, penalty = "identity", solver = "rue", 
                ka_start = 50, ka_a = 1, ka_b = 0.00005) +
               sx(Z = Z1, K = K1, penalty = "gmrf", solver = "lanczos", 
                        ka_start = 50, ka_a = 1, ka_b = 0.00005, linear_constraint = "TRUE") + 
                     sx(Z = Z2, K = K2, penalty = "gmrf", solver = "lanczos", 
                        ka_start = 50, ka_a = 1, ka_b = 0.00005, linear_constraint = "TRUE") + 
                     sx(Z = Z3, K = K3, penalty = "gmrf", solver = "lanczos", 
                        ka_start = 50, ka_a = 1, ka_b = 0.00005, linear_constraint = "TRUE") +
                     sx(x = X2, knots = 1, penalty = "identity", solver = "rue", 
                        ka_start = 50, ka_a = 1, ka_b = 0.00005)
             , 
                  family = "poisson", link = "log",
                  data = df, nIter = 100, burnin=111, intercept = "FALSE") 
})

# #Time difference of 1.488806 hours
# n=4000
# mean<-unlist(out$kappa_mean)
# unlist(parallel::mclapply(out$kappa_results,function(x)mean(log(x))))
# var<-unlist(out$kappa_mean2)/(n-1)
# unlist(parallel::mclapply(out$kappa_results,function(x)var(log(x))))
# skew<-sqrt(n)*unlist(out$kappa_mean3)/unlist(out$kappa_mean2)^(1.5)
# 
# for (k in 1:5){
# #cp <- list(mean=mean[k], var.cov=array(var[k], c(1,1)), gamma1=skew[k])
# #dp <- sn::cp2dp(cp, "SN")
# #d<-x<-seq(0,50,by=0.1)
# #for (i in 1:length(d))d[i]<-dmsn(x[i],dp=dp)
# plot(density(out$kappa_results[[k]]))
# #lines(x,d,col="blue")
# lines(seq(0,50,by=0.1),dlnorm(seq(0,50,by=0.1),mean[k],var[k]))
# }
# 


