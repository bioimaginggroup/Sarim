library(BayGMRF)
data(scotland)

library(Matrix)
X<-as.matrix(scotland.data$X)
Q<-BayGMRF::weightedQ(rep(1,117),scotland.coords)
y<-scotland.data$Counts
data<-list("y"=y)
Z<-diag(52)

library(Sarim)
system.time(mod_poisson <- sarim(y ~ 
                      sx(Z = Z, K = Q, penalty = "gmrf", solver = "lanczos", 
                            ka_start = 50, ka_a = 1, ka_b = 0.00005, linear_constraint = "TRUE") 
                     + 
                       sx(x = X, knots = 1, penalty = "identity", solver = "rue", 
                         ka_start = 50, ka_a = 1, ka_b = 0.00005)
                     , 
                     family = "poisson", link = "log",
                     data = data, nIter = 1000, intercept = "FALSE"))


coeff1<-mod_poisson$coef_results[[1]]
coeff2<-mod_poisson$coef_results[[2]]
co1<-apply(coeff1,1,median)
map(co1,scotland.shape)


true<-rep(0,52)
true[18:40]<-1
true[23:25]<-true[29:32]<-true[35:37]<-2
scotland.sim<-scotland.data[,-3]
scotland.sim$E<-rpois(52,scotland.data$E-1)+1
globalgamma<-sample(c(-.1,0,.1))
gamma.true<-rnorm(52,globalgamma[true+1],.01)
map(gamma.true,scotland.shape)
scotland.sim$Counts<-rpois(52,exp(gamma.true+scotland.data$X*0.2))
map(scotland.sim$Counts,scotland.shape)

library(Matrix)
X<-as.matrix(scotland.data$X)
Q<-BayGMRF::weightedQ(rep(1,117),scotland.coords)
y<-scotland.data$Counts
data<-list("y"=scotland.sim$Counts)
Z<-diag(52)

library(Sarim)
modsim <- sarim(y ~ 
                       sx(Z = Z, K = Q, penalty = "gmrf", solver = "lanczos", 
                          ka_start = 50, ka_a = 1, ka_b = 0.00005, linear_constraint = "TRUE") 
                     + 
                       sx(x = X, knots = 1, penalty = "identity", solver = "rue", 
                          ka_start = 50, ka_a = 1, ka_b = 0.00005)
                     , 
                     family = "poisson", link = "log",
                     data = data, nIter = 3500, intercept = "FALSE") 

coeff1<-modsim$coef_results[[1]]
coeff2<-modsim$coef_results[[2]]
plot(density(coeff2))
co1<-apply(coeff1,1,median)
map(co1,scotland.shape)
plot(co2,gamma.true)
