coeff2<-coda::mcmc.list(
  coda::as.mcmc(as.vector(mod_poisson[[1]]$coef_results[[2]])[-(1:100)]),
  coda::as.mcmc(as.vector(mod_poisson[[2]]$coef_results[[2]])[-(1:100)]),
  coda::as.mcmc(as.vector(mod_poisson[[3]]$coef_results[[2]])[-(1:100)]),
  coda::as.mcmc(as.vector(mod_poisson[[4]]$coef_results[[2]])[-(1:100)]))
plot(coeff2)  
coda::gelman.diag(coeff2)

for (i in 1:52)
{
coeff1<-coda::mcmc.list(
  coda::as.mcmc(as.vector(mod_poisson[[1]]$coef_results[[1]][i,]), start=100),
  coda::as.mcmc(as.vector(mod_poisson[[2]]$coef_results[[1]][i,]), start=100),
  coda::as.mcmc(as.vector(mod_poisson[[3]]$coef_results[[1]][i,]), start=100),
  coda::as.mcmc(as.vector(mod_poisson[[4]]$coef_results[[1]][i,]), start=100))
print(coda::gelman.diag(coeff1))
}