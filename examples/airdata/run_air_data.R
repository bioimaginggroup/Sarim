load("~/software/Sarim/examples/airdata/prep.Rdata")
y=vecdat$value
n<-dim(vecdat)[1]
p<-dim(Q)[1]
Z1<-sparseMatrix(i=1:n,j=vecdat$id,x=rep(1,n),dims=c(n,p))
Q1<-Q+0.5*Diagonal(p)
library(Sarim)
form<- "y ~ sx(Z = Z1, K = Q1, penalty = 'gmrf', solver = 'lanczos', 
              ka_start = 50, ka_a = 1, ka_b = 1)"
df <- data.frame("y" = y)
mf <- stats::model.frame(formula = form,data=df)
erg<-sarim(form, data=df)

res<-erg$coef_results[[1]]
spateff<-apply(res,1,mean)
spateff<-matrix(spateff,c(length(latgrid),length(longrid)))
fields::image.plot(spateff)
plot(erg$kappa_results[[1]][-1])
plot(res[1234,])
