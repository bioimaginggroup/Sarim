psrf<-function(out){
  T<-mean(unlist(lapply(out,function(x)return(unlist(x$iterationcounter)))))
nc<-length(out)
m<-parallel::mclapply(out,function(x)return(c(unlist(x$gamma_mean),
  unlist(x$kappa_mean))))
v<-parallel::mclapply(out,function(x)return(c(unlist(x$gamma_mean2),
  unlist(x$kappa_mean2))))
p<-length(m[[1]])
m<-array(unlist(m),c(p,nc))
v<-array(unlist(v),c(p,nc))
mtotal<-apply(m,1,mean)
mmt<-m-mtotal

# vtotal<-(T-1)* (apply(v,1,sum)+(T*apply(mmt^2,1,sum)/(T-1)) / (T*nc-1))

b = T*(apply(mmt^2,1,sum))/(nc-1)
w = apply(v,1,mean)
var.w <- apply(v, 1, stats::var)/nc
var.b <- (2 * b^2)/(nc - 1)

varvec<-function(i,v,m)
{
  return(var(v[i,],m[i,]))
}
C<-dim(v)[1]
#st1 <- stats::var(t(v), t(m))
st1 <- unlist(mclapply(1:C,varvec,v,m))
#st2 <- stats::var(t(v), t(m^2))
st2 <- unlist(mclapply(1:C,varvec,v,m^2))
cov.wb <- (T/nc) * (st2 - 2 * mtotal * st1)

V <- (T - 1) * w/T + (1 + 1/nc) * b/T
var.V <- ((T - 1)^2 * var.w + (1 + 1/nc)^2 * var.b + 
            2 * (T - 1) * (1 + 1/nc) * cov.wb)/T^2
df.V <- (2 * V^2)/var.V

df.adj <- (df.V + 3)/(df.V + 1)
B.df <- nc - 1
W.df <- (2 * w^2)/var.w
R2.random <- (1 + 1/nc) * (1/T) * (b/w)
R2.fixed <- (T - 1)/T
R2.estimate <- R2.fixed + R2.random
psrf.my <- sqrt(abs(df.adj * R2.estimate))
print(summary(psrf.my))
#plot(hist(psrf.my))
print(which(psrf.my==max(psrf.my)))
return(max(psrf.my))
#return(quantile(psrf.my,.9))
}