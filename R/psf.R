psrf<-function(out){
  T<-sum(unlist(lapply(out,function(x)return(unlist(x$iterationcounter)))))
nc<-length(out)
m<-parallel::mclapply(out,function(x)return(c(unlist(x$gamma_mean),
  unlist(x$kappa_mean))))
v<-parallel::mclapply(out,function(x)return(c(unlist(x$gamma_mean2),
  unlist(x$kappa_mean2))/(T-1)))
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
cov.wb <- (T/nc) * diag(stats::var(t(v), t(m^2)) - 2 * 
                                  mtotal * stats::var(t(v), t(m)))

V <- (T - 1) * w/T + (1 + 1/nc) * b/T
var.V <- ((T - 1)^2 * var.w + (1 + 1/nc)^2 * var.b + 
            2 * (T - 1) * (1 + 1/nc) * cov.wb)/T^2
df.V <- (2 * V^2)/var.V
df.adj <- (df.V + 3)/(df.V + 1)
B.df <- nc - 1
W.df <- (2 * w^2)/var.w
R2.fixed <- (T - 1)/T
R2.random <- (1 + 1/nc) * (1/T) * (b/w)
R2.estimate <- R2.fixed + R2.random
psrf.my <- sqrt(df.adj * R2.estimate)
psrf.my[is.nan(psrf.my)] <- 2^6
#print(summary(psrf.my))
return(max(psrf.my))
}