## QuantileContourOT(y, z, nR, nS, m, ContLength)
## Function for computing 2D empirical center-outward quantile contours
## y,z column vectors; (y,z) 2d sample
## nR, nS control the size of mini batch (k = nR x nS)
## m is the number of mini batches for which quantiles are computed
## ContLength is the number of points at which the 0.2, 0.4 and 0.8 contours are evaluated

QuantileContourOT<-function(y, z, nR, nS, m, ContLength){
  library(transport)
  library(doParallel)
  nc<-min(m,detectCores())
  n.R<-nR
  n.S<-nS
  registerDoParallel(nc)
  n<-n.R*n.S
  d<-2
  ## generate uniform grid
  disc.sphere<-cbind(cos((2*pi)*(0:(n.S-1))/n.S ),sin((2*pi)*(0:(n.S-1))/n.S ))
  U.grid<-NULL
  for (i in 1:n.R) U.grid<-rbind(U.grid,i/(n.R+1)*disc.sphere)
  ## x mesh
  quantiles<-foreach (u=1:m) %dopar%{
    ### select batch
    ord = sample(1:length(y),n)
    Y = y[ord]
    Z = z[ord]
    random.sample<-cbind(Y,Z)
    normas.y<-rep(0,n)
    for (i in 1:n) normas.y[i]<-sqrt(sum(random.sample[i,]^2))
    ysup<-max(normas.y)
    ##############################################
    ### Step 1: find optimal assignment
    ### We use transport (auctionbf algorithm)
    ### from 'transport' R package
    ##############################################
    asignacion<-transport(pp(U.grid),pp(random.sample),p=2,method='auctionbf')
    asignacion<-asignacion[,2]
    xxx<-U.grid/ysup
    yyy<-random.sample[asignacion,]/ysup
    ##############################################
    ### Step 2b: solve linear program (16)
    ###  via Karp's algorithm
    ##############################################
    cij<-list(apply(xxx*yyy,1,sum))
    cij<-do.call(cbind,rep(cij,n))
    cij<-cij-xxx%*%t(yyy)
    ind.diag<-cbind(1:n,1:n)
    cij[ind.diag]<-rep(Inf,n)
    karp <- karp_cpp(cij)
    mu.star <- karp$mu_star
    shortest.distances <- karp$shortest
    psi<-(-shortest.distances+shortest.distances[1])*ysup^2
    e0<-abs(mu.star)*ysup^2
    ##############################################
    ### Step 3: computation of cyclically monotone
    ### interpolation
    ##############################################
    xxx<-U.grid
    yyy<-random.sample[asignacion,]
    T.0<-function(z){
      auxfun<-function(i){return(sum(z*yyy[i,])-psi[i])}
      scores<-sapply(1:n,auxfun)
      indice<-which.max(scores)
      return(yyy[indice,])
    }
    centr<-T.0(c(0,0))
    disc.sphere<-cbind(cos(2*pi*(0:ContLength)/ContLength),sin(2*pi*(0:ContLength)/ContLength))
    quant1<-t(apply(0.2*disc.sphere,1,T.0))
    quant2<-t(apply(0.4*disc.sphere,1,T.0))
    quant3<-t(apply(0.8*disc.sphere,1,T.0))
    return(list(center=centr,quantile1=quant1,quantile2=quant2,quantile3=quant3))
  }
  center<-NULL
  for (a in 1:m) center<-rbind(center,quantiles[[a]]$center)
  return(list(center=center,quantiles=quantiles,y=y,z=z,m=m,ContLength=ContLength))
}


