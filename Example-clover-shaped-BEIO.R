source('KarpCPP.R')
source('CenterOutwardContourOT.R')

## Example 1
set.seed(532342)

n.R<-20
n.S<-20
nn<-n.R*n.S

n<-1000

#### Clover shaped model
rmixture<-function(a){
  aux<-runif(1)
  o1<-c(0,0)
  o2<-10*c(cos(-pi/6),sin(-pi/6))
  o3<-10*c(cos(7*pi/6),sin(7*pi/6))
  o4<-c(0,10)
  result<-(rnorm(2)+(aux<=0.25)*o1+(aux>0.25)*(aux<=0.5)*o2+(aux>0.5)*(aux<=0.75)*o3+(aux>0.75)*o4)/10
  return(result)
}

rotation<-function(x){
  angle<-pi*x/2
  rotmatrix<-matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)),ncol=2)
  return(rotmatrix)
}

random.sample<-t(sapply(1:n,rmixture))%*%rotation(-pi/3)


normas.y<-rep(0,n)
for (i in 1:n) normas.y[i]<-sqrt(sum(random.sample[i,]^2))
ysup<-max(normas.y)
lim.w<-max(abs(random.sample))


normas.y<-rep(0,n)
for (i in 1:n) normas.y[i]<-sqrt(sum(random.sample[i,]^2))
ysup<-max(normas.y)
lim.w<-max(abs(random.sample))



m<-150

start.time<-Sys.time()
fit1<-QuantileContourOT(random.sample[,1],random.sample[,2],n.R,n.S,m, ContLength = 50)
end.time<-Sys.time()
time1<-end.time-start.time
print(time1)


quant1<-quant2<-quant3<-matrix(0,nrow=dim(fit1$quantiles[[1]]$quantile1)[1],ncol=dim(fit1$quantiles[[1]]$quantile1)[2])
for (i in 1:m){
  quant1<-quant1+fit1$quantiles[[i]]$quantile1
  quant2<-quant2+fit1$quantiles[[i]]$quantile2
  quant3<-quant3+fit1$quantiles[[i]]$quantile3
}
quant1<-quant1/m
quant2<-quant2/m
quant3<-quant3/m


png('Fig-clover-shaped-BEIO-alt.png',width=500,height=500)
plot(random.sample,pch=21,col=rgb(0.7,0.7,0.7,0.01),xlab=expression(x[1]),ylab=expression(x[2]))
lines(quant1[,1],quant1[,2],col='darkolivegreen',lwd=2)
lines(quant2[,1],quant2[,2],col='darkolivegreen3',lwd=2)
lines(quant3[,1],quant3[,2],col='darkolivegreen1',lwd=2)
dev.off()
