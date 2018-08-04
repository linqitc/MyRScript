beta<-function(x,y){
  s<-(x*y)*0.5/1.96
  a<-((1-x)/(s*x)^2)-x
  b<-a/x-a
  n=1000
  z=rbeta(n,a,b)
  return(z)
}
gamma<-function(x,y) {
  s<-(x*y)*0.5/1.96
  a=(x/s)^2
  b=a/x
  n=1000
  z=rgamma(n,a,b)
  return(z)
}
