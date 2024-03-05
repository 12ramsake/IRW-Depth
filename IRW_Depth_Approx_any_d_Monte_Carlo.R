getUnivariateDepth<-function(u,data,even){
  
  ##all the points projected onto u
  projections<-data%*%(u)
  ##sample ranks
  n<-nrow(data)
  temp<-rank(projections)

    ##sample ranks
    if(!even)
      maxD<-(n+1)/2
    else
      maxD<-n/2
    
    change<-function(r){
      return(n-r+1)
    }
    temp[temp>maxD]<-change(temp[temp>maxD])
  
  
  return(temp[length(temp)])#returns index of points
}

#this is the one
getUnivariateDepth<-function(u,point,data){
  
  projections<-data%*%(u)
  xp<-sum(u*c(point))
  Fn<-ecdf(projections)
  
  return(min(Fn(xp),1-Fn(xp)))#returns index of points
}


calculateDepth<-function(point,data,N=10000,sim=F,x11){
  
  normalize<-function(y){
    return(y/sqrt(sum(y^2)))
  }
  

    if(!sim){
      d<-ncol(data)
      x11<-replicate(d,rnorm(N))
    }
    #ncol=d
    #turn to unit vectors
    x1<-apply(x11,MARGIN=1,normalize)
 
#  even<-(nrow(data)+1)%%2==0
  depths<-c()
  data<-apply(data,2,as.numeric)
  for(i in 1:N)
    depths[i]<-getUnivariateDepth(matrix(x1[,i],ncol=1),point=point,data=data)
    #depths[i]<-getUnivariateDepth(matrix(x1[,i],ncol=1),data=rbind(data,point),even=even)
  
  return(2*mean(depths))
}







## testing 



testData<-replicate(2,rnorm(20))
calculateDepth(c(1,2),testData,N=1000)
getPointDepth2(c(1,2),testData)

## testing 

data("iris")

c1<-subset(iris,Species=="virginica")
c2<-subset(iris,Species=="versicolor")
c3<-subset(iris,Species=="setosa")

N<-1000
d<-4
x1<-replicate(d,rnorm(N))

virD<-apply(c1[,-5],1,calculateDepth,data=c1[,-5],sim=T,N=N,x11=x1)
verD<-apply(c2[,-5],1,calculateDepth,data=c2[,-5],sim=T,N=N,x11=x1)
setD<-apply(c3[,-5],1,calculateDepth,data=c3[,-5],sim=T,N=N,x11=x1)


ddPlot<-function(x,y){
  
  qqplot(virD,verD,asp=1)
  abline(0,1)
  
}

ddPlot(virD,verD)
ddPlot(virD,setD)
ddPlot(setD,verD)

plot(c1[,1:2],col="red")
points(c2[,1:2],col="blue")
points(c3[,1:3],col="pink")
dim(c1)
dim(c2)
dim(c3)

train<-rbind(sample(1:50,40,replace = F),
         sample(1:50,40,replace = F),
         sample(1:50,40,replace = F)
         )


test<-rbind(c1[-train[1,],-5],c2[-train[2,],-5],c3[-train[3,],-5])
vir<-apply(test,1,calculateDepth,data=c1[train[1,],-5],sim=T,N=N,x11=x1)
ver<-apply(test,1,calculateDepth,data=c2[train[2,],-5],sim=T,N=N,x11=x1)
set<-apply(test,1,calculateDepth,data=c3[train[3,],-5],sim=T,N=N,x11=x1)

rs<-t(apply(cbind(vir,ver,set),1,rank))
pred<-c()
for(i in 1:30){
  for(j in 1:3){
    if(rs[i,j]==3)
      pred[i]<-j
  }
}


err<-(sum(pred[1:10]!=1)+sum(pred[11:20]!=2)+sum(pred[21:30]!=3))
err<-err/30
err


