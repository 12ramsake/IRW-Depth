getUnivariateDepth<-function(u,point,data){
  
  projections<-as.matrix(data)%*%(u)
  xp<-sum(u*point)
  Fn<-ecdf(projections)
  addition<-any(apply(data,1,identical,y=point))
  return(min(Fn(xp),1-Fn(xp)+addition/nrow(data)))#returns index of points
}

getUnivariateDepthSimpl<-function(u,point,data){
  
  projections<-as.matrix(data)%*%(u)
  xp<-sum(u*point)
  Fn<-ecdf(projections)
  return(Fn(xp)*(1-Fn(xp)))#returns index of points
}


findCritPoints<-function(data){
    
    n<-length(data[,1])
    data2<-cbind(data[combn(1:n,2)[1,],],data[combn(1:n,2)[2,],])
    
    cp<-atan((data2[,1]-data2[,3])/(data2[,4]-data2[,2]))
    if(any(data2[,4]==data2[,2])){
      cp[data2[,4]==data2[,2]]=pi/2
      
    }
    print(cp)
    if(any(cp<0))
      cp[cp<0]<-cp[cp<0]+pi
    
    
    #now have n choose 2 possible switching points
    return(cp)
  }
  
findXtraCritPoints<-function(newPoint,data){
    
    rows<-matrix(TRUE,nrow=nrow(data),ncol=1)
    for(i in 1:nrow(data))
      rows[i,1]=sum(newPoint==data[i,])!=2
    
    cp<-atan((data[rows[,1],1]-unlist(newPoint[1]))/(unlist(newPoint[2])-data[rows[,1],2]))
    if(any(cp<0))
      cp[cp<0]<-cp[cp<0]+pi
    
    
    #now have n choose 2 possible switching points
    return(cp)
  }
  
getRanks<-function(x,data,even,depth=T){
    
    ##all the points projected onto u
    projections<-as.matrix(data)%*%t(t(x))
    ##sample ranks
    
    n<-nrow(data)
    temp<-rank(projections)
    
    if(depth){
      ##sample ranks
      if(!even)
        maxD<-(n+1)/2
      else
        maxD<-n/2
      
      change<-function(r){
        return(n-r+1)
      }
      temp[temp>maxD]<-change(temp[temp>maxD])
    }
    
    return(temp)#returns index of points
  }
  
rank2depth<-function(n,temp,even){
    
    
    ##sample ranks
    if(!even)
      maxD<-(n+1)/2
    else
      maxD<-n/2
    
    change<-function(r){
      return(n-r+1)
    }
    temp[temp>maxD]<-change(temp[temp>maxD])
    
    
    return(temp)#returns index of points
  }
  #n is max rank
rank2depth2<-function(n,temp){
    
    depth<-0:(n-1)/(n-1)
    ##sample ranks
    depth<-apply(cbind(depth[temp],1-depth[temp]),1,min)
    
    return(depth)#returns index of points
  }
  
# Use if multiplicities =1
getPointDepth<-Vectorize(function(x,y,data){
  
  N<-nrow(data)
  criticalPoints<-findXtraCritPoints(c(x,y),data)
  #the sort
  order<-order(criticalPoints)
  ##add extra slot
  angles2<-(criticalPoints[order[1]])/2
  ##cartesian
  point<-c(cos(angles2),sin(angles2))
  ##even number of points?
  even<-(N+1)%%2==0
  #get RANKS for first section
  ranks<-getRanks(point,data=rbind(data,c(x,y)),even=even,depth=F)
  ##vector of ranks at each section
  pointRanks<-ranks[N+1]
  
  ##adds the next rank
  incRank<-function(initRanks,swapIndex,depth){
    
    ##if it had larger rank it is on right side, 
    ##then the rank of the new point swaps, 
    ##and is increased by 1
    if(initRanks[swapIndex]>depth[1])
      {depth[length(depth)+1]=depth[length(depth)]+1}
    
    else
      {depth[length(depth)+1]=depth[length(depth)]-1}
    
    return(depth)
  }
  
  #weight for each rank
  for(i in 1:N){pointRanks<-incRank(ranks,order[i],pointRanks)}
  
  
  weights<-rbind(c(0,criticalPoints[order][1]),
                 cbind(criticalPoints[order][1:(N-1)],
                       criticalPoints[order][2:(N)]),
                 c(criticalPoints[order][N],pi))
  w<-(weights[,2]-weights[,1])/pi
  
#  depths<-rank2depth(N+1,pointRanks,even)
  depths<-rank2depth2(N+1,pointRanks)
  depth<-sum(depths*w)
  
  
   return(2*depth)
},vectorize.args = c('x','y'))


# Use if multiplicities >1
getPointDepth3<-Vectorize(function(x,y,data,Cue=F){
  
  N<-nrow(data)
  criticalPoints<-sort(findXtraCritPoints(c(x,y),data))
  #the sort
 
  N2<-length(criticalPoints)
  weights<-rbind(c(0,criticalPoints[1]),
                 cbind(criticalPoints[1:(N2-1)],
                       criticalPoints[2:(N2)]),
                 c(criticalPoints[N2],pi))
  
  depths<-c()
  
  #weight for each rank
  angles<-(weights[,2]-weights[,1])/2+weights[,1]
  
  for(i in 1:nrow(weights)){
    if(!Cue)
      depths[i]<-getUnivariateDepth(matrix(c(cos(angles[i]),sin(angles[i])),ncol=1),c(x,y),as.matrix(data))
    else
      depths[i]<-getUnivariateDepthSimpl(matrix(c(cos(angles[i]),sin(angles[i])),ncol=1),c(x,y),as.matrix(data))*2
    
    }
  
  
  w<-(weights[,2]-weights[,1])/pi
  
  #  depths<-rank2depth(N+1,pointRanks,even)
  #depths<-rank2depth2(N+1,pointRanks)
  depth<-sum(depths*w)
  
  
  return(2*depth)
},vectorize.args = c('x','y'))

# Cue is for cuevas IDD depth, set to F for IRW
getPointDepth2<-function(theta,data,Cue=F){
  if(anyDuplicated(rbind(data,theta))==0&&!Cue)
    return(getPointDepth(theta[1],theta[2],data))
  else
    return(getPointDepth3(theta[1],theta[2],data,Cue))
}

get