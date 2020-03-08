
###########################
# getDist
###########################

library(pdist)
library(plyr)

#computes weights for a given data, standalone
getWeights<-function(xx, out,cv=FALSE, train.snames=NULL){
   #add a check for missing rownames
  if(is.null(rownames(xx))){stop("provide rownames in getDist")}
  #how can we compute faster HRs, apply is faster than for 
  
  if(cv==FALSE){
    #always intersect with out
    inter = intersect(names(out),rownames(xx))
    xx = xx[inter,,drop=FALSE]
    out = out[inter]
    
    #calculate survobj upfront
    ll<-list()
    for(i in 1:ncol(xx)){
      if(length(unique(out))==2){
        ll[[i]]<-summary(glm(out ~ xx[,i], family="binomial"))$coefficients
        }
      
      if(length(unique(out)) > 2){
        ll[[i]]<-sqrt(kruskal.test(xx[,i] ~ as.factor(out))$statistic)
      }
    }
    if(length(unique(out)) == 2){
      lodds<-unlist(lapply(ll, function(x) x[2,1]))}
    
    if(length(unique(out)) > 2){
      lodds<-unlist(ll)}
    
    xx.wt <- t(t(xx) * sqrt(abs(lodds)))
    return(xx.wt)
  }
  
  if(cv==TRUE){
    inter = intersect(names(out), intersect(rownames(xx),train.snames))
    xx.train = xx[inter,,drop=FALSE]
    out = out[inter]
    
    ll<-list()
    for(i in 1:ncol(xx)){
      
      if(length(unique(out)) == 2){
        ll[[i]]<-summary(glm(out ~ xx.train[,i], family="binomial"))$coefficients
      }
      
      if(length(unique(out)) > 2){
        
        ll[[i]]<-sqrt(kruskal.test(xx.train[,i] ~ as.factor(out))$statistic)
      }  
    }
    
    if(length(unique(out)) == 2){
    lodds<-unlist(lapply(ll, function(x) x[2,1]))}
    
    if(length(unique(out)) > 2){
      lodds<-unlist(ll)}
    
    #on features which are columns, same in all and train
    xx.train.wt <- t(t(xx.train) * sqrt(abs(lodds)))
    xx.wt <- t(t(xx) * sqrt(abs(lodds)))
    wt.mat<-list(all=xx.wt, train=xx.train.wt)
    return(wt.mat)
  }
  
}

getUnionDist<-function(rnames,dat, type=NULL){
  #compute distance across union of genes/features

  dist.dat<-list()

  for (i in 1:length(dat)){
    m <- dat[[i]][intersect(rownames(dat[[i]]),rnames),,drop=FALSE]
    #if there are missing samples, add row of NAs for it
    if (length(intersect(rownames(m),rnames)) != length(rnames)){
      m.na<-matrix(NA,nrow=length(setdiff(rnames,rownames(m))),ncol=ncol(m))
      rownames(m.na) = setdiff(rnames,rownames(m))
      m<-rbind(m,m.na)
      m<-m[rnames,]
    }
    if(!(is.null(type))){
      if(type=="mut"){dist.dat[[i]] = dist_wtbinary(m)
                    rownames(dist.dat[[i]]) = colnames(dist.dat[[i]]) = rownames(m)}}
    
    if(is.null(type)){    
      m2 = m*m
      m2.ss = sum(m2, na.rm=T)
      m.tr = m/sqrt(m2.ss)
      dist.dat[[i]] = as.matrix(dist(m.tr, method="euclidean"))
    }
  }
  return(dist.dat)
}


getDist<-function(datasets,out,cv=FALSE,train.snames=NULL,type=NULL){
  # add other checks
  if (is.list(datasets) == F)
    stop("datasets must be a list")
  
  #if out has no events stop 
  if (length(unique(out))==1)
    stop("Need more than one category in outcome")
  
  #convert everything to numeric
  dat<-lapply(datasets, function(x) as.data.frame(aaply(x,1,as.numeric,.drop=FALSE)) )
  rnames<-unlist(lapply(datasets, function(x) rownames(x)))
  rnames<-unique(rnames)

  if(is.null(rnames))
    stop("rowanmes=NULL, add sample names to matrix of datasets list object")

  dat.wt<-lapply(dat, function(x) getWeights(x,out,cv,train.snames))

  if(cv==TRUE){
    dat.all.wt<-lapply(dat.wt, function(x) x$all)
    dat.train.wt<-lapply(dat.wt, function(x) x$train)

    dat.all.dist<-getUnionDist(rnames, dat.all.wt,type)
    dat.train.dist<-getUnionDist(train.snames, dat.train.wt, type)
    dat.dist<-list(train=dat.train.dist, all=dat.all.dist)
    return(dat.dist)
  }

  if(cv==FALSE){
    dat.dist<-getUnionDist(rnames, dat.wt, type)
    return(dat.dist)
  }

}



#take average
combineDist<-function(dist.dat){
  arr <- array(unlist(dist.dat), dim = c(nrow(dist.dat[[1]]),ncol(dist.dat[[1]]),length(dist.dat)))
  #calculate mean distance after removing NA
  combMat <- rowMeans(arr, dim = 2, na.rm = TRUE)
  rownames(combMat) <- rownames(dist.dat[[1]])
  colnames(combMat) <- colnames(dist.dat[[1]])
  #copy the rest of the triangle
  combMat[lower.tri(combMat)] <- t(combMat)[lower.tri(combMat)]
  colna = colSums(is.na(combMat))
  tab.idx = sort(table(unlist(apply(combMat,2,function(x) which(is.na(x))))),decreasing=T)
  idx = as.numeric(names(tab.idx))
  iter=1
  combMatDel = combMat
  #remove incomplete pairwise information
  while(sum(is.na(combMatDel))!=0){
    if (iter==1){del.idx = idx[iter]}
    if (iter!=1){del.idx = c(del.idx, idx[iter])}
    combMatDel = combMat[-del.idx, -del.idx]
    iter = iter + 1
    }

  if(nrow(combMatDel) < 2){stop("only one sample left after overlaping for complete pairs, low sample overlap!")}
  return(combMatDel)
}


##########################
# survClust
###########################

survclust<-function(combine.dist,out,k, cmd.k=NULL){
  if(is.null(names(out)))
    stop("rowanmes of out can't be NULL")

  if(is.null(rownames(combine.dist)))
    stop("rowanmes of combine.dist can't be NULL")

  if (!requireNamespace("pdist", quietly = TRUE)) {
    stop("pdist package needed for this function to work. Please install it.",
         call. = FALSE)
  }

  inter <- intersect(names(out), rownames(combine.dist))

  #run cmdscale
  combine.dist <- combine.dist[inter,inter]

  if(is.null(cmd.k)){cmd.k =nrow(combine.dist)-1 }
  if(!(is.null(cmd.k))){cmd.k =as.numeric(cmd.k) }

  cmd.combine.dist<-cmdscale(combine.dist,k=cmd.k)
  my.k = as.numeric(k)
  #run kmeans with 100 starts
  fit<-kmeans(cmd.combine.dist,my.k,nstart=100)

  #return fit 
  return(fit)
}

#predict test labels on survclust fitted
predict.test.label<-function(all.cmd,fit,k){
  all.cmd = as.matrix(all.cmd)
  train.snames = names(fit$cluster)
  test.snames = setdiff(rownames(all.cmd),train.snames)

  #where row - samples, col - genes
  centroid = matrix(NA, nrow = k, ncol = ncol(all.cmd))
  for (kk in 1:k) {
    #meaning k clust has one sample. #WARNING #check
    if(is.vector(all.cmd[names(fit$cluster)[which(fit$cluster==kk)],]) & ncol(all.cmd) > 1){
      message(paste0("k=",k, " training cluster has one sample, prediction might be inaccurate"))
      centroid[kk, ]=all.cmd[names(fit$cluster)[which(fit$cluster==kk)], ]
    }

    if (!(is.null(dim(all.cmd[names(fit$cluster)[fit$cluster==kk],])))){
      if(ncol(all.cmd)> 1){centroid[kk, ]=apply(all.cmd[names(fit$cluster)[which(fit$cluster==kk)], ], 2, mean)}
    }

    if(ncol(all.cmd)==1){centroid[kk,] = mean(all.cmd[names(fit$cluster)[which(fit$cluster==kk)], ])}
  }

  dist.whole = apply(centroid,1,function(x) as.matrix(pdist(x,all.cmd)))

  #assign the cluster membership
  dist.labels = apply(dist.whole,1,which.min)
  names(dist.labels) = rownames(all.cmd)
  test.labels = dist.labels[test.snames]
  
  #is missing a class label via pdist
  if(length(unique(test.labels)) != k){
    message(paste0("k=", k, " was reduced to ", length(unique(test.labels)), " in test label prediction"))}
    #do we really need to relabel? 
    #ul<-unique(test.labels)
    #nl<-1:length(ul)
    #tt<-rep(NA, length(test.labels))
    #names(tt) = names(test.labels)
    #for(i in 1:length(ul)){
     # tt[which(test.labels==ul[i])] = nl[i]
    #}
    #test.labels = tt}

  return(list(test.labels = test.labels))
}

#Runs cross validation on survclust to attain best k
#calculate sum of squares of each simulated dataset
do.ss.stats<-function(mm,labels){
  ll = unique(labels)
  mm[lower.tri(mm,diag=TRUE)]<-NA
  tss = sum(mm, na.rm=T)
  wss<-rep(NA, length(ll))
  for (i in 1:length(ll)){
    wss[i] = sum(mm[labels==ll[i], labels==ll[i]], na.rm=T)
  }

  tot.withinss = sum(wss)
  within.over.tot = tot.withinss/tss
  return(within.over.tot)
}

get.centroid<-function (mat, labels, f)
{
  ul <- unique(labels)
  if (ncol(mat) == 1) {
    warning("cmd reduces matrix to one eigen value! Noisy data?")
  }
  centroids <- matrix(NA, nrow = length(ul), ncol = ncol(mat))
  for (i in 1:length(ul)) {
    mat.ul <- mat[names(labels)[which(labels == ul[i])], ]
    if (is.vector(mat.ul)) {
      centroids[i, ] = mat.ul
    }
    else {
      centroids[i, ] <- apply(mat.ul, 2, mean)
    }
  }
  rownames(centroids) = paste0("f", f, "_k", ul)
  return(centroids)
}

get.relabel<-function(pattern, olabel, centroid.cluster,kk){
  relabel<-rep(NA, length(olabel))
  names(relabel) = names(olabel)

  for(i in 1:kk){
    kpattern<-paste0(pattern,i)
    idx<-which(names(centroid.cluster)==kpattern)
    #change current label to this
    if(length(idx)!=0){
      change.label = centroid.cluster[idx]
      idx2 = which(olabel == i)
      if(length(idx2)!=0){relabel[idx2] = change.label}
    }
  }
  if(any(is.na(relabel))){warning("there is a NA in relabel, something is wrong with clustering, noisy data or pick lower k?")}
  return(relabel)
}

#this is done to provide meaning to cluster labels across folds.
#we run kmeans on centroid vector of all the folds to determine their closeness.
#random start -10
cv.relabel<-function(mat, train.labels,cv.test.labels, k,fold){

  centroids<-list()

  for(i in 1:length(train.labels)){
    centroids[[i]]<-get.centroid(mat, train.labels[[i]],i)
  }

  centroids.all<-do.call(rbind.data.frame, lapply(centroids, function(x) x))
  #do kmeans on the centroids
  centroids.kmeans<-kmeans(centroids.all,k,nstart=20)
  #print(centroids.kmeans$cluster)
  #centroids cluster labels
  all.cluster<-centroids.kmeans$cluster

  relabel<-rep(NA,nrow(mat))
  names(relabel) = rownames(mat)

  for(i in 1:fold){
    pattern = paste0("f",i,"_k")
    rr<-get.relabel(pattern, cv.test.labels[[i]], all.cluster,k)
    relabel[names(rr)] = rr
  }

  return(relabel)
}

#############################
# perform cross validation 
#############################

cv.survclust<-function(x, out,k,fold, cmd.k=NULL, type=NULL){

  my.k <- as.numeric(k)
  fold <- as.numeric(fold)

  #To get an idea of total samples
  dist.dat<-getDist(x,out, type=type)
  #this for calculating ss on test labels
  combine.dist<-combineDist(dist.dat)
  inter <- intersect(names(out), rownames(combine.dist))
  out = out[inter]
  ll <- seq(1,length(inter))
  
  #we will use this for cluster relabeling of test
  if(is.null(cmd.k)){this.k = nrow(combine.dist)-1}
  if(!(is.null(cmd.k))){this.k = as.numeric(cmd.k)}
  
  combine.dist.cmd<-cmdscale(combine.dist, k=this.k)
  
  folds <- cut(seq(1,length(ll)),breaks=fold,labels=FALSE)
  ll.rand<- sample(ll,length(ll))

  #Perform n fold cross validation
  cv.test.labels<-list()
  #cv.test.rand.index =NA
  survfit<-list()

  for(i in 1:fold){
    #Segement your data by fold using the which() function
    test.idx <- ll.rand[which(folds==i)]
    train.idx <- setdiff(ll,test.idx)
    train.snames = names(out)[train.idx]
    out.train<- out[train.snames]

    #multiply by coxph abs(log(HR))
    distwt<-getDist(x,out,cv=TRUE,train.snames, type=type)
    train.dist.dat<-distwt$train
    all.dist.dat<-distwt$all

    #combine dist
    train.combine.dist<-combineDist(train.dist.dat)
    all.combine.dist<-combineDist(all.dist.dat)
    inter <- intersect(names(out), rownames(all.combine.dist))
    all.combine.dist = all.combine.dist[inter,inter]
    
    if(is.null(cmd.k)){cmd.k.all =nrow(all.combine.dist)-1 }
    if(!(is.null(cmd.k))){cmd.k.all =as.numeric(cmd.k) }
    #as cmd is on dist, and nrow is different for all and training set
    #but multiplied with training HR
    cmd.whole = cmdscale(all.combine.dist,k=cmd.k.all)
    #get training fit labels
    fit=survclust(train.combine.dist,out,my.k,cmd.k)

    #calculate test logrank and concordance
    #we basically predict on whole
    test =predict.test.label(cmd.whole,fit,my.k)
    cv.test.labels[[i]] = test$test.labels
    survfit[[i]]<-fit
  }

  train.labels = lapply(survfit, function(x) x$cluster)
  cv.test.relabels<-cv.relabel(combine.dist.cmd, train.labels,cv.test.labels,my.k,fold)
  min.labels = min(table(cv.test.relabels))
  idx = which(min.labels <=5)
  if (length(idx)!=0){message(paste0("k= ", my.k, " has 5 or few samples in cluster solution"))}

  message(paste0("finished ", fold, " cross validation, total samples-", length(cv.test.relabels)))
  cv.test.relabels = cv.test.relabels[names(out)]
  if (length(unique(cv.test.relabels)) != my.k){warning(paste0("Test labels not equal to chosen k ",my.k)) }

  #if everything collapses after test relabeling
  if (length(unique(cv.test.relabels)) ==1){
                      cv.all.logrank = NA
                      warning("Everything collapsed after test relabel, logrank test is NA")
  }

  cv.test.ss<-do.ss.stats(combine.dist, cv.test.relabels)
  cv.fit = list(cv.labels = cv.test.relabels, cv.ss = cv.test.ss)
  return(cv.fit)

}

cv.voting<-function(fit,dd,kk, cmd.k=NULL, kk.test=TRUE, minlabel.test=TRUE){
  
  cc=combineDist(dd)
  if(is.null(cmd.k)){cmd.mat = cmdscale(cc, nrow(cc)-1)}
  if(!is.null(cmd.k)){cmd.mat = cmdscale(cc, cmd.k)}
  
  fit = fit[[kk-1]] 
  if(nrow(cmd.mat) != length(fit[[1]]$cv.labels)){stop("unequal samples in CV and cmd mat")}
  #remove fits with solutions != k
  ttt =unlist(lapply(fit, function(x) length(unique(x$cv.labels)) ))
  idx = which(ttt < kk)
  
  #test for solutions less than 5
  min.labels = lapply(fit, function(x) min(table(x$cv.labels)))
  idx2 = which(min.labels < 5)
  
  if(length(idx)==0){ kk.test=FALSE}
  if(length(idx2)==0){minlabel.test=FALSE}
  
  if(kk.test==TRUE & minlabel.test==TRUE){idx = unique(c(idx,idx2)); fit = fit[-idx]}
  if(kk.test==TRUE & minlabel.test==FALSE){fit = fit[-idx]}
  if(kk.test==FALSE & minlabel.test==TRUE){fit = fit[-idx2]}
  if(kk.test==FALSE & minlabel.test==FALSE){fit = fit}
  
  message(paste0("performing consensus on ", length(fit), " rounds"))
  if(length(fit)==0){stop("cross validation returned labels not equal to kk, pick another kk OR relax filters")}
  cv.rounds = length(fit)
  
  cmd.mat = cmd.mat[names(fit[[1]]$cv.labels),]

  centroids<-list()
  for (i in 1:cv.rounds){
    centroids[[i]]<-get.centroid(cmd.mat, fit[[i]]$cv.labels,i)
  }
  
  centroids.all<-do.call(rbind.data.frame, lapply(centroids, function(x) x))
  #do kmeans on the centroids
  centroids.kmeans<-kmeans(centroids.all,kk,nstart=100)
  all.cluster<-centroids.kmeans$cluster
  #print(table(all.cluster))
  #print(all.cluster)
  relabels<-list()
  for(i in 1:cv.rounds){
    pattern = paste0("f",i,"_k")
    relabels[[i]]<-get.relabel(pattern, fit[[i]]$cv.labels, all.cluster,kk)
    
  }
  relabels.all<-do.call(rbind.data.frame, lapply(relabels, function(x) x))
  relabels.all = apply(relabels.all, 2, as.numeric)
  colnames(relabels.all) = names(fit[[1]]$cv.labels)
  final.labels<-apply(relabels.all,2,function(x) names(table(x))[which.max(table(x))])
  return(unlist(final.labels))
}

