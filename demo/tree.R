## this uses the opposite of the standard row/col convention... to be fixed
build.kdtree <- function(p, nbucket = 4*nrow(p)){
  nd <- nrow(p)
  n <- ncol(p)

  rbb <- c(paste(c('bbl','bbu'),rep(1:nd, each = 2),sep = '.'))
  rows <- rbb
  ntrd <- length(rows)

  rchi <- c('child.l', 'child.u')
  rowsi <- c(rchi,'nlev','istart')

  ntri <- length(rowsi)

  ##numnode <- 2^(1+ceiling(log(n/nbucket)/log(2))) - 1
  nf <- 2^(1+floor(log(n/nbucket)/log(2)))
  nmax <- 2^(1+ceiling(log(n/nbucket)/log(2))) - 1
  numnode <- min(2*n - (nbucket-1)*nf - 1, nmax)
  
  trd <- matrix(nrow = ntrd, ncol = numnode, dimnames = list(rows, NULL))
  tri <- matrix(data = -1, nrow = ntri, ncol = numnode, dimnames = list(rowsi, NULL))

## a basic kd tree node keeps track of:
## 1. box centre
## 2. edge lengths
## 3. children
## 4. points
##   - if we use a selection function which finds the median and
##     puts points in the array so that they are high-low sorted
##     then every node knows how to find all leaves contained by it
##     quite simply
##     -> sort(x, partial = )

## ROOT STARTS AT LEVEL 1, NODE INDEX 1
## node: current node 
## lev:  level in tree (somewhat redundant)
## nlev: number of poins at this level
## istart: starting index in point array
## nbucket: maximum number of points in a leaf
## nn: total number of nodes in your subtree
  
  build.tree <- function(node = 1, lev = 1, nlev = n, istart = 1, nbucket = nbucket){
    ip <- istart + 0:(nlev-1)
    
    trd[rows,node] <<- sapply(1:nd, function(i) range(p[i,ip]))

    ## set some accounting info

    tri['nlev', node] <<- nlev
    tri['istart', node] <<- istart

    ni <- 1

    if(nlev > nbucket){
      ## next we split the node, cycling through dimensions equally
      d.split <- (lev-1)%%nd + 1

      ## This is more efficiently done in the c backend
      med <- median(p[d.split, ip])
      i.lo <- which(p[d.split,ip] <= med)
      
      p[,ip] <<- cbind(p[,ip[i.lo],drop=FALSE],p[,ip[-i.lo],drop=FALSE])

      cnl2 <- ceiling(nlev/2)

      tri['child.l', node] <<- node + ni

      ni <- ni + build.tree(node = tri['child.l',node], lev = lev+1,
                            nlev = cnl2, istart = istart, nbucket = nbucket)

      tri['child.u', node] <<- node + ni
      
      ni <- ni + build.tree(node = tri['child.u',node], lev = lev+1,
                            nlev = nlev - cnl2, istart = istart + cnl2, nbucket = nbucket)
    }
    return(ni)
  }

  build.tree(nbucket=nbucket)

  tree <- list(p = p, trd = trd, tri = tri, nbucket = nbucket, numnode = numnode, rbb = rbb)
  invisible(tree)
}

## Two problems need be addressed:
## a) orientation of the two rectangles wrt one another
## b) nature of intersection
## 1. No overlap.
## 2. Search box partially contained in node 
## 3. Search box totally contained in node
## 4. Search box partially contains node.
## 5. Search box totally contains node.
## 6. Intersection of search boxes with no contained corners.

## This can be reduced to three states:
## miss, hit-open (partial intersection), and hit-done (the search box envelopes the node)

KD.NOCHILD <- -1

KD.MISS <- 0
KD.HITOPEN <- 1
KD.HITDONE <- 2

## This might look a little mysterious, but it's the result of manually mapping out
## the 6 different scenarios above, so it should look mysterious
boxIntersect <- function(bba, bbb){
  ## bounding boxes - lower.1, upper.1, lower.2, upper.2, ...
  pa <- matrix(bba, byrow = TRUE, ncol = 2)
  pb <- matrix(bbb, byrow = TRUE, ncol = 2)

  ## not bothering to be tricky here - just do two tables of relevant comparisons
  tta <- pa[,1] < pb
  ttb <- pa[,2] < pb

  tt <- cbind(tta,ttb)

  miss <- any(sapply(1:nrow(tt), function(i){ (sum(tt[i,])%%4) == 0 }))

  if(miss)
    return(KD.MISS)

  done <- all(sapply(1:nrow(tt), function(i){ all(xor(tt[i,c(1,2)],tt[i,c(3,4)])) }))
  return(ifelse(done, KD.HITDONE, KD.HITOPEN))
}

boxSearch <- function(tree, box, node = 1){
  ihit <- node

  if((res <- boxIntersect(box, tree$trd[tree$rbb,node])) == KD.MISS)
    return(integer(0))

  is.leaf <- tree$tri['child.l',node] == KD.NOCHILD
  
  if((res == KD.HITOPEN) && (!is.leaf))
    ihit <- c(boxSearch(tree, box, node = tree$tri['child.l',node]),
              boxSearch(tree, box, node = tree$tri['child.u',node]))

  return(ihit)
}

## playing around with the tree 
if(FALSE){
  set.seed(42)

  n <- 67
  nd <- 2
  nbucket <- 8
  
  q <- matrix(runif(nd*n), nrow = nd, ncol = n)
  
  mytree <- build.kdtree(q, nbucket = nbucket)

  attach(mytree)

  op <- par(mfrow = c(2, 2), pty = "s")

  for(i in 1:4){
    plot(p[1,],p[2,],pch=19, type = 'n', xlab = 'x', ylab = 'y', main = paste("level",i))
    ilev <- round(log(n / tri["nlev",])/log(2))+1
    lev.plot <- c(i)
    i.plot <- ilev %in% lev.plot
    colors <- c("blue","green","orange", "red", "purple")
    rect(trd["bbl.1",i.plot], trd["bbl.2",i.plot], trd["bbu.1",i.plot], trd["bbu.2",i.plot],lwd = 4, border = (colors[ilev])[i.plot])
    points(p[1,],p[2,],pch=19)
  }
  
  detach(mytree)

  ## do a box search
  x11()
  plot(q[1,],q[2,], pch = 19, xlab ='x', ylab='y')

  box <- c(0.25,0.5,0.25,0.5)
  nb <- boxSearch(mytree, box)

  rect(mytree$trd[1,nb],mytree$trd[3,nb],mytree$trd[2,nb],mytree$trd[4,nb], border = 'blue', lwd = 3)
  rect(box[1],box[3],box[2],box[4], border = 'red', lwd = 3)
  
  if(FALSE){
    set.seed(42)
    x11()
    while(readline(prompt='continue?') == 'y'){
      bba <- c(range(runif(2)),range(runif(2)))
      bbb <- c(range(runif(2)),range(runif(2)))
      
      res <- boxIntersect(bba,bbb)
      
      plot(NA, main = res, xlab = 'x', ylab = 'y', xlim = c(0,1), ylim = c(0,1), type = 'n')
      rect(bba[1],bba[3],bba[2],bba[4], border = 'red')
      rect(bbb[1],bbb[3],bbb[2],bbb[4], border = 'green')
    }
  }

}
