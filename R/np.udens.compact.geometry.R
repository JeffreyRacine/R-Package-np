# Original physical geometry is used only for observed uniform support masks.
# Continuous quadrature keeps its normalized integration coordinates.
.np_udens_compact_geometry <- function(state,grid) {
  if(!identical(state$kernel,'uniform'))return(NULL)
  x<-as.matrix(state$x);grid<-as.matrix(grid);h<-state[['physical.h',exact=TRUE]]
  if(!is.numeric(h)||length(h)!=ncol(x)||any(!is.finite(h))||any(h<=0)||
     ncol(grid)!=ncol(x)||any(!is.finite(c(x,grid))))
    stop('invalid compact physical observation geometry',call.=FALSE)
  list(x=x,grid=grid,h=h)
}

.np_udens_compact_geometry_slice <- function(physical,rows=NULL,columns=NULL,observations=NULL) {
  if(is.null(physical))return(NULL)
  if(!is.null(rows))physical$grid<-physical$grid[rows,,drop=FALSE]
  if(!is.null(observations))physical$x<-physical$x[observations,,drop=FALSE]
  if(!is.null(columns)) {
    physical$grid<-physical$grid[,columns,drop=FALSE]
    physical$x<-physical$x[,columns,drop=FALSE];physical$h<-physical$h[columns]
  }
  physical
}

# Distinct physical targets must not be collapsed by normalized rounding.
.np_udens_compact_physical_grid <- function(grid,physical) {
  raw<-physical$grid
  keep<-which(!duplicated(as.data.frame(raw)))
  keep<-keep[do.call(order,as.data.frame(raw[keep,,drop=FALSE]))]
  key<-function(z)apply(z,1L,function(v)paste(format(v,digits=17,scientific=TRUE),collapse=':'))
  list(grid=if(is.matrix(grid))grid[keep,,drop=FALSE] else grid[keep],
    physical=.np_udens_compact_geometry_slice(physical,rows=keep),
    restore=match(key(raw),key(raw[keep,,drop=FALSE])))
}

# This small factory forces only observation inputs, never a quadrature frame.
.np_udens_compact_observed_axis <- function(x,h,grid,kernel,physical) {
  force(x);force(h);force(grid);force(kernel);force(physical)
  function(rows) {
    da<-list(list(lv=log(x[,1]),l1v=log1p(-x[,1])),
      list(lv=log(x[,2]),l1v=log1p(-x[,2])))
    .np_udens_compact_normal_q2_axis(kernel,
      .np_udens_compact_geometry_slice(physical,rows=rows))(grid[rows,,drop=FALSE],da,h)
  }
}
