#include <math.h>

#include <R.h>
#include "headers.h"
#include "tree.h"

/*
  inputs
  p: point coordinate array [dims][obs]
  nump: number of points
  ndim: dimensionality of the data
  nbucket: maximum number of particles contained in a leaf - this will likely require some performance tuning.
    the minimum is ndim for non-degenerate leaf shapes. At all levels in the tree, a bounding box is computed
    around each partition of the data. That means between ndim and 2*ndim points lie on the bounding box itself.
    At the moment I think a sensible minimum is 4*ndim, but it's just a gut feeling. 
  
  outputs
  ip: an index array for putting the data in tree order
      the array is of the type array[new_index] = orig_index
      all training data needs to be reordered according to this scheme.
  kdt: the tree structure

*/

void build_kdtree(double ** p, int nump, int ndim, int nbucket, int ** ip, KDT ** kdt){
  KDT * kdx;
  int nodecount;
  int nf = (1 << (1 + (int)floor(log((double)nump/(double)nbucket)/log(2.0))));
  int maxnode = (1 << (1 + (int)ceil(log((double)nump/(double)nbucket)/log(2.0)))) - 1; 
  
  int numnode = MIN(2*nump - (nbucket-1)*nf - 1, maxnode);

  *kdt = (KDT *)malloc(sizeof(KDT));
  if(!(*kdt != NULL))
    error("!(*kdt != NULL)");

  kdx = *kdt;

  kdx->kdn = (KDN *)malloc(numnode*sizeof(KDN));
  if(!(kdx->kdn != NULL))
    error("!(kdx->kdn != NULL)");

  kdx->bb = (double *)malloc(numnode*2*ndim*sizeof(double));
  if(!(kdx->bb != NULL))
    error("!(kdx->bb != NULL)");

  for(int i = 0; i < numnode; i++){
    kdx->kdn[i].bb = kdx->bb+2*i*ndim;

    kdx->kdn[i].childl = -1;
    kdx->kdn[i].childu = -1;
  }

  kdx->numnode = numnode;
  kdx->nbucket = nbucket;
  kdx->ndim = ndim;

  // tree is now fully allocated and ready to be built

  // for the tree to work, points need to be put into tree order
  *ip = (int *)malloc(nump*sizeof(int));
  if(!(*ip != NULL))
    error("!(*ip != NULL)");

  for(int i = 0; i < nump; i++){
    (*ip)[i] = i;
  }

  nodecount = build_tree(p, kdx, *ip, 0, 0, nump, 0);

  if(!(nodecount == numnode))
    error("!(nodecount == numnode)");
}

void kdSelect(double ** p, KDT * kdt, int * ip, int d, int k, int l, int r){
  int ti;
  int i,j;

  //  int tj, tk;

  double v;

  while (r > l) {
    v = p[d][ip[k]];

    ti = ip[r];
    ip[r] = ip[k];
    ip[k] = ti;

    i = l - 1;
    j = r;
    while (1) {
      while (i < j) if (p[d][ip[++i]] >= v) break;
      while (i < j) if (p[d][ip[--j]] <= v) break;
      ti = ip[i];
      ip[i] = ip[j];
      ip[j] = ti;
      
      if (j <= i) break;
    }

    ip[j] = ip[i];
    ip[i] = ip[r];
    ip[r] = ti;
    if (i >= k) r = i - 1;
    if (i <= k) l = i + 1;
  }
}

int build_tree(double ** p, KDT * kdt, int * ip, int node, int d, int nlev, int istart){
  int i,j, ni = 1;
  KDN * const tn = &(kdt->kdn[node]);

  tn->nlev = nlev;
  tn->istart = istart;

  // find bounding box for node
  for(i = 0; i < kdt->ndim; i++){
    tn->bb[2*i] = DBL_MAX;
    tn->bb[2*i+1] = -DBL_MAX;
  }
    
  for(j = istart; j < (istart+nlev); j++){
    for(i = 0; i < kdt->ndim; i++){
      tn->bb[2*i] = MIN(tn->bb[2*i], p[i][ip[j]]);
      tn->bb[2*i+1] = MAX(tn->bb[2*i+1], p[i][ip[j]]);
    }
  }

  if(nlev > kdt->nbucket){
    const int k = istart + nlev/2 - 1;
    const int d_next = (d+1) % kdt->ndim;

    // split on dimension 'd'
    kdSelect(p, kdt, ip, d, k, istart, istart + nlev - 1);

    tn->childl = node+ni;
    ni += build_tree(p, kdt, ip, tn->childl, d_next, k - istart + 1, istart);

    tn->childu = node+ni;
    ni += build_tree(p, kdt, ip, tn->childu, d_next, nlev - nlev/2, k + 1);
  }
  return ni;
}

/*
  inputs
  bbs: the search box
  bbb: the test box
  ndim: dimensionality of the data
  
  outputs
  returns KD_MISS neither boxes touch or contain one another
  returns KD_HITOPEN if box edges intersect or the test box contains the search box
  returns KD_HITDONE if the search box contains the test box

*/

int boxIntersect(double * bbs, double * bbb, int ndim){
  int i, tt[4], miss;
  int hdone = 1;

  // here we build up a 'truth table' of sorts
  // at the same time we test for a total miss and break
  // away early if that's the case
  for(i = 0; i < ndim; i++){
    miss  = tt[0] = bbs[2*i] < bbb[2*i];
    miss += tt[1] = bbs[2*i] < bbb[2*i+1];
    miss += tt[2] = bbs[2*i+1] < bbb[2*i];
    miss += tt[3] = bbs[2*i+1] < bbb[2*i+1];
    if(miss%4 == 0) return(KD_MISS);

    hdone = hdone && (tt[0]^tt[2]) && (tt[1]^tt[3]);
  }
  return(hdone ? KD_HITDONE : KD_HITOPEN);
}

/*
  builds a list of nodes intersected by search box bb
  it uses the nl (node list) structure to aggregate results as the search progresses
*/
void boxSearch(KDT * kdt, int node, double * bb, NL * nl){
  int res = boxIntersect(bb, kdt->kdn[node].bb, kdt->ndim);

  if(res == KD_MISS) return;

  if(nl->n == nl->nalloc){
    nl->node = realloc(nl->node, MAX(10, 2*nl->nalloc)*sizeof(int));
    if(!(nl->node != NULL))
      error("!(nl->node != NULL)");
    
    nl->nalloc = MAX(10, 2*nl->nalloc);
  }

  if((res == KD_HITDONE) || (kdt->kdn[node].childl == KD_NOCHILD))
    nl->node[nl->n++] = node;
  else { // KD_HITOPEN
    boxSearch(kdt, kdt->kdn[node].childl, bb, nl);
    boxSearch(kdt, kdt->kdn[node].childu, bb, nl);
  }
}

void free_kdtree(KDT ** kdt){
  KDT * kdx = *kdt;
  
  free(kdx->kdn);
  free(kdx->bb);
  free(kdx);
  *kdt = NULL;
}
