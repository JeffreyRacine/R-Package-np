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

void build_kdtree(double ** p, int nump, int ndim, int nbucket, int * ip, KDT ** kdt){
  KDT * kdx;
  int nodecount;
  int nf = (1 << (1 + (int)floor(log(MAX((double)nump,(double)nbucket)/(double)nbucket)/log(2.0))));
  int maxnode = (1 << (1 + (int)ceil(log(MAX((double)nump,(double)nbucket)/(double)nbucket)/log(2.0)))) - 1; 
  
  int numnode = MIN(2*MAX(nump,nbucket) - (nbucket-1)*nf - 1, maxnode);

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

  nodecount = build_tree(p, kdx, ip, 0, 0, nump, 0);

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

int boxIntersectPartial(double * bbs, double * bbb, int * restrict idim, int nidim){
  int i, tt[4], miss;
  int hdone = 1;

  // here we build up a 'truth table' of sorts
  // at the same time we test for a total miss and break
  // away early if that's the case
  for(i = 0; i < nidim; i++){
    miss  = tt[0] = bbs[2*idim[i]] < bbb[2*idim[i]];
    miss += tt[1] = bbs[2*idim[i]] < bbb[2*idim[i]+1];
    miss += tt[2] = bbs[2*idim[i]+1] < bbb[2*idim[i]];
    miss += tt[3] = bbs[2*idim[i]+1] < bbb[2*idim[i]+1];
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
    nl->node = (int *)realloc(nl->node, MAX(10, 2*nl->nalloc)*sizeof(int));
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

// purely iterative version
void boxSearchNL(KDT * restrict kdt, NL * restrict search, double * restrict bb, NL * restrict nl, XL * restrict xl){
  NL nls = {.node = NULL, .n = 0, .nalloc = 0};
  mirror_nl(search, &nls);

  while (nls.n > 0){
    const int node = nls.node[nls.n - 1];

    int res = boxIntersect(bb, kdt->kdn[node].bb, kdt->ndim);

    if(res == KD_MISS) {
      nls.n--;
      continue;
    }

    if((res == KD_HITDONE) || (kdt->kdn[node].childl == KD_NOCHILD)){
      if(nl != NULL){
        check_grow_nl(nl);
        nl->node[nl->n++] = node;
      }

      if(xl != NULL){
        merge_end_xl(xl, &kdt->kdn[node]);
      }

      nls.n--;
    }
    else { // KD_HITOPEN
      check_grow_nl(&nls);
      nls.node[nls.n - 1] = kdt->kdn[node].childu;
      nls.node[nls.n++] = kdt->kdn[node].childl;
    }
  }
  clean_nl(&nls);
}

// purely iterative version
void boxSearchNLPartial(KDT * restrict kdt, NL * restrict search, double * restrict bb, NL * restrict nl, XL * restrict xl, int * idim, int nidim){
  NL nls = {.node = NULL, .n = 0, .nalloc = 0};
  mirror_nl(search, &nls);

  while (nls.n > 0){
    const int node = nls.node[nls.n - 1];

    int res;

    if(nidim > 0)
      res = boxIntersectPartial(bb, kdt->kdn[node].bb, idim, nidim);
    else
      res = KD_HITDONE;

    if(res == KD_MISS) {
      nls.n--;
      continue;
    }

    if((res == KD_HITDONE) || (kdt->kdn[node].childl == KD_NOCHILD)){
      if(nl != NULL){
        check_grow_nl(nl);
        nl->node[nl->n++] = node;
      }

      if(xl != NULL){
        merge_end_xl(xl, &kdt->kdn[node]);
      }

      nls.n--;
    }
    else { // KD_HITOPEN
      check_grow_nl(&nls);
      nls.node[nls.n - 1] = kdt->kdn[node].childu;
      nls.node[nls.n++] = kdt->kdn[node].childl;
    }
  }
  clean_nl(&nls);
}

// index-restricted (inclusive) partial search
// NB - 
// 1) use create_fake_nodes to ensure that all indices are relative to idx[0], and are clamped to the range idx[0]..idx[1]
// 2) once you are done with your fake results, call reset_fake_nodes
// 3) call create_fake_nodes before reuse if reset_fake_nodes has been called in the meanwhile

void boxSearchNLPartialIdx(KDT * restrict kdt, NL * restrict search, double * restrict bb, NL * restrict nl, XL * restrict xl, int * idim, int nidim, int * idx){
  int res, tt[4];

  NL nls = {.node = NULL, .n = 0, .nalloc = 0};
  mirror_nl(search, &nls);

  while (nls.n > 0){
    const int node = nls.node[nls.n - 1];
    
    const int il = kdt->kdn[node].istart;
    const int ih = kdt->kdn[node].istart+kdt->kdn[node].nlev-1;

    tt[0] = (il >= idx[0]);
    tt[1] = (il > idx[1]);

    tt[2] = (ih < idx[0]);
    tt[3] = (ih <= idx[1]);

    if(!((tt[0] == tt[1]) && (tt[2] == tt[3]) && (tt[0] != tt[2])))
      if(nidim > 0)
        res = boxIntersectPartial(bb, kdt->kdn[node].bb, idim, nidim);
      else
        res = KD_HITDONE;
    else
      res = KD_MISS;

    if(res == KD_MISS) {
      nls.n--;
      continue;
    }

    if((res == KD_HITDONE) || (kdt->kdn[node].childl == KD_NOCHILD)){
      if(nl != NULL){
        check_grow_nl(nl);
        nl->node[nl->n++] = node;
      }

      if(xl != NULL){
        merge_end_xl_idx(xl, &kdt->kdn[node], idx);
      }

      nls.n--;
    }
    else { // KD_HITOPEN
      check_grow_nl(&nls);
      nls.node[nls.n - 1] = kdt->kdn[node].childu;
      nls.node[nls.n++] = kdt->kdn[node].childl;
    }
  }
  clean_nl(&nls);
}

void check_grow_nl(NL * restrict nl){
  if(nl->n == nl->nalloc){
    nl->node = (int *)realloc(nl->node, MAX(10, 2*nl->nalloc)*sizeof(int));
    if(!(nl->node != NULL))
      error("!(nl->node != NULL)");
    
    nl->nalloc = MAX(10, 2*nl->nalloc);
  }
}

void clean_nl(NL * restrict nl){
  if(nl){
    if(nl->node != NULL)
      free(nl->node);
    nl->node = NULL;
    nl->n = nl->nalloc = 0;    
  }
}

void clean_xl(XL * restrict xl){
  if(xl){
    if(xl->istart != NULL)
      free(xl->istart);

    if(xl->nlev != NULL)
      free(xl->nlev);

    xl->nlev = xl->istart = NULL;
    xl->n = xl->nalloc = 0;    
  }
}

void mirror_nl(NL * restrict nla, NL * restrict nlb){
  if(nla->n > nlb->nalloc){
    nlb->node = (int *)realloc(nlb->node, (1+nla->n)*sizeof(int));
    nlb->nalloc = nla->n + 1;
  }

  for(int i = 0; i < nla->n; i++)
    nlb->node[i] = nla->node[i];
  
  nlb->n = nla->n; 
}

void mirror_xl(XL * restrict xla, XL * restrict xlb){
  if(xla->n > xlb->nalloc){
    xlb->istart = (int *)realloc(xlb->istart, (1+xla->n)*sizeof(int));
    xlb->nlev = (int *)realloc(xlb->nlev, (1+xla->n)*sizeof(int));
    xlb->nalloc = xla->n + 1;
  }

  for(int i = 0; i < xla->n; i++){
    xlb->istart[i] = xla->istart[i];
    xlb->nlev[i] = xla->nlev[i];
  }
  
  xlb->n = xla->n; 
}

void merge_end_xl(XL * restrict xl, KDN * restrict kdn){
  if(xl->n == xl->nalloc){
    xl->istart = (int *)realloc(xl->istart, MAX(10,2*xl->nalloc)*sizeof(int));
    xl->nlev = (int *)realloc(xl->nlev, MAX(10,2*xl->nalloc)*sizeof(int));

    xl->nalloc = MAX(10, 2*xl->nalloc);
  }

  if((xl->n != 0) && (kdn->istart == (xl->istart[xl->n-1] + xl->nlev[xl->n-1]))){
    xl->nlev[xl->n-1] += kdn->nlev;
  } else {
    xl->istart[xl->n] = kdn->istart;
    xl->nlev[xl->n++] = kdn->nlev;
  }
}

void merge_end_xl_idx(XL * restrict xl, KDN * restrict kdn, int * restrict idx){
  if(xl->n == xl->nalloc){
    xl->istart = (int *)realloc(xl->istart, MAX(10,2*xl->nalloc)*sizeof(int));
    xl->nlev = (int *)realloc(xl->nlev, MAX(10,2*xl->nalloc)*sizeof(int));

    xl->nalloc = MAX(10, 2*xl->nalloc);
  }
  
  const int xln = xl->n;
  if((xln != 0) && ((kdn->istart - idx[0]) == (xl->istart[xln-1] + xl->nlev[xln-1]))){
    xl->nlev[xln-1] += kdn->nlev;
    xl->nlev[xln-1] = MIN(idx[1]-idx[0], xl->istart[xln-1] + xl->nlev[xln-1] - 1) - xl->istart[xln-1] + 1;
  } else {
    xl->istart[xln] = MAX(0, kdn->istart-idx[0]);
    xl->nlev[xln] = MIN(idx[1], kdn->istart + kdn->nlev - 1) - idx[0] - xl->istart[xln] + 1;
    xl->n++;
  }
}

void free_kdtree(KDT ** kdt){
  KDT * kdx = *kdt;
  
  free(kdx->kdn);
  free(kdx->bb);
  free(kdx);
  *kdt = NULL;
}
