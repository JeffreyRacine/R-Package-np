#ifndef KDT_INCLUDED
#define KDT_INCLUDED

#define KD_NOCHILD -1

#define KD_MISS 0
#define KD_HITOPEN 1
#define KD_HITDONE 2

//  bb = bounding box, in form: lower.dim1, upper.dim1, lower.dim2 ... and so on

struct kdn {
  double * bb;
  int childl, childu, nlev, istart;
};

typedef struct kdn KDN;

struct kdt {
  KDN * kdn;
  double * bb;
  int ndim, nbucket, numnode;
};

typedef struct kdt KDT;

struct nl {
  int * node, n, nalloc;
};

typedef struct nl NL;

struct xl {
  int * istart, * nlev, n, nalloc;
};

typedef struct xl XL;

void build_kdtree(double ** p, int nump, int ndim, int nbucket, int * ip, KDT ** kdt);
void kdSelect(double ** p, KDT * kdt, int * ip, int d, int k, int l, int r);
int build_tree(double ** p, KDT * kdt, int * ip, int node, int d, int nlev, int istart);
void free_kdtree(KDT ** kdt);
int boxIntersect(double * bbs, double * bbb, int ndim);
void boxSearch(KDT * kdt, int node, double * bb, NL * nl);
void check_grow_nl(NL * nl);
void boxSearchNL(KDT * restrict kdt, NL * restrict search, double * restrict bb, NL * restrict nl, XL * restrict xl);
void clean_nl(NL * restrict nl);
void clean_xl(XL * restrict xl);
void mirror_nl(NL * restrict nla, NL * restrict nlb);
void mirror_xl(XL * restrict xla, XL * restrict xlb);
int boxIntersectPartial(double * bbs, double * bbb, int * restrict idim, int nidim);
void boxSearchNLPartial(KDT * restrict kdt, NL * restrict search, double * restrict bb, NL * restrict nl, XL * restrict xl, int * idim, int nidim);
void boxSearchNLPartialIdx(KDT * restrict kdt, NL * restrict search, double * restrict bb, NL * restrict nl, XL * restrict xl, int * idim, int nidim, int * idx);

void merge_end_xl(XL * restrict xl, KDN * restrict kdn);
void merge_end_xl_idx(XL * restrict xl, KDN * restrict kdn, int * restrict idx);
#endif
