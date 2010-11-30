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

void build_kdtree(double ** p, int nump, int ndim, int nbucket, int ** ip, KDT ** kdt);
void kdSelect(double ** p, KDT * kdt, int * ip, int d, int k, int l, int r);
int build_tree(double ** p, KDT * kdt, int * ip, int node, int d, int nlev, int istart);
void free_kdtree(KDT ** kdt);
int boxIntersect(double * bbs, double * bbb, int ndim);
void boxSearch(KDT * kdt, int node, double * bb, NL * nl);
