/* Conditional-only ownership around the unchanged Numerical Recipes source.
 * The shared nr.c object and every regression caller remain unchanged.
 * Private names keep line-search state and RNG symbols separate. */
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
typedef union NPCSearchBlock NPCSearchBlock;
union NPCSearchBlock {
  max_align_t alignment;
  struct { NPCSearchBlock *previous, *next; } links;
};
typedef struct { NPCSearchBlock *first; } NPCSearchArena;
static NPCSearchArena *npc_search_arena;
static void *npc_search_allocate(size_t bytes) {
  if(!npc_search_arena || bytes>SIZE_MAX-sizeof(NPCSearchBlock))
    error("conditional search allocation context is invalid");
  NPCSearchBlock *block=malloc(sizeof(*block)+bytes);
  if(!block)return NULL;
  block->links.previous=NULL;block->links.next=npc_search_arena->first;
  if(block->links.next)block->links.next->links.previous=block;
  npc_search_arena->first=block;
  return block+1;
}
static void npc_search_release(void *pointer) {
  if(!pointer)return;
  NPCSearchBlock *block=(NPCSearchBlock*)pointer-1;
  if(block->links.previous)block->links.previous->links.next=block->links.next;
  else npc_search_arena->first=block->links.next;
  if(block->links.next)block->links.next->links.previous=block->links.previous;
  free(block);
}
#define sort npc_private_sort
#define powell npc_private_powell
#define nrerror npc_private_nrerror
#define vector npc_private_vector
#define free_vector npc_private_free_vector
#define linmin npc_private_linmin
#define brent npc_private_brent
#define mnbrak npc_private_mnbrak
#define f1dim npc_private_f1dim
#define erfun npc_private_erfun
#define np_reset_nr_rng_state npc_private_np_reset_nr_rng_state
#define ran3 npc_private_ran3
#define gasdev npc_private_gasdev
#define chidev npc_private_chidev
#define ncom npc_private_ncom
#define pcom npc_private_pcom
#define xicom npc_private_xicom
#define nrfunc npc_private_nrfunc
#define iff npc_private_iff
#define malloc npc_search_allocate
#define free npc_search_release
#include "nr.c"
#undef malloc
#undef free
#undef sort
#undef powell
#undef nrerror
#undef vector
#undef free_vector
#undef linmin
#undef brent
#undef mnbrak
#undef f1dim
#undef erfun
#undef np_reset_nr_rng_state
#undef ran3
#undef gasdev
#undef chidev
#undef ncom
#undef pcom
#undef xicom
#undef nrfunc
#undef iff

#include "conditional_search.h"
typedef struct {
  NPCSearchArena arena, *previous;
  int restrict_search, integer, n, itmax, *iterations;
  double *restricted, *point, **directions, ftol, tol, small, *value;
  double (*objective)(double *);
  int saved_ncom;
  double *saved_pcom,*saved_xicom,(*saved_nrfunc)(double *);
} NPCSearchCall;
static SEXP npc_search_execute(void *raw) {
  NPCSearchCall *c=raw;
  npc_private_powell(c->restrict_search,c->integer,c->restricted,c->point,
    c->directions,c->n,c->ftol,c->tol,c->small,c->itmax,c->iterations,
    c->value,c->objective);
  return R_NilValue;
}
static void npc_search_cleanup(void *raw) {
  NPCSearchCall *c=raw;
  while(c->arena.first)npc_search_release(c->arena.first+1);
  npc_search_arena=c->previous;
  npc_private_ncom=c->saved_ncom;npc_private_pcom=c->saved_pcom;
  npc_private_xicom=c->saved_xicom;npc_private_nrfunc=c->saved_nrfunc;
}
void np_conditional_powell(int restrict_search,int integer,double *restricted,
  double *point,double **directions,int n,double ftol,double tol,double small,
  int itmax,int *iterations,double *value,double (*objective)(double *)) {
  NPCSearchCall c={.previous=npc_search_arena,.restrict_search=restrict_search,
    .integer=integer,.n=n,.itmax=itmax,.iterations=iterations,
    .restricted=restricted,.point=point,.directions=directions,
    .ftol=ftol,.tol=tol,.small=small,.value=value,.objective=objective,
    .saved_ncom=npc_private_ncom,.saved_pcom=npc_private_pcom,
    .saved_xicom=npc_private_xicom,.saved_nrfunc=npc_private_nrfunc};
  npc_search_arena=&c.arena;
  R_ExecWithCleanup(npc_search_execute,&c,npc_search_cleanup,&c);
}
