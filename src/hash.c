#include <stdint.h>
#include <stdlib.h>
#include "hash.h"


size_t thfhash0(uint64_t key, size_t n){
  return (key % n);
}


int thcreate_r(size_t nel, struct th_table * tab){
  int i;
  tab->ht = NULL;
  tab->ht = (struct th_entry *)malloc(nel*sizeof(struct th_entry));
  if(tab->ht == NULL) return TH_ERROR;
  tab->size = 0;
  tab->maxsize = nel;

  for(i = 0; i < nel; i++)
    tab->ht[i].vacant = 1;

  return TH_SUCCESS;
}

int thdestroy_r(struct th_table *tab){
  free(tab->ht);
  tab->size = 0;
  tab->maxsize = 0;
  return TH_SUCCESS;
}


int thsearch_r(const struct th_entry * const q, int action, struct th_entry ** ret, struct th_table * tab){
  size_t it, i;

  struct th_entry * ht = NULL;
  const size_t nel = tab->maxsize;

  it = thfhash0(q->key.ukey,tab->maxsize);

  ht = tab->ht;

  if(action == TH_ENTER) {
    if(tab->size == tab->maxsize) return TH_FAILURE;
    if(tab->size == 0){
      ht[it].key.ukey = q->key.ukey;
      ht[it].data = q->data;
      ht[it].vacant = 0;
      tab->size += 1;
      *ret = ht+it;
      return TH_SUCCESS;
    }
  }

  if((action == TH_SEARCH) && tab->size == 0) {
    *ret = NULL;
    return TH_FAILURE;
  }

  for(i = 0; i < tab->maxsize; i++, it++){
    const size_t ot = it % nel;

    if(ht[ot].vacant){
      if(action == TH_ENTER){
        ht[ot].key.ukey = q->key.ukey;
        ht[ot].data = q->data;
        ht[ot].vacant = 0;
        tab->size += 1;
        *ret = ht+ot;
        return TH_SUCCESS;
      } else {
        *ret = NULL;
        return TH_FAILURE;
      }
    }

    if(ht[ot].key.ukey == q->key.ukey){
      if(action == TH_SEARCH){
        *ret = ht + ot;
        return TH_SUCCESS;
      } else if(action == TH_ENTER){
        *ret = NULL;
        return TH_FAILURE;
      }
    }
  }

  // we should only make it here if we somehow search through 
  // the entire table without hitting a vacancy or our key
  *ret = NULL;
  return TH_FAILURE;
}
