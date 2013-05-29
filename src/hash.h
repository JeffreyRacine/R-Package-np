
#define TH_ERROR 0
#define TH_SUCCESS 1
#define TH_FAILURE 2

#define TH_ENTER 0
#define TH_SEARCH 1


struct th_entry {
  union {
    double dkey;
    uint64_t ukey;
  };
  
  int data;
  uint8_t vacant;
};


struct th_table {
  struct th_entry * ht;

  size_t size, maxsize;
};

int thcreate_r(size_t nel, struct th_table * tab);
int thdestroy_r(struct th_table *tab);

int thsearch_r(const struct th_entry * const q, int action, struct th_entry ** ret, struct th_table * tab);
