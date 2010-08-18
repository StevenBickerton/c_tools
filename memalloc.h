/* header for memory allocation code */

#ifndef MEM_LIB
#define MEM_LIB

/* defines */
#define MAX_FILE 128

/* typedefs */
typedef struct {
  int id;
  size_t size;
  char file[MAX_FILE];
  int const line;
  void *ptr;
} MEMBLOCK;

typedef struct {
  int nalloc;
  int id;
  MEMBLOCK **blocks;
} MEMLIST;

/*  */
extern MEMLIST memlist;


/* macros */
#define mem_malloc(S)     p_malloc(__FILE__, __LINE__, &memlist, S)
#define mem_realloc(P, S) p_realloc(__FILE__, __LINE__, &memlist, P, S)
#define mem_free(S)       p_free(&memlist, S);

/* prototypes */
void *q_malloc(size_t size);
void *q_realloc(void *ptr, size_t size);
void *p_malloc(char const *file, int const line, MEMLIST *mem, size_t size);
void *p_realloc(char const *file, int const line, MEMLIST *mem, void *ptr, size_t size);
void p_free(MEMLIST *mem, void *vptr);
void free_leaks(MEMLIST *mem);
int show_leaks(MEMLIST *mem);

#endif	/* MEM_LIB */

