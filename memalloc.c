/*            */
/* Steven Bickerton */
/* Dept. of Astrophysical Sciences, Princeton University */
/* bick@astro.princeton.edu*/
/* Created: Tue Jul  1, 2008  10:00:43 DST */
/* Host: bender.astro.princeton.edu */
/* Working Directory: /Users/bick/working/photoz/eazy-0.99/src  */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <errno.h>
#include <string.h>
#include <stdarg.h>

#include "memalloc.h"

MEMLIST memlist = {0, -1, NULL};


void *q_malloc(size_t size) {
  void *tmp = malloc( size );
  if ( tmp == NULL ) {
    perror("q_malloc() failed\n"); 
    exit(EXIT_FAILURE);
  }
  return tmp;
}
void *q_realloc(void *ptr, size_t size) {
  void *tmp = realloc(ptr, size);
  if ( tmp == NULL ) {
    perror("q_realloc() failed\n"); 
    exit(EXIT_FAILURE);
  }
  return tmp;
}

/* A malloc routine which maintains a list of allocated memory */
void *p_malloc(char const *file, int const line, MEMLIST *mem, size_t size) {
  
  MEMBLOCK *tmp = malloc( sizeof(MEMBLOCK) + size );
  assert (tmp != NULL);
  
  tmp->id = ++mem->id;
  tmp->size = size;
  memcpy(tmp->file, file, MAX_FILE);
  *(int *)&tmp->line = line;
  tmp->ptr = tmp + 1;
  
  if (mem->id >= mem->nalloc) {
    mem->nalloc = 2*mem->nalloc + 1;
    mem->blocks = realloc(mem->blocks, mem->nalloc*sizeof(MEMBLOCK *));
    assert(mem->blocks != NULL);
  }
  mem->blocks[tmp->id] = tmp;
  
  return tmp->ptr;
}

/* A realloc routine which maintains a list of allocated memory */
void *p_realloc(char const *file, int const line, MEMLIST *mem, 
		void *vptr, size_t size) {
  
  MEMBLOCK *ptr = (MEMBLOCK *) vptr - 1;
  assert (ptr != NULL);		/* must be allocated and not yet dealloc'd */
  int id_tmp = ptr->id;
  assert ( id_tmp >= 0  &&  id_tmp <= mem->id ); /* id in valid range */

  MEMBLOCK *tmp = realloc( ptr, sizeof(MEMBLOCK) + size );
  assert (tmp != NULL);
  
  tmp->id = id_tmp; 		/* don't change the id */
  tmp->size = size;		/* reset the size */
  memcpy(tmp->file, file, MAX_FILE); /* ... perhaps new position at realloc */
  *(int *)&tmp->line = line;	/*  ... new line too */
  tmp->ptr = tmp + 1;

  mem->blocks[tmp->id] = tmp;	/*  */

  return tmp->ptr;
}


/* The free() equivalent for p_ez_malloc */
void p_free(MEMLIST *mem, void *vptr) {
  if (!vptr) return;

  MEMBLOCK *ptr = (MEMBLOCK *)vptr - 1;
  int id_tmp = ptr->id;
  
  assert(id_tmp >= 0 && id_tmp < mem->nalloc);
  mem->blocks[id_tmp] = NULL;
  free(ptr);
}



/* a routine to free all memory allocated with p_ez_malloc */
void free_leaks(MEMLIST *mem) {

  for (int i = 0; i<=mem->id; i++) {
    if (mem->blocks[i] != NULL) {
      //printf("Freeing %-6d %8ld %s:%d\n", mem->blocks[i]->id,
      //(int) mem->blocks[i]->size, mem->blocks[i]->file, 
      //mem->blocks[i]->line);
      mem_free(mem->blocks[i]->ptr);
    }
  }
  /*  */
  mem->blocks = NULL;
  free(mem->blocks);
  mem->nalloc = 0;
  mem->id = -1;

}
	

/* Print info about all memory allocated with p_ez_malloc(),
   but not freed with ez_free() */
int show_leaks(MEMLIST *mem) {
  int nleak = 0;
  size_t sum_leak_size = 0;
  for (int i = 0; i <= mem->id; i++) {
    if (mem->blocks[i] != NULL) {
      ++nleak;
      sum_leak_size += mem->blocks[i]->size;
      printf("Block %-6d %8ld %s:%d\n", mem->blocks[i]->id,
	     (long) mem->blocks[i]->size, mem->blocks[i]->file, 
	     mem->blocks[i]->line);
    }
  }
  if (nleak) {
    printf("Total %.3fM leaked\n", (double) sum_leak_size/(1024*1024));
  }
  return nleak;
}


