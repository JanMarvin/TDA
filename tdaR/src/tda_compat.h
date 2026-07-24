#ifndef TDA_COMPAT_H
#define TDA_COMPAT_H

#include <stddef.h>

/* Stable sort with a context pointer, defined in t_sort.c.
   compar signature: int (*)(const void *a, const void *b, void *thunk)

   A merge sort rather than the C library's qsort_r/qsort_s: TDA sorts index
   arrays by data value, and where values tie the order in which the indices
   come out decides which observation is processed first.  glibc's qsort is a
   stable merge sort, so every reference output was produced with ties kept in
   input order; the Windows and BSD library sorts are not stable and reorder
   ties, which shows up as different numbers (e.g. the intermediate F of tied
   events in mple).  Doing the sort here makes every platform agree. */
void tda_qsort_r(void *base, size_t n, size_t sz,
                 int (*cmp)(const void *, const void *, void *), void *thunk);

#endif /* TDA_COMPAT_H */
