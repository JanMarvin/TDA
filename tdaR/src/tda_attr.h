#ifndef TDA_ATTR_H
#define TDA_ATTR_H

#include <stdio.h>

#if defined(__MINGW_PRINTF_FORMAT)
#  define TDA_PRINTF(f, a) __attribute__((format(__MINGW_PRINTF_FORMAT, f, a)))
#elif defined(__GNUC__) || defined(__clang__)
#  define TDA_PRINTF(f, a) __attribute__((format(printf, f, a)))
#else
#  define TDA_PRINTF(f, a)
#endif

#endif
