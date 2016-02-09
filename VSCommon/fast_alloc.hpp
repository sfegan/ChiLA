//-*-mode:c++; mode:font-lock;-*-

/*! \file fast_alloc.hpp

  Fast allocation and free

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/06/2007

  $Id: fast_alloc.hpp,v 1.1 2007/12/04 17:15:04 sfegan Exp $

*/

#include <cassert>

#ifdef USEALLOCA
#include <alloca.h>
#define FASTCALLOC(T,n) static_cast<T*>(alloca(n*sizeof(T)))
#define FASTMALLOC(n) alloca(n)
#define FASTFREE(x) {}
#else
#include <cstdlib>
#define FASTCALLOC(T,n) static_cast<T*>(calloc(n,sizeof(T)))
#define FASTMALLOC(n) malloc(n)
#define FASTFREE(x) free(static_cast<void*>(x))
#endif
