//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAFunctor.hpp
  Generic Functor class.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n        

  \version    0.1
  \date       11/10/2008
*/

#ifndef VSAFUNCTOR_HPP
#define VSAFUNCTOR_HPP

namespace VERITAS
{
  namespace VSAMath
  {
    template< typename T >
    class Functor
    {
    public:
      Functor() {}
      virtual ~Functor() {}

      virtual double operator() (const T& x) const = 0;
    };

    template< typename T1, typename T2, typename T3 >
    class Functor3
    {
    public:
      Functor3() {}
      virtual ~Functor3() {}

      virtual double operator() 
	(const T1& x1, const T2& x2, const T3& x3) const = 0;
    };
  }
}

#endif // VSAFUNCTOR_HPP
