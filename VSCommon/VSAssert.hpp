//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAssert.hpp

  Assert using throw

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       12/04/2007

  $Id: VSAssert.hpp,v 1.3 2007/12/05 01:34:55 sfegan Exp $

*/

#ifndef VSASSERT_HPP
#define VSASSERT_HPP

#include <stdexcept>
#include <cassert>
#include <string>

namespace VERITAS
{

  class VSAssert: public std::logic_error
  {
  public:
    VSAssert(const char* assertion, const char* file, unsigned line,
	     const char* func):
      logic_error("Assertion failed"), m_assertion(assertion),
      m_file(file), m_line(line), m_func(func) { msg(); }

    virtual ~VSAssert() throw();

    virtual const char* what() const throw();

    const std::string& assertion() { return m_assertion; }
    const std::string& file() { return m_file; }
    unsigned line() { return m_line; }
    const std::string& func() { return m_func; }
					    
  private:
    std::string m_message;
    std::string m_assertion;
    std::string m_file;
    unsigned    m_line;
    std::string m_func;

    void msg();
  };
  
  void __vsassert_fail (const char *assertion, const char *file,
			unsigned int line, const char *function)
    throw (VSAssert) __attribute__ ((__noreturn__));
}

#if (defined NDEBUG)||(defined NVSASSERT)

#define vsassert(expr)				\
  assert(expr)

#else

#define __VSASSERT_VOID_CAST static_cast<void>

#define vsassert(expr)							\
  ((expr)								\
   ? __VSASSERT_VOID_CAST (0)						\
   : VERITAS::__vsassert_fail (#expr, __FILE__, __LINE__, __PRETTY_FUNCTION__))

#endif

#endif // ifndef VSASSERT_HPP
