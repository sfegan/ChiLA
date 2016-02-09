//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAssert.hpp

  Assert using throw

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       12/04/2007

  $Id: VSAssert.cpp,v 1.2 2008/07/13 04:31:18 matthew Exp $

*/

#include<sstream>
#include<VSAssert.hpp>

using namespace VERITAS;

VSAssert::~VSAssert() throw()
{
  // nothing to see here
}

const char* VSAssert::what() const throw()
{
  return m_message.c_str();
}

void VSAssert::msg()
{
  std::ostringstream stream;
  stream << m_file << ':' << m_line << ": " << m_func << ": Assertion `"
	 << m_assertion << "' failed.";
  m_message = stream.str();
}

void VERITAS::__vsassert_fail (const char *assertion, const char *file,
		               unsigned int line, const char *function)
    throw (VSAssert)
{
  throw VSAssert(assertion, file, line, function);
}
