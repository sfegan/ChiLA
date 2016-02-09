//-*-mode:c++; mode:font-lock;-*-

/*! \file VSLog.hpp

  Logging to multiple streams

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/03/2006
*/

#ifndef VSLOG_HPP
#define VSLOG_HPP

#include<list>
#include<memory>
#include<iostream>
#include<iomanip>
#include<bits/char_traits.h>

namespace VERITAS
{

  class VSLog
  {
  public:
    void attachStream(std::ostream& stream) { m_streams.push_back(&stream); }
    void attachStream(std::ostream* stream) { m_streams.push_back(stream); }

    typedef std::ostream::char_type Ch;
    typedef std::ostream::traits_type Tr;

    template<typename T> VSLog& operator << (const T& x)
    {
      for(std::list<std::ostream*>::iterator i_ostream = m_streams.begin();
	  i_ostream!=m_streams.end(); i_ostream++)(**i_ostream) << x;
      return *this;
    }

    VSLog& operator << (std::basic_ostream<Ch,Tr>& (*f) (std::basic_ostream<Ch,Tr>&))
    {
      for(std::list<std::ostream*>::iterator i_ostream = m_streams.begin();
	  i_ostream!=m_streams.end(); i_ostream++)f(**i_ostream);
      return *this;
    }

    VSLog& operator << (std::ios_base& (*f) (std::ios_base&))
    {
      for(std::list<std::ostream*>::iterator i_ostream = m_streams.begin();
	  i_ostream!=m_streams.end(); i_ostream++)f(**i_ostream);
      return *this;
    }

    VSLog& operator << (std::basic_ios<Ch,Tr>& (*f) (std::basic_ios<Ch,Tr>&))
    {
      for(std::list<std::ostream*>::iterator i_ostream = m_streams.begin();
	  i_ostream!=m_streams.end(); i_ostream++)f(**i_ostream);
      return *this;
    }

    static VSLog& instance()
    {
      if(s_instance.get() == 0)s_instance.reset(new VSLog);
      return *s_instance.get();
    }

    static VSLog& i() { return instance(); }

  protected:
    VSLog(): m_streams() { /* nothing to see here */ }

  private:
    static std::auto_ptr<VSLog> s_instance;
    std::list<std::ostream*> m_streams;
  };

};

#endif // VSLOG_HPP
