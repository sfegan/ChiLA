//-*-mode:c++; mode:font-lock;-*-

/*! \file VSFileLock.hpp

  File locking utility.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 1.2 $
  \date       09/19/2009

  $Id: VSFileLock.hpp,v 1.2 2009/12/16 00:19:20 matthew Exp $

*/

#ifndef VSFILELOCK_HPP
#define VSFILELOCK_HPP

#include<string>

namespace VERITAS
{
  // ==========================================================================
  // Perform BSD file locking with flock(2)
  // ==========================================================================
  class VSFileLockBSD
  {
  public:

    VSFileLockBSD(const std::string& file);
    ~VSFileLockBSD();

    static bool testLock(const std::string& file);

    static bool acquireLock(const std::string& file, int& lfp, 
			    bool block = true);
    static bool acquireLock(const std::string& file, bool block = true);

    static bool releaseLock(const std::string& file, int lfp);   
    static bool releaseLock(const std::string& file);   

  private:

    std::string m_file;
    int         m_lfp;

    static int s_lfp;
  };

  // ==========================================================================
  // Perform POSIX file locking with fcntl(2)
  // ==========================================================================
  class VSFileLock
  {
  public:
    static bool testLock(const std::string& file);
    static bool acquireLock(const std::string& file);
    static bool releaseLock(const std::string& file);   

  private:
    static int s_lfp;
  };
}

#endif // VSFILELOCK_HPP
