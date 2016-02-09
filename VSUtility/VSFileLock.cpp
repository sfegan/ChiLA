//-*-mode:c++; mode:font-lock;-*-

/*! \file VSFileLock.cpp

  File locking utility.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 1.4 $
  \date       09/19/2009

  $Id: VSFileLock.cpp,v 1.4 2009/12/16 00:19:20 matthew Exp $

*/

#include <fcntl.h>
#include <iostream>
#include <cerrno>
#include <cstring>
#include <sys/file.h>

#include <VSFileLock.hpp>

using namespace VERITAS;

// ============================================================================
// VSFileLockBSD
// ============================================================================
int VSFileLockBSD::s_lfp = 0;

VSFileLockBSD::VSFileLockBSD(const std::string& file):
  m_file(file), m_lfp()
{
  acquireLock(m_file,m_lfp);
}

VSFileLockBSD::~VSFileLockBSD()
{
  releaseLock(m_file,m_lfp);
}

bool VSFileLockBSD::testLock(const std::string& file) 
{
  int lfp = open(file.c_str(),O_RDONLY);
  if(lfp<0) 
    {
      lfp = 0;
      std::cerr << "Failed to open/create file for locking: " 
		<< file << std::endl;
      std::cerr << strerror(errno) << std::endl;
      return false;
    }

  if(flock(lfp,LOCK_EX | LOCK_NB)) 
    {
      std::cerr << "VSFileLockBSD::testLock(): "
		<< "Could not acquire lock on file: " 
		<< file << std::endl;  
      std::cerr << strerror(errno) << std::endl;
      return false;
    }

  close(lfp);
  return true;
}

bool VSFileLockBSD::acquireLock(const std::string& file, bool block) 
{
  return acquireLock(file,s_lfp,block);
}

bool VSFileLockBSD::acquireLock(const std::string& file, int& lfp, bool block) 
{
  //  ,S_IRUSR|S_IWUSR
  lfp = open(file.c_str(),O_RDONLY);
  if(lfp<0) 
    {
      lfp = 0;
      std::cerr << "Failed to open/create file for locking: " 
		<< file << std::endl;
      std::cerr << strerror(errno) << std::endl;
      return false;
    }

  int flag = LOCK_EX;
  if(!block) flag |= LOCK_NB;

  if(flock(lfp,flag)) 
    {
      std::cerr << "VSFileLockBSD::acquireLock(): "
		<< "Could not acquire lock on file: " 
		<< file << std::endl;
      std::cerr << strerror(errno) << std::endl;
      return false;
    }

  return true;
}

bool VSFileLockBSD::releaseLock(const std::string& file) 
{
  return releaseLock(file,s_lfp);
}

bool VSFileLockBSD::releaseLock(const std::string& file, int lfp) 
{
  if(flock(lfp,LOCK_UN)) 
    {
      std::cerr << "VSFileLockBSD::releaseLock(): "
		<< "Could not release lock on file: " 
		<< file << std::endl;
      std::cerr << strerror(errno) << std::endl;
      return false;
    }

  close(lfp);
  return true;
}

// ============================================================================
// VSFileLock
// ============================================================================
int VSFileLock::s_lfp = 0;

bool VSFileLock::testLock(const std::string& file) 
{
  s_lfp = open(file.c_str(),O_WRONLY|O_CREAT,S_IRUSR|S_IWUSR);
  if(s_lfp<0) 
    {
      s_lfp = 0;
      std::cerr << "Failed to open/create lock file: " << file << std::endl;
      std::cerr << strerror(errno) << std::endl;
      return false;
    }

  if(lockf(s_lfp,F_TEST,0)<0) 
    {
      std::cerr << "VSFileLock::testLock(): Could not acquire lock on file: " 
		<< file << std::endl;  
      std::cerr << strerror(errno) << std::endl;
      return false;
    }

  close(s_lfp);
  return true;
}

bool VSFileLock::acquireLock(const std::string& file) 
{
  s_lfp = open(file.c_str(),O_WRONLY|O_CREAT,S_IRUSR|S_IWUSR);
  if(s_lfp<0) 
    {
      s_lfp = 0;
      std::cerr << "Failed to open/create lock file: " << file << std::endl;
      std::cerr << strerror(errno) << std::endl;
      return false;
    }

  if(lockf(s_lfp,F_TLOCK,0)<0) 
    {
      std::cerr << "VSFileLock::acquireLock(): "
		<< "Could not acquire lock on file: " 
		<< file << std::endl;
      std::cerr << strerror(errno) << std::endl;
      return false;
    }

  return true;
}

bool VSFileLock::releaseLock(const std::string& file) 
{
  if(lockf(s_lfp,F_TLOCK,0)<0) 
    {
      std::cerr << "VSFileLock::releaseLock(): "
		<< "Could not release lock on file: " 
		<< file << std::endl;
      std::cerr << strerror(errno) << std::endl;
      return false;
    }

  close(s_lfp);
  return true;
}
