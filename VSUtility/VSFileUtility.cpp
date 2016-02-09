//-*-mode:c++; mode:font-lock;-*-

/*! \file VSFileUtility.cpp

  Miscellaneous utility file functions

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       11/11/2006
*/

// New header for consistency with OAWG package

/**
 * \class VSFileUtility
 * \ingroup common
 * \brief This is a one-line description of this cpp file.
 *
 * Here is a tedious verbose multi-line description of all
 * the details of the code, more than you would
 * ever want to read. Generally, all the important documentation
 * goes in the .cpp files.
 *
 * Original Author: Stephen Fegan
 * $Author: sfegan $
 * $Date: 2007/01/29 02:52:59 $
 * $Revision: 1.5 $
 * $Tag$
 *
 **/

// These deinitions tell the makefile which library the cpp file
// should be included in
// VA_LIBRARY_TAG: libSP24common.a
// VA_LIBRARY_TAG: libSP24commonLite.a

#include <unistd.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <pwd.h>

#include <VSDataConverter.hpp>

#include "VSFileUtility.hpp"

#ifndef _OAWG
using namespace VERITAS;
#endif

bool VSFileUtility::exists(const std::string& filename)
{
  if(filename.empty())return false;
  struct stat statbuf;
  if(stat(filename.c_str(),&statbuf)<0)return false;
  return true;
}

bool VSFileUtility::isFile(const std::string& filename)
{
  if(filename.empty())return false;
  struct stat statbuf;
  if(stat(filename.c_str(),&statbuf)<0)return false;
  return S_ISREG(statbuf.st_mode);
}

bool VSFileUtility::isDirectory(const std::string& filename)
{
  if(filename.empty())return false;
  struct stat statbuf;
  if(stat(filename.c_str(),&statbuf)<0)return false;
  return S_ISDIR(statbuf.st_mode);
}

bool VSFileUtility::isReadable(const std::string& filename)
{
  if(filename.empty())return false;
  return access(filename.c_str(),R_OK)==0;
}

bool VSFileUtility::isWritable(const std::string& filename)
{
  if(filename.empty())return false;
  return access(filename.c_str(),W_OK)==0;
}

bool VSFileUtility::canWriteFile(const std::string& filename)
{
  if(filename.empty())return false;

  // Test whether the file exists and is writable OR does not exist
  // but the directory in which it would be created is writable
  if(exists(filename))
    return isFile(filename)&&isWritable(filename);
  else
    {
      std::string directory = filename;
      unsigned dpos = 0;
      for(unsigned ipos=0;ipos<filename.length();ipos++)
	if(filename[ipos]=='/')dpos=ipos;
      if((dpos==0)&&(filename[dpos]!='/'))directory=".";
      else directory=directory.substr(0,dpos);
      return isDirectory(directory)&&isWritable(directory);
    }
}

void VSFileUtility::expandFilename(std::string& filename)
{ 
  /* Do leading tilde expansion, replace with home directory */
  if((!filename.empty())&&(filename[0]=='~'))
    {
      unsigned ipos = 1;
      while((ipos!=filename.length())&&(filename[ipos]!='/'))ipos++;
      std::string user = filename.substr(1,ipos-1);
      if(user.empty())
	{
	  char* login = getlogin();
	  if(login)user = login;
	  else 
	    {
	      struct passwd* pwd = getpwuid(getuid());
	      if(pwd)user = pwd->pw_name;
	    }
	}

      if(!user.empty())
	{
	  struct passwd* pwd = getpwnam(user.c_str());
	  if(pwd)filename.replace(0,ipos,pwd->pw_dir);
	}
    }
}

unsigned VSFileUtility::extractNumberFromFilename(const std::string& filename)
{
  std::string::const_iterator longest_num_start = filename.begin();

  for(std::string::const_iterator ichar = filename.begin();
      ichar!=filename.end();ichar++)if(*ichar=='/')longest_num_start=ichar;
  if(*longest_num_start=='/')longest_num_start++;
  
  std::string::const_iterator longest_num_end = longest_num_start;
  std::string::const_iterator number_begin = longest_num_start;

  for(std::string::const_iterator ichar = longest_num_start; 
      ichar!=filename.end(); ichar++)
    {
      if(isdigit(*ichar))
	{
	  if(number_begin==filename.end())number_begin=ichar;
	}
      else if(number_begin!=filename.end())
	{
	  if((ichar-number_begin)>(longest_num_end-longest_num_start))
	    longest_num_start = number_begin, longest_num_end = ichar;
	  number_begin = filename.end();	  
	}
    }

  if(number_begin!=filename.end()&&
     ((filename.end()-number_begin)>(longest_num_end-longest_num_start)))
    longest_num_start = number_begin, longest_num_end = filename.end();
  
  unsigned n = 0;

  if(longest_num_end != longest_num_start)
    VSDataConverter::
      fromString(n,std::string(longest_num_start,longest_num_end));

  return n;
}

void  VSFileUtility::
replaceQuestionWithNumber(std::string& filename, unsigned n)
{
  std::string::size_type iquestion = filename.find('?');
  if(iquestion != filename.npos)
    filename.replace(iquestion,1, VSDataConverter::toString(n));
}
