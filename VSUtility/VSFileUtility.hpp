//-*-mode:c++; mode:font-lock;-*-

// This file is shared between the Simulations and OAWG. To avoid a
// fork in development, please coordinate changes with Stephen Fegan.
// email: sfegan@astro.ucla.edu

// Old header for consistency with Simulations package

/*! \file VSFileUtility.hpp

  Miscellaneous utility file functions

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       11/11/2006
*/

// New header for consistency with OAWG package

/**
 *
 * Original Author: Stephen Fegan
 * $Author: sfegan $
 * $Date: 2007/01/29 02:52:59 $
 * $Revision: 1.6 $
 * $Tag$
 *
 **/

#ifndef VSFILEUTILITY_HPP
#define VSFILEUTILITY_HPP

#include<string>

#ifndef _OAWG
namespace VERITAS
{
#endif

  class VSFileUtility
  {
  public:
    //! Return true if filename exists on the filesystem
    static bool exists(const std::string& filename);

    //! Return true if filename exists and corresponds to a regular file
    static bool isFile(const std::string& filename);
    //! Return true if filename exists and corresponds to a directory
    static bool isDirectory(const std::string& filename);

    //! Return true if user has permission to read from filename
    static bool isReadable(const std::string& filename);
    //! Return true if user has permission to write to filename
    static bool isWritable(const std::string& filename);

    //! Return true if filename exists, is regular and writable or
    // if file does not exist but its parent directory is writable
    static bool canWriteFile(const std::string& filename);

    //! Expand the leading tilde in the filename, replacing it with the home
    // directory of the current or specified user
    static void expandFilename(std::string& filename);

    //! Extract the longest positive number from the filename portion
    // of a path/file specification. For example, 
    // "/veritas/data/V12345_v2.cvbf" would return 12345.
    static unsigned extractNumberFromFilename(const std::string& filename);

    //! Replace the first occurance of a question mark '?' with an
    // unsigned integer. For example passing: filename="?.root" and n=12345
    // would result in filename being returned as "12345.root"
    static void replaceQuestionWithNumber(std::string& filename, unsigned n);
  };

#ifndef _OAWG
}
#endif

#endif // VSFILEUTILITY_HPP
