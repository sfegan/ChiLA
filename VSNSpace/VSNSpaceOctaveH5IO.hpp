/*! \file  VSNSpaceOctaveH5IO.hpp

  NSpace IO to Octave H5 file

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       06/16/2006

*/

#include<VSOctaveIO.hpp>
#include<VSNSpace.hpp>

#ifndef VSNSPACEOCTAVEH5IO_HPP
#define VSNSPACEOCTAVEH5IO_HPP

namespace VERITAS 
{

  class VSNSpaceOctaveH5IO: public VSNSpaceIO
  {
  public:
    VSNSpaceOctaveH5IO() { /* nothing to see here */ }
    virtual ~VSNSpaceOctaveH5IO();

    virtual bool readSpace(const std::string& filename,
			   VSNSpace::Space& space);
    virtual bool writeSpace(const std::string& filename,
			    const VSNSpace::Space& space);
    virtual bool readVolume(const std::string& filename, 
			    VSNSpace::Volume& vol);
    virtual bool writeVolume(const std::string& filename,
			     const VSNSpace::Volume& vol, bool sparse = false);
    virtual bool readOrdering(const std::string& filename,
			      VSNSpace::Ordering& ordering);
    virtual bool writeOrdering(const std::string& filename,
			       const VSNSpace::Ordering& ordering);
    virtual bool readHistogram(const std::string& filename,
			       VSNSpace& hist);
    virtual bool writeHistogram(const std::string& filename,
				const VSNSpace& hist);

    bool readSpace(VSOctaveH5ReaderStruct* s, VSNSpace::Space& space);
    static bool isSpace(VSOctaveH5ReaderStruct* s);
    bool writeSpace(VSOctaveH5WriterStruct* s, const VSNSpace::Space& space);

    bool readVolume(VSOctaveH5ReaderStruct* s, VSNSpace::Volume& vol);
    static bool isVolume(VSOctaveH5ReaderStruct* s);
    bool writeVolume(VSOctaveH5WriterStruct* s, const VSNSpace::Volume& vol,
		     bool sparse = false);

    bool readHistogram(VSOctaveH5ReaderStruct* s, VSNSpace& hist);
    static bool isHistogram(VSOctaveH5ReaderStruct* s);
    bool writeHistogram(VSOctaveH5WriterStruct* s, const VSNSpace& hist);
  };

}

#endif // defined VSNSPACEOCTAVEH5IO_HPP
