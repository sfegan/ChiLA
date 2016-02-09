//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventSupplierTextFile.hpp

  Input of event from a text file

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n        

  \version    0.1
  \date       08/14/2006
  \note
*/

#ifndef VSEVENTSUPPLIERTEXTFILE_HPP
#define VSEVENTSUPPLIERTEXTFILE_HPP

#include<vector>
#include<fstream>

#include<VSSimDB.hpp>
#include<VSEventSupplier.hpp>

//! VERITAS namespace
namespace VERITAS 
{

  class VSEventSupplierTextFile: public VSEventSupplier
  {
  public:    
    VSEventSupplierTextFile(const std::string& event_filename, const std::string& pe_filename);
    virtual ~VSEventSupplierTextFile();
    virtual std::vector<SimParam> getSimParam();
    virtual bool getNextEvent(Event& e);
    virtual uint32_t getNumEvents();
    void readEventFile();
  private:

    std::ifstream                  fEventFileHandle;
    std::ifstream                  fPEFileHandle;

    std::string                    fEventFilename;
    std::string                    fPEFilename;
    unsigned                       fCompleteEvents;
    unsigned                       fIEvent;
    std::vector<VSSimDBEventData>  fEvents;

  };

}

#endif // defined VSEVENTSUPPLIERTEXTFILE_HPP
