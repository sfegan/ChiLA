//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventSupplierTextFile.cpp

  Input of event from a text file 

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n        

  \version    0.1
  \date       08/14/2006
  \note
*/

#include <cassert>

#include "VSEventSupplierTextFile.hpp"

using namespace VERITAS;

VSEventSupplierTextFile::
VSEventSupplierTextFile(const std::string& event_filename, const std::string& pe_filename)
  : VSEventSupplier(), fEventFileHandle(), fPEFileHandle(),
    fEventFilename(event_filename),fPEFilename(pe_filename),
    fCompleteEvents(), fIEvent(), fEvents()
{

  fEventFileHandle.open( fEventFilename.c_str(), std::ios::in );
  fPEFileHandle.open( fPEFilename.c_str(), std::ios::in );

  assert( fEventFileHandle.is_open() );
  assert( fPEFileHandle.is_open() );


  readEventFile();
  fCompleteEvents = fEvents.size();
  fIEvent=0;
}

VSEventSupplierTextFile::~VSEventSupplierTextFile()
{

}

std::vector<VSEventSupplier::SimParam> VSEventSupplierTextFile::getSimParam()
{
  std::vector<SimParam> tp_vec;
  return tp_vec;
}

void VSEventSupplierTextFile::readEventFile() {

  VSSimDBEventData event;

  unsigned int numHitTelescopes=0;

  while(fEventFileHandle 
	>> event.fEventID 
	>> event.fTargetZenithRad 
	>> event.fTargetAzimuthRad
	>> event.fPrimaryZenithRad   
	>> event.fPrimaryAzimuthRad
	>> event.fPrimaryCoreEastM
	>> event.fPrimaryCoreNorthM
	>> event.fPrimaryCoreUpASLM
	>> numHitTelescopes 
	>> event.fEventComplete && !fEventFileHandle.eof()) {

    fEvents.push_back(event);
  }


}

bool VSEventSupplierTextFile::getNextEvent(Event& e)
{
  if(fIEvent >= fEvents.size())return false;

  e.fEventID                 = fEvents[fIEvent].fEventID;
  e.fTargetZenithRad         = fEvents[fIEvent].fTargetZenithRad;
  e.fTargetAzimuthRad        = fEvents[fIEvent].fTargetAzimuthRad;
  e.fPrimarySpeciesCORSIKA   = 0;
  e.fPrimaryEnergyTeV        = 0;
  e.fPrimaryZenithRad        = fEvents[fIEvent].fPrimaryZenithRad;
  e.fPrimaryAzimuthRad       = fEvents[fIEvent].fPrimaryAzimuthRad;
  e.fPrimaryCoreEastM        = fEvents[fIEvent].fPrimaryCoreEastM;
  e.fPrimaryCoreNorthM       = fEvents[fIEvent].fPrimaryCoreNorthM;
  e.fPrimaryCoreUpASLM       = fEvents[fIEvent].fPrimaryCoreUpASLM;
  e.fSamplingRadiusM         = 0;
  e.fTableIndex              = 0;
  e.fEventComplete           = fEvents[fIEvent].fEventComplete;

  fIEvent++;

  std::map< unsigned int, std::map< unsigned int, VSEventSupplier::Pixel > > pixels_map;

  unsigned int eventID = e.fEventID;
  unsigned int eventNumber, telescopeID, pixelID;
  unsigned int lastPos = 0;
  float time;
  while ( ( fPEFileHandle >> eventNumber >> telescopeID >> pixelID >> time 
	    && eventNumber <= eventID && !fPEFileHandle.eof() ) ) {

    VSEventSupplier::PE pe;
    pe.fTimeNS = time;

    pixels_map[telescopeID][pixelID].fPEs.push_back(pe);
    pixels_map[telescopeID][pixelID].fPixelID = pixelID;


    lastPos = fPEFileHandle.tellg();

  }

  fPEFileHandle.seekg(lastPos); // return the file stream to its last position

  for(std::map< unsigned int, std::map< unsigned int, VSEventSupplier::Pixel > >::iterator scopeItr = pixels_map.begin();
      scopeItr != pixels_map.end(); ++scopeItr) {
    VSEventSupplier::Scope scope;
    scope.fScopeID = scopeItr->first;

    for(std::map< unsigned int, VSEventSupplier::Pixel >::iterator pixelItr = scopeItr->second.begin();
	pixelItr != scopeItr->second.end(); ++pixelItr) {    
      scope.fPixels.push_back(pixelItr->second);
    }

    e.fScopes.push_back(scope);

  }

  return true;
}

uint32_t VSEventSupplierTextFile::getNumEvents()
{
  return fCompleteEvents;
}
