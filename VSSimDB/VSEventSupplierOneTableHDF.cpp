//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventSupplierOneTableHDF.cpp

  Input of event from single table from simulations database 

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       07/08/2006
  \note
*/

#include <cassert>

#include <VSFileUtility.hpp>
#include <VSHDFEventStore.hpp>

#include "VSEventSupplierOneTableHDF.hpp"

using namespace VERITAS;

VSEventSupplierOneTableHDF::
VSEventSupplierOneTableHDF(const std::string& table_name, VSSimDB* sim_db, 
			   std::string directory,
			   bool check_hdf_event_count)
  : VSEventSupplierOneTableCommon(table_name,sim_db),
    fDirectory(directory), fWorkunitIds(),
    fReader(), fWorkunitFile(), fIDOfOpenWorkunit(), fNEventInWorkunit(), 
    fICompleteEventAtWorkunitStart(),
    fEV(), fSC(), fPX(), fPE(),
    fREventID(), fRTargetZenithRad(), fRTargetAzimuthRad(), 
    fRPrimaryZenithRad(), fRPrimaryAzimuthRad(),
    fRPrimaryCoreEastM(), fRPrimaryCoreNorthM(), fRPrimaryCoreUpASLM(),
    fRNumHitScopes(), fRScopeStart(), fRScopeCount(),
    fRScopeID(), fRScopeZenithRad(), fRScopeAzimuthRad(), 
    fRPixelStart(), fRPixelCount(),
    fRPixelID(), fRPEStart(), fRPECount(),
    fRPixelTimeNS()
{
  initialize(check_hdf_event_count);
}

VSEventSupplierOneTableHDF::
VSEventSupplierOneTableHDF(const std::string& table_name, VSSimDB* sim_db,
			   uint32_t nevent_to_supply, 
			   double ievent_to_start_fraction, 
			   std::string directory,
			   bool check_hdf_event_count)
  : VSEventSupplierOneTableCommon(table_name,sim_db,nevent_to_supply,
				  ievent_to_start_fraction),
    fDirectory(directory), fWorkunitIds(),
    fReader(), fWorkunitFile(), fIDOfOpenWorkunit(), fNEventInWorkunit(), 
    fICompleteEventAtWorkunitStart(),
    fEV(), fSC(), fPX(), fPE(),
    fREventID(), fRTargetZenithRad(), fRTargetAzimuthRad(), 
    fRPrimaryZenithRad(), fRPrimaryAzimuthRad(),
    fRPrimaryCoreEastM(), fRPrimaryCoreNorthM(), fRPrimaryCoreUpASLM(),
    fRNumHitScopes(), fRScopeStart(), fRScopeCount(),
    fRScopeID(), fRScopeZenithRad(), fRScopeAzimuthRad(), 
    fRPixelStart(), fRPixelCount(),
    fRPixelID(), fRPEStart(), fRPECount(),
    fRPixelTimeNS()
{
  initialize(check_hdf_event_count);
}

VSEventSupplierOneTableHDF::~VSEventSupplierOneTableHDF()
{
  close();
}

void VSEventSupplierOneTableHDF::initialize(bool check_hdf_event_count)
{
  unsigned nevents = 0;

  std::vector<unsigned> workunits;
  fSimDB->getWorkunitIDs(fTableName, workunits);
  for(unsigned iworkunit=0;iworkunit<workunits.size();iworkunit++)
    {
      std::string filename = 
	fDirectory
	+ std::string("/")
	+ VSHDFEventStore::nameForWorkunitID(workunits[iworkunit]);
  
      if(!VSFileUtility::isFile(filename))
	{
	  std::cerr << "Warning: " << filename
		    << " does not exist.  Table: " << fTableName
		    << std::endl;
	}
      else if(!VSOctaveH5ReaderStruct::isHDF5(filename))
	{
	  std::cerr << "Warning: " << filename
		    << " is not a valid HDF5 file.  Table: " << fTableName
		    << std::endl;
	}
      else if(check_hdf_event_count)
	{
	  try
	    {
	      const unsigned Ne = 100;
	      fWorkunitIds.push_back(workunits[iworkunit]);
	      fReader        = new VSOctaveH5ReaderStruct(filename);
	      fEV            = fReader->readStruct("event");
	      fREventID      = fEV->readVectorElements<uint32_t>("EventID",Ne);
	      nevents += fREventID->rows();
	      delete fREventID;
	      delete fEV;
	      delete fReader;
	      fREventID = NULL;
	      fEV       = NULL;
	      fReader   = NULL;
	    }
	  catch(const VSOctaveH5Exception& e)
	    {
	      std::cerr << "Exception while opening file: " << filename
			<< " Table: " << fTableName
			<< std::endl;	      
	      delete fREventID;
	      delete fEV;
	      delete fReader; 
	      throw(e);
	    }
	}
      else 
	fWorkunitIds.push_back(workunits[iworkunit]);
    }

  if(check_hdf_event_count && 
     getICompleteEvent() + getNEventToSupply() > nevents)
    setStartEventAndCount(fTableName,nevents,0);
}

bool VSEventSupplierOneTableHDF::getNextEvent(Event& e)
{
  unsigned icomplete_event = getICompleteEvent();
  if(!nextEvent())return false;

  unsigned ievent_at = icomplete_event-fICompleteEventAtWorkunitStart;

  while((!fReader)||(ievent_at >= fNEventInWorkunit))
    {
      if(fReader)
	{
	  fICompleteEventAtWorkunitStart += fNEventInWorkunit;
	  close();
	  fReader = NULL;
	}

      assert(!fWorkunitIds.empty());
      assert(open(fWorkunitIds.front()));
      fWorkunitIds.pop_front();

      ievent_at = icomplete_event-fICompleteEventAtWorkunitStart;
    }

  try
    {
      e.fEventID                 = fREventID->at(ievent_at);
      e.fTargetZenithRad         = fRTargetZenithRad->at(ievent_at);
      e.fTargetAzimuthRad        = fRTargetAzimuthRad->at(ievent_at);
      e.fPrimarySpeciesCORSIKA   = fTableParam.fPrimaryID;
      e.fPrimaryEnergyTeV        = fTableParam.fEnergyTeV;
      e.fPrimaryZenithRad        = fRPrimaryZenithRad->at(ievent_at);
      e.fPrimaryAzimuthRad       = fRPrimaryAzimuthRad->at(ievent_at);
      e.fPrimaryCoreEastM        = fRPrimaryCoreEastM->at(ievent_at);
      e.fPrimaryCoreNorthM       = fRPrimaryCoreNorthM->at(ievent_at);
      e.fPrimaryCoreUpASLM       = fRPrimaryCoreUpASLM->at(ievent_at);
      e.fSamplingRadiusM         = fTableParam.fSamplingRadiusM;
      e.fTableIndex              = 0;
      e.fEventComplete           = true;

      unsigned iscope_start      = fRScopeStart->at(ievent_at);
      unsigned nscope            = fRScopeCount->at(ievent_at);

      e.fScopes.resize(nscope);

      for(unsigned iscope=0;iscope<nscope;iscope++)
	{
	  Scope& s(e.fScopes[iscope]);
	  
	  unsigned iscope_at     = iscope_start + iscope;
	  s.fScopeID             = fRScopeID->at(iscope_at);
	  s.fScopeZenithRad      = fRScopeZenithRad->at(iscope_at);
	  s.fScopeAzimuthRad     = fRScopeAzimuthRad->at(iscope_at);
	  
	  unsigned ipix_start    = fRPixelStart->at(iscope_at);
	  unsigned npix          = fRPixelCount->at(iscope_at);
	  
	  s.fPixels.resize(npix);
	  
	  for(unsigned ipix=0;ipix<npix;ipix++)
	    {
	      Pixel& p(s.fPixels[ipix]);
	      
	      unsigned ipix_at   = ipix_start + ipix;
	      p.fPixelID         = fRPixelID->at(ipix_at);
	      
	      unsigned ipe_start = fRPEStart->at(ipix_at);
	      unsigned npe       = fRPECount->at(ipix_at);
	      
	      p.fPEs.resize(npe);
	      
	      for(unsigned ipe=0;ipe<npe;ipe++)
		p.fPEs[ipe].fTimeNS = fRPixelTimeNS->at(ipe_start+ipe);
	    }
	}
    }
  catch(const VSOctaveH5Exception& e)
    {
      std::cerr << "Exception while reading file: " << fWorkunitFile
		<< " Workunit: " << fIDOfOpenWorkunit
		<< std::endl;	      
      throw(e);
    }

  return true;
}

bool VSEventSupplierOneTableHDF::open(unsigned workunit_run_id)
{
  const unsigned Ne = 200;
  const unsigned Ns = Ne * 3;
  const unsigned Nx = Ns * 10;
  const unsigned Np = Nx * 10;

  std::string filename = 
    fDirectory
    + std::string("/")
    + VSHDFEventStore::nameForWorkunitID(workunit_run_id);
  
  try
    {
      fReader             = new VSOctaveH5ReaderStruct(filename);

      if(!fReader)return false;

      fIDOfOpenWorkunit   = workunit_run_id;
      fWorkunitFile       = filename;

      // VERIFY THAT FILE IS CORRECT WORKUNIT AND IS FROM CORRECT TABLE
      VSOctaveH5ReaderStruct* tps = fReader->readStruct("table_param");

      unsigned h5_workunit_run_id;
      fReader->readScalar("workunit_run_id",h5_workunit_run_id);
      assert(h5_workunit_run_id == workunit_run_id);

      std::string h5_table_name;
      tps->readString("TableName",h5_table_name);
      assert(h5_table_name == fTableParam.fTableName);

      delete tps;

      // READERS -- INFRASTRUCTURE
      fEV                 = fReader->readStruct("event");
      fSC                 = fReader->readStruct("scope");
      fPX                 = fReader->readStruct("pixel");
      fPE                 = fReader->readStruct("pe");

      // READERS -- EVENT
      fREventID           = fEV->readVectorElements<uint32_t>("EventID",Ne);
      fRTargetZenithRad   = 
	fEV->readVectorElements<float>("TargetZenithRad",Ne);
      fRTargetAzimuthRad  = 
	fEV->readVectorElements<float>("TargetAzimuthRad",Ne);
      fRPrimaryZenithRad  = 
	fEV->readVectorElements<float>("PrimaryZenithRad",Ne);
      fRPrimaryAzimuthRad = 
	fEV->readVectorElements<float>("PrimaryAzimuthRad",Ne);
      fRPrimaryCoreEastM  = 
	fEV->readVectorElements<float>("PrimaryCoreEastM",Ne);
      fRPrimaryCoreNorthM = 
	fEV->readVectorElements<float>("PrimaryCoreNorthM",Ne);
      fRPrimaryCoreUpASLM = 
	fEV->readVectorElements<float>("PrimaryCoreUpASLM",Ne);
      fRNumHitScopes      = 
	fEV->readVectorElements<uint16_t>("NumHitScopes",Ne); 
      fRScopeStart        = fEV->readVectorElements<uint32_t>("ScopeStart",Ne);
      fRScopeCount        = fEV->readVectorElements<uint16_t>("ScopeCount",Ne);

      fNEventInWorkunit   = fREventID->rows();

      // READERS -- SCOPE
      fRScopeID           = fSC->readVectorElements<uint16_t>("ScopeID",Ns);
      fRScopeZenithRad    = 
	fSC->readVectorElements<float>("ScopeZenithRad",Ns);
      fRScopeAzimuthRad   = 
	fSC->readVectorElements<float>("ScopeAzimuthRad",Ns);
      fRPixelStart        = fSC->readVectorElements<uint32_t>("PixelStart",Ns);
      fRPixelCount        = fSC->readVectorElements<uint32_t>("PixelCount",Ns);

      // READERS -- PIXEL
      fRPixelID           = fPX->readVectorElements<uint32_t>("PixelID",Nx);
      fRPEStart           = fPX->readVectorElements<uint32_t>("PEStart",Nx);
      fRPECount           = fPX->readVectorElements<uint32_t>("PECount",Nx);

      // READERS -- PE
      fRPixelTimeNS       = fPE->readVectorElements<float>("PixelTimeNS",Np);
    }
  catch(const VSOctaveH5Exception& e)
    {
      std::cerr << "Exception while reading file: " << filename
		<< std::endl;	      
      delete fReader; 
      fReader = NULL;
      throw(e);
    }

  return true;
}

#define DELETEZERO(x) delete x; x=0

void VSEventSupplierOneTableHDF::close()
{
  fIDOfOpenWorkunit   = 0;
  fNEventInWorkunit   = 0;

  DELETEZERO(fREventID);
  DELETEZERO(fRTargetZenithRad);
  DELETEZERO(fRTargetAzimuthRad);
  DELETEZERO(fRPrimaryZenithRad);
  DELETEZERO(fRPrimaryAzimuthRad);
  DELETEZERO(fRPrimaryCoreEastM);
  DELETEZERO(fRPrimaryCoreNorthM);
  DELETEZERO(fRPrimaryCoreUpASLM);
  DELETEZERO(fRNumHitScopes);
  DELETEZERO(fRScopeStart);
  DELETEZERO(fRScopeCount);
  DELETEZERO(fRScopeID);
  DELETEZERO(fRScopeZenithRad);
  DELETEZERO(fRScopeAzimuthRad);
  DELETEZERO(fRPixelStart);
  DELETEZERO(fRPixelCount);
  DELETEZERO(fRPixelID);
  DELETEZERO(fRPEStart);
  DELETEZERO(fRPECount);
  DELETEZERO(fRPixelTimeNS);

  DELETEZERO(fEV);
  DELETEZERO(fSC);
  DELETEZERO(fPX);
  DELETEZERO(fPE);

  DELETEZERO(fReader);
}
