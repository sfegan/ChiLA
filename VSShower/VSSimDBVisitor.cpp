//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimDBVisitor.hpp
  CORSIKA visitor stores ray-traced PEs in the database

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       07/26/2005
*/

#include <VSSimDBVisitor.hpp>

using namespace VERITAS;

static inline float d2r(const float& x)
{
  return x/180.0*M_PI;
}

VSSimDBVisitor::
VSSimDBVisitor(uint32_t workunit_runid, VSSimDBEventStore* simdb,
	       bool store_empty_telescopes,
	       VSOTelescopeArray* array, RandomNumbers* rng,
	       VSTargeting* targeting, VSArrayTracking* tracking,
	       VSPrimaryArrivalDistribution* pad,
	       VSSimDBWavelengthData* quaneff_data,
	       VSSimDBWavelengthData* mirrref_data,
	       VSSimDBWavelengthAltitudeData* atmoabs_data):
  VSRayTracingVisitor(array,rng,targeting,tracking,pad,quaneff_data,
		      mirrref_data,atmoabs_data), 
  fWorkunitRunID(workunit_runid), fDB(simdb),
  fStoreEmptyTelescopes(store_empty_telescopes), fArray(array),
  fDBEventData(), fDBScopeData(), fDBPEData(),
  fQuanEffData(quaneff_data), fMirrRefData(mirrref_data),
  fAtmoAbsData(atmoabs_data),
  fHitScope(false), fHitPixels(), fExtraHitPixels()
{
  // nothing to see here
}

VSSimDBVisitor::~VSSimDBVisitor()
{
  // nothing to see here
}

void VSSimDBVisitor::visitEvent(const VSCORSIKAEvent& event, bool& veto)
{
  VSRayTracingVisitor::visitEvent(event, veto);
  if(!veto)
    {     
      fDBEventData.fEnergyGeV = event.pri_energy*1E3;
    }
}

// Event as seen by sampled array
void VSSimDBVisitor::visitEventUse(const VSCORSIKAEventUse& use, bool& veto)
{
  VSRayTracingVisitor::visitEventUse(use, veto);
  if(!veto)
    {
      fDBEventData.fEventID           = 0;  // Event ID  in by DB call
      fDBEventData.fWorkunitRunID     = fWorkunitRunID;      
      fDBEventData.fTargetZenithRad   = fTargetZenithRad;
      fDBEventData.fTargetAzimuthRad  = fTargetAzimuthRad;
      fDBEventData.fPrimaryZenithRad  = d2r(90.0-fEvent.pri_elevation);
      fDBEventData.fPrimaryAzimuthRad = d2r(fEvent.pri_azimuth);
      fDBEventData.fPrimaryCoreEastM  = -use.pri_ycore;
      fDBEventData.fPrimaryCoreNorthM = use.pri_xcore;
      fDBEventData.fPrimaryCoreUpASLM = fObsLevel/100.0; // stupid units!!
      fDBEventData.fNumHitScopes      = 0;
      fDBEventData.fEventComplete     = false;
      fDB->insertEvent(fDBEventData);

      fDBScopeData.fEventID           = fDBEventData.fEventID;
      fDBPEData.fEventID              = fDBEventData.fEventID;
    }
}

void VSSimDBVisitor::leaveEventUse(bool veto)
{
  if(!veto)
    {
      fDB->updateHitScopes(fDBEventData);
    }
}

void VSSimDBVisitor::visitTelescopeEvent(const VSCORSIKATelescopeEvent& scope,
					 bool& veto)
{
  VSRayTracingVisitor::visitTelescopeEvent(scope, veto);
  if(!veto)
    {
      if(fStoreEmptyTelescopes)
	{
	  const VSOTelescope* _scope = fArray->telescope(scope.scope_num);
	  //fDBScopeData.fEventID         = fDBEventData.fEventID;
	  fDBScopeData.fScopeID           = _scope->id();
	  fDBScopeData.fScopeZenithRad    = M_PI_2-_scope->elevation();
	  fDBScopeData.fScopeAzimuthRad   = _scope->azimuth();
	  fDBScopeData.fNumHitPixels      = 0;

	  fDB->insertScope(fDBScopeData);
	}
      fHitScope=false;
    }
}

void VSSimDBVisitor::leaveTelescopeEvent(bool veto)
{
  if((!veto)&&((fStoreEmptyTelescopes)||(fHitScope)))
    {
      fDB->updateHitPixels(fDBScopeData);
    }
}

void VSSimDBVisitor::visitPE(const VSORayTracer::TraceInfo& info)
{
  if(!fHitScope)
    {
      const VSOTelescope* scope = info.scope;

      fDBEventData.fNumHitScopes++;
      
      if(!fStoreEmptyTelescopes)
	{
	  //fDBScopeData.fEventID         = fDBEventData.fEventID;
	  fDBScopeData.fScopeID           = scope->id();
	  fDBScopeData.fScopeZenithRad    = M_PI_2-scope->elevation();
	  fDBScopeData.fScopeAzimuthRad   = scope->azimuth();
	  fDBScopeData.fNumHitPixels      = 0;

	  fDB->insertScope(fDBScopeData);
	}

      fDBPEData.fScopeID              = fDBScopeData.fScopeID;

      if(fHitPixels.size() < scope->numPixels())
	fHitPixels.resize(scope->numPixels());
      for(unsigned ipix=0;ipix<scope->numPixels();ipix++)
	fHitPixels[ipix]=false;
      fExtraHitPixels.clear();

      fHitScope=true;
    }

  // Hack to allow camera faking
  if((info.pixel->id() < fHitPixels.size()))
    {
      if(fHitPixels[info.pixel->id()] == false)
	{
	  fHitPixels[info.pixel->id()] = true;
	  fDBScopeData.fNumHitPixels++;
	}
    }
  else
    {
      if(fExtraHitPixels.find(info.pixel->id()) == fExtraHitPixels.end())
	{
	  fExtraHitPixels.insert(info.pixel->id());
	  fDBScopeData.fNumHitPixels++;	  
	}
    }

  //fDBPEData.fEventID     = fDBEventData.fEventID;
  //fDBPEData.fScopeID     = fDBScopeData.fScopeID;
  fDBPEData.fPixelID       = info.pixel->id();
  fDBPEData.fPixelTimeNS   = 1e9 * info.fplane_t;
  for(unsigned ipe = 0; ipe < fNPE; ipe++) fDB->insertPE(fDBPEData);
}


#ifdef TEST_MAIN

#include <VSOptions.hpp>
#include <VSDBFactory.hpp>

#include "RandomNumbers.hpp"

int main(int argc, char** argv)
{
  VSOptions options(argc,argv);
  VSDBFactory::configure(&options);

  char *progname = *argv;
  argv++, argc--;

  // --------------------------------------------------------------------------
  // Get the database and table names from the command line arguments
  // --------------------------------------------------------------------------
  if(argc<3)
    {
      std::cerr << "Usage: " << progname << " database optics_id filename" 
		<< std::endl;
      exit(EXIT_FAILURE);
    }

  std::string database(*argv);
  argc--; argv++;

  unsigned optics_id = 0;
  VSDataConverter::fromString(optics_id, *argv);
  argc--; argv++;

  std::string filename(*argv);
  argc--; argv++;

  // --------------------------------------------------------------------------
  // Create the database
  // --------------------------------------------------------------------------
  
  VSDatabase* db = VSDBFactory::getInstance()->createVSDB();
  db->createDatabase(database,VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
  db->useDatabase(database);

  // --------------------------------------------------------------------------
  // Get the definition of the array from the database
  // --------------------------------------------------------------------------

  VSOTelescopeArray array;
  array.readFromDatabase(db, optics_id);
 
  // --------------------------------------------------------------------------
  // Simulation database
  // --------------------------------------------------------------------------
  VSSimDB sim_db(db);
  sim_db.createInfrastructureTables();

  // --------------------------------------------------------------------------
  // Create a table and prepare for INSERTs
  // --------------------------------------------------------------------------
  VSSimDBTableParam table;
  table.fPrimaryID          = 1;
  table.fEnergyGeV          = 10.0;
  table.fZenithMinRad       = d2r(0);
  table.fZenithMaxRad       = d2r(20);
  table.fAzimuthMinRad      = d2r(0);
  table.fAzimuthMaxRad      = d2r(360);
  table.fSamplingRadiusM    = 300.0;
  table.fTargetEventCount   = 1000;
  table.fWorkUnitEventCount = 100;
  table.fOpticsID           = 0;
  table.fTableName          = sim_db.tableName(table);
  sim_db.createDataTables(table);

  sim_db.setInsertTableName(table.fTableName);

  // --------------------------------------------------------------------------
  // Set up the CORSIKA visitor
  // --------------------------------------------------------------------------

  RandomNumbers rng("/tmp/seeds.dat");

  VSParallelArrayTracking manager;
  VSSimDBVisitor visitor(&array, &rng, &manager, &sim_db);
  VSCORSIKAEventDispatcher dispatcher(&visitor);

  // --------------------------------------------------------------------------
  // INSERT the data
  // --------------------------------------------------------------------------
  dispatcher.processFile(filename.c_str());
}

#endif
