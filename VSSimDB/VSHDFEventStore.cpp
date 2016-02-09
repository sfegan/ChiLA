//-*-mode:c++; mode:font-lock;-*-

/*! \file VSHDFEventStore.cpp

  Classes for writing event data to HDF5 files with DB assistance

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       22/10/2006
  \note
*/

#include <iomanip>
#include <sstream>
#include <string>

#include <VSHDFEventStore.hpp>

using namespace VERITAS;

VSHDFEventStore::
VSHDFEventStore(VSSimDB* simdb, uint32_t workunit_run_id,
		const VSSimDBTableParam* table_param,
		VSOctaveH5WriterStruct* writer):
  VSSimDBEventStore(),
  fSimDB(simdb), fWorkunitRunID(workunit_run_id), fWriter(writer),
  fMyWriter(),
  fScopeCount(), fScopeCache(),
  fEV(), fSC(), fPX(), fPE(),
  fWEventID(), fWTargetZenithRad(), fWTargetAzimuthRad(), 
  fWPrimaryZenithRad(), fWPrimaryAzimuthRad(), 
  fWPrimaryCoreEastM(), fWPrimaryCoreNorthM(), fWPrimaryCoreUpASLM(),
  fWNumHitScopes(), fWScopeStart(), 
  fWScopeCount(), fWScopeID(), fWScopeZenithRad(), fWScopeAzimuthRad(), 
  /* fWNumHitPixels(), */ fWPixelStart(), fWPixelCount(),
  fWPixelID(), fWPEStart(), fWPECount(),
  fWPixelTimeNS()
{
  if(fWriter == 0)
    {
      fMyWriter = new VSOctaveH5Writer(nameForWorkunitID(workunit_run_id));
      fWriter = fMyWriter;
    }

  fWriter->writeScalar("workunit_run_id",workunit_run_id);

  // WRITE THE TABLE PARAMETERS
  VSOctaveH5WriterStruct* tps = fWriter->writeStruct("table_param");
  tps->writeScalar("TableID", table_param->fTableID);
  tps->writeScalar("PrimaryID", table_param->fPrimaryID);
  tps->writeScalar("EnergyGeV", table_param->fEnergyGeV);
  tps->writeScalar("ZenithMinRad", table_param->fZenithMinRad);
  tps->writeScalar("ZenithMaxRad", table_param->fZenithMaxRad);
  tps->writeScalar("AzimuthMinRad", table_param->fAzimuthMinRad);
  tps->writeScalar("AzimuthMaxRad", table_param->fAzimuthMaxRad);
  tps->writeScalar("OpticsID", table_param->fOpticsID);
  tps->writeScalar("SamplingRadiusM", table_param->fSamplingRadiusM);
  tps->writeScalar("TargetEventCount", table_param->fTargetEventCount);
  tps->writeScalar("WorkunitEventCount", table_param->fWorkunitEventCount);
  tps->writeScalar("WorkunitPriority", table_param->fWorkunitPriority);
  tps->writeString("TableName", table_param->fTableName);
  delete tps;

  // WRITERS -- INFRASTRUCTURE
  fEV = fWriter->writeStruct("event");
  fSC = fWriter->writeStruct("scope");
  fPX = fWriter->writeStruct("pixel");
  fPE = fWriter->writeStruct("pe");

  // WRITERS -- EVENT
  fWEventID           = fEV->writeExpandableVector<uint32_t>("EventID");
  fWTargetZenithRad   = fEV->writeExpandableVector<float>("TargetZenithRad");
  fWTargetAzimuthRad  = fEV->writeExpandableVector<float>("TargetAzimuthRad");
  fWPrimaryZenithRad  = fEV->writeExpandableVector<float>("PrimaryZenithRad");
  fWPrimaryAzimuthRad = fEV->writeExpandableVector<float>("PrimaryAzimuthRad");
  fWPrimaryCoreEastM  = fEV->writeExpandableVector<float>("PrimaryCoreEastM");
  fWPrimaryCoreNorthM = fEV->writeExpandableVector<float>("PrimaryCoreNorthM");
  fWPrimaryCoreUpASLM = fEV->writeExpandableVector<float>("PrimaryCoreUpASLM");
  fWNumHitScopes      = fEV->writeExpandableVector<uint16_t>("NumHitScopes"); 
  fWScopeStart        = fEV->writeExpandableVector<uint32_t>("ScopeStart");
  fWScopeCount        = fEV->writeExpandableVector<uint16_t>("ScopeCount");

  // WRITERS -- SCOPE
  fWScopeID           = fSC->writeExpandableVector<uint16_t>("ScopeID");
  fWScopeZenithRad    = fSC->writeExpandableVector<float>("ScopeZenithRad");
  fWScopeAzimuthRad   = fSC->writeExpandableVector<float>("ScopeAzimuthRad");
  //fWNumHitPixels      = fSC->writeExpandableVector<uint32_t>("NumHitPixels");
  fWPixelStart        = fSC->writeExpandableVector<uint32_t>("PixelStart");
  fWPixelCount        = fSC->writeExpandableVector<uint32_t>("PixelCount");

  // WRITERS -- PIXEL
  fWPixelID           = fPX->writeExpandableVector<uint32_t>("PixelID");
  fWPEStart           = fPX->writeExpandableVector<uint32_t>("PEStart");
  fWPECount           = fPX->writeExpandableVector<uint32_t>("PECount");

  // WRITERS -- PE
  fWPixelTimeNS       = fPE->writeExpandableVector<float>("PixelTimeNS");
}

VSHDFEventStore::~VSHDFEventStore()
{
  // Not strictly necessary to delete these as they would be deleted
  // by the writer anyway... but no harm!

  // WRITERS -- EVENT
  delete fWEventID;
  delete fWTargetZenithRad;
  delete fWTargetAzimuthRad;
  delete fWPrimaryZenithRad;
  delete fWPrimaryAzimuthRad;
  delete fWPrimaryCoreEastM;
  delete fWPrimaryCoreNorthM;
  delete fWPrimaryCoreUpASLM;
  delete fWNumHitScopes;
  delete fWScopeStart;
  delete fWScopeCount;

  // WRITERS -- SCOPE
  delete fWScopeID;
  delete fWScopeZenithRad;
  delete fWScopeAzimuthRad;
  //delete fWNumHitPixels;
  delete fWPixelStart;
  delete fWPixelCount;
  
  // WRITERS -- PIXEL
  delete fWPixelID;
  delete fWPEStart;
  delete fWPECount;
  
  // WRITERS -- PE
  delete fWPixelTimeNS;

  // WRITERS -- INFRASTRUCTURE
  delete fEV;
  delete fSC;
  delete fPX;
  delete fPE;

  delete fMyWriter;
}

int VSHDFEventStore::insertEvent(VSSimDBEventData& event)
{ 
  // Insert event into database - get consecutive event ID
  int ret = fSimDB->insertEvent(event);
  if(ret != 1)return ret;

  fWEventID->append(event.fEventID);
  fWTargetZenithRad->append(event.fTargetZenithRad);
  fWTargetAzimuthRad->append(event.fTargetAzimuthRad);
  fWPrimaryZenithRad->append(event.fPrimaryZenithRad);
  fWPrimaryAzimuthRad->append(event.fPrimaryAzimuthRad);
  fWPrimaryCoreEastM->append(event.fPrimaryCoreEastM);
  fWPrimaryCoreNorthM->append(event.fPrimaryCoreNorthM);
  fWPrimaryCoreUpASLM->append(event.fPrimaryCoreUpASLM);
  fWScopeStart->append(fWScopeID->rows());
  fScopeCount=0;

  return ret;
}

int VSHDFEventStore::insertScope(const VSSimDBScopeData& scope)
{
  fScopeCount++;
  assert(fScopeCache.empty());

  fWScopeID->append(scope.fScopeID);
  fWScopeZenithRad->append(scope.fScopeZenithRad);
  fWScopeAzimuthRad->append(scope.fScopeAzimuthRad);
  fWPixelStart->append(fWPixelID->rows());  

  return 1;
}

int VSHDFEventStore::insertPE(const VSSimDBPEData& pe)
{
  ChanCache*& chan_cache(fScopeCache[pe.fPixelID]);
  if(!chan_cache)chan_cache = new ChanCache;
  chan_cache->PEs.insert(pe.fPixelTimeNS);

  return 1;
}

int VSHDFEventStore::updateHitScopes(const VSSimDBEventData& event)
{
  // fSimDB->updateHitScopes(event);
  fWNumHitScopes->append(event.fNumHitScopes);
  fWScopeCount->append(fScopeCount);

  return 1;
}

int VSHDFEventStore::updateHitPixels(const VSSimDBScopeData& scope)
{
  unsigned npixel = 0;
  for(std::map<uint32_t, ChanCache*>::iterator ipixel=fScopeCache.begin();
      ipixel!=fScopeCache.end();ipixel++)
    {
      npixel++;
      
      fWPixelID->append(ipixel->first);
      fWPEStart->append(fWPixelTimeNS->rows());

      unsigned npe = 0;
      ChanCache* chan = ipixel->second;
      for(std::multiset<float>::const_iterator ipe = chan->PEs.begin();
	  ipe != chan->PEs.end(); ipe++)
	{
	  fWPixelTimeNS->append(*ipe);
	  npe++;
	}
      delete chan;

      fWPECount->append(npe);
    }
  fScopeCache.clear();

  //fWNumHitPixels->append(scope.fNumHitPixels);
  fWPixelCount->append(npixel);

  return 1;
}

std::string VSHDFEventStore::nameForWorkunitID(unsigned workunit_run_id)
{
  std::ostringstream stream;
  stream << "workunit" << std::setw(10) << std::setfill('0') 
	 << workunit_run_id << ".h5";
  return stream.str();
}
