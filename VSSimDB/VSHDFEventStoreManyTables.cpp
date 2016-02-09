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

#include <VSFileUtility.hpp>
#include <VSAAlgebra.hpp>
#include <VSHDFEventStoreManyTables.hpp>

using namespace VERITAS;
using namespace VSAAlgebra;

VSHDFEventStoreManyTables::
VSHDFEventStoreManyTables(VSSimDB* simdb): 
  VSSimDBEventStore(), fSimDB(simdb), fEventStore(), fTableParam(), 
  fTableID(), fTableName()
{
  std::vector<VSSimDBTableParam> data_tables = fSimDB->getAllDataTables();

  for(std::vector<VSSimDBTableParam>::iterator itr = 
	data_tables.begin(); itr != data_tables.end(); ++itr)
    fTableParam[itr->fTableID] = *itr;
}

VSHDFEventStoreManyTables::~VSHDFEventStoreManyTables()
{
  for(std::map< uint32_t, VSHDFEventStore* >::iterator itr = 
	fEventStore.begin(); itr != fEventStore.end(); ++itr)
    {
      uint32_t workunit_run_id = itr->second->getWorkunitRunID();
      delete itr->second;
      fSimDB->registerWorkunitRunFinish(workunit_run_id);
    }
}

int VSHDFEventStoreManyTables::insertEvent(VSSimDBEventData& event)
{ 
  unsigned table_id = findTable(event);

  if(table_id != fTableID)
    {
      fTableID = table_id;
      fTableName = fTableParam[table_id].fTableName;
      fSimDB->setInsertTableName(fTableName,true);
    }

  if(fEventStore.find(fTableID) == fEventStore.end()) openFile(fTableID);

  event.fWorkunitRunID = fEventStore[fTableID]->getWorkunitRunID();
  event.fEventComplete = true; 
  return fEventStore[fTableID]->insertEvent(event);
}

int VSHDFEventStoreManyTables::insertScope(const VSSimDBScopeData& scope)
{
  return fEventStore[fTableID]->insertScope(scope);
}

int VSHDFEventStoreManyTables::insertPE(const VSSimDBPEData& pe)
{
  return fEventStore[fTableID]->insertPE(pe);
}

int VSHDFEventStoreManyTables::updateHitScopes(const VSSimDBEventData& event)
{
  return fEventStore[fTableID]->updateHitScopes(event);
}

int VSHDFEventStoreManyTables::updateHitPixels(const VSSimDBScopeData& scope)
{
  return fEventStore[fTableID]->updateHitPixels(scope);
}

void VSHDFEventStoreManyTables::openFile(unsigned table_id)
{
  std::string hostname;
  char hostbuffer[1000];
  gethostname(hostbuffer,1000);
  hostname=hostbuffer;
  unsigned jobid=getppid();

  uint32_t workunit_run_id = 
    fSimDB->registerWorkunitRunStart(table_id,hostname,jobid);

  vsassert(!VSFileUtility::
	   isFile(VSHDFEventStore::nameForWorkunitID(workunit_run_id)));

  fEventStore[table_id] = new VSHDFEventStore(fSimDB,workunit_run_id,
					      &fTableParam[table_id]);
}

unsigned VSHDFEventStoreManyTables::findTable(const VSSimDBEventData& event)
{
  unsigned table_id = 0;
  double dth_min = 0;
  double de_min = 0;

  for(std::map< uint32_t, VSSimDBTableParam >::iterator itr = 
	fTableParam.begin(); itr != fTableParam.end(); ++itr)
    {
      double zn = (itr->second.fZenithMaxRad+itr->second.fZenithMinRad)/2.;
      double az = (itr->second.fAzimuthMaxRad+itr->second.fAzimuthMinRad)/2.;

      if(itr->second.fAzimuthMaxRad < itr->second.fAzimuthMinRad)
	az += M_PI;

      Vec3D e1 = Vec3D::makePolar(zn,az);
      Vec3D e2 = 
	Vec3D::makePolar(event.fPrimaryZenithRad,event.fPrimaryAzimuthRad);

      double dth = acos(e1*e2);
      double de = fabs(log10(itr->second.fEnergyGeV) - 
		       log10(event.fEnergyGeV));

      if(itr == fTableParam.begin() || dth < dth_min || de < de_min)
	{
	  dth_min = dth;
	  de_min = de;
	  table_id = itr->second.fTableID;
	}
    }

  return table_id;
}
