//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventSupplierDBTable.cpp

  Input of event from single table from simulations database 

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       07/08/2006
  \note
*/

#include <cassert>

#include "VSEventSupplierDBTable.hpp"

using namespace VERITAS;

VSEventSupplierDBTable::
VSEventSupplierDBTable(const std::string& table_name, VSSimDB* sim_db,
		       uint32_t max_buffer_size)
  : VSEventSupplier(), fSimDB(sim_db), fTableName(table_name),
    fTableParam(), fMaxEventNum(), fEventNumBuffer(), 
    fMaxBufferSize(max_buffer_size),
    fNextEventNum(fEventNumBuffer.end()), fLastEventNum(0), 
    fNCompleteEvent(), fICompleteEvent()
{
  VSSimDBTableParam* table = fSimDB->getDataTableByName(table_name);
  assert(table);

  fMaxEventNum = fSimDB->getMaxEventNumByName(table_name);
  fNCompleteEvent = 
    fSimDB->getCompleteEventCountByName(table_name, fMaxEventNum);

  fTableParam = SimParam(*table,fNCompleteEvent);
  delete table;
}

VSEventSupplierDBTable::~VSEventSupplierDBTable()
{
  // nothing to see here
}

uint32_t VSEventSupplierDBTable::getNextEventNum()
{
  if(fICompleteEvent>=fNCompleteEvent)
    {
      if(!fEventNumBuffer.empty())
	{
	  fEventNumBuffer.clear();
	  fNextEventNum=fEventNumBuffer.end();
	}
      return 0;
    }
  
  fICompleteEvent++;

  if(fNextEventNum==fEventNumBuffer.end())
    {
      fEventNumBuffer.clear();
      while((fEventNumBuffer.empty())&&(fLastEventNum<fMaxEventNum))
	{
	  uint32_t min_event_num = fLastEventNum;
	  uint32_t max_event_num = fLastEventNum+fMaxBufferSize;
	  if(max_event_num>fMaxEventNum)max_event_num=fMaxEventNum;
	  fSimDB->getLimitedCompleteEventNums(fTableName, 
					      min_event_num, max_event_num,
					      fEventNumBuffer);
	  if(fEventNumBuffer.empty())fLastEventNum=max_event_num;
	}
      fNextEventNum=fEventNumBuffer.begin();
    }

  assert(fNextEventNum!=fEventNumBuffer.end());

  fLastEventNum = *fNextEventNum;
  fNextEventNum++;
  return fLastEventNum;
}

std::vector<VSEventSupplier::SimParam> VSEventSupplierDBTable::getSimParam()
{
  std::vector<SimParam> tpvec;
  tpvec.push_back(fTableParam);
  return tpvec;
}

bool VSEventSupplierDBTable::getNextEvent(Event& e)
{
  uint32_t event_id = getNextEventNum();
  if(event_id == 0)return false;

  VSSimDBEventData evd;
  assert(fSimDB->getEventByNum(fTableName, event_id, evd) >= 0);

  e.fEventID                 = event_id;
  e.fTargetZenithRad         = evd.fTargetZenithRad;
  e.fTargetAzimuthRad        = evd.fTargetAzimuthRad;
  e.fPrimarySpeciesCORSIKA   = fTableParam.fPrimaryID;
  e.fPrimaryEnergyTeV        = fTableParam.fEnergyTeV;
  e.fPrimaryZenithRad        = evd.fPrimaryZenithRad;
  e.fPrimaryAzimuthRad       = evd.fPrimaryAzimuthRad;
  e.fPrimaryCoreEastM        = evd.fPrimaryCoreEastM;
  e.fPrimaryCoreNorthM       = evd.fPrimaryCoreNorthM;
  e.fPrimaryCoreUpASLM       = evd.fPrimaryCoreUpASLM;
  e.fSamplingRadiusM         = fTableParam.fSamplingRadiusM;
  e.fTableIndex              = 0;
  e.fEventComplete           = evd.fEventComplete;

  std::vector<VSSimDBScopeData> scopes;
  assert(fSimDB->getAllScopeByNum(fTableName, e.fEventID, scopes) >= 0);

  unsigned nscope = scopes.size();
  e.fScopes.resize(nscope);

  std::vector<std::pair<bool, unsigned> > scope_lookup;
  scope_lookup.reserve(100);

  for(unsigned iscope=0; iscope<nscope; iscope++)
    {
      if(scopes[iscope].fScopeID >= scope_lookup.size())
	scope_lookup.resize(scopes[iscope].fScopeID+1);
      scope_lookup[scopes[iscope].fScopeID].first = true;
      scope_lookup[scopes[iscope].fScopeID].second = iscope;
      e.fScopes[iscope].fScopeID         = scopes[iscope].fScopeID;
      e.fScopes[iscope].fScopeZenithRad  = scopes[iscope].fScopeZenithRad;
      e.fScopes[iscope].fScopeAzimuthRad = scopes[iscope].fScopeAzimuthRad;
      e.fScopes[iscope].fPixels.clear();
      e.fScopes[iscope].fPixels.reserve(1000);
    }

  std::vector<VSSimDBPEData> pes;
  std::vector<std::vector<std::pair<bool, unsigned> > > pix_lookup(nscope);
  assert(fSimDB->getAllPEByNum(fTableName, e.fEventID, pes) >= 0);

  unsigned npe = pes.size();
  for(unsigned ipe=0; ipe<npe; ipe++)
    {
      unsigned iscope = pes[ipe].fScopeID;

      if((iscope >= scope_lookup.size())
	 ||(scope_lookup[iscope].first == false))continue;

      unsigned iscopeslot = scope_lookup[pes[ipe].fScopeID].second;
      unsigned ipix = pes[ipe].fPixelID;
      if(ipix >= pix_lookup[iscopeslot].size())
	pix_lookup[iscopeslot].resize(ipix+1);
      if(pix_lookup[iscopeslot][ipix].first == false)
	{
	  unsigned islot = e.fScopes[iscopeslot].fPixels.size();
	  pix_lookup[iscopeslot][ipix].first = true;
	  pix_lookup[iscopeslot][ipix].second = islot;
	  e.fScopes[iscopeslot].fPixels.resize(islot+1);
	  e.fScopes[iscopeslot].fPixels[islot].fPixelID = ipix;
	  e.fScopes[iscopeslot].fPixels[islot].fPEs.reserve(200);
	}

      unsigned islot = pix_lookup[iscopeslot][ipix].second;

      PE pe;
      pe.fTimeNS = pes[ipe].fPixelTimeNS;

      e.fScopes[iscopeslot].fPixels[islot].fPEs.push_back(pe);
    }
  
  return true;
}

uint32_t VSEventSupplierDBTable::getNumEvents()
{
  return fNCompleteEvent;
}
