//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventSupplierOneTableCommon.cpp

  Input of event from single table from simulations dataCommon 

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       07/08/2006
  \note
*/

#include <cassert>

#include "VSEventSupplierOneTableCommon.hpp"

using namespace VERITAS;

VSEventSupplierOneTableCommon::
VSEventSupplierOneTableCommon(const std::string& table_name, 
			      VSSimDB* sim_db, uint32_t max_buffer_size)
  : VSEventSupplier(), fSimDB(sim_db), fTableName(table_name),
    fTableParam(), fMaxEventNum(), fEventNumBuffer(), 
    fMaxBufferSize(max_buffer_size),
    fNextEventNum(fEventNumBuffer.end()), fLastEventNum(0), 
    fNEventToSupply(), fZEventToSupply(), fICompleteEvent()
{
  fMaxEventNum = fSimDB->getMaxEventNumByName(table_name);
  fNEventToSupply =
    fSimDB->getCompleteEventCountByName(table_name, fMaxEventNum);
  setTableParam(table_name);
}

VSEventSupplierOneTableCommon::
VSEventSupplierOneTableCommon(const std::string& table_name, 
			      VSSimDB* sim_db,
			      uint32_t nevent_to_supply, 
			      double ievent_to_start_fraction,
			      uint32_t max_buffer_size)
  : VSEventSupplier(), fSimDB(sim_db), fTableName(table_name),
    fTableParam(), fMaxEventNum(), fEventNumBuffer(), 
    fMaxBufferSize(max_buffer_size),
    fNextEventNum(fEventNumBuffer.end()), fLastEventNum(0), 
    fNEventToSupply(), fZEventToSupply(), fICompleteEvent()
{
  fMaxEventNum = fSimDB->getMaxEventNumByName(table_name);
  setStartEventAndCount(table_name,nevent_to_supply,ievent_to_start_fraction);
  setTableParam(table_name);
}

VSEventSupplierOneTableCommon::
VSEventSupplierOneTableCommon(const std::string& table_name, 
			      VSSimDB* sim_db)
  : VSEventSupplier(), fSimDB(sim_db), fTableName(table_name),
    fTableParam(), fMaxEventNum(), fEventNumBuffer(), 
    fMaxBufferSize(), fNextEventNum(fEventNumBuffer.end()), fLastEventNum(0), 
    fNEventToSupply(), fZEventToSupply(), fICompleteEvent()
{
  fNEventToSupply = fSimDB->getCompleteEventCountByName(table_name);
  setTableParam(table_name);
}

VSEventSupplierOneTableCommon::
VSEventSupplierOneTableCommon(const std::string& table_name, 
			      VSSimDB* sim_db,
			      uint32_t nevent_to_supply, 
			      double ievent_to_start_fraction)
  : VSEventSupplier(), fSimDB(sim_db), fTableName(table_name),
    fTableParam(), fMaxEventNum(), fEventNumBuffer(), 
    fMaxBufferSize(), fNextEventNum(fEventNumBuffer.end()), fLastEventNum(0), 
    fNEventToSupply(), fZEventToSupply(), fICompleteEvent()
{
  setStartEventAndCount(table_name,nevent_to_supply,ievent_to_start_fraction);
  setTableParam(table_name);
}

void VSEventSupplierOneTableCommon::
setTableParam(const std::string& table_name)
{
  VSSimDBTableParam* table = fSimDB->getDataTableByName(table_name);
  assert(table);
  fTableParam = SimParam(*table,fNEventToSupply);
  delete table;
}

void VSEventSupplierOneTableCommon::
setStartEventAndCount(const std::string& table_name,
		      uint32_t nevent, double fraction)
{
  uint32_t nevent_complete =
    fSimDB->getCompleteEventCountByName(table_name, fMaxEventNum);

  if(nevent < nevent_complete)
    {
      fNEventToSupply = nevent;
      if(fraction<=0.0)
	fICompleteEvent = 0;
      else if(fraction>=1.0)
	fICompleteEvent = nevent_complete-fNEventToSupply;
      else
	{
	  double ievent = double(nevent_complete-fNEventToSupply+1)*fraction;
	  fICompleteEvent = uint32_t(ievent);
#warning TEMPORARY ASSERTION
	  assert(fICompleteEvent <= (nevent_complete-fNEventToSupply));
	}
      fZEventToSupply = fICompleteEvent;
    }
  else
    {
      fNEventToSupply = nevent_complete;
      fZEventToSupply = 0;
    }
}

VSEventSupplierOneTableCommon::~VSEventSupplierOneTableCommon()
{
  // nothing to see here
}

bool VSEventSupplierOneTableCommon::nextEvent()
{
  if(fICompleteEvent>=fNEventToSupply+fZEventToSupply)return false;
  fICompleteEvent += 1;
  return true;
}

uint32_t VSEventSupplierOneTableCommon::getNextEventNum()
{
  if(!nextEvent())
    {
      if(!fEventNumBuffer.empty())
	{
	  fEventNumBuffer.clear();
	  fNextEventNum=fEventNumBuffer.end();
	}
      return 0;
    }

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

std::vector<VSEventSupplier::SimParam> 
VSEventSupplierOneTableCommon::getSimParam()
{
  std::vector<SimParam> tpvec;
  tpvec.push_back(fTableParam);
  return tpvec;
}

uint32_t VSEventSupplierOneTableCommon::getNumEvents()
{
  return fNEventToSupply;
}
