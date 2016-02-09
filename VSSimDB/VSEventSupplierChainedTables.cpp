#include <VSEventSupplierDBTable.hpp>
#include <VSEventSupplierChainedTables.hpp>

using namespace VERITAS;

VSEventSupplierChainedTables::
VSEventSupplierChainedTables(VSSimDB* sim_db, 
			     double max_zn_rad, double min_zn_rad)
  : fSimDB(sim_db), fITable(0), fTables(), fTableParam(), fNEvent()
{
  std::vector<VSSimDBTableParam> tables = fSimDB->getAllDataTables();
  for(std::vector<VSSimDBTableParam>::const_iterator itable=tables.begin();
      itable!=tables.end(); itable++)
    if((itable->fZenithMinRad<=max_zn_rad)
       &&(itable->fZenithMaxRad>=min_zn_rad))
      {
	Table t;
	t.param = *itable;
	t.supplier = new VSEventSupplierDBTable(t.param.fTableName, sim_db);
	uint32_t nevent = t.supplier->getNumEvents();
	fNEvent += nevent;
	fTables.push_back(t);
	fTableParam.push_back(SimParam(*itable,nevent));
      }
}

VSEventSupplierChainedTables::~VSEventSupplierChainedTables()
{
  while(!fTables.empty())
    {
      delete fTables.front().supplier;
      fTables.pop_front();
    }
}

std::vector<VSEventSupplier::SimParam> VSEventSupplierChainedTables::
getSimParam()
{
  return fTableParam;
}

bool VSEventSupplierChainedTables::getNextEvent(Event& e)
{
  while(!fTables.empty())
    {
      if(fTables.front().supplier->getNextEvent(e))
	{
	  e.fTableIndex = fITable;
	  return true;
	}
      delete fTables.front().supplier;
      fTables.pop_front();
      fITable++;
    }
  return false;
}

uint32_t VSEventSupplierChainedTables::getNumEvents()
{
  return fNEvent;
}

