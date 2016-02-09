//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventSupplierManyTableschpp

  Dispatch events from multiple suppliers, possibly in random order

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       02/21/2007
  \note
*/

#include<cassert>

#include<VSEventSupplierManyTables.hpp>

using namespace VERITAS;

VSEventSupplierManyTables::
VSEventSupplierManyTables(VSEventSupplierFactory* factory,
			  const std::list<TableName>& tables,
			  RandomNumbers* rng)
  : VSEventSupplier(),
    fNEvent(), fNEventRemaining(), fSuppliers(), fTableParam(), fRNG(rng)
{
  for(std::list<TableName>::const_iterator itable = tables.begin();
      itable != tables.end(); itable++)
    fSuppliers.push_back(factory->newOneTableSupplier(*itable));
  finishConstruction();
}

VSEventSupplierManyTables::
VSEventSupplierManyTables(VSEventSupplierFactory* factory,
			  const std::list<TableNameAndCount>& tables,
			  RandomNumbers* rng)
  : VSEventSupplier(),
    fNEvent(), fNEventRemaining(), fSuppliers(), fTableParam(), fRNG(rng)
{
  for(std::list<TableNameAndCount>::const_iterator itable = tables.begin();
      itable != tables.end(); itable++)
    {
      double f = 0;
      if(fRNG)f=fRNG->Uniform();
      VSEventSupplier* s =
	factory->newOneTableSupplier(itable->first, itable->second, f);
      fSuppliers.push_back(s);
    }
  finishConstruction();
}

void VSEventSupplierManyTables::finishConstruction()
{
  unsigned itable=0;
  unsigned ntable=fSuppliers.size();
  fTableParam.reserve(ntable);
  for(std::list<Supplier>::iterator isupplier = fSuppliers.begin();
      isupplier != fSuppliers.end(); isupplier++)
    {
      assert(isupplier->supplier);
      isupplier->table_index = itable;
      isupplier->count       = isupplier->supplier->getNumEvents();
      fNEvent += isupplier->count;

      std::vector<VSEventSupplier::SimParam> one_param = 
	isupplier->supplier->getSimParam();
      assert(one_param.size() == 1);
      fTableParam.push_back(one_param[0]);

      itable++;
    }

  // Establish precondition that all suppliers have at least one event
  // available.

  std::list<Supplier>::iterator isupplier = fSuppliers.begin();
  while(isupplier != fSuppliers.end())
    {
      if(isupplier->count == 0)
	{
	  delete isupplier->supplier;
	  std::list<Supplier>::iterator nsupplier = isupplier;
	  nsupplier++;
	  fSuppliers.erase(isupplier);
	  isupplier=nsupplier;
	}
      else isupplier++;
    }

  fNEventRemaining = fNEvent;
}

VSEventSupplierManyTables::~VSEventSupplierManyTables()
{
  for(std::list<Supplier>::iterator isupplier = fSuppliers.begin();
      isupplier != fSuppliers.end(); isupplier++)
    delete isupplier->supplier;
}

std::vector<VSEventSupplier::SimParam> VSEventSupplierManyTables::getSimParam()
{
  return fTableParam;
}

bool VSEventSupplierManyTables::getNextEvent(Event& e)
{
  if(fNEventRemaining==0)return false;

  std::list<Supplier>::iterator isupplier = fSuppliers.begin();
#warning TEMPORARY ASSERTION
  assert(isupplier != fSuppliers.end());

  if(fRNG)
    {
      uint32_t ievent = 
	uint32_t(floor(double(fNEventRemaining)*fRNG->Uniform()));
      uint32_t jevent = isupplier->count;
      while(jevent <= ievent)
	{
#warning TEMPORARY ASSERTION
	  isupplier++;
	  assert(isupplier != fSuppliers.end());
	  jevent += isupplier->count;
	}
    }
  
  assert(isupplier->supplier->getNextEvent(e));
  e.fTableIndex = isupplier->table_index;
  isupplier->count--;
  fNEventRemaining--;

  if(isupplier->count == 0)
    {
      delete isupplier->supplier;
      fSuppliers.erase(isupplier);
    }

  return true;
}

uint32_t VSEventSupplierManyTables::getNumEvents()
{
  return fNEvent;
}

