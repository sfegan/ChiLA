//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventSupplierManyTables.hpp

  Dispatch events from multiple suppliers, possibly in random order

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       02/21/2007
  \note
*/

#ifndef VSEVENTSUPPLIERMANYTABLE_HPP
#define VSEVENTSUPPLIERMANYTABLE_HPP

#include<list>
#include<vector>

#include<RandomNumbers.hpp>

#include<VSEventSupplier.hpp>
#include<VSEventSupplierFactory.hpp>

//! VERITAS namespace
namespace VERITAS 
{

  class VSEventSupplierManyTables: public VSEventSupplier
  {
  public:    

    typedef std::string TableName;
    typedef std::pair<TableName,uint32_t> TableNameAndCount;

    VSEventSupplierManyTables(VSEventSupplierFactory* factory,
			      const std::list<TableName>& tables,
			      RandomNumbers* rng = 0);

    VSEventSupplierManyTables(VSEventSupplierFactory* factory,
			      const std::list<TableNameAndCount>& tables,
			      RandomNumbers* rng = 0);

    virtual ~VSEventSupplierManyTables();
    virtual std::vector<SimParam> getSimParam();
    virtual bool getNextEvent(Event& e);
    virtual uint32_t getNumEvents();

  private:
    void finishConstruction();

    struct Supplier
    {
      Supplier(): supplier(), table_index(), count() { }
      Supplier(VSEventSupplier* _supplier): 
	supplier(_supplier), table_index(), count() { }
      VSEventSupplier* supplier;
      uint32_t table_index;
      uint32_t count;
    };

    unsigned                       fNEvent;
    unsigned                       fNEventRemaining;
    std::list<Supplier>            fSuppliers;
    std::vector<SimParam>          fTableParam;
    RandomNumbers*                 fRNG;
  };

}

#endif // defined VSEVENTSUPPLIERMANYTABLES_HPP
