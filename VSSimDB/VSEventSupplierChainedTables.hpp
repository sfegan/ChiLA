//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventSupplierChainedTables.hpp

  Input of events from database tables in a chained configuration

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       10/26/2006
  \note
*/

#ifndef VSEVENTSUPPLIERCHAINEDTABLE_HPP
#define VSEVENTSUPPLIERCHAINEDTABLE_HPP

#include<cmath>
#include<list>
#include<vector>

#include<VSSimDB.hpp>
#include<VSEventSupplier.hpp>

//! VERITAS namespace
namespace VERITAS 
{

  class VSEventSupplierChainedTables: public VSEventSupplier
  {
  public:    
    VSEventSupplierChainedTables(VSSimDB* sim_db, 
				 double max_zn_rad=M_PI/2, 
				 double min_zn_rad=0);
    virtual ~VSEventSupplierChainedTables();
    virtual std::vector<SimParam> getSimParam();
    virtual bool getNextEvent(Event& e);
    virtual uint32_t getNumEvents();
  private:
    struct Table
    {
      Table(): param(), supplier() { /* nothing to see here */ }
      VSSimDBTableParam param;
      VSEventSupplier* supplier;
    };

    VSSimDB*                       fSimDB;
    unsigned                       fITable;
    std::list<Table>               fTables;
    std::vector<SimParam>          fTableParam;
    unsigned                       fNEvent;
  };

}

#endif // defined VSEVENTSUPPLIERCHAINEDTABLES_HPP
