//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventSupplierOneTableDB.hpp

  Input of event from single table from simulations database 

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       07/08/2006
  \note
*/

#ifndef VSEVENTSUPPLIERONETABLEDB_HPP
#define VSEVENTSUPPLIERONETABLEDB_HPP

#include<vector>

#include<VSSimDB.hpp>
#include<VSEventSupplier.hpp>
#include<VSEventSupplierOneTableCommon.hpp>

//! VERITAS namespace
namespace VERITAS 
{

  class VSEventSupplierOneTableDB: public VSEventSupplierOneTableCommon
  {
  public:    
    VSEventSupplierOneTableDB(const std::string& table_name, VSSimDB* sim_db,
			      uint32_t max_buffer_size = 16*1024);

    VSEventSupplierOneTableDB(const std::string& table_name, VSSimDB* sim_db,
			      uint32_t nevent_to_supply, 
			      double ievent_to_start_fraction,
			      uint32_t max_buffer_size = 16*1024);

    virtual ~VSEventSupplierOneTableDB();
    virtual bool getNextEvent(Event& e);
  };

}

#endif // defined VSEVENTSUPPLIERONETABLEDB_HPP
