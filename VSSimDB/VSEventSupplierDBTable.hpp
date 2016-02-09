//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventSupplierDBTable.hpp

  Input of event from single table from simulations database 

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       07/08/2006
  \note
*/

#ifndef VSEVENTSUPPLIERDBTABLE_HPP
#define VSEVENTSUPPLIERDBTABLE_HPP

#include<vector>

#include<VSSimDB.hpp>
#include<VSEventSupplier.hpp>

//! VERITAS namespace
namespace VERITAS 
{

  class VSEventSupplierDBTable: public VSEventSupplier
  {
  public:    
    VSEventSupplierDBTable(const std::string& table_name, VSSimDB* sim_db,
			   uint32_t max_buffer_size = 16*1024);
    virtual ~VSEventSupplierDBTable();
    virtual std::vector<SimParam> getSimParam();
    virtual bool getNextEvent(Event& e);
    virtual uint32_t getNumEvents();

  protected:
    uint32_t getNextEventNum();

    VSSimDB*                         fSimDB;
    std::string                      fTableName;
    SimParam                         fTableParam;

  private:
    uint32_t                         fMaxEventNum;
    std::vector<uint32_t>            fEventNumBuffer;
    unsigned                         fMaxBufferSize;
    std::vector<uint32_t>::iterator  fNextEventNum;
    uint32_t                         fLastEventNum;

    uint32_t                         fNCompleteEvent;
    uint32_t                         fICompleteEvent;
  };

}

#endif // defined VSEVENTSUPPLIERDBTABLE_HPP
