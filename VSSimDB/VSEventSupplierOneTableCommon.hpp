//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventSupplierOneTableCommon.hpp

  Input of event from single table from simulations dataCommon 

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       07/08/2006
  \note
*/

#ifndef VSEVENTSUPPLIERONETABLECOMMON_HPP
#define VSEVENTSUPPLIERONETABLECOMMON_HPP

#include<vector>

#include<VSSimDB.hpp>
#include<VSEventSupplier.hpp>

//! VERITAS namespace
namespace VERITAS 
{

  class VSEventSupplierOneTableCommon: public VSEventSupplier
  {
  public:    
    VSEventSupplierOneTableCommon(const std::string& table_name, 
				  VSSimDB* sim_db, uint32_t max_buffer_size);

    VSEventSupplierOneTableCommon(const std::string& table_name, 
				  VSSimDB* sim_db,
				  uint32_t nevent_to_supply, 
				  double ievent_to_start_fraction,
				  uint32_t max_buffer_size);

    VSEventSupplierOneTableCommon(const std::string& table_name, 
				  VSSimDB* sim_db);

    VSEventSupplierOneTableCommon(const std::string& table_name, 
				  VSSimDB* sim_db,
				  uint32_t nevent_to_supply, 
				  double ievent_to_start_fraction);

    virtual ~VSEventSupplierOneTableCommon();
    virtual std::vector<SimParam> getSimParam();
    virtual uint32_t getNumEvents();

  protected:
    bool nextEvent();
    uint32_t getNextEventNum();

    uint32_t getICompleteEvent() const { return fICompleteEvent; }
    uint32_t getNEventToSupply() const { return fNEventToSupply; }
    void setStartEventAndCount(const std::string& table_name,
			       uint32_t nevent, double fraction);

    VSSimDB*                         fSimDB;
    std::string                      fTableName;
    SimParam                         fTableParam;

  private:
    VSEventSupplierOneTableCommon(const VSEventSupplierOneTableCommon&);
    VSEventSupplierOneTableCommon& operator=(VSEventSupplierOneTableCommon&);

    void setTableParam(const std::string& table_name);

    uint32_t                         fMaxEventNum;
    std::vector<uint32_t>            fEventNumBuffer;
    unsigned                         fMaxBufferSize;
    std::vector<uint32_t>::iterator  fNextEventNum;
    uint32_t                         fLastEventNum;

    uint32_t                         fNEventToSupply;
    uint32_t                         fZEventToSupply;

    uint32_t                         fICompleteEvent;
  };

}

#endif // defined VSEVENTSUPPLIERONETABLECOMMON_HPP
