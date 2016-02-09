//-*-mode:c++; mode:font-lock;-*-

/*! \file VSChannelReMappingDBTable.hpp

  Class to provide remapping of channels from database optics

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       09/12/2006
  \note
*/

#ifndef VSCHANNELREMAPPINGDBTABLE_HPP
#define VSCHANNELREMAPPINGDBTABLE_HPP

#include<vector>
#include<stdint.h>

#include<VSSimDB.hpp>
#include<VSChannelReMapping.hpp>

//! VERITAS namespace
namespace VERITAS 
{

  class VSChannelReMappingDBTable: public VSChannelReMapping
  {
  public:
    struct Scope
    {
      uint32_t               fNumAnlChannels;
      uint32_t               fNumSimChannels;
      std::vector<uint32_t>  fAnlMapping;
      std::vector<uint32_t>  fSimMapping;
      std::vector<bool>      fIsPresentInAnl;      
    };

    VSChannelReMappingDBTable(const std::string& table_name, VSSimDB* sim_db,
			      double field_of_view_deg);
    virtual ~VSChannelReMappingDBTable();
    virtual unsigned numScopes();
    virtual unsigned numAnlChannels(unsigned iscope);
    virtual unsigned numSimChannels(unsigned iscope);
    virtual unsigned mapAnlToSim(unsigned iscope, unsigned ianalysis);
    virtual unsigned mapSimToAnl(unsigned iscope, unsigned isim);
    virtual bool isPresentInAnl(unsigned iscope, unsigned isim);

  private:
    std::vector< Scope >        fScopes;
    uint32_t                    fNumScopes;
  };

};

#endif // defined VSCHANNELREMAPPINGDBTABLE_HPP
