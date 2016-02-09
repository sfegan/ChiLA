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

#include <VSOTelescopeArray.hpp>
#include "VSChannelReMappingDBTable.hpp"

using namespace VERITAS;
using namespace Physics;

VSChannelReMappingDBTable::
VSChannelReMappingDBTable(const std::string& table_name, VSSimDB* sim_db,
			  double field_of_view_deg)
  : VSChannelReMapping(), fScopes(), fNumScopes()
{
  const double fov = field_of_view_deg/2/180*M_PI;

  // --------------------------------------------------------------------------
  // Get optics ID for the table
  // --------------------------------------------------------------------------

  VSSimDBTableParam* table_param;
  table_param = sim_db->getDataTableByName(table_name);
  assert(table_param != 0);
  unsigned optics_id = table_param->fOpticsID;
  delete table_param;

  // --------------------------------------------------------------------------
  // Get the definition of the array from the database
  // --------------------------------------------------------------------------

  VSOTelescopeArray array;
  array.readFromDatabase(sim_db->db(), optics_id);

  fNumScopes = array.numTelescopes();
  fScopes.resize(fNumScopes);
  
  for(unsigned iscope=0;iscope<fNumScopes;iscope++)
    {
      unsigned npix=array.telescope(iscope)->numPixels();
      fScopes[iscope].fNumAnlChannels = 0;

      for(unsigned ipix=0;ipix<npix; ipix++)
	{	  
	  Vec3D r(array.telescope(iscope)->pixel(ipix)
		  ->incomingSkyVectorAtZenith(1.0));
	  
	  double theta = acos(r * Vec3D(0,0,-1));
	  
	  if(theta <= fov)
	    {
	      fScopes[iscope].fNumSimChannels = ipix+1;
	      fScopes[iscope].fAnlMapping.push_back(ipix);
	      fScopes[iscope].fSimMapping.
		push_back(fScopes[iscope].fNumAnlChannels);
	      fScopes[iscope].fNumAnlChannels++;
	      fScopes[iscope].fIsPresentInAnl.push_back(true);
	    }
	  else
	    {
	      fScopes[iscope].fSimMapping.push_back(0);
	      fScopes[iscope].fIsPresentInAnl.push_back(false);
	    }
	}
    }
}

VSChannelReMappingDBTable::
~VSChannelReMappingDBTable()
{
  // nothing to see here
}

unsigned VSChannelReMappingDBTable::numScopes()
{
  return fNumScopes;
}

unsigned VSChannelReMappingDBTable::numAnlChannels(unsigned iscope)
{
  vsassert(iscope < fNumScopes);
  return fScopes[iscope].fNumAnlChannels;
}

unsigned VSChannelReMappingDBTable::numSimChannels(unsigned iscope)
{
  vsassert(iscope < fNumScopes);
  return fScopes[iscope].fNumSimChannels;
}

unsigned VSChannelReMappingDBTable::
mapAnlToSim(unsigned iscope, unsigned ianalysis)
{
  vsassert(iscope < fNumScopes);
  vsassert(ianalysis < fScopes[iscope].fNumAnlChannels);
  return fScopes[iscope].fAnlMapping[ianalysis];
}

unsigned VSChannelReMappingDBTable::
mapSimToAnl(unsigned iscope, unsigned isim)
{
  vsassert(iscope < fNumScopes);
  vsassert(isim < fScopes[iscope].fNumSimChannels);
  return fScopes[iscope].fSimMapping[isim];
}

bool VSChannelReMappingDBTable::
isPresentInAnl(unsigned iscope, unsigned isim)
{
  vsassert(iscope < fNumScopes);
  vsassert(isim < fScopes[iscope].fNumSimChannels);
  return fScopes[iscope].fIsPresentInAnl[isim];
}
