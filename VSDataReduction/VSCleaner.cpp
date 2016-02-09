//-*-mode:c++; mode:font-lock;-*-

/*! \file VSCleaner.cpp

  Various classes to supply pointing information

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/10/2006

  $Id: VSCleaner.cpp,v 3.1 2007/12/04 18:05:03 sfegan Exp $

*/

#include<vsassert>

#include<WhippleCams.h>
#include<VSCleaner.hpp>
#include<VSDataConverter.hpp>
#include<fast_alloc.hpp>

using namespace VERITAS;

// ============================================================================
// VSCLEANER - BASE CLASS
// ============================================================================

VSCleaner::~VSCleaner()
{
  // nothing to see here
}

// ============================================================================
// VSPICBNDCLEANER - STANDARD WHIPPLE PICTURE BOUNDARY CLEANING
// ============================================================================

VSPicBndCleaner::~VSPicBndCleaner()
{
  // nothing to see here
}

unsigned VSPicBndCleaner::clean(unsigned nchan,
				unsigned* mask, const double* signal)
{
  unsigned nimage = 0;

  vsassert(nchan>=m_nchan);

  CleaningState* state = FASTCALLOC(CleaningState,m_nchan);

  for(unsigned ichan=0;ichan<m_nchan;ichan++)
    state[ichan] = (mask[ichan]==1)?CS_MASKED:CS_UNDETERMINED;
  
  for(unsigned ichan=0;ichan<m_nchan;ichan++)
    {
      if(state[ichan]==CS_MASKED)continue;
      if(signal[ichan] >= m_plevel)
        {
	  if(state[ichan]==CS_UNDETERMINED)nimage++;
          state[ichan]=CS_PIC;
          for(unsigned jneighbor=0;jneighbor<NUM_NEIGHBORS;jneighbor++)
            {
              int m_nchan=m_neighbors[ichan][jneighbor];
	      if(m_nchan==-1)break;
              if(state[m_nchan]==CS_MASKED)continue;
              if((state[m_nchan]==CS_UNDETERMINED)&&(signal[m_nchan]>=m_blevel))
                nimage++, state[m_nchan]=CS_BND;
            }
        }
      else if(state[ichan]==CS_UNDETERMINED)state[ichan]=CS_MASKED;
    }

  for(unsigned ichan=0;ichan<m_nchan;ichan++)
    mask[ichan]=(state[ichan] == CS_MASKED)?1:0;
  for(unsigned ichan=m_nchan;ichan<nchan;ichan++)
    mask[ichan]=1;

  FASTFREE(state);

  return nimage;
}

// ============================================================================
// VSREGIONALCLEANER - REGIONAL (ISLAND) CLEANING METHOD
// ============================================================================

VSRegionalCleaner::~VSRegionalCleaner()
{
  // nothing to see here
}

unsigned VSRegionalCleaner::clean(unsigned nchan,
				  unsigned* mask, const double* signal)
{
  vsassert(nchan>=m_nchan);

  CleaningState* state = FASTCALLOC(CleaningState,m_nchan);

  for(unsigned int ichan=0;ichan<m_nchan;ichan++)
    state[ichan] = (mask[ichan]==1)?CS_MASKED:CS_UNDETERMINED;

  // Go through all the tubes and assign them to a "region" of continuous
  // tubes above a threshold of r_level. The "0" region is a pretend one for
  // channels that are below the r_level threshold. Hence the initialisation
  // of region_size with 1... which skips the "0" entry in the vector

  unsigned* region_size = FASTCALLOC(unsigned,m_nchan);
  unsigned* channel_stack = FASTCALLOC(unsigned,m_nchan);
  unsigned* channel_region = FASTCALLOC(unsigned,m_nchan);

  unsigned iregion = 0;
  region_size[iregion] = 0; // For tubes that are not part of a region

  for(unsigned ichan=0;ichan<m_nchan;ichan++)
    {
      if(state[ichan] != CS_UNDETERMINED)continue;
      unsigned channel_stack_occupancy = 0;
      if(signal[ichan] >= m_rlevel)
	{
	  region_size[++iregion] = 0;
	  channel_stack[channel_stack_occupancy++] = ichan;
	  state[ichan] = CS_REGION;
	}
      else
	{
	  state[ichan] = CS_MASKED;
	  continue;
	}

      while(channel_stack_occupancy)
	{
	  unsigned jchan = channel_stack[--channel_stack_occupancy];
	  region_size[iregion]++;
	  channel_region[jchan]=iregion;

	  if(signal[jchan] >= m_ilevel)state[jchan] = CS_IMMEDIATE;
#warning try assert here
	  else state[jchan] = CS_REGION;
	  
          for(unsigned jneighbor=0;jneighbor<NUM_NEIGHBORS;jneighbor++)
            {
              int m_nchan=m_neighbors[jchan][jneighbor];
	      if(m_nchan==-1)break;
              if(state[m_nchan]==CS_MASKED)continue;
              if(state[m_nchan]==CS_UNDETERMINED)
		{
		  if(signal[m_nchan] >= m_rlevel)
		    channel_stack[channel_stack_occupancy++] = m_nchan,
		      state[m_nchan] = CS_REGION;
		  else
		    state[m_nchan] = CS_MASKED;
		}
	    }
	}
    }

  unsigned nimage = 0;
  for(unsigned ichan=0;ichan<m_nchan;ichan++)
    {
      if((state[ichan] == CS_IMMEDIATE)||
	 ((state[ichan] == CS_REGION)
	  &&(region_size[channel_region[ichan]] >= m_rsize)))
	mask[ichan] = 0, nimage++;
      else
	mask[ichan] = 1;
    }      

  for(unsigned ichan=m_nchan;ichan<nchan;ichan++)
    mask[ichan]=1;
  
  FASTFREE(region_size);
  FASTFREE(channel_stack);
  FASTFREE(channel_region);
  FASTFREE(state);

  return nimage;
}

// ============================================================================
// CLEANERFACTORY
// ============================================================================

VSCleaner* VSCleanerFactory::
getCleaner(const std::vector<std::string>& args,
	   unsigned nchan, const int (*neighbors)[NUM_NEIGHBORS])
{
  if(args.size()==0)return 0;
  else if(args[0] == "picbnd")
    {
      double piclevel = 4.25;
      double bndlevel = 2.25;
      if(args.size() >= 2)VSDataConverter::fromString(piclevel,args[1]);
      if(args.size() >= 3)VSDataConverter::fromString(bndlevel,args[2]);
      return new VSPicBndCleaner(piclevel, bndlevel, nchan, neighbors);
    }
  else if(args[0] == "regional")
    {
      double immlevel = 5.0;
      double reglevel = 2.0;
      unsigned regsize = 3;
      if(args.size() >= 2)VSDataConverter::fromString(reglevel,args[1]);
      if(args.size() >= 3)VSDataConverter::fromString(regsize,args[2]);
      if(args.size() >= 4)VSDataConverter::fromString(immlevel,args[3]);
      return new VSRegionalCleaner(reglevel, regsize, immlevel,
				   nchan, neighbors);

    }
  return 0;
}
