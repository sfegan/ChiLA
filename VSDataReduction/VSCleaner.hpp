//-*-mode:c++; mode:font-lock;-*-

/*! \file VSCleaner.hpp

  Various classes to clean images

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/10/2006

  $Id: VSCleaner.hpp,v 3.1 2009/06/19 04:35:00 matthew Exp $

*/

#ifndef VSCLEANER_HPP
#define VSCLEANER_HPP

#include<vector>
#include <WhippleCams.h>

namespace VERITAS
{

  class VSCleaner
  {
  public:
    VSCleaner() { /* nothing to see here */ }
    virtual ~VSCleaner();
    virtual unsigned clean(unsigned nchan,
			   unsigned* mask, const double* signal) = 0;
  private:
    VSCleaner(const VSCleaner&);
    VSCleaner& operator= (const VSCleaner&);
  };

  class VSPicBndCleaner: public VSCleaner
  {
  public:
    VSPicBndCleaner(double plevel, double blevel,
		    unsigned nchan, const int (*neighbors)[NUM_NEIGHBORS]):
      VSCleaner(), 
      m_nchan(nchan), m_neighbors(neighbors),
      m_plevel(plevel), m_blevel(blevel) { }
    virtual ~VSPicBndCleaner();
    virtual unsigned clean(unsigned nchan,
			   unsigned* mask, const double* signal);
  private:
    VSPicBndCleaner(const VSPicBndCleaner&);
    VSPicBndCleaner& operator= (const VSPicBndCleaner&);

    enum CleaningState { CS_MASKED, CS_UNDETERMINED, CS_PIC, CS_BND };
    unsigned m_nchan;
    const int (*m_neighbors)[NUM_NEIGHBORS];
    double m_plevel;
    double m_blevel;
  };

  class VSRegionalCleaner: public VSCleaner
  {
  public:
    VSRegionalCleaner(double rlevel, unsigned rsize, double ilevel,
		      unsigned nchan, const int (*neighbors)[NUM_NEIGHBORS]):
      VSCleaner(), 
      m_nchan(nchan), m_neighbors(neighbors),
      m_rlevel(rlevel), m_rsize(rsize), m_ilevel(ilevel) { }
    virtual ~VSRegionalCleaner();
    virtual unsigned clean(unsigned nchan,
			   unsigned* mask, const double* signal);
  private:
    VSRegionalCleaner(const VSRegionalCleaner&);
    VSRegionalCleaner& operator= (const VSRegionalCleaner&);

    enum CleaningState { CS_MASKED, CS_UNDETERMINED, CS_REGION, CS_IMMEDIATE };
    unsigned m_nchan;
    const int (*m_neighbors)[NUM_NEIGHBORS];
    double m_rlevel;
    unsigned m_rsize;
    double m_ilevel;
  };


  class VSCleanerFactory
  {
  public:
    static VSCleaner* getCleaner(const std::vector<std::string>& args,
				 unsigned nchan, 
				 const int (*neighbors)[NUM_NEIGHBORS]);
  private:
    VSCleanerFactory(const VSCleanerFactory&);
    VSCleanerFactory& operator= (const VSCleanerFactory&);
  };

} // namespace VERITAS

#endif
