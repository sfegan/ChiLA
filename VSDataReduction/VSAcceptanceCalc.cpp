//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAcceptanceCalc.cpp

  Class for generating a CR acceptance model.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.9 $
  \date       07/29/2007

  $Id: VSAcceptanceCalc.cpp,v 3.9 2010/05/31 02:20:22 matthew Exp $

*/

#include <VSAcceptanceCalc.hpp>
#include <VSAMath.hpp>
#include <VSALinearLeastSquares.hpp>
#include <VSANonlinearFitting.hpp>
#include <VSAQuadrature.hpp>
#include <VSBkgndModel.hpp>

using namespace VERITAS;
using namespace VERITAS::VSAFunction;

// ----------------------------------------------------------------------------
// VSAcceptanceCalc
// ----------------------------------------------------------------------------
VSAcceptanceCalc::VSAcceptanceCalc(const std::string bkgnd_model)
{ 
  std::vector< std::string > bkgnd_param;
  VSDataConverter::fromString(bkgnd_param, bkgnd_model);

  vsassert(bkgnd_param.size() >=1);

  std::string model = bkgnd_param[0];
  bkgnd_param.erase(bkgnd_param.begin());

  std::string param = VSDataConverter::toString(bkgnd_param);

  if(model == "bessel2")
    m_bkgnd_model = new VSAcceptanceModelBessel2(param,
						 bin_size_deg, offset_max);
  else if(model == "poly")
    m_bkgnd_model = new VSAcceptanceModelPoly(param,
					      bin_size_deg, offset_max);
  else
    {
      std::cerr << "VSIntegralAnalysis(): Unknown background model: " 
		<< bkgnd_model << std::endl;
      exit(EXIT_FAILURE);
    }
}

VSAcceptanceCalc::~VSAcceptanceCalc()
{
  delete m_bkgnd_model;
}

VSAcceptanceData* VSIntegralAnalysis::
fitAcceptanceCamera(VSAnalysisStage3Data& data)
{
  VSAMath::Data<VSModelCoord> fit_data;

  const unsigned nptg = data.ptg_xy.size();
  std::vector< VSAAlgebra::Vec2D > ptg_xy = data.ptg_xy;
  std::vector< VSSimple2DHist<double,double> > exposure_hists(nptg);
  std::vector< VSSimple2DHist<double,double> > sky_counts_hists(nptg);

  for(std::vector< VSAnalysisStage3Data::RunData >::const_iterator itr =
	data.run_data().begin(); itr != data.run_data().end(); ++itr)
    {
      exposure_hists[itr->iptg()].merge(itr->sky_exposure_hist);
      sky_counts_hists[itr->iptg()].merge(itr->sky_counts_hist);
    }

  // Construct the data model -------------------------------------------------
  m_bkgnd_model->clear();
  for(unsigned iptg = 0; iptg < nptg; iptg++)
    {
      m_bkgnd_model->setObs(ptg_xy[iptg]);

      double zcount = 0;
      for(VSSimple2DHist<double,double>::iterator itr = 
	    cam_counts_hists[iptg].begin(); itr !=
	    cam_counts_hists[iptg].end(); ++itr)
	{      
	  VSAAlgebra::Vec2D xy(itr->x(),itr->y());

	  if(exposure_hists[iptg].countForIndex(itr->bin())==0) continue;
	  else if(m_exclusion_region.isExcluded(xy)) continue;
	  zcount+=itr->count();
	}

      // Skip this pointing if all bins are empty -----------------------------
      if(zcount == 0) continue;

      // Create the dataset ---------------------------------------------------
      for(VSSimple2DHist<double,double>::iterator itr = 
	    sky_counts_hists[iptg].begin(); itr !=
	    sky_counts_hists[iptg].end(); ++itr)
	{      
	  VSAAlgebra::Vec2D xy(itr->x(),itr->y());

	  if(exposure_hists[iptg].countForIndex(itr->bin())==0) continue;
	  else if(m_exclusion_region.isExcluded(xy)) continue;
	
	  double z = itr->count();
	  VSModelCoord c(itr->x(),itr->y(),iptg);
	  VSAMath::DataPoint<VSModelCoord> p(c,z,sqrt(z));
	  fit_data.insert(p);
	}

      std::cout << iptg << " " << ptg_xy[iptg].x() << " "
		<< ptg_xy[iptg].y() << " " 
		<< fit_data.size() << " "
		<< zcount << std::endl;
    }

  VSAcceptanceData* d = m_bkgnd_model->fit(fit_data);
  return d;
}
