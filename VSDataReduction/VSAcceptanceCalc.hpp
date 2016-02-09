//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAcceptanceCalc.hpp

  Class for generating a CR acceptance model.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.8 $
  \date       07/29/2007

  $Id: VSAcceptanceCalc.hpp,v 3.8 2010/05/31 02:20:22 matthew Exp $

*/

#ifndef VSACCEPTANCECALC_HPP
#define VSACCEPTANCECALC_HPP

#include <vector>
#include <VSResultsData.hpp>
#include <VSAFunction.hpp>
#include <VSDataModel.hpp>

namespace VERITAS
{
  class VSAcceptanceCalc
  {
  public:

    VSAcceptanceCalc(const std::string bkgnd_model);
    ~VSAcceptanceCalc();

    VSAcceptanceData* fitAcceptance(VSAnalysisStage3Data& data);

  protected:
    VSAcceptanceModel*                  m_bkgnd_model;
  };
}

#endif // VSACCEPTANCECALC_HPP
