//-*-mode:c++; mode:font-lock;-*-

/*! \file VSMuonAnalysisData.hpp

  Data structures for muon analysis

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/13/2007

  $Id: VSMuonAnalysisData.hpp,v 3.3 2007/12/04 18:05:03 sfegan Exp $

*/

#define MUON_TEST_NON_UNIFORMITY

#ifndef VSMUONANALYSISDATA_HPP
#define VSMUONANALYSISDATA_HPP

#include<VSOctaveIO.hpp>

namespace VERITAS
{

  struct VSMuonAnalysisDatum
  {
    VSMuonAnalysisDatum():
      ievent(),
      chi2(), x0(), y0(), r0(), 
      c_nimage(), c_signal(), c_rms(), c_cx(), c_cy(), c_xi(), c_U0_r0(),
      r_nimage(), r_signal(), r_rms(), r_cx(), r_cy(), r_xi(), r_U0_r0(),
      r_U0_r0_corr()
#ifdef MUON_TEST_NON_UNIFORMITY
      , nu_xi(), nu_U0_r0(), nu_U0_r0_corr()
#endif
    { /* nothing to see here */ }

    // Administrative parameters ----------------------------------------------
    unsigned ievent;

    // Parameters of fing fit, which is done using cleaned image --------------
    double   chi2;
    double   x0;
    double   y0;
    double   r0;

    // Muon parameters from cleaned image -------------------------------------
    unsigned c_nimage;
    double   c_signal;
    double   c_rms;
    double   c_cx;
    double   c_cy;
    double   c_xi;
    double   c_U0_r0;

    // Muon parameters from ring image ----------------------------------------
    unsigned r_nimage;
    double   r_signal;
    double   r_rms;
    double   r_cx;
    double   r_cy;
    double   r_xi;
    double   r_U0_r0;
    double   r_U0_r0_corr;

#ifdef MUON_TEST_NON_UNIFORMITY
    // Test for non uniform response ------------------------------------------
    double   nu_xi;
    double   nu_U0_r0;
    double   nu_U0_r0_corr;
#endif
    
    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,ievent);
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,chi2); 
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,x0); 
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,y0); 
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,r0); 
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,c_nimage);
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,c_signal);
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,c_rms);
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,c_cx);
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,c_cy);
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,c_xi);
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,c_U0_r0);
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,r_nimage);
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,r_signal);
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,r_rms);
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,r_cx);
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,r_cy);
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,r_xi);
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,r_U0_r0);
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,r_U0_r0_corr);
#ifdef MUON_TEST_NON_UNIFORMITY
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,nu_xi);
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,nu_U0_r0);
      H5_ADDMEMBER(c,VSMuonAnalysisDatum,nu_U0_r0_corr);
#endif
    }
    
  };

  class VSMuonAnalysisData
  {
  public:
    VSMuonAnalysisData(unsigned nscope = 0):
      scope(nscope) { /* nothing to see here */ }

    typedef VSMuonAnalysisDatum Datum;
    typedef std::vector<Datum> Scope;

    std::vector<Scope> scope;

    void addDatum(unsigned iscope, const Datum& datum)
    {
      vsassert(iscope < scope.size());
      scope[iscope].push_back(datum);
    }

    void clear() { scope.clear(); }
    void load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer) const;
  };

}

#endif // defined VSMUONANALYSISDATA_HPP
