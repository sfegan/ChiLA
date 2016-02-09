#ifndef VSREVENTVISITOR_HPP
#define VSREVENTVISITOR_HPP

// ----------------------------------------------------------------------------
// ROOT Includes
// ----------------------------------------------------------------------------
#include <TROOT.h>
#include <TH1F.h>
#include <TRint.h>
#include <TPDF.h>
#include <TStyle.h>
#include <TH2F.h>
#include <TF1.h>
#include <THStack.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TTree.h>


#include <VSEventDataVisitor.hpp>

namespace VERITAS
{
  class VSREventVisitor : public VSEventDataVisitor
  {
  public:

    VSREventVisitor();
    virtual ~VSREventVisitor();

    virtual void getSimSet(VSEventDataReader::MemberSubset& s) const
    { 
      s.insert("energy_tev");
      s.insert("core_east_m");
      s.insert("core_north_m");
      s.insert("core_elevation_asl_m");
      s.insert("primary_azimuth_deg");
      s.insert("primary_zenith_deg");
    }

    virtual void getScopeSet(VSEventDataReader::MemberSubset& s) const
    { 
      s.insert("used_in_reconstruction");
      s.insert("R");
      s.insert("N");
      s.insert("G");
      s.insert("theta1");
      s.insert("fp_disp");
      s.insert("fp_dist");
      s.insert("lambdad");
    }

    virtual void getArraySet(VSEventDataReader::MemberSubset& s) const
    { 
      // s.insert("event_num");
      // s.insert("scope");
      // s.insert("theta1");
    }
    
    virtual void visitEvent(const VSEventArrayDatum& event);
    virtual void leaveEvent();

  private:

    VSEventArrayDatum m_evt;

    TTree* m_event_tree;
    TFile* m_file;
  };

}

#endif // VSREVENTVISITOR_HPP
