#include <VSREventVisitor.hpp>

using namespace VERITAS;

VSREventVisitor::VSREventVisitor()
{
  m_file = new TFile("test.root","recreate");
  m_event_tree = new TTree("t1","tree");

  //  m_event_tree->Branch("run",&run_id,"run_id/I");
  m_event_tree->Branch("rec_mask",
		       &m_evt.used_in_reconstruction_mask,"rec_mask/I");
  m_event_tree->Branch("msc_width",&m_evt.msc_width,"msc_width/D");
  m_event_tree->Branch("msc_length",&m_evt.msc_length,"msc_length/D");
  m_event_tree->Branch("msc_disp",&m_evt.msc_disp,"msc_disp/D");
  m_event_tree->Branch("N2",&m_evt.N2,"N2/D");
  m_event_tree->Branch("R",&m_evt.R,"R/D");
  m_event_tree->Branch("Rx",&m_evt.Rx,"Rx/D");
  m_event_tree->Branch("Ry",&m_evt.Ry,"Ry/D");
  m_event_tree->Branch("zn",&m_evt.zn,"zn/D");
  m_event_tree->Branch("az",&m_evt.az,"az/D");
  m_event_tree->Branch("fovx_derotated",
		       &m_evt.mean_derotated_fov_x,"fovx_derotated/D");
  m_event_tree->Branch("fovy_derotated",
		       &m_evt.mean_derotated_fov_y,"fovy_derotated/D");
  m_event_tree->Branch("fovx",&m_evt.mean_fov_x,"fovx/D");
  m_event_tree->Branch("fovy",&m_evt.mean_fov_y,"fovy/D");
  m_event_tree->Branch("array_zn",&m_evt.mean_array_zn,"array_zn/D");
  m_event_tree->Branch("array_az",&m_evt.mean_array_az,"array_az/D");
  m_event_tree->Branch("theta0",&m_evt.theta0,"theta0/D");
}

VSREventVisitor::~VSREventVisitor()
{
  m_event_tree->Write();

  delete m_file;
}

void VSREventVisitor::visitEvent(const VSEventArrayDatum& event)
{
  m_evt = event;
}

void VSREventVisitor::leaveEvent()
{
  m_event_tree->Fill();
}
