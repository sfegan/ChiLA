//-*-mode:c++; mode:font-lock;-*-

/*! \file dump_corsika.cpp

  Dump CORSIKA output file to console

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    1.0
  \date       03/02/2005
*/

#include <iostream>
#include <vector>
#include <cmath>

#include <VSCORSIKAEvent.hpp>
#include <VSNSpace.hpp>
#include <VSGenSpace.hpp>

using namespace VERITAS;

class VSDensityHistGenerator: public VSCORSIKAEventVisitor
{
public:
  VSDensityHistGenerator();
  virtual  ~VSDensityHistGenerator();

  virtual void visitTelescopeSpec(const VSCORSIKATelescopeSpec& scopespec);
  virtual void visitEvent(const VSCORSIKAEvent& event, bool& veto);
  virtual void visitEventUse(const VSCORSIKAEventUse& use, bool& veto);
  virtual void visitTelescopeEvent(const VSCORSIKATelescopeEvent& scope,
				   bool& veto);
  virtual void visitPhotonBunch(const VSCORSIKAPhotonBunch& bunch);
private:
  unsigned                            m_nevent;
  VSNSpace::Space                     m_space;
  VSNSpace*                           m_event_density;
  VSNSpace*                           m_total_density;
  VSNSpace*                           m_square_density;
  VSGenSpace<OneSidedIntervalWeight>* m_density;
  unsigned                            m_iscope;
  std::vector<double>                 m_sx;
  std::vector<double>                 m_sy;
};

VSDensityHistGenerator::VSDensityHistGenerator():
  VSCORSIKAEventVisitor(), 
  m_nevent(), m_space(1), 
  m_event_density(), m_total_density(), m_square_density(), m_density()
{
  m_space.axes[0] = VSNSpace::Axis(0, 500*500, 15*15, 0, "Distance squared");
  m_event_density = new VSNSpace(m_space);
  m_total_density = new VSNSpace(m_space);
  m_square_density = new VSNSpace(m_space);
  m_density = new VSGenSpace<OneSidedIntervalWeight>(m_space);
}

VSDensityHistGenerator::~VSDensityHistGenerator()
{
  const double darea = M_PI/m_space.axes[0].bin_factor;
  unsigned npoint = m_space.size();
  for(unsigned ipoint=0; ipoint<npoint; ipoint++)
    {
      VSNSpace::Point p;
      m_space.midPointOfIndexUnchecked(ipoint,p);
      double med = m_density->getWeightUnchecked(ipoint)/darea;
      double mean = 
	m_total_density->getWeightUnchecked(ipoint)/(darea*m_nevent);
      double var = 
	m_square_density->getWeightUnchecked(ipoint)/(darea*darea*m_nevent) -
	mean*mean;

      std::cout << p.x[0] << ' ' << sqrt(p.x[0]) << ' '
		<< med << ' ' << mean << ' ' << sqrt(var) << '\n';
    }
  delete m_density;
  delete m_total_density;
  delete m_square_density;
  delete m_event_density;
}

void VSDensityHistGenerator::
visitTelescopeSpec(const VSCORSIKATelescopeSpec& scopespec)
{
  if(scopespec.scope_num >= m_sx.size())
    {
      m_sx.resize(scopespec.scope_num+1);
      m_sy.resize(scopespec.scope_num+1);
    }

  m_sx[scopespec.scope_num] = scopespec.scope_x;
  m_sy[scopespec.scope_num] = scopespec.scope_y;
}

void VSDensityHistGenerator::
visitEvent(const VSCORSIKAEvent& event, bool& veto)
{
  std::cerr << event.event_num << std::endl;
}

void VSDensityHistGenerator::
visitEventUse(const VSCORSIKAEventUse& use, bool& veto)
{
  m_nevent++;
  *m_density += *m_event_density;
  *m_total_density += *m_event_density;
  *m_event_density *= *m_event_density;
  *m_square_density += *m_event_density;
  m_event_density->clear();
}

void VSDensityHistGenerator::
visitTelescopeEvent(const VSCORSIKATelescopeEvent& scope, bool& veto)
{
  m_iscope = scope.scope_num;
}

void VSDensityHistGenerator::
visitPhotonBunch(const VSCORSIKAPhotonBunch& bunch)
{
  VSNSpace::Point p = m_space.point();
  double ix = bunch.ph_rel_y + m_sx[m_iscope];
  double iy = bunch.ph_rel_x + m_sy[m_iscope];
  p.x[0] = ix*ix + iy*iy;
  m_event_density->accumulate(p,fabs(bunch.ph_count));
}

int main(int argc, char**argv)
{
  //char* progname = *argv;
  argc--,argv++;

  VSDensityHistGenerator visitor;
  VSCORSIKAEventDispatcher dispatcher(&visitor);

  while(argc)
    {
      dispatcher.processFile(*argv);
      argc--,argv++;
    }
}
