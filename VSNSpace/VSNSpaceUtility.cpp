/*! \file VSNSpaceUtility.cpp

  Utility class for creating nspace filters.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/08/2007

*/

#include "VSNSpaceUtility.hpp"
#include "VSNSpaceOctaveH5IO.hpp"

using namespace VERITAS;

// ============================================================================
// INDEX SINK
// ============================================================================

VSNSpaceIndexSink::~VSNSpaceIndexSink()
{
  // nothing to see here
}

VSNSpaceOctaveIndexSink::
VSNSpaceOctaveIndexSink(const std::string& filename):
  VSNSpaceIndexSink(), m_writer(), m_inspace(), m_index()
{
  const char* fn = "OctaveIndexSink::OctaveIndexSink: ";
  m_writer = 
    new VSOctaveH5Writer(filename,true,"# written by OctaveIndexSink");

  if(m_writer==0)
    {
      std::cerr << fn << "could not write to " << filename << std::endl;
      return;
    }

  m_inspace = m_writer->writeExpandableVector<bool>("in_space");
  if(m_inspace == 0)
    {
      delete m_writer;
      m_writer = 0;
      std::cerr << fn << "could not write vector \"in_space\"" << std::endl;
      return;
    }

  m_index = m_writer->writeExpandableVector<VSNSpace::Index>("index");
  if(m_index == 0)
    {
      delete m_writer;
      m_writer = 0;
      m_inspace = 0;
      std::cerr << fn << "could not write vector \"index\"" << std::endl;
      return;
    }
}

VSNSpaceOctaveIndexSink::~VSNSpaceOctaveIndexSink()
{
  delete m_writer;
}

bool VSNSpaceOctaveIndexSink::putIndex(bool inspace, VSNSpace::Index index)
{
  const char* fn = "VSNSpaceOctaveIndexSink::putIndex: ";
  if(!m_inspace->append(inspace))
    {
      std::cerr << fn << "could not write \"in_space\" entry" << std::endl;
      return false;
    }
  if(!m_index->append(index))
    {
      std::cerr << fn << "could not write \"index\" entry" << std::endl;
      return false;
    }
  return true;
}

// ============================================================================
// EVENT PASS SINK
// ============================================================================

VSNSpaceEventPassedSink::~VSNSpaceEventPassedSink()
{
  // nothing to see here
}

VSNSpaceOctaveEventPassedSink::
VSNSpaceOctaveEventPassedSink(const std::string& filename):
  VSNSpaceEventPassedSink(), m_writer(), m_passed()
{
  const char* fn = 
    "VSNSpaceOctaveEventPassedSink::VSNSpaceOctaveEventPassedSink: ";
  m_writer = 
    new VSOctaveH5Writer(filename,true,"# written by OctaveEventPassedSink");

  if(m_writer==0)
    {
      std::cerr << fn << "could not write to " << filename << std::endl;
      return;
    }

  m_passed = m_writer->writeExpandableVector<bool>("passed");
  if(m_passed == 0)
    {
      delete m_writer;
      m_writer = 0;
      std::cerr << fn << "could not write vector \"passed\"" << std::endl;
      return;
    }
}

VSNSpaceOctaveEventPassedSink::~VSNSpaceOctaveEventPassedSink()
{
  delete m_writer;
}

bool VSNSpaceOctaveEventPassedSink::putEvent(bool passed)
{
  const char* fn = "OctaveEventPassedSink::putEvent: ";
  if(!m_passed->append(passed))
    {
      std::cerr << fn << "could not write \"passed\" entry" << std::endl;
      return false;
    }
  return true;
}

// ============================================================================
// NSPACE UTILITY
// ============================================================================

VSNSpaceUtility::Options 
VSNSpaceUtility::s_default_options = VSNSpaceUtility::Options();

VSNSpaceUtility::VSNSpaceUtility(const Options& opt):
  m_options(opt)
{
  
}

VSNSpaceUtility::~VSNSpaceUtility()
{

}

void VSNSpaceUtility::loadHist(const std::string& file,
			       VSNSpace& hist)
{
  const char* fn = "VSNSpaceUtility::loadHist: ";
  VSNSpaceIO* io = new VSNSpaceOctaveH5IO;
  std::cerr << fn << "loading histogram " << file << std::endl;
  if(!io->readHistogram(file,hist))
    {
      std::cerr 
	<< file << ": could not open histogram" << std::endl;
      exit(EXIT_FAILURE);
    }
  delete io;
}

void VSNSpaceUtility::createSpace(VSNSpaceUtility::SpaceDefinition& space_def,
				  VSNSpace::Space& space)
{
  const char* fn = "VSNSpaceUtility::createSpace: ";
  unsigned ndim = space_def.size();
  std::cerr << fn << "space has " << ndim << " dimensions" << std::endl;

  space.resize(ndim);
  for(unsigned idim=0;idim<ndim;idim++)
    {
      std::string name = space_def[idim].first;
      VSNSpace::Coord lo = space_def[idim].second;
      VSNSpace::Coord hi = space_def[idim].third;
      VSNSpace::Index nbins = space_def[idim].fourth;

      if(lo == hi)
	{
	  std::cerr 
	    << fn << "axis " << idim << ": lower bound equals upper bound"
	    << std::endl;
	  exit(EXIT_FAILURE);
	}

      if(nbins == 0)
	{
	  std::cerr 
	    << fn << "axis " << idim << ": number of bins is zero"
	    << std::endl;
	  exit(EXIT_FAILURE);
	}
      
      if(hi > lo)
	space.axes[idim] = VSNSpace::Axis(lo, hi, nbins, name);
      else
	space.axes[idim] = VSNSpace::Axis(hi, lo, nbins, name);

      std::cerr << fn << "axis " << idim << ": " << space.axes[idim]
		<< std::endl;
    }
}

void VSNSpaceUtility::printSpace(VSNSpace::Space& space)
{
  std::cout << "Dimensions:            " << space.ndim << std::endl
	    << "Total space size:      " << space.size() << std::endl
	    << "Space definition:      ";
  
  for(std::vector<VSNSpace::Axis>::const_iterator iaxis = space.axes.begin();
      iaxis != space.axes.end(); iaxis++)
    {
      if(iaxis != space.axes.begin())
	std::cout << std::endl << "                       ";
      std::cout << *iaxis;
    }
  std::cout << std::endl << std::endl;
}

void VSNSpaceUtility::printHist(VSNSpace& hist)
{
  VSNSpace::Space space = hist.space();

  // Histogram statistics -----------------------------------------------------
  VSNSpace::Weight total_weight = hist.totalWeight();
  VSNSpace::Weight min_weight = hist.minWeight();
  VSNSpace::Weight max_weight = hist.maxWeight();

  VSNSpace::Point moments1(space.ndim);
  hist.integrate(VSNSpaceMoments1Functional(),moments1);

  std::cout << "Maximum cell weight:   " << max_weight << std::endl
	    << "Minimum cell weight:   " << min_weight << std::endl
	    << "Total weight:          " << total_weight << std::endl
	    << "Weight per cell:       " << total_weight/double(space.size())
	    << std::endl;
  
  for(unsigned idim=0;idim<moments1.ndim;idim++)
    std::cout << "Centroid " << idim << ":            " 
	      << moments1.x[idim]/total_weight << std::endl;

  // Occupied cells -----------------------------------------------------------
  VSNSpace::Volume occ_volume;
  hist.volumeAbove(0, occ_volume);
  unsigned nocc = occ_volume.countCellsInVolume();

  std::cout << "Occupied cells:        " << nocc << std::endl
	    << "Weight per occ. cell:  " << total_weight/double(nocc)
	    << std::endl;

  dumpHist(hist);
}

void VSNSpaceUtility::dumpHist(VSNSpace& hist)
{  
  VSNSpace::Space space = hist.space();
  
  // Write 1D hist to standard output -----------------------------------------
  if(space.ndim == 1)
    {
      for(unsigned index = 0; index<hist.space().axes[0].nbin; index++)
	{
	  VSNSpace::Point p;
	  hist.space().midPointOfIndexUnchecked(index,p);
	  std::cout << p.x[0] << ' '
		    << hist.getWeightUnchecked(index) << std::endl;
	}
    }
  else
    {
      VSNSpace::Cell c(space.ndim);

      std::cout << "x = [";
      c.i[0]=0;
      for(c.i[1]=0; c.i[1]<hist.space().axes[1].nbin; c.i[1]++)
	{
	  VSNSpace::Point p;
	  hist.space().midPointOfCellUnchecked(c,p);
	  std::cout << ' ' << p.x[1];
	}
      std::cout << " ]" << std::endl;

      std::cout << "y = [";
      c.i[1]=0;
      for(c.i[0]=0; c.i[0]<hist.space().axes[0].nbin; c.i[0]++)
	{
	  VSNSpace::Point p;
	  hist.space().midPointOfCellUnchecked(c,p);
	  std::cout << ' ' << p.x[0];
	}
      std::cout << " ]" << std::endl;
      
      for(c.i[0]=0; c.i[0]<hist.space().axes[0].nbin; c.i[0]++)
	{
	  for(c.i[1]=0; c.i[1]<hist.space().axes[1].nbin; c.i[1]++)
	    {	      
	      VSNSpace::Weight w = -1; // initialize any value
	      hist.getWeight(c,w);
	      if(c.i[1] != 0)std::cout << ' ';
	      std::cout << w;
	    }
	  std::cout << std::endl;
	}
    }
}

void VSNSpaceUtility::createOrdering(const VSNSpace& on,
				     const VSNSpace& off,
				     VSNSpace::Ordering& ordering)
{
  // Mode dependant processing ------------------------------------------------
  VSNSpace roff(off);
  roff*=m_options.ratio;

  VSNSpace excess(on);
  if(!m_options.sig_pure)excess -= roff;

  roff*=m_options.ratio;
  VSNSpace off_var(roff);

  VSNSpace variance(off_var);

  if(m_options.ordering_mode == "significance_simple")
    variance += on;

  // Produce ordering ---------------------------------------------------------
  ordering.m_space = on.space();

  std::cout << "Mode:      " << m_options.ordering_mode << std::endl;
  std::cout << "Threshold: " << m_options.threshold << std::endl;

  if(!VSNSpace::cumulativeOrderingSimple(ordering.m_index, excess, 
					 variance, m_options.threshold))
    {
      std::cerr << "ordering failed" << std::endl;
      exit(EXIT_FAILURE);
    }

  // Prepare output vectors ---------------------------------------------------
  unsigned nordering = ordering.m_index.size();
  std::vector<VSNSpace::Weight> cumsum_on(nordering);
  std::vector<VSNSpace::Weight> cumsum_off(nordering);
  std::vector<VSNSpace::Weight> cumsum_excess(nordering);
  std::vector<VSNSpace::Weight> cumsum_Q(nordering);
  std::vector<VSNSpace::Weight> cumsum_sigma(nordering);

  VSNSpace::Weight sum_on = 0;
  VSNSpace::Weight sum_off = 0;
  VSNSpace::Weight sum_excess = 0;
  VSNSpace::Weight sum_Q = 0;
  VSNSpace::Weight sum_variance = 0;  

  for(unsigned iordering = 0; iordering<nordering; iordering++)
    {
      unsigned index = ordering.m_index[iordering];
      sum_on       += on.getWeightUnchecked(index);
      sum_off      += off.getWeightUnchecked(index);
      sum_excess   += excess.getWeightUnchecked(index);
      sum_Q        += off_var.getWeightUnchecked(index);
      sum_variance += on.getWeightUnchecked(index)
	+ off_var.getWeightUnchecked(index);
      
      cumsum_on[iordering]       = sum_on;
      cumsum_off[iordering]      = sum_off;
      cumsum_excess[iordering]   = sum_excess;
      cumsum_Q[iordering]        = sum_excess/sqrt(sum_Q);
      cumsum_sigma[iordering]    = sum_excess/sqrt(sum_variance);
    }

  ordering.m_counts_on = cumsum_on;
  ordering.m_counts_off = cumsum_off;
  ordering.m_excess = cumsum_excess;
  ordering.m_Q = cumsum_Q;
  ordering.m_sigma = cumsum_sigma;
}

void VSNSpaceUtility::createFilter(VSNSpace::Ordering& ordering,
				   VSNSpace::Volume& filter)
{
  if(m_options.ncell)
    {
      vsassert(m_options.ncell < ordering.m_index.size());
      ordering.m_index.erase(ordering.m_index.begin()+m_options.ncell,
			     ordering.m_index.end());      
      filter.loadSparse(ordering.m_space,ordering.m_index);
    }
  else
    {
      double excess_max = *std::max_element(ordering.m_excess.begin(),
					    ordering.m_excess.end());

      unsigned ncell = 0;
      for(std::vector< double >::iterator itr = ordering.m_excess.begin(); 
	  itr != ordering.m_excess.end(); itr++)
	{

	  if(*itr > m_options.fraction*excess_max)
	    break;

	  ncell++;
	}
     
      vsassert(ncell < ordering.m_index.size());
      ordering.m_index.erase(ordering.m_index.begin()+ncell,
			     ordering.m_index.end());

      filter.loadSparse(ordering.m_space,ordering.m_index);
    }
}

#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSNSpaceUtility::configure(VSOptions& options, 
				const std::string& profile, 
				const std::string& opt_prefix)
{
  options.findWithValue(OPTNAME(opt_prefix,"ordering_mode"),
			s_default_options.ordering_mode,	
			"Set the ordering mode. Valid modes are: "
			"\"Q_simple\", optimize on Q-value in Gaussian "
			"approximation, \"significance_simple\", on "
			"significance in Gaussian approximation or "
			"\"rate\" on excess.");

  options.findWithValue(OPTNAME(opt_prefix,"ratio"),
			s_default_options.ratio,			
			"Set the ratio of ON to OFF counts.");

  options.findWithValue(OPTNAME(opt_prefix,"threshold"), 
			s_default_options.threshold,
			"Set the threshold counts for histogram in ordering "
			"and when producing filter.");

  options.findBoolValue(OPTNAME(opt_prefix,"sig_pure"), 
			s_default_options.sig_pure, true,
			"ON histogram contains pure (simulated) signal with "
			"no background.");

  options.findWithValue(OPTNAME(opt_prefix,"ncell"), 
			s_default_options.ncell,
			"Set the number of cells to be included in the "
			"filter.");

  options.findWithValue(OPTNAME(opt_prefix,"fraction"), 
			s_default_options.fraction,
			"Set the fraction of the maximum filter.");
}
