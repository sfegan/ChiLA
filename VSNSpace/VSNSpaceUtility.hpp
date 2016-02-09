/*! \file VSNSpaceUtility.hpp

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

#ifndef VSNSPACEUTILITY_HPP
#define VSNSPACEUTILITY_HPP

#include <vector>
#include <string>

#include <VSDataConverter.hpp>
#include <VSNSpace.hpp>
#include <VSOctaveH5Reader.hpp>
#include <VSNSpaceDataSource.hpp>
#include <VSOptions.hpp>

namespace VERITAS
{
  // ==========================================================================
  // INDEX SINK
  // ==========================================================================

  class VSNSpaceIndexSink
  {
  public:
    virtual ~VSNSpaceIndexSink();
    virtual bool putIndex(bool inspace, VSNSpace::Index index) = 0;
  };

  class VSNSpaceOctaveIndexSink: public VSNSpaceIndexSink
  {
  public:
    VSNSpaceOctaveIndexSink(const std::string& filename);
    virtual ~VSNSpaceOctaveIndexSink();
    virtual bool putIndex(bool inspace, VSNSpace::Index index);
  private:
    VSOctaveH5Writer*                            m_writer;
    VSOctaveH5WriterVector<bool>*                m_inspace;
    VSOctaveH5WriterVector<VSNSpace::Index>*     m_index;
  };

  // ==========================================================================
  // EVENT PASS SINK
  // ==========================================================================

  class VSNSpaceEventPassedSink
  {
  public:
    virtual ~VSNSpaceEventPassedSink();
    virtual bool putEvent(bool passed) = 0;
  };

  class VSNSpaceOctaveEventPassedSink: public VSNSpaceEventPassedSink
  {
  public:
    VSNSpaceOctaveEventPassedSink(const std::string& filename);
    virtual ~VSNSpaceOctaveEventPassedSink();
    virtual bool putEvent(bool passed);
  private:
    VSOctaveH5Writer*                            m_writer;
    VSOctaveH5WriterVector<bool>*                m_passed;
  };

  // ==========================================================================
  // MOMENTS FUNCTIONAL
  // ==========================================================================

  class VSNSpaceMoments1Functional
  {
  public:
    void operator() (const VSNSpace::Point& p, VSNSpace::Point& x) const 
    { x=p; }
  };

  // ==========================================================================
  // NSPACE UTILITY
  // ==========================================================================

  class VSNSpaceUtility
  {
  public:

    struct Options
    {
      Options():
	ordering_mode("Q_simple"),
	ncell(0),
	fraction(1.0),
	ratio(1.0),
	threshold(0),
	sig_pure(true)
      { }

      std::string              ordering_mode;
      unsigned                 ncell;
      double                   fraction;
      double                   ratio;
      double                   threshold;
      bool                     sig_pure;
    };

    typedef 
    quad<std::string, VSNSpace::Coord, VSNSpace::Coord, VSNSpace::Index> 
    AxisDefinition;

    typedef std::vector< AxisDefinition > SpaceDefinition;
    VSNSpaceUtility(const Options& opt = s_default_options);
    ~VSNSpaceUtility();

    static void createSpace(SpaceDefinition& space_def,
			    VSNSpace::Space& space);

    static void printSpace(VSNSpace::Space& space);

    static void printHist(VSNSpace& hist);

    static void dumpHist(VSNSpace& hist);

    static void loadHist(const std::string& file, VSNSpace& hist);

    static void projectHist(VSNSpace& hist,
			    VSNSpace::Volume& volume,
			    const std::set<unsigned>& dims);

    void createOrdering(const VSNSpace& on, const VSNSpace& off,
			VSNSpace::Ordering& ordering);

    void createFilter(VSNSpace::Ordering& ordering,
		      VSNSpace::Volume& filter);

    static void configure(VSOptions& options, 
			  const std::string& profile="",
			  const std::string& opt_prefix="");

  private:
    Options                        m_options;

    static Options                 s_default_options;
  };

}

#endif // VSNSPACEUTILITY_HPP
