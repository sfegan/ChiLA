//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventDataVisitor.hpp

  Base class visitor for stage2 event and simulation data.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       05/14/2007

  $Id: VSEventDataVisitor.hpp,v 3.23 2010/06/20 00:52:15 matthew Exp $

*/

#ifndef VSEVENTDATAVISITOR_HPP
#define VSEVENTDATAVISITOR_HPP

#include <vector>

#include <SphericalCoords.h>
#include <VSAAlgebra.hpp>
#include <VSEventData.hpp>
#include <VSSimulationData.hpp>
#include <VSMergedCalibrationData.hpp>
#include <VSDiagnosticsData.hpp>
#include <VSAnalysisStage1.hpp>
#include <VSCutsEvaluator.hpp>
#include <VSScaledParameterLibrary.hpp>
#include <VSResultsSimData.hpp>
#include <VSSimEnergyWeightCalc.hpp>
#include <VSSourceInjector.hpp>
#include <VSSpectrumCalc.hpp>

namespace VERITAS
{
  /////////////////////////////////////////////////////////////////////////////
  //! Base visitor class for processing stage2 data.  Classes which
  //! perform analysis on stage2 data files should inherit from this
  //! class.  Visitor classes are intended to be used in conjunction
  //! with VSEventDataDispatcher.
  /////////////////////////////////////////////////////////////////////////////
  class VSEventDataVisitor
  {
  public:
    VSEventDataVisitor();
    virtual ~VSEventDataVisitor();

    //! Return the set of simulation parameters which should be loaded
    //! to VSArraySimulationDatum.
    virtual void getSimSet(VSEventDataReader::MemberSubset& s) const { }
    virtual void getArraySet(VSEventDataReader::MemberSubset& s) const { }
    virtual void getScopeSet(VSEventDataReader::MemberSubset& s) const { }

    //! Process data structures associated with the run.  This method
    //! is called first when processing a data file.
    //! @param stage1 Stage1 analysis data structure.
    //! @param obs Observations data structure.
    virtual void visitRun(const VSAnalysisStage1Data& stage1,
			  const VSTargetTable::Observation& obs,
			  const VSArrayMergedCalibrationData& cal);
    virtual void leaveRun();

    virtual void visitDiagnostics(const VSArrayDiagnosticsData& diag);
    virtual void leaveDiagnostics();

    virtual void visitEvent(const VSEventArrayDatum& event);
    virtual void leaveEvent();

    virtual void visitScope(const VSEventScopeDatum& scope, unsigned iscope);
    virtual void leaveScope();

    virtual void visitSimEvent(const VSArraySimulationDatum& sim);
    virtual void leaveSimEvent();

    virtual void visitSimHeader(const VSHeaderSimulationDatum& header);
    virtual void leaveSimHeader();
  };

  /////////////////////////////////////////////////////////////////////////////
  //! Class responsible for reading stage2 data files and dispatching
  //! visitor classes inheriting from VSEventDataVisitor.
  /////////////////////////////////////////////////////////////////////////////
  class VSEventDataDispatcher
  {
  public:
    VSEventDataDispatcher();
    VSEventDataDispatcher(VSEventDataVisitor* visitor);

    //! Set the visitor class which will be dispatched when processing
    //! stage2 files.
    void setVisitor(VSEventDataVisitor* visitor) { m_visitor = visitor; }

    //! Process the supplied list of stage2 files.
    //! @param file_list std::list of paths to stage2 files. 
    //! @param nevents Number of events to process.  All events are
    //! processed by default or when nevents = 0.
    //! @param verbose Print detailed information about each run as it
    //! is opened.
    void processFiles(const std::list<std::string>& file_list,
		      unsigned nevents = 0, bool verbose = false);

    //! Process the internal list of stage2 files.
    //! @param nevents Number of events to process.  All events are processed
    //! by default or when nevents = 0. 
    //! @param verbose Print detailed information about each run as it
    //! is opened.
    void processFiles(unsigned nevents = 0, bool verbose = false);

    void loadList(const std::list<std::string>& file_list);

    const std::vector<std::string>& getFiles() const
    { return m_file_list; }

  private:
    std::vector<std::string> m_file_list;

    VSEventDataVisitor*    m_visitor;
  };

}

#endif // VSEVENTDATAVISITOR_HPP
