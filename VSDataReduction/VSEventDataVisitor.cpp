//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventDataVisitor.cpp
  

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       03/17/2005
*/

#include "VSEventDataVisitor.hpp"
#include "VSStage3SimCalc.hpp"
#include <VSFileUtility.hpp>
#include <VSLineTokenizer.hpp>
#include <VSFileLock.hpp>
#include <VSSimpleStat.hpp>

using namespace VERITAS;

// ============================================================================
// VSEventDataVisitor
// ============================================================================

VSEventDataVisitor::VSEventDataVisitor()
{

}

VSEventDataVisitor::~VSEventDataVisitor()
{

}

void VSEventDataVisitor::visitRun(const VSAnalysisStage1Data& stage1,
				  const VSTargetTable::Observation& obs,
				  const VSArrayMergedCalibrationData& cal)
{

}

void VSEventDataVisitor::leaveRun()
{

}

void VSEventDataVisitor::visitDiagnostics(const VSArrayDiagnosticsData& diag)
{

}

void VSEventDataVisitor::leaveDiagnostics()
{

}

void VSEventDataVisitor::visitEvent(const VSEventArrayDatum& event)
{

}

void VSEventDataVisitor::leaveEvent()
{

}

void VSEventDataVisitor::visitScope(const VSEventScopeDatum& scope,
				    unsigned iscope)
{

}

void VSEventDataVisitor::leaveScope()
{

}

void VSEventDataVisitor::visitSimEvent(const VSArraySimulationDatum& sim)
{

}

void VSEventDataVisitor::leaveSimEvent()
{

}

void VSEventDataVisitor::visitSimHeader(const VSHeaderSimulationDatum& header)
{

}

void VSEventDataVisitor::leaveSimHeader()
{

}

// ============================================================================
// VSEventDataDispatcher
// ============================================================================
VSEventDataDispatcher::VSEventDataDispatcher():
  m_file_list(), m_visitor()
{
  
}

VSEventDataDispatcher::VSEventDataDispatcher(VSEventDataVisitor* visitor):
  m_file_list(), m_visitor(visitor)
{
  
}

void VSEventDataDispatcher::loadList(const std::list<std::string>& file_list)
{
  for(std::list<std::string>::const_iterator itr = file_list.begin();
      itr != file_list.end(); ++itr)
    {
      std::string file = *itr;

      VSFileUtility::expandFilename(file);
      if(!VSFileUtility::isFile(file)) 
	{
	  std::cerr << std::string(__PRETTY_FUNCTION__) << ": "
		    << file << " is not a valid file." << std::endl;
	  exit(EXIT_FAILURE);
	}
      else if(!VSOctaveH5ReaderStruct::isHDF5(file))
	{
	  std::ifstream fileStream(itr->c_str());
	  VSLineTokenizer tokenizer;

	  while( !fileStream.eof() )
	    {
	      VSTokenList tokens;
	      std::string tmp;
	      getline(fileStream,tmp);
	      tokenizer.tokenize(tmp,tokens);
	      if(tokens.size() == 0)continue;

	      if(tokens[0].string().find('#') == std::string::npos)
		{
		  std::string h5_file = tokens[0].string();		  
		  VSFileUtility::expandFilename(h5_file);
		  if(!VSFileUtility::isFile(h5_file))
		    {
		      std::cerr 
			<< std::string(__PRETTY_FUNCTION__) << ": "
			<< h5_file 
			<< " is not a valid file." << std::endl;
		      exit(EXIT_FAILURE);
		    } 
		  else if(!VSOctaveH5ReaderStruct::isHDF5(h5_file))
		    {
		      std::cerr 
			<< std::string(__PRETTY_FUNCTION__) << ": "
			<< h5_file 
			<< " is not a valid HDF5 file." << std::endl;
		      exit(EXIT_FAILURE);
		    }
		  else
		    m_file_list.push_back(h5_file);	      
		}
	    }
	}
      else
	m_file_list.push_back(file);
    }
}

void VSEventDataDispatcher::
processFiles(const std::list<std::string>& file_list, unsigned nevents, 
	     bool verbose)
{
  loadList(file_list);
  processFiles(nevents,verbose);
}

void VSEventDataDispatcher::processFiles(unsigned nevents, bool verbose)
{


  if(!m_file_list.size())
    {
      std::cerr 
	<< std::string(__PRETTY_FUNCTION__) << ": "
	<< "Empty input file list." << std::endl;
      exit(EXIT_FAILURE);
    }

  // --------------------------------------------------------------------------
  // Loop over stage2 files
  // --------------------------------------------------------------------------
  for(std::vector<std::string>::iterator itr = m_file_list.begin();
      itr != m_file_list.end(); ++itr)
    {
      std::string stage2_file = *itr;
      VSFileUtility::expandFilename(stage2_file);
     
      try
	{
	  VSFileLockBSD::acquireLock(stage2_file);

	  VSOctaveH5Reader reader(stage2_file);

	  VSAnalysisStage1Data         stage1;
	  VSTargetTable::Observation   obs;
	  VSArrayMergedCalibrationData cal;
	  VSArrayDiagnosticsData       diagnostics;
	  VSHeaderSimulationDatum      sim_header;
	  VSArraySimulationReader*     sim_reader = NULL;

	  // ------------------------------------------------------------------
	  // Obtain the subset of parameters for each event that will
	  // be loaded
	  VSEventDataReader::MemberSubset sim_set;
	  VSEventDataReader::MemberSubset scope_set;
	  VSEventDataReader::MemberSubset array_set;

	  m_visitor->getSimSet(sim_set);
	  m_visitor->getScopeSet(scope_set);
	  m_visitor->getArraySet(array_set);

	  std::vector<unsigned> event_to_sim;    

	  stage1.load(reader.readStruct("stage1"));
	  obs.load(reader.readStruct("observation"));
	  cal.load(&reader);
	  //diagnostics.load(reader.readStruct("diagnostics"));

	  if(verbose)
	    {
	      std::cout 
		<< std::string(79,'-') << std::endl
		<< "File..........: " << stage2_file << std::endl
		<< "Run...........: " 
		<< stage1.run_info.run_number << std::endl
		<< "Date..........: " 
		<< stage1.run_info.first_event_time.getString() 
		<< std::endl
		<< "Source Name...: " << std::setw(15) 
		<< obs.name   << std::endl
		<< "Source RA.....: " << std::setw(15)
		<< obs.src_radec_J2000.longitude().hmsString(1) 
		<< std::endl
		<< "Source Dec....: " << std::setw(15)
		<< obs.src_radec_J2000.latitude().dmsString(1) 
		<< std::endl
		<< "Mode..........: " << std::setw(15)
		<< obs.mode_string << std::endl
		<< "Wobble Offset.: " << std::setw(15)
		<< obs.wobble_theta_rad*180./M_PI << std::endl
		<< "Wobble Angle..: " << std::setw(15)
		<< obs.wobble_phi_rad*180./M_PI   << std::endl
		<< "Mean El.......: " << std::setw(15)
		<< 90-stage1.run_info.zn_mean_deg << std::endl
		<< "Mean Az.......: " << std::setw(15)
		<< stage1.run_info.az_mean_deg  << std::endl
		<< "Scaled Dev....: " << std::setw(15)
		<< cal.mean_scaled_dev << std::endl;
	    }	  

	  // Determine whether this is a simulation file ----------------------
	  bool is_sim = false;
	  if(reader.isStruct("stage1.sim_info")) is_sim = true;

	  if(is_sim)
	    {
	      sim_header.load(reader.readStruct("sim_header"));
	      sim_reader = 
		new VSArraySimulationReader(reader.readStruct("sim_event"),
					    sim_set);
	      reader.readVector("event_to_sim",event_to_sim);
	      m_visitor->visitSimHeader(sim_header);
	    }

	  m_visitor->visitRun(stage1,obs,cal);


	  VSEventDataReader ed(reader.readStruct("events"),
			       array_set,scope_set);
	  unsigned nrow = ed.rows();
	  
	  std::cerr << stage2_file << ": " << nrow << std::endl;

	  // Loop on events ---------------------------------------------------
	  for(unsigned irow=0;irow<nrow;irow++)
	    {
	      if(nevents && irow == nevents) break;

	      VSEventArrayDatum                evt;
	      VSArraySimulationDatum           sim;
	      ed.element(evt,irow);
	      if(sim_reader) sim_reader->element(sim,event_to_sim[irow]);

	      if(is_sim) m_visitor->visitSimEvent(sim);
	      m_visitor->visitEvent(evt);

	      // Loop on telescopes -------------------------------------------
	      const unsigned nscope = evt.scope.size();
	      for(unsigned iscope = 0; iscope < nscope; iscope++)
		{
		  if(!evt.scope[iscope]) continue;
		  m_visitor->visitScope(*evt.scope[iscope],iscope);
		}

	      m_visitor->leaveEvent();
	      m_visitor->leaveSimEvent();
	    }

	  m_visitor->leaveRun();

	  if(is_sim) m_visitor->leaveSimHeader();	    

	  delete sim_reader;

	  VSFileLockBSD::releaseLock(stage2_file);
	}
      catch(const VSOctaveH5Exception& e)
	{
	  std::cerr << "Caught instance of VSOctaveH5Exception" << std::endl
		    << e.message() << std::endl
		    << "Skipping file: " << stage2_file << std::endl;
	}
      catch(const std::exception& x)
	{
	  std::cerr << "Caught instance of " << x.what() << std::endl
		    << "Skipping file: " << stage2_file << std::endl;
	}
      catch(...)
	{
	  std::cerr << "Caught some exception" << std::endl
		    << "Skipping file: " << stage2_file << std::endl;
	}
    }
}

