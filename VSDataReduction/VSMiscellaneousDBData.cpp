//-*-mode:c++; mode:font-lock;-*-

/*! \file VSMiscellaneousDBData.hpp

  Various data from the run that comes from the database but is not
  strictly necessary for the analaysis. (Happy birthday Eoin)

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       02/24/2007 

  $Id: VSMiscellaneousDBData.cpp,v 3.8 2008/04/14 22:29:59 sfegan Exp $

*/


#include <VSMiscellaneousDBData.hpp>
#include <VSChannelMap.hpp>

using namespace VERITAS;

bool VSMiscellaneousDBData::Scope::
getCorrectionParameters(const VSTime& time, 
			CorrectionParametersDatum& cpd) const
{
  for(std::vector<CorrectionParametersDatum>::const_iterator icp = 
	correction_parameters.begin(); icp != correction_parameters.end();
      icp++)
    if((time>=icp->db_start_time)
       &&((icp->db_end_time_is_null)||(time<=icp->db_end_time)))
      {
	cpd = *icp;
	return true;
      }
  return false;
}

bool VSMiscellaneousDBData::Scope::
getNextCorrectionParameters(const VSTime& time,
			    CorrectionParametersDatum& cpd) const
{
  std::vector<CorrectionParametersDatum>::const_iterator ifind = 
    correction_parameters.end();
  for(std::vector<CorrectionParametersDatum>::const_iterator icp = 
	correction_parameters.begin(); icp != correction_parameters.end();
      icp++)
    if((time>=icp->db_start_time)
       &&((icp->db_end_time_is_null)||(time<=icp->db_end_time)))
      {
	ifind = icp;
	break;
      }

  if(ifind != correction_parameters.end())ifind++;
  if(ifind == correction_parameters.end())return false;
  
  cpd = *ifind;
  return true;
}

void VSMiscellaneousDBData::clear()
{
  scope.clear();
  fir.clear();
  target_table.clear();
}

void VSMiscellaneousDBData::
fillFromDB(const VSTime& start_time, const VSTime& stop_time,
	   const std::vector<unsigned>& telescopes)
{
  VSChannelMap channel_map(start_time);
  VSCentralizedDBAccess* db = VSCentralizedDBAccess::getInstance();

  clear();

  // --------------------------------------------------------------------------
  // Scope information
  // --------------------------------------------------------------------------
	   
  unsigned nscope=0;
  for(unsigned iscopeindex=0;iscopeindex<telescopes.size();iscopeindex++)
    if(telescopes[iscopeindex]>=nscope)nscope=telescopes[iscopeindex]+1;

  if(nscope)
    {
      scope.resize(nscope);

      for(unsigned iscopeindex=0;iscopeindex<telescopes.size();iscopeindex++)
	{
	  unsigned iscope = telescopes[iscopeindex];

	  scope[iscope].has_scope = true;

	  // Tracking Targets -------------------------------------------------

	  db->getTrackingTarget(iscope, start_time, stop_time, 
				scope[iscope].tracking_targets);

	  // Correction Parameters --------------------------------------------

	  db->getCorrectionParameters(iscope, 
				      scope[iscope].correction_parameters);

	  // HV Status --------------------------------------------------------

	  std::vector<VSCentralizedDBAccess::HVStatusDatum>  hv_status_data;
	  db->getHVStatus(iscope, start_time, stop_time, hv_status_data);

	  if(hv_status_data.size())
	    {
	      unsigned ndata = hv_status_data.size();

	      unsigned nchan = 0;
	      for(unsigned idatum = 0; idatum<ndata; idatum++)
		if(hv_status_data[idatum].ichan >= nchan)
		  nchan = hv_status_data[idatum].ichan+1;

	      unsigned nmeas = 0;
	      unsigned nmeaschan = 0;
	      VATime last_meas_time = hv_status_data.front().timestamp;
	      for(unsigned idatum = 0; idatum<ndata; idatum++)
		{
		  if(hv_status_data[idatum].timestamp != last_meas_time)
		    {
		      if(nmeaschan == nchan)nmeas++;
		      last_meas_time = hv_status_data[idatum].timestamp;
		      nmeaschan = 0;
		    }
		  nmeaschan++;
		}
	      if(nmeaschan == nchan)nmeas++;

	      scope[iscope].hv_status.resize(nmeas,nchan);
	  
	      std::vector<float> temp_current(nchan);
	      std::vector<float> temp_voltage(nchan);
	      unsigned imeas = 0;
	      nmeaschan = 0;
	      last_meas_time = hv_status_data.front().timestamp;
	      for(unsigned idatum = 0; idatum<ndata; idatum++)
		{
		  if(hv_status_data[idatum].timestamp != last_meas_time)
		    {
		      if(nmeaschan == nchan)
			{
			  for(unsigned jchan=0;jchan<nchan;jchan++)
			    {
			      scope[iscope].hv_status[imeas].chan[jchan].
				voltage = temp_voltage[jchan];
			      scope[iscope].hv_status[imeas].chan[jchan].
				current = temp_current[jchan];
			    }
			  scope[iscope].hv_status[imeas].timestamp = 
			    last_meas_time;
			  imeas++;
			}
		      last_meas_time = hv_status_data[idatum].timestamp;
		      nmeaschan = 0;
		    }

		  unsigned ichan = hv_status_data[idatum].ichan;
		  nmeaschan++;
		  temp_voltage[ichan] = hv_status_data[idatum].voltage;
		  temp_current[ichan] = hv_status_data[idatum].current;
		}
	  
	      if(nmeaschan == nchan)
		{
		  for(unsigned jchan=0;jchan<nchan;jchan++)
		    {
		      scope[iscope].hv_status[imeas].chan[jchan].voltage = 
			temp_voltage[jchan];
		      scope[iscope].hv_status[imeas].chan[jchan].current = 
			temp_current[jchan];
		    }
		  scope[iscope].hv_status[imeas].timestamp = last_meas_time;
		}
	    }

	  // L1 Rate ----------------------------------------------------------

	  std::vector<VSCentralizedDBAccess::L1RateDatum>  l1_rate_data;
	  db->getL1Rate(iscope, start_time, stop_time, l1_rate_data);

	  if(l1_rate_data.size())
	    {
	      unsigned ndata = l1_rate_data.size();
	      unsigned nchan = 0;
	      unsigned nmeas = 1;
	      VATime last_meas_time = l1_rate_data.front().timestamp;
	      for(unsigned idatum = 0; idatum<ndata; idatum++)
		{
		  if(l1_rate_data[idatum].ichan >= nchan)
		    nchan = l1_rate_data[idatum].ichan+1;
		  if(l1_rate_data[idatum].timestamp != last_meas_time)
		    nmeas++, last_meas_time = l1_rate_data[idatum].timestamp;
		}

	      scope[iscope].l1_rate.resize(nmeas,nchan);
	  
	      unsigned imeas = 0;
	      last_meas_time = l1_rate_data.front().timestamp;
	      scope[iscope].l1_rate[imeas].timestamp = last_meas_time;
	      for(unsigned idatum = 0; idatum<ndata; idatum++)
		{
		  unsigned ichan = l1_rate_data[idatum].ichan;
		  scope[iscope].l1_rate[imeas].rate[ichan] =
		    l1_rate_data[idatum].rate;
		  if(l1_rate_data[idatum].timestamp != last_meas_time)
		    imeas++,
		      last_meas_time = l1_rate_data[idatum].timestamp,
		      scope[iscope].l1_rate[imeas].timestamp = last_meas_time;
		}
	    }

	  // L3 Telescope -----------------------------------------------------

	  db->getL3Scope(iscope, start_time, stop_time, 
			 scope[iscope].l3_scope);
	}
    }

  // --------------------------------------------------------------------------
  // FIR
  // --------------------------------------------------------------------------

  const unsigned nfir = channel_map.nfir();
  if(nfir)
    {
      fir.resize(nfir);
      for(unsigned ifir=0;ifir<nfir;ifir++)
	db->getFIR(ifir, start_time, stop_time, fir[ifir]);
    }

  // --------------------------------------------------------------------------
  // TARGET TABLE
  // --------------------------------------------------------------------------
  
  db->getTargetTable(target_table);
}

void VSMiscellaneousDBData::
load(VSOctaveH5ReaderStruct* reader)
{
  clear();

  // --------------------------------------------------------------------------
  // Read scope array
  // --------------------------------------------------------------------------

  if(reader->isValid("t"))
    {
      VSOctaveH5ReaderCellVector* c = reader->readCellVector("t");
      vsassert(c);
      unsigned nscope = c->dimensions();

      scope.resize(nscope);

      for(unsigned iscope=0; iscope<nscope; iscope++)
	if(!c->isEmpty(iscope))
	  {
	    scope[iscope].has_scope = true;

	    VSOctaveH5ReaderStruct* s = c->readStruct(iscope);
	    vsassert(s);

	    if(s->isValid("tracking_targets"))
	      {
		VSOctaveH5ReaderCellVector* cc = 
		  s->readCellVector("tracking_targets");
		vsassert(cc);
		unsigned ntar = cc->dimensions();
		scope[iscope].tracking_targets.resize(ntar);
		for(unsigned itar=0;itar<ntar;itar++)
		  cc->readComposite(itar,scope[iscope].tracking_targets[itar]);
		delete cc;
	      }
	    
	    if(s->isValid("correction_parameters"))
	      s->readCompositeVector("correction_parameters",
				     scope[iscope].correction_parameters);
	    
	    if(s->isValid("hv_status"))
	      {
		VSOctaveH5ReaderStruct* ss = s->readStruct("hv_status");
		vsassert(ss);
		unsigned nchan;
		unsigned nmeas;
		std::vector<uint32_t> t_mjd;
		std::vector<uint64_t> t_dns;
		float* current;
		float* voltage;
		ss->readVector("timestamp_mjd",t_mjd);
		ss->readVector("timestamp_dns",t_dns);
		ss->readMatrix("current",nmeas,nchan,current);
		ss->readMatrix("voltage",nmeas,nchan,voltage);
		scope[iscope].hv_status.resize(nmeas,nchan);	    
		for(unsigned imeas=0;imeas<nmeas;imeas++)
		  {
		    OneHVMeasScope& M(scope[iscope].hv_status[imeas]);
		    M.timestamp.setFromMJDIntAndNS(t_mjd[imeas],t_dns[imeas]);
		    for(unsigned ichan=0;ichan<nchan;ichan++)
		      {
			M.chan[ichan].current = current[imeas*nchan+ichan];
			M.chan[ichan].voltage = voltage[imeas*nchan+ichan];
		      }
		  }
		delete[] voltage;
		delete[] current;
		delete ss;	    
	      }

	    if(s->isValid("l1_rate"))
	      {
		VSOctaveH5ReaderStruct* ss = s->readStruct("l1_rate");
		vsassert(ss);
		unsigned nchan;
		unsigned nmeas;
		std::vector<uint32_t> t_mjd;
		std::vector<uint64_t> t_dns;
		float* rate;
		ss->readVector("timestamp_mjd",t_mjd);
		ss->readVector("timestamp_dns",t_dns);
		ss->readMatrix("rate",nmeas,nchan,rate);
		scope[iscope].l1_rate.resize(nmeas,nchan);	    
		for(unsigned imeas=0;imeas<nmeas;imeas++)
		  {
		    OneL1RateMeasScope& M(scope[iscope].l1_rate[imeas]);
		    M.timestamp.setFromMJDIntAndNS(t_mjd[imeas],t_dns[imeas]);
		    for(unsigned ichan=0;ichan<nchan;ichan++)
		      M.rate[ichan] = rate[imeas*nchan+ichan];
		  }
		delete[] rate;
		delete ss;	    
	      }

	    if(s->isValid("l3_scope"))
	      s->readCompositeVector("l3_scope", scope[iscope].l3_scope);

	    if(s->isValid("vpm_stars"))
	      s->readCompositeVector("vpm_stars", scope[iscope].vpm_stars);

	    if(s->isValid("vpm_leds"))
	      s->readCompositeVector("vpm_leds", scope[iscope].vpm_leds);

	    delete s;
	  }  
      delete c;
    }

  // --------------------------------------------------------------------------
  // FIR measurements
  // --------------------------------------------------------------------------

  if(reader->isValid("fir"))
    {
      VSOctaveH5ReaderCellVector* c = reader->readCellVector("fir");
      vsassert(c);
      unsigned nfir = c->dimensions();
      fir.resize(nfir);
      for(unsigned ifir=0;ifir<nfir;ifir++)
	c->readCompositeVector(ifir,fir[ifir]);
    }      

  // --------------------------------------------------------------------------
  // Target table
  // --------------------------------------------------------------------------

  if(reader->isValid("target_table"))
    reader->readCompositeVector("target_table",target_table);
}

void VSMiscellaneousDBData::
save(VSOctaveH5WriterStruct* writer) const
{
  // --------------------------------------------------------------------------
  // Write scope array
  // --------------------------------------------------------------------------

  const unsigned nscope = scope.size();
  if(nscope)
    {
      VSOctaveH5WriterCellVector* c = writer->writeCellVector("t",nscope);

      for(unsigned iscope=0; iscope<nscope; iscope++)
	if(scope[iscope].has_scope)
	  {
	    VSOctaveH5WriterStruct* s = c->writeStruct(iscope);

	    if(!scope[iscope].tracking_targets.empty())
	      {
		unsigned ntar = scope[iscope].tracking_targets.size();
		VSOctaveH5WriterCellVector* cc = 
		  s->writeCellVector("tracking_targets",ntar);
		for(unsigned itar=0;itar<ntar;itar++)
		  cc->writeComposite(itar,scope[iscope].
				     tracking_targets[itar]);
		delete cc;
	      }

	    if(!scope[iscope].correction_parameters.empty())
	      s->writeCompositeVector("correction_parameters",
				      scope[iscope].correction_parameters);

	    if(!scope[iscope].hv_status.empty())
	      {
		unsigned nmeas = scope[iscope].hv_status.size();
		unsigned nchan = scope[iscope].hv_status[0].chan.size();
		std::vector<uint32_t> t_mjd(nmeas);
		std::vector<uint64_t> t_dns(nmeas);
		float* current = new float[nmeas*nchan];
		float* voltage = new float[nmeas*nchan];
		for(unsigned imeas=0;imeas<nmeas;imeas++)
		  {
		    const OneHVMeasScope& M(scope[iscope].hv_status[imeas]);
		    t_mjd[imeas] = M.timestamp.getMJDInt();
		    t_dns[imeas] = M.timestamp.getDayNS();
		    for(unsigned ichan=0;ichan<nchan;ichan++)
		      {
			current[imeas*nchan+ichan] = M.chan[ichan].current;
			voltage[imeas*nchan+ichan] = M.chan[ichan].voltage;
		      }
		  }
		VSOctaveH5WriterStruct* ss = s->writeStruct("hv_status");
		ss->writeVector("timestamp_mjd",t_mjd);
		ss->writeVector("timestamp_dns",t_dns);
		ss->writeMatrix("current",nmeas,nchan,current);
		ss->writeMatrix("voltage",nmeas,nchan,voltage);
		delete ss;	    
		delete[] voltage;
		delete[] current;
	      }

	    if(!scope[iscope].l1_rate.empty())
	      {
		unsigned nmeas = scope[iscope].l1_rate.size();
		unsigned nchan = scope[iscope].l1_rate[0].rate.size();
		std::vector<uint32_t> t_mjd(nmeas);
		std::vector<uint64_t> t_dns(nmeas);
		float* rate = new float[nmeas*nchan];
		for(unsigned imeas=0;imeas<nmeas;imeas++)
		  {
		    const OneL1RateMeasScope& M(scope[iscope].l1_rate[imeas]);
		    t_mjd[imeas] = M.timestamp.getMJDInt();
		    t_dns[imeas] = M.timestamp.getDayNS();
		    for(unsigned ichan=0;ichan<nchan;ichan++)
		      rate[imeas*nchan+ichan] = M.rate[ichan];
		  }
		VSOctaveH5WriterStruct* ss = s->writeStruct("l1_rate");
		ss->writeVector("timestamp_mjd",t_mjd);
		ss->writeVector("timestamp_dns",t_dns);
		ss->writeMatrix("rate",nmeas,nchan,rate);
		delete ss;	    
		delete[] rate;
	      }

	    if(!scope[iscope].l3_scope.empty())
	      s->writeCompositeVector("l3_scope", scope[iscope].l3_scope);

	    if(!scope[iscope].vpm_stars.empty())
	      s->writeCompositeVector("vpm_stars", scope[iscope].vpm_stars);

	    if(!scope[iscope].vpm_leds.empty())
	      s->writeCompositeVector("vpm_leds", scope[iscope].vpm_leds);

	    delete s;
	  }

      delete c;
    }

  // --------------------------------------------------------------------------
  // FIR measurements
  // --------------------------------------------------------------------------

  unsigned nfir = fir.size();
  if(nfir)
    {
      VSOctaveH5WriterCellVector* c = writer->writeCellVector("fir",nfir);
      vsassert(c);
      for(unsigned ifir=0;ifir<nfir;ifir++)
	c->writeCompositeVector(ifir,fir[ifir]);
    }

  // --------------------------------------------------------------------------
  // Target table
  // --------------------------------------------------------------------------

  if(!target_table.empty())
    writer->writeCompositeVector("target_table",target_table);
}

