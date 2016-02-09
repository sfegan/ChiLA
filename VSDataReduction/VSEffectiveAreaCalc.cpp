#include <VSEffectiveAreaCalc.hpp>
#include <VSANonlinearFitting.hpp>
#include <VSALinearLeastSquares.hpp>
#include <VSNSpaceOctaveH5IO.hpp>
#include <VSALocalRegression.hpp>

using namespace VERITAS;


// ============================================================================
// VSEffectiveAreaCalc
// ============================================================================
VSPolyFn::VSPolyFn(unsigned n1, unsigned n2, double x0): 
  VSAFunction::ParamFn<double>(n1+n2),
  m_x0(x0), m_n(n1+n2), m_p1(n1), m_p2(n2)
{
  
}

void VSPolyFn::operator() (const double& x, VSAAlgebra::VecND& v) const
{
  v.resize(m_n);
  v.clear();

  if(x < m_x0)
    {
      VSAAlgebra::VecND v1;
      m_p1(x,v1);

      for(unsigned ip = 0; ip < v1.ndim(); ip++)
	v[ip] = v1[ip];
    }
  else
    {
      VSAAlgebra::VecND f1;
      VSAAlgebra::VecND f2;
      VSAAlgebra::VecND f3;

      VSAAlgebra::VecND df1;
      VSAAlgebra::VecND df2;
      

      m_p1(m_x0,f1);
      m_p2(m_x0,f2);
      m_p1.dydx(m_x0,df1);
      m_p2.dydx(m_x0,df2);

      m_p2(x,f3);

      double dx = x-m_x0;

      const unsigned nd1 = f1.ndim();
      for(unsigned ip = 0; ip < nd1; ip++) v[ip] = f1[ip] + dx*df1[ip];

      const unsigned nd2 = f2.ndim();
      for(unsigned ip = 2; ip < nd2; ip++)
	v[nd1+ip-2] = f3[ip] - f2[ip] - dx*df2[ip];
    }
  
}


double VSPolyFn::val(const double& x) const
{
  VSAAlgebra::VecND v;
  (*this)(x,v);
  return VSAFunction::ParamFn<double>::param()*v;
}

double VSPolyFn::val(const double& x, 
		     const VSAAlgebra::VecND& a) const
{
  VSAAlgebra::VecND v;
  (*this)(x,v);
  return a * v;
}

void VSPolyFn::dyda(const double& x, 
		    VSAAlgebra::VecND& dyda) const
{
  (*this)(x,dyda);
}

void VSPolyFn::dyda(const double& x, 
		    const VSAAlgebra::VecND& a,
		    VSAAlgebra::VecND& dyda) const
{
  (*this)(x,dyda);
}

void VSEffectiveAreaCalc::Data::save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeScalar("offset_deg",offset_deg);

  egymc_total_hist.save(writer->writeStruct("egymc_total_hist"));
  egymc_selected_hist.save(writer->writeStruct("egymc_selected_hist"));
  egymc_area_hist.save(writer->writeStruct("egymc_area_hist"));
  egymc_effarea_hist.save(writer->writeStruct("egymc_effarea_hist"));
  egymc_log_effarea_hist.save(writer->writeStruct("egymc_log_effarea_hist"));

  egymc_effarea_fit_hist.save(writer->writeStruct("egymc_effarea_fit_hist"));
 //egymc_effarea_fit2_hist.save(writer->writeStruct("egymc_effarea_fit2_hist"));
  egymc_log_effarea_fit_hist.
    save(writer->writeStruct("egymc_log_effarea_fit_hist"));  



  const unsigned ndim = effarea_fit_cov.ndim();
  double * cov = new double[ndim*ndim];
  std::vector<double> p;
  std::vector<double> err;
  
  for(unsigned idim1 = 0; idim1 < ndim; idim1++)
    {
      p.push_back(effarea_fit_param(idim1));
      err.push_back(sqrt(effarea_fit_cov(idim1,idim1)));
      for(unsigned idim2 = 0; idim2 < ndim; idim2++)
	cov[idim2+idim1*ndim] = effarea_fit_cov(idim1,idim2);
    }


  writer->writeVector("effarea_param",p);
  writer->writeVector("effarea_param_err",err);
  writer->writeMatrix("effarea_param_cov",ndim,ndim,cov);

  delete [] cov;
}

VSEffectiveAreaCalc::VSEffectiveAreaCalc():
  m_effarea_param_hist(),
  m_log10_egybin(),
  m_log10_egylo(),
  m_log10_egyhi(),
  m_negy(),
  m_data(), m_data_ptr()
{

}

VSEffectiveAreaCalc::~VSEffectiveAreaCalc()
{
  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++) delete m_data[idata];
}

void VSEffectiveAreaCalc::loadSimInfo(VSSimInfoData* sim_info)
{
  m_data_ptr = NULL;

  for(std::vector<Data*>::iterator itr = m_data.begin();
      itr != m_data.end(); ++itr)
    {
      double doff = fabs(sim_info->wobble_theta_deg-(*itr)->offset_deg);
      
      if(doff < 0.05) 
	{
	  m_data_ptr = *itr;
	  break;
	}
      else if(sim_info->wobble_theta_deg < (*itr)->offset_deg)
	{
	  m_data_ptr = new Data(sim_info->wobble_theta_deg,
				m_log10_egybin,m_log10_egylo,m_log10_egyhi);
	  m_data.insert(itr,m_data_ptr);
	  break;
	}
    }

  if(m_data_ptr == NULL)
    {
      m_data_ptr = new Data(sim_info->wobble_theta_deg,
			    m_log10_egybin,m_log10_egylo,m_log10_egyhi);
      m_data.push_back(m_data_ptr);
    }
}

void VSEffectiveAreaCalc::loadHeader(const VSHeaderSimulationDatum& sim_header)
{
  vsassert(m_data_ptr != NULL);

  for(std::vector< VSTableSimulationDatum >::const_iterator itr = 
	sim_header.tables.begin(); itr != sim_header.tables.end(); ++itr)
    {
      double area = std::pow(itr->sampling_radius_m,2)*M_PI;
      double log10_egy = std::log10(itr->energy_tev);
      
      if(m_data_ptr->egymc_area_hist.countForVal(log10_egy) == 0)
	m_data_ptr->egymc_area_hist.accumulate(log10_egy,area,0);

      m_data_ptr->egymc_total_hist.accumulate(log10_egy,itr->event_count,0);
    }
}

void VSEffectiveAreaCalc::setEnergyBinning(double egybin, 
					   double egylo, double egyhi)
{
  if(m_log10_egybin == 0)
    {
      m_log10_egybin = egybin;
      m_log10_egylo = egylo;
      m_log10_egyhi = egyhi;
      m_negy = lround((egyhi-egylo)/egybin);

      m_effarea_offset_hist = 
	VSSimple2DHist<double,double>(egybin,egylo,egyhi,0.04,0,1.8);
    }
}

void VSEffectiveAreaCalc::accumulate(double energy_tev)
{
  m_data_ptr->egymc_selected_hist.accumulate(std::log10(energy_tev));
}

void VSEffectiveAreaCalc::calcEffarea()
{
  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    calcEffarea(m_data[idata]);
}

void VSEffectiveAreaCalc::calcEffarea(Data* data)
{
  for(VSLimitedErrorsHist<double, double>::iterator itr = 
	data->egymc_total_hist.begin(); itr != 
	data->egymc_total_hist.end(); ++itr)
    {
      double log10_egy = itr->center();  

      double ntotal = itr->count();
      double nselected = 
	data->egymc_selected_hist.countForVal(itr->center());
      double area = data->egymc_area_hist.countForVal(itr->center());

      double eff = nselected/ntotal;
      double eff_var = eff*(1-eff)/ntotal;
      double effarea = eff*area;
      double effarea_var = eff_var*area*area;

      data->egymc_effarea_hist.accumulate(log10_egy,effarea,effarea_var);

      if(effarea > 0)
	{
	  data->egymc_log_effarea_hist.
	    accumulate(log10_egy,std::log(effarea),
		       effarea_var/std::pow(effarea,2));
	}
    }
}

void VSEffectiveAreaCalc::fit()
{
  std::vector< std::vector< std::pair<double,double> > > offset_spline_data;

  const unsigned negy = m_effarea_offset_hist.nXBins();
  offset_spline_data.resize(negy);

  std::vector< VSAFunction::Spline > offset_spline;
  VSAMath::Data<VSAAlgebra::Vec2D> offset_fit_data;

  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    {
      std::cout << std::string(__PRETTY_FUNCTION__)  
		<< ": Fitting Effective Area Model " 
		<< " Offset  " << std::setw(15) << m_data[idata]->offset_deg
		<< std::endl;

      fit(m_data[idata]);

      m_effarea_param_hist.
	resize(m_data[idata]->log_effarea_fit_param.ndim(),
	       VSLimitedErrorsHist<double,double>(0.05,0.,3));

      offset_spline.resize(m_data[idata]->effarea_fit_param.ndim());

      for(unsigned ip = 0; ip < m_data[idata]->effarea_fit_param.ndim(); ip++)
	{
	  offset_spline[ip].setPoint(m_data[idata]->offset_deg,
				     m_data[idata]->effarea_fit_param[ip]);
	}		     

      for(unsigned ip = 0; ip < m_data[idata]->log_effarea_fit_param.ndim(); 
	  ip++)
	{
	  m_effarea_param_hist[ip].
	    accumulate(m_data[idata]->offset_deg,
		       m_data[idata]->log_effarea_fit_param[ip],
		       m_data[idata]->log_effarea_fit_cov(ip,ip));
	}

      for(unsigned iegy = 0; iegy < negy; iegy++)
	{	  
	  double log10_egy = m_effarea_offset_hist.xBinToCenter(iegy);
	  double effarea = 
	    m_data[idata]->egymc_log_effarea_fit_hist.countForVal(log10_egy);
	  offset_spline_data[iegy].
	    push_back(std::make_pair(m_data[idata]->offset_deg,effarea));

	  m_effarea_offset_hist.accumulate(log10_egy,
					   m_data[idata]->offset_deg,
					   effarea);
	}

      for(unsigned iegy = 0; iegy < negy; iegy++)
	{	  
	  double log10_egy = m_effarea_offset_hist.xBinToCenter(iegy);
	  double effarea = 
	    m_data[idata]->egymc_log_effarea_fit_hist.countForVal(log10_egy);
	  double effarea_var = 
	    m_data[idata]->egymc_log_effarea_fit_hist.varForVal(log10_egy);

	  if(effarea_var == 0) effarea_var = 1;

	  if(effarea == 0) continue;

	  VSAAlgebra::Vec2D cnd(log10_egy,m_data[idata]->offset_deg);
	  
	  VSAMath::DataPoint<VSAAlgebra::Vec2D> dp(cnd,effarea,
						   sqrt(effarea_var));

	  offset_fit_data.insert(dp);
	}
    }

  // VSAMath::LocalRegression2D<VSAAlgebra::Vec2D> 
  //   lr(offset_fit_data,0.25,0.5,1,5);

  for(unsigned iegy = 0; iegy < negy; iegy++)
    {
      double log10_egy = m_effarea_offset_hist.xBinToCenter(iegy);
      vsassert(iegy < offset_spline_data.size());

      VSAFunction::Spline spline(offset_spline_data[iegy]);
      for(unsigned iy = 0; iy < m_effarea_offset_hist.nYBins(); iy++)
	{
	  double offset = m_effarea_offset_hist.yBinToCenter(iy);
	  //double log_effarea = lr.val(VSAAlgebra::Vec2D(log10_egy,offset));
	  m_effarea_offset_hist.setBin(iegy,iy,std::exp(spline.val(offset)));
	}
    }

  // for(unsigned idata = 0; idata < ndata; idata++)
  //   {
  //     double offset = m_data[idata]->offset_deg;

  //     for(unsigned iegy = 0; iegy < negy; iegy++)
  // 	{
  // 	  double log10_egy = m_effarea_offset_hist.xBinToCenter(iegy);
  // 	  double log_effarea = lr.val(VSAAlgebra::Vec2D(log10_egy,offset));
  // 	  m_data[idata]->egymc_effarea_fit2_hist.
  // 	    accumulate(log10_egy,std::exp(log_effarea),0.0);
  // 	}
  //   }

  m_effarea_offset_nspace = VSNSpace(m_effarea_offset_hist);  
}

void VSEffectiveAreaCalc::fit(Data* data)
{
  VSAMath::Data<double> fit_data;
  VSAMath::Data<double> log_fit_data;

  VSLimitedErrorsHist<double,double>& hlog = data->egymc_log_effarea_hist;

  double emin = hlog.hiLimit();
  double emax = hlog.loLimit();

  for(VSLimitedErrorsHist<double,double>::iterator itr = hlog.begin();
      itr != hlog.end(); ++itr)
    {
      double ncount = data->egymc_selected_hist.countForVal(itr->center());

      if(itr->center() < emin && ncount > 5) emin = itr->center();
      if(itr->center() > emax && ncount > 1) emax = itr->center();
      if(ncount == 0) continue;

      VSAMath::DataPoint<double> p(itr->center(),itr->count(),itr->err());
      log_fit_data.insert(p);
    }

  // --------------------------------------------------------------------------
  // Perform a preliminary fit to the effective area which will be
  // used to find the break energy.
  // --------------------------------------------------------------------------
  VSAMath::LocalRegression1D<double> lr(log_fit_data,0.35,1,5);

  // --------------------------------------------------------------------------
  // Find the break energy identified as the point at which the
  // the logarithmic slope > 2
  // --------------------------------------------------------------------------
  double x0 = 0;
  const double dx = 0.0001;
  for(double x = emin; x < 1.5; x += 0.01)
    {
      if(x > emax) continue;

      double y1 = lr.val(x-dx);
      double y2 = lr.val(x);

      double dydx = (y2-y1)/dx;

      if(dydx < 2)
	{
	  x0 = x;
	  break;
	}
    }

  VSPolyFn polyfn(5,4,x0);

  VSAMath::Fitsvd<VSPolyFn,double> svdfit(log_fit_data,polyfn);
  svdfit.fit();

  data->log_effarea_fit_param = svdfit.param();
  data->log_effarea_fit_cov = svdfit.cov();

  double y1 = polyfn.val(emin,data->log_effarea_fit_param);
  double y2 = polyfn.val(emin+0.0001,data->log_effarea_fit_param);
  double dydx = (y2-y1)/0.0001;

  for(VSLimitedErrorsHist<double,double>::iterator itr = 
	data->egymc_effarea_hist.begin();
      itr != data->egymc_effarea_hist.end(); ++itr)
    {     
      double v = 
	polyfn.val(itr->center(),data->log_effarea_fit_param);
      //      data->egymc_log_effarea_fit_hist.accumulate(itr->center(),v,0);

      VSAAlgebra::VecND dyda;
      polyfn.dyda(itr->center(),dyda);

      double v_var = dyda*svdfit.cov()*dyda;

      if(itr->center() > emin)
	{
	  data->egymc_log_effarea_fit_hist.
	    accumulate(itr->center(),v,sqrt(v_var));
	  data->egymc_effarea_fit_hist.accumulate(itr->center(),std::exp(v),0);
	}
      else if(itr->center() <= emin)
	{
	  v = dydx*itr->center() + y1 - dydx*emin;
	  data->egymc_log_effarea_fit_hist.accumulate(itr->center(),v,0);
	  data->egymc_effarea_fit_hist.accumulate(itr->center(),std::exp(v),0);
	}
    }

}

void VSEffectiveAreaCalc::save(VSOctaveH5WriterStruct* writer) const
{
  const unsigned ndata = m_data.size();

  VSOctaveH5WriterCellVector* wc = writer->writeCellVector("offset", ndata);
  vsassert(wc);

  for(unsigned idata = 0; idata < ndata; idata++)
    {      
      VSOctaveH5WriterStruct* ws = wc->writeStruct(idata);
      vsassert(ws); 
      m_data[idata]->save(ws);
      delete ws;
    }
  delete wc;

  m_effarea_offset_hist.save(writer->writeStruct("effarea_offset_hist"));

  const unsigned np = m_effarea_param_hist.size();
  
  wc = writer->writeCellVector("effarea_param_hist", np);
  for(unsigned ip = 0; ip < np; ip++)
    {
      VSOctaveH5WriterStruct* ws = wc->writeStruct(ip);
      m_effarea_param_hist[ip].save(ws);
      delete ws;
    }
  delete wc;
}

// ============================================================================
// VSEnergyKernelFn
// ============================================================================
VSEnergyKernelFn::VSEnergyKernelFn(double dx, unsigned negy,
				   double dloge, double logemin): 
  VSAFunction::ParamFn<VSACoord::CoordND>(), 
  m_pol(),
  m_ebin(dloge,logemin)
{
  m_pol.push_back(new VSAFunction::Poly(2)); 
  m_pol.push_back(new VSAFunction::Poly(2)); 
  m_pol.push_back(new VSAFunction::Poly(1)); 
  m_pol.push_back(new VSAFunction::Poly(2)); 
  m_pol.push_back(new VSAFunction::Poly(0));



  unsigned nparm = 0;
  for(unsigned ipol = 0; ipol < m_pol.size(); ipol++)
    nparm += m_pol[ipol]->nparm();

  setNParam(nparm);

  const unsigned np1 = m_pol[0]->nparm();
  const unsigned np2 = m_pol[1]->nparm();
  const unsigned np3 = m_pol[2]->nparm();
  const unsigned np4 = m_pol[3]->nparm();

  
  setParam(0,0.06);
  setParam(np1+1,1);
  setParam(np1+np2,0.05);
  setParam(np1+np2+np3+1,1);
  setParam(np1+np2+np3+np4,0.8);

  for(unsigned iegy = 0; iegy < negy; iegy++)
    m_fn.push_back(VSEnergyResponseFn(dx));
}

VSEnergyKernelFn::~VSEnergyKernelFn()
{
  for(unsigned ipol = 0; ipol < m_pol.size(); ipol++)
    delete m_pol[ipol];
}

void VSEnergyKernelFn::setNorm(unsigned iegy, double counts)
{
  vsassert(iegy < m_fn.size());
  m_fn[iegy].setNorm(counts);
}

double VSEnergyKernelFn::val(const VSACoord::CoordND& x) const 
{ 
  return val(x,param());
}

double VSEnergyKernelFn::val(const VSACoord::CoordND& x, 
			     const VSAAlgebra::VecND& a) const 
{ 
  unsigned iegy = m_ebin.valToBin(x[1]);
  VSAAlgebra::VecND afn(5);  
  coeff(x[1],a,afn);

  double v = m_fn[iegy].val(x,afn);

  return v;
  //return m_fn[iegy].val(x,afn);
}

double VSEnergyKernelFn::val(const VSACoord::CoordND& x, 
			     const VSAAlgebra::VecND& a,
			     double counts) const 
{ 
  unsigned iegy = m_ebin.valToBin(x[1]);  
  VSAAlgebra::VecND afn(5);
  coeff(x[1],a,afn);

  VSACoord::CoordND x2(1);
  x2[0] = x[0];
  return m_fn[iegy].val(x2,afn,counts);
}

void VSEnergyKernelFn::dyda(const VSACoord::CoordND& x, 
			    VSAAlgebra::VecND& dyda) const 
{
  VSEnergyKernelFn::dyda(x,param(),dyda);
}

void VSEnergyKernelFn::dyda(const VSACoord::CoordND& x, 
			    const VSAAlgebra::VecND& a, 
			    VSAAlgebra::VecND& dyda) const 
{
  dyda.resize(nparm());

  VSAAlgebra::VecND afn(5);
  coeff(x[1],a,afn);

  unsigned iegy = m_ebin.valToBin(x[1]);
  VSAAlgebra::VecND dfdn;
  m_fn[iegy].dyda(x,afn,dfdn);

  unsigned ip1 = 0;
  for(unsigned ipol = 0; ipol < m_pol.size(); ipol++)
    {
      VSAAlgebra::VecND fpol;
      m_pol[ipol]->dyda(x[1],fpol);

      for(unsigned ip2 = 0; ip2 < fpol.ndim(); ip2++)
	{
	  dyda(ip1) = dfdn(ipol)*fpol(ip2);
	  ip1++;
	}
    }
}

VSEnergyKernelFn::VSEnergyKernelFn(const VSEnergyKernelFn& o):
  VSAFunction::ParamFn<VSACoord::CoordND>(o)
{
  for(unsigned ipol = 0; ipol < o.m_pol.size(); ipol++)
    m_pol.push_back(o.m_pol[ipol]->clone());

  m_fn = o.m_fn;
  m_ebin = o.m_ebin;
}

// VSEnergyKernelFn& VSEnergyKernelFn::operator= (const VSEnergyKernelFn& o)
// {
  
// }

void VSEnergyKernelFn::coeff(double x, 
			     const VSAAlgebra::VecND& a,
			     VSAAlgebra::VecND& afn) const
{
  afn.resize(5);
  unsigned ip = 0;
  for(unsigned ipol = 0; ipol < m_pol.size(); ipol++)
    {
      afn(ipol) = m_pol[ipol]->val(x,a.subVector(ip,m_pol[ipol]->nparm()));
      ip += m_pol[ipol]->nparm();
    }
}

void VSEnergyKernelFn::coeff_var(double x, 
				 const VSAAlgebra::VecND& a,
				 const VSAAlgebra::MatrixND& cov,
				 VSAAlgebra::VecND& var) const
{
  var.resize(5);

  unsigned ip = 0;
  for(unsigned ipol = 0; ipol < m_pol.size(); ipol++)
    {
      const unsigned np = m_pol[ipol]->nparm();
      VSAAlgebra::VecND ap = a.subVector(ip,np);
      VSAAlgebra::VecND dyda;
      m_pol[ipol]->dyda(x,ap,dyda);
      VSAAlgebra::MatrixND c = cov.subMatrix(ip,ip,np,np);
      var(ipol) = dyda*(c*dyda);
      ip += m_pol[ipol]->nparm();
    }
}

// ============================================================================
// VSEnergyKernelCalc
// ============================================================================
void VSEnergyKernelCalc::Data::save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeScalar("offset_deg",offset_deg);

  egymc_egy_hist.save(writer->writeStruct("egymc_egy_hist"));
  egymc_egyerr_hist.save(writer->writeStruct("egymc_egyerr_hist"));
  egymc_biasdev_hist.save(writer->writeStruct("egymc_biasdev_hist"));
  egymc_mse_hist.save(writer->writeStruct("egymc_mse_hist"));
  egymc_bias_hist.save(writer->writeStruct("egymc_bias_hist"));
  egymc_dev_hist.save(writer->writeStruct("egymc_dev_hist"));
  egymc_res_hist.save(writer->writeStruct("egymc_res_hist"));
  kernel_hist.save(writer->writeStruct("kernel_hist"));
  kernel_fit_hist.save(writer->writeStruct("kernel_fit_hist"));
  kernel_fit_counts_hist.save(writer->writeStruct("kernel_fit_counts_hist"));

  param_hist[0].save(writer->writeStruct("param0_hist"));
  param_hist[1].save(writer->writeStruct("param1_hist"));
  param_hist[2].save(writer->writeStruct("param2_hist"));
  param_hist[3].save(writer->writeStruct("param3_hist"));
  param_hist[4].save(writer->writeStruct("param4_hist"));

  param_fit_hist[0].save(writer->writeStruct("param0_fit_hist"));
  param_fit_hist[1].save(writer->writeStruct("param1_fit_hist"));
  param_fit_hist[2].save(writer->writeStruct("param2_fit_hist"));
  param_fit_hist[3].save(writer->writeStruct("param3_fit_hist"));
  param_fit_hist[4].save(writer->writeStruct("param4_fit_hist"));

  VSNSpaceOctaveH5IO io;
  io.writeHistogram(writer->writeStruct("kernel_nspace"),kernel_nspace);

  writer->writeStructCellVector("egy_hist",egy_hist);
  writer->writeStructCellVector("egy_fit_hist",egy_fit_hist);
  writer->writeStructCellVector("egy_fit2_hist",egy_fit2_hist);
}

VSEnergyKernelCalc::VSEnergyKernelCalc():
  m_log10_egybin(),
  m_log10_egylo(),
  m_log10_egyhi(),
  m_negy(),
  m_data(), m_data_ptr()
{

}

VSEnergyKernelCalc::~VSEnergyKernelCalc()
{
  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++) delete m_data[idata];
}

void VSEnergyKernelCalc::loadSimInfo(VSSimInfoData* sim_info)
{
  m_data_ptr = NULL;

  for(std::vector<Data*>::iterator itr = m_data.begin();
      itr != m_data.end(); ++itr)
    {
      double doff = fabs(sim_info->wobble_theta_deg-(*itr)->offset_deg);
      
      if(doff < 0.05) 
	{
	  m_data_ptr = *itr;
	  break;
	}
      else if(sim_info->wobble_theta_deg < (*itr)->offset_deg)
	{
	  m_data_ptr = new Data(sim_info->wobble_theta_deg,
				m_log10_egybin,m_log10_egylo,m_log10_egyhi);
	  m_data.insert(itr,m_data_ptr);
	  break;
	}
    }

  if(m_data_ptr == NULL)
    {
      m_data_ptr = new Data(sim_info->wobble_theta_deg,
			    m_log10_egybin,m_log10_egylo,m_log10_egyhi);
      m_data.push_back(m_data_ptr);
    }
}

void VSEnergyKernelCalc::loadHeader(const VSHeaderSimulationDatum& sim_header)
{
  vsassert(m_data_ptr != NULL);

//   for(std::vector< VSTableSimulationDatum >::const_iterator itr = 
// 	sim_header.tables.begin(); itr != sim_header.tables.end(); ++itr)
//     {
//       double area = std::pow(itr->sampling_radius_m,2)*M_PI;
//       double log10_egy = std::log10(itr->energy_tev);
      
//       if(m_data_ptr->egymc_area_hist.countForVal(log10_egy) == 0)
// 	m_data_ptr->egymc_area_hist.accumulate(log10_egy,area,0);

//       m_data_ptr->egymc_total_hist.accumulate(log10_egy,itr->event_count,0);
//     }
}

void VSEnergyKernelCalc::setEnergyBinning(double egybin, 
					   double egylo, double egyhi)
{
  if(m_log10_egybin == 0)
    {
      m_log10_egybin = egybin;
      m_log10_egylo = egylo;
      m_log10_egyhi = egyhi;
      m_negy = lround((egyhi-egylo)/egybin);

      m_sigma1_offset_hist = 
	VSSimple2DHist<double,double>(egybin,egylo,egyhi,0.04,0,1.8);
      m_bias1_offset_hist = 
	VSSimple2DHist<double,double>(egybin,egylo,egyhi,0.04,0,1.8);
      m_sigma2_offset_hist = 
	VSSimple2DHist<double,double>(egybin,egylo,egyhi,0.04,0,1.8);
      m_bias2_offset_hist = 
	VSSimple2DHist<double,double>(egybin,egylo,egyhi,0.04,0,1.8);
      m_alpha_offset_hist = 
	VSSimple2DHist<double,double>(egybin,egylo,egyhi,0.04,0,1.8);
    }
}

void VSEnergyKernelCalc::accumulate(double emc_tev, double log10_erec_tev)
{
  m_data_ptr->egymc_egy_hist.
    accumulate(std::log10(emc_tev),log10_erec_tev);
  m_data_ptr->egymc_egyerr_hist.
    accumulate(std::log10(emc_tev),log10_erec_tev-std::log10(emc_tev));
}

void VSEnergyKernelCalc::calcKernel()
{
  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    calcKernel(m_data[idata]);
}

void VSEnergyKernelCalc::calcKernel(Data* data)
{
  const unsigned negy = data->egymc_egy_hist.nXBins();
  data->egy_hist.resize(negy);
  data->egy_fit_hist.resize(negy);
  data->egy_fit2_hist.resize(negy);

  for(unsigned iegymc = 0; iegymc < negy; iegymc++)
    {
      data->egy_hist[iegymc] = 
	VSLimitedErrorsHist<double,double>(data->kernel_hist.yBinSize(),
					   data->kernel_hist.yLoLimit(),
					   data->kernel_hist.yHiLimit());
      data->egy_fit_hist[iegymc] = 
	VSLimitedErrorsHist<double,double>(data->kernel_hist.yBinSize(),
					   data->kernel_hist.yLoLimit(),
					   data->kernel_hist.yHiLimit());
      
      data->egy_fit2_hist[iegymc] = 
	VSLimitedErrorsHist<double,double>(data->kernel_hist.yBinSize(),
					   data->kernel_hist.yLoLimit(),
					   data->kernel_hist.yHiLimit());

      for(unsigned iegy = 0; iegy < data->egymc_egy_hist.nYBins(); iegy++)
	{
	  double x = data->egymc_egy_hist.yBinToCenter(iegy);
	  data->egy_hist[iegymc].
	    accumulate(x,data->egymc_egy_hist.count(iegymc,iegy));
	}
    }

  data->kernel_hist = data->egymc_egy_hist;
  data->kernel_hist.normalizeY();
  data->kernel_nspace = VSNSpace(data->kernel_hist);  
}

void VSEnergyKernelCalc::fit()
{
  const unsigned negy = m_sigma1_offset_hist.nXBins();
  std::vector<VSAFunction::Spline> bias1_spline(negy);
  std::vector<VSAFunction::Spline> bias2_spline(negy);
  std::vector<VSAFunction::Spline> alpha_spline(negy);

  VSAMath::Data<VSAAlgebra::Vec2D> sigma1_xys;
  VSAMath::Data<VSAAlgebra::Vec2D> sigma2_xys;

  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    {
      std::cout << std::string(__PRETTY_FUNCTION__)  
		<< ": Fitting Energy Response Fn " 
		<< " Offset  " << std::setw(15) << m_data[idata]->offset_deg
		<< std::endl;

      fit(m_data[idata]);

      for(unsigned iegy = 0; iegy < negy; iegy++)
	{
	  double log10_egy = m_sigma1_offset_hist.xBinToCenter(iegy);
	  double sigma1 = 
	    m_data[idata]->param_fit_hist[0].countForVal(log10_egy);
	  double sigma1_err = 
	    sqrt(m_data[idata]->param_fit_hist[0].varForVal(log10_egy));
	  double bias1 = 
	    m_data[idata]->param_fit_hist[1].countForVal(log10_egy);
	  double sigma2 = 
	    m_data[idata]->param_fit_hist[2].countForVal(log10_egy);
	  double sigma2_err = 
	    sqrt(m_data[idata]->param_fit_hist[2].varForVal(log10_egy));
	  double bias2 = 
	    m_data[idata]->param_fit_hist[3].countForVal(log10_egy);
	  double alpha =
	    m_data[idata]->param_fit_hist[4].countForVal(log10_egy);

	  VSAAlgebra::Vec2D xy(log10_egy,m_data[idata]->offset_deg);

	  // std::cout << std::setw(10) << iegy
	  // 	    << std::setw(15) << m_data[idata]->offset_deg
	  // 	    << std::setw(15) << sigma1
	  // 	    << std::setw(15) << sigma1_err
	  // 	    << std::setw(15) << sigma2
	  // 	    << std::setw(15) << sigma2_err
	  // 	    << std::endl;

	  if(sigma1_err > 1E-5 && std::isfinite(sigma1_err))
	    {
	      sigma1_xys.
		insert(VSAMath::
		       DataPoint<VSAAlgebra::Vec2D>(xy,sigma1,sigma1_err));
	    }

	  if(sigma2_err > 1E-5 && std::isfinite(sigma2_err))
	    {
	      sigma2_xys.
		insert(VSAMath::
		       DataPoint<VSAAlgebra::Vec2D>(xy,sigma2,sigma2_err));
	    }

	  bias1_spline[iegy].setPoint(m_data[idata]->offset_deg,bias1);
	  bias2_spline[iegy].setPoint(m_data[idata]->offset_deg,bias2);
	  alpha_spline[iegy].setPoint(m_data[idata]->offset_deg,alpha);
	}

    }  

  VSAMath::LocalRegression2D<VSAAlgebra::Vec2D> sigma1_lr(sigma1_xys,0.7,0.7,1);
  VSAMath::LocalRegression2D<VSAAlgebra::Vec2D> sigma2_lr(sigma2_xys,0.7,0.7,1);

  for(unsigned iegy = 0; iegy < negy; iegy++)
    {
      double log10_egy = m_sigma1_offset_hist.xBinToCenter(iegy);

      bias1_spline[iegy].spline();
      bias2_spline[iegy].spline();
      alpha_spline[iegy].spline();

      for(unsigned iy = 0; iy < m_sigma1_offset_hist.nYBins(); iy++)
	{
	  double offset = m_sigma1_offset_hist.yBinToCenter(iy);

	  VSAAlgebra::Vec2D xy(log10_egy,offset);

	  try
	    {
	      m_sigma1_offset_hist.setBin(iegy,iy,sigma1_lr.val(xy));
	      m_sigma2_offset_hist.setBin(iegy,iy,sigma2_lr.val(xy));
	    }
	  catch(const std::exception& e)
	    {
	      std::cerr << std::string(__PRETTY_FUNCTION__) + ": " 
			<< " E = " << std::setw(20) << log10_egy 
			<< " Offset = " << std::setw(20) << offset << std::endl
			<< e.what() << std::endl;
	    }

	  m_bias1_offset_hist.setBin(iegy,iy,bias1_spline[iegy].val(offset));
	  m_bias2_offset_hist.setBin(iegy,iy,bias2_spline[iegy].val(offset));
	  m_alpha_offset_hist.setBin(iegy,iy,alpha_spline[iegy].val(offset));
	}
    }

  m_sigma1_offset_nspace = VSNSpace(m_sigma1_offset_hist);  
  m_bias1_offset_nspace = VSNSpace(m_bias1_offset_hist);  
  m_sigma2_offset_nspace = VSNSpace(m_sigma2_offset_hist);  
  m_bias2_offset_nspace = VSNSpace(m_bias2_offset_hist);  
  m_alpha_offset_nspace = VSNSpace(m_alpha_offset_hist); 
}

void VSEnergyKernelCalc::fit(Data* data)
{
  const unsigned negy = data->egy_hist.size();
  const double dloge = data->egymc_egy_hist.xBinSize();

  std::vector< double > counts(negy);
  std::vector< double > mean(negy);
  std::vector< double > rms(negy);

  VSEnergyKernelFn kfn(data->egymc_egy_hist.yBinSize(),negy,dloge,
		       data->egymc_egy_hist.xLoLimit());

  VSAMath::Data<VSACoord::CoordND> fit_data2;

  for(unsigned iegymc = 0; iegymc < negy; iegymc++)
    {
      double egy = data->egymc_egy_hist.binToCenterX(iegymc);
	    

      VSAMath::Data<VSACoord::CoordND> fit_data;

      double s = 0;
      double emin = data->egy_hist[iegymc].hiLimit();
      double emax = data->egy_hist[iegymc].loLimit();

      VSSimpleStat2<double,double> stat;

      for(VSLimitedErrorsHist<double,double>::iterator itr = 
	    data->egy_hist[iegymc].begin(); itr != 
	    data->egy_hist[iegymc].end(); ++itr)
	{
	  VSACoord::CoordND c(2);
	  c[0] = itr->center();
	  c[1] = egy;
	  VSAMath::DataPoint<VSACoord::CoordND> dp(c,itr->count(),itr->err());

	  stat.accumulate(itr->center(),itr->count());

	  if(itr->count() > 0 && itr->center() < emin) emin = itr->center();
	  if(itr->count() > 0 && itr->center() > emax) emax = itr->center();
	  s+= itr->count();
	}  

      emin = stat.mean() - 3*stat.dev();
      emax = stat.mean() + 3*stat.dev();

      counts[iegymc] = s;
      mean[iegymc] = stat.mean();
      rms[iegymc] = stat.dev();


      double bias = mean[iegymc]-egy;
      double mse = std::pow(bias,2)+stat.var();
      double mse_var = 4*std::pow(bias,2)*stat.mean_var() + stat.var_var();
      double res = std::pow(10,mean[iegymc]-egy)*log(10)*stat.dev();
      double res_var = 
	std::pow(res,2)*(stat.dev_var()/stat.var() +
			 std::pow(log(10),2)*stat.mean_var());

      data->egymc_bias_hist.accumulate(egy,bias,stat.mean_var());
      data->egymc_dev_hist.accumulate(egy,rms[iegymc],stat.dev_var());
      data->egymc_res_hist.accumulate(egy,res,res_var);
      data->egymc_biasdev_hist.accumulate(egy,bias,stat.var()); 
      data->egymc_mse_hist.accumulate(egy,mse,mse_var); 

      kfn.setNorm(iegymc,s);

      if(s < 10) continue;

      for(VSLimitedErrorsHist<double,double>::iterator itr = 
	    data->egy_hist[iegymc].begin(); 
	  itr != data->egy_hist[iegymc].end(); ++itr)
	{
	  VSACoord::CoordND c(2);
	  c[0] = itr->center();
	  c[1] = egy;
	  VSAMath::DataPoint<VSACoord::CoordND> dp(c,itr->count(),itr->err());
	  if(itr->center() > emin - 0.25 && itr->center() < emax + 0.25) 
	    {
	      fit_data.insert(dp);
	      fit_data2.insert(dp);
	    }
	}        

      VSEnergyResponseFn fn(data->egy_hist[iegymc].binSize(),s);

      typedef VSAFunction::LnPoissonLikelihood<VSEnergyResponseFn,
	VSACoord::CoordND> LnLFn;
      VSAMath::NLFitter< LnLFn >* fitter = 
	VSAMath::NLFitterFactory::createLnL(&fn,fit_data);
      
      VSAAlgebra::VecND p(5);

      p(0) = stat.dev();
      p(1) = stat.mean();
      p(2) = stat.dev();
      p(3) = stat.mean();
      p(4) = 0.8;

      try
	{
	  fitter->setLoBound(4,0);
	  fitter->setHiBound(4,1);
	  fitter->hold(2,0.10);
	  fitter->hold(4,0.8);
	  fitter->initialize(p);
	  fitter->fit();
	  fitter->free(2);
	  fitter->initialize(fitter->param());
	  fitter->fit();
	}
      catch(const std::exception& e)
	{
	  std::cerr << e.what() << std::endl;
	}
      
      p = fitter->param();
      VSAAlgebra::MatrixND cov = fitter->cov();

#ifdef DEBUG
      std::cout 
	<< std::setw(5) << iegymc
	<< std::setw(15) << egy 
	<< std::setw(15) << s
	<< std::setw(15) << p(0) 
	<< std::setw(15) << p(1) 
	<< std::setw(15) << p(2) 
	<< std::setw(15) << p(3) 
	<< std::setw(15) << p(4) 
	<< std::setw(15) << fitter->chi2()/(double)fitter->ndf() << std::endl;
#endif

      delete fitter;


      for(unsigned i = 0; i < data->param_hist.size(); i++)
	{
	  if(std::isfinite(cov(i,i)))
	    {
	      data->param_hist[i].accumulate(egy,p(i),cov(i,i));
	    }
	}

      data->egy_fit_hist[iegymc].fill(0);

      for(VSLimitedErrorsHist<double,double>::iterator itr = 
	    data->egy_fit_hist[iegymc].begin(); itr !=
	    data->egy_fit_hist[iegymc].end(); ++itr)
	{
	  VSACoord::CoordND c(2);
	  c[0] = itr->center();
	  c[1] = egy;
	  data->egy_fit_hist[iegymc].setBin(itr->bin(),fn.val(c,p),0);
	}

    }

  for(unsigned i = 0; i < data->param_fit_hist.size(); i++)
    data->param_fit_hist[i].fill(0.);

  VSAAlgebra::VecND plast;
  for(unsigned iegymc = 0; iegymc < negy; iegymc++)
    {
      double egy = data->egymc_egy_hist.binToCenterX(iegymc);
      VSAMath::Data<VSACoord::CoordND> fd;

#ifdef DEBUG
      std::cout << std::string(79,'-') << std::endl;
      std::cout << iegymc << " " << egy << std::endl;
#endif

      double s = counts[iegymc];

      int iegylo = iegymc;
      int iegyhi = iegymc;

      VSAMath::Data<double> mean_data;
      VSAMath::Data<double> rms_data;

      while(s < 2000 || iegyhi-iegylo < 8)
	{
	  if(iegylo > 0) s += counts[--iegylo];
	  if(iegyhi+1 < (int)negy) s += counts[++iegyhi];

	  if(iegylo == 0 && iegyhi+1 == (int)negy) break;
	}
     
      for(int iegy = iegylo; iegy <= iegyhi; iegy++)
	{
	  double egy2 = data->egymc_egy_hist.binToCenterX(iegy);

	  if(counts[iegy]>1 && rms[iegy] > 0)
	    {
	      double mean_err = rms[iegy]/sqrt(counts[iegy]);
	      double rms_err = sqrt(2.)*mean_err;
	      
	      mean_data.insert(VSAMath::DataPoint<double>(egy2,
							  mean[iegy],
							  mean_err));

	      rms_data.insert(VSAMath::DataPoint<double>(egy2,
							 rms[iegy],
							 rms_err));
	    }
	}

      VSAAlgebra::VecND p_mean;
      VSAMath::PolyFit::fit(2, mean_data, p_mean);
	      
      VSAAlgebra::VecND p_rms;
      VSAMath::PolyFit::fit(1, rms_data, p_rms);

      VSSimpleStat2<double,double> stat;
      for(VSAMath::Data<VSACoord::CoordND>::iterator itr = fit_data2.begin();
	  itr != fit_data2.end(); ++itr)
	{
	  int iegy = data->egymc_egy_hist.xValToBin(itr->x[1]);

	  if(iegy >= iegylo && iegy <= iegyhi)
	    {
	      fd.insert(*itr);
	      stat.accumulate(itr->x[1],itr->y);
	    }
	}


      // double mean_egy = stat.mean();
      // std::cout << "KERNEL FITTING " 
      // 		<< std::setw(15) << egy 
      // 		<< std::setw(15) << mean_egy 
      // 		<< std::setw(10) << iegylo 
      // 		<< std::setw(10) << iegyhi 
      // 		<< std::setw(15) << s
      // 		<< std::endl;

      typedef VSAFunction::LnPoissonLikelihood<VSEnergyKernelFn,
	VSACoord::CoordND> LnLFn;
      VSAMath::NLFitter< LnLFn >* fitter = 
	VSAMath::NLFitterFactory::createLnL(&kfn,fd);

      try
	{

	  // PARAMETERS
	  // 0-2: Sigma1
	  // 3-5: Mean1
	  // 6-7: Sigma2
	  // 8-10: Mean2
	  // 11: Alpha


	  if(plast.ndim() == 0)
	    {
	      fitter->hold(0,p_rms[0]);
	      fitter->hold(1,p_rms[1]);
	      fitter->hold(2,0);    
	      fitter->hold(3,p_mean[0]);    
	      fitter->hold(4,p_mean[1]);    
	      fitter->hold(5,p_mean[2]);
	      fitter->hold(6,0.1);
	      fitter->hold(7,0.);
	      fitter->hold(8,p_mean[0]);
	      fitter->hold(9,p_mean[1]);
	      fitter->hold(10,p_mean[2]);
	      fitter->hold(11,0.8);	  

	      fitter->setLoBound(11,0.6);
	      fitter->setHiBound(11,1);

	      // --------------------------------------------------------------
	      // Fit the first Gaussian

	      // std::cout << "-----------------------------------------------" 
	      // 	    << std::endl;
	      // for(unsigned ip = 0; ip < fitter->nparam(); ip++)
	      //   std::cout << std::setw(5) << ip
	      // 	      << std::setw(15) << fitter->param(ip)
	      // 	      << std::setw(15) << fitter->err(ip)
	      // 	      << std::endl;
	  
	      fitter->free(3);
	      fitter->free(4);
	      fitter->free(5);
	      fitter->initialize();
	      fitter->fit();
	  
	  
	      // std::cout << "-----------------------------------------------" 
	      // 	    << std::endl;
	      // for(unsigned ip = 0; ip < fitter->nparam(); ip++)
	      //   std::cout << std::setw(5) << ip
	      // 	      << std::setw(15) << fitter->param(ip)
	      // 	      << std::setw(15) << fitter->err(ip)
	      // 	      << std::endl;


	      fitter->free(0);
	      fitter->free(1);
	      //	      fitter->free(2);
	      fitter->initialize(fitter->param());
	      fitter->fit();
      
	      // ---------------------------------------------------------------
	      // Fit the second Gaussian

	      // std::cout << "-----------------------------------------------" 
	      // 	    << std::endl;
	      // for(unsigned ip = 0; ip < fitter->nparam(); ip++)
	      //   std::cout << std::setw(5) << ip
	      // 	      << std::setw(15) << fitter->param(ip)
	      // 	      << std::setw(15) << fitter->err(ip)
	      // 	      << std::endl;
	  
	      fitter->free(8);
	      fitter->free(9);
	      fitter->free(10);
	      fitter->set(8,fitter->param(3));
	      fitter->set(9,fitter->param(4));
	      fitter->set(10,fitter->param(5));
	      fitter->initialize(fitter->param());
	      fitter->fit();
	  
	      // std::cout << "-----------------------------------------------" 
	      // 	    << std::endl;
	      // for(unsigned ip = 0; ip < fitter->nparam(); ip++)
	      //   std::cout << std::setw(5) << ip
	      // 	      << std::setw(15) << fitter->param(ip)
	      // 	      << std::setw(15) << fitter->err(ip)
	      // 	      << std::endl;
	  
	      fitter->free(6);
	      fitter->free(7);
	      fitter->initialize(fitter->param());	
	      fitter->fit();
	    }
	  else
	    {
	      //	      fitter->free();
	      fitter->hold(0);
	      fitter->hold(1);
	      fitter->hold(2);
	      fitter->hold(6);
	      fitter->hold(7);
	      fitter->hold(8);
	      fitter->hold(9);
	      fitter->hold(10);
	      fitter->hold(11,0.8);
	      fitter->initialize(plast);	      
	      fitter->fit();

	      // std::cout << "-----------------------------------------------" 
	      // 		<< std::endl;
	      // for(unsigned ip = 0; ip < fitter->nparam(); ip++)
	      // 	std::cout << std::setw(5) << ip
	      // 		  << std::setw(15) << fitter->param(ip)
	      // 		  << std::setw(15) << fitter->err(ip)
	      // 		  << std::endl;


	      fitter->free();
	      fitter->hold(2);
	      fitter->hold(6);
	      fitter->hold(7);
	      fitter->hold(8);
	      fitter->hold(9);
	      fitter->hold(10);
	      fitter->hold(11,0.8);
	      fitter->initialize(fitter->param());	      
	      fitter->fit();

	      // std::cout << "-----------------------------------------------" 
	      // 		<< std::endl;
	      // for(unsigned ip = 0; ip < fitter->nparam(); ip++)
	      // 	std::cout << std::setw(5) << ip
	      // 		  << std::setw(15) << fitter->param(ip)
	      // 		  << std::setw(15) << fitter->err(ip)
	      // 		  << std::endl;

	      fitter->free();
	      fitter->hold(2);
	      fitter->hold(6);
	      fitter->hold(7);
	      fitter->hold(10);
	      fitter->hold(11,0.8);
	      fitter->initialize(fitter->param());	      
	      fitter->fit();

	      // std::cout << "-----------------------------------------------" 
	      // 		<< std::endl;
	      // for(unsigned ip = 0; ip < fitter->nparam(); ip++)
	      // 	std::cout << std::setw(5) << ip
	      // 		  << std::setw(15) << fitter->param(ip)
	      // 		  << std::setw(15) << fitter->err(ip)
	      // 		  << std::endl;


	      fitter->free();
	      fitter->hold(2);	      
	      fitter->hold(11,0.8);
	      fitter->initialize(fitter->param());
	      fitter->fit();
	    }
	 
	  // fitter->free(11);
	  // fitter->initialize(fitter->param());	
	  // fitter->fit();
	  
	}
      catch(const std::exception& e)
	{
	  std::cerr << e.what() << std::endl;
	}

      VSAAlgebra::VecND param = fitter->param();
      VSAAlgebra::MatrixND cov = fitter->cov();

      param(0) = fabs(param(0));
      param(2) = fabs(param(2));

      plast = param;

      delete fitter;

      VSAAlgebra::VecND coeff(5);
      VSAAlgebra::VecND coeff_var(5);



      kfn.coeff(egy,param,coeff);
      kfn.coeff_var(egy,param,cov,coeff_var);
      for(unsigned i = 0; i < coeff.ndim(); i++)
	{	  
	  if(i==0 || i==2) coeff(i) = fabs(coeff(i));
	  data->param_fit_hist[i].accumulate(egy,coeff(i),coeff_var(i));
	}

#ifdef DEBUG
      std::cout << "-----------------------------------------------" 
		<< std::endl;
      for(unsigned ip = 0; ip < param.ndim(); ip++)
	std::cout << std::setw(5) << ip
		  << std::setw(15) << param(ip)
		  << std::setw(15) << sqrt(cov(ip,ip))
		  << std::endl;

      std::cout << "-----------------------------------------------" 
		<< std::endl;

      for(unsigned i = 0; i < coeff.ndim(); i++)
	std::cout << std::setw(15) << coeff(i) 
		  << std::setw(15) << sqrt(coeff_var(i)) 
		  << std::endl;
#endif
    }

  unsigned iegymin = 0;
  while(counts[iegymin] < 10) iegymin++;
  for(unsigned iegymc = 0; iegymc < iegymin; iegymc++)
    {
      for(unsigned i = 0; i < data->param_fit_hist.size(); i++)
  	{
  	  double c = data->param_fit_hist[i].count(iegymin);
  	  double var = data->param_fit_hist[i].var(iegymin);
  	  data->param_fit_hist[i].setBin(iegymc,c,var);
  	}
    }


  for(unsigned iegymc = 0; iegymc < negy; iegymc++)
    {
      double egy = data->egymc_egy_hist.binToCenterX(iegymc);

      VSAAlgebra::VecND coeff(5);
      

      for(unsigned i = 0; i < coeff.ndim(); i++)
	coeff(i) = data->param_fit_hist[i].count(iegymc);
      
      VSEnergyResponseFn fn(data->egy_hist[iegymc].binSize(),counts[iegymc]);

      for(VSLimitedErrorsHist<double,double>::iterator itr = 
	    data->egy_fit2_hist[iegymc].begin(); itr !=
	    data->egy_fit2_hist[iegymc].end(); ++itr)
	{
 	  VSACoord::CoordND c(2);
 	  c[0] = itr->center();
 	  c[1] = egy;
	  data->egy_fit2_hist[iegymc].setBin(itr->bin(),
					     fn.val(c,coeff),0);

	  data->kernel_fit_hist.setBin(iegymc,itr->bin(),
				       fn.val(c,coeff,1));
	  data->kernel_fit_counts_hist.setBin(iegymc,itr->bin(),
					      fn.val(c,coeff));
	}
    }

}

void VSEnergyKernelCalc::save(VSOctaveH5WriterStruct* writer) const
{
  const unsigned ndata = m_data.size();

  VSOctaveH5WriterCellVector* wc = writer->writeCellVector("offset", ndata);
  vsassert(wc);

  for(unsigned idata = 0; idata < ndata; idata++)
    {      
      VSOctaveH5WriterStruct* ws = wc->writeStruct(idata);
      vsassert(ws); 
      m_data[idata]->save(ws);
      delete ws;
    }
  delete wc;

  m_sigma1_offset_hist.save(writer->writeStruct("sigma1_offset_hist"));
  m_bias1_offset_hist.save(writer->writeStruct("bias1_offset_hist"));
  m_sigma2_offset_hist.save(writer->writeStruct("sigma2_offset_hist"));
  m_bias2_offset_hist.save(writer->writeStruct("bias2_offset_hist"));
  m_alpha_offset_hist.save(writer->writeStruct("alpha_offset_hist"));
}

