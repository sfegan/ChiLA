#include <fstream>

#include <VSALinearLeastSquares.hpp>
#include <VSANonlinearFitting.hpp>
#include <VSAQuadrature.hpp>
#include <VSBkgndModel.hpp>
#include <VSAFunction.hpp>

using namespace VERITAS;
using namespace VERITAS::VSAFunction;

// ============================================================================
// VSAcceptanceModel
// ============================================================================
VSAcceptanceModel::VSAcceptanceModel(double bin_size_deg, double offset_max): 
  VSAFunction::ParamFn<VSModelCoord>(),
  m_bin_size(bin_size_deg),
  m_domega(std::pow(bin_size_deg,2)), 
  m_offset_max(offset_max), m_obs_xy()
{ 
  
}

void VSAcceptanceModel::setObs(const VSAAlgebra::Vec2D& xy)
{
  m_obs_xy.push_back(xy);

  VSAAlgebra::VecND p = param();
  setNParam(nparm()+1);

  for(unsigned ip = 0; ip < p.ndim(); ip++)
    setParam(ip,p[ip]);

  setParam(nparm()-1,0.0);
}

// ============================================================================
// VSAcceptanceModelFixed
// ============================================================================
VSAcceptanceModelFixed::
VSAcceptanceModelFixed(VSAcceptanceModel* model,
		       const VSSimple2DHist<double,double>& counts_hist):
  VSAcceptanceModel(*model)
{
  const unsigned nobs = model->nobs();
  for(unsigned iobs = 0; iobs < nobs; iobs++)
    {
      VSSimple2DHist<double,double> h = counts_hist;
      h.clear();
      
      for(VSSimple2DHist<double,double>::iterator itr = h.begin(); itr !=
	    h.end(); ++itr)
	{
	  VSModelCoord c(itr->x(), itr->y(), iobs);
	  h.setBin(itr->bin(),model->val(c));
	}

      m_fn_val.push_back(h);
    }
}

double VSAcceptanceModelFixed::val(const VSAAlgebra::Vec2D& xy) const
{ 
  const unsigned nobs = m_fn_val.size();
  double val = 0;
  for(unsigned iobs = 0; iobs < nobs; iobs++)
    {
      VSModelCoord c(xy,iobs);
      val += VSAcceptanceModelFixed::val(c);
    }
  return val;
}

double VSAcceptanceModelFixed::val(const VSModelCoord& x) const
{
  return m_fn_val[x.iobs()].countForVal(x.x(),x.y());
}
    
double VSAcceptanceModelFixed::val(const VSModelCoord& x, 
				   const VSAAlgebra::VecND& a) const
{
  return m_fn_val[x.iobs()].countForVal(x.x(),x.y());
}

void VSAcceptanceModelFixed::dyda(const VSModelCoord& x, 
				  VSAAlgebra::VecND& dyda) const
{

}

void VSAcceptanceModelFixed::dyda(const VSModelCoord& x, 
				  const VSAAlgebra::VecND& a,
				  VSAAlgebra::VecND& dyda) const
{
  dyda.resize(nparm());
  dyda.clear();
}

double VSAcceptanceModelFixed::getFoVCounts(const VSAAlgebra::Vec2D& xy) 
  const 
{ 
  return 0;
}

double VSAcceptanceModelFixed::getFoVAcceptance(const VSAAlgebra::Vec2D& xy) 
  const 
{ 
  return 0;
}

void VSAcceptanceModelFixed::integrate(const VSAAlgebra::VecND& a,
				       const VSAAlgebra::MatrixND& cov,
				       double& integral,
				       double& integral_err) const
{
  integral = 0;
  integral_err = 0;
}

double VSAcceptanceModelFixed::integrate(const VSAAlgebra::Vec2D& xy, 
					 double R1, double R2) const
{
  return 0;
}

// ============================================================================
// VSAcceptanceModelPoly
// ============================================================================
VSAcceptanceModelPoly::VSAcceptanceModelPoly(const std::string& param,
					     double bin_size_deg, 
					     double offset_max): 
  VSAcceptanceModel(bin_size_deg,offset_max), m_fn() 
{ 

}    

double VSAcceptanceModelPoly::val(const VSModelCoord& x) const
{
  return val(x,m_fn.param());
}

double VSAcceptanceModelPoly::val(const VSModelCoord& x, 
				  const VSAAlgebra::VecND& a) const
{
  double val = 0;
  const unsigned nobs = m_obs_xy.size();
  for(unsigned iobs = 0; iobs < nobs; iobs++)
    { 
      VSAAlgebra::Vec2D xy(x.x(),x.y());
      xy -= m_obs_xy[iobs];
      val += m_fn.val(VSACoord::Coord2D(xy),a);
    }
  return val;
}

void VSAcceptanceModelPoly::dyda(const VSModelCoord& x, 
				 VSAAlgebra::VecND& dyda) const
{
  VSAcceptanceModelPoly::dyda(x,param(),dyda);
}

void VSAcceptanceModelPoly::dyda(const VSModelCoord& x, 
			   const VSAAlgebra::VecND& a,
			   VSAAlgebra::VecND& dyda) const
{
  dyda.resize(nparm());
  dyda.clear();
  
  const unsigned nobs = m_obs_xy.size();
  for(unsigned iobs = 0; iobs < nobs; iobs++)
    { 
      VSAAlgebra::Vec2D xy(x.x(),x.y());
      xy -= m_obs_xy[iobs];
      VSAAlgebra::VecND tmp;
      m_fn.dyda(VSACoord::Coord2D(xy),a,tmp);
      dyda += tmp;
    }
}

double VSAcceptanceModelPoly::getFoVAcceptance(const VSAAlgebra::Vec2D& xy) 
  const
{ 
  return m_fn.val(VSModelCoord(xy,0));
}

// void VSAcceptanceModelPoly::initialize(const VSAcceptanceData& data)
// {
//   const unsigned max_ncoeff = 6; 
//   const double dr2 = 0.1;
//   const double domega = M_PI*dr2;

//   VSLimitedErrorsHist<double,double> hist(dr2,0,2.0);

//   VSAMath::Data<double> data_tmp;

//   for(VSAMath::Data<VSModelCoord>::const_iterator itr = data.begin();
//       itr != data.end(); ++itr)
//     {
//       double r2 = std::pow(itr->x.r(),2);
//       hist.accumulate(r2,itr->y,itr->y);
//     }

//   for(VSLimitedErrorsHist<double,double>::iterator itr = 
//  	hist.begin(); itr != hist.end(); ++itr)
//     {
//       vsassert(std::isfinite(itr->count()));
//       if(itr->count() == 0) continue;

//       VSAMath::DataPoint<double> p(itr->val(),itr->count(),itr->err());
//       data_tmp.insert(p);
//     }

//   VSAAlgebra::VecND param;
//   double chi2 = 0;
//   double ndf = 0;

//   for(unsigned ncoeff = 1; ncoeff < max_ncoeff; ncoeff++)
//     {
//       m_fn.set(ncoeff);
//       ndf = data.size()-ncoeff-1;
//       try
// 	{
// 	  chi2 = VSAMath::PolyFit::fit(ncoeff,data_tmp,param);
// 	}
//       catch(const std::string& s)
// 	{
// 	  std::cout << s << std::endl;
// 	  continue;
// 	}

//       if(chi2/ndf < 1.5)
// 	break;
//     }

//   for(unsigned ip = 0; ip < param.size(); ip++)
//     param[ip]=param[ip]/domega*m_domega;

//   m_fn.setParam(param);
//}

void VSAcceptanceModelPoly::integrate(const VSAAlgebra::VecND& a,
				      const VSAAlgebra::MatrixND& cov,
				      double& integral,
				      double& integral_err) const
{
  integral = integrate(VSAAlgebra::Vec2D(0,0),m_offset_max);
  integral_err = 0;
}

double VSAcceptanceModelPoly::integrate(const VSAAlgebra::Vec2D& xy, 
					double R1, double R2) const
{
  return integrate(xy,R2) - integrate(xy,R1);
}

double VSAcceptanceModelPoly::integrate(const VSAAlgebra::Vec2D& xy, 
					double R) const
{
  double r = xy.norm();

  vsassert(m_fn.nparm() <= 7);

  const double d2 = std::pow(r,2);
  const double d4 = std::pow(r,4);
  const double d6 = std::pow(r,6);
  const double d8 = std::pow(r,8);
  const double d10 = std::pow(r,10);
  const double d12 = std::pow(r,12);

  const double R2 = std::pow(R,2);
  const double R4 = std::pow(R,4);
  const double R6 = std::pow(R,6);
  const double R8 = std::pow(R,8);
  const double R10 = std::pow(R,10);
  const double R12 = std::pow(R,12);

  const double b0 = 0.5;
  const double b1 = (1./2.)*d2 + (1./4.)*R2;
  const double b2 = (1./2.)*d4 + d2*R2 + (1./6.)*R4;
  const double b3 = (1./2.)*d6 + (18./8.)*d4*R2 + (12./8.)*d2*R4 + (1./8.)*R6;
  const double b4 = (1./2.)*d8 + (40./10.)*d6*R2 + (60./10.)*d4*R4 + 
    (20./10.)*d2*R6 + (1./10.)*R8;
  const double b5 = (1./2.)*d10 + (75./12.)*d8*R2 + (200./12.)*d6*R4 + 
    (150./12.)*d4*R6 + (30./12.)*d2*R8 + (1./12.)*R10;
  const double b6 = (1./2.)*d12 + (126./14.)*d10*R2 + (525./14.)*d8*R4 + 
    (700./14.)*d6*R6 + (315./14.)*d4*R8 + (42./14.)*d2*R10 + (1./14.)*R12;
  
  const double b[7] = { b0, b1, b2, b3, b4, b5, b6 };
  double sum = 0;

  for(unsigned icoeff = 0; icoeff < m_fn.nparm(); icoeff++)
    sum += M_PI*R2*m_fn.param(icoeff)*b[icoeff];
    
  return sum;
}

// ============================================================================
// VSAcceptanceModelBessel2
// ============================================================================
VSAcceptanceDataBessel2::
VSAcceptanceDataBessel2(double offset_max,
			const VSAAlgebra::VecND& p,
			const VSAAlgebra::MatrixND& cov):
  VSAcceptanceData(offset_max,"bessel2",p,cov), 
  c_param(),
  nparam_shape(),
  m0_param(),
  mn_param(10),
  m0_param_hist(1,0,10),
  mx_param_hist(),
  my_param_hist(),
  mr_param_hist(),
  mphi_param_hist()
{ 
  const unsigned m = 2;

  for(unsigned i = 0; i < m; i++)
    {
      mx_param_hist.push_back(VSLimitedErrorsHist<double, double>(1,0,10));
      my_param_hist.push_back(VSLimitedErrorsHist<double, double>(1,0,10));
      mr_param_hist.push_back(VSLimitedErrorsHist<double, double>(1,0,10));
      mphi_param_hist.push_back(VSLimitedErrorsHist<double, double>(1,0,10));
    }
}

bool VSAcceptanceDataBessel2::load(VSOctaveH5ReaderStruct* reader)
{
  return true;
}

void VSAcceptanceDataBessel2::save(VSOctaveH5WriterStruct* writer) const
{
  VSAcceptanceData::save(writer);

  writer->writeScalar("c_param",c_param);
  writer->writeScalar("nparam_shape",nparam_shape);

  writer->writeCompositeVector("m0_param",m0_param);

  VSOctaveH5WriterCellVector* wc = 
    writer->writeCellVector("mn_param", mn_param.size());
  for(unsigned ip = 0; ip < mn_param.size(); ip++)
    {
      wc->writeCompositeVector(ip,mn_param[ip]);
    }
  delete wc;

  m0_param_hist.save(writer->writeStruct("m0_param_hist"));

  const unsigned m = mx_param_hist.size();
  for(unsigned im = 0; im < m; im++)
    {
      std::ostringstream os;
      os << "m" << im+1;
      mx_param_hist[im].save(writer->writeStruct(os.str() + "x_param_hist"));
      my_param_hist[im].save(writer->writeStruct(os.str() + "y_param_hist"));
      mr_param_hist[im].save(writer->writeStruct(os.str() + "r_param_hist"));
      mphi_param_hist[im].
	save(writer->writeStruct(os.str() + "phi_param_hist"));
    }
}

VSAcceptanceModelBessel2::VSAcceptanceModelBessel2(const std::string& param,
						   double bin_size_deg, 
						   double offset_max):
  VSAcceptanceModel(bin_size_deg,offset_max), m_fn(offset_max)
{  
  std::vector< std::pair< unsigned, unsigned > > mn;
  VSDataConverter::fromString(mn, param);

  if(mn.size() == 0)
    {
      mn.push_back(std::make_pair(0,0));
      mn.push_back(std::make_pair(0,1));
      mn.push_back(std::make_pair(0,2));
      mn.push_back(std::make_pair(0,3));
    }

  m_fn.set(mn);
  setNParam(m_fn.nparm());
  for(unsigned ip = 0; ip < nparm(); ip++) setParam(ip,0);
}

void VSAcceptanceModelBessel2::fillData(VSAcceptanceDataBessel2* data)
{
  const unsigned nm = data->mx_param_hist.size();
  std::vector< std::vector< double > > mcov(nm);

  data->c_param = m_fn.getC();
  data->nparam_shape = m_fn.nparm();

  std::vector< std::map<unsigned,double> > mx_param(10);
  std::vector< std::map<unsigned,double> > mx_param_cov(10);
  std::vector< std::map<unsigned,double> > my_param(10);
  std::vector< std::map<unsigned,double> > my_param_cov(10);
  std::vector< std::map<unsigned,double> > mxy_param_cov(10);

  for(FourierBesselSeries2::const_iterator itr = m_fn.begin(); itr !=
	m_fn.end(); ++itr)
    {
      unsigned absm = std::abs(itr->m); 
      double param_norm = data->param(itr->ip)/data->c_param;
      double cov_norm = data->param_cov(itr->ip,itr->ip)/
	std::pow(data->c_param,2);

      if(itr->m==0)
	{
	  VSAcceptanceDataBessel2::Param p;
	  p.c = data->param(itr->ip);
	  p.c_err = data->param_err(itr->ip);
	  p.n = itr->n;
	  data->m0_param.push_back(p);
	}
      else if(itr->m > 0)
	{
	  mx_param[absm-1][itr->n] = param_norm;
	  mx_param_cov[absm-1][itr->n] = cov_norm;

	  data->mx_param_hist[absm-1].setBin(itr->n,param_norm,cov_norm);
	  for(FourierBesselSeries2::const_iterator itr2 = m_fn.begin(); 
	      itr2 != m_fn.end(); ++itr2)
	    if(itr2->m == -itr->m) 
	      mxy_param_cov[absm-1][itr->n] = 
		data->param_cov(itr->ip,itr2->ip)/std::pow(data->c_param,2);
	}
      else if(itr->m < 0)
	{
	  my_param[absm-1][itr->n] = param_norm;
	  my_param_cov[absm-1][itr->n] = cov_norm;

	  data->my_param_hist[absm-1].setBin(itr->n,param_norm,cov_norm);
	}       
    }

  for(unsigned m = 0; m < nm; m++)
    {

      for(std::map<unsigned,double>::iterator itr = mx_param[m].begin();
	  itr != mx_param[m].end(); ++itr)
	{
	  unsigned n = itr->first;

	  vsassert(my_param[m].find(itr->first) != my_param[m].end());
	  double mx = mx_param[m][n];
	  double mx_var = mx_param_cov[m][n];
	  double my = my_param[m][n];
	  double my_var = my_param_cov[m][n];	  

	  double mphi = atan2(my,mx);
	  double mr = sqrt(std::pow(mx,2)+std::pow(my,2));
	  double mr_var = 0;
	  double mphi_var = 0;

	  if(mr > 0)
	    {
	      mr_var = std::pow(mx/mr,2)*mx_var + std::pow(my/mr,2)*my_var +
		2*mxy_param_cov[m][n]*mx*my/std::pow(mr,2);
	      mphi_var =
		std::pow(mx/(mr*mr),2)*mx_var + 
		std::pow(my/(mr*mr),2)*my_var -
		2*mx*my/std::pow(mr,4)*mxy_param_cov[m][n];
	    }

	  mphi *= 180/M_PI;
	  mphi_var *= std::pow(180/M_PI,2);

	  VSAcceptanceDataBessel2::AzimuthalParam azp;

	  azp.n = n;
	  azp.cx = mx;
	  azp.cx_err = sqrt(mx_var);
	  azp.cy = sqrt(my);
	  azp.cy_err = sqrt(my_var);
	  azp.cr = mr;
	  azp.cr_err = sqrt(mr_var);
	  azp.cphi = mphi;
	  azp.cphi_err = sqrt(mphi_var);

	  data->mn_param[m].push_back(azp);

	  data->mr_param_hist[m].accumulate(0,mr,mr_var);
	  data->mphi_param_hist[m].accumulate(0,mphi,mphi_var);
	}
      
    }
  


//   for(unsigned m = 0; m < nm; m++)
//     {
//       for(VSLimitedErrorsHist<double, double>::iterator itr = 
// 	    data->mx_param_hist[m].begin(); 
// 	  itr != data->mx_param_hist[m].end(); ++itr)
// 	{
// 	  double mx = data->mx_param_hist[m].count(itr->bin());
// 	  double mx_var = data->mx_param_hist[m].var(itr->bin());
// 	  double my = data->my_param_hist[m].count(itr->bin());
// 	  double my_var = data->my_param_hist[m].var(itr->bin());
// 	  double mphi = atan2(my,mx);
// 	  double mr = sqrt(std::pow(mx,2)+std::pow(my,2));
// 	  double mr_var = 
// 	    std::pow(mx/mr,2)*mx_var + std::pow(my/mr,2)*my_var +
// 	    2*mcov[m][itr->bin()]*mx*my/std::pow(mr,2);

// 	  double mphi_var =
// 	    std::pow(mx/(mr*mr),2)*mx_var + 
// 	    std::pow(my/(mr*mr),2)*my_var -
// 	    2*mx*my/std::pow(mr,4)*mcov[m][itr->bin()];

// 	  mphi *= 180/M_PI;
// 	  mphi_var *= std::pow(180/M_PI,2);

// 	  data->mr_param_hist[m].setBin(itr->bin(),mr,mr_var);
// 	  data->mphi_param_hist[m].setBin(itr->bin(),mphi,mphi_var);

// 	  if(m == 0)
// 	    {
// 	      m1r_param.push_back(mr);
// 	      m1r_param_err.push_back(sqrt(mr_var));
// 	      m1phi_param.push_back(mphi);
// 	      m1phi_param_err.push_back(sqrt(mphi_var));
// 	    }
// 	  else if(m == 1)
// 	    {
// 	      m2r_param.push_back(mr);
// 	      m2r_param_err.push_back(sqrt(mr_var));
// 	      m2phi_param.push_back(mphi);
// 	      m2phi_param_err.push_back(sqrt(mphi_var));
// 	    }
// 	}
//     }

}

// ============================================================================
// Evaluation Methods
// ============================================================================
double VSAcceptanceModelBessel2::val(const VSAAlgebra::Vec2D& xy) const
{ 
  const unsigned nobs = m_obs_xy.size();

  double val = 0;
  for(unsigned iobs = 0; iobs < nobs; iobs++)
    {
      VSModelCoord c(xy,iobs);
      val += VSAcceptanceModelBessel2::val(c);
    }
  return val;
}

void VSAcceptanceModelBessel2::val(const VSAAlgebra::Vec2D& xy, 
				   const VSAAlgebra::VecND& a,
				   const VSAAlgebra::MatrixND& cov, 
				   double& val, double& err) const
{
  const unsigned nobs = m_obs_xy.size();

  val = 0;
  err = 0;

  VSAAlgebra::VecND dyda(nparm());
  for(unsigned iobs = 0; iobs < nobs; iobs++)
    {
      VSModelCoord c(xy,iobs);
      val += VSAcceptanceModelBessel2::val(c,a);

      VSAAlgebra::VecND tmp;
      VSAcceptanceModelBessel2::dyda(c,a,tmp);
      dyda += tmp;
    }

  err = sqrt(dyda*(cov*dyda));
}

double VSAcceptanceModelBessel2::val(const VSModelCoord& x) const
{
  VSAAlgebra::Vec2D xy(x.x(),x.y());
  xy -= m_obs_xy[x.iobs()];

  if(xy.norm() > m_offset_max) return 0;

  double val = (m_fn.val(VSACoord::Coord2D(xy)))*param(x.iobs()+m_fn.nparm());
  return val*m_domega;
}
    
double VSAcceptanceModelBessel2::val(const VSModelCoord& x, 
				     const VSAAlgebra::VecND& a) const
{
  VSAAlgebra::Vec2D xy(x.x(),x.y());
  xy -= m_obs_xy[x.iobs()];

  if(xy.norm() > m_offset_max) return 0;

  VSAAlgebra::VecND fna = a.subVector(0,m_fn.nparm());

  double val = (m_fn.val(VSACoord::Coord2D(xy),fna))*a(x.iobs()+m_fn.nparm());
  return val*m_domega;
}

void VSAcceptanceModelBessel2::dyda(const VSModelCoord& x, 
				    VSAAlgebra::VecND& dyda) const
{
  dyda.resize(nparm());
  dyda.clear();

  VSAAlgebra::Vec2D xy(x.x(),x.y());
  xy -= m_obs_xy[x.iobs()];

  if(xy.norm() > m_offset_max) return;

  VSAAlgebra::VecND tmp;
  m_fn.dyda(VSACoord::Coord2D(xy),tmp);

  tmp *= param(x.iobs()+m_fn.nparm());

  for(unsigned ip = 0; ip < tmp.ndim(); ip++) dyda(ip) += tmp(ip);

  dyda(m_fn.nparm()+x.iobs()) = m_fn.val(VSACoord::Coord2D(xy));
  dyda *= m_domega;
}

void VSAcceptanceModelBessel2::dyda(const VSModelCoord& x, 
				    const VSAAlgebra::VecND& a,
				    VSAAlgebra::VecND& dyda) const
{
  dyda.resize(nparm());
  dyda.clear();

  VSAAlgebra::Vec2D xy(x.x(),x.y());
  xy -= m_obs_xy[x.iobs()];

  if(xy.norm() > m_offset_max) return;

  VSAAlgebra::VecND tmp;
  m_fn.dyda(VSACoord::Coord2D(xy),a.subVector(0,m_fn.nparm()),tmp);

  tmp *= a(x.iobs()+m_fn.nparm());

  for(unsigned ip = 0; ip < tmp.ndim(); ip++) dyda(ip) += tmp(ip);

  dyda(m_fn.nparm()+x.iobs()) = 
    m_fn.val(VSACoord::Coord2D(xy),a.subVector(0,m_fn.nparm()));
  dyda *= m_domega;
}

double VSAcceptanceModelBessel2::getFoVCounts(const VSAAlgebra::Vec2D& xy) 
  const 
{ 
  if(xy.norm() > m_offset_max) return 0;
  const unsigned nobs = m_obs_xy.size();
  VSAAlgebra::VecND a = param().subVector(0,m_fn.nparm());

  double val = 0;
  for(unsigned iobs = 0; iobs < nobs; iobs++)
    val += m_fn.val(VSACoord::Coord2D(xy),a)*param(iobs+m_fn.nparm());

  return val;
}

double VSAcceptanceModelBessel2::getFoVAcceptance(const VSAAlgebra::Vec2D& xy) 
  const 
{ 
  if(xy.norm() > m_offset_max) return 0;

  const VSAAlgebra::VecND a = param().subVector(0,m_fn.nparm());
  return m_fn.val(VSACoord::Coord2D(xy),a);
}

double VSAcceptanceModelBessel2::getFoVAcceptance(double r2) 
  const 
{ 
  if(sqrt(r2) > m_offset_max) return 0;

  VSACoord::Polar<VSACoord::Coord2D> rphi1(sqrt(r2),0);
  VSACoord::Polar<VSACoord::Coord2D> rphi2(sqrt(r2),2*M_PI);

  VSAFunction::FourierBesselSeries2 fn(m_fn);
  const VSAAlgebra::VecND a = param().subVector(0,m_fn.nparm());
  fn.setParam(a);

  return 0.5*VSAMath::Quadrature::integrate1D(fn,rphi1,rphi2,1,1E-4);
}

VSAcceptanceData* 
VSAcceptanceModelBessel2::fit(const VSAMath::Data<VSModelCoord>& data)
{
  // --------------------------------------------------------------------------
  // Generate first pass estimates of fit parameters
  VSAAlgebra::VecND fit_param(nparm(),0.);
  calc(data,fit_param);  
  
  VSAcceptanceDataBessel2* d = new VSAcceptanceDataBessel2(m_offset_max);

  d->param      = VSAAlgebra::VecND(nparm(),0.);
  d->param_err  = VSAAlgebra::VecND(nparm(),0.);
  d->param_cov  = VSAAlgebra::MatrixND(nparm(),nparm(),0.);

  // --------------------------------------------------------------------------
  // Instantiate the non-linear fitter
  typedef VSAMath::NLFitter< LnLFn > LnLFitter;
  LnLFitter* fitter = VSAMath::NLFitterFactory::createLnL(this,data);

  fitter->initialize(fit_param);
  fitter->fit();
  
  fit_param = fitter->param();
  VSAAlgebra::MatrixND  fit_cov = fitter->cov();

  d->param     = fit_param;
  d->param_err = fitter->err();
  d->param_cov = fit_cov;
  d->chi2      = fitter->chi2();

  setParam(fit_param);
  fillData(d);

  d->acceptance_param.resize(m_fn.nparm());
  d->acceptance_param_err.resize(m_fn.nparm());
  d->acceptance_param_cov.resize(m_fn.nparm(),m_fn.nparm());

  for(unsigned ip1 = 0; ip1 < m_fn.nparm(); ip1++)
    {
      d->acceptance_param(ip1) = d->param(ip1);
      d->acceptance_param_err(ip1) = d->param_err(ip1);

      for(unsigned ip2 = 0; ip2 < m_fn.nparm(); ip2++)
	d->acceptance_param_cov(ip1,ip2) = d->param_cov(ip1,ip2);
    }

  delete fitter;

  return d;
}

void VSAcceptanceModelBessel2::integrate(const VSAAlgebra::VecND& a,
					 const VSAAlgebra::MatrixND& cov,
					 double& integral,
					 double& integral_err) const
{
  const double rmax = 2.5;
  const double xlo = -rmax;
  const double xhi = rmax;
  const double ylo = -rmax;
  const double yhi = rmax;
  
  integral = 0;
  const unsigned nobs = m_obs_xy.size();
  VSAAlgebra::VecND b(a.ndim());

  for(unsigned iobs = 0; iobs < nobs; iobs++)
    {
      for(double x = xlo; x < xhi; x+= m_bin_size)
	for(double y = ylo; y < yhi; y+= m_bin_size)
	  {
	    VSModelCoord c(x,y,iobs);
	    integral += val(c);
	    
	    VSAAlgebra::VecND dyda;
	    VSAcceptanceModelBessel2::dyda(c,a,dyda);

	    b += dyda;	  
	  }
    }

  VSAAlgebra::VecND b2 = cov*b;
  integral_err = sqrt(b*b2);
}

double VSAcceptanceModelBessel2::integrate(const VSAAlgebra::Vec2D& xy, 
					   double R1,double R2) const
{
  const unsigned nobs = m_obs_xy.size();

  double z = 0;
  for(unsigned iobs = 0; iobs < nobs; iobs++)
    {
      VSACoord::Polar<VSModelCoord> rphi1(xy,R1,0);
      VSACoord::Polar<VSModelCoord> rphi2(xy,R2,2*M_PI);

      rphi1.setObs(iobs);
      rphi2.setObs(iobs);

      try
	{
	  z += VSAMath::Quadrature::integrate(*this,rphi1,rphi2,1E-3);
	}
      catch(const std::exception& e)
	{
	  std::cerr 
	    << std::string(__PRETTY_FUNCTION__) 
	    << ": Caught exception." << std::endl
	    << "X = " << xy.x() << " Y = " << xy.y() << std::endl
	    << "R1 = " << R1 << " R2 = " << R2 << " IOBS = " << iobs
	    << std::endl
	    << e.what() << std::endl;
	}
    }

  return z/m_domega;
}

double VSAcceptanceModelBessel2::integrateFoV(double R1,double R2) const
{
  double z = 0;

  const unsigned nobs = m_obs_xy.size();
  for(unsigned iobs = 0; iobs < nobs; iobs++)
    {
      const VSAAlgebra::Vec2D& xy = m_obs_xy[iobs];
      VSACoord::Polar<VSModelCoord> rphi1(xy,R1,0);
      VSACoord::Polar<VSModelCoord> rphi2(xy,R2,2*M_PI);

      rphi1.setObs(iobs);
      rphi2.setObs(iobs);

      try
	{
	  z += VSAMath::Quadrature::integrate(*this,rphi1,rphi2,1E-3);
	}
      catch(const std::exception& e)
	{
	  std::cerr 
	    << std::string(__PRETTY_FUNCTION__) 
	    << ": Caught exception." << std::endl
	    << "X = " << xy.x() << " Y = " << xy.y() << std::endl
	    << "R1 = " << R1 << " R2 = " << R2 << " IOBS = " << iobs
	    << std::endl
	    << e.what() << std::endl;
	}
    }

  return z/m_domega;
}

void VSAcceptanceModelBessel2::calc(const VSAMath::Data<VSModelCoord>& data,
				    VSAAlgebra::VecND& p)
{
  std::cout << "VSAcceptanceModelBessel2::calc()" << std::endl;
  vsassert(data.size());

  p.resize(nparm());

  unsigned nr = 0;
  for(FourierBesselSeries2::const_iterator itr = m_fn.begin(); itr !=
	m_fn.end(); ++itr)
    if(itr->m==0) nr++;

  vsassert(nr > 0);

  const unsigned naz = 1;
  const unsigned nobs = m_obs_xy.size();

  VSAFunction::FourierBesselSeries bessel_fn(naz,nr,m_offset_max);
  const unsigned np = bessel_fn.nparm();

  std::vector< VSAAlgebra::MatrixND > t(nobs,VSAAlgebra::MatrixND(np,np));
  std::vector< VSAAlgebra::MatrixND > nt(nobs,VSAAlgebra::MatrixND(np,np));
  std::vector< VSAAlgebra::VecND > s(nobs,  VSAAlgebra::VecND(np));
  std::vector< VSAAlgebra::VecND > ns(nobs,  VSAAlgebra::VecND(np));

  std::vector< double > nsum(nobs,0);
  std::vector< double > nbin(nobs,0);

  for(unsigned idata = 0; idata < data.size(); idata++)
    {
      unsigned iobs = data[idata].x.iobs();

      VSAAlgebra::VecND v;
      VSAAlgebra::Vec2D xy(data[idata].x.x(),data[idata].x.y());
      xy -= m_obs_xy[iobs];

      bessel_fn(VSACoord::Coord2D(xy),v);
      
      double y = data[idata].y;

      nbin[iobs]++;
      nsum[iobs] += y;
      s[iobs] += v;
      ns[iobs] += y*v;

      for(unsigned ix = 0; ix < np; ix++)
	for(unsigned iy = 0; iy < np; iy++)
	  {
	    t[iobs](ix,iy) += v[ix]*v[iy];
	    nt[iobs](ix,iy) += v[ix]*v[iy]*y;
	  }
    }

  VSAAlgebra::MatrixND q(np,np);
  VSAAlgebra::VecND b(np);
  VSAAlgebra::VecND c(np);

  for(unsigned iobs = 0; iobs < nobs; iobs++)
    {
      if(nbin[iobs] == 0) continue;
      b += ns[iobs] - (nsum[iobs]/nbin[iobs])*s[iobs];
      for(unsigned ix = 0; ix < np; ix++)
	for(unsigned iy = 0; iy < np; iy++)
	  q(ix,iy)+=nsum[iobs]/nbin[iobs]*t[iobs](ix,iy)+nt[iobs](ix,iy)-
	    2*nsum[iobs]/(nbin[iobs]*nbin[iobs])*s[iobs](ix)*s[iobs](iy);
    }

  VSAAlgebra::SVD svd(q);
  svd.solve(b,c);

  double epsilon = m_domega/(M_PI*std::pow(m_offset_max,2));
  double c00 = 1;
  for(unsigned ip = 0; ip < c.ndim(); ip++)
    {
      c00 += c(ip)*4/VSAMath::numeric_constants<double>::bessel_zero(0,ip) +
	c(ip)*c(ip);
    }

  c00 = 1/sqrt(c00);

  for(unsigned ip = 0; ip < c.ndim(); ip++)
    {
      std::cout << "INPUT " 
		<< std::setw(10) << ip 
		<< std::setw(15) << c(ip)*c00 << std::endl;

      p(ip) = c(ip)*c00;
    }


  std::cout << "c00: " << c00 << std::endl;

  for(unsigned iobs = 0; iobs < nobs; iobs++)
    {
      double s2 = nbin[iobs] + 2*c*s[iobs] + c*(t[iobs]*c);
      double bnorm = nsum[iobs]/(s2*epsilon*c00*c00);

      if(std::isfinite(bnorm) && bnorm > 0) p(m_fn.nparm()+iobs) = bnorm;
      else p(m_fn.nparm()+iobs) = 1.0;

      std::cout << "B" << iobs << ": " 
		<< std::setw(15) << p(m_fn.nparm()+iobs) << std::endl;
    }
}
