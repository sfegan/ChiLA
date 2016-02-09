#include <VSLTLibraryVisitor.hpp>
#include <SphericalCoords.h>
#include <VSANonlinearFitting.hpp>
#include <VSScaledParameterLibrary.hpp>
#include <VSNSpaceOctaveH5IO.hpp>
#include <VSALinearLeastSquares.hpp>
#include <VSAFunction.hpp>
#include <VSALocalRegression.hpp>

using namespace VERITAS;
using namespace SEphem;

// ============================================================================
// VSEnergyLTLibraryVisitor
// ============================================================================

VSEnergyLTLibraryVisitor::Fn::Fn():
  VSAFunction::ParamFn<VSAAlgebra::Vec2D>(14)
{

}

void VSEnergyLTLibraryVisitor::Fn::operator() 
  (const VSAAlgebra::Vec2D& x, VSAAlgebra::VecND& v) const
{
  v.resize(nparm());

  v[0] = 1;
  v[1] = x.x();
  v[2] = std::pow(x.x(),2);
  v[3] = std::pow(x.x(),3);
  v[4] = std::pow(x.x(),4);
  v[5] = x.y();
  v[6] = std::pow(x.y(),2);
  v[7] = std::pow(x.y(),3);
  v[8] = std::pow(x.y(),4);
  v[9] = x.x()*x.y();
  v[10] = std::pow(x.x(),2)*x.y();
  v[11] = std::pow(x.x(),3)*x.y();
  v[12] = x.x()*std::pow(x.y(),2);
  v[13] = x.x()*std::pow(x.y(),3);
}

double VSEnergyLTLibraryVisitor::Fn::val(const VSAAlgebra::Vec2D& x, 
					 const VSAAlgebra::VecND& a) const
{
  VSAAlgebra::VecND v;
  (*this)(x,v);
  return a * v;
}

// ----------------------------------------------------------------------------
// VSEnergyLTLibraryVisitor::Spaces
// ----------------------------------------------------------------------------
VSEnergyLTLibraryVisitor::Spaces::
Spaces(const VSNSpace::Space& energy_def, const VSNSpace::Space& size_def, 
       unsigned _min_count, const std::string& _comment_base):
  min_count(_min_count), comment_base(_comment_base),
  energy_space(energy_def), 
  n_energy(energy_def), 
  log10_energy_mean(energy_def),  log10_energy_rms(energy_def),
  log10_energy_med(energy_def),   log10_energy_i68(energy_def),
  n_size(size_def), 
  log10_size_mean(size_def),  log10_size_rms(size_def),
  log10_size_med(size_def),   log10_size_i68(size_def)
{ 

}

void VSEnergyLTLibraryVisitor::Spaces::
accumulate(const VSNSpace::Point& p_energy, double log10_energy,
	   const VSNSpace::Point& p_size, double log10_size, double weight)
{
  try
    {
      n_energy[p_energy]++;
      log10_energy_mean[p_energy].accumulate(log10_energy,weight);
      log10_energy_rms[p_energy].accumulate(log10_energy,weight);
      log10_energy_med[p_energy].accumulate(log10_energy,weight);
      log10_energy_i68[p_energy].accumulate(log10_energy,weight);

      n_size[p_size]++;
      log10_size_mean[p_size].accumulate(log10_size,weight);
      log10_size_rms[p_size].accumulate(log10_size,weight);
      log10_size_med[p_size].accumulate(log10_size,weight);
      log10_size_i68[p_size].accumulate(log10_size,weight);
    }
  catch(std::out_of_range&)
    {
      // nothing to see here
    }
}

void VSEnergyLTLibraryVisitor::Spaces::
calculateMomentTables(const VSNSpace& n,
		      VSNSpace& exp, VSNSpace& exp_var, 
		      VSNSpace& rms, VSNSpace& rms_var)
{
  VSNSpace::Volume mask = n>min_count;	
  exp.clearOutsideVolume(mask);
  rms.clearOutsideVolume(mask);

  exp_var = rms;
  rms_var = rms;

  exp_var *= exp_var;
  exp_var /= n;

  exp_var.clearOutsideVolume(mask);

  rms_var *= rms;
  rms_var *= 0.5;
  rms_var /= n;

  rms_var.clearOutsideVolume(mask);
}

void VSEnergyLTLibraryVisitor::Spaces::
calculateMedianTables(const VSNSpace& n,
		      VSNSpace& exp, VSNSpace& exp_var, 
		      VSNSpace& rms, VSNSpace& rms_var)
{
  VSNSpace::Volume mask = n>min_count;	
  exp.clearOutsideVolume(mask);
  rms.clearOutsideVolume(mask);

  exp_var = rms;
  rms_var = rms;

  exp_var *= exp_var;
  exp_var *= 1.57;
  exp_var /= n;

  exp_var.clearOutsideVolume(mask);

  rms_var *= rms;
  rms_var *= 0.5;
  rms_var /= n;

  rms_var.clearOutsideVolume(mask);
}

void VSEnergyLTLibraryVisitor::Spaces::
extrapolate(VSNSpace& exp, VSNSpace& exp_var)
{
 
  VSNSpace::Cell c(exp.space().cell());
  const unsigned ncell = exp.space().size();
  for(unsigned icell = 0; icell < ncell; icell++)
    {   
      exp.space().cellOfIndexUnchecked(icell,c);      
      if(c.i[1] != 0) continue;

      VSAMath::Data<double> xys;
      for(c.i[1]=0;c.i[1]<exp.axis(1).nbin;c.i[1]++)
	{
	  double logS = exp.axis(1).midCoordUnchecked(c.i[1]);

	  if(exp_var[c] > 0)
	    {
	      VSAMath::DataPoint<double> dp(logS,exp[c],sqrt(exp_var[c]));
	      xys.insert(dp);
	    }
	}

      if(xys.size() <= 1) continue;
      const unsigned np = std::min((int)3,(int)(xys.size()-1));

      VSAAlgebra::VecND param;
      VSAMath::PolyFit::fit(np, xys, param);
      
      for(c.i[1]=0;c.i[1]<exp.axis(1).nbin;c.i[1]++)
	{
	  double logS = exp.space().axes[1].midCoordUnchecked(c.i[1]);
	  if(exp[c] > 0) exp[c] = VSAMath::PolyFit::val(param,logS);
	}
    }
}

void VSEnergyLTLibraryVisitor::Spaces::
extrapolate2(VSNSpace& exp, VSNSpace& exp_var)
{
  VSNSpace::Cell c(exp.space().cell());
  
  VSAMath::Data<VSAAlgebra::Vec2D> exp_xys;
  // --------------------------------------------------------------------------
  // Loop over Impact Distance
  // --------------------------------------------------------------------------
  const unsigned nR = exp.space().axes[0].nbin;
  for(unsigned ir=0; ir<nR; ir++)
    {      
      c.i[0] = ir;
      double R = exp.space().axes[0].midCoordUnchecked(c.i[0]);

      const unsigned nS = exp.space().axes[1].nbin;      
      for(c.i[1]=0;c.i[1]<nS;c.i[1]++)
	{
	  double logS = exp.space().axes[1].midCoordUnchecked(c.i[1]);

	  if(exp_var[c] > 0)
	    {
	      VSAMath::DataPoint<VSAAlgebra::Vec2D> 
		dp(VSAAlgebra::Vec2D(logS,R),exp[c],sqrt(exp_var[c]));
	      exp_xys.insert(dp);
	    }
	}
    }

  //  VSAFunction::Poly2D<VSAAlgebra::Vec2D> fn(4);
  Fn fn;

  //  VSAMath::Fitsvd<VSAFunction::Poly2D<VSAAlgebra::Vec2D>,VSAAlgebra::Vec2D> 
  VSAMath::Fitsvd<Fn,VSAAlgebra::Vec2D> 
    svdfit(exp_xys,fn);
  svdfit.fit();

  VSAAlgebra::VecND param = svdfit.param();

  for(unsigned ir=0; ir<nR; ir++)
    {      
      c.i[0] = ir;
      double R = exp.space().axes[0].midCoordUnchecked(c.i[0]);

      const unsigned nS = exp.space().axes[1].nbin;      
      for(c.i[1]=0;c.i[1]<nS;c.i[1]++)
	{
	  double logS = exp.space().axes[1].midCoordUnchecked(c.i[1]);
	  VSAAlgebra::Vec2D x(logS,R);
	  exp[c] = fn.val(x,param);
	}
    }
}

void VSEnergyLTLibraryVisitor::Spaces::smooth(VSNSpace& exp, VSNSpace& exp_var,
					      double dx1, double dx2,
					      unsigned norder)
{
  VSNSpace::Cell c(exp.space().cell());

  const unsigned ncell = exp.space().size();
  for(unsigned icell = 0; icell < ncell; icell++)
    {
      exp.space().cellOfIndexUnchecked(icell,c);      
      if(c.i[0] != 0 || c.i[1] != 0) continue;

      VSAMath::Data<VSAAlgebra::Vec2D> xys;
      for(c.i[0]=0; c.i[0] < exp.axis(0).nbin; c.i[0]++)
	{      
	  double x0 = exp.axis(0).midCoordUnchecked(c.i[0]);	  
	  for(c.i[1]=0;c.i[1]<exp.axis(1).nbin;c.i[1]++)
	    {
	      double x1 = exp.axis(1).midCoordUnchecked(c.i[1]);	      
	      if(exp_var[c] > 0)
		{
		  VSAMath::DataPoint<VSAAlgebra::Vec2D> 
		    dp(VSAAlgebra::Vec2D(x0,x1),exp[c],sqrt(exp_var[c]));
		  xys.insert(dp);
		}
	    }
	}
      
      VSAMath::LocalRegression2D<VSAAlgebra::Vec2D> lr2d(xys,dx1,dx2,norder);
      
      for(c.i[0]=0; c.i[0] < exp.axis(0).nbin; c.i[0]++)
	{      
	  double x0 = exp.axis(0).midCoordUnchecked(c.i[0]);
	  for(c.i[1]=0;c.i[1]<exp.axis(1).nbin;c.i[1]++)
	    {
	      double x1 = exp.axis(1).midCoordUnchecked(c.i[1]);
	      VSAAlgebra::Vec2D x(x0,x1);
	      
	      double val, var;
	      lr2d.val(x,val,var);
	      
	      exp[c] = val;
	      exp_var[c] = var;
	    }
	}
    }
}

void VSEnergyLTLibraryVisitor::Spaces::
calculateMomentEnergyTables(VSNSpace& exp, VSNSpace& exp_var, 
			    VSNSpace& rms, VSNSpace& rms_var)
{
  VSNSpace::Volume mask = n_energy>min_count;	
  log10_energy_mean.nspace(exp);
  exp.clearOutsideVolume(mask);
  log10_energy_rms.nspace(rms);
  rms.clearOutsideVolume(mask);
  exp.setComment(comment_base+" log10 energy mean");
  rms.setComment(comment_base+" log10 energy standard deviation");
}
      
void VSEnergyLTLibraryVisitor::Spaces::
calculateMedianEnergyTables(VSNSpace& exp, VSNSpace& exp_var, 
			    VSNSpace& rms, VSNSpace& rms_var)
{
  log10_energy_med.nspace(exp);
  log10_energy_i68.nspace(rms);

  calculateMedianTables(n_energy,exp,exp_var,rms,rms_var);
  exp.setComment(comment_base+" log10 energy median");
  rms.setComment(comment_base+" log10 energy 68% containment interval");
}
      
void VSEnergyLTLibraryVisitor::Spaces::
calculateMomentSizeTables(VSNSpace& exp, VSNSpace& exp_var, 
			  VSNSpace& rms, VSNSpace& rms_var)
{
  VSNSpace::Volume mask = n_size>min_count;	
  log10_size_mean.nspace(exp);
  exp.clearOutsideVolume(mask);
  log10_size_rms.nspace(rms);
  rms.clearOutsideVolume(mask);

  calculateMomentTables(n_size,exp,exp_var,rms,rms_var);
  exp.setComment(comment_base+" log10 size mean");
  rms.setComment(comment_base+" log10 size standard deviation");
}
      
void VSEnergyLTLibraryVisitor::Spaces::
calculateMedianSizeTables(VSNSpace& exp, VSNSpace& exp_var, 
			  VSNSpace& rms, VSNSpace& rms_var)
{
  log10_size_med.nspace(exp);
  log10_size_i68.nspace(rms);

  calculateMedianTables(n_size,exp,exp_var,rms,rms_var);
  exp.setComment(comment_base+" log10 size median");
  rms.setComment(comment_base+" log10 size 68% containment interval");
}

void VSEnergyLTLibraryVisitor::Spaces::fit(const VSAMath::Data<double>& xys,
					   double& fit_chi2,
					   double& fit_rchi2,
					   VSAAlgebra::VecND& fit_param,
					   VSAAlgebra::MatrixND& fit_cov,
					   double poly_rchi2_max_change,
					   unsigned poly_max_order)
{
  const unsigned ndata = xys.size();

  unsigned nparam = ndata;
  if(nparam>poly_max_order)nparam=poly_max_order;

  fit_chi2 = std::numeric_limits<double>::infinity();
  fit_rchi2 = fit_chi2;
  double fit_min_rchi2 = fit_chi2;
  
  // std::cout << "===========================================" << std::endl;
  // std::cout << "SIZE " << xys.size() << std::endl;
  // for(unsigned i = 0; i < xys.size(); i++)
  //   std::cout << std::setw(5) << i
  // 	      << std::setw(15) << xys[i].x
  // 	      << std::setw(15) << xys[i].y
  // 	      << std::setw(15) << xys[i].sigma
  // 	      << std::endl;

  // ----------------------------------------------------------------------
  // Continue increasing the number of parameters in the fit
  // until the difference in chi2 reaches the threshold
  // ----------------------------------------------------------------------
  for(unsigned jparam=0;jparam<nparam-1;jparam++)
    {
      unsigned iparam = nparam-jparam-1;
      unsigned ndf = ndata-iparam-1;

      // std::cout << "-----------------------------------" << std::endl;
      //      std::cout << "NPARAM " << iparam << std::endl;

      VSAAlgebra::VecND test_param;
      VSAAlgebra::MatrixND test_cov;

      
      VSAMath::Fitcsvd<VSAFunction::Poly,double> 
	fitcsvd(xys,VSAFunction::Poly(iparam));

      fitcsvd.addMonotonicConstraints(-1.5,0.0625,57);

      fitcsvd.fit();

      // double chi2 = VSAMath::PolyFit::fit(iparam, xys, 
      // 					  test_param, &test_cov);

      
      // for(unsigned ip = 0; ip < test_param.ndim(); ip++)
      // 	{
      // 	  std::cout << ip << " " 
      // 		    << std::setw(20) << test_param(ip) 
      // 		    << std::setw(20) << sqrt(test_cov(ip,ip)) 
      // 		    << std::setw(20) << fitcsvd.param()[ip] 
      // 		    << std::setw(20) << sqrt(fitcsvd.cov()(ip,ip))
      // 		    << std::endl;
      // 	}


      //  std::cout << fitcsvd.param() << std::endl;
      //  std::cout << "CHI2 " << fitcsvd.chi2() << " " << chi2 << std::endl;

      test_param = fitcsvd.param();
      test_cov = fitcsvd.cov();
      double chi2 = fitcsvd.chi2();
      

      double rchi2 = 0;
      if(ndf > 0) rchi2 = chi2/double(ndf);
      
      // std::cout << iparam << " " << chi2 << " " << rchi2 << " "
      // 		<< fit_min_rchi2 << " " << poly_rchi2_max_change*fit_min_rchi2
      // 		<< std::endl;

      if(rchi2 < fit_min_rchi2)fit_min_rchi2 = rchi2;
      if(poly_rchi2_max_change*fit_min_rchi2 < rchi2 || iparam == 0)break;
      
      fit_param  = test_param;
      fit_cov  = test_cov;
      fit_chi2   = chi2;
      fit_rchi2  = rchi2;
    }
}

void VSEnergyLTLibraryVisitor::Spaces::
size2Energy(const VSNSpace& s_exp, const VSNSpace& s_exp_var, 
	    const VSNSpace& s_rms, const VSNSpace& s_rms_var, 
	    const VSNSpace& s_n, 
	    VSNSpace& e_exp, VSNSpace& e_exp_var, 
	    VSNSpace& e_rms, VSNSpace& e_rms_var, 
	    std::vector<SizeFitData>& size_fit_data,
	    double poly_rchi2_max_change,
	    unsigned poly_max_order)
{
  e_exp = VSNSpace(energy_space, 
		   std::string("log10 energy from: ")+s_exp.comment());
  e_exp_var = VSNSpace(energy_space, 
		       std::string("log10 energy from: ")+s_exp.comment());
  e_rms = VSNSpace(energy_space, 
		   std::string("log10 energy rms from: ")+s_exp.comment());
  e_rms_var = VSNSpace(energy_space, 
		       std::string("log10 energy rms from: ")+s_exp.comment());

  const VSNSpace::Space& s_space(s_exp.space());
  const unsigned nE = s_space.axes[1].nbin;
  const unsigned ncell = s_space.size();
  
  VSNSpace::Cell c(s_space.cell());
  
  VSNSpace s_rms_smooth = s_rms;
  VSNSpace s_rms_var_smooth = s_rms_var;
  VSNSpace s_log_rms_smooth = s_rms;
  VSNSpace s_log_rms_var_smooth = s_rms_var;

  for(unsigned icell = 0; icell < ncell; icell++)
    {
      s_space.cellOfIndexUnchecked(icell, c);
      if(s_rms_var[c] > 0)
	{
	  s_log_rms_smooth[c] = std::log(s_rms[c]);
	  s_log_rms_var_smooth[c] *= std::pow(s_rms[c],-2);
	}
    }

  smooth(s_log_rms_smooth,s_log_rms_var_smooth,
	 s_log_rms_smooth.axis(0).bin_size*3,
	 s_log_rms_smooth.axis(1).bin_size*4,1);

  for(unsigned icell = 0; icell < ncell; icell++)
    {
      s_space.cellOfIndexUnchecked(icell, c);
      s_rms_smooth[c] = std::exp(s_log_rms_smooth[c]);
      s_rms_var_smooth[c] = 
	std::pow(s_rms_smooth[c],2)*s_log_rms_var_smooth[c];
    }

  for(unsigned icell = 0; icell < ncell; icell++)
    {
      s_space.cellOfIndexUnchecked(icell, c);
      if(c.i[1]!=0)continue;      

      // double R = s_space.axes[0].midCoordUnchecked(c.i[0]);
      // std::cout << icell << " " << R << std::endl;

      SizeFitData sfd;
      sfd.R = s_space.axes[0].midCoordUnchecked(c.i[0]);
      if(c.ndim == 3) sfd.disp = s_space.axes[2].midCoordUnchecked(c.i[2]);

      sfd.size_exp_hist = 
	VSLimitedErrorsHist<double,double>(s_space.axes[1].bin_size,
					   s_space.axes[1].lo_bound,
					   s_space.axes[1].hi_bound);
      
      sfd.size_rms_hist = sfd.size_exp_hist;
      sfd.size_fit_exp_hist = sfd.size_exp_hist;
      sfd.size_fit_rms_hist = sfd.size_exp_hist;

      // -----------------------------------------------------------------------
      // In each distance bin fit a polynomial to Size vs. Energy
      // -----------------------------------------------------------------------
      VSAMath::Data<double> exp_xys;
      VSAMath::Data<double> rms_xys;
      double xmin = std::numeric_limits<double>::infinity();
      double xmax = -xmin;
      for(c.i[1]=0;c.i[1]<nE;c.i[1]++)
	{
	  VSNSpace::Point p(2);
	  s_space.midPointOfCellUnchecked(c, p);

	  sfd.size_exp_hist.accumulate(p.x[1],s_exp[c],s_exp_var[c]);
	  sfd.size_rms_hist.accumulate(p.x[1],s_rms[c],s_rms_var[c]);

	  if(s_rms[c] > 0)
	    {
	      VSAMath::DataPoint<double> exp_datum;
	      VSAMath::DataPoint<double> rms_datum;
	      if(p.x[1] > xmax)xmax=p.x[1];
	      if(p.x[1] < xmin)xmin=p.x[1];
	      exp_datum.x = p.x[1];
	      exp_datum.y = s_exp[c];
	      exp_datum.sigma = s_rms[c]/sqrt(s_n[c]);
	      exp_xys.insert(exp_datum);

	      rms_datum.x = p.x[1];
	      rms_datum.y = s_rms[c];
	      rms_datum.sigma = sqrt(s_rms_var[c]);
	      rms_xys.insert(rms_datum);
	    }
	}

      if(exp_xys.size() <= 1)continue;

      VSAAlgebra::VecND exp_fit_param;
      VSAAlgebra::MatrixND exp_fit_cov;
      double exp_fit_chi2;
      double exp_fit_rchi2;
      fit(exp_xys,exp_fit_chi2,exp_fit_rchi2,exp_fit_param,exp_fit_cov,
	  poly_rchi2_max_change,poly_max_order);

      sfd.exp_fit_param = exp_fit_param;
      sfd.exp_fit_nparam = exp_fit_param.ndim();
      sfd.exp_fit_chi2 = exp_fit_chi2;
      sfd.exp_fit_rchi2 = exp_fit_rchi2;

      for(c.i[1]=0;c.i[1]<nE;c.i[1]++)
	{
	  VSNSpace::Point p(2);
	  s_space.midPointOfCellUnchecked(c, p);

	  sfd.size_fit_exp_hist.
	    accumulate(p.x[1],VSAMath::PolyFit::val(exp_fit_param,p.x[1]),0);

	  sfd.size_fit_rms_hist.
	    accumulate(p.x[1],s_rms_smooth[p],s_rms_var_smooth[p]);

	  // sfd.size_fit_rms_hist.
	  //   accumulate(p.x[1],VSAMath::PolyFit::val(rms_fit_param,p.x[1]),0);
	}

      // ----------------------------------------------------------------------
      // Find first point on polynomial fit to Size vs. Energy with
      // positive slope
      // ----------------------------------------------------------------------
      VSAAlgebra::VecND dfit_param;
      VSAMath::PolyFit::differentiate(exp_fit_param,dfit_param);

      double x;
      double dx = 0.001;

      for(x=xmin;x<=xmax;x+=dx)
	{
	  double dydx = VSAMath::PolyFit::val(dfit_param,x);
	  if(dydx>0)break;
	}

      if(x>xmax)continue; // slope is ALWAYS negative


      double y = VSAMath::PolyFit::val(exp_fit_param,x);
      VSAAlgebra::VecND dlogS_da1;
      VSAMath::PolyFit::dyda(exp_fit_param,x,dlogS_da1);

      sfd.exp_hist = 
	VSLimitedErrorsHist<double,double>(energy_space.axes[1].bin_size,
					   energy_space.axes[1].lo_bound,
					   energy_space.axes[1].hi_bound);
      
      sfd.rms_hist = 
	VSLimitedErrorsHist<double,double>(energy_space.axes[1].bin_size,
					   energy_space.axes[1].lo_bound,
					   energy_space.axes[1].hi_bound);

      // -----------------------------------------------------------------------
      // Loop over Image Size
      // -----------------------------------------------------------------------
      unsigned nS = energy_space.axes[1].nbin;
      for(c.i[1]=0;c.i[1]<nS;c.i[1]++)
	{
	  double logS = energy_space.axes[1].midCoordUnchecked(c.i[1]);
	    
	  if(y>logS)continue;
	  double yold = y;
	  while((x<=xmax)&&(y<=logS))
	    x += dx, yold = y, y = VSAMath::PolyFit::val(exp_fit_param,x);
	  if(x>xmax)break;
	    
	  double logE = (logS-y)/(y-yold)*dx + x;

	  VSAAlgebra::VecND dlogS_da2;
	  VSAMath::PolyFit::dyda(exp_fit_param,x,dlogS_da2);

	  VSAAlgebra::VecND dlogE_da = 
	    -dlogS_da1/(y-yold)*dx - 
	    (dlogS_da1-dlogS_da2)*(logS-y)/std::pow(y-yold,2)*dx;

	  double logE_var = dlogE_da*(exp_fit_cov*dlogE_da);

	  double dlogS_dlogE = VSAMath::PolyFit::val(dfit_param,logE);

	  double logE_rms = 0;
	  double logE_rms_var = 0;
	    
	  VSNSpace::Point p;
	  s_space.midPointOfCellUnchecked(c, p);
	  p.x[1] = logE;

	  double logS_rms = s_rms_smooth[p];
	  double logS_rms_var = s_rms_var_smooth[p];

	  // std::cout << p.x[0] << " " << p.x[1] << " " << logS_rms
	  // 	    << std::endl;

	  // double logS_rms = VSAMath::PolyFit::val(rms_fit_param,logE);
	  // double logS_rms_var = VSAMath::PolyFit::var(rms_fit_param,
	  // 					       rms_fit_cov,logE);

	  // double logS_rms = 0;
	  // double logS_rms_var = 0;

	  // s_rms.interpolateWeight(p, logS_rms);
	  // s_rms_var.interpolateWeight(p, logS_rms_var);

	  logE_rms = logS_rms/dlogS_dlogE;
	  logE_rms_var = logS_rms_var/std::pow(dlogS_dlogE,2);

	  // std::cout << std::setw(20) << R
	  // 	    << std::setw(20) << logS
	  // 	    << std::setw(20) << logE_rms
	  // 	    << std::setw(20) << logE_rms_var
	  // 	    << std::setw(20) << logS_rms
	  // 	    << std::setw(20) << logS_rms2
	  // 	    << std::setw(20) << logS_rms_var
	  // 	    << std::setw(20) << logS_rms2_var
	  // 	    << std::endl;

	  //	  std::cout << logS << ' ' << logE << ' ' << logE_rms << '\n';

	  sfd.exp_hist.accumulate(logS,logE,logE_var); 
	  sfd.rms_hist.accumulate(logS,logE_rms,logE_rms_var);

	  e_exp[c] = logE;
	  e_exp_var[c] = logE_var;
	  e_rms[c] = logE_rms;
	  e_rms_var[c] = logE_rms_var;
	}

      size_fit_data.push_back(sfd);
    }
}

// ----------------------------------------------------------------------------
// VSEnergyLTLibraryVisitor::SizeFitData
// ----------------------------------------------------------------------------
void VSEnergyLTLibraryVisitor::SizeFitData::
save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeScalar("ndata",ndata);
  writer->writeScalar("R",R);
  writer->writeScalar("disp",disp);
  writer->writeScalar("exp_fit_nparam",exp_fit_nparam);
  writer->writeScalar("exp_fit_chi2",exp_fit_chi2);
  writer->writeScalar("exp_fit_rchi2",exp_fit_rchi2);  
  exp_fit_param.save(writer,"exp_fit_param");
  size_fit_exp_hist.save(writer->writeStruct("size_fit_exp_hist"));
  size_fit_rms_hist.save(writer->writeStruct("size_fit_rms_hist"));
  size_exp_hist.save(writer->writeStruct("size_exp_hist"));
  size_rms_hist.save(writer->writeStruct("size_rms_hist"));
  exp_hist.save(writer->writeStruct("exp_hist"));
  rms_hist.save(writer->writeStruct("rms_hist"));
}

// ----------------------------------------------------------------------------
// VSEnergyLTLibraryVisitor::LTData
// ----------------------------------------------------------------------------
void VSEnergyLTLibraryVisitor::LTData::
save(VSOctaveH5WriterStruct* writer) const
{
  VSNSpaceOctaveH5IO io;

  io.writeHistogram(writer->writeStruct("size_n"),m_size_n);
  io.writeHistogram(writer->writeStruct("size_exp"),m_size_exp);
  io.writeHistogram(writer->writeStruct("size_exp_var"),m_size_exp_var);
  io.writeHistogram(writer->writeStruct("size_rms"),m_size_rms);
  io.writeHistogram(writer->writeStruct("size_rms_var"),m_size_rms_var);
  io.writeHistogram(writer->writeStruct("energy_n"),m_energy_n);
  io.writeHistogram(writer->writeStruct("energy_mask"),m_energy_mask);
  io.writeHistogram(writer->writeStruct("energy_exp"),m_energy_exp);
  io.writeHistogram(writer->writeStruct("energy_exp_var"),m_energy_exp_var);
  io.writeHistogram(writer->writeStruct("energy_rms"),m_energy_rms);
  io.writeHistogram(writer->writeStruct("energy_rms_var"),m_energy_rms_var);

  writer->writeStructCellVector("energy_fit_exp_hist",m_energy_fit_exp_hist);
  writer->writeStructCellVector("energy_fit_rms_hist",m_energy_fit_rms_hist);
  writer->writeStructCellVector("energy_exp_hist",m_energy_exp_hist);
  writer->writeStructCellVector("energy_rms_hist",m_energy_rms_hist);


  const unsigned nfit_data = m_size_fit_data.size();
  VSOctaveH5WriterCellVector* wc = 
    writer->writeCellVector("polyfit",nfit_data);
  for(unsigned ifit = 0; ifit < nfit_data; ifit++)
    {
      VSOctaveH5WriterStruct* ws = wc->writeStruct(ifit);
      m_size_fit_data[ifit].save(ws);
      delete ws;
    }

  delete wc;
}

// ----------------------------------------------------------------------------
// VSEnergyLTLibraryVisitor::Data
// ----------------------------------------------------------------------------
VSEnergyLTLibraryVisitor::
Data::Data(const VSNSpace::Space& energy_space, 
	   const VSNSpace::Space& size_space, 
	   unsigned min_count,
	   double zn, double az, double ped): 
  zenith_deg(zn), azimuth_deg(az), ped_dev(ped), m_azel(), 
  m_array_lt(), m_scope_lt(),
  m_sp_all(energy_space,size_space,min_count,"array"), m_sp_tel()
{ 
  m_azel = SEphem::SphericalCoords::makeDeg(zn,az);
}

VSEnergyLTLibraryVisitor::Data::~Data()
{
  const unsigned nscope = m_sp_tel.size();
  for(unsigned iscope = 0; iscope < nscope; iscope++)
    {
      if(m_sp_tel[iscope]) delete m_sp_tel[iscope];
    }
}

// ----------------------------------------------------------------------------
// VSEnergyLTLibraryVisitor
// ----------------------------------------------------------------------------
VSEnergyLTLibraryVisitor::Options 
VSEnergyLTLibraryVisitor::s_default_options = 
VSEnergyLTLibraryVisitor::Options();

VSEnergyLTLibraryVisitor::
VSEnergyLTLibraryVisitor(const Options& opt):
  m_cuts(), m_egywt_calc(), m_ndim(),
  m_energy_space(), m_size_space(),
  m_zn_deg(), m_az_deg(), m_ped_dev(), m_offset_deg(),
  m_log10_egybin(), m_log10_egylo(), m_log10_egyhi(),
  m_data(), m_data_ptr(), m_wt(1), m_is_selected(),
  m_sim(), m_evt(), m_sim_header(),
  m_options(opt)
{
  if(opt.space == "N_R")
    {
      m_ndim = 2;
      m_energy_space = VSNSpace::Space(2);

      m_energy_space.axes[0] = VSNSpace::Axis(-10.0, 800.0, 20.0, 0, "R");
      m_energy_space.axes[1] = VSNSpace::Axis(  1.0,   6.0,  0.1, 0, "N");
    }
  else if(opt.space == "N_R_disp")
    {
      m_ndim = 3;
      m_energy_space = VSNSpace::Space(3);
      m_energy_space.axes[0] = VSNSpace::Axis(-20.0, 800.0, 40.0, 0, "R");
      m_energy_space.axes[1] = VSNSpace::Axis(  1.0,   6.0,  0.2, 0, "N");
      m_energy_space.axes[2] = VSNSpace::Axis( -0.1,   3.6,  0.2, 0, "fp_disp");
    }
  else if(opt.space == "lambdad_R")
    {
      m_ndim = 2;
      m_energy_space = VSNSpace::Space(2);

      m_energy_space.axes[0] = VSNSpace::Axis(-10.0, 800.0, 20.0, 0, "R");
      m_energy_space.axes[1] = VSNSpace::Axis(  6.0,  13.0,  0.1, 0, "lambdad");
    }
  else if(opt.space == "lambdad_G")
    {
      m_ndim = 2;
      m_energy_space = VSNSpace::Space(2);

      m_energy_space.axes[0] = VSNSpace::Axis( 50.0,1200.0, 50.0, 0, "G");
      m_energy_space.axes[1] = VSNSpace::Axis(  6.0,  13.0,  0.2, 0, "lambdad");
    }
  else if(opt.space == "lambdad_G_eta")
    {
      m_ndim = 3;
      m_energy_space = VSNSpace::Space(3);

      m_energy_space.axes[0] = VSNSpace::Axis( 50.0,1200.0, 50.0, 0, "G");
      m_energy_space.axes[1] = VSNSpace::Axis(  6.0,  13.0,  0.2, 0, "lambdad");
      m_energy_space.axes[2] = VSNSpace::Axis(  0.0,   5.0,  0.2, 0, "eta");
    }
  else
    {
      std::cerr << "Unrecognized space type: " << opt.space
		<< std::endl;
      exit(EXIT_FAILURE);
    }


  m_egywt_calc = new VSSimEnergyWeightCalc;
  m_cuts = new VSCutsEvaluator;
}

VSEnergyLTLibraryVisitor::~VSEnergyLTLibraryVisitor()
{
  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++) delete m_data[idata];

  delete m_egywt_calc;
  delete m_cuts;
}

void VSEnergyLTLibraryVisitor::
visitRun(const VSAnalysisStage1Data& stage1,
	 const VSTargetTable::Observation& obs,
	 const VSArrayMergedCalibrationData& cal)
{
  // Determine Offset and Pointing --------------------------------------------
  m_zn_deg = stage1.run_info.zn_mean_deg;
  m_az_deg = stage1.run_info.az_mean_deg;
  m_ped_dev = cal.mean_scaled_dev;
  m_offset_deg = stage1.sim_info->wobble_theta_deg;
  
  std::cout << "OFFSET  " << std::setw(20) << m_offset_deg 
	    << " ZN     " << std::setw(20) << m_zn_deg 
	    << " AZ     " << std::setw(20) << m_az_deg 
	    << " PEDDEV " << std::setw(20) << m_ped_dev << std::endl;

  m_data_ptr = NULL;

  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    {
      SphericalCoords azel =
	SphericalCoords::makeDeg(stage1.run_info.zn_mean_deg,
				 stage1.run_info.az_mean_deg);

      double dzn = fabs(m_zn_deg - m_data[idata]->zenith_deg);
      if(dzn < 0.5) m_zn_deg = m_data[idata]->zenith_deg;

      double dtheta = azel.separation(m_data[idata]->m_azel).deg();      
      if(dtheta < 0.5) m_az_deg = m_data[idata]->azimuth_deg;

      double dped = fabs(cal.mean_scaled_dev - m_data[idata]->ped_dev);

      if(dtheta < 0.5 && dped < 0.1)
	{
	  m_data_ptr = m_data[idata]; 
	  break;
	}
    }

  if(m_data_ptr == NULL)
    {
      m_data_ptr = new Data(m_energy_space,m_size_space,m_options.min_count,
			    m_zn_deg,m_az_deg,m_ped_dev);
      m_data.push_back(m_data_ptr);
    }

  m_scope_pos.resize(stage1.sim_info->scope_positions.size());
  for(unsigned iscope=0;iscope<m_scope_pos.size();iscope++)
    m_scope_pos[iscope] = 
      VSAAlgebra::Vec3D(stage1.sim_info->scope_positions[iscope].first,
			stage1.sim_info->scope_positions[iscope].second,
			stage1.sim_info->scope_positions[iscope].third);

}

void VSEnergyLTLibraryVisitor::leaveRun()
{

}

void VSEnergyLTLibraryVisitor::visitEvent(const VSEventArrayDatum& evt)
{
  m_evt = evt;

  if(m_cuts->isSelected(m_evt) && m_evt.theta1 < m_options.theta_cut &&
     std::isfinite(m_evt.theta1)) 
    m_is_selected = true;
  else m_is_selected = false;

  const unsigned nscope = evt.scope.size();

  if(nscope > m_data_ptr->m_sp_tel.size()) 
    m_data_ptr->m_sp_tel.resize(nscope,0);

  for(unsigned iscope=0;iscope<nscope;iscope++)
    if((!m_options.no_scope_tables)&&(evt.scope[iscope])&&
       !m_data_ptr->m_sp_tel[iscope])
      m_data_ptr->m_sp_tel[iscope] = 
	new Spaces(m_energy_space,m_size_space,m_options.min_count,
		   "T"+VSDataConverter::toString(iscope+1));

  if(m_options.no_reconstructed_impact)
    {
      m_sim_core.set(m_sim.core_east_m, 
		     m_sim.core_north_m, 
		     m_sim.core_elevation_asl_m);
      
      double szn = sin(m_sim.primary_zenith_deg*M_PI/180.0);
      m_sim_azel.set(szn*sin(m_sim.primary_azimuth_deg*M_PI/180.0),
		     szn*cos(m_sim.primary_azimuth_deg*M_PI/180.0),
		     cos(m_sim.primary_zenith_deg*M_PI/180.0));
    }
}

void VSEnergyLTLibraryVisitor::leaveEvent()
{

}

void VSEnergyLTLibraryVisitor::
visitScope(const VSEventScopeDatum& sd, unsigned iscope)
{
  // --------------------------------------------------------------------------
  // Apply selection and various quality cuts
  // --------------------------------------------------------------------------
  if(!m_is_selected) return;
  else if(sd.N <= 0  || (!sd.used_in_reconstruction && 
			 !m_options.no_reconstruction_cut)) return;
  else if(sd.fp_dist > m_options.dist_cut) return;

  // --------------------------------------------------------------------------
  // Set the distance with respect to this telescope (dr_perp) to
  // either its true or reconstructed value.
  // --------------------------------------------------------------------------
  double dr_perp = sd.R;
  if(m_options.no_reconstructed_impact)
    {
      VSAAlgebra::Vec3D dri = m_sim_core-m_scope_pos[iscope];
      double dri_e = dri*m_sim_azel;
      dr_perp = sqrt(dri*dri - dri_e*dri_e);
    }
  
  double log10_lambda = log10(sd.lambdad);
  double log10_size = log10(sd.N);
  double log10_energy = log10(m_sim.energy_tev);
    
  VSNSpace::Point p_size(m_ndim);
  VSNSpace::Point p_energy(m_ndim);
  if(m_options.space == "N_R")
    {
      p_energy.x[0] = dr_perp;
      p_energy.x[1] = log10_size;
      p_size.x[0] = dr_perp;
    }
  else if(m_options.space == "N_R_disp")
    {
      p_energy.x[0] = dr_perp;
      p_energy.x[1] = log10_size;
      p_energy.x[2] = sd.fp_disp;

      p_size.x[0] = dr_perp;
      p_size.x[2] = sd.fp_disp;
    }
  else if(m_options.space == "lambdad_R")
    {
      p_energy.x[0] = dr_perp;
      p_energy.x[1] = log10_lambda;

      p_size.x[0] = dr_perp;
    }
  else if(m_options.space == "lambdad_G")
    {
      p_energy.x[0] = sd.G;
      p_energy.x[1] = log10_lambda;
      p_size.x[0] = sd.G;
    }
  else if(m_options.space == "lambdad_eta")
    {
      double eta = -std::log(tan(sd.theta1/2));

      p_energy.x[0] = eta;
      p_energy.x[1] = log10_lambda;
      p_size.x[0] = eta;
    }
  else if(m_options.space == "lambdad_G_eta")
    {
      double eta = -std::log(tan(sd.theta1/2));

      p_energy.x[0] = sd.G;
      p_energy.x[1] = log10_lambda;
      p_energy.x[2] = eta;

      p_size.x[0] = sd.G;
      p_size.x[2] = eta;
    }

  p_size.x[1] = log10_energy;

  m_data_ptr->m_sp_all.accumulate(p_energy, log10_energy, 
				  p_size, p_energy.x[1] ,m_wt);
  
  vsassert(iscope < m_data_ptr->m_sp_tel.size());
  
  if(m_data_ptr->m_sp_tel[iscope])
    m_data_ptr->m_sp_tel[iscope]->
      accumulate(p_energy, log10_energy, p_size, p_energy.x[1] ,m_wt);
}

void VSEnergyLTLibraryVisitor::leaveScope()
{

}

void VSEnergyLTLibraryVisitor::
visitSimEvent(const VSArraySimulationDatum& sim)
{  
  m_wt = m_egywt_calc->calcWeight(sim.energy_tev);
  m_sim = sim;
}

void VSEnergyLTLibraryVisitor::leaveSimEvent()
{

}

void VSEnergyLTLibraryVisitor::
visitSimHeader(const VSHeaderSimulationDatum& header)
{ 
  m_size_space = VSNSpace::Space(m_ndim);
  
  double log10_egylo = log10(header.tables.front().energy_tev);
  double log10_egyhi = log10(header.tables.back().energy_tev);
  std::set< double > egy_set;
  unsigned negy = 0;

  for(std::vector< VSTableSimulationDatum >::const_iterator itr = 
	header.tables.begin(); itr != header.tables.end(); ++itr)
    {
      double log10_egy = std::log10(itr->energy_tev);

      if(log10_egy < log10_egylo) log10_egylo = log10_egy;
      if(log10_egy > log10_egyhi) log10_egyhi = log10_egy;


      std::set< double >::iterator itr = egy_set.begin();

      for(; itr != egy_set.end(); ++itr)
	if(fabs(log10_egy-*itr) < 0.001) break;

      if(itr == egy_set.end())
	{
	  egy_set.insert(log10_egy);
	  negy++;
	}
    }

  if(negy > 1)
    {
      m_log10_egybin = (log10_egyhi-log10_egylo)/(double)(negy-1);
      m_log10_egylo  = log10_egylo - m_log10_egybin/2.;
      m_log10_egyhi  = log10_egyhi + m_log10_egybin/2.;
    }

  for(unsigned idim=0;idim<m_ndim;idim++)if(idim!=1)
    m_size_space.axes[idim] = m_energy_space.axes[idim];
  m_size_space.axes[1] = 
    VSNSpace::Axis(m_log10_egylo, m_log10_egyhi,
		   m_log10_egybin, 0, "MC energy [log10 E/TeV]");

  m_egywt_calc->calcWeighting(header);
}

void VSEnergyLTLibraryVisitor::leaveSimHeader()
{

}

void VSEnergyLTLibraryVisitor::createLibrary()
{
  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    {
      std::cout << __PRETTY_FUNCTION__  
		<< ": Generating Energy LT " 
		<< " Zn  " << std::setw(15) << m_data[idata]->zenith_deg
		<< " Az  " << std::setw(15) << m_data[idata]->azimuth_deg
		<< " Ped " << std::setw(15) << m_data[idata]->ped_dev
		<< std::endl;

      createLibrary(m_data[idata]->m_sp_all,m_data[idata]->m_array_lt);

      const unsigned nscope = m_data[idata]->m_sp_tel.size();
      m_data[idata]->m_scope_lt.resize(nscope);

      for(unsigned iscope = 0; iscope < nscope; iscope++)
	{
	  if(!m_data[idata]->m_sp_tel[iscope]) continue;
	  createLibrary(*m_data[idata]->m_sp_tel[iscope],
			m_data[idata]->m_scope_lt[iscope]);
	}

    }
}

void VSEnergyLTLibraryVisitor::createLibrary(Spaces& sp, LTData& data)
{
  // --------------------------------------------------------------------------
  // Generate lookup table of size 
  // --------------------------------------------------------------------------
  if(m_options.median)
    sp.calculateMedianSizeTables(data.m_size_exp, data.m_size_exp_var, 
				 data.m_size_rms, data.m_size_rms_var);
  else
    sp.calculateMomentSizeTables(data.m_size_exp, data.m_size_exp_var, 
				 data.m_size_rms, data.m_size_rms_var);
      
  // --------------------------------------------------------------------------
  // Generate lookup table of energy 
  // --------------------------------------------------------------------------
  if(m_options.use_size_tables)
    {
      sp.size2Energy(data.m_size_exp, data.m_size_exp_var, 
		     data.m_size_rms, data.m_size_rms_var, sp.n_size,
		     data.m_energy_exp, data.m_energy_exp_var, 
		     data.m_energy_rms, data.m_energy_rms_var,
		     data.m_size_fit_data,
		     m_options.poly_rchi2_max_change,
		     m_options.poly_max_order);
    }
  else
    {
      data.m_energy_n = sp.n_energy;

      if(m_options.median)
	sp.calculateMedianEnergyTables(data.m_energy_exp, 
				       data.m_energy_exp_var,
				       data.m_energy_rms,
				       data.m_energy_rms_var);
      else
	sp.calculateMomentEnergyTables(data.m_energy_exp, 
				       data.m_energy_exp_var,
				       data.m_energy_rms,
				       data.m_energy_rms_var);
    }

  VSNSpace::Volume v = data.m_energy_exp_var>0;
  v.expandFromEdge();
  data.m_energy_mask = VSNSpace(v);
  data.m_energy_mask.setComment(m_options.space);
  v.expandFromEdge();

  VSNSpace log_rms = data.m_energy_rms;
  VSNSpace log_rms_var = data.m_energy_rms_var;

  VSNSpace::Cell c(data.m_energy_exp.space().cell());
  const unsigned ncell = data.m_energy_exp.space().size();

  for(unsigned icell = 0; icell < ncell; icell++)
    {
      data.m_energy_exp.space().cellOfIndexUnchecked(icell,c);      
      if(c.i[1] != 0) continue;

      VSLimitedErrorsHist<double,double> 
	hexp(data.m_energy_exp.space().axes[1].bin_size,
	     data.m_energy_exp.space().axes[1].lo_bound,
	     data.m_energy_exp.space().axes[1].hi_bound);

      VSLimitedErrorsHist<double,double> 
	hrms(data.m_energy_exp.space().axes[1].bin_size,
	     data.m_energy_exp.space().axes[1].lo_bound,
	     data.m_energy_exp.space().axes[1].hi_bound);

      for(c.i[1] = 0; c.i[1] < data.m_energy_exp.space().axes[1].nbin; c.i[1]++)
	{
	  if(data.m_energy_rms[c] > 0)
	    {
	      double rms = data.m_energy_rms[c];
	      log_rms[c] = std::log(rms);
	      log_rms_var[c] *= std::pow(rms,-2);
	    }

	  double logS = 
	    data.m_energy_exp.space().axes[1].midCoordUnchecked(c.i[1]);  
	  hexp.accumulate(logS,data.m_energy_exp[c],data.m_energy_exp_var[c]);
	  hrms.accumulate(logS,data.m_energy_rms[c],data.m_energy_rms_var[c]);
	}

      data.m_energy_exp_hist.push_back(hexp);
      data.m_energy_rms_hist.push_back(hrms);
    }

  if(!m_options.no_extrapolation)
    {
      sp.extrapolate(data.m_energy_exp, data.m_energy_exp_var); 
      sp.smooth(log_rms, log_rms_var,
		log_rms.axis(0).bin_size*4,log_rms.axis(1).bin_size*5,1);

      VSNSpace energy_exp_smooth = data.m_energy_exp;
      VSNSpace energy_exp_var_smooth = data.m_energy_exp_var;

      sp.smooth(energy_exp_smooth, energy_exp_var_smooth,
		log_rms.axis(0).bin_size*3,log_rms.axis(1).bin_size*4,1);

      for(unsigned icell = 0; icell < ncell; icell++)
	{
	  data.m_energy_exp.space().cellOfIndexUnchecked(icell,c); 
	  data.m_energy_rms[c] = std::exp(log_rms[c]);
	  if(data.m_energy_exp_var[c] == 0)
	    data.m_energy_exp[c] = energy_exp_smooth[c];
	}

      data.m_energy_exp.clearOutsideVolume(v);
      data.m_energy_rms.clearOutsideVolume(v);
    }


  for(unsigned icell = 0; icell < ncell; icell++)
    {
      data.m_energy_exp.space().cellOfIndexUnchecked(icell,c);      
      if(c.i[1] != 0) continue;

      VSLimitedErrorsHist<double,double> 
	hexp(data.m_energy_exp.space().axes[1].bin_size,
	     data.m_energy_exp.space().axes[1].lo_bound,
	     data.m_energy_exp.space().axes[1].hi_bound);

      VSLimitedErrorsHist<double,double> 
	hrms(data.m_energy_exp.space().axes[1].bin_size,
	     data.m_energy_exp.space().axes[1].lo_bound,
	     data.m_energy_exp.space().axes[1].hi_bound);

      for(c.i[1] = 0; c.i[1] < data.m_energy_exp.space().axes[1].nbin; c.i[1]++)
	{
	  double logS = 
	    data.m_energy_exp.space().axes[1].midCoordUnchecked(c.i[1]);

	  hexp.accumulate(logS,data.m_energy_exp[c],0.);
	  hrms.accumulate(logS,data.m_energy_rms[c],0.);
	}

      data.m_energy_fit_exp_hist.push_back(hexp);
      data.m_energy_fit_rms_hist.push_back(hrms);
    }

}

void VSEnergyLTLibraryVisitor::save(VSOctaveH5WriterStruct* writer) const
{
  std::cout << "Writing Library File..." << std::endl;

  VSScaledParameterLibraryWriter wlibrary(writer);

  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    {      
      wlibrary.write(m_data[idata]->m_array_lt.m_energy_exp,
		     m_data[idata]->m_array_lt.m_energy_mask,
		     VSScaledParameterLibraryWriter::SPS_ENERGY_EXPECTED,
		     false, m_data[idata]->zenith_deg, 
		     m_data[idata]->zenith_deg,
		     false, m_data[idata]->azimuth_deg, 
		     m_data[idata]->azimuth_deg,
		     true,0,
		     false, m_data[idata]->ped_dev,m_data[idata]->ped_dev);

      wlibrary.write(m_data[idata]->m_array_lt.m_energy_rms,
		     m_data[idata]->m_array_lt.m_energy_mask,
		     VSScaledParameterLibraryWriter::SPS_ENERGY_RMS,
		     false, m_data[idata]->zenith_deg, 
		     m_data[idata]->zenith_deg,
		     false, m_data[idata]->azimuth_deg, 
		     m_data[idata]->azimuth_deg,
		     true,0,
		     false, m_data[idata]->ped_dev,m_data[idata]->ped_dev);


      const unsigned nscope =  m_data[idata]->m_sp_tel.size();
      for(unsigned iscope = 0; iscope < nscope; iscope++)
	{
	  if(!m_data[idata]->m_sp_tel[iscope]) continue;

	  wlibrary.write(m_data[idata]->m_scope_lt[iscope].m_energy_exp,
			 VSScaledParameterLibraryWriter::SPS_ENERGY_EXPECTED,
			 false, m_data[idata]->zenith_deg, 
			 m_data[idata]->zenith_deg,
			 false, m_data[idata]->azimuth_deg, 
			 m_data[idata]->azimuth_deg,
			 false,iscope,
			 false, m_data[idata]->ped_dev,m_data[idata]->ped_dev);

	  wlibrary.write(m_data[idata]->m_scope_lt[iscope].m_energy_rms,
			 VSScaledParameterLibraryWriter::SPS_ENERGY_RMS,
			 false, m_data[idata]->zenith_deg, 
			 m_data[idata]->zenith_deg,
			 false, m_data[idata]->azimuth_deg, 
			 m_data[idata]->azimuth_deg,
			 false,iscope,
			 false, m_data[idata]->ped_dev,m_data[idata]->ped_dev);
	}      
    }

  if(!m_options.no_diagnostics)
    {
      VSOctaveH5WriterCellVector* wc = writer->writeCellVector("data", ndata);
      vsassert(wc);

      for(unsigned idata = 0; idata < ndata; idata++)
	{ 
	  VSOctaveH5WriterStruct* ws = wc->writeStruct(idata);
	  vsassert(ws);  

	  ws->writeScalar("zenith_deg",m_data[idata]->zenith_deg);
	  ws->writeScalar("azimuth_deg",m_data[idata]->azimuth_deg);
	  ws->writeScalar("ped_dev",m_data[idata]->ped_dev);

	  m_data[idata]->m_array_lt.save(ws);      

	  delete ws;
	}

      delete wc;      
    }
}

void VSEnergyLTLibraryVisitor::configure(VSOptions& options,
					 const std::string& profile, 
					 const std::string& opt_prefix)
{
  std::vector< std::string > sim_energy_weight;
  sim_energy_weight.push_back("powerlaw");
  sim_energy_weight.push_back("2.5");

  VSSimEnergyWeightCalc::getDefaultOptions().sim_energy_weight =
    sim_energy_weight;

  VSSimEnergyWeightCalc::configure(options);
  VSCutsEvaluator::configure(options);

  options.findWithValue("min_events", s_default_options.min_count,
			"Minimum number of events per cell.");  

  options.findBoolValue("use_size_tables", 
			s_default_options.use_size_tables, true,
			"Accumulate sized based tables (i.e. Gernot's "
			"method A) and extrapolate energy tables from them.");

  options.findWithValue("poly_max_order", s_default_options.poly_max_order,
			"Set largest degree polynomial which can be fitted "
			"to the size-based lookup tables during inversion.");

  options.findWithValue("poly_rchi2_max_change", 
			s_default_options.poly_rchi2_max_change,
			"Set amount to which reduced chi-squared can be "
			"increased from its minimum when fitting size-based "
			"tables with polynomials during inversion.");

  options.findBoolValue("median", s_default_options.median, true,
			"Calculate tables with median and 68% containment "
			"intervals rather than mean and variance.");

  options.findBoolValue("no_scope_tables", 
			s_default_options.no_scope_tables, true,
			"Do not generate the per-telescope lookup tables.");

  options.findWithValue("space", s_default_options.space,
   			"Define the set of parameters to be used for the "
			"energy lookup table.  Options are: \"N_R\" (Size,"
			"Impact Distance), \"N_R_disp\" (Size,Impact Distance,"
			"Image Displacemet).");  

  options.findBoolValue("no_reconstruction_cut", 
			s_default_options.no_reconstruction_cut, true,
			"Do not apply reconstruction cuts to events being "
			"entered in the lookup tables.");

  options.findBoolValue("no_reconstructed_impact", 
			s_default_options.no_reconstructed_impact, true,
			"Do not use reconstructed impact distance, instead "
			"used actual impact distance of simulated particle "
			"on the ground");

  options.findWithValue("theta_cut", s_default_options.theta_cut,
			"Cut on the maximum separation in degrees between "
			"the true and reconstructed primary direction.");

  options.findWithValue("dist_cut", s_default_options.dist_cut,
			"Cut on the maximum distance of the image centroid "
			"from the center of the FOV.");

  // options.findBoolValue("no_extrapolation", 
  // 			s_default_options.no_extrapolation, true,
  // 			"Do not extrapolate.");

  options.findBoolValue("no_diagnostics", 
			s_default_options.no_diagnostics, true,
			"Do not save extra diagnostic data "
			"associated with lookup table generation.");

  options.findWithValue("min_events", 
			s_default_options.min_count,
			"Minimum number of events per cell.");  
}

// ============================================================================
// VSScaledParameterLibraryVisitor
// ============================================================================
#define CALC(M1,M2)					\
  VSNSpace::Volume mask = n>min_count;			\
  M1.nspace(exp);					\
  exp.clearOutsideVolume(mask);				\
  M2.nspace(rms);					\
  rms.clearOutsideVolume(mask);

// ----------------------------------------------------------------------------
// VSScaledParameterLibraryVisitor::Spaces
// ----------------------------------------------------------------------------
VSScaledParameterLibraryVisitor::
Spaces::Spaces(const VSNSpace::Space& def, 
	       unsigned _min_count,
	       const std::string& _comment_base):
  min_count(_min_count), comment_base(_comment_base),
  n(def), 
  width_mean(def),  width_rms(def),
  length_mean(def), length_rms(def),
  disp_mean(def),   disp_rms(def),
  width_med(def),   width_i68(def),
  length_med(def),  length_i68(def),
  disp_med(def),    disp_i68(def)
{ /* nothing to see here */ }

void VSScaledParameterLibraryVisitor::
Spaces::accumulate(const VSNSpace::Point& p, 
		   double width, double length, double disp, 
		   double weight)
{
  try
    {
      n[p]++;
      
      width_mean[p].accumulate(width,weight);
      width_rms[p].accumulate(width,weight);
      length_mean[p].accumulate(length,weight);
      length_rms[p].accumulate(length,weight);
      disp_mean[p].accumulate(disp,weight);
      disp_rms[p].accumulate(disp,weight);
      
      width_med[p].accumulate(width,weight);
      width_i68[p].accumulate(width,weight);
      length_med[p].accumulate(length,weight);
      length_i68[p].accumulate(length,weight);
      disp_med[p].accumulate(disp,weight);
      disp_i68[p].accumulate(disp,weight);
    }
  catch(std::out_of_range&)
    {
      // nothing to see here
    }
}

void VSScaledParameterLibraryVisitor::
Spaces::calculateMomentWidthTables(VSNSpace& exp, VSNSpace& exp_var,
				   VSNSpace& rms, VSNSpace& rms_var)
{
  CALC(width_mean, width_rms);
  exp.setComment(comment_base+" width mean");
  rms.setComment(comment_base+" width standard deviation");
}

void VSScaledParameterLibraryVisitor::
Spaces::calculateMomentLengthTables(VSNSpace& exp, VSNSpace& exp_var,
				    VSNSpace& rms, VSNSpace& rms_var)
{
  CALC(length_mean, length_rms);
  exp.setComment(comment_base+" length mean");
  rms.setComment(comment_base+" length standard deviation");
}

void VSScaledParameterLibraryVisitor::
Spaces::calculateMomentDispTables(VSNSpace& exp, VSNSpace& exp_var,
				  VSNSpace& rms, VSNSpace& rms_var)
{
  CALC(disp_mean, disp_rms);
  exp.setComment(comment_base+" disp mean");
  rms.setComment(comment_base+" disp standard deviation");
}

void VSScaledParameterLibraryVisitor::
Spaces::calculateMedianWidthTables(VSNSpace& exp, VSNSpace& exp_var,
				   VSNSpace& rms, VSNSpace& rms_var)
{
  CALC(width_med, width_i68);

  width_rms.nspace(exp_var);
  width_rms.nspace(rms_var);

  exp_var *= exp_var;
  exp_var *= 1.57;
  exp_var /= n;

  rms_var *= rms_var;
  rms_var *= 0.5;
  rms_var /= n;
  
  exp_var.clearOutsideVolume(n>min_count);
  rms_var.clearOutsideVolume(n>min_count);

  exp.setComment(comment_base+" width median");
  rms.setComment(comment_base+" width 68% containment interval");
}

void VSScaledParameterLibraryVisitor::
Spaces::calculateMedianLengthTables(VSNSpace& exp, VSNSpace& exp_var,
				    VSNSpace& rms, VSNSpace& rms_var)
{
  CALC(length_med, length_i68);

  length_rms.nspace(exp_var);
  length_rms.nspace(rms_var);

  exp_var *= exp_var;
  exp_var *= 1.57;
  exp_var /= n;

  rms_var *= rms_var;
  rms_var *= 0.5;
  rms_var /= n;

  exp_var.clearOutsideVolume(n>min_count);
  rms_var.clearOutsideVolume(n>min_count);

  exp.setComment(comment_base+" length median");
  rms.setComment(comment_base+" length 68% containment interval");
}

void VSScaledParameterLibraryVisitor::
Spaces::calculateMedianDispTables(VSNSpace& exp, VSNSpace& exp_var,
				  VSNSpace& rms, VSNSpace& rms_var)
{
  CALC(disp_med, disp_i68);

  disp_rms.nspace(exp_var);
  disp_rms.nspace(rms_var);

  exp_var *= exp_var;
  exp_var *= 1.57;
  exp_var /= n;

  rms_var *= rms_var;
  rms_var *= 0.5;
  rms_var /= n;

  exp_var.clearOutsideVolume(n>min_count);
  rms_var.clearOutsideVolume(n>min_count);

  exp.setComment(comment_base+" disp median");
  rms.setComment(comment_base+" disp 68% containment interval");
}

// ----------------------------------------------------------------------------
// VSScaledParameterLibraryVisitor::Data 
// ----------------------------------------------------------------------------
VSScaledParameterLibraryVisitor::LTData::LTData(): 
  m_width_exp(), m_width_exp_var(), m_width_rms(), m_width_rms_var(),  
  m_width_fit_exp(), m_width_fit_rms(), 
  m_length_exp(), m_length_exp_var(), m_length_rms(), m_length_rms_var(),
  m_length_fit_exp(), m_length_fit_rms(), 
  m_disp_exp(), m_disp_exp_var(), m_disp_rms(), m_disp_rms_var(),
  m_disp_fit_exp(), m_disp_fit_rms(), 
  m_width_fit_exp_hist(), m_width_fit_rms_hist(),
  m_width_exp_hist(), m_width_rms_hist(),
  m_length_fit_exp_hist(), m_length_fit_rms_hist(),
  m_length_exp_hist(), m_length_rms_hist(),
  m_disp_fit_exp_hist(), m_disp_fit_rms_hist(),
  m_disp_exp_hist(), m_disp_rms_hist()
{ }

void VSScaledParameterLibraryVisitor::
LTData::save(VSOctaveH5WriterStruct* ws) const
{

  VSNSpaceOctaveH5IO io;

  io.writeHistogram(ws->writeStruct("width_mask"),m_width_mask);
  io.writeHistogram(ws->writeStruct("width_exp"),m_width_exp);
  io.writeHistogram(ws->writeStruct("width_exp_var"),m_width_exp_var);
  io.writeHistogram(ws->writeStruct("width_rms"),m_width_rms);
  io.writeHistogram(ws->writeStruct("width_rms_var"),m_width_rms_var);
  io.writeHistogram(ws->writeStruct("width_fit_exp"),m_width_fit_exp);
  io.writeHistogram(ws->writeStruct("width_fit_rms"),m_width_fit_rms);

  io.writeHistogram(ws->writeStruct("length_mask"),m_length_mask);
  io.writeHistogram(ws->writeStruct("length_exp"),m_length_exp);
  io.writeHistogram(ws->writeStruct("length_exp_var"),m_length_exp_var);
  io.writeHistogram(ws->writeStruct("length_rms"),m_length_rms);
  io.writeHistogram(ws->writeStruct("length_rms_var"),m_length_rms_var);
  io.writeHistogram(ws->writeStruct("length_fit_exp"),m_length_fit_exp);
  io.writeHistogram(ws->writeStruct("length_fit_rms"),m_length_fit_rms);

  io.writeHistogram(ws->writeStruct("disp_mask"),m_disp_mask);
  io.writeHistogram(ws->writeStruct("disp_exp"),m_disp_exp);
  io.writeHistogram(ws->writeStruct("disp_exp_var"),m_disp_exp_var);
  io.writeHistogram(ws->writeStruct("disp_rms"),m_disp_rms);
  io.writeHistogram(ws->writeStruct("disp_rms_var"),m_disp_rms_var);
  io.writeHistogram(ws->writeStruct("disp_fit_exp"),m_disp_fit_exp);
  io.writeHistogram(ws->writeStruct("disp_fit_rms"),m_disp_fit_rms);

  ws->writeStructCellVector("width_exp_hist",m_width_exp_hist);
  ws->writeStructCellVector("width_fit_exp_hist",m_width_fit_exp_hist);
  ws->writeStructCellVector("width_rms_hist",m_width_rms_hist);
  ws->writeStructCellVector("width_fit_rms_hist",m_width_fit_rms_hist);

  ws->writeStructCellVector("length_exp_hist",m_length_exp_hist);
  ws->writeStructCellVector("length_fit_exp_hist",m_length_fit_exp_hist);
  ws->writeStructCellVector("length_rms_hist",m_length_rms_hist);
  ws->writeStructCellVector("length_fit_rms_hist",m_length_fit_rms_hist);
  
  ws->writeStructCellVector("disp_exp_hist",m_disp_exp_hist);
  ws->writeStructCellVector("disp_fit_exp_hist",m_disp_fit_exp_hist);
  ws->writeStructCellVector("disp_rms_hist",m_disp_rms_hist);
  ws->writeStructCellVector("disp_fit_rms_hist",m_disp_fit_rms_hist);
}

VSScaledParameterLibraryVisitor::
Data::Data(const VSNSpace::Space& space, unsigned min_count,
	   double zn, double az, double ped): 
  zenith_deg(zn), azimuth_deg(az), ped_dev(ped), m_azel(), 
  m_array_lt(), m_scope_lt(),
  m_sp_all(space,min_count,"array"), m_sp_tel()
{ 
  m_azel = SEphem::SphericalCoords::makeDeg(zn,az);
}

VSScaledParameterLibraryVisitor::Data::~Data()
{
  const unsigned nscope = m_sp_tel.size();
  for(unsigned iscope = 0; iscope < nscope; iscope++)
    {
      if(m_sp_tel[iscope]) delete m_sp_tel[iscope];
    }
}

// ----------------------------------------------------------------------------
// VSScaledParameterLibraryVisitor
// ----------------------------------------------------------------------------
VSScaledParameterLibraryVisitor::Options 
VSScaledParameterLibraryVisitor::s_default_options = 
VSScaledParameterLibraryVisitor::Options();

VSScaledParameterLibraryVisitor::
VSScaledParameterLibraryVisitor(const Options& opt):
  m_cuts(), m_egywt_calc(),
  m_space(), m_zn_deg(),
  m_az_deg(), m_ped_dev(), m_offset_deg(),
  m_data(), m_data_ptr(), m_wt(1), m_is_selected(),
  m_sim(), m_evt(), m_sim_header(),
  m_options(opt)
{
  m_space = VSScaledParameterLibraryWriter::defaultSpace(m_options.ndim);

  m_egywt_calc = new VSSimEnergyWeightCalc;
  m_cuts = new VSCutsEvaluator;
}

VSScaledParameterLibraryVisitor::~VSScaledParameterLibraryVisitor()
{
  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++) delete m_data[idata];

  delete m_egywt_calc;
  delete m_cuts;
}

void VSScaledParameterLibraryVisitor::
visitRun(const VSAnalysisStage1Data& stage1,
	 const VSTargetTable::Observation& obs,
	 const VSArrayMergedCalibrationData& cal)
{
  // Determine Offset and Pointing --------------------------------------------
  m_zn_deg = stage1.run_info.zn_mean_deg;
  m_az_deg = stage1.run_info.az_mean_deg;
  m_ped_dev = cal.mean_scaled_dev;
  m_offset_deg = stage1.sim_info->wobble_theta_deg;
  
  std::cout << "OFFSET  " << std::setw(20) << m_offset_deg 
	    << " ZN     " << std::setw(20) << m_zn_deg 
	    << " AZ     " << std::setw(20) << m_az_deg 
	    << " PEDDEV " << std::setw(20) << m_ped_dev << std::endl;

  m_data_ptr = NULL;

  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    {
      SphericalCoords azel =
	SphericalCoords::makeDeg(stage1.run_info.zn_mean_deg,
				 stage1.run_info.az_mean_deg);

      double dzn = fabs(m_zn_deg - m_data[idata]->zenith_deg);
      if(dzn < 0.5) m_zn_deg = m_data[idata]->zenith_deg;

      double dtheta = azel.separation(m_data[idata]->m_azel).deg();
      
      if(dtheta < 0.5) m_az_deg = m_data[idata]->azimuth_deg;

      double dped =
	fabs(cal.mean_scaled_dev - m_data[idata]->ped_dev);
      
      // std::cout << stage1.run_info.az_mean_deg << " "
      // 		<< m_data[idata]->azimuth_deg << " "
      // 		<< stage1.run_info.zn_mean_deg << " "
      // 		<< m_data[idata]->zenith_deg << " " << dtheta 
      // 		<< std::endl;

      if(dtheta < 0.5 && dped < 0.1)
	{
	  m_data_ptr = m_data[idata]; 
	  break;
	}
    }

  if(m_data_ptr == NULL)
    {
      m_data_ptr = new Data(m_space,m_options.min_count,
			    m_zn_deg,m_az_deg,m_ped_dev);
      m_data.push_back(m_data_ptr);
    }
}

void VSScaledParameterLibraryVisitor::leaveRun()
{

}

void VSScaledParameterLibraryVisitor::visitEvent(const VSEventArrayDatum& evt)
{
  m_evt = evt;

  if(m_cuts->isSelected(m_evt) &&
     m_evt.theta1 < m_options.theta_cut) m_is_selected = true;
  //if(m_cuts->isSelected(m_evt)) m_is_selected = true;
  else m_is_selected = false;

  const unsigned nscope = evt.scope.size();

  if(nscope > m_data_ptr->m_sp_tel.size()) 
    m_data_ptr->m_sp_tel.resize(nscope,0);

  for(unsigned iscope=0;iscope<nscope;iscope++)
    if((!m_options.no_scope_tables)&&(evt.scope[iscope])&&
       !m_data_ptr->m_sp_tel[iscope])
      m_data_ptr->m_sp_tel[iscope] = 
	new Spaces(m_space,m_options.min_count,
		   "T"+VSDataConverter::toString(iscope+1));
}

void VSScaledParameterLibraryVisitor::leaveEvent()
{

}

void VSScaledParameterLibraryVisitor::
visitScope(const VSEventScopeDatum& sd, unsigned iscope)
{
  if(!m_is_selected) return;

  VSNSpace::Point p;
  VSScaledParameterLibraryWriter::point(p,sd.R,sd.N,sd.fp_disp,
					m_options.ndim);
		    
  double fp_width  = sd.intrinsic_width;
  double fp_length = sd.intrinsic_length;
  double fp_disp   = sd.fp_disp;

  if(sd.used_in_reconstruction)
    m_data_ptr->m_sp_all.accumulate(p,fp_width,fp_length,fp_disp,m_wt);

  vsassert(iscope < m_data_ptr->m_sp_tel.size());

  if(m_data_ptr->m_sp_tel[iscope] && sd.used_in_reconstruction)
    m_data_ptr->m_sp_tel[iscope]->
      accumulate(p,fp_width,fp_length,fp_disp,m_wt);
}

void VSScaledParameterLibraryVisitor::leaveScope()
{

}

void VSScaledParameterLibraryVisitor::
visitSimEvent(const VSArraySimulationDatum& sim)
{  
  m_wt = m_egywt_calc->calcWeight(sim.energy_tev);
}

void VSScaledParameterLibraryVisitor::leaveSimEvent()
{

}

void VSScaledParameterLibraryVisitor::
visitSimHeader(const VSHeaderSimulationDatum& header)
{ 
  m_egywt_calc->calcWeighting(header);
}

void VSScaledParameterLibraryVisitor::leaveSimHeader()
{

}

void VSScaledParameterLibraryVisitor::createLibrary()
{
  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    {
      std::cout << __PRETTY_FUNCTION__  
		<< ": Generating Scaled Parameter LT " 
		<< " Zn  " << std::setw(15) << m_data[idata]->zenith_deg
		<< " Az  " << std::setw(15) << m_data[idata]->azimuth_deg
		<< " Ped " << std::setw(15) << m_data[idata]->ped_dev
		<< std::endl;

      createLibrary(m_data[idata]->m_sp_all,m_data[idata]->m_array_lt);

      const unsigned nscope = m_data[idata]->m_sp_tel.size();
      m_data[idata]->m_scope_lt.resize(nscope);

      for(unsigned iscope = 0; iscope < nscope; iscope++)
	{
	  if(!m_data[idata]->m_sp_tel[iscope]) continue;
	  createLibrary(*m_data[idata]->m_sp_tel[iscope],
			m_data[idata]->m_scope_lt[iscope]);
	}

      fit(m_data[idata]->m_array_lt);
    }
}

void VSScaledParameterLibraryVisitor::fit(LTData& data)
{
  fit(data.m_width_exp, data.m_width_exp_var, data.m_width_fit_exp,4);
  smooth(data.m_width_rms, data.m_width_rms_var, data.m_width_fit_rms);

  VSNSpace::Volume width_mask = data.m_width_exp_var>0;
  width_mask.expandFromEdge();
  data.m_width_mask = VSNSpace(width_mask);
  width_mask.expandFromEdge();
  data.m_width_fit_exp.clearOutsideVolume(width_mask);
  data.m_width_fit_rms.clearOutsideVolume(width_mask);

  fit(data.m_length_exp,data.m_length_exp_var,data.m_length_fit_exp,3);
  smooth(data.m_length_rms,data.m_length_rms_var,data.m_length_fit_rms);

  VSNSpace::Volume length_mask = data.m_length_exp_var>0;
  length_mask.expandFromEdge();
  data.m_length_mask = VSNSpace(length_mask);
  length_mask.expandFromEdge();
  data.m_length_fit_exp.clearOutsideVolume(length_mask);
  data.m_length_fit_rms.clearOutsideVolume(length_mask);


  fit(data.m_disp_exp,data.m_disp_exp_var,data.m_disp_fit_exp,3);
  smooth(data.m_disp_rms,data.m_disp_rms_var,data.m_disp_fit_rms);

  VSNSpace::Volume disp_mask = data.m_disp_exp_var>0;
  disp_mask.expandFromEdge();
  data.m_disp_mask = VSNSpace(disp_mask);
  disp_mask.expandFromEdge();
  data.m_disp_fit_exp.clearOutsideVolume(disp_mask);
  data.m_disp_fit_rms.clearOutsideVolume(disp_mask);


  const unsigned nbin = data.m_width_exp.axis(0).nbin;
  std::set<unsigned> dims;
  dims.insert(0);
  
  data.m_width_fit_exp_hist.resize(nbin);
  data.m_width_fit_rms_hist.resize(nbin);
  data.m_width_exp_hist.resize(nbin);
  data.m_width_rms_hist.resize(nbin);
  
  data.m_length_fit_exp_hist.resize(nbin);
  data.m_length_fit_rms_hist.resize(nbin);
  data.m_length_exp_hist.resize(nbin);
  data.m_length_rms_hist.resize(nbin);

  data.m_disp_fit_exp_hist.resize(nbin);
  data.m_disp_fit_rms_hist.resize(nbin);
  data.m_disp_exp_hist.resize(nbin);
  data.m_disp_rms_hist.resize(nbin);

  for(unsigned i = 0; i < nbin; i++)
    {
      std::vector<unsigned> indices;
      indices.push_back(i);
      data.m_width_exp.
	hist(dims,indices,data.m_width_exp_hist[i].contentsHist());
      data.m_width_exp_var.
	hist(dims,indices,data.m_width_exp_hist[i].varianceHist());

      data.m_width_fit_exp.
	hist(dims,indices,data.m_width_fit_exp_hist[i].contentsHist());
      data.m_width_fit_exp_hist[i].varianceHist() = 
	data.m_width_fit_exp_hist[i].contentsHist();
      data.m_width_fit_exp_hist[i].varianceHist().clear();
      data.m_width_fit_exp_hist[i].varianceHist().fill(0.);

      data.m_width_rms.
	hist(dims,indices,data.m_width_rms_hist[i].contentsHist());
      data.m_width_rms_var.
	hist(dims,indices,data.m_width_rms_hist[i].varianceHist());

      data.m_width_fit_rms.
	hist(dims,indices,data.m_width_fit_rms_hist[i].contentsHist());
      data.m_width_fit_rms_hist[i].varianceHist() = 
	data.m_width_fit_rms_hist[i].contentsHist();
      data.m_width_fit_rms_hist[i].varianceHist().clear();
      data.m_width_fit_rms_hist[i].varianceHist().fill(0.);

      data.m_length_exp.
	hist(dims,indices,data.m_length_exp_hist[i].contentsHist());
      data.m_length_exp_var.
	hist(dims,indices,data.m_length_exp_hist[i].varianceHist());

      data.m_length_fit_exp.
	hist(dims,indices,data.m_length_fit_exp_hist[i].contentsHist());
      data.m_length_fit_exp_hist[i].varianceHist() = 
	data.m_length_fit_exp_hist[i].contentsHist();
      data.m_length_fit_exp_hist[i].varianceHist().clear();
      data.m_length_fit_exp_hist[i].varianceHist().fill(0.);

      data.m_length_rms.
	hist(dims,indices,data.m_length_rms_hist[i].contentsHist());
      data.m_length_rms_var.
	hist(dims,indices,data.m_length_rms_hist[i].varianceHist());

      data.m_length_fit_rms.
	hist(dims,indices,data.m_length_fit_rms_hist[i].contentsHist());
      data.m_length_fit_rms_hist[i].varianceHist() = 
	data.m_length_fit_rms_hist[i].contentsHist();
      data.m_length_fit_rms_hist[i].varianceHist().clear();
      data.m_length_fit_rms_hist[i].varianceHist().fill(0.);

      data.m_disp_exp.
	hist(dims,indices,data.m_disp_exp_hist[i].contentsHist());
      data.m_disp_exp_var.
	hist(dims,indices,data.m_disp_exp_hist[i].varianceHist());

      data.m_disp_fit_exp.
	hist(dims,indices,data.m_disp_fit_exp_hist[i].contentsHist());
      data.m_disp_fit_exp_hist[i].varianceHist() = 
	data.m_disp_fit_exp_hist[i].contentsHist();
      data.m_disp_fit_exp_hist[i].varianceHist().clear();
      data.m_disp_fit_exp_hist[i].varianceHist().fill(0.);
    }

  data.m_width_exp = data.m_width_fit_exp;
  data.m_width_rms = data.m_width_fit_rms;
  data.m_length_exp = data.m_length_fit_exp;
  data.m_length_rms = data.m_length_fit_rms;
  data.m_disp_exp = data.m_disp_fit_exp;
  data.m_disp_rms = data.m_disp_fit_rms;
}

void VSScaledParameterLibraryVisitor::createLibrary(Spaces& sp, LTData& data)
{
  if(m_options.median)
    {
      sp.calculateMedianWidthTables(data.m_width_exp, 
				    data.m_width_exp_var, 
				    data.m_width_rms,
				    data.m_width_rms_var);
      sp.calculateMedianLengthTables(data.m_length_exp, 
				     data.m_length_exp_var, 
				     data.m_length_rms,
				     data.m_length_rms_var);
      sp.calculateMedianDispTables(data.m_disp_exp, 
				   data.m_disp_exp_var, 
				   data.m_disp_rms,
				   data.m_disp_rms_var);      
    }
  else
    {
      sp.calculateMomentWidthTables(data.m_width_exp, 
				    data.m_width_exp_var, 
				    data.m_width_rms,
				    data.m_width_rms_var);
      sp.calculateMomentLengthTables(data.m_length_exp, 
				     data.m_length_exp_var, 
				     data.m_length_rms,
				     data.m_length_rms_var);
      sp.calculateMomentDispTables(data.m_disp_exp, 
				   data.m_disp_exp_var, 
				   data.m_disp_rms,
				   data.m_disp_rms_var);   
    }
}

void VSScaledParameterLibraryVisitor::
smooth(VSNSpace& exp, VSNSpace& var, VSNSpace& fit)
{
  fit = exp;

  VSAMath::Data<VSAAlgebra::Vec2D> xys;

  VSNSpace::Cell c(2);
  // Loop on Impact Parameter -------------------------------------------------
  for(c.i[0] = 0; c.i[0] < exp.axis(0).nbin; c.i[0]++)
    {
      double x0 = exp.axis(0).midCoordUnchecked(c.i[0]);	  
      for(c.i[1] = 0; c.i[1] < exp.axis(1).nbin; c.i[1]++)
	{
	  double x1 = exp.axis(1).midCoordUnchecked(c.i[1]);
	  if(var[c] <= 0) continue;
	  
	  double log_exp = std::log(exp[c]);
	  double log_exp_err = sqrt(var[c])/exp[c];

	  VSAMath::DataPoint<VSAAlgebra::Vec2D> 
	    dp(VSAAlgebra::Vec2D(x0,x1),log_exp,log_exp_err);
	  xys.insert(dp);
	}
    }

  const double dx1 = 100.;
  const double dx2 = 1.0;

  VSAMath::LocalRegression2D<VSAAlgebra::Vec2D> lr(xys,dx1,dx2,1);

  // Loop on Impact Parameter -------------------------------------------------
  for(c.i[0] = 0; c.i[0] < exp.axis(0).nbin; c.i[0]++)
    {
      double x0 = exp.axis(0).midCoordUnchecked(c.i[0]);	  
      for(c.i[1] = 0; c.i[1] < exp.axis(1).nbin; c.i[1]++)
	{
	  double x1 = exp.axis(1).midCoordUnchecked(c.i[1]);
	  VSAAlgebra::Vec2D x(x0,x1);	  
	  fit[c] = std::exp(lr.val(x));
	}
    }
}

void VSScaledParameterLibraryVisitor::
fit(VSNSpace& exp, VSNSpace& var, VSNSpace& fit, unsigned npoly)
{
  VSNSpace fit_var = exp;
  fit_var.clear();

  fit = exp;
  VSAMath::Data<VSAAlgebra::Vec2D> fit_xys;
  VSNSpace::Cell c(2);
  // Loop on Impact Parameter -------------------------------------------------
  for(c.i[0] = 0; c.i[0] < exp.axis(0).nbin; c.i[0]++)
    {
      //      double R = exp.axis(0).midCoordUnchecked(c.i[0]);
      VSAMath::Data<double> fit_size_data;

      for(c.i[1] = 0; c.i[1] < exp.axis(1).nbin; c.i[1]++)
	{
	  VSNSpace::Point p;
	  exp.space().midPointOfCellUnchecked(c,p);

	  VSAMath::DataPoint<double> dp(p.x[1],exp[c],sqrt(var[c]));

	  if(exp[c] > 0) fit_size_data.insert(dp);

	  // if(exp[c] > 0)
	  //   {
	  //     std::cout << p.x[1] << " " << exp[c] << " " << sqrt(var[c])
	  // 		<< std::endl;
	  //   }
	}

      if(fit_size_data.size() < 3) continue;

      unsigned nparam = std::min(npoly+1,fit_size_data.size()-1);

      //      std::cout << R << " " << nparam << std::endl;

      VSAAlgebra::MatrixND  cov;
      VSAAlgebra::VecND     param;
      VSAFunction::Poly fn(nparam-1);

      VSAMath::Fitcsvd<VSAFunction::Poly,double> 
	fitcsvd(fit_size_data,VSAFunction::Poly(nparam-1));

      fitcsvd.addMonotonicConstraints(2,0.2,20);

      try
	{
	  fitcsvd.fit();
	}
      catch(const std::exception& e)
	{
	  std::cerr << e.what() << std::endl;
	  continue;
	}

      param = fitcsvd.param();
      cov = fitcsvd.cov();

      // for(unsigned i = 0; i < param.ndim(); i++)
      // 	{
      // 	  std::cout << std::setw(5) << i
      // 		    << std::setw(15) << param(i)
      // 		    << std::setw(15) << sqrt(cov(i,i))
      // 		    << std::setw(15) << cov(i,i)
      // 		    << std::endl;
      // 	  //if(cov(i,i) < 0) exit(1);
      // 	}

      for(c.i[1] = 0; c.i[1] < exp.axis(1).nbin; c.i[1]++)
	{
	  double x0 = exp.axis(0).midCoordUnchecked(c.i[0]);
	  double x1 = exp.axis(1).midCoordUnchecked(c.i[1]);

	  double val = fn.val(x1,param);
	  double err = sqrt(VSAMath::PolyFit::var(param,cov,x1));

	  // std::cout << x0 << " " << x1 << " " 
	  // 	    << val << " " << err << std::endl;

	  if(var[c] > 0)
	    {
	      fit[c] = val;	  
	      double log_exp = std::log(val);
	      double log_exp_err = err/val;

	      VSAMath::DataPoint<VSAAlgebra::Vec2D> 
		dp(VSAAlgebra::Vec2D(x0,x1),log_exp,log_exp_err);
	      fit_xys.insert(dp);
	    }
	}

    }

  const double dx1 = 100.;
  const double dx2 = 1.0;

  VSAMath::LocalRegression2D<VSAAlgebra::Vec2D> lr(fit_xys,dx1,dx2,1);

  for(c.i[0] = 0; c.i[0] < exp.axis(0).nbin; c.i[0]++)
    {
      for(c.i[1] = 0; c.i[1] < exp.axis(1).nbin; c.i[1]++)
	{
	  if(fit[c] > 0) continue;
	  double x0 = exp.axis(0).midCoordUnchecked(c.i[0]);
	  double x1 = exp.axis(1).midCoordUnchecked(c.i[1]);

	  fit[c] = std::exp(lr.val(VSAAlgebra::Vec2D(x0,x1)));
	}
    }
}

void VSScaledParameterLibraryVisitor::save(VSOctaveH5WriterStruct* writer) 
  const
{
  std::cout << "Writing Library File..." << std::endl;

  VSScaledParameterLibraryWriter wlibrary(writer);

  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    {   
      wlibrary.write(m_data[idata]->m_array_lt.m_width_exp,
		     m_data[idata]->m_array_lt.m_width_mask,
		     VSScaledParameterLibraryWriter::SPS_WIDTH_EXPECTED,
		     false, m_data[idata]->zenith_deg, 
		     m_data[idata]->zenith_deg,
		     false, m_data[idata]->azimuth_deg, 
		     m_data[idata]->azimuth_deg,
		     true,0,
		     false, m_data[idata]->ped_dev,m_data[idata]->ped_dev);

      wlibrary.write(m_data[idata]->m_array_lt.m_width_rms,
		     m_data[idata]->m_array_lt.m_width_mask,
		     VSScaledParameterLibraryWriter::SPS_WIDTH_RMS,
		     false, m_data[idata]->zenith_deg, 
		     m_data[idata]->zenith_deg,
		     false, m_data[idata]->azimuth_deg, 
		     m_data[idata]->azimuth_deg,
		     true,0,
		     false, m_data[idata]->ped_dev,m_data[idata]->ped_dev);

      wlibrary.write(m_data[idata]->m_array_lt.m_length_exp,
		     m_data[idata]->m_array_lt.m_length_mask,
		     VSScaledParameterLibraryWriter::SPS_LENGTH_EXPECTED,
		     false, m_data[idata]->zenith_deg, 
		     m_data[idata]->zenith_deg,
		     false, m_data[idata]->azimuth_deg, 
		     m_data[idata]->azimuth_deg,
		     true,0,
		     false, m_data[idata]->ped_dev,m_data[idata]->ped_dev);

      wlibrary.write(m_data[idata]->m_array_lt.m_length_rms,
		     m_data[idata]->m_array_lt.m_length_mask,
		     VSScaledParameterLibraryWriter::SPS_LENGTH_RMS,
		     false, m_data[idata]->zenith_deg, 
		     m_data[idata]->zenith_deg,
		     false, m_data[idata]->azimuth_deg, 
		     m_data[idata]->azimuth_deg,
		     true,0,
		     false, m_data[idata]->ped_dev,m_data[idata]->ped_dev);
      
      wlibrary.write(m_data[idata]->m_array_lt.m_disp_exp,
		     m_data[idata]->m_array_lt.m_disp_mask,
		     VSScaledParameterLibraryWriter::SPS_DISP_EXPECTED,
		     false, m_data[idata]->zenith_deg, 
		     m_data[idata]->zenith_deg,
		     false, m_data[idata]->azimuth_deg, 
		     m_data[idata]->azimuth_deg,
		     true,0,
		     false, m_data[idata]->ped_dev,m_data[idata]->ped_dev);

      wlibrary.write(m_data[idata]->m_array_lt.m_disp_rms,
		     m_data[idata]->m_array_lt.m_disp_mask,
		     VSScaledParameterLibraryWriter::SPS_DISP_RMS,
		     false, m_data[idata]->zenith_deg, 
		     m_data[idata]->zenith_deg,
		     false, m_data[idata]->azimuth_deg, 
		     m_data[idata]->azimuth_deg,
		     true,0,
		     false, m_data[idata]->ped_dev,m_data[idata]->ped_dev);

      
      const unsigned nscope =  m_data[idata]->m_sp_tel.size();
      for(unsigned iscope = 0; iscope < nscope; iscope++)
	{
	  if(!m_data[idata]->m_sp_tel[iscope]) continue;

	  wlibrary.write(m_data[idata]->m_scope_lt[iscope].m_width_exp,
			 m_data[idata]->m_scope_lt[iscope].m_width_mask,
			 VSScaledParameterLibraryWriter::SPS_WIDTH_EXPECTED,
			 false, m_data[idata]->zenith_deg, 
			 m_data[idata]->zenith_deg,
			 false, m_data[idata]->azimuth_deg, 
			 m_data[idata]->azimuth_deg,
			 false,iscope,
			 false, m_data[idata]->ped_dev,m_data[idata]->ped_dev);
	  
	  wlibrary.write(m_data[idata]->m_scope_lt[iscope].m_width_rms,
			 m_data[idata]->m_scope_lt[iscope].m_width_mask,
			 VSScaledParameterLibraryWriter::SPS_WIDTH_RMS,
			 false, m_data[idata]->zenith_deg, 
			 m_data[idata]->zenith_deg,
			 false, m_data[idata]->azimuth_deg, 
			 m_data[idata]->azimuth_deg,
			 false,iscope,
			 false, m_data[idata]->ped_dev,m_data[idata]->ped_dev);

	  wlibrary.write(m_data[idata]->m_scope_lt[iscope].m_length_exp,
			 m_data[idata]->m_scope_lt[iscope].m_length_mask,
			 VSScaledParameterLibraryWriter::SPS_LENGTH_EXPECTED,
			 false, m_data[idata]->zenith_deg, 
			 m_data[idata]->zenith_deg,
			 false, m_data[idata]->azimuth_deg, 
			 m_data[idata]->azimuth_deg,
			 false,iscope,
			 false, m_data[idata]->ped_dev,m_data[idata]->ped_dev);

	  wlibrary.write(m_data[idata]->m_scope_lt[iscope].m_length_rms,
			 m_data[idata]->m_scope_lt[iscope].m_length_mask,
			 VSScaledParameterLibraryWriter::SPS_LENGTH_RMS,
			 false, m_data[idata]->zenith_deg, 
			 m_data[idata]->zenith_deg,
			 false, m_data[idata]->azimuth_deg, 
			 m_data[idata]->azimuth_deg,
			 false,iscope,
			 false, m_data[idata]->ped_dev,m_data[idata]->ped_dev);
      
	  wlibrary.write(m_data[idata]->m_scope_lt[iscope].m_disp_exp,
			 m_data[idata]->m_scope_lt[iscope].m_disp_mask,
			 VSScaledParameterLibraryWriter::SPS_DISP_EXPECTED,
			 false, m_data[idata]->zenith_deg, 
			 m_data[idata]->zenith_deg,
			 false, m_data[idata]->azimuth_deg, 
			 m_data[idata]->azimuth_deg,
			 false,iscope,
			 false, m_data[idata]->ped_dev,m_data[idata]->ped_dev);

	  wlibrary.write(m_data[idata]->m_scope_lt[iscope].m_disp_rms,
			 m_data[idata]->m_scope_lt[iscope].m_disp_mask,
			 VSScaledParameterLibraryWriter::SPS_DISP_RMS,
			 false, m_data[idata]->zenith_deg, 
			 m_data[idata]->zenith_deg,
			 false, m_data[idata]->azimuth_deg, 
			 m_data[idata]->azimuth_deg,
			 false,iscope,
			 false, m_data[idata]->ped_dev,m_data[idata]->ped_dev);
	}


    }


  if(!m_options.no_diagnostics)
    {
      VSOctaveH5WriterCellVector* wc = writer->writeCellVector("data", ndata);
      vsassert(wc);

      for(unsigned idata = 0; idata < ndata; idata++)
	{ 
	  VSOctaveH5WriterStruct* ws = wc->writeStruct(idata);
	  vsassert(ws);  

	  ws->writeScalar("zenith_deg",m_data[idata]->zenith_deg);
	  ws->writeScalar("azimuth_deg",m_data[idata]->azimuth_deg);
	  ws->writeScalar("ped_dev",m_data[idata]->ped_dev);

	  m_data[idata]->m_array_lt.save(ws);      

	  delete ws;
	}

      delete wc;      
    }

}

#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSScaledParameterLibraryVisitor::configure(VSOptions& options,
						const std::string& profile, 
						const std::string& opt_prefix)
{
  std::vector< std::string > sim_energy_weight;
  sim_energy_weight.push_back("powerlaw");
  sim_energy_weight.push_back("2.5");

  VSSimEnergyWeightCalc::getDefaultOptions().sim_energy_weight =
    sim_energy_weight;

  VSSimEnergyWeightCalc::configure(options);
  VSCutsEvaluator::configure(options);

  options.findWithValue(OPTNAME(opt_prefix,"min_events"), 
			s_default_options.min_count,
			"Minimum number of events per cell.");  

  options.findBoolValue(OPTNAME(opt_prefix,"median"), 
			s_default_options.median, true,
			"Calculate tables with median and 68% containment "
			"intervals rather than mean and variance.");

  options.findWithValue("theta_cut", s_default_options.theta_cut,
			"Cut on the maximum separation in degrees between "
			"the true and reconstructed primary direction.");

  options.findBoolValue(OPTNAME(opt_prefix,"no_scope_tables"), 
			s_default_options.no_scope_tables, 
			true,
			"Do not generate the per-telescope lookup tables.");

  options.findBoolValue("no_diagnostics", 
			s_default_options.no_diagnostics, true,
			"Do not save extra diagnostic data "
			"associated with lookup table generation.");

  options.findWithValue(OPTNAME(opt_prefix,"ndim"), s_default_options.ndim,
			"Set the number of dimensions in the SP tables. "
			"Should be 2 or 3.");  
}

// ============================================================================
// VSEffectiveAreaLibraryVisitor
// ============================================================================
VSEffectiveAreaLibraryVisitor::Options 
VSEffectiveAreaLibraryVisitor::s_default_options = 
VSEffectiveAreaLibraryVisitor::Options();

VSEffectiveAreaLibraryVisitor::Data::~Data()
{
  delete m_sp_calc;
  delete m_egy_calc;
}

VSEffectiveAreaLibraryVisitor::
VSEffectiveAreaLibraryVisitor(const Options& opt):
  VSEventDataVisitor(), m_cuts(), 
  m_data(), m_data_ptr(),
  m_sim(), m_evt(), m_sim_header(),
  m_options(opt)
{
  m_cuts = new VSCutsEvaluator;
}
 
VSEffectiveAreaLibraryVisitor::~VSEffectiveAreaLibraryVisitor()
{
  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++) delete m_data[idata];

  delete m_cuts;
}

void VSEffectiveAreaLibraryVisitor::
visitRun(const VSAnalysisStage1Data& stage1,
	 const VSTargetTable::Observation& obs,
	 const VSArrayMergedCalibrationData& cal)
{
  // Determine Offset and Pointing --------------------------------------------
  m_zn_deg = stage1.run_info.zn_mean_deg;
  m_az_deg = stage1.run_info.az_mean_deg;
  m_ped_dev = cal.mean_scaled_dev;
  m_offset_deg = stage1.sim_info->wobble_theta_deg;
  
  std::cout << "OFFSET  " << std::setw(20) << m_offset_deg 
	    << " ZN     " << std::setw(20) << m_zn_deg 
	    << " AZ     " << std::setw(20) << m_az_deg 
	    << " PEDDEV " << std::setw(20) << m_ped_dev << std::endl;

  // Create EffectiveAreaCalc -------------------------------------------------
  m_data_ptr = NULL;

  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    {
      SphericalCoords azel =
	SphericalCoords::makeDeg(stage1.run_info.zn_mean_deg,
				 stage1.run_info.az_mean_deg);

      double dzn = fabs(m_zn_deg - m_data[idata]->zenith_deg);
      if(dzn < 0.5) m_zn_deg = m_data[idata]->zenith_deg;

      double dtheta = azel.separation(m_data[idata]->m_azel).deg();
      
      if(dtheta < 0.5) m_az_deg = m_data[idata]->azimuth_deg;

      double dped =
	fabs(cal.mean_scaled_dev - m_data[idata]->ped_dev);
      
      if(dtheta < 0.5 && dped < 0.1)
	{
	  m_data_ptr = m_data[idata]; 
	  break;
	}
    }

  if(m_data_ptr == NULL)
    {
      m_data_ptr = new Data(m_zn_deg,m_az_deg,m_ped_dev);
      m_data.push_back(m_data_ptr);
    }
  
  m_data_ptr->m_effarea_calc.setEnergyBinning(m_log10_egybin,
					      m_log10_egylo,
					      m_log10_egyhi);
  m_data_ptr->m_effarea_trigger_calc.setEnergyBinning(m_log10_egybin,
					      m_log10_egylo,
					      m_log10_egyhi);
  m_data_ptr->m_kernel_calc.setEnergyBinning(m_log10_egybin,
					     m_log10_egylo,
					     m_log10_egyhi);
  m_data_ptr->m_psf_calc.setEnergyBinning(m_log10_egybin,
					  m_log10_egylo,
					  m_log10_egyhi);

  m_data_ptr->m_effarea_calc.loadSimInfo(stage1.sim_info);
  m_data_ptr->m_effarea_trigger_calc.loadSimInfo(stage1.sim_info);
  m_data_ptr->m_kernel_calc.loadSimInfo(stage1.sim_info);
  m_data_ptr->m_psf_calc.loadSimInfo(stage1.sim_info);

  m_data_ptr->m_effarea_calc.loadHeader(m_sim_header);
  m_data_ptr->m_effarea_trigger_calc.loadHeader(m_sim_header);


  m_data_ptr->m_sp_calc = new VSScaledParameterCalc;
  m_data_ptr->m_sp_calc->load(stage1.run_info.nchan,
			      stage1.run_info.zn_mean_deg,
			      stage1.run_info.az_mean_deg,cal.mean_scaled_dev);

  m_data_ptr->m_egy_calc = new VSEnergyCalcLT;
  m_data_ptr->m_egy_calc->load(stage1.run_info.nchan,
			       stage1.run_info.zn_mean_deg,
			       stage1.run_info.az_mean_deg,
			       cal.mean_scaled_dev);
}

void VSEffectiveAreaLibraryVisitor::leaveRun()
{

}

void VSEffectiveAreaLibraryVisitor::visitEvent(const VSEventArrayDatum& evt)
{
  m_evt = evt;
  m_data_ptr->m_sp_calc->calcSP(m_evt);
  bool has_mlt_energy = m_data_ptr->m_egy_calc->calcEnergy(m_evt);


  if(!std::isfinite(m_evt.mlt_log10_energy) || 
     !std::isfinite(m_evt.theta1) || 
     m_evt.theta1 > m_options.theta_cut) return;

  if(m_cuts->isSelected(m_evt)) 
    {  
      m_data_ptr->m_effarea_calc.accumulate(m_sim.energy_tev);

      if(m_evt.theta1 < m_options.kernel_theta_cut && has_mlt_energy)
	m_data_ptr->m_kernel_calc.accumulate(m_sim.energy_tev,
					     m_evt.mlt_log10_energy);

      m_data_ptr->m_psf_calc.accumulate(m_sim.energy_tev,
					std::pow(m_evt.theta1,2));
    }
}

void VSEffectiveAreaLibraryVisitor::leaveEvent()
{

}

void VSEffectiveAreaLibraryVisitor::
visitSimEvent(const VSArraySimulationDatum& sim)
{  
  m_sim = sim;
  m_data_ptr->m_effarea_trigger_calc.accumulate(m_sim.energy_tev);
}

void VSEffectiveAreaLibraryVisitor::leaveSimEvent()
{

}

void VSEffectiveAreaLibraryVisitor::
visitSimHeader(const VSHeaderSimulationDatum& header)
{ 
  m_sim_header = header;

  double log10_egylo = log10(header.tables.front().energy_tev);
  double log10_egyhi = log10(header.tables.back().energy_tev);
  std::set< double > egy_set;
  unsigned negy = 0;

  for(std::vector< VSTableSimulationDatum >::const_iterator itr = 
	header.tables.begin(); itr != header.tables.end(); ++itr)
    {
      double log10_egy = std::log10(itr->energy_tev);

      if(log10_egy < log10_egylo) log10_egylo = log10_egy;
      if(log10_egy > log10_egyhi) log10_egyhi = log10_egy;


      std::set< double >::iterator itr = egy_set.begin();

      for(; itr != egy_set.end(); ++itr)
	if(fabs(log10_egy-*itr) < 0.001) break;

      if(itr == egy_set.end())
	{
	  egy_set.insert(log10_egy);
	  negy++;
	}
    }

  if(negy > 1)
    {
      m_log10_egybin = (log10_egyhi-log10_egylo)/(double)(negy-1);
      m_log10_egylo  = log10_egylo - m_log10_egybin/2.;
      m_log10_egyhi  = log10_egyhi + m_log10_egybin/2.;
    }
}

void VSEffectiveAreaLibraryVisitor::leaveSimHeader()
{

}

void VSEffectiveAreaLibraryVisitor::createLibrary()
{
  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    {
      std::cout << __PRETTY_FUNCTION__  
		<< ": Generating Instrument Response Fn " 
		<< " Zn  " << std::setw(15) << m_data[idata]->zenith_deg
		<< " Az  " << std::setw(15) << m_data[idata]->azimuth_deg
		<< " Ped " << std::setw(15) << m_data[idata]->ped_dev
		<< std::endl;

      m_data[idata]->m_effarea_calc.calcEffarea();
      m_data[idata]->m_effarea_trigger_calc.calcEffarea();
      m_data[idata]->m_effarea_calc.fit();
      m_data[idata]->m_effarea_trigger_calc.fit();
      m_data[idata]->m_kernel_calc.calcKernel();
      m_data[idata]->m_kernel_calc.fit();
      m_data[idata]->m_psf_calc.fit();
    }
}

void VSEffectiveAreaLibraryVisitor::save(VSOctaveH5WriterStruct* writer) const
{
  std::cout << "Writing Library File..." << std::endl;

  typedef VSScaledParameterLibraryWriter SPL;

  SPL wlibrary(writer);

  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    {
      double zn = m_data[idata]->zenith_deg;
      double az = m_data[idata]->azimuth_deg;
      double ped = m_data[idata]->ped_dev;

      wlibrary.write(m_data[idata]->m_effarea_calc.getEffarea(),
		     SPL::SPS_EFFECTIVE_AREA,
		     false, zn, zn, false, az, az, true, 0, false, ped, ped);

      wlibrary.write(m_data[idata]->m_kernel_calc.getSigma1(),
		     SPL::SPS_KERNEL_SIGMA1,
		     false, zn, zn, false, az, az, true, 0, false, ped, ped);
  
      wlibrary.write(m_data[idata]->m_kernel_calc.getBias1(),
		     SPL::SPS_KERNEL_BIAS1,
		     false, zn, zn, false, az, az, true, 0, false, ped, ped);

      wlibrary.write(m_data[idata]->m_kernel_calc.getSigma2(),
		     SPL::SPS_KERNEL_SIGMA2,
		     false, zn, zn, false, az, az, true, 0, false, ped, ped);

      wlibrary.write(m_data[idata]->m_kernel_calc.getBias2(),
		     SPL::SPS_KERNEL_BIAS2,
		     false, zn, zn, false, az, az, true, 0, false, ped, ped);

      wlibrary.write(m_data[idata]->m_kernel_calc.getAlpha(),
		     SPL::SPS_KERNEL_ALPHA,
		     false, zn, zn, false, az, az, true, 0, false, ped, ped);

      wlibrary.write(m_data[idata]->m_psf_calc.getSigma1(),
		     VSScaledParameterLibraryWriter::SPS_PSF_SIGMA1,
		     false, zn, zn, false, az, az, true, 0, false, ped, ped);

      wlibrary.write(m_data[idata]->m_psf_calc.getSigma2(),
		     VSScaledParameterLibraryWriter::SPS_PSF_SIGMA2,
		     false, zn, zn, false, az, az, true, 0, false, ped, ped);

      wlibrary.write(m_data[idata]->m_psf_calc.getAlpha(),
		     VSScaledParameterLibraryWriter::SPS_PSF_ALPHA,
		     false, zn, zn, false, az, az, true, 0, false, ped, ped);
    }

  if(!m_options.no_diagnostics)
    {
      VSOctaveH5WriterCellVector* wc = writer->writeCellVector("data", ndata);
      vsassert(wc);

      for(unsigned idata = 0; idata < ndata; idata++)
	{
	  VSOctaveH5WriterStruct* ws = wc->writeStruct(idata);
	  vsassert(ws);  
      
	  ws->writeScalar("zenith_deg",m_data[idata]->zenith_deg);
	  ws->writeScalar("azimuth_deg",m_data[idata]->azimuth_deg);
	  ws->writeScalar("ped_dev",m_data[idata]->ped_dev);

	  m_data[idata]->m_effarea_calc.save(ws->writeStruct("effarea"));
	  m_data[idata]->
	    m_effarea_trigger_calc.save(ws->writeStruct("effarea_trigger"));
	  m_data[idata]->m_kernel_calc.save(ws->writeStruct("kernel"));
	  m_data[idata]->m_psf_calc.save(ws->writeStruct("psf"));
	  
	  m_data[idata]->m_sp_calc->save(ws->writeStruct("sp_calc"));
	  m_data[idata]->m_egy_calc->save(ws->writeStruct("egy_calc"));

	  delete ws;
	}

      delete wc;
    }
}

void VSEffectiveAreaLibraryVisitor::configure(VSOptions& options,
					      const std::string& profile, 
					      const std::string& opt_prefix)
{
  VSScaledParameterCalc::configure(options,profile,opt_prefix);
  VSEnergyCalcLT::configure(options,profile,opt_prefix);
  VSCutsEvaluator::configure(options);

  options.findWithValue(OPTNAME(opt_prefix,"kernel_theta_cut"), 
			s_default_options.kernel_theta_cut,
			"Cut on the maximum separation in degrees "
			"between the true and reconstructed primary "
			"direction used to select the sample of events "
			"with which to model the energy response function. "
			"This option should be close to the theta_cut which "
			"will be used for spectral analysis.");  

  options.findWithValue(OPTNAME(opt_prefix,"theta_cut"), 
			s_default_options.theta_cut,
			"Quality cut on the maximum separation in degrees "
			"between the true and reconstructed primary "
			"direction.  This cut is intended to exclude badly "
			"reconstructed events that fall out of the main psf "
			"distribution.");  

  options.findBoolValue(OPTNAME(opt_prefix,"no_diagnostics"), 
			s_default_options.no_diagnostics, true,
			"Do not save extra diagnostic data "
			"associated with lookup table generation.");
}
