//-*-mode:c++; mode:font-lock;-*-

/*! \file VSReconstruction.cpp
  Implementation of VVV "in the sky" chi-2 shower core reconstruction
  algorithm

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/16/2005
*/

#include<cmath>
#include<cfloat>
#include<cassert>

#include <Vec3D.hpp>

#include "VSWeightedStat.hpp"

#include "VSReconstruction.hpp"

using namespace VERITAS;
using namespace Physics;

VSReconstruction::
VSReconstruction(const ScopeSpecs& scopes, double curvature_distance):
  fScopes(), fCurvatureDistance(curvature_distance),
  ec_x(), ec_y(), ec_z(), rc_x(), rc_y(), rc_z(), chi2(), 
  nscopes(),  nrays(), weight(),
  scope_array_size(), scope_nray(), ray_weight(), 
  es_x(), es_y(), es_z(), scope_id(), r_x(), r_y(), r_z(), 
  scope_ray_weight(), scope_weighting(), 
  Q_scope_xx(), Q_scope_yy(), Q_scope_zz(), 
  Q_scope_xy(), Q_scope_xz(), Q_scope_yz(), 
  ray_array_size(), 
  es_dot_ec(), one_minus_es_dot_ec_squared(), q_tau_x(), q_tau_y(), q_tau_z(), 
  q_p_x(), q_p_y(), q_p_z(), Q_xx(), Q_yy(), Q_zz(), Q_xy(), Q_xz(), Q_yz()
{
  setupArray(scopes);
}

VSReconstruction::~VSReconstruction()
{
  deleteChi2Workspace();
  deleteArray();
}

Vec3D VSReconstruction::findRotationVector(const Axis& optical_axis) const
{
  double oa_x = optical_axis.cx;
  double oa_y = optical_axis.cy;

  if((fabs(oa_x)<=DBL_EPSILON)&&(fabs(oa_y)<=DBL_EPSILON))
    return Physics::Vec3D(0,0,0);

  double oa_z = -sqrt(1.0-oa_x*oa_x-oa_y*oa_y);

#if 0  
  std::cout << oa_x << ' ' << oa_y << ' ' << oa_z << std::endl;
#endif

  Physics::Vec3D oa_rot(0,0,-1); 
  Vec3D oa(oa_x,oa_y,oa_z);
  double cos_theta = oa_rot*oa;
  oa_rot = oa_rot^oa;
  double sin_theta = oa_rot.Norm();
  oa_rot *= atan2(sin_theta,cos_theta)/sin_theta;
  return oa_rot;
}

double VSReconstruction::
distanceToAxis(const Axis& optical_axis, 
	       double ec_x, double ec_y, double ec_z) const
{
  double oa_x = optical_axis.cx;
  double oa_y = optical_axis.cy;
  double oa_z = -sqrt(1.0-oa_x*oa_x-oa_y*oa_y);

  Vec3D ec(ec_x,ec_y,ec_z);
  Vec3D oa(oa_x,oa_y,oa_z);
  double cos_theta = ec*oa;
  double sin_theta = Vec3D(ec^oa).Norm();
  return atan2(sin_theta,cos_theta);
}

template<typename T> static void delete_and_zero(T*&x)
{
  delete[] x;
  x=0;
}

void VSReconstruction::deleteArray()
{
  for(unsigned iscope=0; iscope<scope_array_size; iscope++)
    {
      delete_and_zero(ray_weight[iscope]);
      delete_and_zero(es_x[iscope]);
      delete_and_zero(es_y[iscope]);
      delete_and_zero(es_z[iscope]);
    }
  
  delete_and_zero(scope_nray);
  delete_and_zero(ray_weight);
  delete_and_zero(es_x);
  delete_and_zero(es_y);
  delete_and_zero(es_z);
  delete_and_zero(scope_id);
  delete_and_zero(r_x);
  delete_and_zero(r_y);
  delete_and_zero(r_z);
  delete_and_zero(scope_ray_weight);
  delete_and_zero(scope_weighting);
  delete_and_zero(Q_scope_xx);
  delete_and_zero(Q_scope_yy);
  delete_and_zero(Q_scope_zz);
  delete_and_zero(Q_scope_xy);
  delete_and_zero(Q_scope_xz);
  delete_and_zero(Q_scope_yz);
  
  scope_array_size = 0;
}

void VSReconstruction::setupArray(const ScopeSpecs& scopes)
{
  fScopes=scopes;

  if(scopes.size() > scope_array_size)
    {
      deleteArray();

      scope_array_size = scopes.size();

      scope_nray = new unsigned[scope_array_size];
      ray_weight = new double*[scope_array_size];
      es_x = new double*[scope_array_size];
      es_y = new double*[scope_array_size];
      es_z = new double*[scope_array_size];
      scope_id = new unsigned[scope_array_size];
      r_x = new double[scope_array_size];
      r_y = new double[scope_array_size];
      r_z = new double[scope_array_size];
      scope_ray_weight = new double[scope_array_size];
      scope_weighting = new double[scope_array_size];
      Q_scope_xx = new double[scope_array_size];
      Q_scope_yy = new double[scope_array_size];
      Q_scope_zz = new double[scope_array_size];
      Q_scope_xy = new double[scope_array_size];
      Q_scope_xz = new double[scope_array_size];
      Q_scope_yz = new double[scope_array_size];

      if(ray_array_size)
	for(unsigned iscope=0; iscope<scope_array_size; iscope++)
	  {
	    ray_weight[iscope] = new double[ray_array_size];
	    es_x[iscope] = new double[ray_array_size];
	    es_y[iscope] = new double[ray_array_size];
	    es_z[iscope] = new double[ray_array_size];
	  }
      else
	for(unsigned iscope=0; iscope<scope_array_size; iscope++)
	  {
	    ray_weight[iscope] = 0;
	    es_x[iscope] = 0;
	    es_y[iscope] = 0;
	    es_z[iscope] = 0;
	  }
    }
}

void VSReconstruction::deleteChi2Workspace()
{
  for(unsigned iscope=0; iscope<scope_array_size; iscope++)
    {
      delete_and_zero(ray_weight[iscope]);
      delete_and_zero(es_x[iscope]);
      delete_and_zero(es_y[iscope]);
      delete_and_zero(es_z[iscope]);
    }

  delete_and_zero(es_dot_ec);
  delete_and_zero(one_minus_es_dot_ec_squared);
  delete_and_zero(q_tau_x);
  delete_and_zero(q_tau_y);
  delete_and_zero(q_tau_z);
  delete_and_zero(q_p_x);
  delete_and_zero(q_p_y);
  delete_and_zero(q_p_z);
  delete_and_zero(Q_xx);
  delete_and_zero(Q_yy);
  delete_and_zero(Q_zz);
  delete_and_zero(Q_xy);
  delete_and_zero(Q_xz);
  delete_and_zero(Q_yz);

  ray_array_size = 0;
}

void VSReconstruction::setupChi2Workspace(const AllScopeRays& rays)
{
  if(rays.size() != fScopes.size())
    {
      std::cerr << "rays.size()=" << rays.size() << std::endl
		<< "fScopes.size()=" << fScopes.size() << std::endl;
      assert(rays.size() == fScopes.size());
    }

  nscopes = 0;
  nrays = 0;
  unsigned max_scope_nrays = 0;

  for(unsigned iscope=0; iscope<rays.size(); iscope++)
    {
      const unsigned nray = rays[iscope].size();
      if(nray)
	{
	  nscopes++;  
	  nrays+=nray;
	  if(nray>max_scope_nrays)max_scope_nrays=nray;
	}
    }

  if(ray_array_size < max_scope_nrays)
    {
      deleteChi2Workspace();

      ray_array_size = max_scope_nrays;

      for(unsigned iscope=0; iscope<scope_array_size; iscope++)
	{
	  ray_weight[iscope] = new double[ray_array_size];
	  es_x[iscope] = new double[ray_array_size];
	  es_y[iscope] = new double[ray_array_size];
	  es_z[iscope] = new double[ray_array_size];
	}
      
      es_dot_ec = new double[ray_array_size];
      one_minus_es_dot_ec_squared = new double[ray_array_size];
      q_tau_x = new double[ray_array_size];
      q_tau_y = new double[ray_array_size];
      q_tau_z = new double[ray_array_size];
      q_p_x = new double[ray_array_size];
      q_p_y = new double[ray_array_size];
      q_p_z = new double[ray_array_size];
      Q_xx = new double[ray_array_size];
      Q_yy = new double[ray_array_size];
      Q_zz = new double[ray_array_size];
      Q_xy = new double[ray_array_size];
      Q_xz = new double[ray_array_size];
      Q_yz = new double[ray_array_size];
    }

  weight = 0;
  for(unsigned iscope=0, jscope=0; iscope<rays.size(); iscope++)
    {
      const unsigned nray = rays[iscope].size();
      if(nray)
	{
	  scope_ray_weight[jscope] = 0;
	  scope_nray[jscope] = nray;
	  for(unsigned iray=0; iray<nray; iray++)
	    {
	      const double cx = rays[iscope][iray].cx;
	      const double cy = rays[iscope][iray].cy;
	      ray_weight[jscope][iray] = rays[iscope][iray].weight;
	      es_x[jscope][iray] = cx;
	      es_y[jscope][iray] = cy;
	      es_z[jscope][iray] = -sqrt(1.0-cx*cx-cy*cy);
	      scope_ray_weight[jscope] += rays[iscope][iray].weight;
	    }
	  scope_id[jscope] = iscope;
	  r_x[jscope] = fScopes[iscope].scope_x;
	  r_y[jscope] = fScopes[iscope].scope_y;
	  r_z[jscope] = fScopes[iscope].scope_z;
	  weight += scope_ray_weight[jscope];
	  jscope++;
	}
    }
}

class NullMapOutputIterator
{
public:
  NullMapOutputIterator& operator++() { return *this; }
  NullMapOutputIterator& operator++(int) { return *this; }
  NullMapOutputIterator& operator*() { return *this; }
  NullMapOutputIterator& operator=(const VSReconstruction::Chi2MapElement&) { return *this; }
};

void VSReconstruction::calculateCurvature(MinimizationResults& res)
{
  Physics::Vec3D oa_rot(findRotationVector(Axis(res.ec_x,res.ec_y)));

  double curv[9];

  unsigned i=0;
  for(int x_bin = -1; x_bin<=1; x_bin++)
    {
      double x = double(x_bin)*fCurvatureDistance;
      for(int y_bin = -1; y_bin<=1; y_bin++)
	{	
	  if((y_bin==0)&&(x_bin==0)){ i++; continue; }
	  double y = double(y_bin)*fCurvatureDistance;
	  double r = sqrt(x*x+y*y);

	  Physics::Vec3D ec(0,0,1);
	  ec.Rotate(Physics::Vec3D(r*M_PI/180.0,0,0));
	  ec.Rotate(Physics::Vec3D(0,0,atan2(y,x)));
	  ec.Rotate(oa_rot);

	  // Transfer to calculateChi2() variables
	  ec_x = ec.x;
	  ec_y = ec.y;
	  ec_z = ec.z;

	  calculateChi2ForHypothesis(true);

	  curv[i]=chi2;
#if 0
	  std::cerr << i << ' ' << curv[i] << std::endl;
#endif
	  i++;
	}
    }

  double norm = 1.0/(res.chi2*fCurvatureDistance*fCurvatureDistance);
  double chi2_pp_xy = (curv[0] - curv[2] - curv[6] + curv[8])*norm/4.0;
  double chi2_pp_xx = (curv[3] + curv[5] - 2.0*res.chi2)*norm;
  double chi2_pp_yy = (curv[1] + curv[7] - 2.0*res.chi2)*norm;

#if 0
  std::cerr << norm << std::endl;
  std::cerr << "xy " << chi2_pp_xy << std::endl;
  std::cerr << "xx " << chi2_pp_xx << std::endl;
  std::cerr << "yy " << chi2_pp_yy << std::endl;
#endif

#define SQR(x) ((x)*(x)) // Bad! Uses x twice

  double diff2 = SQR(chi2_pp_xx-chi2_pp_yy)+4.0*SQR(chi2_pp_xy);

  if(diff2<DBL_EPSILON)
    {
      res.chi2_curv_l = 1.0/sqrt((chi2_pp_xx+chi2_pp_yy)/2.0);
      res.chi2_curv_w = 1.0/sqrt((chi2_pp_xx+chi2_pp_yy)/2.0);
      res.chi2_curv_ang = 0;
    }
  else
    {
      double diff = sqrt(diff2);
      res.chi2_curv_l = 1.0/sqrt((chi2_pp_xx+chi2_pp_yy-diff)/2.0);
      res.chi2_curv_w = 1.0/sqrt((chi2_pp_xx+chi2_pp_yy+diff)/2.0);
      res.chi2_curv_ang = 
	atan2(-2.0*chi2_pp_xy, chi2_pp_xx-chi2_pp_yy)/2.0;
    }
}

void VSReconstruction::
calculatePProfile(MinimizationResults& res)
{
  // Assume workspace has been configured by "setupChi2Workspace()"
  
  res.ps.clear();
  res.ps.reserve(res.nrays);
  
  VSWeightedStat<double> ps_stats;
  res.scope_id.resize(nscopes);
  res.scope_ray_weight.resize(nscopes);
  res.scope_theta_stat.resize(nscopes);
  res.scope_ps_stat.resize(nscopes);
  for(unsigned iscope=0; iscope<nscopes;iscope++)
    {
      VSWeightedStat<double> scope_theta_stats;
      VSWeightedStat<double> scope_ps_stats;

      const unsigned nray = scope_nray[iscope];
            
      const double *const ws_p = ray_weight[iscope];
      const double *const es_x_p = es_x[iscope];
      const double *const es_y_p = es_y[iscope];
      const double *const es_z_p = es_z[iscope];

      for(unsigned iray=0; iray<nray; iray++)
	{
	  const double my_es_dot_ec = 
	    res.ec_x*es_x_p[iray] 
	    + res.ec_y*es_y_p[iray] 
	    + res.ec_z*es_z_p[iray];
	  
	  const double my_one_minus_es_dot_ec_squared =
	    1.0 - my_es_dot_ec*my_es_dot_ec;
      
#if 0
	  if(my_one_minus_es_dot_ec_squared < DBL_EPSILON)
	    std::cout << iscope << '/' << iray << ' '
		      << res.ec_x << " * " << es_x_p[iray] << " + "
		      << res.ec_y << " * " << es_y_p[iray] << " + "
		      << res.ec_z << " * " << es_z_p[iray] << " = "
		      << my_es_dot_ec << ' ' 
		      << my_one_minus_es_dot_ec_squared << ' '
		      << 1.0/my_one_minus_es_dot_ec_squared 
		      <<std::endl;
#endif
      
	  if(my_one_minus_es_dot_ec_squared >= DBL_EPSILON)
	    {
	      const double proj_x = 
		(res.ec_x - my_es_dot_ec * es_x_p[iray])/
		my_one_minus_es_dot_ec_squared;
	      const double proj_y = 
		(res.ec_y - my_es_dot_ec * es_y_p[iray])/
		my_one_minus_es_dot_ec_squared;
	      const double proj_z = 
		(res.ec_z - my_es_dot_ec * es_z_p[iray])/
		my_one_minus_es_dot_ec_squared;

	      const double del_r_x = r_x[iscope] - res.rc_x;
	      const double del_r_y = r_y[iscope] - res.rc_y;
	      const double del_r_z = r_z[iscope] - res.rc_z;

	      const double ps = 
		del_r_x*proj_x + del_r_y*proj_y + del_r_z*proj_z;

	      ps_stats.accumulate(ps,ws_p[iray]);
	      scope_ps_stats.accumulate(ps,ws_p[iray]);
	      scope_theta_stats.accumulate(acos(my_es_dot_ec),ws_p[iray]);
	      res.ps.push_back(std::make_pair(ps,ws_p[iray]));
	    }
	}
      if(nray)
	{
	  res.scope_id[iscope]         = scope_id[iscope];
	  res.scope_ray_weight[iscope] = scope_ray_weight[iscope];
	  res.scope_theta_stat[iscope] = 
	    MinimizationResults::Stat(scope_theta_stats.mean(), 
				      scope_theta_stats.var());
	  res.scope_ps_stat[iscope]    = 
	    MinimizationResults::Stat(scope_ps_stats.mean(), 
				      scope_ps_stats.var());
	}
    }

  res.ps_mean = ps_stats.mean();
  res.ps_var = ps_stats.var();
}

void VSReconstruction::
gridSearch(const AllScopeRays& rays, const Axis& optical_axis,
	   std::vector<GridSpecs> grid, MinimizationResults& res)
{
  NullMapOutputIterator out;
  gridSearch(rays, optical_axis, grid, res, out);
}

void VSReconstruction::
calculateChi2ForAxis(const AllScopeRays& rays, 
		     const MinimizationResults& axis, 
		     MinimizationResults& res)
{
  setupChi2Workspace(rays);

  // Transfer to calculateChi2() variables
  ec_x = axis.ec_x;
  ec_y = axis.ec_y;
  ec_z = axis.ec_z;
  rc_x = axis.rc_x;
  rc_y = axis.rc_y;
  rc_z = axis.rc_z;
  
  calculateChi2ForHypothesis(false);

  res.nrays=nrays;
  res.chi2 = chi2;
  res.ec_x = ec_x;
  res.ec_y = ec_y;
  res.ec_z = ec_z;
  res.rc_x = rc_x;
  res.rc_y = rc_y;
  res.rc_z = rc_z;

  // The curvature parameters make no sense in the context when 
  // chi^2 is not necessarily a minimum

  res.chi2_curv_l = 0.0;
  res.chi2_curv_w = 0.0;
  res.chi2_curv_ang = 0.0; 
}
